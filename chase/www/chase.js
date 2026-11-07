/* Author: Flavio Ballante, Chemical Biology Consortium Sweden (CBCS-KI), Karolinska Institutet - 2026 */
/* Acknowledgments: thanks to Dr. Francesco Massai, Chemical Biology Consortium Sweden (CBCS-KI), Karolinska Institutet, */
/* for valuable scientific input and feedback during the development of CHASE. */
/* CHASE browser glue: plotly rendering, structure tooltips, click delegation, JPEG export. */
(function () {
  'use strict';
  // bound: plots that already have their hover and click handlers.
  // mouse: where the mouse pointer was last, so we know where to put the tooltip.
  // nonce: goes up with every hover. If the server's tooltip arrives after the mouse has moved on to
  // another point, its number no longer matches and we ignore it.
  const bound = {};
  let mouse = { x: 0, y: 0 };
  let nonce = 0;
  let tip = null;

  // ---- structure tooltip: one floating box, created the first time it is needed
  function tipEl() {
    if (!tip) { tip = document.createElement('div'); tip.className = 'chase-tip'; tip.hidden = true; document.body.appendChild(tip); }
    return tip;
  }
  // Put the tooltip next to the mouse pointer. If it would stick out of the window, show it on the
  // other side of the pointer instead.
  function place() {
    const t = tipEl(); const r = t.getBoundingClientRect();
    let x = mouse.x + 18, y = mouse.y + 14;
    if (x + r.width > window.innerWidth - 8) x = mouse.x - r.width - 18;
    if (y + r.height > window.innerHeight - 8) y = Math.max(8, window.innerHeight - r.height - 8);
    t.style.left = Math.max(8, x) + 'px'; t.style.top = y + 'px';
  }
  function hideTip() { nonce++; tipEl().hidden = true; }
  document.addEventListener('mousemove', (e) => { mouse = { x: e.clientX, y: e.clientY }; if (tip && !tip.hidden) place(); }, { passive: true });

  // When the user hovers over or clicks a plot point, tell the server which point it was
  // (its customdata holds a well or compound number). The histogram has no customdata and is ignored.
  function bind(el, id) {
    el.on('plotly_hover', (ev) => {
      const p = ev.points && ev.points[0];
      if (!p || p.customdata === undefined || p.customdata === null || p.customdata < 0 || p.data.type === 'histogram') return;
      nonce++;
      Shiny.setInputValue('chase_hover', { plot: id, idx: p.customdata, nonce: nonce }, { priority: 'event' });
    });
    el.on('plotly_unhover', hideTip);
    el.on('plotly_click', (ev) => {
      const p = ev.points && ev.points[0];
      if (!p || p.customdata === undefined || p.customdata === null || p.customdata < 0 || p.data.type === 'histogram') return;
      Shiny.setInputValue('chase_click', { plot: id, idx: p.customdata }, { priority: 'event' });
    });
    el.addEventListener('mouseleave', hideTip);
  }

  // Messages from the server (sent with session.send_custom_message in app.py).
  function init() {
    // chase-plot: draw a figure into the box with id m.id, or show the text m.empty if there is nothing to plot.
    // Plotly.react updates the plot in place, so the user's zoom is kept.
    Shiny.addCustomMessageHandler('chase-plot', (m) => {
      const el = document.getElementById(m.id);
      if (!el) return;
      if (!m.fig) {
        if (bound[m.id]) { Plotly.purge(el); delete bound[m.id]; }
        el.innerHTML = m.empty ? '<div class="empty">' + m.empty + '</div>' : '';
        return;
      }
      if (!bound[m.id]) el.innerHTML = '';
      Plotly.react(el, m.fig.data, m.fig.layout, m.config || { responsive: true, displaylogo: false });
      if (!bound[m.id]) { bound[m.id] = true; bind(el, m.id); }
    });
    // chase-tip: the tooltip HTML for the point under the mouse. We only show it if the mouse is still on that point.
    Shiny.addCustomMessageHandler('chase-tip', (m) => {
      if (m.nonce !== nonce) return;
      const t = tipEl(); t.innerHTML = m.html; t.hidden = false; place();
    });
    Shiny.addCustomMessageHandler('chase-scroll', (m) => {
      const el = document.getElementById(m.id); if (el) el.scrollIntoView({ behavior: 'smooth', block: m.block || 'nearest' });
    });
    // chase-store: remember a value in this browser (used for 'Your name').
    Shiny.addCustomMessageHandler('chase-store', (m) => { try { localStorage.setItem(m.key, m.value); } catch (e) { /* storage unavailable */ } });
  }
  // Once Shiny is connected, send the name remembered from last time to the server.
  function restoreUser() {
    try { const u = localStorage.getItem('chase.user'); if (u) Shiny.setInputValue('chase_saved_user', u); } catch (e) { /* storage unavailable */ }
  }
  if (window.jQuery) jQuery(document).one('shiny:connected', restoreUser);
  // Set up the message handlers now if Shiny is ready, otherwise as soon as it connects.
  if (window.Shiny && Shiny.addCustomMessageHandler) init(); else if (window.jQuery) jQuery(document).one('shiny:connected', init);

  // Any element with data-chase="<action>" becomes the Shiny input chase_<action> when clicked.
  // The element's other data-* attributes are sent along with it.
  document.addEventListener('click', (e) => {
    const t = e.target.closest('[data-chase]');
    if (!t || t.disabled) return;
    // _n is the current time, so every click is a new value. Without it, clicking the same button twice
    // would only reach the server once.
    e.preventDefault();
    const payload = Object.assign({}, t.dataset); delete payload.chase; payload._n = Date.now();
    Shiny.setInputValue('chase_' + t.dataset.chase, payload, { priority: 'event' });
  });
  // The triage note is sent when the user leaves the note box after changing it.
  document.addEventListener('change', (e) => {
    if (e.target.id === 'tri_note') Shiny.setInputValue('chase_note', { key: e.target.dataset.key, note: e.target.value }, { priority: 'event' });
  });

  // ---- plot colours (each browser remembers its own)
  function readColors() { const o = {}; document.querySelectorAll('input.chase-color').forEach((i) => { o[i.dataset.key] = i.value; }); return o; }
  function writeColors(o) { document.querySelectorAll('input.chase-color').forEach((i) => { if (o[i.dataset.key]) i.value = o[i.dataset.key]; }); }
  function sendColors() { const o = readColors(); Shiny.setInputValue('chase_colors', o); try { localStorage.setItem('chase.colors', JSON.stringify(o)); } catch (e) { /* storage unavailable */ } }
  // A colour picker sends many events while it is being dragged. We wait until it has been still for
  // 120 ms and only then send the colours.
  let colorT;
  document.addEventListener('input', (e) => { if (e.target.classList && e.target.classList.contains('chase-color')) { clearTimeout(colorT); colorT = setTimeout(sendColors, 120); } });
  function restoreColors() {
    let saved = null; try { saved = JSON.parse(localStorage.getItem('chase.colors') || 'null'); } catch (e) { saved = null; }
    if (saved) { writeColors(saved); Shiny.setInputValue('chase_colors', readColors()); }
  }
  // Shiny announces that it is connected with a jQuery event, not a normal browser event,
  // so we have to listen through jQuery.
  if (window.jQuery) jQuery(document).one('shiny:connected', restoreColors);
  // chase-colors: the server sets the colours (after a reset, or when a session file is opened).
  function onColors(m) { writeColors(m.colors || {}); sendColors(); }
  if (window.Shiny && Shiny.addCustomMessageHandler) Shiny.addCustomMessageHandler('chase-colors', onColors);
  else if (window.jQuery) jQuery(document).one('shiny:connected', () => Shiny.addCustomMessageHandler('chase-colors', onColors));
  // Show or hide a panel (the colour pickers), and tell screen readers whether it is open.
  window.chaseToggle = function (id, btn) { const el = document.getElementById(id); if (!el) return; el.hidden = !el.hidden; if (btn) btn.setAttribute('aria-expanded', String(!el.hidden)); };

  // ---- unsaved triage decisions: warn before the tab is closed. Downloading the board counts as saving it.
  let unsaved = 0;
  function onDirty(m) { unsaved = m.n || 0; }
  if (window.Shiny && Shiny.addCustomMessageHandler) Shiny.addCustomMessageHandler('chase-dirty', onDirty);
  else if (window.jQuery) jQuery(document).one('shiny:connected', () => Shiny.addCustomMessageHandler('chase-dirty', onDirty));
  // If there are unsaved decisions, the browser shows its own 'Leave site?' question.
  window.addEventListener('beforeunload', (e) => { if (unsaved > 0) { e.preventDefault(); e.returnValue = ''; } });
  document.addEventListener('click', (e) => {
    // This listener runs in the 'capture' phase (the 'true' at the end), so it runs before the download starts.
    if (e.target.closest('#export_triage, #save_session')) { unsaved = 0; Shiny.setInputValue('chase_saved', Date.now(), { priority: 'event' }); }
  }, true);

  // ---- save a plot as JPEG at three times the screen resolution
  // (at least 900 x 520 pixels for the main plot, 760 x 460 for a plate map)
  window.chaseJPEG = function (id, filename) {
    const el = document.getElementById(id);
    if (!el || !el.data) return;
    Plotly.downloadImage(el, { format: 'jpeg', filename: filename, width: Math.max(900, el.clientWidth), height: Math.max(520, el.clientHeight + 40), scale: 3 });
  };
  window.chaseHeatJPEG = function (id, filename) {
    const el = document.getElementById(id);
    if (!el || !el.data) return;
    Plotly.downloadImage(el, { format: 'jpeg', filename: filename, width: Math.max(760, el.clientWidth), height: Math.max(460, el.clientHeight), scale: 3 });
  };
})();
