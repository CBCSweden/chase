"""The 'How to use' page.

Author: Flavio Ballante, Chemical Biology Consortium Sweden (CBCS-KI), Karolinska Institutet - 2026
Acknowledgments: thanks to Dr. Francesco Massai, Chemical Biology Consortium Sweden (CBCS-KI), Karolinska Institutet,
for valuable scientific input and feedback during the development of CHASE.
"""
from __future__ import annotations

import importlib.metadata as md
import platform
import sys

from . import __version__, __date__

# The full name with the letters of the acronym in bold. It is shown in the top bar and on the help page.
ACRO = '<span class="acro"><b>C</b>hemotype-<b>H</b>it <b>A</b>nalysis of <b>S</b>creening <b>E</b>xperiments</span>'


def versions() -> list[tuple[str, str, str]]:
    """Return (component, version, purpose) for CHASE and every library it uses, read from the installed packages."""
    # The versions are read from the installed packages when CHASE starts, so the list always matches
    # what is really installed.
    def v(p):
        """Installed version of package ``p``, or 'not installed'."""
        try:
            return md.version(p)
        except md.PackageNotFoundError:
            return "not installed"
    return [
        ("CHASE", f"{__version__} ({__date__})", "Chemotype-Hit Analysis of Screening Experiments"),
        ("Python", sys.version.split()[0], platform.platform(terse=True)),
        ("RDKit", v("rdkit"), "SMILES parsing, Morgan fingerprints (r = 3, 2048 bits), Murcko scaffolds, depiction"),
        ("Shiny for Python", v("shiny"), "Web interface"),
        ("Plotly", v("plotly"), "Plots (plotly.js served locally from the plotly package)"),
        ("pandas", v("pandas"), "Tables"),
        ("NumPy", v("numpy"), "Numerics"),
        ("SciPy", v("scipy"), "Hypergeometric test, correlations"),
        ("openpyxl", v("openpyxl"), "Reading .xlsx files"),
        ("PAINS filters", "RDKit FilterCatalog PAINS A/B/C (480 patterns)", "Baell & Holloway, J. Med. Chem. 2010, 53, 2719"),
        ("Covalent alerts", f"{len(_cov()[0])} SMARTS", _cov()[2]),
    ]


def _cov():
    """Return the loaded covalent alerts (imported lazily to keep start-up fast)."""
    from .core.alerts import covalent_alerts
    return covalent_alerts()


def versions_text() -> str:
    """Return the software versions as plain text, one per line, for copying into a methods section."""
    return "\n".join(f"{a}: {b}" for a, b, _ in versions())


# The text of the help page, as HTML. It covers the input files, the steps, the plot, chemotypes,
# triage, cliffs, plates, PAINS and covalent flags, privacy, exports and the methods.
# The software versions, the acknowledgments, the references and the credits are added below it in app.py.
BODY = f"""
<h2>How to use CHASE</h2>
<section><p><b>CHASE</b> ({ACRO}) turns raw single-concentration screening reads into a <b>% response</b> for each compound. It plots them by compound number,
sorts compounds into activity classes and groups them into chemotypes, so biochemists and chemists can see which structural series are active and decide together
what to follow up. This is the Python edition: it runs locally from a conda environment, needs no internet connection, and uses RDKit for all chemistry.</p></section>

<section><h3>Input files</h3>
<h4>Screening data (.xlsx, .xls, .csv, .tsv)</h4>
<p>One row per well with at least these four columns. Header names need not match exactly; reassign them under <b>Columns</b> in step 1.</p>
<table><thead><tr><th>Column</th><th>Example</th><th>Notes</th></tr></thead><tbody>
<tr><td>Raw Data</td><td><code>1031</code></td><td>Numeric readout (signal, RFU, counts).</td></tr>
<tr><td>Plate Number</td><td><code>1</code></td><td>Used for per-plate normalization and QC.</td></tr>
<tr><td>384 well ID</td><td><code>A08</code></td><td>Plate maps read IDs such as A01, A1 or AF48 (96, 384 and 1536 wells).</td></tr>
<tr><td>Batch ID</td><td><code>310480c04</code></td><td>Compound batch, or the control name (e.g. <code>DMSO</code>, <code>CTRL</code>).</td></tr>
</tbody></table>
<ul><li>Rows without a numeric readout or a plate are skipped, so summary rows such as <code>Average</code> at the bottom of a sheet are ignored.</li>
<li>Wells with a signal but no Batch ID are treated as empty and left out.</li>
<li>Format the Batch ID column as text in Excel; otherwise IDs such as <code>274E06</code> become numbers and no longer match the structures file. <b>Data checks</b> flags IDs that look affected.</li></ul>
<h4>Structures (.csv, .tsv, .smi, .xlsx)</h4>
<p>One row per compound: an ID column matching the Batch IDs and a SMILES column, optionally a display name. CHASE picks the ID column that matches the most Batch IDs; matching ignores case and also tries IDs without leading zeros.</p></section>

<section><h3>Step by step</h3><ol>
<li><b>Load files</b> in step 1 (or <b>Load example data</b>). Check the column assignments.</li>
<li><b>Controls &amp; normalization:</b> choose the positive control (100%) and negative control (0%) Batch IDs, per-plate or all-plate normalization and mean or median, then press <b>Normalize &amp; plot</b>.</li>
<li><b>Activity classes:</b> k-means (k = 2-4), fixed thresholds or a robust z cut-off. Changes apply immediately.</li>
<li><b>Chemotypes:</b> Butina clustering at a Tanimoto cut-off, or Bemis-Murcko scaffolds.</li></ol></section>

<section><h3>The plot</h3><ul>
<li>Each point is a sample well: x is its compound number (plate order, then file order), y its % response, colour its compound's class. Hover for the structure; click to pin the compound.</li>
<li>Drag to zoom, double-click to reset. Shaded bands mark alternate plates; dashed lines are class borders; open diamonds and squares are negative and positive control wells.</li>
<li>The narrow panel is a histogram of compound % response on a log count scale.</li>
<li><b>Z' per plate:</b> the strip above the plot shows one bar per plate, spanning that plate's compounds: green for Z' >= 0.5, amber for 0-0.5, red below 0 (values below -1 are drawn at -1). The dotted line marks 0.5; hover a bar for the exact Z' and S/B. Untick <b>Z' per plate</b> to hide it.</li>
<li><b>Colors</b> opens colour pickers for each activity class and both control types. The colours apply to the plots, chips and chemotype strips, are remembered in your browser, and are saved in session files. <b>Reset colours</b> restores the defaults.</li>
<li><b>Save JPEG</b> exports the current view (with your colours) at 3x resolution. <b>Save data (.csv)</b> writes every point of the plot: compound number, plot position, plate, well, well type, Batch ID, name, raw signal, % response, the compound's mean, class, chemotype, plate Z' and SMILES.</li></ul></section>

<section><h3>Chemotypes</h3><ul>
<li><b>Active</b> lists chemotypes with at least one non-inactive compound; <b>Inactive</b> those where every member is inactive.</li>
<li>Cards show the representative structure, hits, hit rate, maximum % response, a strip of member responses (-20 to 120%), and for active chemotypes the enrichment factor and p-value. A green <b>enriched</b> chip marks q &lt; 0.05.</li>
<li>Click a card to list its members and ring them on the plot.</li></ul></section>

<section><h3>Triage board</h3>
<p>Pin a compound and choose <b>Confirm</b>, <b>Resynthesize</b>, <b>Buy analogues</b> or <b>Discard</b>, then add a note (saved when you leave the text box). Click the active decision again to clear it. In a chemotype's panel, one click applies a decision to its undecided members, or to all of them when every member already has one.</p>
<ul><li><b>Nothing is stored in a database.</b> Decisions live only in the open browser tab; CHASE writes nothing to disk. Keep them with <b>Save session</b> (everything) or <b>Export triage (.csv)</b> (the board only), and bring them back with <b>Open session</b> or <b>Import triage (.csv)</b>. The status line above the board shows how many changes are unsaved, and the browser warns you before closing a tab with unsaved decisions.</li>
<li>Each decision records <b>your name</b> (set it above the board; the browser remembers it) and the time. The <b>campaign</b> name (by default the screening file name) labels the exported files.</li>
<li>To work on a board together, pass the session file or the triage CSV: importing adds the decisions to the board, and for a compound that is already decided the newer decision wins. <b>Clear board</b> empties the board after asking.</li>
<li>The triage CSV has the columns campaign, batch_id, compound_name, decision, note, pct_response, activity_class, chemotype, decided_by and decided_at; only batch_id and decision are required when you import a file prepared elsewhere.</ul></section>

<section><h3>PAINS and covalent flags</h3>
<p>Every compound with a structure is checked against two sets of substructure alerts:</p>
<ul><li><b>PAINS</b>: the 480 pan-assay interference patterns of Baell &amp; Holloway,<sup>1</sup> as implemented in RDKit's FilterCatalog (families A, B and C). These substructures (rhodanines, catechols, quinones, Mannich bases ...) often give false positives through aggregation, redox cycling, metal chelation or reactivity.</li>
<li><b>Covalent warheads</b>: electrophiles that can react with protein nucleophiles: acrylamides and other Michael acceptors, haloacetamides, epoxides, aziridines, sulfonyl fluorides, aldehydes, boronic acids, beta-lactams, activated esters, SNAr-prone arenes and more (31 SMARTS). The list is a plain text file you can edit: copy <code>chase/core/covalent_alerts.tsv</code> to <code>~/.chase/covalent_alerts.tsv</code> (or set <code>CHASE_COVALENT</code> to a file) and CHASE uses your version.</li></ul>
<p>Flags appear as red <b>PAINS</b> and violet <b>Covalent</b> chips in tooltips, the compound panel, chemotype members, cards, cliffs and the triage board. On the plot, PAINS hits are marked with an x, covalent hits with an open triangle; the selector in the plot header shows the marks on hits only (default), on all compounds, or not at all. The compound panel names each alert and highlights the matching atoms in the structure. The <b>Flagged hits</b> tile counts hits with any flag. A flag is a prompt to look more closely, not a verdict: some approved drugs carry PAINS substructures, and covalent warheads may be exactly what a campaign is looking for.</p></section>

<section><h3>Activity cliffs</h3>
<p>Pairs of very similar compounds (Tanimoto at or above the cut-off) whose % response differs by at least the chosen number of points, ranked by SALI. Pairs flagged <b>same compound, other batch</b> are two batches of one structure; a large difference there usually means a batch problem rather than SAR.</p></section>

<section><h3>Plates &amp; replicates</h3><ul>
<li><b>Plate maps:</b> every plate in its layout, coloured violet (below 0) - neutral (0) - orange (100+). Click a plate to open it large; hover a well for its structure. The header compares edge and interior wells.</li>
<li><b>Replicates:</b> replicate 1 against replicate 2 (first two wells in plate order) with Pearson r, Spearman rho and the median difference; ringed points differ by more than the threshold.</li>
<li><b>QC table:</b> controls, CVs, S/B and Z' per plate.</li></ul></section>

<section><h3>Sessions</h3>
<p><b>Save session</b> writes one .json file with the screening rows, structures, all settings, the pinned selection and the triage decisions. <b>Open session</b> restores everything and reruns the analysis; the session's decisions join the board (newer decision wins per compound). Session files are interchangeable with the HTML version of CHASE. A session file contains your structures, so share it like the original data.</p></section>

<section><h3>Privacy: what CHASE keeps</h3>
<p>CHASE stores nothing on the computer or server it runs on. Your data exist only in memory while your browser tab is open:</p>
<ul><li><b>Uploaded files</b> are read and their temporary copies deleted at once. On Linux, Shiny's upload folder is placed in RAM (<code>/dev/shm</code>), so uploads never touch the disk.</li>
<li><b>Data, results, structure drawings and the triage board</b> belong to your session alone. They are never shared with other users of the same server, and they are discarded when you close the tab.</li>
<li><b>No database and no log of your data.</b> Requests are not logged.</li>
<li><b>What leaves CHASE</b> is only what you download: session files, CSV exports and JPEGs. Your browser remembers your plot colours and your name in its own local storage, on your computer.</li></ul></section>

<section><h3>Methods</h3>

<h4>Normalization (% response)</h4>
<div class="formula">% response = 100 * (raw - neg) / (pos - neg)</div>
<p><b>raw</b> is the signal of a well, <b>pos</b> and <b>neg</b> are the signal of the positive and negative controls (their mean or median, as chosen under <i>Control statistic</i>). With <i>Per plate</i> the controls of the same plate are used, with <i>All plates</i> the controls of the whole screen. A well that reads like the negative control gets 0%, one that reads like the positive control gets 100%. This works whether the positive control gives a higher or a lower signal than the negative one. A plate that is missing one of its controls is normalized with the all-plate values and is marked in the QC table.</p>

<h4>Plate quality (Z', S/B, CV)</h4>
<div class="formula">Z'  = 1 - 3 * (SDpos + SDneg) / |mean pos - mean neg|
S/B = mean pos / mean neg
CV  = 100 * SD / mean</div>
<p>All three are calculated from each plate's own control wells (SD is the sample standard deviation).</p>
<ul><li><b>Z'</b> tells how well the assay separates the two controls: it compares the gap between their means with their spread. 1 would be a perfect assay with no noise. 0.5 or more is excellent, between 0 and 0.5 the controls are close or noisy, and below 0 they overlap, so the plate can't reliably tell actives from inactives.</li>
<li><b>S/B</b> (signal to background) is how many times stronger the positive control's signal is than the negative control's.</li>
<li><b>CV</b> (coefficient of variation) is the spread of a control's wells as a percentage of their mean. Low values mean the control wells agree with each other.</li></ul>

<h4>Activity classes</h4>
<div class="formula">k-means:  borders halfway between the k cluster centres
Robust z: border = median + n * 1.4826 * MAD</div>
<ul><li><b>k-means</b> groups the compound means (the average % response of each compound's wells) into k clusters along the % response axis. CHASE runs it from several starting points and keeps the solution with the tightest clusters. A class border sits halfway between two neighbouring cluster centres.</li>
<li><b>Robust z</b> uses one border. The median and the MAD (median absolute deviation: the median distance of the values from their median) describe the bulk of inactive compounds without being pulled up by the few actives. Multiplied by 1.4826, the MAD equals the standard deviation for normally distributed data, so the border is "n standard deviations above the median".</li></ul>

<h4>Chemotypes</h4>
<div class="formula">Tanimoto = bits set in both fingerprints / bits set in either fingerprint</div>
<ul><li><b>Fingerprints:</b> RDKit Morgan fingerprints, radius 3 (similar to ECFP6), 2048 bits. Each bit stands for a small atom environment in the molecule.</li>
<li><b>Similarity:</b> Tanimoto similarity runs from 0 (nothing in common) to 1 (identical fingerprints). Taylor-Butina clustering, as described by Taylor<sup>2</sup> and Butina,<sup>3</sup> compares all pairs. The compound with the most neighbours at or above the cut-off becomes the centre of the first chemotype and takes those neighbours with it, then the next compound with the most remaining neighbours, and so on. This is done for up to 20,000 compounds; larger sets are grouped by scaffold.</li>
<li><b>Scaffold:</b> the molecular frameworks of Bemis &amp; Murcko,<sup>4</sup> computed with RDKit: ring systems and the linkers between them, without side chains and without stereochemistry. When chemotypes are made by scaffold, all compounds with the same scaffold form one chemotype.</li></ul>

<h4>Enrichment of hits in a chemotype</h4>
<div class="formula">EF = (hits in chemotype / n) / (hits in library / N)
p  = P(X >= hits in chemotype),  X ~ Hypergeometric(N, hits in library, n)
q  = Benjamini-Hochberg adjusted p</div>
<p>"Hits" are compounds in any class above Inactive. <b>n</b> is the number of compounds in the chemotype and <b>N</b> the number of compounds with a readable structure.</p>
<ul><li><b>EF</b> (enrichment factor) compares the chemotype's hit rate with the hit rate of the whole library. EF = 3 means hits are three times more common in this chemotype than on average.</li>
<li><b>p</b> is the chance of getting at least this many hits if n compounds were picked at random from the library. A small p means the hits are unlikely to cluster in this chemotype by chance.</li>
<li><b>q</b> corrects p for testing many chemotypes at once (Benjamini-Hochberg false discovery rate). Chemotypes with q &lt; 0.05 are marked <i>enriched</i>.</li></ul>

<h4>Activity cliffs</h4>
<div class="formula">SALI(a,b) = |Response_a - Response_b| / (1 - Tanimoto(a,b))</div>
<p><b>Response_a</b> and <b>Response_b</b> are the % responses of compounds a and b, and <b>Tanimoto(a,b)</b> is their fingerprint similarity. SALI, the structure-activity landscape index of Guha &amp; Van Drie,<sup>5</sup> is large when two very similar compounds differ a lot in activity, which points to a small structural change with a big effect. Pairs with identical fingerprints (Tanimoto 1) would divide by zero, so 1 - Tanimoto is never allowed below 0.01. Only pairs at or above both thresholds (Tanimoto and the difference in % response) are listed, sorted by SALI.</p>
<p>Guha and Van Drie defined SALI with potency values such as pIC50. CHASE uses the % response of a single-concentration screen instead, so its SALI values are useful for ranking pairs within one screen but can't be compared with SALI values reported in the literature.</p>

<h4>Edge effect</h4>
<div class="formula">Edge Delta = mean % (outer ring of sample wells) - mean % (interior sample wells)</div>
<p>The outer rows and columns of a plate can behave differently, for example through evaporation or temperature differences. A Delta far from 0 means the edge wells read systematically higher (positive) or lower (negative) than the rest of the plate.</p></section>

<section><h3>Export</h3><p><b>Export results (.csv)</b>: one row per compound with Batch ID, name, SMILES, wells, mean and SD % response, class and chemotype statistics. <b>Export triage (.csv)</b>: the board's decisions, re-importable with <b>Import triage (.csv)</b>. Plots and plate maps save as JPEG; the whole analysis saves as a session file.</p></section>
"""

# Acknowledgments, shown at the bottom of the help page just above the references.
ACK = ('<section><h3>Acknowledgments</h3><p>Thanks to Dr. Francesco Massai, Chemical Biology Consortium Sweden (CBCS-KI), '
       'Karolinska Institutet, for valuable scientific input and feedback during the development of CHASE.</p></section>')
# References in Journal of Medicinal Chemistry style, numbered in the order they are first cited on the page.
REFS = """<section class="refs"><h3>References</h3>
<p>(1) Baell, J. B.; Holloway, G. A. New Substructure Filters for Removal of Pan Assay Interference Compounds (PAINS) from Screening Libraries and for Their Exclusion in Bioassays. <i>J. Med. Chem.</i> <b>2010</b>, <i>53</i>, 2719-2740. DOI: <a href="https://doi.org/10.1021/jm901137j" target="_blank" rel="noopener">10.1021/jm901137j</a></p>
<p>(2) Taylor, R. Simulation Analysis of Experimental Design Strategies for Screening Random Compounds as Potential New Drugs and Agrochemicals. <i>J. Chem. Inf. Comput. Sci.</i> <b>1995</b>, <i>35</i>, 59-67. DOI: <a href="https://doi.org/10.1021/ci00023a009" target="_blank" rel="noopener">10.1021/ci00023a009</a></p>
<p>(3) Butina, D. Unsupervised Data Base Clustering Based on Daylight's Fingerprint and Tanimoto Similarity: A Fast and Automated Way To Cluster Small and Large Data Sets. <i>J. Chem. Inf. Comput. Sci.</i> <b>1999</b>, <i>39</i>, 747-750. DOI: <a href="https://doi.org/10.1021/ci9803381" target="_blank" rel="noopener">10.1021/ci9803381</a></p>
<p>(4) Bemis, G. W.; Murcko, M. A. The Properties of Known Drugs. 1. Molecular Frameworks. <i>J. Med. Chem.</i> <b>1996</b>, <i>39</i>, 2887-2893. DOI: <a href="https://doi.org/10.1021/jm9602928" target="_blank" rel="noopener">10.1021/jm9602928</a></p>
<p>(5) Guha, R.; Van Drie, J. H. Structure-Activity Landscape Index: Identifying and Quantifying Activity Cliffs. <i>J. Chem. Inf. Model.</i> <b>2008</b>, <i>48</i>, 646-658. DOI: <a href="https://doi.org/10.1021/ci7004093" target="_blank" rel="noopener">10.1021/ci7004093</a></p></section>"""
FOOT = '<footer class="foot">Flavio Ballante, Chemical Biology Consortium Sweden (CBCS-KI), Karolinska Institutet - 2026</footer>'
