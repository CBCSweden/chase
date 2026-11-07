"""Activity classes from compound % response.

Author: Flavio Ballante, Chemical Biology Consortium Sweden (CBCS-KI), Karolinska Institutet - 2026
Acknowledgments: thanks to Dr. Francesco Massai, Chemical Biology Consortium Sweden (CBCS-KI), Karolinska Institutet,
for valuable scientific input and feedback during the development of CHASE.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

# The activity classes for 2, 3 or 4 classes. Each entry is (label shown to the user, short name used
# for the colours and the CSS).
CLASS_SETS = {
    2: [("Inactive", "inactive"), ("Active", "active")],
    3: [("Inactive", "inactive"), ("Weak", "weak"), ("Active", "active")],
    4: [("Inactive", "inactive"), ("Weak", "weak"), ("Moderate", "moderate"), ("Active", "active")],
}


@dataclass
class Classes:
    """Activity classes of all compounds: labels, class borders, class index per compound, counts and display colours."""
    labels: list[tuple[str, str]]
    bounds: list[float]
    cls: np.ndarray          # activity class of each compound (0 = inactive)
    counts: list[int]
    hint: str
    colors: dict = field(default_factory=dict)   # colours chosen by the user, if any (class names plus 'pos' and 'neg')

    @property
    def top(self) -> str:
        """Name of the most active class (e.g. ``"Active"``)."""
        return self.labels[-1][0]


def kmeans1d(vals: np.ndarray, k: int) -> np.ndarray:
    """Exact-enough 1-D k-means: best of several Lloyd runs (quantile, even-spaced and k-means++ starts)."""
    v = np.sort(np.asarray(vals, dtype=float))
    n = len(v)
    # With fewer values than clusters there is nothing to cluster: every value becomes its own centre.
    if n <= k:
        return v.copy()
    # Running sums of the sorted values. With them, the mean and the spread of any cluster can be worked
    # out from two numbers, instead of adding up all its values again in every round.
    pre = np.concatenate([[0.0], np.cumsum(v)])
    pre2 = np.concatenate([[0.0], np.cumsum(v * v)])

    def cuts_for(c):
        """Indices into the sorted values where each cluster starts, given sorted centres."""
        return [0] + [int(np.searchsorted(v, (c[j - 1] + c[j]) / 2, side="left")) for j in range(1, k)] + [n]

    def lloyd(c0):
        """Run Lloyd iterations from the starting centres ``c0``; return (centres, within-cluster sum of squares)."""
        c = np.sort(np.asarray(c0, dtype=float))
        for _ in range(200):
            cuts = cuts_for(c)
            # The new centre of each cluster is the mean of the values in it. An empty cluster keeps its old centre.
            nc = np.array([(pre[cuts[j + 1]] - pre[cuts[j]]) / (cuts[j + 1] - cuts[j]) if cuts[j + 1] > cuts[j] else c[j] for j in range(k)])
            nc.sort()
            if np.allclose(nc, c, atol=1e-9):
                c = nc
                break
            c = nc
        cuts = cuts_for(c)
        # How tight the clusters are: the sum of squared distances to the cluster centres, written as
        # sum(x^2) - sum(x)^2/n for each cluster and taken from the running sums.
        sse = 0.0
        for j in range(k):
            a, b = cuts[j], cuts[j + 1]
            if b > a:
                s = pre[b] - pre[a]
                sse += pre2[b] - pre2[a] - s * s / (b - a)
        return c, sse

    # We try several starting points and keep the best result. Screening data are mostly inactive
    # compounds, and a single start can easily miss a small group of very active ones.
    # The starts are: evenly spaced quantiles, evenly spaced values, and 12 random k-means++ starts
    # (with a fixed seed, so the result is the same every time).
    inits = [v[np.floor(((np.arange(k) + 0.5) / k) * (n - 1)).astype(int)], v[0] + (v[-1] - v[0]) * (np.arange(k) + 0.5) / k]
    rng = np.random.default_rng(7)
    for _ in range(12):
        c = [v[rng.integers(n)]]
        while len(c) < k:
            # k-means++ start: pick the next centre at random, favouring values far from the centres already chosen.
            d = np.min((v[:, None] - np.array(c)[None, :]) ** 2, axis=1)
            tot = d.sum()
            c.append(v[rng.choice(n, p=d / tot)] if tot > 0 else v[rng.integers(n)])
        inits.append(np.array(c))
    best = min((lloyd(c) for c in inits), key=lambda r: r[1])
    return best[0]


def assign(means: np.ndarray, method: str = "kmeans", k: int = 3, thr_a: float | None = 50, thr_w: float | None = 20,
           rob_n: float = 3.0) -> Classes:
    """Split compounds into activity classes from their mean % response.

    Args:
        means: Mean % response per compound.
        method: ``"kmeans"``, ``"threshold"`` or ``"robust"``.
        k: Number of classes for k-means (2-4).
        thr_a: Active threshold (% response) for the threshold method.
        thr_w: Weak threshold for the threshold method; None gives two classes.
        rob_n: Multiple of the robust SD (1.4826 * MAD) for the robust method.

    Returns:
        A :class:`Classes` object; class 0 is always Inactive.
    """
    means = np.asarray(means, dtype=float)
    if not len(means):
        return Classes(CLASS_SETS[2], [], np.zeros(0, int), [0, 0], "")
    if method == "kmeans":
        # Never ask for more classes than there are different values.
        k = min(k, max(2, len(np.unique(np.round(means, 3)))))
        labels = CLASS_SETS[k]
        centers = kmeans1d(means, k)
        bounds = [(centers[j - 1] + centers[j]) / 2 for j in range(1, k)]
        hint = "Cluster centres: " + " | ".join(f"{c:.1f}%" for c in centers) + ". Class borders sit halfway between centres."
    elif method == "threshold":
        if thr_w is not None and thr_a is not None and thr_w < thr_a:
            labels, bounds = CLASS_SETS[3], [thr_w, thr_a]
        else:
            labels, bounds = CLASS_SETS[2], [thr_a if thr_a is not None else 50.0]
        hint = "Leave 'Weak' empty for a two-class split."
    else:
        # Robust z: the cut-off is the median plus n robust standard deviations. The robust SD is
        # 1.4826 * MAD (median absolute deviation), which equals the normal SD for normally distributed data
        # but is hardly affected by the few strong actives.
        med = float(np.median(means))
        mad = float(np.median(np.abs(means - med))) * 1.4826
        labels, bounds = CLASS_SETS[2], [med + rob_n * mad]
        hint = f"Median {med:.2f}%, robust SD {mad:.2f}% -> cut-off {bounds[0]:.1f}%."
    # A compound's class is the number of class borders at or below its mean (so 0 = Inactive).
    cls = np.zeros(len(means), dtype=int)
    for j, b in enumerate(bounds):
        cls[means >= b] = j + 1
    counts = [int((cls == i).sum()) for i in range(len(labels))]
    return Classes(labels, [float(b) for b in bounds], cls, counts, hint)
