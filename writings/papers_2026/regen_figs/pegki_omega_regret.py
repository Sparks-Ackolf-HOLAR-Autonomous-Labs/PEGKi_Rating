#!/usr/bin/env python3
"""Regenerate the PEGKi omega_regret figure at true print size.

Panels: one, the 8 x 8 heatmap of omega_regret(A->B) for all 56
donor-receiver routes (pair arm):
  - only cells with omega_regret >= 0.03 carry printed values, in
    plain Kendall-tau units (no x100 encoding), so the routes that
    matter stand out from the 56 cells;
  - the Experiment receiver column, where the costliest routes point,
    is outlined so the laboratory reading is immediate.
The Spearman +0.52 concordance with the interim dispersion omega is
reported in the body text (Section "omega").
All numbers come from shared/output/omega_regret/omega_regret.json
(omega_regret.pair[donor][receiver].omega_regret).

Size and fonts: single-column RSC figure (3.27 in), fonts per
paper_fig_style, short source names so labels stay horizontal-readable.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from paper_fig_style import INK, apply_style, save  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.patches as mpatches  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
ART = ROOT / "shared/output/omega_regret/omega_regret.json"
OUT = ROOT / "writings/papers_2026/pegki_digital_discovery/figures/omega_regret.png"

# short reader-facing names so the axis labels stay legible in-column
SHORT_NAMES = {
    "expt_zhuo": "Experiment",
    "jarvis_mbj": "TB-mBJ",
    "jarvis_optb88vdw": "OptB88vdW",
    "mp_gga": "MP GGA",
    "mp_gga_expsynth": "MP GGA (ICSD)",
    "mp_r2scan": "r$^2$SCAN",
    "snumat_gga": "SNUMAT GGA",
    "snumat_hse": "SNUMAT HSE",
}
# heatmap order: experiment first, then grouped by database of origin
HEAT_ORDER = [
    "expt_zhuo", "jarvis_mbj", "jarvis_optb88vdw", "mp_gga",
    "mp_gga_expsynth", "mp_r2scan", "snumat_gga", "snumat_hse",
]
ANNOTATE_MIN = 0.03  # print values only where the route costs at least 0.03 in tau


def main():
    apply_style()
    data = json.load(open(ART))
    pair = data["omega_regret"]["pair"]

    n = len(HEAT_ORDER)
    M = np.full((n, n), np.nan)
    for i, don in enumerate(HEAT_ORDER):
        for j, rec in enumerate(HEAT_ORDER):
            if rec in pair[don]:
                M[i, j] = pair[don][rec]["omega_regret"]
    vmax = np.nanmax(M)

    fig, ax = plt.subplots(figsize=(3.27, 3.05))
    Mm = np.ma.masked_invalid(M)
    cmap = plt.get_cmap("Reds").copy()
    cmap.set_bad("#E8E8E8")
    im = ax.imshow(Mm, cmap=cmap, vmin=0.0, vmax=vmax, aspect="equal")
    # thin white grid between cells
    ax.set_xticks(np.arange(-0.5, n, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.7)
    ax.tick_params(which="minor", length=0)

    ax.set_xticks(range(n))
    ax.set_xticklabels([SHORT_NAMES[s] for s in HEAT_ORDER],
                       rotation=45, ha="right", rotation_mode="anchor",
                       fontsize=6.5)
    ax.set_yticks(range(n))
    ax.set_yticklabels([SHORT_NAMES[s] for s in HEAT_ORDER], fontsize=6.5)
    ax.set_xlabel("receiver data source $B$", fontsize=7.5)
    ax.set_ylabel("donor data source $A$", fontsize=7.5)
    ax.grid(which="major", visible=False)

    for i in range(n):
        for j in range(n):
            v = M[i, j]
            if not np.isnan(v) and v >= ANNOTATE_MIN:
                # ".08" not "0.08": the leading zero costs a third of the
                # cell width and carries no information
                ax.text(j, i, f"{v:.2f}".lstrip("0"), ha="center",
                        va="center", fontsize=6,
                        color="white" if v > 0.55 * vmax else INK)
    # outline the Experiment receiver column: the costly routes point here
    jx = HEAT_ORDER.index("expt_zhuo")
    ax.add_patch(mpatches.Rectangle(
        (jx - 0.5, -0.5), 1.0, n, fill=False, edgecolor=INK,
        linewidth=1.1, zorder=4))

    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("early harm of blind reuse,\n"
                 "$\\omega_{\\mathrm{regret}}$ (Kendall $\\tau$ units)",
                 fontsize=7)
    cb.ax.tick_params(labelsize=6.5)
    save(fig, OUT)


if __name__ == "__main__":
    main()
