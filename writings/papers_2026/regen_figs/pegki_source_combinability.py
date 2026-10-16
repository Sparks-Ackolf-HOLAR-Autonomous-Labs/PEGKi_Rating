#!/usr/bin/env python3
"""PEGKi figure: multi-source set operations decide which data sources can be
combined.

Data: shared/output/source_combinability/source_combinability.json.
(a) Lower-triangle matrix over the eight sources: Kendall tau of the gaps on
    the compounds both sources report (A cap B), with the shared count below.
(b) Decision map, one point per pair: x = share of shared compounds whose gaps
    agree within 0.3 eV after one rigid scissor shift (can the values be
    merged?), y = Kendall tau on A cap B (can the ordering be combined?);
    90% bootstrap intervals as bars, colour by relation of the pair. The MP GGA /
    ICSD-matched pair (same calculations) is omitted from (b).
(c) Average-linkage tree on 1 - tau with the grouping cut.
Size and fonts: paper_fig_style, RSC full width at print size.
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from paper_fig_style import (  # noqa: E402
    RSC_FULL, MP_BLUE, SNUMAT_PURPLE, EXPT_GOLD, GREY, INK, apply_style, save)
import matplotlib.pyplot as plt  # noqa: E402
from scipy.cluster.hierarchy import dendrogram, set_link_color_palette  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
ART = ROOT / "shared/output/source_combinability/source_combinability.json"
OUT = ROOT / "writings/papers_2026/pegki_digital_discovery/figures/source_combinability.png"

SHORT = {"expt_zhuo": "Experiment", "snumat_hse": "SNUMAT HSE06",
         "jarvis_mbj": "JARVIS TB-mBJ", "mp_r2scan": "MP r$^2$SCAN",
         "snumat_gga": "SNUMAT GGA", "jarvis_optb88vdw": "JARVIS OptB88vdW",
         "mp_gga": "MP GGA", "mp_gga_expsynth": "MP GGA (ICSD)"}
GGA = {"snumat_gga", "jarvis_optb88vdw", "mp_gga", "mp_gga_expsynth"}


def relation(a, b):
    if "expt_zhuo" in (a, b):
        return "with experiment"
    if a in GGA and b in GGA:
        return "both GGA level"
    return "different DFT levels"


# pair type by marker shape and neutral shades (the database palette is kept
# for sources, not pair types)
COL = {"both GGA level": INK, "different DFT levels": "#777777",
       "with experiment": "#AAAAAA"}
MK = {"both GGA level": "o", "different DFT levels": "s", "with experiment": "^"}
# pairs discussed in the text
LABEL = {frozenset({"snumat_gga", "snumat_hse"}): ("SNUMAT GGA/HSE06", (3, 6), "left"),
         frozenset({"mp_gga", "expt_zhuo"}): ("MP GGA/Exp.", (6, -2), "left"),
         frozenset({"jarvis_optb88vdw", "jarvis_mbj"}): ("OptB88vdW/TB-mBJ", (4, -10), "left"),
         frozenset({"mp_gga", "snumat_gga"}): ("MP/SNUMAT GGA", (4, -9), "left")}


def main():
    apply_style()
    d = json.load(open(ART, encoding="utf-8"))
    order, P = d["order"], d["pairs"]
    n = len(order)

    def get(a, b):
        return P.get(f"{a}|{b}") or P.get(f"{b}|{a}")

    fig = plt.figure(figsize=(RSC_FULL, 3.0))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1.0, 0.85])
    ax0, ax1, ax2 = (fig.add_subplot(gs[0, i]) for i in range(3))

    # (a) tau matrix, lower triangle
    M = np.full((n, n), np.nan)
    for i in range(n):
        for j in range(i):
            r = get(order[i], order[j])
            if "tau" in r:
                M[i, j] = r["tau"]
    im = ax0.imshow(M, cmap="Blues", vmin=0.5, vmax=1.0)
    for i in range(n):
        for j in range(i):
            r = get(order[i], order[j])
            ax0.text(j, i, f"{r['tau']:.2f}", ha="center", va="center", fontsize=6.5,
                     color="white" if r["tau"] > 0.72 else INK)
            # shared-compound counts are given in the caption (0.5k to 24k)
    ax0.set_xticks(range(n))
    ax0.set_xticklabels([SHORT[s] for s in order], rotation=55, ha="right", fontsize=6)
    ax0.set_yticks(range(n))
    ax0.set_yticklabels([SHORT[s] for s in order], fontsize=6)
    ax0.grid(False)
    for sp in ax0.spines.values():
        sp.set_visible(False)
    ax0.set_xlim(-0.5, n - 1.5)
    ax0.set_ylim(n - 0.5, 0.5)
    cb = fig.colorbar(im, ax=ax0, shrink=0.7, pad=0.02)
    cb.set_label("Kendall $\\tau$", fontsize=6.5)
    cb.ax.tick_params(labelsize=6)
    ax0.set_title("(a) ordering on $A\\cap B$", fontsize=8, loc="left")

    # (b) decision map
    for r in P.values():
        # MP GGA (ICSD-matched) is a subset of the same MP GGA calculations
        if "tau" not in r or {r["a"], r["b"]} == {"mp_gga", "mp_gga_expsynth"}:
            continue
        rel = relation(r["a"], r["b"])
        x, y = r["value_agree_frac"], r["tau"]
        xe = np.array([[x - r["value_agree_ci90"][0]], [r["value_agree_ci90"][1] - x]])
        ye = np.array([[y - r["tau_ci90"][0]], [r["tau_ci90"][1] - y]])
        ax1.errorbar(x, y, xerr=xe, yerr=ye, fmt=MK[rel], ms=3.4, color=COL[rel],
                     ecolor=COL[rel], elinewidth=0.6, capsize=0, label=rel)
        lab_ = LABEL.get(frozenset({r["a"], r["b"]}))
        if lab_:
            ax1.annotate(lab_[0], (x, y), xytext=lab_[1], textcoords="offset points", fontsize=6,
                         color=INK, ha=lab_[2])
    h, lab = ax1.get_legend_handles_labels()
    uniq = dict(zip(lab, h))
    ax1.legend(uniq.values(), uniq.keys(), loc="lower right", fontsize=6,
               handletextpad=0.2, borderaxespad=0.2)
    ax1.set_xlim(0.2, 1.0)
    ax1.set_ylim(0.5, 1.0)
    ax1.set_xlabel("fraction of shared compounds\nwithin 0.3 eV after scissor shift", fontsize=7)
    ax1.set_ylabel("Kendall $\\tau$ on shared compounds", fontsize=7)
    ax1.set_title("(b) values versus ordering", fontsize=8, loc="left")

    # (c) tree on 1 - tau
    Z = np.array(d["linkage"])
    cut = d["config"]["group_cut"]
    set_link_color_palette([MP_BLUE, SNUMAT_PURPLE])
    dendrogram(Z, ax=ax2, orientation="right", labels=[SHORT[s] for s in order],
               color_threshold=cut, above_threshold_color=GREY, leaf_font_size=6)
    ax2.axvline(cut, color=INK, lw=0.7, ls=":")
    ax2.set_xlabel("1 $-$ Kendall $\\tau$", fontsize=7)
    ax2.grid(False)
    ax2.spines["left"].set_visible(False)
    ax2.set_title("(c) pooled groups", fontsize=8, loc="left")
    save(fig, OUT)


if __name__ == "__main__":
    main()
