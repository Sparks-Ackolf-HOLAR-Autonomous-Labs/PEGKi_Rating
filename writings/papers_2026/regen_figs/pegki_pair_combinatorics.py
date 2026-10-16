#!/usr/bin/env python3
"""PEGKi figure: the pair structure behind the hypergraph and k-sweep results.

Panels: full width, three panels sharing one idea (class pairs of a
dominance table). (a) observed per-round residual of the whole-round fits
against the value predicted by counting near-tie pairs, one marker per data
source, identity line. (b) mean final Kendall tau of the 50-seed k-sweep
(budget 240) against k. (c) evaluations needed to see every pair of 8
classes at least once (Schonheim covering bound) against k. Panels (b) and
(c) share the k axis so the two saturations line up. From
shared/output/pair_combinatorics/pair_combinatorics.json.

Size and fonts: RSC full width (6.73 in), paper_fig_style fonts, source
colours.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from paper_fig_style import (  # noqa: E402
    RSC_FULL, SOURCE_COLOURS, SOURCE_NAMES, GREY, INK, MP_BLUE, apply_style, save,
)
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
ART = ROOT / "shared/output/pair_combinatorics/pair_combinatorics.json"
OUT = ROOT / "writings/papers_2026/pegki_digital_discovery/figures/pair_combinatorics.png"

SHORT = {"expt_zhuo": "Expt", "jarvis_mbj": "TB-mBJ", "jarvis_optb88vdw": "OptB88vdW",
         "mp_gga": "MP GGA", "mp_gga_expsynth": "ICSD-m.", "mp_r2scan": "r$^2$SCAN",
         "snumat_gga": "SN GGA", "snumat_hse": "HSE06"}
OFFS = {"expt_zhuo": (4, -2), "jarvis_mbj": (-34, 3), "jarvis_optb88vdw": (-30, 4),
        "mp_gga": (-32, -3), "mp_gga_expsynth": (4, -6), "mp_r2scan": (-40, -2),
        "snumat_gga": (4, -4), "snumat_hse": (5, -7)}


def main():
    apply_style()
    d = json.load(open(ART, encoding="utf-8"))
    fig, axes = plt.subplots(1, 3, figsize=(RSC_FULL, 2.35))

    ax = axes[0]
    ident = d["round_residual_identity"]
    lo, hi = 0.20, 0.40
    ax.plot([lo, hi], [lo, hi], "--", color=GREY, lw=0.9)
    for r in ident:
        s = r["source"]
        ax.plot(r["predicted_vs_majority"], r["observed_edge_residual"], "o",
                color=SOURCE_COLOURS[s], ms=4.5)
        ax.annotate(SHORT[s], (r["predicted_vs_majority"], r["observed_edge_residual"]),
                    xytext=OFFS[s], textcoords="offset points", fontsize=6.5, color=INK)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel("Predicted from pair table")
    ax.set_ylabel("Observed per-round residual")
    rs = d["round_residual_identity_summary"]
    ax.set_title(f"(a) One round, Pearson {rs['pearson']:+.2f}", loc="left")

    corr = d["k_sweep_correlations"]
    ks = sorted(int(k[1:]) for k in corr)
    ax = axes[1]
    ax.plot(ks, [corr[f"k{k}"]["mean_final_tau"] for k in ks], "-o", color=MP_BLUE, ms=3.5)
    ax.set_xticks(ks)
    ax.set_xlabel("Classes per query event, k")
    ax.set_ylabel("Final Kendall $\\tau$ (budget 240)")
    ax.set_title("(b) Measured k-sweep", loc="left")

    ax = axes[2]
    cov = d["covering"]["8"]
    kk = [r["k"] for r in cov if r["k"] <= 7]
    ev = [r["evaluations"] for r in cov if r["k"] <= 7]
    ax.bar(kk, ev, color=[MP_BLUE if k in ks else GREY for k in kk], width=0.65)
    for k, e in zip(kk, ev):
        ax.text(k, e + 1, str(e), ha="center", fontsize=6.5, color=INK)
    ax.set_xticks(kk)
    ax.set_xlabel("Classes per query event, k")
    ax.set_ylabel("Evaluations to see all 28 pairs")
    ax.set_title("(c) Pair covering, n = 8", loc="left")
    save(fig, OUT)


if __name__ == "__main__":
    main()
