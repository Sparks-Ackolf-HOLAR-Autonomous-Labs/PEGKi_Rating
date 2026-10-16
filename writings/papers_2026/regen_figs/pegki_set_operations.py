#!/usr/bin/env python3
"""PEGKi figure: set-operation studies on the band-gap archive.

Panels: one column-width panel per fidelity step (three rows). Each row
splits the receiver's shared compounds into the region where one rigid
scissor shift of the donor reproduces the receiver gap within 0.3 eV
(agreement) and the region where it does not (disagreement). Markers: the
transfer score from the donor that defines the split, and from the best
base donor for that region. From
shared/output/set_operation_study/set_operation_study.json.

Size and fonts: RSC column (3.27 in), paper_fig_style fonts.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from paper_fig_style import (  # noqa: E402
    RSC_COL, SOURCE_NAMES, SOURCE_COLOURS, GREY, INK, apply_style, save,
)
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
ART = ROOT / "shared/output/set_operation_study/set_operation_study.json"
OUT = ROOT / "writings/papers_2026/pegki_digital_discovery/figures/set_operations.png"


TITLES = {("mp_gga", "expt_zhuo"): "MP GGA to experiment",
          ("snumat_gga", "snumat_hse"): "SNUMAT GGA to HSE06",
          ("jarvis_optb88vdw", "jarvis_mbj"): "JARVIS OptB88vdW to TB-mBJ"}


def main():
    apply_style()
    d = json.load(open(ART, encoding="utf-8"))
    h1 = d["h1_ownership"]
    h2 = {h["receiver_study"]: h for h in d["h2_rerouting"]}
    fig, axes = plt.subplots(len(h1), 1, figsize=(RSC_COL, 3.9), sharex=True)
    for ax, h in zip(axes, h1):
        a, x = h["donor"], h["receiver"]
        rows = [("agree", h["agree"], h["delta_into_agree"], h["records_agree"]),
                ("disagree", h["disagree"], h["delta_into_disagree"], h["records_disagree"])]
        for yi, (kind, study, dval, nrec) in enumerate(rows):
            best = h2[study]
            bval = dval + best["best_minus_defining"]
            ax.plot([dval, bval], [yi, yi], color=GREY, lw=0.8, zorder=1)
            ax.scatter([dval], [yi], s=26, color=SOURCE_COLOURS[a], zorder=3,
                       edgecolor="white", linewidth=0.5)
            ax.scatter([bval], [yi], s=26, marker="D", zorder=3,
                       color=SOURCE_COLOURS[best["best_donor"]],
                       edgecolor="white", linewidth=0.5)
            ax.annotate(SOURCE_NAMES[best["best_donor"]], (bval, yi),
                        xytext=(4, 0), textcoords="offset points",
                        va="center", fontsize=6.5, color=INK)
        ax.set_yticks([0, 1])
        ax.set_yticklabels([f"agree\n({h['records_agree']:,})",
                            f"disagree\n({h['records_disagree']:,})"], fontsize=6.5)
        ax.set_ylim(-0.6, 1.6)
        ax.invert_yaxis()
        ax.set_title(f"{TITLES[(a, x)]}, scissor {h['scissor_shift_eV']:+.2f} eV",
                     fontsize=7.5, loc="left")
        ax.grid(axis="y", visible=False)
    axes[-1].set_xlabel(r"Transfer score $\Delta$ into the region")
    axes[-1].set_xlim(0.45, 1.65)
    handles = [plt.Line2D([], [], marker="o", ls="", color=GREY, label="donor that defines the split"),
               plt.Line2D([], [], marker="D", ls="", color=GREY, label="best of the eight donors")]
    axes[0].legend(handles=handles, loc="lower right", fontsize=6.5,
                   bbox_to_anchor=(1.0, 1.18), ncol=2, handletextpad=0.2,
                   columnspacing=0.8)
    save(fig, OUT)


if __name__ == "__main__":
    main()
