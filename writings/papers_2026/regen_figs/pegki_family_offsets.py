#!/usr/bin/env python3
"""Regenerate the PEGKi family-resolved offsets figure at true print size.

Panels: two, matching the caption of Fig. \\ref{fig:families}.
Left: per-source offsets relative to mp_gga, refitted within one anion
family at a time, with 200-resample bootstrap intervals; grey diamonds
are the global (all-family) fit. Right: the three ladder spacings per
family with bootstrap intervals plus the global fit. All numbers come
from shared/output/family_offsets/family_offsets.json:
  families[fam].offsets_ref / offset_ref_ci  (left panel)
  global.offsets_ref / offset_ref_ci          (left panel, diamonds)
  families[fam].spacings / spacing_ci         (right panel)
  global.spacings                              (right panel, diamonds)

Size and fonts: RSC full-width figure (6.73 in), fonts per paper_fig_style
(7-8 pt at print size). Labels use "data source".
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from paper_fig_style import (  # noqa: E402
    RSC_FULL, SOURCE_NAMES, GREY, INK, apply_style, save,
)
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
ART = ROOT / "shared/output/family_offsets/family_offsets.json"
OUT = ROOT / "writings/papers_2026/pegki_digital_discovery/figures/family_offsets.png"

# CVD-safe family palette (Okabe-Ito)
FAMILY_COLOURS = {
    "oxide": "#0072B2",
    "halide": "#D55E00",
    "chalcogenide": "#009E73",
    "nitride": "#E69F00",
}
FAMILY_MARKERS = {"oxide": "o", "halide": "s", "chalcogenide": "^", "nitride": "D"}
FAMILIES = ["oxide", "halide", "chalcogenide", "nitride"]

SPACING_ORDER = ["hse_minus_gga_snumat", "mbj_minus_optb88_jarvis", "expt_minus_gga_mp"]
SPACING_LABELS = {
    "hse_minus_gga_snumat": "HSE06\nvs GGA\n(SNUMAT)",
    "mbj_minus_optb88_jarvis": "TB-mBJ vs\nOptB88vdW\n(JARVIS)",
    "expt_minus_gga_mp": "Experiment\nvs GGA\n(Zhuo / MP)",
}


def errbar(ax, x, val, ci, colour, marker, ms=3.4):
    ax.errorbar([x], [val], yerr=[[val - ci[0]], [ci[1] - val]],
                fmt=marker, color=colour, ms=ms, mfc=colour, mec="white",
                mew=0.4, elinewidth=1.0, capsize=1.6, capthick=1.0, zorder=3)


def main():
    apply_style()
    data = json.load(open(ART))
    glob = data["global"]
    fams = data["families"]

    # Sources ordered by the global ladder, descending (top of ladder left)
    sources = sorted(glob["offsets_ref"], key=lambda s: -glob["offsets_ref"][s])

    fig, (axL, axR) = plt.subplots(
        1, 2, figsize=(RSC_FULL, 2.95), width_ratios=[2.35, 1.0])

    # Left: per-source offsets relative to mp_gga, by family
    dodge = {"oxide": -0.27, "halide": -0.09, "chalcogenide": 0.09, "nitride": 0.27}
    for i, src in enumerate(sources):
        axL.errorbar([i], [glob["offsets_ref"][src]],
                     yerr=[[glob["offsets_ref"][src] - glob["offset_ref_ci"][src][0]],
                           [glob["offset_ref_ci"][src][1] - glob["offsets_ref"][src]]],
                     fmt="D", color=GREY, ms=4.6, mfc=GREY, mec="white", mew=0.4,
                     elinewidth=1.0, capsize=0, zorder=2)
        for fam in FAMILIES:
            errbar(axL, i + dodge[fam], fams[fam]["offsets_ref"][src],
                   fams[fam]["offset_ref_ci"][src],
                   FAMILY_COLOURS[fam], FAMILY_MARKERS[fam])
    axL.set_xticks(range(len(sources)))
    axL.set_xticklabels([SOURCE_NAMES[s] for s in sources],
                        rotation=32, ha="right", rotation_mode="anchor")
    axL.set_ylabel("offset vs MP GGA (log units)")
    axL.set_title("Per-source offsets refitted within one anion family", fontsize=8.5)
    axL.grid(axis="x", visible=False)
    axL.axhline(0, color=INK, linewidth=0.7)

    # shared legend (family markers + global diamond), inside left panel
    handles = [plt.Line2D([], [], color=GREY, marker="D", ls="none", ms=4.6,
                          mec="white", mew=0.4, label="global fit (all families)")]
    for fam in FAMILIES:
        handles.append(plt.Line2D([], [], color=FAMILY_COLOURS[fam],
                                  marker=FAMILY_MARKERS[fam], ls="none", ms=3.4,
                                  mec="white", mew=0.4, label=fam))
    axL.legend(handles=handles, loc="upper right", ncol=1, handletextpad=0.3,
               borderaxespad=0.2, labelspacing=0.25)

    # Right: ladder spacings per family
    for j, sp in enumerate(SPACING_ORDER):
        axR.errorbar([j], [glob["spacings"][sp]], fmt="D", color=GREY, ms=4.6,
                     mfc=GREY, mec="white", mew=0.4, zorder=2)
        for fam in FAMILIES:
            errbar(axR, j + dodge[fam] * 0.9, fams[fam]["spacings"][sp],
                   fams[fam]["spacing_ci"][sp],
                   FAMILY_COLOURS[fam], FAMILY_MARKERS[fam])
    axR.set_xticks(range(len(SPACING_ORDER)))
    axR.set_xticklabels([SPACING_LABELS[s] for s in SPACING_ORDER], fontsize=6.5)
    axR.set_xlim(-0.6, len(SPACING_ORDER) - 0.4)
    axR.set_ylabel("ladder spacing (log units)")
    axR.set_title("Ladder spacings by family", fontsize=8.5)
    axR.grid(axis="x", visible=False)

    save(fig, OUT)


if __name__ == "__main__":
    main()
