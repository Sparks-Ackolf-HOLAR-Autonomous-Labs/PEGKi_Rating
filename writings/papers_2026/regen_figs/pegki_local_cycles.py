#!/usr/bin/env python3
"""Regenerate the PEGKi campaign-loop figure (Fig. 10) with uncertainty.

Observation: loops rise, then fade as evaluations accumulate; the bands
show the variation so loose cycles are visible, over long campaigns.
Panels: three, full RSC width, from
shared/output/campaign_loop_bands/campaign_loop_bands.json (800
campaigns per query size, budgets to 9,600 evaluations):
  (a, b) share of campaigns showing a cycle, pairs and 4-sample batches
      in separate panels so the two gradients do not overlap. Solid line: a triangle counts
      as a cycle when the evidence for it is at least even (evidence
      level 0.5). Shaded gradient: the same share when looser (outward
      to 0.05, loops the data cannot rule out) or stricter (inward to
      0.95) evidence is required. Thin dotted line: the yes/no majority
      count.
  (c) loop share of each campaign's comparison record: median line,
      25th to 75th percentile (dark band) and 5th to 95th percentile
      (light band) over campaigns; grey strip, the class-level loop
      share of the full archive (0.006 to 0.043).

Size and fonts: RSC full width, fonts per paper_fig_style, two query
sizes in a blue/vermillion pair, log budget axis.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from paper_fig_style import GREY, INK, RSC_FULL, apply_style, save  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
ART = ROOT / "shared/output/campaign_loop_bands/campaign_loop_bands.json"
DBC = ROOT / "shared/output/decision_body_cycles/decision_body_cycles.json"
OUT = ROOT / "writings/papers_2026/pegki_digital_discovery/figures/local_cycles.png"
SERIES = [("2", "pairs", "#0072B2"), ("4", "4-sample batches", "#D55E00")]
BUDGET = 240


def main():
    apply_style()
    d = json.load(open(ART))
    levels = [float(x) for x in d["config"]["levels"]]
    fig, axs = plt.subplots(1, 3, figsize=(RSC_FULL, 2.85), sharex=True,
                            gridspec_kw={"width_ratios": [1, 1, 1.1]})

    # (a), (b) share with a cycle, graded by evidence level, one panel per size
    for ax, (k, lab, col), tag in zip(axs[:2], SERIES, "ab"):
        rows = d["by_k"][k]
        bs = np.array([int(b) for b in rows])
        cur = {e: np.array([rows[str(b)]["share_cycle_at_level"][str(e)] for b in bs])
               for e in levels}
        # nested fills: darkest next to the 0.5 line, fading outward
        lo_side = [e for e in levels if e >= 0.5]          # stricter, lower curves
        hi_side = [e for e in levels if e <= 0.5][::-1]    # looser, upper curves
        for side in (lo_side, hi_side):
            for i in range(len(side) - 1):
                a = 0.40 * (1 - i / (len(side) - 1)) + 0.06
                ax.fill_between(bs, cur[side[i]], cur[side[i + 1]], color=col,
                                alpha=a, lw=0)
        ax.plot(bs, cur[0.5], color=col, lw=1.6)
        pt = np.array([rows[str(b)]["share_point_cycle"]["mean"] for b in bs])
        ax.plot(bs, pt, color=INK, lw=0.8, ls=":")
        ax.axvline(BUDGET, color=INK, lw=0.7, ls="--")
        ax.set_ylim(0, 1.02)
        ax.set_title(f"({tag}) Cycles, {lab}", loc="left")
        ax.set_xlabel("Evaluations spent")
    axs[0].set_ylabel("Share of campaigns with a cycle")
    axs[1].tick_params(labelleft=False)
    axs[0].text(BUDGET * 1.12, 0.03, "campaign\nbudget", fontsize=7, color=INK)

    # (c) loop share over campaigns
    bx = axs[2]
    try:
        mc = json.load(open(DBC))["material_classes"]
        vals = [v["hodge"]["loop"] for s, v in mc.items() if not s.startswith("_")]
        bx.axhspan(min(vals), max(vals), color=GREY, alpha=0.45, lw=0)
    except (OSError, KeyError):
        pass
    for k, lab, col in SERIES:
        rows = d["by_k"][k]
        bs = np.array([int(b) for b in rows])
        q = {p: np.array([rows[str(b)]["loop_share"][p] for b in bs])
             for p in ("p05", "p25", "p50", "p75", "p95")}
        bx.fill_between(bs, q["p05"], q["p95"], color=col, alpha=0.12, lw=0)
        bx.fill_between(bs, q["p25"], q["p75"], color=col, alpha=0.30, lw=0)
        bx.plot(bs, q["p50"], color=col, lw=1.6)
    bx.axvline(BUDGET, color=INK, lw=0.7, ls="--")
    bx.set_ylim(0, 0.8)
    bx.set_ylabel("Loop share of the win table")
    bx.set_xlabel("Evaluations spent")
    bx.set_title("(c) How loose the loops are", loc="left")
    for a in axs:
        a.set_xscale("log")
        a.set_xlim(8, 9600)

    h = [Line2D([], [], color=c, lw=1.6, label=l) for _, l, c in SERIES]
    h += [Patch(facecolor=GREY, alpha=0.5, lw=0,
                label="(a, b) looser (up) or stricter (down) evidence"),
          Line2D([], [], color=INK, lw=0.8, ls=":", label="(a, b) yes/no majority count"),
          Patch(facecolor=GREY, alpha=0.35, lw=0,
                label="(c) middle half / 90% of campaigns"),
          Patch(facecolor=GREY, alpha=0.8, lw=0,
                label="(c) class-level loop share, full archive")]
    fig.legend(handles=h, loc="outside lower center", ncol=3, handlelength=1.6,
               columnspacing=1.2)
    save(fig, OUT)


if __name__ == "__main__":
    main()
