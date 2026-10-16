#!/usr/bin/env python3
"""PEGKi figure: the fidelity ladder, rung position and rung spacing.

Panels: two, sharing one row order (sources sorted by offset).
(a) pooled Plackett-Luce offset delta per data source, from
    shared/output/hypergraph_trial_bandgap/hypergraph_trial_bandgap.json
    (pooled.channel_offsets_sorted).
(b) learned margin scale a_c from the non-additive pooled model, with
    100-resample bootstrap intervals, from
    shared/output/nonadditive_pooled/nonadditive_pooled.json
    (models.scale.margin_scales, scale_bootstrap_ci).
Takeaway: the rank-only fit puts the sources in the known order (a), and
the learned model shows how widely each source spreads the classes (b).

Size and fonts: RSC full width (6.73 in), paper_fig_style fonts.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from paper_fig_style import (  # noqa: E402
    RSC_FULL, SOURCE_COLOURS, SOURCE_NAMES, INK, GREY, apply_style, save,
)
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
ART_O = ROOT / "shared/output/hypergraph_trial_bandgap/hypergraph_trial_bandgap.json"
ART_S = ROOT / "shared/output/nonadditive_pooled/nonadditive_pooled.json"
OUT = ROOT / "writings/papers_2026/pegki_digital_discovery/figures/ladder.png"


def main():
    apply_style()
    offsets = json.load(open(ART_O, encoding="utf-8"))["pooled"]["channel_offsets_sorted"]
    s = json.load(open(ART_S, encoding="utf-8"))
    scales = s["models"]["scale"]["margin_scales"]
    ci = s["scale_bootstrap_ci"]
    order = sorted(offsets, key=lambda c: offsets[c])
    ys = list(range(len(order)))

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(RSC_FULL, 2.6), sharey=True,
                                   width_ratios=[1.15, 1.0])
    vals = [offsets[c] for c in order]
    axL.barh(ys, vals, color=[SOURCE_COLOURS[c] for c in order], height=0.62)
    axL.axvline(0, color=INK, lw=0.8)
    for y, v in zip(ys, vals):
        axL.annotate(f"{v:+.2f}", xy=(v, y), xytext=(3 if v >= 0 else -3, 0),
                     textcoords="offset points", va="center",
                     ha="left" if v >= 0 else "right", fontsize=7)
    axL.set_xlim(min(vals) - 0.8, max(vals) + 0.8)
    axL.set_yticks(ys)
    axL.set_yticklabels([SOURCE_NAMES[c] for c in order])
    axL.set_xlabel("offset $\\delta$ (log units; higher = wider gaps)")
    axL.set_title("(a) Rung position", loc="left")
    axL.grid(axis="y", visible=False)

    for y, c in zip(ys, order):
        axR.barh(y, scales[c], color=SOURCE_COLOURS[c], height=0.62)
        axR.plot(ci[c], [y, y], color=INK, lw=1.0)
    axR.axvline(1.0, color=GREY, ls="--", lw=1.0)
    axR.text(1.02, len(order) - 0.45, "additive model", color=GREY, fontsize=7,
             va="center")
    axR.set_xlim(0, 1.95)
    axR.set_xlabel("margin scale $a$ (spread of the classes)")
    axR.set_title("(b) Rung spacing", loc="left")
    axR.grid(axis="y", visible=False)
    save(fig, OUT)


if __name__ == "__main__":
    main()
