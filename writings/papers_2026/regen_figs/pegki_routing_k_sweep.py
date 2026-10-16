#!/usr/bin/env python3
"""PEGKi figure: routing rules across query arity k.

Panels: one column-width panel. x = k (2 = pairwise arm, 3-6 = k-way arms
of matched-seed paired campaigns); y = realised early transfer gain on
held-out receivers (leave one receiver out). Lines: maximin linear program
(identical to uniform within 0.001 at every k), routing learned from the
linear-program cold start, best measured donor. Bands: 90% percentile
intervals over 200 resamples of the 50 matched seeds.
From shared/output/routing_k_sweep_50seed/routing_k_sweep.json (means) and
routing_k_sweep_bootstrap.json (bands; shared/scripts/routing_k_sweep_bootstrap.py).

Size and fonts: RSC column (3.27 in), paper_fig_style fonts.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from paper_fig_style import (  # noqa: E402
    RSC_COL, MP_BLUE, JARVIS_GREEN, GREY, apply_style, save,
)
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
ART = ROOT / "shared/output/routing_k_sweep_50seed/routing_k_sweep.json"
BOOT = ROOT / "shared/output/routing_k_sweep_50seed/routing_k_sweep_bootstrap.json"
OUT = ROOT / "writings/papers_2026/pegki_digital_discovery/figures/routing_k_sweep.png"

LINES = [("oracle", "best measured donor", GREY, "--"),
         ("learned_from_lp_fixed_score", "routing learned from the LP", JARVIS_GREEN, "-"),
         ("lp_fixed", "maximin LP and uniform routing", MP_BLUE, "-")]


def main():
    apply_style()
    d = json.load(open(ART, encoding="utf-8"))
    b = json.load(open(BOOT, encoding="utf-8"))["by_k"]
    ks = sorted(int(k) for k in d["by_k"])
    fig, ax = plt.subplots(figsize=(RSC_COL, 2.3))
    for key, label, col, ls in LINES:
        y = [d["by_k"][str(k)]["routing_mean"][key] for k in ks]
        lo = [b[str(k)]["lo"][key] for k in ks]
        hi = [b[str(k)]["hi"][key] for k in ks]
        ax.fill_between(ks, lo, hi, color=col, alpha=0.18, lw=0)
        ax.plot(ks, y, ls, color=col, marker="o", ms=3.5, label=label)
    k2 = d["by_k"]["2"]["routing_mean"]
    gap = k2["learned_from_lp_fixed_score"] - k2["lp_fixed"]
    ax.annotate(f"+{gap:.2f}", xy=(2, (k2["learned_from_lp_fixed_score"] + k2["lp_fixed"]) / 2),
                xytext=(5, 0), textcoords="offset points", va="center", fontsize=7,
                color=JARVIS_GREEN)
    ax.set_xticks(ks)
    ax.set_xticklabels([f"{k}\n(pairs)" if k == 2 else str(k) for k in ks])
    ax.set_xlabel("Candidates ranked per query event, k")
    ax.set_ylabel("Early transfer gain (Kendall $\\tau$)")
    ax.set_ylim(0.08, 0.50)
    ax.legend(loc="upper right")
    save(fig, OUT)


if __name__ == "__main__":
    main()
