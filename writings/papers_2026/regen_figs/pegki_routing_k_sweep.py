#!/usr/bin/env python3
"""PEGKi figure: routing rules across query arity k.

Panels: one column-width panel. x = k (2 = pairwise arm, 3-6 = k-way arms
of matched-seed paired campaigns); y = realised early transfer gain on
held-out receivers (leave one receiver out). Lines: maximin linear program
(identical to uniform), routing learned from the linear-program cold start,
best measured donor. From shared/output/routing_k_sweep_50seed/routing_k_sweep.json.

Size and fonts: RSC column (3.27 in), paper_fig_style fonts.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from paper_fig_style import (  # noqa: E402
    RSC_COL, MP_BLUE, JARVIS_GREEN, GREY, INK, apply_style, save,
)
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
ART = ROOT / "shared/output/routing_k_sweep_50seed/routing_k_sweep.json"
OUT = ROOT / "writings/papers_2026/pegki_digital_discovery/figures/routing_k_sweep.png"

LINES = [("oracle", "best measured donor", GREY, "--"),
         ("learned_from_lp_fixed_score", "learned from the LP", JARVIS_GREEN, "-"),
         ("lp_fixed", "maximin LP (= uniform)", MP_BLUE, "-")]


def main():
    apply_style()
    d = json.load(open(ART, encoding="utf-8"))
    ks = sorted(int(k) for k in d["by_k"])
    fig, ax = plt.subplots(figsize=(RSC_COL, 2.3))
    for key, label, col, ls in LINES:
        y = [d["by_k"][str(k)]["routing_mean"][key] for k in ks]
        ax.plot(ks, y, ls, color=col, marker="o", ms=3.5, label=label)
    ax.set_xticks(ks)
    ax.set_xticklabels([f"{k}\n(pairs)" if k == 2 else str(k) for k in ks])
    ax.set_xlabel("Candidates ranked per query event, k")
    ax.set_ylabel("Early transfer gain (Kendall $\\tau$)")
    ax.legend(loc="upper right")
    save(fig, OUT)


if __name__ == "__main__":
    main()
