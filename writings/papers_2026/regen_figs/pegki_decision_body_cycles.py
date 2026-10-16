#!/usr/bin/env python3
"""Regenerate the PEGKi decision-body loop figure at true print size.

Loops are measured on a graded scale, not as yes/no cycles, and among
the decision bodies (learners, acquisition rules, reuse strategies), not
only among material classes. Panels: two, full RSC width:
  (a) loop share (curl plus harmonic energy of the comparison flow) with
      its interval, for every decision-body tournament and every
      material-class tournament; for the reuse strategies an open marker
      shows the loop share left inside a single receiver, so the reader
      sees how much of the loop comes from pooling environments.
  (b) commitment regret in Kendall tau: the loss from committing to the
      body that won in one context and using it in another, across
      environments versus across campaign stages.
All numbers come from
shared/output/decision_body_cycles/decision_body_cycles.json.

Size and fonts: RSC full width (6.73 in), fonts per paper_fig_style.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from paper_fig_style import (  # noqa: E402
    ACCENT_RED, GREY, INK, MP_BLUE, RSC_FULL, SOURCE_NAMES, apply_style, save)
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
ART = ROOT / "shared/output/decision_body_cycles/decision_body_cycles.json"
OUT = ROOT / "writings/papers_2026/pegki_digital_discovery/figures/decision_body_cycles.png"

# Main-text version: the first-pass tournaments (8 properties, 8 seed
# means) are not plotted; they are reported in Annex S4.9 text only.
BODY_ROWS = [
    ("ML models, all leaderboard tasks", "ML models, 43 tasks"),
    ("Reuse strategies, pairwise queries", "Reuse strategies, pairwise"),
    ("Reuse strategies, 4-way events", "Reuse strategies, 4-way"),
    ("Reuse strategies, 6-way events", "Reuse strategies, 6-way"),
    ("Acquisition rules, long horizon (exp 27b)", "Acquisition rules, 800 campaigns"),
]
CLASS_ORDER = ["expt_zhuo", "mp_gga", "mp_gga_expsynth", "mp_r2scan",
               "jarvis_optb88vdw", "jarvis_mbj", "snumat_gga", "snumat_hse"]
REGRET_ROWS = [
    ("Reuse strategies, pairwise queries", "receiver", "Reuse,\npairwise"),
    ("Reuse strategies, 4-way events", "receiver", "Reuse,\n4-way"),
    ("Reuse strategies, 6-way events", "receiver", "Reuse,\n6-way"),
    ("Acquisition rules, long horizon (exp 27b)", "held-out source", "Acquisition\nrules"),
]
BODY_COLOUR = ACCENT_RED
CLASS_COLOUR = MP_BLUE


def main():
    apply_style()
    d = json.load(open(ART))
    db, mc = d["decision_bodies"], d["material_classes"]
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(RSC_FULL, 3.1),
                                 gridspec_kw={"width_ratios": [1.35, 1.0]})

    # (a) loop share
    labels, y = [], 0
    ys = []
    for key, lab in BODY_ROWS:
        h = db[key]
        v, lo, hi = h["hodge"]["loop"], *h["loop_share_ci90"]
        ax.errorbar(v, y, xerr=[[v - min(lo, v)], [max(hi, v) - v]], fmt="o",
                    color=BODY_COLOUR, ms=4, capsize=2, lw=1)
        ctx = h["contexts"].get("receiver")
        if ctx:
            w = [x for x in ctx["within_context_loop_share"].values()]
            ax.plot(np.mean(w), y, "o", mfc="white", mec=BODY_COLOUR, ms=4, mew=1)
        labels.append(lab)
        ys.append(y)
        y += 1
    y += 0.6
    for src in CLASS_ORDER:
        h = mc[src]
        v, lo, hi = h["hodge"]["loop"], *h["loop_share_ci90"]
        ax.errorbar(v, y, xerr=[[v - min(lo, v)], [max(hi, v) - v]], fmt="s",
                    color=CLASS_COLOUR, ms=3.5, capsize=2, lw=1)
        labels.append("Classes, " + SOURCE_NAMES[src].replace(" (ICSD-matched)", " (ICSD)"))
        ys.append(y)
        y += 1
    ax.set_yticks(ys)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlim(0, 0.24)
    ax.set_xlabel("Loop share of the comparison record")
    ax.set_title("(a) Loops: decision bodies vs material classes", loc="left")
    ax.plot([], [], "o", color=BODY_COLOUR, ms=4, label="decision bodies, pooled")
    ax.plot([], [], "o", mfc="white", mec=BODY_COLOUR, ms=4, label="within one receiver")
    ax.plot([], [], "s", color=CLASS_COLOUR, ms=3.5, label="material classes")
    ax.legend(loc="lower right", handletextpad=0.3)
    ax.grid(axis="y", visible=False)

    # (b) commitment regret
    xs = np.arange(len(REGRET_ROWS))
    wbar = 0.36
    for i, (key, env, lab) in enumerate(REGRET_ROWS):
        ce = db[key]["contexts"][env]["commitment_regret"]
        cs = db[key]["contexts"].get("stage", {}).get("commitment_regret")
        bx.bar(i - wbar / 2, ce["mean"], wbar, color=BODY_COLOUR)
        bx.plot([i - wbar / 2] * 2, [ce["mean"], ce["max"]], color=INK, lw=0.8)
        bx.plot(i - wbar / 2, ce["max"], "_", color=INK, ms=6)
        if cs:
            bx.bar(i + wbar / 2, cs["mean"], wbar, color=GREY)
            bx.plot([i + wbar / 2] * 2, [cs["mean"], cs["max"]], color=INK, lw=0.8)
            bx.plot(i + wbar / 2, cs["max"], "_", color=INK, ms=6)
    bx.set_xticks(xs)
    bx.set_xticklabels([r[2] for r in REGRET_ROWS])
    bx.set_ylabel(r"Commitment regret (Kendall $\tau$)")
    bx.set_title("(b) Cost of committing", loc="left")
    handles = [Patch(color=BODY_COLOUR, label="across environments"),
               Patch(color=GREY, label="across budget stages"),
               Line2D([], [], marker="_", color=INK, ls="none", ms=6,
                      label="worst pair of contexts")]
    bx.set_ylim(0, 0.62)
    bx.legend(handles=handles, loc="upper right")
    bx.grid(axis="x", visible=False)
    save(fig, OUT)


if __name__ == "__main__":
    main()
