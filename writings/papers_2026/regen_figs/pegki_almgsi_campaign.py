#!/usr/bin/env python3
"""PEGKi figure: outcome of the computed Al-Mg-Si discovery campaigns.

Panels (a-c): share of receiver compute needed to find eight of the ten most
stable compounds for the PBE, HSE06 and r2SCAN receivers, per donor, with
fixed trust (open) and learned trust (filled); no-donor baseline in grey.
Whiskers: 5-95% range over the 20 runs with MD-noise draws (arms without
an MD donor have no draw-to-draw spread). Panel (d): PRISMS given each
assumed fidelity order of the three potentials against the donor ordering
chosen by the transfer score, PBE and HSE06 receivers, 95% bootstrap
intervals over runs.

Data: compute/almgsi_campaign/results/campaigns/{qe_pbe_all,qe_hse06_all,
gpaw_r2scan_all}/summary.json and results/analysis/prisms_comparison/
prisms_comparison_{qe_pbe,qe_hse06}.json. No calculation is rerun.
Size and fonts: RSC full width, paper_fig_style fonts.
"""

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from paper_fig_style import (  # noqa: E402
    RSC_FULL, ALMGSI_COLOURS, GREY, INK, apply_style, save,
)
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
LIVE = ROOT / "compute" / "almgsi_campaign" / "results"
OUT = ROOT / "writings/papers_2026/pegki_digital_discovery/figures/almgsi_campaign.png"

DONORS = [("mace_relax", "MACE-MP-0 0 K"), ("archive_mp", "MP archive"),
          ("lammps_meam_md", "MEAM MD"), ("gpumd_nep89_md", "NEP89 MD")]
RECEIVERS = [("qe_pbe_all", "(a) PBE receiver", None),
             ("qe_hse06_all", "(b) HSE06 receiver", ("qe_pbe", "PBE")),
             ("gpaw_r2scan_all", "(c) r$^2$SCAN receiver", ("qe_pbe", "PBE"))]


def load(p):
    return json.loads(Path(p).read_text(encoding="utf-8"))


def dot(ax, y, m, col, filled):
    v = 100 * m["mean"]
    lo, hi = 100 * m["p05"], 100 * m["p95"]
    ax.errorbar([v], [y], xerr=[[max(v - lo, 0)], [max(hi - v, 0)]], fmt="o", ms=4.2,
                mfc=col if filled else "white", mec=col, mew=1.1,
                ecolor=col, elinewidth=0.9, capsize=1.5)


def main():
    apply_style()
    fig = plt.figure(figsize=(RSC_FULL, 4.3))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 0.95])
    a0 = fig.add_subplot(gs[0, 0])
    axes = [a0, fig.add_subplot(gs[0, 1], sharey=a0), fig.add_subplot(gs[0, 2], sharey=a0)]
    rows = [("surrogate", "No donor")] + DONORS + [("qe_pbe", "PBE")]
    for ax, (tag, title, extra) in zip(axes, RECEIVERS):
        arms = load(LIVE / f"campaigns/{tag}/summary.json")["arms"]
        for i, (key, _lab) in enumerate(rows):
            col = ALMGSI_COLOURS.get(key, GREY)
            if f"fixed:{key}|k4" not in arms and key != "surrogate":
                continue
            if key == "surrogate":
                dot(ax, i, arms["surrogate|k4"]["cost_frac_to_8of10"], col, True)
                continue
            dot(ax, i - 0.17, arms[f"fixed:{key}|k4"]["cost_frac_to_8of10"], col, False)
            dot(ax, i + 0.17, arms[f"learned:{key}|k4"]["cost_frac_to_8of10"], col, True)
        ax.set_yticks(range(len(rows)))
        ax.set_yticklabels([lab for _, lab in rows])
        ax.set_ylim(len(rows) - 0.5, -0.6)
        ax.set_xlim(0, 100)
        ax.set_xticks([0, 50, 100])
        ax.set_title(title, loc="left", fontsize=8)
        if ax is not a0:
            ax.tick_params(labelleft=False)
        ax.grid(axis="y", visible=False)
        base = 100 * arms["surrogate|k4"]["cost_frac_to_8of10"]["mean"]
        ax.axvline(base, color=GREY, lw=0.6, ls=":")
    axes[1].set_xlabel("Receiver compute to find 8 of the 10 most stable (%)")
    h = [plt.Line2D([], [], marker="o", ls="", mfc="white", mec=INK, label="fixed trust"),
         plt.Line2D([], [], marker="o", ls="", mfc=INK, mec=INK, label="learned trust")]
    a0.legend(handles=h, loc="lower right", fontsize=7, handletextpad=0.2)

    # (d) PRISMS by assumed source order against the scored donor
    ax = fig.add_subplot(gs[1, :])
    rows = [("prisms_3fid:cost", "PRISMS, order by cost per query"),
            ("prisms_3fid:class", "PRISMS, order by model class"),
            ("prisms_3fid:measured", "PRISMS, order by measured agreement"),
            ("prisms_single:mace_relax", "PRISMS, top-scored potential only"),
            ("surrogate", "No donor"),
            ("fixed:mace_relax", "Top-scored donor ordering, fixed trust")]
    for j, (rec, lab, off, mk) in enumerate([("qe_pbe", "PBE receiver", -0.15, "o"),
                                            ("qe_hse06", "HSE06 receiver", 0.15, "s")]):
        arms = load(LIVE / f"analysis/prisms_comparison/prisms_comparison_{rec}.json")[
            "receivers"][rec]["arms"]
        col = INK if j == 0 else ALMGSI_COLOURS["qe_hse06"]
        for i, (key, _) in enumerate(rows):
            a = arms[key]
            v = 100 * a["cost_frac_to_8of10"]
            lo, hi = (100 * np.array(a["c8_ci"])) if a.get("c8_ci") else (v, v)
            ax.errorbar([v], [i + off], xerr=[[max(v - lo, 0)], [max(hi - v, 0)]], fmt=mk, ms=4,
                        mfc=col if j == 0 else "white", mec=col, ecolor=col,
                        elinewidth=0.9, capsize=1.5, label=lab if i == 0 else None)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([lab for _, lab in rows])
    ax.set_ylim(len(rows) - 0.5, -0.6)
    ax.set_xlim(0, 100)
    ax.set_xlabel("Receiver compute to find 8 of the 10 most stable (%)")
    ax.set_title("(d) Assumed fidelity order against the scored donor (165 compounds)",
                 loc="left", fontsize=8)
    ax.grid(axis="y", visible=False)
    ax.legend(loc="lower right", fontsize=7)
    save(fig, OUT)


if __name__ == "__main__":
    main()
