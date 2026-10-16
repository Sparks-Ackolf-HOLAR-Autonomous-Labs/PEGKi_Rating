#!/usr/bin/env python3
"""Regenerate the PEGKi bridge-budget figure at true print size.

Panels: three, one per source pair, matching the caption of
Fig. \\ref{fig:bridge}: offset sign stability and shared-set tau power
versus the number of shared bridge compositions n, from 300 seeded
subsamples per size. Vertical marks at n=3 (the Al-Mg-Si overlap) and
at the artifact's n_star per pair (8, 8, 12). All numbers come from
shared/output/bridge_budget/bridge_budget.json:
  pairs[*].curve[*].{n, sign_stability, tau_power}
  pairs[*].{n_star, n_shared, tau_full}
  config.{sign_stability, tau_power, resamples}   (threshold lines)

Size and fonts: RSC full-width figure (6.73 in), fonts per paper_fig_style.
Labels use "data source".
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from paper_fig_style import (  # noqa: E402
    MP_BLUE, GREY, apply_style, save,
)
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
ART = ROOT / "shared/output/bridge_budget/bridge_budget.json"
OUT = ROOT / "writings/papers_2026/pegki_digital_discovery/figures/bridge_budget.png"

VERMILLION = "#D55E00"
GREEN = "#0F7B3F"

PAIR_TITLES = {
    ("jarvis_optb88vdw", "jarvis_mbj"): "JARVIS OptB88vdW vs TB-mBJ\n(same-database fidelity pair)",
    ("snumat_gga", "snumat_hse"): "SNUMAT GGA vs HSE06\n(same-database GGA/HSE pair)",
    ("mp_gga", "expt_zhuo"): "MP GGA vs Experiment (Zhuo)\n(DFT vs experiment)",
}


def main():
    apply_style()
    data = json.load(open(ART))
    cfg = data["config"]

    fig, axes = plt.subplots(1, 3, figsize=(6.73, 2.5), sharey=True)

    for ax, pr in zip(axes, data["pairs"]):
        ns = [c["n"] for c in pr["curve"]]
        sign = [c["sign_stability"] for c in pr["curve"]]
        power = [c["tau_power"] for c in pr["curve"]]
        ax.set_xscale("log")
        ax.axhline(cfg["sign_stability"], color=MP_BLUE, lw=0.8, ls=":",
                   alpha=0.8)
        ax.axhline(cfg["tau_power"], color=VERMILLION, lw=0.8, ls=":",
                   alpha=0.8)
        ax.plot(ns, sign, marker="o", ms=3.2, color=MP_BLUE,
                label="offset sign stability")
        ax.plot(ns, power, marker="s", ms=3.2, color=VERMILLION,
                label="shared-set $\\tau$ power ($p<0.05$)")
        ax.axvline(3, color=GREY, lw=0.9, ls="--")
        ax.axvline(pr["n_star"], color=GREEN, lw=0.9, ls="--")
        ax.set_xlim(2.4, 235)
        ax.text(2.67, 0.5, "$n=3$ (Al-Mg-Si)", fontsize=6.5,
                color="#555555", ha="center", va="center", rotation=90)
        ax.text(pr["n_star"] * 1.09, 0.5, f"$n^{{*}}={pr['n_star']}$",
                fontsize=6.5, color=GREEN, ha="left", va="center")
        key = (pr["channel_a"], pr["channel_b"])
        ax.set_title(PAIR_TITLES[key] + "\n"
                     f"{pr['n_shared']} shared, full-set "
                     f"$\\tau = {pr['tau_full']:.2f}$", fontsize=7.5)
        ax.set_xlabel("bridge compositions $n$ (log)")
        ax.set_xticks([3, 10, 30, 100])
        ax.set_xticklabels(["3", "10", "30", "100"])
        ax.set_ylim(-0.05, 1.09)
    axes[0].set_ylabel(f"fraction of {cfg['resamples']} subsamples")
    axes[0].legend(loc="lower right", fontsize=6.5, borderaxespad=0.4)

    save(fig, OUT)


if __name__ == "__main__":
    main()
