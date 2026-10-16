#!/usr/bin/env python3
"""PEGKi figure: solute-cluster binding energies of each computed source
against PBE (Al-Mg-Si).

Data: compute/almgsi_campaign/results/analysis/cluster_transfer.json.
One panel per source; x = PBE binding energy, y = source binding energy, both
in eV (positive = attractive). Dotted lines at zero mark attraction versus
repulsion; panel titles carry Kendall tau against PBE with its 90% bootstrap
interval and tau on the ordered compounds for contrast. All panels share one
y axis so the size of each source's binding energies compares directly with
PBE; the dashed line is y = x.
Size and fonts: paper_fig_style, RSC full width at print size.
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from paper_fig_style import RSC_FULL, MP_BLUE, INK, GREY, ACCENT_RED, apply_style, save  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "compute" / "almgsi_campaign" / "results" / "analysis" / "cluster_transfer.json"
OUT = ROOT / "writings" / "papers_2026" / "pegki_digital_discovery" / "figures" / "cluster_transfer.png"
PANELS = ["mace_relax", "mace_md", "gpumd_nep89_md", "lammps_meam_md"]


def main():
    apply_style()
    d = json.loads(DATA.read_text())
    pbe = d["values"]["qe_pbe"]
    fig, axes = plt.subplots(1, 4, figsize=(RSC_FULL, 2.2), sharey=True)
    for ax, s in zip(axes.flat, PANELS):
        r, v = d["sources"][s], d["values"][s]
        c = sorted(set(v) & set(pbe))
        x = np.array([pbe[k]["E_b_eV"] for k in c])
        y = np.array([v[k]["E_b_eV"] for k in c])
        col = ACCENT_RED if r["kendall_tau"] < 0 else MP_BLUE
        ax.scatter(x, y, s=12, color=col, lw=0, zorder=3)
        ax.axhline(0, color=GREY, lw=0.6, ls=":")
        ax.axvline(0, color=GREY, lw=0.6, ls=":")
        ax.plot([-0.01, 0.14], [-0.01, 0.14], color=GREY, lw=0.7, ls="--", zorder=2)
        ax.set_title(f"{r['label']}\nclusters $\\tau$ {r['kendall_tau']:+.2f}\n"
                     f"compounds $\\tau$ {r['compound_tau_vs_pbe']:+.2f}", fontsize=7)
        ax.tick_params(labelsize=7)
    for ax in axes:
        ax.set_xlabel("PBE binding (eV)", fontsize=7.5)
    for ax in axes[:1]:
        ax.set_ylabel("source binding (eV)", fontsize=7.5)
    fig.tight_layout(pad=0.3, h_pad=0.6, w_pad=0.6)
    save(fig, OUT)


if __name__ == "__main__":
    main()
