#!/usr/bin/env python3
"""PEGKi figure: the observed Al-Mg-Si precipitation sequence as a receiver.

Data: compute/almgsi_campaign/results/analysis/lab_sequence_study.json
(energy per solute atom of every staged item, per computed source). One small
panel per source; x = stage of the observed sequence (Edwards et al. 1998),
y = energy per solute atom relative to fcc Al, hcp Mg, diamond Si. Points are
items, the line joins stage medians, the band spans the stage minimum to
maximum. A source that reproduces the sequence descends from left to right.
All panels share one y axis so the energy scales compare directly.
Panel title carries the agreement (-Kendall tau-b) and its 90% bootstrap CI.
Size and fonts: paper_fig_style, RSC full width at print size.
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from paper_fig_style import RSC_FULL, MP_BLUE, INK, GREY, apply_style, save  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "compute" / "almgsi_campaign" / "results" / "analysis" / "lab_sequence_study.json"
OUT = ROOT / "writings" / "papers_2026" / "pegki_digital_discovery" / "figures" / "lab_sequence.png"
sys.path.insert(0, str(ROOT / "compute" / "almgsi_campaign"))
from lab_sequence_study import STAGE, NSTAGE  # noqa: E402

PANELS = ["qe_pbe", "mace_relax", "mace_md", "archive_mp", "gpumd_nep89_md", "lammps_meam_md"]
TICKS = ["solid\nsol.", "Mg, Si\nclust.", "Mg-Si\nco-cl.", "$\\beta''$",
         "$\\beta'$\nU1, U2", "$\\beta$"]


def main():
    apply_style()
    d = json.loads(DATA.read_text())
    fig, axes = plt.subplots(2, 3, figsize=(RSC_FULL, 3.9), sharex=True, sharey=True)
    for ax, s in zip(axes.flat, PANELS):
        es, r = d["values"].get(s, {}), d["sources"].get(s, {})
        xs, ys = [], []
        for cid, v in es.items():
            if cid in STAGE:
                xs.append(STAGE[cid]); ys.append(v)
        xs, ys = np.array(xs), np.array(ys)
        st = sorted(set(xs.tolist()))
        lo = [ys[xs == k].min() for k in st]
        hi = [ys[xs == k].max() for k in st]
        med = [np.median(ys[xs == k]) for k in st]
        ax.fill_between(st, lo, hi, color=MP_BLUE, alpha=0.15, lw=0)
        ax.plot(st, med, color=MP_BLUE, lw=1.6)
        jit = (np.random.default_rng(3).random(len(xs)) - 0.5) * 0.18
        ax.scatter(xs + jit, ys, s=10, color=INK, zorder=3, lw=0)
        ax.axhline(0, color=GREY, lw=0.6, ls=":")
        if "agreement" in r:
            ci = r["ci90"]
            ttl = f"{r['label']}\nagreement {r['agreement']:+.2f} [{ci[0]:+.2f}, {ci[1]:+.2f}]"
        else:
            ttl = r.get("label", s)
        ax.set_title(ttl, fontsize=7.5)
        ax.set_xlim(-0.5, NSTAGE - 0.5)
        ax.set_xticks(range(NSTAGE))
        ax.set_xticklabels(TICKS, fontsize=6.0)
    for ax in axes[:, 0]:
        ax.set_ylabel("energy per solute atom (eV)")
    fig.tight_layout(h_pad=0.6, w_pad=0.6)
    save(fig, OUT)


if __name__ == "__main__":
    main()
