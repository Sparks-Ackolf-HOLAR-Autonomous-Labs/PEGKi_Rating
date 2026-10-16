#!/usr/bin/env python
"""CALPHAD as a receiver: the equilibrium (convex-hull) energy per composition.

CALPHAD describes phases, not atomic orderings, so it returns one value per
composition: the Gibbs energy of the stable phase assemblage. The matching
quantity from any computed source is its lower convex hull over the pool,
evaluated at each composition: the lowest energy that a mixture of the
source's computed structures reaches there. This script computes that hull
energy for every source at the 27 distinct pool compositions and reports
Kendall tau of each source against CALPHAD (COST 507) at 300, 453 and 803 K.
MD sources at temperature are compared with CALPHAD at the same temperature;
0 K sources are compared at every temperature.

Output: results/analysis/hull_receiver.{json,md}
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.optimize import linprog
from scipy.stats import kendalltau

import receiver_extension_study as rx

HERE = Path(__file__).resolve().parent
RES = HERE / "results"
ZERO_K = ("mace_relax", "qe_pbe", "qe_hse06", "archive_mp", "qe_r2scan")


def hull_energy(points: dict[str, tuple[np.ndarray, float]], x: np.ndarray) -> float:
    """Lowest energy of a mixture of the source's structures (plus elements at 0) at fractions x."""
    X = [p[0] for p in points.values()] + [np.eye(3)[i] for i in range(3)]
    E = [p[1] for p in points.values()] + [0.0, 0.0, 0.0]
    A = np.vstack([np.array(X).T, np.ones(len(X))])
    r = linprog(E, A_eq=A, b_eq=np.append(x, 1.0), bounds=(0, None), method="highs")
    return float(r.fun) if r.success else float("nan")


def main():
    src, idx = rx.load_src()
    ts = rx.ts
    frac = {}
    for s in src.values():
        for c in s:
            if c not in frac:
                frac[c] = ts.fractions_noase(c)
    comps = {}
    for c, x in frac.items():
        comps.setdefault(tuple(np.round(x, 6)), c)
    out, lines = {}, ["# CALPHAD (COST 507) as receiver: hull energy per composition", "",
                      "| T | Source | tau vs CALPHAD | n compositions |", "|---|---|---|---|"]
    for T in (300, 453, 803):
        cal = src.get(f"calphad_{T}")
        if not cal:
            continue
        ref = {k: cal[c]["value"] for k, c in comps.items() if c in cal}
        sources = list(ZERO_K) + ([f"mace_md_{T}" if T != 300 else "mace_md_300"]
                                  + ([f"lammps_meam_md_{T}", f"gpumd_nep89_md_{T}"] if T != 300
                                     else ["lammps_meam_md", "gpumd_nep89_md"])) + [f"qha_{T}"]
        out[T] = {}
        for s in sources:
            if s not in src:
                continue
            pts = {c: (frac[c], v["value"]) for c, v in src[s].items()}
            h = {k: hull_energy(pts, np.array(k)) for k in ref}
            ks = [k for k in ref if np.isfinite(h[k])]
            if len(ks) < 5:
                continue
            t = kendalltau([h[k] for k in ks], [ref[k] for k in ks]).statistic
            out[T][s] = {"tau": float(t), "n": len(ks)}
            lines.append(f"| {T} K | {rx.LABEL.get(s, s)} | {t:+.2f} | {len(ks)} |")
    (RES / "analysis" / "hull_receiver.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    (RES / "analysis" / "hull_receiver.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
