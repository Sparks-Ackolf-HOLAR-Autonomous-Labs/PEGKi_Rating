#!/usr/bin/env python
"""Finite-size check of the vacancy-solute binding energies.

Collects results/analysis/vacancy_binding{,_n3,_n4,_n5}.json (vacancy_study.py
--analyse at VAC_N = 2, 3, 4, 5: 32, 108, 256, 500 sites) and reports, per
source, the binding energy of each complex at every size, the change from the
32-site cell to the largest cell, and the Kendall tau of the binding ordering
between sizes and against PBE at each size. In the 32-site cell a solute in the
second-neighbour shell of the vacancy sits at the same distance from the
vacancy and from one of its periodic images, so the 2nn items are expected to
move most. Output: results/analysis/vacancy_size_study.{json,md}
"""
import json
from pathlib import Path

from scipy.stats import kendalltau

HERE = Path(__file__).resolve().parent
A = HERE / "results" / "analysis"
SIZES = {32: "", 108: "_n3", 256: "_n4", 500: "_n5"}
LABEL = {"mace_relax": "MACE-MP-0", "lammps_meam_relax": "MEAM", "gpumd_nep89_relax": "NEP89", "qe_pbe": "QE PBE"}


def main():
    data = {}
    for n, suf in SIZES.items():
        p = A / f"vacancy_binding{suf}.json"
        if p.exists():
            data[n] = json.loads(p.read_text())["sources"]
    items = list(next(iter(data.values()))["mace_relax"]["binding_eV"])
    out, lines = {}, ["# Vacancy binding against supercell size (eV, positive = bound)", ""]
    for s, lab in LABEL.items():
        sizes = [n for n in data if data[n].get(s, {}).get("binding_eV")]
        if not sizes:
            continue
        b = {n: data[n][s]["binding_eV"] for n in sizes}
        out[s] = {"binding_eV": {str(n): b[n] for n in sizes}}
        lines += [f"## {lab}", "", "| Item | " + " | ".join(f"{n} sites" for n in sizes) + " |",
                  "|---|" + "---|" * len(sizes)]
        for it in items:
            lines.append(f"| {it} | " + " | ".join(f"{b[n][it]:+.3f}" if it in b[n] else "n/a" for n in sizes) + " |")
        if len(sizes) > 1:
            lo, hi = sizes[0], sizes[-1]
            common = [i for i in items if i in b[lo] and i in b[hi]]
            d = {i: b[hi][i] - b[lo][i] for i in common}
            t = float(kendalltau([b[lo][i] for i in common], [b[hi][i] for i in common]).statistic)
            out[s].update({"change_small_to_large": d, "tau_small_vs_large": t, "sizes": sizes})
            big = max(d, key=lambda i: abs(d[i]))
            lines += ["", f"{lo} to {hi} sites: largest change {big} {d[big]:+.3f} eV; "
                      f"Kendall tau of the ordering between the two sizes {t:+.2f}."]
        lines.append("")
    # donor orderings against PBE at each size where PBE exists
    taus = {}
    for n in data:
        pb = data[n].get("qe_pbe", {}).get("binding_eV", {})
        if len(pb) < 4:
            continue
        for s in ("mace_relax", "lammps_meam_relax", "gpumd_nep89_relax"):
            d = data[n][s]["binding_eV"]
            common = [i for i in items if i in d and i in pb]
            taus[f"{s}@{n}"] = float(kendalltau([d[i] for i in common], [pb[i] for i in common]).statistic)
    if taus:
        lines += ["## Donor orderings against QE PBE at the same size", ""]
        lines += [f"- {LABEL[k.split('@')[0]]}, {k.split('@')[1]} sites: tau {v:+.2f}" for k, v in taus.items()]
    out["tau_vs_pbe"] = taus
    (A / "vacancy_size_study.json").write_text(json.dumps(out, indent=1))
    (A / "vacancy_size_study.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
