#!/usr/bin/env python
"""Archived MP r2SCAN as a third DFT level on the pool compounds MP holds.

Compares MP r2SCAN formation energies (fetch_mp_r2scan.py) with the computed
sources (QE PBE, QE HSE06, MACE-MP-0, MEAM and NEP89 relaxations) and with MP
GGA/GGA+U, on the compounds each pair shares: Kendall tau of the orderings and
mean absolute difference. Output: results/analysis/mp_r2scan_study.{json,md}
"""
import csv
import json
from pathlib import Path

import numpy as np
from scipy.stats import kendalltau

HERE = Path(__file__).resolve().parent
RES = HERE / "results"
TABLES = {"QE PBE": "qe_pbe", "QE HSE06": "qe_hse06", "MACE-MP-0 0 K": "mace_relax",
          "MEAM 0 K": "lammps_meam_relax", "NEP89 0 K": "gpumd_nep89_relax"}
JANAF_MG2SI = -77.8 / 96.485 / 3   # eV/atom, calorimetric (NIST-JANAF), as lab_receivers_study.py


def table(name):
    out = {}
    with open(RES / f"table_{name}.csv") as f:
        for r in csv.DictReader(f):
            if r["status"] == "ok" and r["property"] == "formation_energy" and r["value"]:
                out[r["cid"]] = float(r["value"])
    return out


def main():
    mp = json.loads((RES / "mp_r2scan" / "mp_r2scan_formation.json").read_text())["records"]
    src = {"MP r2SCAN": {c: r["r2scan_ef"] for c, r in mp.items() if "r2scan_ef" in r},
           "MP GGA/GGA+U": {c: r["gga_ef"] for c, r in mp.items() if "gga_ef" in r}}
    src.update({k: table(v) for k, v in TABLES.items()})
    names = list(src)
    rows, res = [], {}
    for rec in ("MP r2SCAN", "QE PBE", "QE HSE06"):
        for don in names:
            if don == rec:
                continue
            common = sorted(set(src[rec]) & set(src[don]) & set(src["MP r2SCAN"]))
            common = [c for c in common if c.split("_")[-1] not in ("Al", "Mg", "Si")]  # drop elements
            if len(common) < 5:
                continue
            a = np.array([src[don][c] for c in common]); b = np.array([src[rec][c] for c in common])
            t = float(kendalltau(a, b).statistic)
            res[f"{don}->{rec}"] = {"tau": t, "mae_eV_atom": float(np.mean(np.abs(a - b))), "n": len(common)}
            rows.append(f"| {don} | {rec} | {len(common)} | {t:+.2f} | {np.mean(np.abs(a - b)):.3f} |")
    mg2si = {k: src[k].get("mp_1367_Mg2Si") for k in names}
    lines = ["# MP r2SCAN (archived) against the computed sources", "",
             f"{len(src['MP r2SCAN'])} pool compounds carry an MP r2SCAN formation energy; "
             "comparisons use compounds (elements dropped) that both sources and MP r2SCAN hold.", "",
             "| Donor | Receiver | n | Kendall tau | MAD (eV/atom) |", "|---|---|---|---|---|", *rows, "",
             f"Mg2Si (mp-1367) formation energy, eV/atom; calorimetric {JANAF_MG2SI:+.3f}:", ""]
    lines += [f"- {k}: {v:+.3f}" for k, v in mg2si.items() if v is not None]
    (RES / "analysis").mkdir(exist_ok=True)
    (RES / "analysis" / "mp_r2scan_study.json").write_text(json.dumps({"pairs": res, "mg2si": mg2si,
                                                                      "mg2si_calorimetric": JANAF_MG2SI}, indent=1))
    (RES / "analysis" / "mp_r2scan_study.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
