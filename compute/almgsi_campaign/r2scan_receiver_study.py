#!/usr/bin/env python
"""r2SCAN (GPAW, PAW) as a third DFT receiver on the Al-Mg-Si pool.

Reads the formation energies of every source from results/table_<source>.csv
(aggregate_sources.py) and compares each donor with the GPAW r2SCAN receiver on
the compounds both hold: Kendall tau, mean absolute difference, and how many of
the receiver's ten most stable compounds are among the donor's ten. Also checks
GPAW r2SCAN against the archived MP r2SCAN (VASP, PAW) on shared compounds, the
GPAW PBE code check against QE PBE, and Mg2Si against the calorimetric value.
Output: results/analysis/r2scan_receiver_study.{json,md}
"""
import csv
import json
from pathlib import Path

import numpy as np
from scipy.stats import kendalltau

HERE = Path(__file__).resolve().parent
RES = HERE / "results"
REC = "gpaw_r2scan"
DONORS = {"qe_pbe": "QE PBE", "qe_hse06": "QE HSE06", "mace_relax": "MACE-MP-0, 0 K",
          "mace_md": "MACE-MP-0 MD, 300 K", "gpumd_nep89_md": "NEP89 MD, 300 K",
          "lammps_meam_md": "MEAM MD, 300 K"}
JANAF_MG2SI = -77.822 / 96.485 / 3


def table(path):
    out = {}
    if not path.exists():
        return out
    with open(path) as f:
        for r in csv.DictReader(f):
            if r["status"] == "ok" and r["property"] == "formation_energy" and r["value"]:
                out[r["cid"]] = float(r["value"])
    return out


def elemental(cid):
    return cid.split("_")[-1] in ("Al", "Mg", "Si")


def compare(d, r):
    common = sorted(c for c in set(d) & set(r) if not elemental(c))
    a = np.array([d[c] for c in common]); b = np.array([r[c] for c in common])
    top_r = set(sorted(common, key=r.get)[:10]); top_d = set(sorted(common, key=d.get)[:10])
    return {"n": len(common), "tau": float(kendalltau(a, b).statistic),
            "mad_eV_atom": float(np.mean(np.abs(a - b))), "top10_overlap": len(top_r & top_d)}


def main():
    rec = table(RES / f"table_{REC}.csv")
    res = {"n_receiver": len(rec), "donors": {}}
    for s, lab in DONORS.items():
        d = table(RES / f"table_{s}.csv")
        if d:
            res["donors"][s] = {"label": lab, **compare(d, rec)}
    mp = json.loads((RES / "mp_r2scan" / "mp_r2scan_formation.json").read_text())["records"]
    mpr = {c: v["r2scan_ef"] for c, v in mp.items() if "r2scan_ef" in v}
    res["mp_r2scan_vs_gpaw"] = compare(mpr, rec)
    gp = json.loads((RES / "analysis" / "gpaw_r2scan_formation.json").read_text()).get("gpaw_pbe", {})
    gpbe = gp.get("formation_eV_atom", {})
    if gpbe:
        res["gpaw_pbe_vs_qe_pbe"] = compare(gpbe, table(RES / "table_qe_pbe.csv"))
    res["mg2si"] = {"gpaw_r2scan": rec.get("mp_1367_Mg2Si"), "mp_r2scan": mpr.get("mp_1367_Mg2Si"),
                    "calorimetric": JANAF_MG2SI}
    archived = {c for c, v in mp.items()}
    top = sorted((c for c in rec if not elemental(c)), key=rec.get)[:10]
    res["receiver_top10"] = top
    res["top10_not_archived"] = [c for c in top if c not in archived]
    lines = ["# GPAW r2SCAN as a third DFT receiver", "",
             f"{len(rec)} pool compounds with a GPAW r2SCAN formation energy (PAW, 700 eV, libxc r2SCAN).", "",
             "| Donor | n | Kendall tau | MAD (eV/atom) | top-10 overlap |", "|---|---|---|---|---|"]
    for s, v in res["donors"].items():
        lines.append(f"| {v['label']} | {v['n']} | {v['tau']:+.2f} | {v['mad_eV_atom']:.3f} | {v['top10_overlap']} |")
    m = res["mp_r2scan_vs_gpaw"]
    lines += ["", f"MP r2SCAN (VASP) against GPAW r2SCAN: tau {m['tau']:+.2f}, MAD {m['mad_eV_atom']:.3f} eV/atom (n = {m['n']})."]
    if gpbe:
        g = res["gpaw_pbe_vs_qe_pbe"]
        lines.append(f"GPAW PBE against QE PBE (code check): tau {g['tau']:+.2f}, MAD {g['mad_eV_atom']:.3f} eV/atom (n = {g['n']}).")
    mg = res["mg2si"]
    lines += ["", f"Mg2Si (eV/atom): GPAW r2SCAN {mg['gpaw_r2scan']:+.3f}, MP r2SCAN {mg['mp_r2scan']:+.3f}, "
              f"calorimetric {mg['calorimetric']:+.3f}.", "",
              f"Ten most stable r2SCAN compounds: {', '.join(top)}; not in the MP archive: {len(res['top10_not_archived'])}."]
    (RES / "analysis" / "r2scan_receiver_study.json").write_text(json.dumps(res, indent=1))
    (RES / "analysis" / "r2scan_receiver_study.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
