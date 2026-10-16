#!/usr/bin/env python
"""Vacancy-solute and vacancy-pair binding in fcc Al: a kinetics receiver.

Natural ageing of Al-Mg-Si alloys is controlled by vacancies: quenched-in
vacancies bind to Mg and Si solutes and to small solute clusters, and the
strength of that binding sets how fast clusters form at room temperature.
This script computes the binding energies of a vacancy to single Mg and Si
solutes (first and second neighbour) and to nearest-neighbour Mg-Si, Mg-Mg
and Si-Si pairs, with the same sources as the rest of the campaign:

  mace_relax         MACE-MP-0 (small), positions relaxed at fixed cell
  lammps_meam_relax  MEAM (Jelinek et al. 2012), positions relaxed at fixed cell
  gpumd_nep89_relax  NEP89, positions relaxed at fixed cell (GPUMD FIRE)
  qe_pbe             Quantum ESPRESSO PBE, positions relaxed at fixed cell

Cells are NxNxN conventional fcc Al (4N^3 sites; N = 2, 32 sites, by default;
VAC_N = 3 or 4 selects the 108- and 256-site size check) at each source's own 0 K
lattice parameter (QE: a_PBE(Al) from results/qe_pbe_al_lattice.json). The
binding energy is positive when the complex is bound:

  E_b(X-V)    = E(Al31 X) + E(Al31 V) - E(Al30 X V) - E(Al32)
  E_b(XY-V)   = E(Al30 X Y) + E(Al31 V) - E(Al29 X Y V) - E(Al32)

Usage:
  python vacancy_study.py --make            write the 14 cells to candidates_vac/
  VAC_N=3 python vacancy_study.py ...       same for 3x3x3 cells (candidates_vac_n3/, results/vacancy_n3/)
  python vacancy_study.py --source mace_relax|lammps_meam_relax|gpumd_nep89_relax
  python vacancy_study.py --source qe_pbe   (QE_NP ranks; resumable; long)
  DFT runs on CHPC guest partitions, not the DGX: qe/vacancy_qe.sbatch
  python vacancy_study.py --analyse         results/analysis/vacancy_binding.{json,md}
Per-cell results: results/vacancy/<source>/<cid>.json
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from itertools import combinations
from pathlib import Path

import numpy as np
from ase.build import bulk
from ase.io import read, write

HERE = Path(__file__).resolve().parent
ROOT = Path(os.environ.get("LIVE_ROOT", HERE))
# supercell size: N x N x N conventional cells. The cell ids keep their 32-site
# names at every N (vac_Al32 is the perfect host, vac_V the bare vacancy).
N_REP = int(os.environ.get("VAC_N", "2"))
SUFFIX = "" if N_REP == 2 else f"_n{N_REP}"
CELLS = ROOT / f"candidates_vac{SUFFIX}"
OUT = ROOT / "results" / f"vacancy{SUFFIX}"
SOURCES = ("mace_relax", "lammps_meam_relax", "gpumd_nep89_relax", "qe_pbe")
LABEL = {"mace_relax": "MACE-MP-0", "lammps_meam_relax": "MEAM", "gpumd_nep89_relax": "NEP89",
         "qe_pbe": "QE PBE"}
ITEMS = {  # binding item: (complex, isolated parts)
    "Mg-V 1nn": ("vac_MgV_1nn", ["vac_Mg", "vac_V"]),
    "Mg-V 2nn": ("vac_MgV_2nn", ["vac_Mg", "vac_V"]),
    "Si-V 1nn": ("vac_SiV_1nn", ["vac_Si", "vac_V"]),
    "Si-V 2nn": ("vac_SiV_2nn", ["vac_Si", "vac_V"]),
    "MgSi-V": ("vac_MgSiV", ["vac_MgSi", "vac_V"]),
    "MgMg-V": ("vac_MgMgV", ["vac_MgMg", "vac_V"]),
    "SiSi-V": ("vac_SiSiV", ["vac_SiSi", "vac_V"]),
}


# ----------------------------------------------------------------- cells

def host(a: float = 4.04):
    return bulk("Al", "fcc", a=a, cubic=True).repeat((N_REP, N_REP, N_REP))


def neighbours(atoms, i):
    d = atoms.get_distances(i, range(len(atoms)), mic=True)
    a = atoms.cell[0, 0] / N_REP
    n1 = [j for j in range(len(atoms)) if abs(d[j] - a / np.sqrt(2)) < 0.05]
    n2 = [j for j in range(len(atoms)) if abs(d[j] - a) < 0.05]
    return n1, n2


def make_cells():
    h = host()
    n1, n2 = neighbours(h, 0)
    # a pair (p, q) of first neighbours of site 0 that are first neighbours of each other
    p, q = next((p, q) for p, q in combinations(n1, 2)
                if abs(h.get_distance(p, q, mic=True) - h.cell[0, 0] / N_REP / np.sqrt(2)) < 0.05)
    spec = {  # cid: (vacancy site or None, {site: element})
        "vac_Al32": (None, {}), "vac_V": (0, {}),
        "vac_Mg": (None, {0: "Mg"}), "vac_Si": (None, {0: "Si"}),
        "vac_MgV_1nn": (0, {n1[0]: "Mg"}), "vac_MgV_2nn": (0, {n2[0]: "Mg"}),
        "vac_SiV_1nn": (0, {n1[0]: "Si"}), "vac_SiV_2nn": (0, {n2[0]: "Si"}),
        "vac_MgSi": (None, {p: "Mg", q: "Si"}), "vac_MgMg": (None, {p: "Mg", q: "Mg"}),
        "vac_SiSi": (None, {p: "Si", q: "Si"}),
        "vac_MgSiV": (0, {p: "Mg", q: "Si"}), "vac_MgMgV": (0, {p: "Mg", q: "Mg"}),
        "vac_SiSiV": (0, {p: "Si", q: "Si"}),
    }
    CELLS.mkdir(exist_ok=True)
    for cid, (v, sub) in spec.items():
        a = h.copy()
        for i, el in sub.items():
            a[i].symbol = el
        if v is not None:
            del a[v]
        write(CELLS / f"{cid}.xyz", a)
    print(f"wrote {len(spec)} cells to {CELLS}")


def cell_at(cid, a_lat):
    at = read(CELLS / f"{cid}.xyz")
    at.set_cell(np.eye(3) * N_REP * a_lat, scale_atoms=True)
    return at


# ----------------------------------------------------------------- sources

def a_source(source):
    if source == "qe_pbe":
        return json.loads((ROOT / "results" / "qe_pbe_al_lattice.json").read_text())["a_A"]
    r = json.loads((ROOT / "results" / ("mace_relax" if source == "mace_relax" else source)
                    / "mp_134_Al.json").read_text())
    return (4 * r["volume_per_atom"]) ** (1 / 3)


def save(source, cid, rec):
    (OUT / source).mkdir(parents=True, exist_ok=True)
    (OUT / source / f"{cid}.json").write_text(json.dumps(rec, indent=1), encoding="utf-8")


def run_mace(cids, a):
    sys.path.insert(0, str(HERE))
    import dgx_sources as ds
    from ase.optimize import BFGS
    calc, _ = ds.get_mace("float64")
    for cid in cids:
        t0 = time.perf_counter()
        at = cell_at(cid, a)
        at.calc = calc
        ok = bool(BFGS(at, logfile=None).run(fmax=0.01, steps=500))
        save("mace_relax", cid, {"cid": cid, "status": "ok" if ok else "unconverged",
                                 "energy": float(at.get_potential_energy()), "natoms": len(at),
                                 "a_A": a, "wall_s": time.perf_counter() - t0, "device": "cuda",
                                 "settings": "MACE-MP-0 small float64, BFGS fmax 0.01 eV/A, fixed cell"})
        print("mace", cid, ok, flush=True)


def run_meam(cids, a):
    sys.path.insert(0, str(HERE))
    import dgx_sources as ds
    for cid in cids:
        t0 = time.perf_counter()
        wd = ROOT / f"work_vac{SUFFIX}" / "lammps" / cid
        wd.mkdir(parents=True, exist_ok=True)
        ds._write_lmp_data(cell_at(cid, a), wd / "data.lmp", velocities=False)
        lmp = ds._lmp_new()
        ds._lmp_setup(lmp, wd / "data.lmp")
        lmp.command("min_style cg")
        for _ in range(3):
            lmp.command("minimize 0 1e-8 20000 200000")
        e = lmp.get_thermo("pe")
        lmp.close()
        save("lammps_meam_relax", cid, {"cid": cid, "status": "ok", "energy": float(e), "a_A": a,
                                        "wall_s": time.perf_counter() - t0, "device": "cpu:1",
                                        "settings": "MEAM Jelinek 2012, cg minimize 0 1e-8 x3, fixed box"})
        print("meam", cid, e, flush=True)


def run_nep(cids, a):
    sys.path.insert(0, str(HERE))
    import dgx_sources as ds
    for cid in cids:
        t0 = time.perf_counter()
        wd = ROOT / f"work_vac{SUFFIX}" / "gpumd" / cid
        if wd.exists():
            import shutil
            shutil.rmtree(wd)
        wd.mkdir(parents=True)
        ds.write_gpumd_model(cell_at(cid, a), wd / "model.xyz", with_vel=False)
        (wd / "run.in").write_text(f"potential {ds.NEP89_FILE}\nminimize fire {ds.RELAX_FTOL:g} 50000\n"
                                   "ensemble nve\ntime_step 0\ndump_xyz 1 relaxed.xyz precision double force\n"
                                   "run 1\n", encoding="utf-8")
        ds.run_gpumd(wd)
        r = read(str(wd / "relaxed.xyz"), index=-1, format="extxyz")
        save("gpumd_nep89_relax", cid, {"cid": cid, "status": "ok", "energy": float(r.get_potential_energy()),
                                        "a_A": a, "wall_s": time.perf_counter() - t0, "device": "cuda",
                                        "settings": f"NEP89, GPUMD minimize fire {ds.RELAX_FTOL:g}, fixed box"})
        print("nep", cid, flush=True)


def run_qe(cids, a):
    sys.path.insert(0, str(HERE / "qe"))
    import qe_driver as qd
    work = Path(os.environ.get("QE_WORK", ROOT / "work_vac" / "qe"))
    if os.environ.get("VAC_REVERSE"):   # a second worker walks the list from the other end
        cids = cids[::-1]
    for cid in cids:
        if (OUT / "qe_pbe" / f"{cid}.json").exists() or (work / cid / "pbe").exists():
            continue
        at = cell_at(cid, a)
        rec, out, final = qd.pbe_relax(cid, at, fixed_cell=True, work=work)
        rec["source"] = "qe_pbe_vacancy"
        rec["energy"] = (rec["energy_per_atom"] * rec["natoms"]) if rec.get("energy_per_atom") is not None else None
        rec["a_A"] = a
        save("qe_pbe", cid, rec)
        print("qe", cid, rec["status"], f"{rec['wall_s']:.0f}s", flush=True)


# ----------------------------------------------------------------- analysis

def analyse():
    from scipy.stats import kendalltau
    res = {}
    for s in SOURCES:
        e = {}
        for p in (OUT / s).glob("*.json") if (OUT / s).exists() else []:
            r = json.loads(p.read_text())
            if r.get("status") == "ok" and r.get("energy") is not None:
                e[r["cid"]] = r["energy"]
        cost = sum(json.loads(p.read_text()).get("core_s") or json.loads(p.read_text()).get("wall_s") or 0
                   for p in (OUT / s).glob("*.json")) if (OUT / s).exists() else 0
        b = {}
        for item, (cx, parts) in ITEMS.items():
            need = [cx, *parts, "vac_Al32"]
            if all(c in e for c in need):
                b[item] = sum(e[c] for c in parts) - e[cx] - e["vac_Al32"]
        res[s] = {"binding_eV": b, "n_cells": len(e), "cost_s": cost}
    lines = [f"# Vacancy binding in fcc Al ({4 * N_REP ** 3}-site cells, positive = bound)", "",
             "| Item | " + " | ".join(LABEL[s] for s in SOURCES) + " |",
             "|---|" + "---|" * len(SOURCES)]
    for item in ITEMS:
        lines.append(f"| {item} | " + " | ".join(
            f"{res[s]['binding_eV'][item]:+.3f}" if item in res[s]["binding_eV"] else "n/a" for s in SOURCES) + " |")
    lines += ["", "Kendall tau of the binding ordering against each receiver (items both sources have):", ""]
    taus = {}
    for r_ in SOURCES:
        for d in SOURCES:
            if d == r_:
                continue
            common = [i for i in ITEMS if i in res[r_]["binding_eV"] and i in res[d]["binding_eV"]]
            if len(common) >= 4:
                t = kendalltau([res[d]["binding_eV"][i] for i in common],
                               [res[r_]["binding_eV"][i] for i in common]).statistic
                taus[f"{d}->{r_}"] = {"tau": float(t), "n": len(common)}
                lines.append(f"- {LABEL[d]} -> {LABEL[r_]}: tau = {t:+.2f} (n = {len(common)})")
    out = {"sources": res, "tau": taus}
    (ROOT / "results" / "analysis").mkdir(parents=True, exist_ok=True)
    (ROOT / "results" / "analysis" / f"vacancy_binding{SUFFIX}.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    (ROOT / "results" / "analysis" / f"vacancy_binding{SUFFIX}.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--make", action="store_true")
    ap.add_argument("--source", choices=SOURCES)
    ap.add_argument("--analyse", action="store_true")
    ap.add_argument("--cids", help="comma-separated subset of cells (one per CHPC array task)")
    args = ap.parse_args()
    if args.make:
        make_cells()
    if args.source:
        cids = sorted(p.stem for p in CELLS.glob("*.xyz"))
        if args.cids:
            cids = [c for c in cids if c in args.cids.split(",")]
        a = a_source(args.source)
        {"mace_relax": run_mace, "lammps_meam_relax": run_meam, "gpumd_nep89_relax": run_nep,
         "qe_pbe": run_qe}[args.source](cids, a)
    if args.analyse:
        analyse()


if __name__ == "__main__":
    main()
