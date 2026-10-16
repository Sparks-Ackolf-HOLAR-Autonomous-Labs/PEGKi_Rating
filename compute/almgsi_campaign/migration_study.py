#!/usr/bin/env python
"""Vacancy migration of Si, Mg and Al in fcc Al: how fast Si moves through a contact.

Companion to interface_study.py (does Si leave the Si substrate into an
Al-Mg-Si contact?). Thermodynamics says whether Si dissolves; diffusion
says how far it travels during the contact anneal. Solutes in Al move by
exchanging with a nearest-neighbour vacancy, so the activation energy
for solute X is

    Q_X = E_f(V) + E_m(X) - E_b(X-V 1nn)

with E_f(V) the vacancy formation energy, E_b the vacancy-solute binding
(both from vacancy_study.py, 108-site cells) and E_m the migration
barrier computed here. E_m is the energy of the symmetric saddle: the
jumping atom is held at the midpoint between its site and the vacant
site, all other atoms relax at fixed cell; the initial state is the
relaxed X-V 1nn cell of vacancy_study.py (Al: the bare vacancy cell).
For a direct exchange in fcc the midpoint lies on the mirror plane
between the two sites, so the constrained midpoint is the saddle of the
symmetric path.

Usage:
  python migration_study.py --make                 write mig_{Al,Si,Mg} cells (108 sites)
  python migration_study.py --source mace_relax|lammps_meam_relax
  python migration_study.py --source qe_pbe        (DFT: CHPC guest, qe/migration_qe.sbatch)
  python migration_study.py --analyse              results/analysis/migration.{json,md}
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("VAC_N", "3")
import numpy as np  # noqa: E402
from ase.constraints import FixAtoms  # noqa: E402
from ase.io import read, write  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import vacancy_study as vs  # noqa: E402

N_REP = vs.N_REP
MIG = vs.ROOT / f"candidates_mig_n{N_REP}"
OUT = vs.ROOT / "results" / f"migration_n{N_REP}"
INIT = {"mig_Al": "vac_V", "mig_Si": "vac_SiV_1nn", "mig_Mg": "vac_MgV_1nn"}
KB = 8.617333262e-5


def make():
    h = vs.host()
    n1, _ = vs.neighbours(h, 0)
    j = n1[0]
    MIG.mkdir(exist_ok=True)
    meta = {}
    for cid, el in (("mig_Al", "Al"), ("mig_Si", "Si"), ("mig_Mg", "Mg")):
        a = h.copy()
        a[j].symbol = el
        d = a.get_distance(j, 0, mic=True, vector=True)  # from j to the vacant site 0
        a.positions[j] = a.positions[j] + 0.5 * d
        del a[0]
        mover = j - 1  # index shift after deleting site 0
        write(MIG / f"{cid}.xyz", a)
        meta[cid] = {"mover": mover, "element": el}
    (MIG / "meta.json").write_text(json.dumps(meta, indent=1), encoding="utf-8")
    print("wrote", list(meta), "to", MIG)


def cell_at(cid, a_lat):
    at = read(MIG / f"{cid}.xyz")
    at.set_cell(np.eye(3) * N_REP * a_lat, scale_atoms=True)
    m = json.loads((MIG / "meta.json").read_text())[cid]["mover"]
    at.set_constraint(FixAtoms(indices=[m]))
    return at


def run(source, only=None):
    vs.CELLS = MIG
    vs.OUT = OUT
    vs.cell_at = cell_at
    cids = [c for c in INIT if not only or c == only]
    a = vs.a_source(source)
    if source == "mace_relax":
        vs.run_mace(cids, a)
    elif source == "qe_pbe":
        vs.run_qe(cids, a)
    elif source == "lammps_meam_relax":
        import dgx_sources as ds
        meta = json.loads((MIG / "meta.json").read_text())
        for cid in cids:
            t0 = time.perf_counter()
            wd = vs.ROOT / f"work_mig_n{N_REP}" / "lammps" / cid
            wd.mkdir(parents=True, exist_ok=True)
            ds._write_lmp_data(cell_at(cid, a), wd / "data.lmp", velocities=False)
            lmp = ds._lmp_new()
            ds._lmp_setup(lmp, wd / "data.lmp")
            lmp.command(f"group mover id {meta[cid]['mover'] + 1}")
            lmp.command("fix hold mover setforce 0.0 0.0 0.0")
            lmp.command("min_style cg")
            for _ in range(3):
                lmp.command("minimize 0 1e-8 20000 200000")
            e = lmp.get_thermo("pe")
            lmp.close()
            vs.save("lammps_meam_relax", cid, {"cid": cid, "status": "ok", "energy": float(e), "a_A": a,
                                               "wall_s": time.perf_counter() - t0, "device": "cpu:1",
                                               "settings": "MEAM Jelinek 2012, mover fixed, cg x3, fixed box"})
            print("meam", cid, e, flush=True)
    else:
        raise SystemExit(f"unknown source {source}")


def energy(path):
    if not path.exists():
        return None
    r = json.loads(path.read_text())
    if r.get("status") != "ok":
        return None
    if r.get("energy") is not None:
        return r["energy"]
    return r["energy_per_atom"] * r["natoms"] if r.get("energy_per_atom") is not None else None


def analyse():
    vac = json.loads((vs.ROOT / "results" / "analysis" / f"vacancy_binding{vs.SUFFIX}.json").read_text()) \
        if (vs.ROOT / "results" / "analysis" / f"vacancy_binding{vs.SUFFIX}.json").exists() else None
    rows = {}
    for src in ("qe_pbe", "mace_relax", "lammps_meam_relax"):
        vdir = vs.ROOT / "results" / f"vacancy{vs.SUFFIX}" / src
        e_host, e_v = energy(vdir / "vac_Al32.json"), energy(vdir / "vac_V.json")
        if e_host is None or e_v is None:
            continue
        n = 4 * N_REP ** 3
        ef_v = e_v - (n - 1) / n * e_host
        r = {"E_f_V": ef_v}
        for cid, init in INIT.items():
            el = cid.split("_")[1]
            e_mid, e_ini = energy(OUT / src / f"{cid}.json"), energy(vdir / f"{init}.json")
            if e_mid is None or e_ini is None:
                continue
            em = e_mid - e_ini
            r[f"E_m_{el}"] = em
            if el == "Al":
                r["Q_Al"] = ef_v + em
            else:
                e_x = energy(vdir / f"vac_{el}.json")
                eb = e_x + e_v - e_ini - e_host  # positive = bound, as in vacancy_study.py
                r[f"E_b_{el}V"] = eb
                r[f"Q_{el}"] = ef_v + em - eb
        rows[src] = r
    for src, r in rows.items():
        if "Q_Si" in r and "Q_Al" in r:
            r["D_Si_over_D_Al_723K"] = math.exp(-(r["Q_Si"] - r["Q_Al"]) / (KB * 723.0))
    (vs.ROOT / "results" / "analysis").mkdir(exist_ok=True)
    (vs.ROOT / "results" / "analysis" / "migration.json").write_text(json.dumps(rows, indent=1))
    keys = ["E_f_V", "E_m_Al", "E_m_Si", "E_m_Mg", "E_b_SiV", "E_b_MgV", "Q_Al", "Q_Si", "Q_Mg",
            "D_Si_over_D_Al_723K"]
    L = ["# Vacancy migration in fcc Al (108-site cells, eV)", "",
         "| source | " + " | ".join(keys) + " |", "|---" * (len(keys) + 1) + "|"]
    for src, r in rows.items():
        L.append(f"| {src} | " + " | ".join("" if r.get(k) is None else f"{r[k]:+.3f}" for k in keys) + " |")
    (vs.ROOT / "results" / "analysis" / "migration.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--make", action="store_true")
    ap.add_argument("--source")
    ap.add_argument("--cid")
    ap.add_argument("--analyse", action="store_true")
    a = ap.parse_args()
    if a.make:
        make()
    if a.source:
        run(a.source, a.cid)
    if a.analyse:
        analyse()


if __name__ == "__main__":
    main()
