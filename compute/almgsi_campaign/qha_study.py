#!/usr/bin/env python
"""Quasi-harmonic formation free energies of the Al-Mg-Si pool (MACE-MP-0).

The MD formation energies at 453 K and 803 K contain thermal expansion and
anharmonic relaxation but no vibrational entropy. This receiver adds it: for
each formation-energy candidate, the quasi-harmonic Helmholtz free energy

    F(V, T) = E_0(V) + F_vib(V, T)

is evaluated on five isotropic volumes (0.97 to 1.03 of the 0 K volume,
positions relaxed at each fixed cell), phonons from finite displacements with
phonopy in a supercell at least 9 A across, and F(T) = min over V of a cubic
fit at each temperature. The formation free energy per atom is taken against
fcc Al, hcp Mg and diamond Si treated the same way at the same T.

A candidate with an imaginary mode below -0.3 THz at the 0 K volume is
dynamically unstable in this model; it gets status "unstable" and no value
(no value is filled in).

Usage: python qha_study.py [--cids a,b] [--temps 300,453,803]
Output: results/qha/raw/<cid>.json (per candidate), results/qha/qha_T<T>.json
(live_values-style tables, eV/atom), results/qha/QHA.md
"""
from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.filters import FrechetCellFilter
from ase.io import read
from ase.optimize import BFGS

HERE = Path(__file__).resolve().parent
OUT = HERE / "results" / "qha"
STRAINS = (0.97, 0.985, 1.0, 1.015, 1.03)
TEMPS = (300, 453, 803)
REFS = {"Al": "mp_134_Al", "Mg": "mp_153_Mg", "Si": "mp_149_Si"}
MESH = 20.0          # phonopy automatic q mesh length (A)
IMAG_THz = -0.3


def to_phonopy(a: Atoms):
    from phonopy.structure.atoms import PhonopyAtoms
    return PhonopyAtoms(symbols=a.get_chemical_symbols(), cell=a.cell.array, scaled_positions=a.get_scaled_positions())


def sc_matrix(a: Atoms, length: float = 9.0):
    lens = a.cell.lengths()
    return np.diag([max(1, math.ceil(length / L)) for L in lens])


def relax(a: Atoms, calc, cell: bool):
    a = a.copy()
    a.calc = calc
    obj = FrechetCellFilter(a) if cell else a
    BFGS(obj, logfile=None).run(fmax=0.005, steps=800)
    return a


def phonon_free_energy(a: Atoms, calc, temps):
    from phonopy import Phonopy
    ph = Phonopy(to_phonopy(a), supercell_matrix=sc_matrix(a), primitive_matrix=None)
    ph.generate_displacements(distance=0.01)
    forces = []
    for sc in ph.supercells_with_displacements:
        s = Atoms(sc.symbols, cell=sc.cell, scaled_positions=sc.scaled_positions, pbc=True)
        s.calc = calc
        forces.append(s.get_forces())
    ph.forces = np.array(forces)
    ph.produce_force_constants()
    ph.run_mesh(MESH)
    fmin = float(np.min(ph.mesh.frequencies))
    ph.run_thermal_properties(temperatures=list(temps))
    tp = ph.get_thermal_properties_dict()
    # phonopy returns kJ/mol per unit cell; convert to eV per cell
    fe = np.array(tp["free_energy"]) / 96.485332
    return fe, fmin, len(ph.supercells_with_displacements)


def run_candidate(cid: str, calc, temps) -> dict:
    t0 = time.perf_counter()
    a0 = relax(read(HERE / "candidates" / f"{cid}.xyz"), calc, cell=True)
    n = len(a0)
    rec = {"cid": cid, "natoms": n, "formula": a0.get_chemical_formula(), "temps": list(temps)}
    vols, E, F = [], [], []
    for s in STRAINS:
        a = a0.copy()
        a.set_cell(a0.cell.array * s ** (1 / 3), scale_atoms=True)
        a = relax(a, calc, cell=False)
        fe, fmin, ndisp = phonon_free_energy(a, calc, temps)
        if s == 1.0:
            rec["min_freq_THz"] = fmin
            rec["n_displacements"] = ndisp
        vols.append(a.get_volume() / n)
        E.append(a.get_potential_energy() / n)
        F.append(fe / n)
    F = np.array(F)
    rec["volumes"], rec["E0"] = vols, E
    if rec["min_freq_THz"] < IMAG_THz:
        rec["status"] = "unstable"
    else:
        g = []
        for j, T in enumerate(temps):
            c = np.polyfit(vols, np.array(E) + F[:, j], 3)
            vv = np.linspace(min(vols), max(vols), 401)
            g.append(float(np.polyval(c, vv).min()))
        rec["F_eV_atom"] = g
        rec["status"] = "ok"
    rec["wall_s"] = time.perf_counter() - t0
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cids", default="")
    ap.add_argument("--temps", default=",".join(map(str, TEMPS)))
    args = ap.parse_args()
    temps = [int(t) for t in args.temps.split(",")]
    import dgx_sources as ds
    calc, _ = ds.get_mace("float64")
    idx = json.loads((HERE / "candidates_index.json").read_text())
    cids = args.cids.split(",") if args.cids else \
        list(REFS.values()) + [r["cid"] for r in idx if r["property"] == "formation_energy" and r["cid"] not in REFS.values()]
    (OUT / "raw").mkdir(parents=True, exist_ok=True)
    for cid in cids:
        p = OUT / "raw" / f"{cid}.json"
        if p.exists():
            continue
        try:
            rec = run_candidate(cid, calc, temps)
        except Exception as exc:  # noqa: BLE001
            rec = {"cid": cid, "status": "failed", "error": f"{type(exc).__name__}: {exc}"}
        p.write_text(json.dumps(rec, indent=1), encoding="utf-8")
        print(cid, rec["status"], f"{rec.get('wall_s', 0):.1f}s", rec.get("min_freq_THz"), flush=True)
    tabulate(temps)


def tabulate(temps):
    raw = {p.stem: json.loads(p.read_text()) for p in (OUT / "raw").glob("*.json")}
    ref = {el: raw.get(c) for el, c in REFS.items()}
    lines = ["# Quasi-harmonic formation free energies (MACE-MP-0)", "",
             f"Candidates computed: {len(raw)}; ok: {sum(r['status'] == 'ok' for r in raw.values())}; "
             f"dynamically unstable: {sum(r['status'] == 'unstable' for r in raw.values())}; "
             f"failed: {sum(r['status'] == 'failed' for r in raw.values())}", ""]
    if not all(r and r["status"] == "ok" for r in ref.values()):
        lines.append("Elemental references incomplete; no formation free energies.")
    for j, T in enumerate(temps):
        tab = {}
        for cid, r in raw.items():
            if r.get("status") != "ok" or not all(x and x["status"] == "ok" for x in ref.values()):
                continue
            from collections import Counter
            cnt = Counter(read(HERE / "candidates" / f"{cid}.xyz").get_chemical_symbols())
            n = sum(cnt.values())
            g = r["F_eV_atom"][j] - sum(cnt[el] / n * ref[el]["F_eV_atom"][j] for el in cnt)
            tab[cid] = {"property": "formation_energy", "value": g, "sem": 0.0, "wall_s": r["wall_s"],
                        "status": "ok", "device": "cuda"}
        (OUT / f"qha_T{T}.json").write_text(json.dumps(tab, indent=1), encoding="utf-8")
        lines.append(f"- T = {T} K: {len(tab)} formation free energies")
    (OUT / "QHA.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
