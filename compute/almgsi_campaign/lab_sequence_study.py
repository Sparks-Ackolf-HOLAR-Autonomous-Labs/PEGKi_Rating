#!/usr/bin/env python
"""Laboratory receiver for the computed Al-Mg-Si campaign: the observed
precipitation sequence of 6xxx alloys.

Experiment (electron microscopy, atom-probe tomography, calorimetry) observes
that a supersaturated Al-Mg-Si solid solution decomposes in a fixed order
(Edwards et al., Acta Mater. 46 (1998) 3893):

  stage 0  isolated Mg and Si solutes (supersaturated solid solution)
  stage 1  separate clusters of Mg and of Si (Mg-Mg, Si-Si pairs)
  stage 2  Mg-Si co-clusters of 2 to 4 atoms
  stage 3  beta''  (Mg5Si6, Zandbergen et al. 1997)
  stage 4  beta' (Mg9Si5, Vissers et al. 2007), U1 (MgAl2Si2), U2 (MgAlSi)
           (Andersen et al. 2007)
  stage 5  beta    (Mg2Si, the equilibrium phase)

A source that describes the alloy correctly should place each later stage
lower in energy per solute atom (Mg + Si) than the earlier ones, since each
transformation in the sequence releases energy. This is an ordering measured
in the laboratory; it carries no energy value, so it is a receiver of exactly
the kind PEGKi is built for.

Energy per solute atom, relative to fcc Al, hcp Mg and diamond Si computed by
the same source at the same temperature:
  phases:   E_sol = E_f(per atom) * N / N_solute
  clusters: E_sol = [E(cell) - N_Al * E(Al108)/108 - sum_X N_X mu_X] / N_solute
            (cluster cells and the Al108 host share one lattice per source)

Metrics per source: agreement A = -Kendall tau-b between E_sol and stage
(+1 = the sequence exactly, 0 = unrelated, -1 = reversed), with a bootstrap
over items (500 resamples); the number of cross-stage pairs in the wrong
order; whether beta is the lowest; and the stage means.
The Al-Mg edge (gamma-Mg17Al12, epsilon-Mg23Al30) is reported separately:
both are observed phases (Murray 1982) and should lie below the elements.

Output: results/analysis/lab_sequence_study.json and .md
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from ase.io import read
from scipy.stats import kendalltau, spearmanr

HERE = Path(__file__).resolve().parent
RES = HERE / "results"
OUT = RES / "analysis"
REFS = {"Al": "mp_134_Al", "Mg": "mp_153_Mg", "Si": "mp_149_Si"}
HOST = "clu_Al108"
STAGE = {"clu_Mg1": 0, "clu_Si1": 0,
         "clu_MgMg_1nn": 1, "clu_SiSi_1nn": 1, "clu_SiSi_2nn": 1,
         "clu_MgSi_1nn": 2, "clu_MgSi_2nn": 2, "clu_Mg2Si_tri": 2, "clu_MgSi2_tri": 2,
         "clu_Mg2Si2_tet": 2, "clu_Mg3Si_tet": 2, "clu_MgSi3_tet": 2,
         "mp_568306_Mg5Si6": 3,
         "mp_542854_Mg9Si5": 4, "mp_3833_Mg(AlSi)2": 4, "mp_4499_MgAlSi": 4,
         "mp_1367_Mg2Si": 5}
NSTAGE = 6
PHASE = {"mp_568306_Mg5Si6": "beta''", "mp_542854_Mg9Si5": "beta'", "mp_3833_Mg(AlSi)2": "U1",
         "mp_4499_MgAlSi": "U2", "mp_1367_Mg2Si": "beta", "mp_2151_Mg17Al12": "gamma",
         "mp_17659_Mg23Al30": "epsilon"}
AL_MG_EDGE = ("mp_2151_Mg17Al12", "mp_17659_Mg23Al30")
SOURCES = {"qe_pbe": ("qe_pbe", "qe_pbe_dgx"), "qe_hse06": ("qe_hse06", "qe_hse06_dgx"),
           "mace_relax": ("mace_relax",), "mace_md": ("mace_md",),
           "gpumd_nep89_md": ("gpumd_nep89_md",), "lammps_meam_md": ("lammps_meam_md",),
           "gpumd_nep89_relax": ("gpumd_nep89_relax",), "lammps_meam_relax": ("lammps_meam_relax",)}
LABEL = {"qe_pbe": "QE PBE, 0 K", "qe_hse06": "QE HSE06, 0 K", "mace_relax": "MACE-MP-0, 0 K",
         "mace_md": "MACE-MP-0 MD, 300 K", "gpumd_nep89_md": "NEP89 MD, 300 K",
         "lammps_meam_md": "MEAM MD, 300 K", "gpumd_nep89_relax": "NEP89, 0 K",
         "lammps_meam_relax": "MEAM, 0 K", "archive_mp": "MP archive"}


def counts(cid):
    a = read(str(HERE / "candidates" / f"{cid}.xyz"), format="extxyz")
    c = {}
    for s in a.get_chemical_symbols():
        c[s] = c.get(s, 0) + 1
    return c


def load_records(dirs):
    recs = {}
    for d in dirs:
        p = RES / d
        if not p.exists():
            continue
        for f in sorted(p.glob("*.json")):
            if f.name.startswith("_") or "." in f.stem:
                continue
            r = json.loads(f.read_text(encoding="utf-8"))
            if isinstance(r, dict) and r.get("status") == "ok" and "energy_per_atom" in r:
                recs.setdefault(r["cid"], r)
    return recs


def check_cost(recs):
    """Compute charged to the staged items: wall seconds and core-seconds (device cpu:N)."""
    wall = core = 0.0
    for cid in STAGE:
        r = recs.get(cid)
        if r is None:
            continue
        w = float(r.get("wall_s") or 0.0)
        dev = str(r.get("device") or "")
        n = int(dev.split(":")[1]) if dev.startswith("cpu:") and dev.split(":")[1].isdigit() else 1
        wall += w
        core += w * n
    return {"wall_s": wall, "core_s": core}


def e_solute(recs):
    """Energy per solute atom (eV) for every staged item and the Al-Mg edge."""
    if not all(c in recs for c in REFS.values()):
        return {}
    mu = {el: recs[c]["energy_per_atom"] for el, c in REFS.items()}
    out = {}
    for cid in list(STAGE) + list(AL_MG_EDGE):
        if cid not in recs:
            continue
        n = counts(cid)
        N = sum(n.values())
        nsol = n.get("Mg", 0) + n.get("Si", 0)
        E = recs[cid]["energy_per_atom"] * N
        if cid.startswith("clu_"):
            if HOST not in recs:
                continue
            mu_al = recs[HOST]["energy_per_atom"]
            val = (E - n.get("Al", 0) * mu_al - n.get("Mg", 0) * mu["Mg"] - n.get("Si", 0) * mu["Si"]) / nsol
        else:
            val = (E - sum(n[el] * mu[el] for el in n)) / (nsol if cid in STAGE else N)
        out[cid] = float(val)
    return out


def mp_archive():
    idx = {d["cid"]: d for d in json.loads((HERE / "candidates_index_observed.json").read_text())}
    out = {}
    for cid in list(STAGE) + list(AL_MG_EDGE):
        d = idx.get(cid)
        if d is None or d.get("mp_ef") is None:
            continue
        n = counts(cid)
        N = sum(n.values())
        nsol = n.get("Mg", 0) + n.get("Si", 0)
        out[cid] = float(d["mp_ef"] * N / (nsol if cid in STAGE else N))
    return out


def agreement(es, n_boot=500, seed=7):
    items = [c for c in STAGE if c in es]
    if len(items) < 4 or len({STAGE[c] for c in items}) < 2:
        return None
    x = np.array([es[c] for c in items])
    y = np.array([STAGE[c] for c in items])
    a = -kendalltau(x, y).statistic
    rng = np.random.default_rng(seed)
    bs = []
    for _ in range(n_boot):
        i = rng.integers(0, len(items), len(items))
        if len(set(y[i])) < 2:
            continue
        t = kendalltau(x[i], y[i]).statistic
        if math.isfinite(t):
            bs.append(-t)
    wrong = tot = 0
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            if y[i] == y[j]:
                continue
            tot += 1
            later, earlier = (i, j) if y[i] > y[j] else (j, i)
            wrong += x[later] > x[earlier]
    stage_mean = {int(s): float(np.mean([es[c] for c in items if STAGE[c] == s]))
                  for s in sorted(set(y.tolist()))}
    phases = [c for c in items if not c.startswith("clu_")]
    beta_lowest = ("mp_1367_Mg2Si" in es and len(phases) > 1
                   and min(phases, key=lambda c: es[c]) == "mp_1367_Mg2Si")
    return {"n_items": len(items), "agreement": float(a),
            "ci90": [float(np.percentile(bs, 5)), float(np.percentile(bs, 95))] if bs else None,
            "wrong_pairs": int(wrong), "cross_stage_pairs": int(tot),
            "beta_lowest_phase": bool(beta_lowest), "stage_mean_eV_per_solute": stage_mean}


def main():
    res = {"sources": {}, "values": {}, "cost": {}}
    for s, dirs in SOURCES.items():
        recs = load_records(dirs)
        es = e_solute(recs)
        if es:
            res["values"][s] = es
            res["cost"][s] = check_cost(recs)
    res["values"]["archive_mp"] = mp_archive()
    for s, es in res["values"].items():
        r = agreement(es) or {}
        r["label"] = LABEL[s]
        r["check_cost"] = res["cost"].get(s)
        r["al_mg_edge_eV_per_atom"] = {PHASE[c]: es[c] for c in AL_MG_EDGE if c in es}
        r["phases_eV_per_solute"] = {PHASE[c]: es[c] for c in PHASE if c in es and c in STAGE}
        res["sources"][s] = r
    # Does agreement with PBE (a calculation) predict agreement with the laboratory sequence?
    sa = json.loads((OUT / "source_agreement.json").read_text()) if (OUT / "source_agreement.json").exists() else {}
    res["note"] = ("agreement = -Kendall tau-b(E per solute, stage); +1 reproduces the observed "
                   "sequence. Stage 0-2 items are cells in a 108-site Al host; stages 3-5 are bulk phases.")
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "lab_sequence_study.json").write_text(json.dumps(res, indent=1))
    lines = ["# Observed precipitation sequence as a receiver (Al-Mg-Si)", "",
             "| Source | items | agreement with the sequence (90% CI) | wrong-order pairs | beta lowest phase | "
             "stage means, eV per solute (0 / 1 / 2 / 3 / 4 / 5) |", "|---|---|---|---|---|---|"]
    for s, r in res["sources"].items():
        if "agreement" not in r:
            continue
        sm = " / ".join(f"{r['stage_mean_eV_per_solute'].get(k, float('nan')):+.3f}" for k in range(NSTAGE))
        ci = r["ci90"]
        lines.append(f"| {r['label']} | {r['n_items']} | {r['agreement']:+.2f} "
                     f"[{ci[0]:+.2f}, {ci[1]:+.2f}] | {r['wrong_pairs']}/{r['cross_stage_pairs']} | "
                     f"{'yes' if r['beta_lowest_phase'] else 'no'} | {sm} |")
    lines += ["", "## Phases, eV per solute atom", "",
              "| Source | " + " | ".join(["beta''", "beta'", "U1", "U2", "beta"]) + " | gamma (eV/atom) | epsilon (eV/atom) |",
              "|---" * 8 + "|"]
    for s, r in res["sources"].items():
        p = r["phases_eV_per_solute"]
        e = r["al_mg_edge_eV_per_atom"]
        cells = [f"{p[k]:+.3f}" if k in p else "n/a" for k in ["beta''", "beta'", "U1", "U2", "beta"]]
        cells += [f"{e[k]:+.3f}" if k in e else "n/a" for k in ["gamma", "epsilon"]]
        lines.append(f"| {r['label']} | " + " | ".join(cells) + " |")
    (OUT / "lab_sequence_study.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
