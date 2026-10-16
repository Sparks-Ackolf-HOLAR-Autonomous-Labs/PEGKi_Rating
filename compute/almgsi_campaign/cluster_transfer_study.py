#!/usr/bin/env python
"""Solute-cluster binding energies as a transfer question (Al-Mg-Si).

Early clustering in a 6xxx alloy is governed by whether Mg and Si solutes
attract each other in the Al matrix. Each computed source gives a binding
energy for the same ten solute clusters in a 108-site Al cell
(aggregate_sources.binding; positive = attractive). PBE is the receiver.

Per source, against PBE on the clusters both hold:
  Kendall tau-b of the binding energies (bootstrap 90% CI over clusters),
  sign agreement (attractive versus repulsive),
  mean absolute difference (eV),
and, for contrast, the same source's Kendall tau against PBE on the ordered
compounds (results/analysis/source_agreement.json), so a source that agrees on
compounds but not on clusters is visible.

Output: results/analysis/cluster_transfer.json and .md
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import kendalltau

import aggregate_sources as ag

HERE = Path(__file__).resolve().parent
RES = HERE / "results"
OUT = RES / "analysis"
SOURCES = {"mace_relax": "MACE-MP-0, 0 K", "mace_md": "MACE-MP-0 MD, 300 K",
           "lammps_meam_md": "MEAM MD, 300 K", "gpumd_nep89_md": "NEP89 MD, 300 K"}
RECEIVER = "qe_pbe"


def load(src):
    recs = {}
    for p in sorted((RES / src).glob("*.json")):
        if p.name.startswith("_") or "." in p.stem:
            continue
        r = json.loads(p.read_text(encoding="utf-8"))
        if isinstance(r, dict) and "cid" in r and "energy_per_atom" in r:
            recs[r["cid"]] = r
    return recs


def bindings(src):
    recs = load(src)
    out = {}
    for cid, rec in recs.items():
        if not cid.startswith("clu_") or cid in (ag.HOST, "clu_Mg1", "clu_Si1"):
            continue
        comp, _ = ag.composition(cid)
        v, s, st = ag.binding(rec, recs, comp)
        if st == "ok":
            out[cid] = {"E_b_eV": v, "sem": s}
    return out


def compound_tau():
    p = OUT / "source_agreement.json"
    if not p.exists():
        return {}
    d = json.loads(p.read_text()).get(RECEIVER, {})
    return {s: v.get("kendall_tau") for s, v in d.items() if isinstance(v, dict)}


def main(n_boot=2000, seed=11):
    rec = bindings(RECEIVER)
    ctau = compound_tau()
    res = {"receiver": RECEIVER, "values": {RECEIVER: rec}, "sources": {}}
    rng = np.random.default_rng(seed)
    for s, lab in SOURCES.items():
        b = bindings(s)
        res["values"][s] = b
        common = sorted(set(b) & set(rec))
        if len(common) < 4:
            continue
        x = np.array([b[c]["E_b_eV"] for c in common])
        y = np.array([rec[c]["E_b_eV"] for c in common])
        tau = kendalltau(x, y).statistic
        bs = []
        for _ in range(n_boot):
            i = rng.integers(0, len(common), len(common))
            if len(set(i.tolist())) < 3:
                continue
            t = kendalltau(x[i], y[i]).statistic
            if math.isfinite(t):
                bs.append(t)
        res["sources"][s] = {
            "label": lab, "n": len(common), "kendall_tau": float(tau),
            "ci90": [float(np.percentile(bs, 5)), float(np.percentile(bs, 95))],
            "sign_agreement": float(np.mean(np.sign(x) == np.sign(y))),
            "mad_eV": float(np.mean(np.abs(x - y))),
            "n_attractive": int(np.sum(x > 0)),
            "compound_tau_vs_pbe": ctau.get(s),
        }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "cluster_transfer.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    lines = ["# Solute-cluster binding energies against PBE (Al-Mg-Si)", "",
             f"PBE clusters: {len(rec)}; attractive at PBE: {sum(v['E_b_eV'] > 0 for v in rec.values())}", "",
             "| Source | n | Kendall tau vs PBE (90% CI) | sign agreement | MAD (eV) | attractive | tau vs PBE on compounds |",
             "|---|---|---|---|---|---|---|"]
    for s, r in res["sources"].items():
        ct = r["compound_tau_vs_pbe"]
        lines.append(f"| {r['label']} | {r['n']} | {r['kendall_tau']:+.2f} [{r['ci90'][0]:+.2f}, {r['ci90'][1]:+.2f}] | "
                     f"{r['sign_agreement']:.2f} | {r['mad_eV']:.3f} | {r['n_attractive']}/{r['n']} | "
                     f"{'n/a' if ct is None else f'{ct:+.2f}'} |")
    lines += ["", "## Binding energies (eV, positive = attractive)", "",
              "| Cluster | " + " | ".join(["PBE"] + [SOURCES[s] for s in SOURCES]) + " |",
              "|---" * (len(SOURCES) + 2) + "|"]
    allc = sorted(set().union(*[set(v) for v in res["values"].values()]))
    for c in allc:
        cells = []
        for s in [RECEIVER] + list(SOURCES):
            v = res["values"][s].get(c)
            cells.append("n/a" if v is None else f"{v['E_b_eV']:+.3f}")
        lines.append(f"| {c[4:]} | " + " | ".join(cells) + " |")
    (OUT / "cluster_transfer.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
