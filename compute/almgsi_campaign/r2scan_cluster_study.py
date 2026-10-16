#!/usr/bin/env python
"""Solute-cluster binding at r2SCAN against PBE (Al-Mg-Si).

GPAW r2SCAN single points (PAW, 700 eV, libxc r2SCAN) on the QE PBE relaxed
108-site cluster cells give a second DFT binding energy for the nine Mg/Si
clusters. Reported: per-cluster binding at both levels, Kendall tau-b between
them (bootstrap 90% CI over clusters), sign agreement, and each potential's
tau against r2SCAN next to its tau against PBE.

Caveat recorded in the output: r2SCAN values are at the PBE geometry, so
relaxation at r2SCAN is not included.

Output: results/analysis/r2scan_cluster.json and .md
"""
from __future__ import annotations

import json

import numpy as np
from scipy.stats import kendalltau

import cluster_transfer_study as ct

R2 = "gpaw_r2scan"


def tau_ci(a, b, rng, n_boot=2000):
    a, b = np.asarray(a), np.asarray(b)
    t = kendalltau(a, b)[0]
    bs = []
    for _ in range(n_boot):
        i = rng.integers(0, len(a), len(a))
        if len(set(a[i])) > 1 and len(set(b[i])) > 1:
            bs.append(kendalltau(a[i], b[i])[0])
    return float(t), [float(np.percentile(bs, 5)), float(np.percentile(bs, 95))]


def main(seed=11):
    rng = np.random.default_rng(seed)
    pbe, r2 = ct.bindings(ct.RECEIVER), ct.bindings(R2)
    com = sorted(set(pbe) & set(r2))
    a = [pbe[k]["E_b_eV"] for k in com]
    b = [r2[k]["E_b_eV"] for k in com]
    t, ci = tau_ci(a, b, rng)
    res = {"clusters": {k: {"pbe": pbe[k]["E_b_eV"], "r2scan": r2[k]["E_b_eV"]} for k in com},
           "pbe_vs_r2scan": {"n": len(com), "kendall_tau": t, "ci90": ci,
                             "sign_agree": int(sum((x > 0) == (y > 0) for x, y in zip(a, b))),
                             "mad_eV": float(np.mean(np.abs(np.subtract(a, b))))},
           "sources": {}, "geometry": "QE PBE relaxed; no r2SCAN relaxation"}
    for s, lab in ct.SOURCES.items():
        x = ct.bindings(s)
        row = {"label": lab}
        for name, ref in (("vs_pbe", pbe), ("vs_r2scan", r2)):
            c = sorted(set(x) & set(ref))
            if len(c) >= 4:
                tt, cc = tau_ci([x[k]["E_b_eV"] for k in c], [ref[k]["E_b_eV"] for k in c], rng)
                row[name] = {"n": len(c), "kendall_tau": tt, "ci90": cc}
        res["sources"][s] = row
    (ct.OUT / "r2scan_cluster.json").write_text(json.dumps(res, indent=2))
    p = res["pbe_vs_r2scan"]
    lines = ["# Solute-cluster binding: r2SCAN against PBE", "",
             "GPAW r2SCAN single points at the QE PBE geometry (108-site Al cells). "
             "Binding in eV, positive = attractive.", "",
             "| Cluster | PBE | r2SCAN |", "|---|---|---|"]
    lines += [f"| {k[4:]} | {v['pbe']:+.3f} | {v['r2scan']:+.3f} |" for k, v in res["clusters"].items()]
    lines += ["", f"PBE vs r2SCAN: n = {p['n']}, tau {p['kendall_tau']:+.2f} "
              f"[{p['ci90'][0]:+.2f}, {p['ci90'][1]:+.2f}], sign agreement {p['sign_agree']}/{p['n']}, "
              f"mean abs. difference {p['mad_eV']:.3f} eV.", "",
              "| Source | tau vs PBE | tau vs r2SCAN |", "|---|---|---|"]
    for s, row in res["sources"].items():
        f = lambda d: f"{d['kendall_tau']:+.2f} [{d['ci90'][0]:+.2f}, {d['ci90'][1]:+.2f}]" if d else "n/a"
        lines.append(f"| {row['label']} | {f(row.get('vs_pbe'))} | {f(row.get('vs_r2scan'))} |")
    (ct.OUT / "r2scan_cluster.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
