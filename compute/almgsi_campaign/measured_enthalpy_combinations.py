#!/usr/bin/env python
"""Measured formation-enthalpy receiver over every combination of measurements.

Three compounds have measured formation enthalpies in the sources read
(results/lab_receivers/LITERATURE_CANDIDATES.md, section A1), each from two
independent measurements:

  beta  (Al140Mg89, "Al3Mg2"): Brown and Pratt 1970 (-2.47 kJ/mol-atoms) and
        Predel and Huelse 1978 (-1.22), both as tabulated in Mezbahul-Islam
        et al. 2014, Table 4 (doi:10.1155/2014/704283);
  gamma (Al12Mg17):            Brown and Pratt 1970 (-3.80) and Predel and
        Huelse 1978 (-1.01), same table;
  Mg2Si:                       NIST-JANAF Mg-034 (-25.94 kJ/mol-atoms) and the
        solution-calorimetry value of Feufel et al. 1997 as quoted by Jund et
        al. 2013 (about -22.2 kJ/mol-atoms).

That gives 2 x 2 x 2 = 8 measured receivers. For each receiver and each
computed source the script reports the signed error per compound, the mean
absolute error, Kendall tau over the three compounds, whether the beta-gamma
order agrees, and the binding scale (least-squares slope of computed against
measured through the origin: the fraction of the measured binding the source
reproduces; the origin is fixed because both scales are relative to the pure
elements).

Mapping of compounds to computed values:
  gamma -> mp_2151_Mg17Al12 (the gamma structure itself);
  Mg2Si -> mp_1367_Mg2Si (antifluorite, as in lab_receivers_study.py);
  beta  -> the 1168-atom beta cell is not in the pool. Primary rule: the
           source's own Al-Mg convex hull evaluated at x_Mg = 89/229 (the best
           the source can say at that composition from the Al-Mg structures it
           computed). Alternative rule: the epsilon phase (mp_17659_Mg23Al30,
           x_Mg = 0.434), the observed low-temperature phase nearest in
           composition. Both are reported.

Run with a Python that has ase and scipy (<venv>/bin/python).
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np
from scipy.stats import kendalltau

import lab_receivers_study as L

HERE = Path(__file__).resolve().parent
RES = HERE / "results"
OUT = RES / "lab_receivers"
KJ = L.KJ_PER_EV
X_BETA = 89 / 229

MEAS = {
    "beta": {"BP1970": -2.47, "PH1978": -1.22},
    "gamma": {"BP1970": -3.80, "PH1978": -1.01},
    "Mg2Si": {"JANAF": -77.822 / 3, "Feufel1997": -22.2},
}
MEAS_REF = {
    "BP1970": "Brown and Pratt 1970, doi:10.1007/bf03037810, via Mezbahul-Islam et al. 2014 Table 4",
    "PH1978": "Predel and Huelse 1978, Z. Metallkd. 69, 661, via Mezbahul-Islam et al. 2014 Table 4",
    "JANAF": "NIST-JANAF table Mg-034, doi:10.18434/T42S31",
    "Feufel1997": "Feufel et al. 1997, doi:10.1016/s0925-8388(96)02655-2, as quoted by Jund et al. 2013, doi:10.1088/0953-8984/25/3/035403",
}
SOURCES = L.DONORS + ["gpaw_r2scan", "archive_mp"]
LABEL = dict(L.LABEL, gpaw_r2scan="GPAW r2SCAN, 0 K", archive_mp="MP archive")
GAMMA, MG2SI, EPS = "mp_2151_Mg17Al12", "mp_1367_Mg2Si", "mp_17659_Mg23Al30"


def source_values():
    """{source: {cid: E_f eV/atom}} using live_values first, then per-source records, then the MP index."""
    lv = json.loads((RES / "live_values.json").read_text())
    vals = {}
    for s in L.DONORS + ["gpaw_r2scan"]:
        d = {c: r["value"] for c, r in lv.get(s, {}).items() if r and r.get("value") is not None}
        if s in L.DONORS:
            for c, v in L.formation_per_atom(s).items():
                d.setdefault(c, v)
        vals[s] = d
    mp = {}
    for f in ("candidates_index.json", "candidates_index_observed.json"):
        for r in json.loads((HERE / f).read_text()):
            if r.get("mp_ef") is not None:
                mp[r["cid"]] = r["mp_ef"]
    vals["archive_mp"] = mp
    return vals


def mg_fraction(cid):
    n = L.counts(cid)
    if set(n) - {"Al", "Mg"}:
        return None
    return n.get("Mg", 0) / sum(n.values())


def hull_at(points, x):
    """Lower convex hull of (x, E) points (endpoints at E = 0 included) evaluated at x."""
    pts = sorted(set(points) | {(0.0, 0.0), (1.0, 0.0)})
    hull = []
    for p in pts:
        while len(hull) >= 2:
            (x1, y1), (x2, y2) = hull[-2], hull[-1]
            if (x2 - x1) * (p[1] - y1) - (y2 - y1) * (p[0] - x1) <= 0:
                hull.pop()
            else:
                break
        hull.append(p)
    for (x1, y1), (x2, y2) in zip(hull, hull[1:]):
        if x1 <= x <= x2:
            return y1 + (y2 - y1) * (x - x1) / (x2 - x1), [(x1, y1), (x2, y2)]
    raise ValueError(x)


def computed(vals):
    out = {}
    xcache = {}
    for s, d in vals.items():
        if GAMMA not in d or MG2SI not in d:
            continue
        pts, cids = [], {}
        for c, v in d.items():
            if c.startswith("clu_"):
                continue
            if c not in xcache:
                try:
                    xcache[c] = mg_fraction(c)
                except Exception:
                    xcache[c] = None
            x = xcache[c]
            if x is None or x in (0.0, 1.0):
                continue
            pts.append((x, v))
            cids[(x, v)] = c
        hb, facet = hull_at(pts, X_BETA)
        out[s] = {"gamma": d[GAMMA], "Mg2Si": d[MG2SI],
                  "beta_hull": hb, "beta_hull_facet": [[x, y, cids.get((x, y), "element")] for x, y in facet],
                  "beta_epsilon": d.get(EPS), "n_al_mg_structures": len(pts)}
    return out


def metrics(c, m):
    """c, m: dicts over beta, gamma, Mg2Si in eV/atom."""
    keys = ["beta", "gamma", "Mg2Si"]
    cv, mv = np.array([c[k] for k in keys]), np.array([m[k] for k in keys])
    t = kendalltau(cv, mv).statistic
    return {"error": {k: float(c[k] - m[k]) for k in keys},
            "mae": float(np.mean(np.abs(cv - mv))),
            "kendall_tau": float(t),
            "beta_gamma_order_agrees": bool(np.sign(c["gamma"] - c["beta"]) == np.sign(m["gamma"] - m["beta"])),
            "binding_scale": float(cv @ mv / (mv @ mv)),
            "beta_gamma_step_computed": float(c["gamma"] - c["beta"]),
            "beta_gamma_step_measured": float(m["gamma"] - m["beta"])}


def main():
    vals = source_values()
    comp = computed(vals)
    combos = list(itertools.product(MEAS["beta"], MEAS["gamma"], MEAS["Mg2Si"]))
    res = {"measured_kJ_per_mol_atoms": MEAS, "measured_refs": MEAS_REF, "x_mg_beta": X_BETA,
           "computed_eV_per_atom": comp, "combinations": {}, "summary": {}}
    for rule in ("beta_hull", "beta_epsilon"):
        res["combinations"][rule] = {}
        for b, g, s in combos:
            m = {"beta": MEAS["beta"][b] / KJ, "gamma": MEAS["gamma"][g] / KJ, "Mg2Si": MEAS["Mg2Si"][s] / KJ}
            name = f"beta:{b}|gamma:{g}|Mg2Si:{s}"
            block = {"measured_eV_per_atom": m, "sources": {}}
            for src, cv in comp.items():
                if cv[rule] is None:
                    continue
                block["sources"][src] = metrics({"beta": cv[rule], "gamma": cv["gamma"], "Mg2Si": cv["Mg2Si"]}, m)
            res["combinations"][rule][name] = block
        summ = {}
        for src in comp:
            rows = [blk["sources"][src] for blk in res["combinations"][rule].values() if src in blk["sources"]]
            if not rows:
                continue
            arr = lambda k: [r[k] for r in rows]
            summ[src] = {"label": LABEL.get(src, src), "n_combinations": len(rows),
                         "mae_min": min(arr("mae")), "mae_max": max(arr("mae")),
                         "tau_values": sorted(set(round(x, 3) for x in arr("kendall_tau"))),
                         "tau_mean": float(np.mean(arr("kendall_tau"))),
                         "order_agrees_in": int(sum(arr("beta_gamma_order_agrees"))),
                         "scale_min": min(arr("binding_scale")), "scale_max": max(arr("binding_scale"))}
        res["summary"][rule] = summ
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "measured_enthalpy_combinations.json").write_text(json.dumps(res, indent=1))
    write_md(res)


def write_md(res):
    L_ = ["# Measured formation-enthalpy receiver: all beta, gamma and Mg2Si combinations", "",
          "Script: `../../measured_enthalpy_combinations.py`. Data: `measured_enthalpy_combinations.json`.", "",
          "## Measurements (kJ/mol of atoms; eV/atom)", "", "| Compound | Measurement | kJ/mol-atoms | eV/atom | Source |",
          "|---|---|---|---|---|"]
    for ph, d in res["measured_kJ_per_mol_atoms"].items():
        for k, v in d.items():
            L_.append(f"| {ph} | {k} | {v:.2f} | {v / KJ:.4f} | {res['measured_refs'][k]} |")
    L_ += ["", "Brown and Pratt place gamma below beta; Predel and Huelse place beta below gamma. "
           "Mg2Si is the most strongly bound compound in every combination.", "",
           "## Computed values (eV/atom)", "",
           "| Source | gamma | Mg2Si | beta, hull at x_Mg = 0.389 | beta, epsilon proxy |", "|---|---|---|---|---|"]
    for s, c in res["computed_eV_per_atom"].items():
        e = "" if c["beta_epsilon"] is None else f"{c['beta_epsilon']:+.4f}"
        L_.append(f"| {LABEL.get(s, s)} | {c['gamma']:+.4f} | {c['Mg2Si']:+.4f} | {c['beta_hull']:+.4f} | {e} |")
    for rule, title in (("beta_hull", "beta from each source's Al-Mg hull (primary)"),
                        ("beta_epsilon", "beta mapped to the epsilon phase (alternative)")):
        L_ += ["", f"## Summary over the 8 combinations, {title}", "",
               "| Source | MAE range, eV/atom | Kendall tau values | beta-gamma order agrees | binding scale range |",
               "|---|---|---|---|---|"]
        for s, d in res["summary"][rule].items():
            L_.append(f"| {d['label']} | {d['mae_min']:.3f} to {d['mae_max']:.3f} | "
                      f"{', '.join(f'{t:+.2f}' for t in d['tau_values'])} | {d['order_agrees_in']} of {d['n_combinations']} | "
                      f"{d['scale_min']:.2f} to {d['scale_max']:.2f} |")
    L_ += ["", "Binding scale: least-squares slope of computed against measured through the origin "
           "(1 = the source reproduces the measured binding; below 1 = underbinding).", "",
           "## Per combination (primary rule): Kendall tau / MAE eV/atom", ""]
    blocks = res["combinations"]["beta_hull"]
    srcs = list(res["summary"]["beta_hull"])
    L_.append("| Combination | " + " | ".join(LABEL.get(s, s) for s in srcs) + " |")
    L_.append("|---|" + "---|" * len(srcs))
    for name, blk in blocks.items():
        L_.append(f"| {name} | " + " | ".join(
            f"{blk['sources'][s]['kendall_tau']:+.2f} / {blk['sources'][s]['mae']:.3f}" if s in blk["sources"] else ""
            for s in srcs) + " |")
    (OUT / "MEASURED_ENTHALPY_COMBINATIONS.md").write_text("\n".join(L_) + "\n")


if __name__ == "__main__":
    main()
