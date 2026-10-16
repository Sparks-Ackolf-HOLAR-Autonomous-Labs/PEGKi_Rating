#!/usr/bin/env python
"""Measured age-hardening receiver from Poznak, Thole and Sanders, Metals 8 (2018) 309.

Six alloys from one laboratory and one processing route (two solute levels,
Mg + Si about 0.7 and 1.4 wt%, each at Mg/Si about 0.3, 0.7 and 2.4 by weight),
water quenched, naturally aged 0.5 h, aged at 175 C. Values digitised from
Figs. 3 and 10 (results/lab_receivers/poznak2018_ageing.json).

Donor prediction per alloy and source (precipitation driving force):
  for each observed precipitate P in {beta'' Mg5Si6, beta' Mg9Si5, beta Mg2Si}
  the alloy can form n_P = min(x_Mg / a, x_Si / b) formula units per alloy atom
  (minority solute limits); the released energy per alloy atom is
  n_P (a + b) E_f(P), with E_f the source's formation energy per atom.
  The prediction is the most negative value over P. More negative should mean
  harder: agreement with peak hardness = -tau.
The earlier rule (formation energy of the most stable compound at the closest
Mg/Si ratio, lab_receivers_study.py) is reported alongside.

Run with a Python that has ase and scipy (<venv>/bin/python).
"""
from __future__ import annotations

import json
import math
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.stats import kendalltau

from lab_receivers_study import DONORS, formation_per_atom, mg_si_items
from lab_sequence_study import LABEL

HERE = Path(__file__).resolve().parent
RES = HERE / "results"
OUT = RES / "lab_receivers"
M = {"Al": 26.9815, "Mg": 24.305, "Si": 28.085}
PHASES = {"mp_568306_Mg5Si6": (5, 6), "mp_542854_Mg9Si5": (9, 5), "mp_1367_Mg2Si": (2, 1)}


def at_frac(mg_wt, si_wt):
    n = {"Mg": mg_wt / M["Mg"], "Si": si_wt / M["Si"], "Al": (100 - mg_wt - si_wt) / M["Al"]}
    s = sum(n.values())
    return n["Mg"] / s, n["Si"] / s


def driving_force(ef, lv, a):
    xm, xs = at_frac(a["Mg_wt"], a["Si_wt"])
    best, phase = None, None
    for cid, (p, q) in PHASES.items():
        e = (lv.get(cid) or {}).get("value")
        if e is None:
            e = ef.get(cid)
        if e is None:
            continue
        d = min(xm / p, xs / q) * (p + q) * e
        if best is None or d < best:
            best, phase = d, cid
    return best, phase


def closest_ratio(items, a):
    xm, xs = at_frac(a["Mg_wt"], a["Si_wt"])
    r = xm / xs
    dmin = min(abs(math.log(q / r)) for q, _ in items.values())
    near = {c: e for c, (q, e) in items.items() if abs(abs(math.log(q / r)) - dmin) < 1e-9}
    c = min(near, key=near.get)
    return near[c], c


def tau_ci(x, y, nboot=2000, seed=0):
    """Kendall tau-b with a 90% interval from resampling the alloys."""
    x, y = np.asarray(x), np.asarray(y)
    t = kendalltau(x, y).statistic
    rng = np.random.default_rng(seed)
    bs = []
    for _ in range(nboot):
        i = rng.integers(0, len(x), len(x))
        if len(set(i)) < 3:
            continue
        v = kendalltau(x[i], y[i]).statistic
        if math.isfinite(v):
            bs.append(v)
    return float(t), [float(np.percentile(bs, 5)), float(np.percentile(bs, 95))]


def main():
    data = json.loads((OUT / "poznak2018_ageing.json").read_text())
    alloys = data["alloys"]
    lv = json.loads((RES / "live_values.json").read_text())
    hv = [a["peak_HV_NA05"] for a in alloys]
    loss = [a["peak_HV_NA05"] - a["peak_HV_NA4"] for a in alloys]
    out = {"receiver": data["citation"], "doi": data["doi"], "n_alloys": len(alloys),
           "measured": {a["name"]: {"peak_HV": a["peak_HV_NA05"], "NA_loss_HV": round(l, 1),
                                    "time_to_peak_h": a["time_to_peak_h"],
                                    "x_Mg": round(at_frac(a["Mg_wt"], a["Si_wt"])[0], 5),
                                    "x_Si": round(at_frac(a["Mg_wt"], a["Si_wt"])[1], 5)}
                        for a, l in zip(alloys, loss)},
           "donors": {}}
    for src in DONORS:
        ef = formation_per_atom(src)
        lvs = lv.get(src, {})
        df = [driving_force(ef, lvs, a) for a in alloys]
        if any(d is None for d, _ in df):
            continue
        pred = [d for d, _ in df]
        t, ci = tau_ci(pred, hv)
        # within each solute level: pairs of alloys with the same Mg + Si
        same = [(i, j) for i, j in combinations(range(6), 2) if alloys[i]["name"][:2] == alloys[j]["name"][:2]]
        conc = sum(1 for i, j in same if (pred[i] - pred[j]) * (hv[i] - hv[j]) < 0)
        items = mg_si_items(ef, lvs)
        old = [closest_ratio(items, a) for a in alloys] if items else None
        t_old = kendalltau([e for e, _ in old], hv).statistic if old else None
        out["donors"][src] = {
            "label": LABEL.get(src, src),
            "driving_force_eV_per_alloy_atom": {a["name"]: round(d, 6) for a, (d, _) in zip(alloys, df)},
            "phase_used": {a["name"]: p for a, (_, p) in zip(alloys, df)},
            "agreement_peak_HV": -t, "agreement_peak_HV_ci90": [-ci[1], -ci[0]],
            "within_level_pairs_agreeing": f"{conc}/{len(same)}",
            "closest_ratio_rule_agreement_peak_HV": None if t_old is None else -float(t_old)}
    # composition-only baseline: beta'' (Mg5Si6) formula units the alloy can form, no source energy
    base = [-min(at_frac(a["Mg_wt"], a["Si_wt"])[0] / 5, at_frac(a["Mg_wt"], a["Si_wt"])[1] / 6) * 11 for a in alloys]
    tb, cib = tau_ci(base, hv)
    out["composition_only_beta_dprime"] = {"agreement_peak_HV": -tb, "agreement_peak_HV_ci90": [-cib[1], -cib[0]],
                                           "note": "ordering by Mg5Si6 solute available, identical for every source "
                                                   "with a negative beta'' formation energy"}
    out["phase_energies_eV_per_atom"] = {}
    for src in DONORS:
        ef = formation_per_atom(src)
        out["phase_energies_eV_per_atom"][src] = {cid: ((lv.get(src, {}).get(cid) or {}).get("value") or ef.get(cid))
                                                  for cid in PHASES}
    (OUT / "poznak_receiver.json").write_text(json.dumps(out, indent=1))
    print(json.dumps({k: {kk: v[kk] for kk in ("agreement_peak_HV", "agreement_peak_HV_ci90",
                      "within_level_pairs_agreeing", "closest_ratio_rule_agreement_peak_HV")}
                      for k, v in out["donors"].items()}, indent=1))
    print(json.dumps(out["measured"], indent=1))


if __name__ == "__main__":
    main()
