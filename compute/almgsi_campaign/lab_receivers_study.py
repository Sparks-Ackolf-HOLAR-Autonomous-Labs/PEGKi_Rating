#!/usr/bin/env python
"""Two laboratory receivers for the computed Al-Mg-Si campaign.

Receiver 1, measured formation enthalpies (results/lab_receivers/measured_enthalpy.json):
  open machine-readable compilations were searched for compounds made only of
  Al, Mg and Si (raw/compilation_hits.json). Only Mg2Si is present (NIST-JANAF
  table Mg-034, carried into the Kingsbury/Wang compilation distributed with
  matminer). Mapping rule: the compilation's likely MP id (mp-1367) equals the
  candidate's MP id; the same cid is also the lowest-PBE Mg2Si candidate.
  Metrics: signed and absolute error per source (n = 1, so no rank statistic).

Receiver 2, age hardening (results/lab_receivers/age_hardening.json):
  each alloy is an item, ordered by the measured outcome. A donor's prediction
  for an alloy uses the formation energy (eV/atom) of the most stable
  Mg-Si-containing candidate whose atomic Mg/Si ratio is closest to the
  alloy's (observed phases beta'' and beta' included through the per-source
  record directories). Kendall tau-b of each donor's prediction against the
  measured outcome is reported with n.

Run with a Python that has ase and scipy (<venv>/bin/python).
"""
from __future__ import annotations

import json
import math
from collections import Counter
from pathlib import Path

from ase.io import read
from scipy.stats import kendalltau

from lab_sequence_study import REFS, SOURCES, LABEL, load_records

HERE = Path(__file__).resolve().parent
RES = HERE / "results"
OUT = RES / "lab_receivers"
KJ_PER_EV = 96.485332  # kJ/mol per eV/atom (CODATA F/1000)
DONORS = ["qe_pbe", "qe_hse06", "mace_relax", "mace_md", "lammps_meam_md", "gpumd_nep89_md",
          "lammps_meam_relax", "gpumd_nep89_relax"]


def counts(cid):
    return Counter(read(str(HERE / "candidates" / f"{cid}.xyz"), format="extxyz").get_chemical_symbols())


def formation_per_atom(src):
    """{cid: E_f eV/atom} for every candidate with a record under the source's directories."""
    recs = load_records(SOURCES[src])
    if not all(c in recs for c in REFS.values()):
        return {}
    mu = {el: recs[c]["energy_per_atom"] for el, c in REFS.items()}
    out = {}
    for cid, r in recs.items():
        if cid.startswith("clu_"):
            continue
        try:
            n = counts(cid)
        except Exception:
            continue
        N = sum(n.values())
        out[cid] = (r["energy_per_atom"] * N - sum(k * mu[el] for el, k in n.items())) / N
    return out


def receiver1(ef):
    lv = json.loads((RES / "live_values.json").read_text())
    # JANAF Mg-034, 298.15 K: delta-f H = -77.822 kJ/mol of Mg2Si (3 atoms); raw/janaf_Mg-034_Mg2Si.txt
    dh_kj_fu = -77.822
    val = dh_kj_fu / 3 / KJ_PER_EV
    meas = {"mp_1367_Mg2Si": {
        "property": "formation_energy", "value": val, "sem": None,
        "unit": "eV/atom", "temperature_K": 298.15, "formula": "Mg2Si",
        "source_ref": "NIST-JANAF Thermochemical Tables, table Mg-034 (doi:10.18434/T42S31); "
                      "as carried in the Wang/Kingsbury compilation (matminer expt_formation_enthalpy_kingsbury, "
                      "doi:10.1038/s41598-021-94550-5), value -0.26886 eV/atom, likely_mpid mp-1367",
        "method": "critically evaluated standard enthalpy of formation (JANAF compilation; uncertainty not given "
                  "in the machine-readable table or the compilation)",
        "mapping": "compilation likely_mpid mp-1367 equals candidate mpid; also the lowest-PBE Mg2Si candidate",
        "status": "ok"}}
    rows = {}
    for src in DONORS + ["qe_hse06_qhalf"]:
        v = (lv.get(src, {}).get("mp_1367_Mg2Si") or {}).get("value")
        if v is None:
            v = ef.get(src, {}).get("mp_1367_Mg2Si")
        if v is None:
            continue
        rows[src] = {"value": v, "error": v - val, "abs_error": abs(v - val), "n": 1}
    idx = {d["cid"]: d for d in json.loads((HERE / "candidates_index.json").read_text())}
    mp = idx["mp_1367_Mg2Si"]["mp_ef"]
    rows["archive_mp"] = {"value": mp, "error": mp - val, "abs_error": abs(mp - val), "n": 1}
    # lowest-energy Mg2Si polymorph per source (alternative mapping rule)
    alt = {}
    for src in DONORS:
        cand = {c: (lv.get(src, {}).get(c) or {}).get("value") for c in
                ("mp_1367_Mg2Si", "mp_1101851_Mg2Si", "mp_1074652_Mg2Si")}
        cand = {c: x for c, x in cand.items() if x is not None}
        if cand:
            c = min(cand, key=cand.get)
            alt[src] = {"cid": c, "value": cand[c], "abs_error": abs(cand[c] - val)}
    return meas, {"n_compounds": 1, "kendall_tau": None,
                  "note": "n = 1: no rank statistic is defined; error per source only",
                  "per_source": rows, "lowest_polymorph_rule": alt}


def mg_si_items(ef_src, lv_src):
    """Mg-Si-containing candidates with their atomic Mg/Si ratio and this source's E_f."""
    out = {}
    cids = set(lv_src) | set(ef_src)
    for cid in cids:
        if cid.startswith("clu_"):
            continue
        v = (lv_src.get(cid) or {}).get("value")
        if v is None:
            v = ef_src.get(cid)
        if v is None:
            continue
        try:
            n = counts(cid)
        except Exception:
            continue
        if n.get("Mg", 0) and n.get("Si", 0):
            out[cid] = (n["Mg"] / n["Si"], v)
    return out


def receiver2(ef):
    p = OUT / "age_hardening.json"
    if not p.exists():
        return None
    ah = json.loads(p.read_text())
    lv = json.loads((RES / "live_values.json").read_text())
    sign = {"peak_hardness_HV": -1, "time_to_peak_h": 1, "dsc_beta_dprime_peak_C": 1}
    res = {"mapping_rule": ah.get("mapping_rule"), "agreement_sign": sign, "subsets": {}}
    items = {src: mg_si_items(ef.get(src, {}), lv.get(src, {})) for src in DONORS}
    for subset, alloys in (("sn_free", [a for a in ah["alloys"] if a["sn_free"]]), ("all_six", ah["alloys"])):
        res["subsets"][subset] = {}
        for outcome in ah["outcome_keys"]:
            meas = {a["name"]: a["outcomes"][outcome]["value"] for a in alloys}
            names = list(meas)
            block = {"n": len(names), "measured": meas, "donors": {}}
            for src in DONORS:
                it = items[src]
                if not it:
                    continue
                pred, picked = {}, {}
                for a in alloys:
                    r = a["mg_si_atomic_ratio"]
                    dmin = min(abs(math.log(q / r)) for q, _ in it.values())
                    near = {c: e for c, (q, e) in it.items() if abs(abs(math.log(q / r)) - dmin) < 1e-9}
                    c = min(near, key=near.get)
                    pred[a["name"]], picked[a["name"]] = near[c], c
                t = kendalltau([pred[n] for n in names], [meas[n] for n in names])
                tau = float(t.statistic) if math.isfinite(t.statistic) else None
                block["donors"][src] = {"label": LABEL.get(src, src), "prediction_eV_per_atom": pred,
                                        "picked_cid": picked, "kendall_tau_b": tau,
                                        "agreement": None if tau is None else sign[outcome] * tau,
                                        "p_value": float(t.pvalue) if math.isfinite(t.pvalue) else None,
                                        "n": len(names)}
            res["subsets"][subset][outcome] = block
    return res


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    ef = {s: formation_per_atom(s) for s in DONORS}
    meas, r1 = receiver1(ef)
    (OUT / "measured_enthalpy.json").write_text(json.dumps(meas, indent=1))
    r2 = receiver2(ef)
    (OUT / "lab_receivers_metrics.json").write_text(json.dumps({"measured_enthalpy": r1, "age_hardening": r2}, indent=1))
    print(json.dumps(r1, indent=1))
    if r2:
        for sub, blocks in r2["subsets"].items():
            for o, b in blocks.items():
                print(sub, o, "n =", b["n"], b["measured"])
                for s, d in b["donors"].items():
                    print(f"  {s:18s} tau={d['kendall_tau_b']}  agreement={d['agreement']}  picks={d['picked_cid']}  pred={ {k: round(v, 3) for k, v in d['prediction_eV_per_atom'].items()} }")


if __name__ == "__main__":
    main()
