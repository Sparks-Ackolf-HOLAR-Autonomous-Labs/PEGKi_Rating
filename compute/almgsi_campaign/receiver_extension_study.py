#!/usr/bin/env python
"""PEGKi transfer scores and measured campaigns for the added Al-Mg-Si receivers.

Receivers that were available but unused (MACE-MP-0 MD at 300 K, the MP
archive) and the added receivers (quasi-harmonic free energies,
CALPHAD, r2SCAN) are scored the same way as the DFT receivers in
dft_receiver_study.py: Plackett-Luce ratings of composition classes per
source, Delta(D -> R) = sigma + rho - omega for every donor D, and the
measured early share of R's ten most stable compounds found by the
learned-trust campaign with donor D (campaign.py, mean over k = 1 and 4,
20 seeds). Small receivers without campaigns (measured enthalpies, age
hardening, vacancy binding) are reported by Kendall tau only, in their own
scripts.

Output: results/analysis/receiver_extension.{json,md}
"""
from __future__ import annotations

import json
import math
from pathlib import Path

from scipy.stats import spearmanr

import temperature_study as ts

HERE = Path(__file__).resolve().parent
RES = HERE / "results"

# source name -> (file relative to HERE, key inside file or None)
EXTRA = {
    "mace_md_300": ("results_T/T300/live_values.json", "mace_md"),
    **{f"{e}_{T}": (f"results_T/T{T}/live_values.json", e)
       for T in (453, 803) for e in ("mace_md", "lammps_meam_md", "gpumd_nep89_md")},
    **{f"qha_{T}": (f"results/qha/qha_T{T}.json", None) for T in (300, 453, 803)},
    **{f"calphad_{T}": (f"results/calphad/calphad_T{T}.json", None) for T in (300, 453, 803)},
    "qe_r2scan": ("results/qe_r2scan_table.json", None),
}
BASE = ("mace_relax", "gpumd_nep89_md", "lammps_meam_md", "qe_pbe", "qe_hse06")
LABEL = {"mace_relax": "MACE-MP-0, 0 K", "gpumd_nep89_md": "NEP89 MD, 300 K", "lammps_meam_md": "MEAM MD, 300 K",
         "qe_pbe": "QE PBE", "qe_hse06": "QE HSE06", "archive_mp": "MP archive", "mace_md_300": "MACE-MP-0 MD, 300 K",
         "qe_r2scan": "QE r2SCAN"}
for T in (300, 453, 803):
    LABEL[f"qha_{T}"] = f"MACE-MP-0 quasi-harmonic, {T} K"
    LABEL[f"calphad_{T}"] = f"CALPHAD, {T} K"
    for e, n in (("mace_md", "MACE-MP-0 MD"), ("lammps_meam_md", "MEAM MD"), ("gpumd_nep89_md", "NEP89 MD")):
        LABEL.setdefault(f"{e}_{T}", f"{n}, {T} K")
# receiver -> campaign tag
RECEIVERS = {"mace_md_300": "mace_md_300_all", "archive_mp": "archive_mp_all",
             "qha_453": "qha_453_all", "qha_803": "qha_803_all",
             "calphad_453": "calphad_453_all", "calphad_803": "calphad_803_all",
             "qe_r2scan": "qe_r2scan_all"}


def clean(rows, idx):
    t = {}
    for cid, r in rows.items():
        if cid.startswith("_") or cid not in idx or idx[cid]["property"] != "formation_energy":
            continue
        if r.get("status", "ok") != "ok" or r.get("value") is None or not math.isfinite(r["value"]):
            continue
        t[cid] = {"value": float(r["value"]), "sem": float(r.get("sem") or 0.0),
                  "cost": float(r.get("core_s") or r.get("wall_s") or 0.0), "property": "formation_energy"}
    return t


def load_src():
    idx = ts.camp.load_index()
    lv = json.loads((RES / "live_values.json").read_text(encoding="utf-8"))
    src = {s: clean(lv.get(s, {}), idx) for s in BASE}
    src["archive_mp"] = {c: {"value": float(r["mp_ef"]), "sem": 0.0, "cost": 0.0, "property": "formation_energy"}
                         for c, r in idx.items() if r.get("mp_ef") is not None and r["property"] == "formation_energy"}
    for name, (path, key) in EXTRA.items():
        p = HERE / path
        if p.exists():
            d = json.loads(p.read_text(encoding="utf-8"))
            src[name] = clean(d[key] if key else d, idx)
    return {s: t for s, t in src.items() if t}, idx


def main():
    src, idx = load_src()
    spaces, ltw = ts.pegki_spaces(src, idx)
    have = {R: tag for R, tag in RECEIVERS.items() if (RES / "campaigns" / tag / "summary.json").exists()}
    D = ts.pegki_deltas(spaces, ltw, list(have))
    out = {"receivers": {}, "n_items": {s: len(t) for s, t in src.items()}}
    lines = ["# Added Al-Mg-Si receivers: PEGKi transfer scores and measured campaigns", "",
             "| Receiver | Donor | tau (pool) | Delta | learned-trust early share | fixed-trust cost to 8 of 10 (k=4) |",
             "|---|---|---|---|---|---|"]
    for R, tag in have.items():
        summ = json.loads((RES / "campaigns" / tag / "summary.json").read_text(encoding="utf-8"))
        arms = summ["arms"]
        rows = {}
        for d, ag in summ["agreement"].items():
            # campaign donor names for MD at temperature -> source names used here
            alias = d.replace("meam_md_", "lammps_meam_md_").replace("nep89_md_", "gpumd_nep89_md_")
            key = f"{alias}->{R}"
            if f"learned:{d}|k1" not in arms:
                continue
            early = 0.5 * (arms[f"learned:{d}|k1"]["early_recall_area"]["mean"]
                           + arms[f"learned:{d}|k4"]["early_recall_area"]["mean"])
            c8 = arms[f"fixed:{d}|k4"]["cost_frac_to_8of10"]["mean"]
            rows[d] = {**D.get(key, {}), "tau_full": ag["tau_full"], "learned_early_share": early,
                       "fixed_cost_frac_to_8of10_k4": c8}
            dl = D[key]["delta"] if key in D else float("nan")
            lines.append(f"| {LABEL.get(R, R)} | {LABEL.get(alias, alias)} | {ag['tau_full']:+.2f} | {dl:+.3f} | "
                         f"{early:.3f} | {c8:.3f} |")
        scored = [d for d in rows if "delta" in rows[d]]
        rs = (spearmanr([rows[d]["delta"] for d in scored], [rows[d]["learned_early_share"] for d in scored]).statistic
              if len(scored) >= 3 else float("nan"))
        top_d = max(scored, key=lambda d: rows[d]["delta"]) if scored else None
        top_m = max(rows, key=lambda d: rows[d]["learned_early_share"])
        cold = 0.5 * (arms["surrogate|k1"]["early_recall_area"]["mean"] + arms["surrogate|k4"]["early_recall_area"]["mean"])
        out["receivers"][R] = {"donors": rows, "spearman_delta_vs_early": float(rs), "top_by_delta": top_d,
                               "top_measured": top_m, "no_donor_early_share": cold, "n_pool": summ["n_pool"],
                               "shapley": summ.get("shapley", {})}
        lines += ["", f"{LABEL.get(R, R)} (n = {summ['n_pool']}): no-donor early share {cold:.3f}; "
                      f"Spearman(Delta, measured) = {rs:+.2f}; top by Delta = {LABEL.get(top_d, top_d)}; "
                      f"top measured = {LABEL.get(top_m, top_m)}", ""]
        for mode, sh in summ.get("shapley", {}).items():
            lines.append(f"- Shapley ({mode}): " + ", ".join(f"{LABEL.get(k.replace('meam_md_', 'lammps_meam_md_').replace('nep89_md_', 'gpumd_nep89_md_'), k)} {v:+.3f}" for k, v in
                                                            sorted(sh["shapley"].items(), key=lambda x: -x[1]))
                         + f"; all donors {sh['v_all']:.3f}, best coalition {sh['v_best']:.3f}")
        lines.append("")
    (RES / "analysis" / "receiver_extension.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    (RES / "analysis" / "receiver_extension.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
