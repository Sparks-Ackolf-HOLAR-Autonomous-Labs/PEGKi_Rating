"""Property-family receivers: PEGKi transfer scores vs paired campaigns.

TENKi experiment 28 used the bulk-modulus (JARVIS, AFLOW, MP), Vickers
hardness (Chen-Niu, JARVIS) and Al-Mg-Si formation-energy (AFLOW, MP)
channels only as zero-shot targets for the band-gap-trained acquisition
policy. Here every one of them becomes a RECEIVER of donor guidance:

  bulk_modulus      3 sources -> 6 ordered donor->receiver pairs
  formation_energy  2 sources -> 2 ordered pairs (lower is better)
  hardness          1 source (exp 28 skipped exp_hardness: 0 classes of
                    >= 30 records) -> single-source family, used only as a
                    receiver of the three bulk-modulus donors, restricted to
                    compositions (reduced formulas) the hardness channel
                    also contains. Motivation: the Chen-Niu model computes
                    Vickers hardness from the shear and bulk moduli
                    (Pugh ratio), so a bulk-modulus class ordering is a
                    physically plausible prior for hardness classes.

Per ordered pair A -> B:
  (1) PEGKi: Plackett-Luce class ratings per channel (fit_space of
      set_operation_study.py, the same fit as compute/almgsi_campaign/
      temperature_study.py), then Delta = sigma + rho - omega with the
      transfer_equation.py components (learned_transfer_weights W0 = 1,1,1).
  (2) TENKi: paired campaigns from experiment 26 (_job imported, not
      reimplemented): receiver B is the ground truth, donor A's class
      ordering is an anchored prior with FIXED trust (lambda = 8, "full")
      and LEARNED trust (prequential lambda, "selective"); matched seeds
      across cold/fixed/learned; pairwise and k = 4 arms at equal cost;
      budget 240 evaluations, checkpoints every 24 (exp 26 defaults).
      Measured outcome T_early = mean over the first half of the
      checkpoints of tau(condition) - tau(cold), averaged over seeds.
  (3) Spearman(Delta, T_early) over all ordered pairs, and per receiver
      whether the top-Delta donor is the best measured donor.

Classes are exp 28's: anion classes from exp 27's classify_formula for the
CSV channels (>= 30 records per class), chemical subsystems for formation
energy (>= 10 records). Exp 26 hard-codes higher_is_better=True, so the
formation-energy values are passed NEGATED (-dHf: more stable = higher);
the script asserts that this reproduces exp 28's lower-is-better
reference ranking. No invented observations: every campaign event draws
one archived record per participating class.

Usage (from shared/, shared uv env):
  UV_PROJECT_ENVIRONMENT=<venv> uv run --no-sync \
      python scripts/property_family_receivers.py [--seeds 20] [--workers 16]
Outputs: output/property_family_receivers/{results.json,results.md,
figure.png} plus per-(receiver, seed) checkpoints. ASCII only.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
SHARED = HERE.parent
SO = SHARED.parent
EXP_DIR = SO / "tenki" / "TENKi" / "experiments"
sys.path.insert(0, str(SHARED))
os.environ.setdefault("PEGKI_ROOT", str(SHARED))


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


exp26 = _load("exp26_paired_campaign", EXP_DIR / "26_paired_campaign.py")
exp28 = _load("exp28_cross_property_policy",
              EXP_DIR / "28_cross_property_policy.py")
setop = _load("set_operation_study", HERE / "set_operation_study.py")
ltw = setop.ltw

LABEL = {"bulk_jarvis": "JARVIS bulk modulus",
         "bulk_aflow": "AFLOW bulk modulus",
         "bulk_mp": "MP bulk modulus",
         "hardness_chen_niu": "JARVIS Chen-Niu hardness",
         "fe_aflow_dHf": "AFLOW formation energy",
         "fe_mp_dHf": "MP formation energy"}
BULK = ("bulk_jarvis", "bulk_aflow", "bulk_mp")
FE = ("fe_aflow_dHf", "fe_mp_dHf")
HARD = "hardness_chen_niu"
CSV_SPEC = {"bulk_jarvis": ("dft_3d_elastic.csv", "formula", "bulk_gpa"),
            "bulk_aflow": ("aflow_elastic.csv", "compound", "bulk_gpa"),
            "bulk_mp": ("mp_elastic.csv", "compound", "bulk_gpa"),
            HARD: ("dft_3d_elastic.csv", "formula", "chen_niu_vickers")}
SHARED_TAG = "@hardness_compositions"


def lab(name):
    base, _, tag = name.partition("@")
    return LABEL.get(base, base) + (" (shared compositions)" if tag else "")


# ------------------------------------------------------------------ data

def reduced(formula):
    try:
        from pymatgen.core import Composition
        return Composition(formula).reduced_formula
    except Exception:
        return None


def load_csv_records(name, jdir):
    """[(reduced formula, class, value)] with exp 28's filters (value > 0)."""
    fname, fcol, vcol = CSV_SPEC[name]
    out = []
    with open(jdir / fname, "r", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            v = (row.get(vcol) or "").strip()
            f = row.get(fcol) or ""
            if not v or not f:
                continue
            try:
                val = float(v)
            except ValueError:
                continue
            if val <= 0.0:
                continue
            out.append((reduced(f), exp28.classify_formula(f), val))
    return out


def classes_from(records, min_size, keep=None):
    raw = {}
    for f, c, v in records:
        if keep is not None and f not in keep:
            continue
        raw.setdefault(c, []).append(v)
    return {c: np.asarray(sorted(v), dtype=float)
            for c, v in raw.items() if len(v) >= min_size}


def build_channels(min_class_size, fe_min_class_size):
    class A:
        pass
    a = A()
    a.min_class_size, a.fe_min_class_size = min_class_size, fe_min_class_size
    a.jarvis_dir, a.tracka_jsonl = None, None
    evaluated, skipped = exp28.build_target_channels(a)
    ch, meta = {}, {}
    for name, spec in evaluated.items():
        hib = spec["higher_is_better"]
        cls = spec["classes"]
        # exp 26 is higher-is-better only: negate lower-is-better channels
        conv = {c: np.asarray(sorted(v if hib else -np.asarray(v)), dtype=float)
                for c, v in cls.items()}
        assert exp26.reference_ranking(conv) == exp28.reference_ranking(cls, hib), name
        ch[name] = conv
        meta[name] = {"family": spec["family"], "higher_is_better": hib,
                      "class_sizes": spec["class_sizes"],
                      "reference_ranking": exp28.reference_ranking(cls, hib)}
    return ch, meta, skipped


def hardness_donors(ch, min_class_size):
    """Bulk donors restricted to compositions present in the hardness channel."""
    jdir = SHARED / "data" / "jarvis"
    hard = load_csv_records(HARD, jdir)
    # consistency with exp 28's hardness classes
    chk = classes_from(hard, min_class_size)
    assert {c: len(v) for c, v in chk.items()} == \
        {c: len(v) for c, v in ch[HARD].items()}, "hardness classes differ from exp 28"
    keep = {f for f, _, _ in hard if f}
    out, info = {}, {}
    for b in BULK:
        recs = load_csv_records(b, jdir)
        full = classes_from(recs, min_class_size)
        assert {c: len(v) for c, v in full.items()} == \
            {c: len(v) for c, v in ch[b].items()}, f"{b} classes differ from exp 28"
        sub = classes_from(recs, min_class_size, keep)
        name = b + SHARED_TAG
        out[name] = sub
        info[name] = {"n_records_full": int(sum(len(v) for v in full.values())),
                      "n_records_shared": int(sum(len(v) for v in sub.values())),
                      "n_compositions_shared": len({f for f, _, _ in recs if f in keep}),
                      "class_sizes": {c: int(len(v)) for c, v in sub.items()}}
    return out, info, len(keep)


# ------------------------------------------------------------------ PEGKi

def pegki(spaces_in, pairs, n_events, n_boot, seed):
    spaces = {}
    for name, cls in spaces_in.items():
        ks, _ = setop.fit_space(name, {c: list(v) for c, v in cls.items()},
                                n_events, n_boot, seed)
        if ks is not None:
            spaces[name] = ks
    D = {}
    for a, b in pairs:
        s = float(ltw.compute_structural_compatibility(spaces[a], spaces[b]))
        r = float(ltw.compute_mapping_efficiency(spaces[a], spaces[b]))
        o = float(ltw.compute_transfer_overhead(spaces[a], spaces[b]))
        D[(a, b)] = {"sigma": s, "rho": r, "omega": o,
                     "delta": float(np.dot(ltw.W0, [s, r, -o]))}
    ratings = {n: dict(sorted(ks.policy_mus.items(), key=lambda kv: -kv[1]))
               for n, ks in spaces.items()}
    return D, ratings


# ------------------------------------------------------------------ main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=20)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--kway", type=int, default=4)
    ap.add_argument("--budget", type=int, default=240)
    ap.add_argument("--budget-step", type=int, default=24)
    ap.add_argument("--lam", type=float, default=8.0)
    ap.add_argument("--min-class-size", type=int, default=30)
    ap.add_argument("--fe-min-class-size", type=int, default=10)
    ap.add_argument("--pl-events", type=int, default=400)
    ap.add_argument("--pl-boot", type=int, default=20)
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--prior-shapes", nargs="+", default=["even", "margin"],
                    choices=list(exp26.PRIOR_SHAPES))
    ap.add_argument("--output-dir", default=str(SHARED / "output" / "property_family_receivers"))
    args = ap.parse_args()
    out = Path(args.output_dir)
    ck_dir = out / "checkpoints"
    ck_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    ch, meta, skipped = build_channels(args.min_class_size, args.fe_min_class_size)
    hdon, hinfo, n_hard_comp = hardness_donors(ch, args.min_class_size)
    print({n: len(c) for n, c in ch.items()}, flush=True)
    print({n: v["n_records_shared"] for n, v in hinfo.items()}, flush=True)

    # receiver -> {channel name: classes} (receiver + its donors)
    groups = {b: {x: ch[x] for x in BULK} for b in BULK}
    groups.update({f: {x: ch[x] for x in FE} for f in FE if f in ch})
    groups[HARD] = {HARD: ch[HARD], **hdon}
    pairs = [(d, r) for r, g in groups.items() for d in g if d != r]

    D, ratings = pegki({**ch, **hdon}, pairs, args.pl_events, args.pl_boot, 11)

    checkpoints = list(range(args.budget_step, args.budget + 1, args.budget_step))
    half = max(1, len(checkpoints) // 2)
    shapes = args.prior_shapes
    jobs = []
    for shape in shapes:
        for r, g in groups.items():
            for s in range(args.seeds):
                seed = args.seed + 1000 * s
                if (ck_dir / f"{shape}_{r}_s{seed}.json").exists():
                    continue
                jobs.append((r, seed, g, args.kway, args.budget, checkpoints, args.lam, shape))
    print(f"[plan] {len(jobs)} (prior shape, receiver, seed) jobs", flush=True)

    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        for job, res in zip(jobs, ex.map(exp26._job, jobs)):
            p = ck_dir / f"{job[7]}_{res['receiver']}_s{res['seed']}.json"
            tmp = p.with_suffix(".tmp")
            tmp.write_text(json.dumps(res), encoding="utf-8")
            os.replace(tmp, p)
            print(f"[done] {job[7]} {res['receiver']} seed {res['seed']} "
                  f"({time.time() - t0:.0f} s)", flush=True)

    # ---- aggregate (matched seeds)
    boot = np.random.default_rng(args.seed)
    TRUSTS = ("full", "selective")
    CONDS = [f"{shape}|{t}|{a}" for shape in shapes for t in TRUSTS for a in exp26.ARMS]
    SHAPE_LABEL = {"even": "order prior", "margin": "margin prior"}
    TRUST_LABEL = {"full": "fixed trust", "selective": "learned trust"}

    def cond_label(c):
        sh, t, a = c.split("|")
        return f"{SHAPE_LABEL[sh]}, {TRUST_LABEL[t]}, " + ("pairwise" if a == "pair" else f"k = {args.kway}")

    cells, combined, cold_ref = {}, {}, {}
    for shape in shapes:
        for r in groups:
            rows = [json.loads((ck_dir / f"{shape}_{r}_s{args.seed + 1000 * s}.json").read_text(encoding="utf-8"))
                    for s in range(args.seeds)]
            for d in list(groups[r]) + [exp26.COMBINED]:
                if d == r:
                    continue
                for cond in ("full", "selective", "combined"):
                    for arm in exp26.ARMS:
                        key = f"{cond}|{arm}|{d}"
                        if key not in rows[0]["curves"]:
                            continue
                        early = np.array([np.mean(x["curves"][key][:half])
                                          - np.mean(x["curves"][f"cold|{arm}"][:half]) for x in rows])
                        final = np.array([x["curves"][key][-1] - x["curves"][f"cold|{arm}"][-1] for x in rows])
                        v = {"T_early": float(early.mean()), "T_early_ci": exp26.bootstrap_ci(early, boot),
                             "T_final": float(final.mean())}
                        if d == exp26.COMBINED:
                            combined.setdefault(r, {})[f"{shape}|{arm}"] = v
                        else:
                            cells.setdefault((d, r), {})[f"{shape}|{cond}|{arm}"] = v
            if shape == shapes[0]:
                cold_ref[r] = {"cold_early_tau": {a: float(np.mean([np.mean(x["curves"][f"cold|{a}"][:half])
                                                                    for x in rows])) for a in exp26.ARMS},
                               "cold_final_tau": {a: float(np.mean([x["curves"][f"cold|{a}"][-1] for x in rows]))
                                                  for a in exp26.ARMS}}
                meta[r].update(cold_ref[r])

    pair_rows = []
    for (d, r) in pairs:
        row = {"donor": d, "receiver": r, "family": meta[r]["family"],
               "cross_property": r == HARD, **D[(d, r)]}
        for c in CONDS:
            row[c] = cells[(d, r)][c]
        pair_rows.append(row)

    def corr(rows, c):
        if len(rows) < 3:
            return None
        v = spearmanr([x["delta"] for x in rows], [x[c]["T_early"] for x in rows]).statistic
        return None if not np.isfinite(v) else float(v)

    spearman = {c: {"all_pairs": corr(pair_rows, c),
                    "bulk_modulus_pairs": corr([x for x in pair_rows if x["receiver"] in BULK], c),
                    "hardness_receiver": corr([x for x in pair_rows if x["receiver"] == HARD], c)}
                for c in CONDS}
    routing = {}
    for r in groups:
        rr = [x for x in pair_rows if x["receiver"] == r]
        top = max(rr, key=lambda x: x["delta"])["donor"]
        pick = next(x for x in rr if x["donor"] == top)
        routing[r] = {}
        for c in CONDS:
            bestv = max(x[c]["T_early"] for x in rr)
            best = [x["donor"] for x in rr if x[c]["T_early"] >= bestv - 1e-9]
            routing[r][c] = {"top_delta_donor": top, "best_measured_donors": best,
                             "n_donors": len(rr), "match": (top in best) if len(rr) > 1 else None,
                             "tie_all_donors": len(rr) > 1 and len(best) == len(rr),
                             "pick_T_early": pick[c]["T_early"], "best_T_early": bestv,
                             "regret": bestv - pick[c]["T_early"]}
    multi = [r for r in routing if routing[r][CONDS[0]]["n_donors"] > 1]
    match_rate = {c: {"strict": float(np.mean([routing[r][c]["match"] and not routing[r][c]["tie_all_donors"]
                                               for r in multi])),
                      "ties_as_match": float(np.mean([routing[r][c]["match"] for r in multi])),
                      "n_receivers": len(multi)} for c in CONDS}

    res = {"config": vars(args), "checkpoints": checkpoints,
           "early_window_checkpoints": checkpoints[:half],
           "channels": meta, "skipped_channels": {n: {k: v for k, v in s.items() if k != "classes"}
                                                 for n, s in skipped.items()},
           "hardness_shared_compositions": {"n_hardness_compositions": n_hard_comp, "donors": hinfo},
           "class_ratings_PL": ratings,
           "condition_labels": {c: cond_label(c) for c in CONDS},
           "pairs": pair_rows,
           "combined_donor": combined,
           "spearman_delta_vs_T_early": spearman,
           "routing_top_delta_vs_best_measured": routing,
           "routing_match_rate": match_rate,
           "elapsed_s": round(time.time() - t0, 1)}
    (out / "results.json").write_text(json.dumps(res, indent=1), encoding="utf-8")

    # ---- markdown
    f2 = lambda v: "n/a" if v is None else f"{v:+.2f}"
    L = ["# Property-family receivers: PEGKi Delta vs measured paired-campaign gain", "",
         f"Receivers: bulk modulus (JARVIS, AFLOW, MP), Al-Mg-Si formation energy (AFLOW, MP), "
         f"and JARVIS Chen-Niu Vickers hardness (single source; donors are the three bulk-modulus "
         f"channels restricted to the {n_hard_comp} hardness compositions). Classes as in TENKi "
         f"experiment 28. Paired campaigns from experiment 26 (`_job`, imported): {args.seeds} matched "
         f"seeds, budget {args.budget}, checkpoints every {args.budget_step}, fixed trust lambda = "
         f"{args.lam}, learned trust prequential. T_early = mean Kendall tau gain over the cold start "
         f"in the first {half} checkpoints, with 95% bootstrap CI over seeds. Order prior = donor class "
         f"order only (exp 26 default); margin prior = donor median margins standardised to the same "
         f"spread (exp 26 --prior-shape margin).", "",
         "## PEGKi components", "", "| Donor -> receiver | sigma | rho | omega | Delta |", "|---|---|---|---|---|"]
    order = sorted(pair_rows, key=lambda x: (x["family"], x["receiver"], -x["delta"]))
    for x in order:
        L.append(f"| {lab(x['donor'])} -> {lab(x['receiver'])} | {x['sigma']:+.3f} | {x['rho']:+.3f} | "
                 f"{x['omega']:+.3f} | {x['delta']:+.3f} |")
    for shape in shapes:
        cs = [c for c in CONDS if c.startswith(shape + "|")]
        L += ["", f"## Measured T_early, {SHAPE_LABEL[shape]}", "",
              "| Donor -> receiver | Delta | " + " | ".join(cond_label(c).split(", ", 1)[1] for c in cs) + " |",
              "|---|---|" + "---|" * len(cs)]
        for x in order:
            L.append(f"| {lab(x['donor'])} -> {lab(x['receiver'])} | {x['delta']:+.3f} | "
                     + " | ".join(f"{x[c]['T_early']:+.3f} [{x[c]['T_early_ci'][0]:+.3f}, {x[c]['T_early_ci'][1]:+.3f}]"
                                  for c in cs) + " |")
        L.append("| all donors pooled (learned trust) -> receiver, pairwise / k-way | | "
                 + "; ".join(f"{lab(r)} {combined[r][shape + '|pair']['T_early']:+.3f} / "
                             f"{combined[r][shape + '|kway']['T_early']:+.3f}" for r in groups) + " |")
    L += ["", "## Spearman(Delta, T_early)", "",
          f"| Condition | all pairs (n = {len(pair_rows)}) | bulk-modulus pairs (n = 6) | hardness receiver (n = 3) |",
          "|---|---|---|---|"]
    for c in CONDS:
        s = spearman[c]
        L.append(f"| {cond_label(c)} | {f2(s['all_pairs'])} | {f2(s['bulk_modulus_pairs'])} | "
                 f"{f2(s['hardness_receiver'])} |")
    L += ["", "n/a: all donors gave the same T_early (constant input to the rank correlation).", "",
          "## Is the top-Delta donor the best measured donor?", "",
          "| Receiver | top-Delta donor | " + " | ".join(cond_label(c) for c in CONDS) + " |",
          "|---|---|" + "---|" * len(CONDS)]
    for r in groups:
        cellf = lambda q: ("single donor" if q["n_donors"] == 1 else
                           "tie (all donors equal)" if q["tie_all_donors"] else "yes" if q["match"] else
                           f"no ({', '.join(lab(b) for b in q['best_measured_donors'])}; regret {q['regret']:.3f})")
        L.append(f"| {lab(r)} | {lab(routing[r][CONDS[0]]['top_delta_donor'])} | "
                 + " | ".join(cellf(routing[r][c]) for c in CONDS) + " |")
    L.append(f"| match rate over the {len(multi)} multi-donor receivers, strict (ties as match) | | "
             + " | ".join(f"{match_rate[c]['strict']:.2f} ({match_rate[c]['ties_as_match']:.2f})" for c in CONDS) + " |")
    L += ["", "## Cold-start reference (no donor)", "",
          "| Receiver | classes | cold early tau (pair, k-way) | cold final tau (pair, k-way) |", "|---|---|---|---|"]
    for r in groups:
        m = meta[r]
        L.append(f"| {lab(r)} | {len(m['class_sizes'])} | {m['cold_early_tau']['pair']:+.3f}, "
                 f"{m['cold_early_tau']['kway']:+.3f} | {m['cold_final_tau']['pair']:+.3f}, {m['cold_final_tau']['kway']:+.3f} |")
    L += ["", "Hardness donors on shared compositions: "
          + "; ".join(f"{lab(n)} {v['n_records_shared']} of {v['n_records_full']} records" for n, v in hinfo.items())
          + ". The JARVIS bulk-modulus donor and the hardness receiver come from the same JARVIS entries "
            "(Chen-Niu hardness is computed from the JARVIS bulk and shear moduli), so that pair is not an "
            "independent source.", "",
          "Formation-energy values enter the campaigns negated (more stable = higher), reproducing the "
          "lower-is-better reference ranking of experiment 28 (asserted at runtime).",
          f"Wall time {res['elapsed_s']:.0f} s."]
    (out / "results.md").write_text("\n".join(L) + "\n", encoding="utf-8")

    # ---- figure (true print size, 3.27 in column)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(SO / "writings" / "papers_2026"))
    try:
        import paper_fig_style as pfs
        pfs.apply_style()
        cols = {"bulk_modulus": pfs.MP_BLUE, "formation_energy": pfs.EXPT_GOLD, "hardness": pfs.ACCENT_RED}
    except Exception:
        cols = {"bulk_modulus": "#2E5FA3", "formation_energy": "#B8860B", "hardness": "#C4452C"}
        plt.rcParams.update({"font.size": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
                             "legend.fontsize": 7, "axes.titlesize": 9, "figure.constrained_layout.use": True})
    fam_label = {"bulk_modulus": "bulk modulus", "formation_energy": "formation energy",
                 "hardness": "hardness, bulk-modulus donors"}
    show = [c for c in (f"even|selective|kway", f"margin|selective|kway") if c in CONDS]
    fig, axes = plt.subplots(len(show), 1, figsize=(3.27, 1.95 * len(show) + 0.3), sharex=True, squeeze=False)
    axes = axes[:, 0]
    for i, (ax, c) in enumerate(zip(axes, show)):
        for fam in ("bulk_modulus", "formation_energy", "hardness"):
            rr = [x for x in pair_rows if x["family"] == fam]
            if not rr:
                continue
            xs = np.array([x["delta"] for x in rr])
            ys = np.array([x[c]["T_early"] for x in rr])
            lo = ys - np.array([x[c]["T_early_ci"][0] for x in rr])
            hi = np.array([x[c]["T_early_ci"][1] for x in rr]) - ys
            ax.errorbar(xs, ys, yerr=[lo, hi], fmt="o", ms=3.5, lw=0.8, capsize=1.5,
                        color=cols[fam], label=fam_label[fam])
        ax.axhline(0, color="#8C8C8C", lw=0.6)
        s = spearman[c]["all_pairs"]
        ax.set_title(f"({'ab'[i]}) {SHAPE_LABEL[c.split('|')[0]]}, learned trust, k = {args.kway}",
                     loc="left")
        ax.text(0.98, 0.04, f"Spearman {f2(s)} (n = {len(pair_rows)})", transform=ax.transAxes,
                ha="right", va="bottom", fontsize=7)
        ax.set_ylabel("early tau gain")
    axes[-1].set_xlabel(r"PEGKi transfer score $\Delta$ (donor $\rightarrow$ receiver)")
    axes[0].legend(loc="center left", bbox_to_anchor=(0.0, 0.5), handletextpad=0.3, borderaxespad=0.2)
    fig.savefig(out / "figure.png", dpi=600, facecolor="white")
    plt.close(fig)
    print("\n".join(L))


if __name__ == "__main__":
    main()
