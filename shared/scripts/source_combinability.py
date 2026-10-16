"""
Multi-source set operations: which band-gap data sources can be combined?

Before any rating is fitted, set operations on the archived records answer a
practical question: can two (or more) sources be pooled, and at what level?

For every pair of sources (A, B) on the version-locked archive:

  * A cap B   : compounds (reduced formulas) both sources report; one value per
                compound and source (median over polymorph records);
  * coverage  : |A cap B| / |A cup B|, the share of the union that can be
                checked at all;
  * values    : after one rigid scissor shift (median of gap_B - gap_A on
                A cap B), the share of shared compounds whose gaps agree within
                TOL eV. A high share means the values themselves can be merged;
  * ordering  : Kendall tau-b of the gaps on A cap B. A high tau with a low
                value share means the ordering can be combined while the
                values cannot;
  * margins   : IQR(B) / IQR(A) on A cap B, how much B stretches the spacing
                between compounds relative to A.

Multi-source step: sources are grouped by average-linkage clustering on
(1 - tau); for each group the intersection over all members (the compounds on
which the whole group can be checked together) and the pairwise minimum tau
are reported.

Uncertainty: 200 bootstrap resamples of the shared compounds per pair (90%
percentile intervals). Every number is computed from archived records; no
value is created.

Usage (from shared/):
    uv run python scripts/source_combinability.py --universe <jsonl>
ASCII only.
"""

import argparse
import hashlib
import json
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.stats import kendalltau

MIN_GAP_EV = 0.01  # mirrors scripts/hypergraph_trial_bandgap.py
SOURCES = ["expt_zhuo", "snumat_hse", "jarvis_mbj", "mp_r2scan",
           "snumat_gga", "jarvis_optb88vdw", "mp_gga", "mp_gga_expsynth"]
LABEL = {
    "mp_gga": "MP GGA", "mp_r2scan": "MP r2SCAN",
    "mp_gga_expsynth": "MP GGA (ICSD-matched)", "expt_zhuo": "Experiment",
    "jarvis_optb88vdw": "JARVIS OptB88vdW", "jarvis_mbj": "JARVIS TB-mBJ",
    "snumat_gga": "SNUMAT GGA", "snumat_hse": "SNUMAT HSE06"}


def load(path):
    vals = {}
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            r = json.loads(line)
            f = (r.get("x") or {}).get("formula")
            g = (r.get("y") or {}).get("bandgap_eV")
            s = r.get("s")
            if not f or g is None or s not in LABEL or float(g) <= MIN_GAP_EV:
                continue
            vals.setdefault(s, {}).setdefault(f, []).append(float(g))
    return {s: {f: float(np.median(v)) for f, v in d.items()} for s, d in vals.items()}


def iqr(x):
    q = np.percentile(x, [25, 75])
    return float(q[1] - q[0])


def pair_stats(a, b, tol):
    x, y = np.asarray(a), np.asarray(b)
    tau = kendalltau(x, y).statistic
    d = y - x
    shift = float(np.median(d))
    frac = float(np.mean(np.abs(d - shift) <= tol))
    scale = iqr(y) / iqr(x) if iqr(x) > 0 else float("nan")
    return float(tau), shift, frac, float(scale)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--universe", default="output/bandgap_archive_build/bandgap_universe.jsonl")
    ap.add_argument("--output-dir", default="output/source_combinability")
    ap.add_argument("--tol", type=float, default=0.3)
    ap.add_argument("--min-shared", type=int, default=30)
    ap.add_argument("--n-boot", type=int, default=200)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--group-cut", type=float, default=0.2,
                    help="cut height on 1 - tau for the grouping")
    args = ap.parse_args()

    root = Path(__file__).resolve().parent.parent
    uni = (root / args.universe).resolve()
    sha = hashlib.sha256(uni.read_bytes()).hexdigest()[:16]
    V = load(uni)
    rng = np.random.default_rng(args.seed)

    pairs = {}
    for a, b in combinations(SOURCES, 2):
        shared = sorted(set(V[a]) & set(V[b]))
        union = set(V[a]) | set(V[b])
        rec = {"a": a, "b": b, "n_a": len(V[a]), "n_b": len(V[b]),
               "n_shared": len(shared), "n_union": len(union),
               "coverage": len(shared) / len(union) if union else 0.0,
               "subset": set(V[a]) <= set(V[b]) or set(V[b]) <= set(V[a])}
        if len(shared) >= args.min_shared:
            xa = np.array([V[a][f] for f in shared])
            xb = np.array([V[b][f] for f in shared])
            tau, shift, frac, scale = pair_stats(xa, xb, args.tol)
            boots = []
            for _ in range(args.n_boot):
                i = rng.integers(0, len(shared), len(shared))
                boots.append(pair_stats(xa[i], xb[i], args.tol))
            boots = np.array(boots)
            ci = lambda k: [float(np.nanpercentile(boots[:, k], 5)),
                            float(np.nanpercentile(boots[:, k], 95))]
            rec.update({"tau": tau, "tau_ci90": ci(0),
                        "scissor_shift_eV": shift, "scissor_shift_ci90": ci(1),
                        "value_agree_frac": frac, "value_agree_ci90": ci(2),
                        "margin_scale_b_over_a": scale, "margin_scale_ci90": ci(3)})
        pairs[f"{a}|{b}"] = rec

    # Multi-source grouping on 1 - tau (unassessable pairs get distance 1)
    n = len(SOURCES)
    D = np.ones((n, n))
    np.fill_diagonal(D, 0.0)
    for i, j in combinations(range(n), 2):
        r = pairs[f"{SOURCES[i]}|{SOURCES[j]}"]
        if "tau" in r:
            D[i, j] = D[j, i] = 1.0 - r["tau"]
    cond = D[np.triu_indices(n, 1)]
    Z = linkage(cond, method="average")
    lab = fcluster(Z, t=args.group_cut, criterion="distance")
    groups = []
    for g in sorted(set(lab)):
        mem = [SOURCES[i] for i in range(n) if lab[i] == g]
        core = set.intersection(*[set(V[s]) for s in mem])
        taus = [pairs.get(f"{p}|{q}", pairs.get(f"{q}|{p}")).get("tau")
                for p, q in combinations(mem, 2)]
        taus = [t for t in taus if t is not None]
        groups.append({"members": mem, "n_core": len(core),
                       "n_union": len(set.union(*[set(V[s]) for s in mem])),
                       "min_pair_tau": min(taus) if taus else None})
    all_core = set.intersection(*[set(V[s]) for s in SOURCES])

    # Cut sweep: every distinct grouping the dendrogram allows, with the
    # range of cut heights (on 1 - tau) that produces it.
    heights = sorted(set(round(float(h), 6) for h in Z[:, 2]))
    edges = [0.0] + heights + [max(1.0, heights[-1] + 1e-3)]
    sweep = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        lab_c = fcluster(Z, t=(lo + hi) / 2, criterion="distance")
        sweep.append({"tau_from": 1.0 - hi, "tau_to": 1.0 - lo,
                      "groups": [[SOURCES[i] for i in range(n) if lab_c[i] == g]
                                 for g in sorted(set(lab_c)) if (lab_c == g).sum() > 1]})

    res = {"config": vars(args), "universe_sha256_16": sha,
           "n_compounds": {s: len(V[s]) for s in SOURCES},
           "pairs": pairs, "linkage": Z.tolist(), "order": SOURCES,
           "groups": groups, "cut_sweep": sweep, "n_core_all_sources": len(all_core)}
    out = (root / args.output_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    (out / "source_combinability.json").write_text(json.dumps(res, indent=1), encoding="utf-8")

    L = ["# Multi-source set operations: which data sources can be combined", "",
         f"Universe sha256[:16] = {sha}; tolerance {args.tol} eV; "
         f"{args.n_boot} bootstrap resamples; pairs with fewer than "
         f"{args.min_shared} shared compounds are not assessable.", "",
         "| A | B | shared | coverage | scissor (eV) | values within tol | Kendall tau | margin scale B/A |",
         "|---|---|---|---|---|---|---|---|"]
    for r in pairs.values():
        if "tau" in r:
            L.append(f"| {LABEL[r['a']]} | {LABEL[r['b']]} | {r['n_shared']} | "
                     f"{r['coverage']:.3f} | {r['scissor_shift_eV']:+.2f} | "
                     f"{r['value_agree_frac']:.2f} [{r['value_agree_ci90'][0]:.2f}, {r['value_agree_ci90'][1]:.2f}] | "
                     f"{r['tau']:+.2f} [{r['tau_ci90'][0]:+.2f}, {r['tau_ci90'][1]:+.2f}] | "
                     f"{r['margin_scale_b_over_a']:.2f} |")
        else:
            L.append(f"| {LABEL[r['a']]} | {LABEL[r['b']]} | {r['n_shared']} | "
                     f"{r['coverage']:.3f} | n/a | n/a | n/a | n/a |")
    L += ["", f"Groups at 1 - tau <= {args.group_cut} (average linkage):", ""]
    for g in groups:
        mt = "n/a" if g["min_pair_tau"] is None else f"{g['min_pair_tau']:+.2f}"
        L.append(f"- {', '.join(LABEL[s] for s in g['members'])}: common compounds "
                 f"{g['n_core']}, union {g['n_union']}, min pairwise tau {mt}")
    L += ["", "Cut sweep (each row: the tau range of the cut that gives this "
          "grouping; singletons omitted):", "",
          "| cut tau range | groups with two or more sources |", "|---|---|"]
    for c in sweep:
        gs = "; ".join("{" + ", ".join(LABEL[s] for s in g) + "}" for g in c["groups"]) or "none"
        L.append(f"| {c['tau_from']:+.3f} to {c['tau_to']:+.3f} | {gs} |")
    L += ["", f"Compounds reported by all eight sources: {len(all_core)}"]
    (out / "source_combinability.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))


if __name__ == "__main__":
    main()
