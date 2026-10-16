"""
Pair combinatorics behind the set-operation studies and the hypergraph layer

Both families of results in the PEGKi draft are statements about one finite
object: for a data source (or a set-operation region of it) with n material
classes, the pairwise dominance table

    P_ij = Pr(gap of a random record of class i > gap of a random record of j)

over the C(n,2) unordered class pairs. P is a weighted tournament. This script
computes, from the version-locked archive only, the combinatorial quantities
that the draft's results reduce to, and checks each reduction against the
committed result artifacts.

  (1) Resolvability r = mean over pairs of |2 P_ij - 1| and the near-tie count
      (pairs with |2 P_ij - 1| < 0.1). A comparison of two classes drawn from a
      near-tie pair is close to a coin toss and carries almost no information
      about their order.
  (2) Linearity of expectation over pairs gives an exact identity for one
      full round (one record per class): the expected fraction of pairs that
      the round orders against a consensus order c is
          E[residual] = (1/C(n,2)) sum_{i<j} Pr(round disagrees with c on {i,j})
                      = (1 - mean_{i<j} s^c_ij (2 P_ij - 1)) / 2,
      with s^c_ij = +1 if c puts i above j. Checked against the per-hyperedge
      residuals of shared/output/hypergraph_trial_bandgap (300 full rounds per
      source).
  (3) Directed 3-cycles of the majority tournament (i beats j beats k beats i
      in the P > 1/2 sense, the non-transitive dice pattern) and the Hodge
      baseline: on the complete comparison graph K_n the cycle space has
      dimension C(n-1,2), so isotropic noise puts a fraction (n-2)/n of its
      energy into the cyclic part. The measured cyclic fractions (0.02-0.07)
      are compared with that baseline.
  (4) Order ceiling: pairs where the majority direction of P disagrees with
      the median order used as the reference ranking. No amount of
      comparison data can resolve such a pair toward the reference, so they
      bound the attainable Kendall tau. Checked against the 50-seed k-sweep
      final tau (shared/output/bandgap_cluster_run1/hyperedge_budget_k*).
  (5) Covering designs: the minimum number of k-subsets of n classes that
      covers every class pair at least once, C(n,k,2), bounded below by the
      Schonheim bound; k times that number is the evaluation cost of seeing
      every pair once. Its shrinking returns with k mirror the saturation of
      the k-sweep between k = 4 and k = 6.
  (6) Set-operation regions: r of each agreement, disagreement, and coverage
      region, and the pair-level concordance predictor
          C(A -> X) = mean_{i<j} s^A_ij (2 P^X_ij - 1),
      with s^A the donor's median order, correlated with the committed
      transfer scores Delta(A -> X) of shared/output/set_operation_study.

Every quantity is a count or an average over archived records; nothing is
fitted except where an existing artifact is read back. ASCII only.

Usage (from shared/):
    uv run python scripts/pair_combinatorics.py
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import argparse
import hashlib
import importlib.util
import itertools
import json
import math

import numpy as np


def _load(name):
    path = Path(__file__).resolve().parent / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


sos = _load("set_operation_study")
LABEL = sos.LABEL


def dominance(classes):
    """Exact P_ij from the empirical distributions (ties count one half)."""
    names = sorted(classes)
    arrs = {c: np.sort(np.asarray(classes[c], float)) for c in names}
    P = {}
    for a, b in itertools.combinations(names, 2):
        x, y = arrs[a], arrs[b]
        less = np.searchsorted(y, x, side="left").sum()    # y < x
        leq = np.searchsorted(y, x, side="right").sum()    # y <= x
        p = (less + 0.5 * (leq - less)) / (len(x) * len(y))
        P[(a, b)] = float(p)
        P[(b, a)] = float(1.0 - p)
    return names, P


def median_sign(classes, a, b):
    ma, mb = np.median(classes[a]), np.median(classes[b])
    return 0 if ma == mb else (1 if ma > mb else -1)


def summarise(classes):
    names, P = dominance(classes)
    pairs = list(itertools.combinations(names, 2))
    m = [abs(2 * P[p] - 1) for p in pairs]
    near = [p for p, v in zip(pairs, m) if v < 0.1]
    # majority tournament: i -> j if P_ij > 1/2
    beats = {(a, b) for a in names for b in names if a != b and P[(a, b)] > 0.5}
    cyc = [t for t in itertools.combinations(names, 3)
           if ((t[0], t[1]) in beats and (t[1], t[2]) in beats and (t[2], t[0]) in beats)
           or ((t[1], t[0]) in beats and (t[2], t[1]) in beats and (t[0], t[2]) in beats)]
    disc = [p for p in pairs
            if median_sign(classes, *p) * (1 if P[p] > 0.5 else -1) < 0
            and P[p] != 0.5]
    n = len(names)
    npairs = len(pairs)
    # expected single-round residual against the median (reference) order
    conc_med = np.mean([median_sign(classes, *p) * (2 * P[p] - 1) for p in pairs])
    return {
        "n_classes": n, "n_pairs": npairs,
        "resolvability": float(np.mean(m)),
        "near_tie_pairs": len(near),
        "near_tie_list": [f"{a}~{b}" for a, b in near],
        "majority_3cycles": len(cyc),
        "median_vs_majority_discordant_pairs": len(disc),
        "order_ceiling_tau": float(1 - 2 * len(disc) / npairs) if npairs else None,
        "expected_round_residual_vs_median": float((1 - conc_med) / 2),
        "expected_round_residual_vs_majority": float((1 - np.mean(m)) / 2),
        "hodge_isotropic_cyclic_fraction": float((n - 2) / n) if n else None,
        "median_spread_eV": float(max(np.median(v) for v in classes.values())
                                  - min(np.median(v) for v in classes.values())),
    }


def concordance(donor_classes, recv_classes):
    common = sorted(set(donor_classes) & set(recv_classes))
    if len(common) < 3:
        return None
    _, P = dominance({c: recv_classes[c] for c in common})
    vals = [median_sign(donor_classes, a, b) * (2 * P[(a, b)] - 1)
            for a, b in itertools.combinations(common, 2)]
    return float(np.mean(vals))


def schoenheim(n, k):
    """Schonheim lower bound on the covering number C(n,k,2)."""
    if k >= n:
        return 1
    return math.ceil(n / k * math.ceil((n - 1) / (k - 1)))


def spearman(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    rx = np.argsort(np.argsort(x)).astype(float)
    ry = np.argsort(np.argsort(y)).astype(float)
    return float(np.corrcoef(rx, ry)[0, 1])


def pearson(x, y):
    return float(np.corrcoef(np.asarray(x, float), np.asarray(y, float))[0, 1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--universe", default="output/bandgap_archive_build/bandgap_universe.jsonl")
    ap.add_argument("--output-dir", default="output/pair_combinatorics")
    args = ap.parse_args()
    root = Path(__file__).resolve().parent.parent
    universe = (root / args.universe).resolve()
    sha = hashlib.sha256(universe.read_bytes()).hexdigest()[:16]
    out = (root / args.output_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)

    so = json.loads((root / "output/set_operation_study/set_operation_study.json").read_text())
    cfg = so["config"]
    recs = sos.load_records(universe)
    base, derived, _meta = sos.build_studies(recs, cfg["tol"], cfg["min_class_size"],
                                             cfg["min_derived_class_size"])
    studies = {**base, **derived}
    res = {"universe_sha256_16": sha, "sources": {}, "regions": {}}
    for name in sorted(studies):
        if len(studies[name]) < 3:
            continue
        s = summarise(studies[name])
        s["label"] = sos.pretty(name)
        (res["regions"] if "|" in name else res["sources"])[name] = s

    # (2) round-residual identity against the committed hypergraph trial
    ht = json.loads((root / "output/hypergraph_trial_bandgap/hypergraph_trial_bandgap.json").read_text())
    within = ht.get("channels") or {}
    ident = []
    for src, s in res["sources"].items():
        rec = within.get(src) if isinstance(within, dict) else None
        if not rec:
            continue
        obs = rec.get("mean_edge_residual")
        cyc = rec.get("pairwise_cyclic_energy_fraction")
        ident.append({"source": src, "observed_edge_residual": obs,
                      "predicted_vs_majority": s["expected_round_residual_vs_majority"],
                      "predicted_vs_median": s["expected_round_residual_vs_median"],
                      "observed_cyclic_fraction": cyc,
                      "isotropic_cyclic_fraction": s["hodge_isotropic_cyclic_fraction"],
                      "majority_3cycles": s["majority_3cycles"]})
    res["round_residual_identity"] = ident
    err = [abs(r["observed_edge_residual"] - r["predicted_vs_majority"]) for r in ident]
    res["round_residual_identity_summary"] = {
        "max_abs_error": float(max(err)), "mean_abs_error": float(np.mean(err)),
        "pearson": pearson([r["observed_edge_residual"] for r in ident],
                           [r["predicted_vs_majority"] for r in ident]),
        "sources_with_majority_3cycles": int(sum(r["majority_3cycles"] > 0 for r in ident))}

    # (4) order ceiling and resolvability against the k-sweep
    ks = {}
    for k in (3, 4, 6):
        p = root / f"output/bandgap_cluster_run1/hyperedge_budget_k{k}/hyperedge_budget.json"
        if p.exists():
            d = json.loads(p.read_text())["per_channel_final"]
            ks[k] = d["kway_random"]
            ks[2] = d["pair_random"]
    sweep = []
    for src, s in res["sources"].items():
        row = {"source": src, "resolvability": s["resolvability"],
               "near_tie_pairs": s["near_tie_pairs"],
               "order_ceiling_tau": s["order_ceiling_tau"], "n_classes": s["n_classes"]}
        for k in sorted(ks):
            row[f"final_tau_k{k}"] = ks[k].get(src)
        sweep.append(row)
    res["k_sweep_rows"] = sweep
    corr = {}
    for k in sorted(ks):
        rows = [r for r in sweep if r.get(f"final_tau_k{k}") is not None]
        corr[f"k{k}"] = {
            "spearman_resolvability": spearman([r["resolvability"] for r in rows],
                                               [r[f"final_tau_k{k}"] for r in rows]),
            "spearman_near_ties": spearman([-r["near_tie_pairs"] for r in rows],
                                           [r[f"final_tau_k{k}"] for r in rows]),
            "mean_final_tau": float(np.mean([r[f"final_tau_k{k}"] for r in rows])),
            "mean_gap_to_ceiling": float(np.mean([r["order_ceiling_tau"] - r[f"final_tau_k{k}"]
                                                  for r in rows]))}
    res["k_sweep_correlations"] = corr

    # (5) covering designs for the typical class count
    ns = sorted({s["n_classes"] for s in res["sources"].values()})
    res["covering"] = {str(n): [{"k": k, "schoenheim_rounds": schoenheim(n, k),
                                 "evaluations": k * schoenheim(n, k),
                                 "pair_observations_per_evaluation": (k - 1) / 2}
                                for k in range(2, n + 1)] for n in ns}

    # (6) set-operation regions: H1 resolvability and the concordance predictor
    h1 = []
    for h in so["h1_ownership"]:
        ra = res["regions"].get(h["agree"], {})
        rd = res["regions"].get(h["disagree"], {})
        h1.append({"donor": h["donor"], "receiver": h["receiver"],
                   "delta_difference": h["difference"],
                   "r_agree": ra.get("resolvability"), "r_disagree": rd.get("resolvability"),
                   "near_ties_agree": ra.get("near_tie_pairs"),
                   "near_ties_disagree": rd.get("near_tie_pairs"),
                   "n_classes": ra.get("n_classes")})
    res["h1_resolvability"] = h1
    # rebuild the rated spaces exactly as set_operation_study.py does, so each
    # committed Delta can be split into its sigma, rho, omega components
    ltw = sos.ltw
    spaces = {}
    for n in sorted(studies):
        ks_, _sd = sos.fit_space(n, studies[n], cfg["n_events"], cfg["n_boot"], cfg["seed"])
        if ks_ is not None and len(ks_.policy_mus) >= 3:
            spaces[n] = ks_
    pts = []
    for key, v in so["delta"].items():
        a, b = key.split("->")
        if a not in spaces or b not in spaces:
            continue
        c = concordance(studies[a], studies[b])
        if c is None:
            continue
        sig = ltw.compute_structural_compatibility(spaces[a], spaces[b])
        rho = ltw.compute_mapping_efficiency(spaces[a], spaces[b])
        om = ltw.compute_transfer_overhead(spaces[a], spaces[b])
        pts.append({"donor": a, "receiver": b, "delta": v, "concordance": c,
                    "sigma": sig, "rho": rho, "omega": om,
                    "delta_recomputed_error": abs(sig + rho - om - v),
                    "receiver_resolvability": summarise(studies[b])["resolvability"]})
    comp = {}
    for term in ("sigma", "rho", "omega", "delta"):
        comp[term] = {
            "spearman_vs_concordance": spearman([p["concordance"] for p in pts],
                                                [p[term] for p in pts]),
            "spearman_vs_receiver_resolvability": spearman(
                [p["receiver_resolvability"] for p in pts], [p[term] for p in pts]),
            "mean": float(np.mean([p[term] for p in pts])),
            "sd": float(np.std([p[term] for p in pts]))}
    res["delta_vs_concordance"] = {
        "n_routes": len(pts),
        "max_recompute_error": float(max(p["delta_recomputed_error"] for p in pts)),
        "pearson": pearson([p["concordance"] for p in pts], [p["delta"] for p in pts]),
        "spearman": spearman([p["concordance"] for p in pts], [p["delta"] for p in pts]),
        "spearman_receiver_resolvability": spearman(
            [p["receiver_resolvability"] for p in pts], [p["delta"] for p in pts]),
        "components": comp}
    for h in h1:
        for side in ("agree", "disagree"):
            reg = next(x["agree" if side == "agree" else "disagree"]
                       for x in so["h1_ownership"]
                       if x["donor"] == h["donor"] and x["receiver"] == h["receiver"])
            p_ = next(p for p in pts if p["donor"] == h["donor"] and p["receiver"] == reg)
            for term in ("sigma", "rho", "omega", "concordance"):
                h[f"{term}_{side}"] = p_[term]
    res["delta_vs_concordance_points"] = pts
    (out / "pair_combinatorics.json").write_text(json.dumps(res, indent=2))

    L = ["# Pair combinatorics behind the set-operation and hypergraph results", "",
         f"Archive sha256[:16] = `{sha}`. P_ij = Pr(random record of class i has the larger gap).", "",
         "## Base data sources", "",
         "| source | classes | pairs | resolvability r | near-tie pairs | majority 3-cycles | median-vs-majority discordant pairs | order ceiling tau | median spread (eV) |",
         "|---|---|---|---|---|---|---|---|---|"]
    for src, s in res["sources"].items():
        L.append(f"| {s['label']} | {s['n_classes']} | {s['n_pairs']} | {s['resolvability']:.3f} | "
                 f"{s['near_tie_pairs']} | {s['majority_3cycles']} | "
                 f"{s['median_vs_majority_discordant_pairs']} | {s['order_ceiling_tau']:+.3f} | "
                 f"{s['median_spread_eV']:.2f} |")
    L += ["", "## One full round: predicted versus observed residual", "",
          "| source | observed edge residual | predicted (majority order) | predicted (median order) | observed cyclic fraction | isotropic-noise cyclic fraction | majority 3-cycles |",
          "|---|---|---|---|---|---|---|"]
    for r in ident:
        L.append(f"| {LABEL.get(r['source'], r['source'])} | {r['observed_edge_residual']:.3f} | "
                 f"{r['predicted_vs_majority']:.3f} | {r['predicted_vs_median']:.3f} | "
                 f"{r['observed_cyclic_fraction'] if r['observed_cyclic_fraction'] is None else format(r['observed_cyclic_fraction'], '.3f')} | "
                 f"{r['isotropic_cyclic_fraction']:.3f} | {r['majority_3cycles']} |")
    rs = res["round_residual_identity_summary"]
    L += ["", f"Identity check: max |observed - predicted| = {rs['max_abs_error']:.3f}, "
          f"mean {rs['mean_abs_error']:.3f}, Pearson {rs['pearson']:+.3f}; sources with a "
          f"directed 3-cycle: {rs['sources_with_majority_3cycles']}."]
    L += ["", "## k-sweep (50 seeds, budget 240) against the pair structure", "",
          "| k | mean final tau | mean gap to order ceiling | Spearman(r, tau) over 8 sources | Spearman(-near ties, tau) |",
          "|---|---|---|---|---|"]
    for k, c in corr.items():
        L.append(f"| {k[1:]} | {c['mean_final_tau']:+.3f} | {c['mean_gap_to_ceiling']:+.3f} | "
                 f"{c['spearman_resolvability']:+.2f} | {c['spearman_near_ties']:+.2f} |")
    L += ["", "## Covering designs (Schonheim lower bound)", ""]
    for n, rows in res["covering"].items():
        L.append(f"n = {n} classes: " + "; ".join(
            f"k={r['k']}: {r['schoenheim_rounds']} rounds, {r['evaluations']} evaluations"
            for r in rows))
    L += ["", "## Set-operation regions", "",
          "| donor to receiver | Delta(agree) - Delta(disagree) | r agree | r disagree | near ties agree / disagree | classes |",
          "|---|---|---|---|---|---|"]
    for h in h1:
        L.append(f"| {LABEL[h['donor']]} to {LABEL[h['receiver']]} | {h['delta_difference']:+.3f} | "
                 f"{h['r_agree']:.3f} | {h['r_disagree']:.3f} | {h['near_ties_agree']} / "
                 f"{h['near_ties_disagree']} | {h['n_classes']} |")
    L += ["", "Components into agree / disagree regions:", ""]
    for h in h1:
        L.append(f"- {LABEL[h['donor']]} to {LABEL[h['receiver']]}: " + "; ".join(
            f"{t} {h[t + '_agree']:.3f} / {h[t + '_disagree']:.3f}"
            for t in ("concordance", "sigma", "rho", "omega")))
    dc = res["delta_vs_concordance"]
    L += ["", "| term | Spearman vs pair concordance | Spearman vs receiver resolvability | mean | sd |",
          "|---|---|---|---|---|"]
    for t, c in dc["components"].items():
        L.append(f"| {t} | {c['spearman_vs_concordance']:+.3f} | "
                 f"{c['spearman_vs_receiver_resolvability']:+.3f} | {c['mean']:.3f} | {c['sd']:.3f} |")
    L += ["", f"Pair-concordance predictor versus committed Delta over {dc['n_routes']} routes: "
          f"Pearson {dc['pearson']:+.3f}, Spearman {dc['spearman']:+.3f}; "
          f"receiver resolvability alone, Spearman {dc['spearman_receiver_resolvability']:+.3f}."]
    (out / "pair_combinatorics.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
