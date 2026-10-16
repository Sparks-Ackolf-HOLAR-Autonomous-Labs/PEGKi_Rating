"""
Learned transfer-score weights and learned routing, with the fixed choices
kept as the cold start

Two fixed choices in the PEGKi transfer layer are the ones reviewers question
most:

  1. The transfer score Delta(A->B) = alpha*sigma + gamma*rho - delta*omega
     uses fixed weights alpha = gamma = delta = 1
     (src/evaluation/phase3/transfer_equation.py).
  2. Routing uses the Nash maximin linear program
     (src/evaluation/phase3/game_theory.py), which assumes the receiver is an
     adversary and returns one mixed strategy for all receivers.

Neither is a law of nature. Here both are kept as the DEFAULT (cold start)
and are replaced by versions learned from measured campaign outcomes whenever
those outcomes exist:

  * Components sigma, rho, omega are computed with the committed
    transfer_equation functions on class-level knowledge spaces built from the
    version-locked band-gap archive: each source's material classes are rated
    by Plackett-Luce on that source's within-source ranking events (the same
    events as scripts/hypergraph_trial_bandgap.py); rating uncertainty is the
    event-bootstrap SD; the per-class error is the L1 dispersion around the
    class median (the class-as-policy summary of materials_summaries.py).
  * The target is the measured early transfer gain T_early(A->B) of the
    50-seed paired campaign protocol (TENKi exp 26, CHPC stage 2, fixed-trust
    full reuse), available for all 56 ordered source pairs, in both query
    arms (pairwise and 4-way).
  * Learned weights: ridge regression of T_early on (sigma, rho, -omega) with
    the penalty centred on the calibrated cold start c*(1, 1, 1):
        min sum (T - b - x.v)^2 + lam * ||v - c*(1,1,1)||^2
    lam = infinity is exactly the fixed-weight score (after a linear
    calibration b + c*Delta that does not change any ordering); lam is chosen
    by inner leave-one-receiver-out validation.
  * Learned routing with the LP as cold start:
        p(A | B) = (1 - eta) * p_LP(A) + eta * softmax_A(Delta(A->B) / tau)
    eta = 0 is exactly the LP. (eta, tau) are chosen on training receivers,
    on the fixed score and on the learned score; compared with the LP alone,
    uniform routing, and the oracle (best measured donor).
  * Evaluation is leave-one-receiver-out: every reported number for receiver
    B comes from weights, lam and tau fitted without any cell whose receiver
    is B.

Usage (from shared/):
    uv run python scripts/learned_transfer_weights.py
Outputs: learned_transfer_weights.json / .md in --output-dir. ASCII only.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import argparse
import hashlib
import importlib.util
import json

import numpy as np
from scipy.stats import spearmanr

from src.rating.hypergraph_rating import PlackettLuceModel  # noqa: E402
from src.evaluation.phase3.knowledge_spaces import (  # noqa: E402
    KnowledgeSpace, KnowledgeSpaceType)
from src.evaluation.phase3.transfer_equation import (  # noqa: E402
    compute_structural_compatibility, compute_mapping_efficiency,
    compute_transfer_overhead)
from src.evaluation.phase3.game_theory import compute_nash_routing  # noqa: E402


def _load_trial():
    path = Path(__file__).resolve().parent / "hypergraph_trial_bandgap.py"
    spec = importlib.util.spec_from_file_location("hypergraph_trial_bandgap", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


trial = _load_trial()
LAM_GRID = (0.0, 0.01, 0.1, 1.0, 10.0, 100.0, float("inf"))
TAU_GRID = (0.01, 0.03, 0.1, 0.3, 1.0, 3.0)
ETA_GRID = (0.0, 0.25, 0.5, 0.75, 1.0)
W0 = np.array([1.0, 1.0, 1.0])  # (alpha, gamma, delta) on (sigma, rho, -omega)


def build_spaces(universe, n_events, min_class_size, n_boot, seed):
    per_channel, _formulas, _n_seen, _n_used = trial.load_universe(universe)
    ch_edges = trial.within_channel_edges(per_channel, n_events, min_class_size)
    spaces, info = {}, {}
    for ch in sorted(ch_edges):
        edges = ch_edges[ch]
        m = PlackettLuceModel().fit(edges)
        theta = dict(m.log_strengths)
        rng = np.random.default_rng([seed, len(ch)])
        reps = []
        for _ in range(n_boot):
            samp = [edges[j] for j in rng.integers(0, len(edges), len(edges))]
            reps.append(PlackettLuceModel().fit(samp).log_strengths)
        sd = {c: float(np.std([r.get(c, np.nan) for r in reps])) for c in theta}
        classes = per_channel[ch]
        maes = {c: float(np.mean(np.abs(np.asarray(classes[c]) - np.median(classes[c]))))
                for c in theta}
        ranks = {c: r + 1 for r, c in enumerate(sorted(theta, key=lambda k: -theta[k]))}
        spaces[ch] = KnowledgeSpace(
            space_type=KnowledgeSpaceType.K_H, study_type=ch,
            policy_mus=dict(theta), policy_sigmas=sd, policy_ranks=ranks,
            policy_maes=maes, policy_mae_stds={}, policy_stat_mus={},
            policy_bc_mus={},
            policy_s_transitive={c: {"bandgap": v} for c, v in theta.items()})
        info[ch] = {"classes": sorted(theta), "theta": theta, "theta_sd": sd,
                    "class_l1_dispersion": maes, "n_events": len(edges)}
    return spaces, info


def components(spaces):
    comp = {}
    for a in spaces:
        for b in spaces:
            if a != b:
                comp[(a, b)] = np.array([
                    compute_structural_compatibility(spaces[a], spaces[b]),
                    compute_mapping_efficiency(spaces[a], spaces[b]),
                    -compute_transfer_overhead(spaces[a], spaces[b])])
    return comp


def fit_weights(X, y, lam):
    """Ridge toward the calibrated cold start; returns (b, v)."""
    d0 = X @ W0
    A = np.column_stack([np.ones_like(d0), d0])
    b0, c0 = np.linalg.lstsq(A, y, rcond=None)[0]
    if np.isinf(lam):
        return float(b0), c0 * W0
    target = c0 * W0
    # centre for the intercept, then ridge on v around target
    xm, ym = X.mean(axis=0), y.mean()
    Xc, yc = X - xm, y - ym
    v = np.linalg.solve(Xc.T @ Xc + lam * np.eye(3), Xc.T @ yc + lam * target)
    b = float(ym - xm @ v)
    return b, v


def lor_cv_lam(cells, X, y, receivers):
    """Inner leave-one-receiver-out validation error per lam."""
    errs = {}
    for lam in LAM_GRID:
        se = []
        for r in receivers:
            tr = np.array([c[1] != r for c in cells])
            b, v = fit_weights(X[tr], y[tr], lam)
            pred = b + X[~tr] @ v
            se.extend((pred - y[~tr]) ** 2)
        errs[lam] = float(np.mean(se))
    best = min(errs, key=errs.get)
    return best, errs


def lp_strategy(names, D):
    info = {"study_types": names, "matrix": D.tolist()}
    s = compute_nash_routing(info)["strategy"]
    return np.array([s.get(n, 0.0) for n in names])


def route_value(p, donors_T):
    p = np.maximum(p, 0)
    if p.sum() <= 0:
        p = np.ones_like(p)
    p = p / p.sum()
    return float(p @ donors_T)


def softmax_route(scores, tau):
    z = scores / tau
    z = z - z.max()
    e = np.exp(z)
    return e / e.sum()


def blend(p_lp, scores, eta, tau):
    p_lp = np.maximum(np.asarray(p_lp, float), 0)
    p_lp = p_lp / p_lp.sum() if p_lp.sum() > 0 else np.ones_like(p_lp) / len(p_lp)
    return (1 - eta) * p_lp + eta * softmax_route(np.asarray(scores, float), tau)


def evaluate(arm, comp, T, names):
    cells = [(a, b) for (a, b) in comp if b in T.get(a, {})]
    X = np.array([comp[c] for c in cells])
    y = np.array([T[c[0]][c[1]]["T_early"] for c in cells])
    receivers = sorted({c[1] for c in cells})
    idx = {n: i for i, n in enumerate(names)}
    out = {"per_receiver": {}, "arm": arm}
    pred_cold_all, pred_learn_all, y_all, rec_all = [], [], [], []
    route = {k: [] for k in ("lp_fixed", "lp_learned",
                             "learned_from_lp_fixed_score",
                             "learned_from_lp_learned_score",
                             "uniform", "oracle")}
    lam_chosen, tau_chosen = {}, {}
    for r in receivers:
        tr = np.array([c[1] != r for c in cells])
        cells_tr = [c for c, t in zip(cells, tr) if t]
        rec_tr = [x for x in receivers if x != r]
        lam, _ = lor_cv_lam(cells_tr, X[tr], y[tr], rec_tr)
        b, v = fit_weights(X[tr], y[tr], lam)
        bc, vc = fit_weights(X[tr], y[tr], float("inf"))
        lam_chosen[r] = lam
        # learned routing with the LP as cold start:
        #   p = (1 - eta) * p_LP + eta * softmax(score / tau)
        # eta = 0 is exactly the LP; (eta, tau) chosen on training receivers
        n = len(names)
        Dfix, Dlearn = np.zeros((n, n)), np.zeros((n, n))
        for (a, bb), x in comp.items():
            Dfix[idx[a], idx[bb]] = x @ W0
            Dlearn[idx[a], idx[bb]] = b + x @ v
        p_lp_fix, p_lp_learn = lp_strategy(names, Dfix), lp_strategy(names, Dlearn)

        def choose(score_fn, p_lp):
            best = (0.0, TAU_GRID[0], -np.inf)
            for eta in ETA_GRID:
                for tau in TAU_GRID:
                    vals = []
                    for r2 in rec_tr:
                        don2 = [c for c in cells_tr if c[1] == r2]
                        sc = np.array([score_fn(c) for c in don2])
                        p = blend(p_lp[[idx[c[0]] for c in don2]], sc, eta, tau)
                        vals.append(route_value(p, np.array(
                            [T[c[0]][r2]["T_early"] for c in don2])))
                    if np.mean(vals) > best[2] + 1e-12:
                        best = (eta, tau, float(np.mean(vals)))
            return best

        eta_l, best_tau, _ = choose(lambda c: b + comp[c] @ v, p_lp_learn)
        eta_c, best_tau_c, _ = choose(lambda c: comp[c] @ W0, p_lp_fix)
        tau_chosen[r] = best_tau
        don = [c for c in cells if c[1] == r]
        Ty = np.array([T[c[0]][r]["T_early"] for c in don])
        d_cold = np.array([comp[c] @ W0 for c in don])
        d_learn = np.array([b + comp[c] @ v for c in don])
        pred_cold_all.extend(bc + np.array([comp[c] @ vc for c in don]))
        pred_learn_all.extend(d_learn)
        y_all.extend(Ty)
        rec_all.extend([r] * len(don))
        dn = [idx[c[0]] for c in don]
        route["lp_fixed"].append(route_value(p_lp_fix[dn], Ty))
        route["lp_learned"].append(route_value(p_lp_learn[dn], Ty))
        route["learned_from_lp_learned_score"].append(route_value(
            blend(p_lp_learn[dn], d_learn, eta_l, best_tau), Ty))
        route["learned_from_lp_fixed_score"].append(route_value(
            blend(p_lp_fix[dn], d_cold, eta_c, best_tau_c), Ty))
        route["uniform"].append(float(Ty.mean()))
        route["oracle"].append(float(Ty.max()))
        rho_c = spearmanr(d_cold, Ty).correlation
        rho_l = spearmanr(d_learn, Ty).correlation
        out["per_receiver"][r] = {
            "lam": lam, "eta_fixed_score": eta_c, "tau_fixed_score": best_tau_c,
            "eta_learned_score": eta_l, "tau_learned_score": best_tau, "weights_learned": (v / np.abs(v).sum()).tolist(),
            "spearman_fixed": float(rho_c), "spearman_learned": float(rho_l),
            "best_donor_measured": don[int(np.argmax(Ty))][0],
            "top_donor_fixed": don[int(np.argmax(d_cold))][0],
            "top_donor_learned": don[int(np.argmax(d_learn))][0],
            "routing": {k: route[k][-1] for k in route}}
    y_all = np.array(y_all)
    pc, pl = np.array(pred_cold_all), np.array(pred_learn_all)
    out["mse_fixed_calibrated"] = float(np.mean((pc - y_all) ** 2))
    out["mse_learned"] = float(np.mean((pl - y_all) ** 2))
    out["sign_acc_fixed_raw"] = float(np.mean((np.array([comp[(c0, c1)] @ W0 for c0, c1 in
                                       [(a, b) for a, b in cells if True]]) > 0) ==
                                       (np.array([T[a][b]["T_early"] for a, b in cells]) > 0)))
    out["sign_acc_fixed_calibrated"] = float(np.mean((pc > 0) == (y_all > 0)))
    out["sign_acc_learned"] = float(np.mean((pl > 0) == (y_all > 0)))
    out["spearman_pooled_fixed"] = float(spearmanr(pc, y_all).correlation)
    out["spearman_pooled_learned"] = float(spearmanr(pl, y_all).correlation)
    out["mean_within_receiver_spearman_fixed"] = float(np.nanmean(
        [v["spearman_fixed"] for v in out["per_receiver"].values()]))
    out["mean_within_receiver_spearman_learned"] = float(np.nanmean(
        [v["spearman_learned"] for v in out["per_receiver"].values()]))
    out["top_donor_hits_fixed"] = int(sum(v["top_donor_fixed"] == v["best_donor_measured"]
                                          for v in out["per_receiver"].values()))
    out["top_donor_hits_learned"] = int(sum(v["top_donor_learned"] == v["best_donor_measured"]
                                            for v in out["per_receiver"].values()))
    out["routing_mean"] = {k: float(np.mean(v)) for k, v in route.items()}
    out["n_receivers"] = len(receivers)
    # weights on all cells, lam by full leave-one-receiver-out
    lam_all, errs = lor_cv_lam(cells, X, y, receivers)
    b, v = fit_weights(X, y, lam_all)
    bc, vc = fit_weights(X, y, float("inf"))
    out["all_cells"] = {"lam": lam_all, "cv_mse_by_lam": {str(k): e for k, e in errs.items()},
                        "intercept": b, "v": v.tolist(),
                        "weights_normalised_alpha_gamma_delta": (v / np.abs(v).sum()).tolist(),
                        "cold_start_calibration": {"b": bc, "c": float(vc[0])}}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--universe", default="output/bandgap_archive_build/bandgap_universe.jsonl")
    ap.add_argument("--paired", default="output/bandgap_cluster_run2/paired_campaign_k4/paired_campaign.json")
    ap.add_argument("--n-events", type=int, default=300)
    ap.add_argument("--min-class-size", type=int, default=30)
    ap.add_argument("--n-boot", type=int, default=30)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--output-dir", default="output/learned_transfer_weights")
    args = ap.parse_args()
    root = Path(__file__).resolve().parent.parent
    universe = (root / args.universe).resolve()
    sha = hashlib.sha256(universe.read_bytes()).hexdigest()[:16]
    out = (root / args.output_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)

    spaces, info = build_spaces(universe, args.n_events, args.min_class_size,
                                args.n_boot, args.seed)
    names = sorted(spaces)
    comp = components(spaces)
    paired = json.load(open(root / args.paired, encoding="utf-8"))
    res = {"config": vars(args), "universe_sha256_16": sha, "sources": names,
           "spaces": info,
           "components": {f"{a}->{b}": {"sigma": float(x[0]), "rho": float(x[1]),
                                        "omega": float(-x[2]), "delta_fixed": float(x @ W0)}
                          for (a, b), x in comp.items()},
           "arms": {}}
    for arm in ("pair", "kway"):
        T = paired["transfer_matrices"]["full"][arm]
        r = evaluate(arm, comp, T, names)
        res["arms"][arm] = r
        print(f"[{arm}] within-receiver Spearman fixed {r['mean_within_receiver_spearman_fixed']:+.3f} "
              f"learned {r['mean_within_receiver_spearman_learned']:+.3f}; "
              f"sign acc raw {r['sign_acc_fixed_raw']:.3f} cal {r['sign_acc_fixed_calibrated']:.3f} "
              f"learned {r['sign_acc_learned']:.3f}; mse {r['mse_fixed_calibrated']:.4f} -> {r['mse_learned']:.4f}")
        print(f"[{arm}] routing " + ", ".join(f"{k} {v:+.4f}" for k, v in r["routing_mean"].items()))
        print(f"[{arm}] all-cell weights {np.round(r['all_cells']['weights_normalised_alpha_gamma_delta'], 3)} lam {r['all_cells']['lam']}")
    (out / "learned_transfer_weights.json").write_text(json.dumps(res, indent=2))

    L = ["# Learned transfer-score weights and routing (fixed choices as cold start)", "",
         f"Archive sha256[:16] = `{sha}`. Target: measured early transfer gain T_early "
         "of the 50-seed paired campaign (TENKi exp 26, CHPC, fixed-trust full reuse). "
         "All numbers are leave-one-receiver-out.", ""]
    for arm, r in res["arms"].items():
        L += [f"## Query arm: {arm}", "",
              "| quantity | fixed weights (1,1,1) | learned weights |", "|---|---|---|",
              f"| mean within-receiver Spearman (donor ordering) | {r['mean_within_receiver_spearman_fixed']:+.3f} | {r['mean_within_receiver_spearman_learned']:+.3f} |",
              f"| pooled Spearman over 56 cells | {r['spearman_pooled_fixed']:+.3f} | {r['spearman_pooled_learned']:+.3f} |",
              f"| sign accuracy (transfer helps or hurts) | {r['sign_acc_fixed_raw']:.3f} raw, {r['sign_acc_fixed_calibrated']:.3f} calibrated | {r['sign_acc_learned']:.3f} |",
              f"| mean squared error of T_early | {r['mse_fixed_calibrated']:.4f} (calibrated) | {r['mse_learned']:.4f} |",
              f"| receivers whose best donor is ranked first | {r['top_donor_hits_fixed']}/{r['n_receivers']} | {r['top_donor_hits_learned']}/{r['n_receivers']} |",
              "", "Routing value (mean measured T_early obtained by the routing rule on held-out receivers):", "",
              "| rule | value |", "|---|---|"]
        for k, v in r["routing_mean"].items():
            L.append(f"| {k} | {v:+.4f} |")
        a = r["all_cells"]
        L += ["", f"Weights fitted on all 56 cells (lam = {a['lam']} by leave-one-receiver-out): "
              f"normalised (alpha, gamma, delta) = ({a['weights_normalised_alpha_gamma_delta'][0]:+.3f}, "
              f"{a['weights_normalised_alpha_gamma_delta'][1]:+.3f}, {a['weights_normalised_alpha_gamma_delta'][2]:+.3f}); "
              "the cold start is (+0.333, +0.333, +0.333).", ""]
    (out / "learned_transfer_weights.md").write_text("\n".join(L) + "\n")
    print(f"[done] {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
