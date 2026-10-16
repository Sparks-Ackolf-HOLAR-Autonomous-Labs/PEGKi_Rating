"""
Four checks of the PEGKi transfer score on the band-gap archive

The design choices of the transfer score Delta(A->B) = sigma + rho - omega
and of the maximin routing are tested directly on materials data:

  1. Weight robustness: each weight is swept over {0.25, 0.5, 1, 2, 3} with
     the other two at 1; the donor ordering of every receiver is compared
     with the (1, 1, 1) ordering (Spearman).
  2. The overhead term omega: the score and the maximin routing are
     recomputed without omega; the allocation and the measured early gain of
     the routing are compared, together with the agreement of the score with
     the measured gain.
  3. Agreement versus transfer: a consensus-medoid rule picks, for each
     receiver, the donor whose class ordering agrees most with the other
     sources (receiver excluded). Its measured gain is compared with the
     donor ranked first by Delta and with the best measured donor.
  4. Noise: Gaussian noise of 0.1, 0.25 and 0.5 eV is added to every
     archived gap, the class ratings and the Delta matrix are rebuilt, and
     the donor ordering (mean outgoing Delta) and each source's class
     ordering are compared with the noise-free ones.
  5. Whole rounds versus split pairs: within each source, the Plackett-Luce
     fit on whole k-way rounds is compared (Kendall tau of the class
     orderings) with the fit on the same rounds split into all pairs.

Measured gains are the early transfer gains T_early of the 50-seed paired
campaigns (TENKi exp 26, CHPC stage 2, fixed-trust full reuse), in the
pairwise and the 4-way query arms. The class-level knowledge spaces are the
ones of scripts/learned_transfer_weights.py.

Usage (from shared/):
    uv run python scripts/transfer_score_checks.py [--universe PATH]
Outputs: transfer_score_checks.json / .md in --output-dir. ASCII only.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import argparse
import hashlib
import importlib.util
import json

import numpy as np
from scipy.stats import spearmanr, kendalltau

from src.rating.hypergraph_rating import PlackettLuceModel, Hyperedge  # noqa: E402
from src.evaluation.phase3.knowledge_spaces import (  # noqa: E402
    KnowledgeSpace, KnowledgeSpaceType)


def _load(name):
    path = Path(__file__).resolve().parent / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ltw = _load("learned_transfer_weights")
trial = ltw.trial
LAMBDAS = (0.25, 0.5, 1.0, 2.0, 3.0)
NOISE_EV = (0.1, 0.25, 0.5)
COMP_NAMES = ("sigma", "rho", "omega")


def spaces_from(per_channel, n_events, min_class_size, n_boot, seed):
    ch_edges = trial.within_channel_edges(per_channel, n_events, min_class_size)
    spaces, theta_all = {}, {}
    for ch in sorted(ch_edges):
        edges = ch_edges[ch]
        theta = dict(PlackettLuceModel().fit(edges).log_strengths)
        rng = np.random.default_rng([seed, len(ch)])
        reps = []
        for _ in range(n_boot):
            samp = [edges[j] for j in rng.integers(0, len(edges), len(edges))]
            reps.append(PlackettLuceModel().fit(samp).log_strengths)
        sd = {c: float(np.std([r.get(c, np.nan) for r in reps])) for c in theta}
        cl = per_channel[ch]
        maes = {c: float(np.mean(np.abs(np.asarray(cl[c]) - np.median(cl[c]))))
                for c in theta}
        ranks = {c: r + 1 for r, c in enumerate(sorted(theta, key=lambda k: -theta[k]))}
        spaces[ch] = KnowledgeSpace(
            space_type=KnowledgeSpaceType.K_H, study_type=ch,
            policy_mus=dict(theta), policy_sigmas=sd, policy_ranks=ranks,
            policy_maes=maes, policy_mae_stds={}, policy_stat_mus={},
            policy_bc_mus={},
            policy_s_transitive={c: {"bandgap": v} for c, v in theta.items()})
        theta_all[ch] = theta
    return spaces, theta_all


def delta_matrix(comp, names, w):
    idx = {n: i for i, n in enumerate(names)}
    D = np.zeros((len(names), len(names)))
    for (a, b), x in comp.items():
        D[idx[a], idx[b]] = float(x @ w)
    return D


def class_tau(t1, t2):
    common = sorted(set(t1) & set(t2))
    if len(common) < 3:
        return np.nan
    return float(kendalltau([t1[c] for c in common], [t2[c] for c in common]).correlation)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--universe", default="output/bandgap_archive_build/bandgap_universe.jsonl")
    ap.add_argument("--paired", default="output/bandgap_cluster_run2/paired_campaign_k4/paired_campaign.json")
    ap.add_argument("--n-events", type=int, default=300)
    ap.add_argument("--min-class-size", type=int, default=30)
    ap.add_argument("--n-boot", type=int, default=30)
    ap.add_argument("--noise-boot", type=int, default=10)
    ap.add_argument("--noise-reps", type=int, default=3)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--output-dir", default="output/transfer_score_checks")
    args = ap.parse_args()
    root = Path(__file__).resolve().parent.parent
    universe = (root / args.universe).resolve()
    sha = hashlib.sha256(universe.read_bytes()).hexdigest()[:16]
    out = (root / args.output_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)

    per_channel, _f, _s, _u = trial.load_universe(universe)
    spaces, theta = spaces_from(per_channel, args.n_events, args.min_class_size,
                                args.n_boot, args.seed)
    names = sorted(spaces)
    comp = ltw.components(spaces)
    paired = json.load(open(root / args.paired, encoding="utf-8"))
    res = {"config": vars(args), "universe_sha256_16": sha, "sources": names}
    W0 = np.array([1.0, 1.0, 1.0])

    # 1. weight robustness ------------------------------------------------
    sweep = {}
    for k, cname in enumerate(COMP_NAMES):
        rows = []
        for lam in LAMBDAS:
            w = W0.copy()
            w[k] = lam
            rhos = []
            for r in names:
                don = [a for a in names if a != r]
                d0 = [comp[(a, r)] @ W0 for a in don]
                d1 = [comp[(a, r)] @ w for a in don]
                rhos.append(spearmanr(d0, d1).correlation)
            rows.append({"lambda": lam, "min_receiver_spearman": float(np.min(rhos)),
                         "mean_receiver_spearman": float(np.mean(rhos))})
        sweep[cname] = rows
    res["weight_sweep"] = sweep

    # 2-3. omega ablation and agreement versus transfer, per query arm -----
    res["arms"] = {}
    for arm in ("pair", "kway"):
        T = paired["transfer_matrices"]["full"][arm]
        arm_res = {}
        for label, w in (("with_omega", W0), ("without_omega", np.array([1.0, 1.0, 0.0]))):
            D = delta_matrix(comp, names, w)
            p = ltw.lp_strategy(names, D)
            vals, rhos = [], []
            for j, r in enumerate(names):
                dn = [i for i, a in enumerate(names) if a != r]
                Ty = np.array([T[names[i]][r]["T_early"] for i in dn])
                vals.append(ltw.route_value(p[dn], Ty))
                rhos.append(spearmanr(D[dn, j], Ty).correlation)
            pn = p / p.sum()
            ent = float(-(pn[pn > 0] * np.log(pn[pn > 0])).sum() / np.log(len(pn)))
            arm_res[label] = {"lp_strategy": dict(zip(names, pn.tolist())),
                              "lp_max_weight": float(pn.max()),
                              "lp_normalised_entropy": ent,
                              "routing_value": float(np.mean(vals)),
                              "mean_receiver_spearman_with_T": float(np.nanmean(rhos))}
        uni = []
        med, top, orc, agree_rho = [], [], [], []
        per_r = {}
        for r in names:
            don = [a for a in names if a != r]
            Ty = np.array([T[a][r]["T_early"] for a in don])
            uni.append(float(Ty.mean()))
            agree = np.array([np.nanmean([class_tau(theta[a], theta[s])
                                          for s in names if s not in (a, r)]) for a in don])
            dfix = np.array([comp[(a, r)] @ W0 for a in don])
            med.append(float(Ty[int(np.argmax(agree))]))
            top.append(float(Ty[int(np.argmax(dfix))]))
            orc.append(float(Ty.max()))
            agree_rho.append(spearmanr(agree, Ty).correlation)
            per_r[r] = {"medoid_donor": don[int(np.argmax(agree))],
                        "delta_top_donor": don[int(np.argmax(dfix))],
                        "best_measured_donor": don[int(np.argmax(Ty))]}
        arm_res["agreement_vs_transfer"] = {
            "uniform": float(np.mean(uni)), "consensus_medoid": float(np.mean(med)),
            "delta_top_donor": float(np.mean(top)), "best_measured": float(np.mean(orc)),
            "mean_receiver_spearman_agreement_with_T": float(np.nanmean(agree_rho)),
            "per_receiver": per_r}
        res["arms"][arm] = arm_res

    # 4. noise ---------------------------------------------------------------
    base_sp, base_th = spaces_from(per_channel, args.n_events, args.min_class_size,
                                   args.noise_boot, args.seed)
    base_comp = ltw.components(base_sp)
    base_out = np.array([np.mean([base_comp[(a, b)] @ W0 for b in names if b != a])
                         for a in names])
    noise = {}
    for s in NOISE_EV:
        donor_rho, class_t = [], {n: [] for n in names}
        for rep in range(args.noise_reps):
            rng = np.random.default_rng([args.seed, int(s * 1000), rep])
            noisy = {}
            for ch in sorted(per_channel):
                noisy[ch] = {}
                for cl in sorted(per_channel[ch]):
                    v = np.asarray(per_channel[ch][cl]) + rng.normal(0.0, s, len(per_channel[ch][cl]))
                    noisy[ch][cl] = sorted(v.tolist())
            sp, th = spaces_from(noisy, args.n_events, args.min_class_size,
                                 args.noise_boot, args.seed)
            cp = ltw.components(sp)
            outg = np.array([np.mean([cp[(a, b)] @ W0 for b in names if b != a]) for a in names])
            donor_rho.append(float(spearmanr(base_out, outg).correlation))
            for n in names:
                class_t[n].append(class_tau(base_th[n], th[n]))
        noise[str(s)] = {"donor_ordering_spearman_mean": float(np.mean(donor_rho)),
                         "donor_ordering_spearman_min": float(np.min(donor_rho)),
                         "class_ordering_tau_mean": {n: float(np.mean(v)) for n, v in class_t.items()}}
        print(f"[noise {s}] donor rho {np.mean(donor_rho):+.3f}")
    res["noise"] = {"noise_free_mean_outgoing_delta": dict(zip(names, base_out.tolist())),
                    "levels": noise}

    # 5. whole rounds versus split pairs ------------------------------------
    ch_edges = trial.within_channel_edges(per_channel, args.n_events, args.min_class_size)
    wr = {}
    for ch in sorted(ch_edges):
        edges = ch_edges[ch]
        whole = PlackettLuceModel().fit(edges).log_strengths
        pairs = [Hyperedge(ranking=(e.ranking[i], e.ranking[j]), channel=e.channel)
                 for e in edges for i in range(len(e.ranking))
                 for j in range(i + 1, len(e.ranking))]
        split = PlackettLuceModel().fit(pairs).log_strengths
        wr[ch] = class_tau(whole, split)
    res["whole_vs_split_tau"] = wr

    (out / "transfer_score_checks.json").write_text(json.dumps(res, indent=2))

    L = ["# Transfer-score checks on the band-gap archive", "",
         f"Archive sha256[:16] = `{sha}`. Measured gains: 50-seed paired campaigns "
         "(TENKi exp 26, CHPC, fixed-trust full reuse).", "",
         "## 1. Weight sweep (donor ordering per receiver versus the (1,1,1) ordering)", "",
         "| component | lambda | min Spearman | mean Spearman |", "|---|---|---|---|"]
    for c, rows in sweep.items():
        for row in rows:
            L.append(f"| {c} | {row['lambda']} | {row['min_receiver_spearman']:+.3f} | {row['mean_receiver_spearman']:+.3f} |")
    for arm, a in res["arms"].items():
        L += ["", f"## 2. Overhead term, {arm} arm", "",
              "| score | LP max weight | LP entropy (1 = uniform) | routing value | mean Spearman with gain |",
              "|---|---|---|---|---|"]
        for lab in ("with_omega", "without_omega"):
            x = a[lab]
            L.append(f"| {lab} | {x['lp_max_weight']:.3f} | {x['lp_normalised_entropy']:.3f} | "
                     f"{x['routing_value']:+.4f} | {x['mean_receiver_spearman_with_T']:+.3f} |")
        g = a["agreement_vs_transfer"]
        L += ["", f"## 3. Agreement versus transfer, {arm} arm", "",
              "| rule | mean measured early gain |", "|---|---|",
              f"| uniform over donors | {g['uniform']:+.4f} |",
              f"| consensus medoid | {g['consensus_medoid']:+.4f} |",
              f"| donor ranked first by Delta | {g['delta_top_donor']:+.4f} |",
              f"| best measured donor | {g['best_measured']:+.4f} |", "",
              f"Mean within-receiver Spearman of agreement with measured gain: "
              f"{g['mean_receiver_spearman_agreement_with_T']:+.3f}"]
    L += ["", "## 4. Noise", "",
          "| noise (eV) | donor-ordering Spearman mean | min | " + " | ".join(names) + " |",
          "|---|---|---|" + "---|" * len(names)]
    for s, x in noise.items():
        L.append(f"| {s} | {x['donor_ordering_spearman_mean']:+.3f} | {x['donor_ordering_spearman_min']:+.3f} | "
                 + " | ".join(f"{x['class_ordering_tau_mean'][n]:+.2f}" for n in names) + " |")
    L += ["", "## 5. Whole rounds versus split pairs (Kendall tau of class orderings)", "",
          "| " + " | ".join(names) + " |", "|" + "---|" * len(names),
          "| " + " | ".join(f"{wr.get(n, float('nan')):+.3f}" for n in names) + " |"]
    (out / "transfer_score_checks.md").write_text("\n".join(L) + "\n")
    print(f"[done] {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
