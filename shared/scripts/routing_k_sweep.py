"""
Routing rules across query arity k (k = 2, 3, 4, 5, 6)

learned_transfer_weights.py evaluates the routing rules (maximin linear
program, uniform, routing learned from the linear-program cold start, oracle)
on the measured early transfer gain of the paired campaigns, for the pairwise
arm (k = 2) and one k-way arm. Here the same evaluation is repeated for every
k-way arm that has a paired-campaign run. All runs use matched seeds, so the
pairwise arm is identical across runs; this is checked, and k = 2 is reported
once.

Components sigma, rho, omega are computed once (they do not depend on k); only
the measured target changes with k. Evaluation is leave-one-receiver-out.

Usage (from shared/):
    uv run python scripts/routing_k_sweep.py --universe <jsonl> \
        --paired 3=<dir>/paired_campaign.json 4=... 5=... 6=...
Outputs routing_k_sweep.json / .md / .png in --output-dir. ASCII only.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import argparse
import hashlib
import importlib.util
import json

import numpy as np


def _load(name):
    path = Path(__file__).resolve().parent / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ltw = _load("learned_transfer_weights")
RULES = [("lp_fixed", "maximin linear program"),
         ("uniform", "uniform over donors"),
         ("learned_from_lp_fixed_score", "learned from LP, default score"),
         ("learned_from_lp_learned_score", "learned from LP, learned score"),
         ("oracle", "best measured donor")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--universe", default="output/bandgap_archive_build/bandgap_universe.jsonl")
    ap.add_argument("--paired", nargs="+", required=True,
                    help="k=path pairs to paired_campaign.json files")
    ap.add_argument("--n-events", type=int, default=300)
    ap.add_argument("--min-class-size", type=int, default=30)
    ap.add_argument("--n-boot", type=int, default=30)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--output-dir", default="output/routing_k_sweep")
    args = ap.parse_args()
    root = Path(__file__).resolve().parent.parent
    universe = (root / args.universe).resolve()
    sha = hashlib.sha256(universe.read_bytes()).hexdigest()[:16]
    out = (root / args.output_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)

    spaces, _info = ltw.build_spaces(universe, args.n_events, args.min_class_size,
                                     args.n_boot, args.seed)
    names = sorted(spaces)
    comp = ltw.components(spaces)

    runs = {}
    for item in args.paired:
        k, p = item.split("=", 1)
        runs[int(k)] = json.load(open(Path(p).expanduser(), encoding="utf-8"))
    ks = sorted(runs)
    pair_ref = runs[ks[0]]["transfer_matrices"]["full"]["pair"]
    same = all(abs(runs[k]["transfer_matrices"]["full"]["pair"][a][b]["T_early"]
                   - pair_ref[a][b]["T_early"]) < 1e-12
               for k in ks for a in pair_ref for b in pair_ref[a])
    n_seeds = {k: runs[k]["config"].get("n_seeds") for k in ks}

    res = {"universe_sha256_16": sha, "n_seeds": n_seeds,
           "pairwise_arm_identical_across_runs": bool(same), "by_k": {}}
    arms = [(2, pair_ref)] + [(k, runs[k]["transfer_matrices"]["full"]["kway"]) for k in ks]
    for k, T in arms:
        r = ltw.evaluate("pair" if k == 2 else f"k{k}", comp, T, names)
        # mean measured early gain of every donor cell, for scale
        allT = [T[a][b]["T_early"] for a in T for b in T[a]]
        res["by_k"][k] = {"routing_mean": r["routing_mean"],
                          "mean_cell_gain": float(np.mean(allT)),
                          "spearman_fixed": r["mean_within_receiver_spearman_fixed"],
                          "spearman_learned": r["mean_within_receiver_spearman_learned"],
                          "top_donor_hits_fixed": r["top_donor_hits_fixed"],
                          "per_receiver": {rc: {"best_donor": v["best_donor_measured"],
                                                "top_donor_fixed": v["top_donor_fixed"],
                                                "eta_fixed_score": v["eta_fixed_score"]}
                                           for rc, v in r["per_receiver"].items()}}
        print(f"[k={k}] " + ", ".join(f"{a} {r['routing_mean'][a]:+.3f}" for a, _ in RULES))
    (out / "routing_k_sweep.json").write_text(json.dumps(res, indent=2))

    L = ["# Routing rules across query arity k", "",
         f"Archive sha256[:16] = `{sha}`; seeds per run {n_seeds}; pairwise arm identical "
         f"across runs: {same}. Values: mean measured early transfer gain on held-out "
         "receivers (leave one receiver out).", "",
         "| k | " + " | ".join(lbl for _, lbl in RULES) + " | learned minus LP |",
         "|---|" + "---|" * (len(RULES) + 1)]
    for k in sorted(res["by_k"]):
        rm = res["by_k"][k]["routing_mean"]
        L.append(f"| {k} | " + " | ".join(f"{rm[a]:+.3f}" for a, _ in RULES)
                 + f" | {rm['learned_from_lp_fixed_score'] - rm['lp_fixed']:+.3f} |")
    L += ["", "Best measured donor per receiver:", ""]
    for k in sorted(res["by_k"]):
        pr = res["by_k"][k]["per_receiver"]
        L.append(f"- k={k}: " + ", ".join(f"{rc}: {v['best_donor']}" for rc, v in sorted(pr.items())))
    (out / "routing_k_sweep.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
