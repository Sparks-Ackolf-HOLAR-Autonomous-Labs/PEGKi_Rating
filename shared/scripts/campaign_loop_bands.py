"""
campaign_loop_bands.py
----------------------
Loops seen during a campaign, with their uncertainty.

The earlier analysis (local_cycles.py, level L3) replayed random pairwise
and 4-way campaigns to 2,400 evaluations and reported, per budget, the
share of campaigns whose win table contains a directed 3-cycle. That is
one yes/no line per query format. Learning campaigns run for far more
rounds, and the line should carry its variation so that loose cycles are
visible, not only cleanly closed ones.

This script therefore
  * replays far longer campaigns (to 9,600 evaluations, i.e. 4,800
    pairwise rounds or 2,400 four-way rounds) with 100 campaigns per
    data source and query format (800 campaigns per format), for
    k = 2, 3, 4 and 6;
  * replaces the yes/no cycle test by an evidence level. For every
    observed class pair the Jeffreys posterior Beta(W_ij + 1/2,
    W_ji + 1/2) gives q_ij = Pr(class i truly beats class j more often
    than not). A fully observed triangle is cyclic with probability
        c = q_ab q_bc q_ca + (1 - q_ab)(1 - q_bc)(1 - q_ca)
    (edges treated as independent). A campaign "shows a cycle at
    evidence level e" if some triangle has c >= e. e = 0.5 is close to
    the old point test; e = 0.05 counts loose loops that the data cannot
    rule out; e = 0.95 counts only tightly closed ones. The shaded
    region of the figure is this family of curves over e.
  * records per campaign: expected number of cyclic triangles (sum of
    c), the Hodge loop share (curl plus harmonic energy fraction of the
    log-odds flow on the observed comparison graph) and its global part,
    so every quantity has a distribution over campaigns (5/25/50/75/95
    percentiles) and a bootstrap interval of the mean.

Population truth for reference: the class tournaments of every source
have loop share 0.006 to 0.043 (decision_body_cycles.py), so a loop that
persists at large budgets would point to structure, and one that fades
is sampling.

Data: the version-locked band-gap archive (sha256[:16] recorded); one
evaluation = one archived record; no invented observations.

ASCII only. Usage (from shared/):
    uv run python scripts/campaign_loop_bands.py [--campaigns 100]
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import argparse
import hashlib
import importlib.util
import itertools
import json
from multiprocessing import Pool

import numpy as np
from scipy.special import betainc


def _load(name):
    path = Path(__file__).resolve().parent / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


lc = _load("local_cycles")

BUDGETS = (8, 16, 24, 32, 48, 64, 96, 128, 176, 240, 360, 480, 720, 960,
           1440, 1920, 2880, 3840, 5760, 7680, 9600)
LEVELS = (0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
KS = (2, 3, 4, 6)
PCTS = (5, 25, 50, 75, 95)


def graded(W, n):
    """Per-checkpoint graded loop measures on the win table W."""
    h = lc.hodge_observed(W, n)
    if h is None:
        return None
    N = W + W.T
    with np.errstate(invalid="ignore"):
        # Pr(p_ij > 1/2) under Beta(W_ij + 1/2, W_ji + 1/2)
        q = np.where(N > 0, 1.0 - betainc(W + 0.5, W.T + 0.5, 0.5), 0.5)
    cmax, csum = 0.0, 0.0
    for a, b, c in itertools.combinations(range(n), 3):
        if N[a, b] > 0 and N[b, c] > 0 and N[a, c] > 0:
            p = q[a, b] * q[b, c] * q[c, a] + q[b, a] * q[c, b] * q[a, c]
            csum += p
            cmax = max(cmax, p)
    return {"loop": h["curl"] + h["harm"], "harm": h["harm"],
            "exp_cycles": csum, "max_cycle_prob": cmax,
            "point_cycle": float(h["n_cyc"] > 0)}


def campaign(task):
    classes, k, seed = task
    rng = np.random.default_rng(seed)
    names = sorted(classes)
    n = len(names)
    arrs = [np.asarray(classes[c], float) for c in names]
    W = np.zeros((n, n))
    spent, bi, out = 0, 0, []
    kk = min(k, n)
    while bi < len(BUDGETS):
        pick = rng.choice(n, size=kk, replace=False)
        draws = sorted(((arrs[c][rng.integers(len(arrs[c]))] + 1e-9 * rng.random(), c)
                        for c in pick), reverse=True)
        for a in range(kk):
            for b in range(a + 1, kk):
                W[draws[a][1], draws[b][1]] += 1
        spent += kk
        while bi < len(BUDGETS) and spent >= BUDGETS[bi]:
            out.append(graded(W, n))
            bi += 1
    return out


def summarise(vals, rng, n_boot=500):
    v = np.asarray(vals, float)
    boots = [v[rng.integers(0, len(v), len(v))].mean() for _ in range(n_boot)]
    return {"mean": float(v.mean()),
            "mean_ci90": [float(np.percentile(boots, 5)), float(np.percentile(boots, 95))],
            **{f"p{p:02d}": float(np.percentile(v, p)) for p in PCTS}}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--universe", default="output/bandgap_archive_build/bandgap_universe.jsonl")
    ap.add_argument("--output-dir", default="output/campaign_loop_bands")
    ap.add_argument("--campaigns", type=int, default=100, help="per source and k")
    ap.add_argument("--workers", type=int, default=16)
    args = ap.parse_args()
    root = Path(__file__).resolve().parent.parent
    universe = Path(args.universe)
    if not universe.is_absolute():
        universe = (root / universe).resolve()
    sha = hashlib.sha256(universe.read_bytes()).hexdigest()[:16]
    out = (root / args.output_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    recs = lc.load(universe)
    sources = list(lc.LABEL)
    whole = {}
    for s in sources:
        d = {}
        for src, f, c, g, cs, es in recs:
            if src == s:
                d.setdefault(c, []).append(g)
        whole[s] = {c: sorted(v) for c, v in d.items() if len(v) >= lc.MIN_CLASS}
    rng = np.random.default_rng(0)
    res = {"archive_sha256_16": sha,
           "config": {"budgets": BUDGETS, "levels": LEVELS, "ks": KS,
                      "campaigns_per_source": args.campaigns,
                      "sources": sources, "min_class": lc.MIN_CLASS},
           "by_k": {}}
    with Pool(args.workers) as pool:
        for k in KS:
            tasks = [(whole[s], k, 7919 * k + 100003 * si + c)
                     for si, s in enumerate(sources) for c in range(args.campaigns)]
            runs = pool.map(campaign, tasks, chunksize=4)
            rows = {}
            for bi, b in enumerate(BUDGETS):
                g = [r[bi] for r in runs if r[bi] is not None]
                mx = np.array([x["max_cycle_prob"] for x in g])
                rows[str(b)] = {
                    "n_campaigns": len(g),
                    "share_point_cycle": summarise([x["point_cycle"] for x in g], rng),
                    "share_cycle_at_level": {str(e): float((mx >= e).mean()) for e in LEVELS},
                    "expected_cycles": summarise([x["exp_cycles"] for x in g], rng),
                    "loop_share": summarise([x["loop"] for x in g], rng),
                    "global_loop_share": summarise([x["harm"] for x in g], rng),
                }
            res["by_k"][str(k)] = rows
            print("k", k, {b: round(rows[str(b)]["share_cycle_at_level"]["0.5"], 3)
                           for b in (240, 2400 if 2400 in BUDGETS else 1920, 9600)},
                  flush=True)
    (out / "campaign_loop_bands.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    L = ["# Loops seen during a campaign, with uncertainty", "",
         f"Archive sha256[:16] `{sha}`; {args.campaigns} campaigns per source and "
         f"query size ({args.campaigns * len(sources)} per size); budgets to {BUDGETS[-1]} "
         "evaluations.", "",
         "Share of campaigns with a cycle at evidence level e (loose e = 0.05, "
         "middle e = 0.5, strict e = 0.95), and the median loop share "
         "(5th to 95th percentile over campaigns).", "",
         "| k | budget | e = 0.05 | e = 0.5 | e = 0.95 | expected cycles, mean | loop share, median (p05 to p95) |",
         "|---|---|---|---|---|---|---|"]
    for k in KS:
        for b in (24, 96, 240, 960, 2880, 9600):
            r = res["by_k"][str(k)][str(b)]
            lv, ls = r["share_cycle_at_level"], r["loop_share"]
            L.append(f"| {k} | {b} | {lv['0.05']:.2f} | {lv['0.5']:.2f} | {lv['0.95']:.2f} | "
                     f"{r['expected_cycles']['mean']:.2f} | {ls['p50']:.3f} "
                     f"({ls['p05']:.3f} to {ls['p95']:.3f}) |")
    (out / "campaign_loop_bands.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
