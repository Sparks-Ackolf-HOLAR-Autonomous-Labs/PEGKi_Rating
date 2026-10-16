"""
Does the pooled cross-source model have to be additive? (learned departure)

The pooled model behind the fidelity-ladder result is additive by design:
log w[source, class] = theta_class + delta_source. That was a modelling
choice, not a finding. This experiment keeps the additive model as the cold
start and lets the version-locked band-gap archive decide how far to move
from it, using src/rating/interaction_pooled.InteractionPooledPlackettLuce:

    log w[c, i] = a_c * theta_i + delta_c + U_c . V_i

Models compared (same events, same split):
    additive     a_c = 1, no interaction (the committed pooled model)
    scale        learned per-source margin scale a_c
    scale+r1     scale plus a rank-1 source x class interaction
    scale+r2     scale plus a rank-2 interaction

The penalty that pulls each model back to the additive cold start is chosen
on a validation split (grid below), the chosen model is refitted on
train+validation, and all models are scored on held-out test events.
Reported: held-out log-likelihood per event against the additive model, with
event-bootstrap CIs, split by event type (within-source class rankings versus
cross-source bridge rankings of one compound), the learned margin scales, and
the class-dependent ladder spacings implied by the chosen model next to the
separately fitted family-resolved spacings (scripts/family_offsets.py).

Events are built exactly as in scripts/hypergraph_trial_bandgap.py (imported).

Usage (from shared/):
    uv run python scripts/nonadditive_pooled.py
    uv run python scripts/nonadditive_pooled.py --smoke
Outputs: nonadditive_pooled.json / .md / .png in --output-dir. ASCII only.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import argparse
import hashlib
import importlib.util
import json
import time

import numpy as np

from src.rating.hypergraph_rating import Hyperedge, PooledPlackettLuce  # noqa: E402
from src.rating.interaction_pooled import InteractionPooledPlackettLuce  # noqa: E402


def _load_trial():
    path = Path(__file__).resolve().parent / "hypergraph_trial_bandgap.py"
    spec = importlib.util.spec_from_file_location("hypergraph_trial_bandgap", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


trial = _load_trial()
SEP = PooledPlackettLuce.SEP
MODELS = {
    "additive": dict(rank=0, fit_scale=False),
    "scale": dict(rank=0, fit_scale=True),
    "scale+r1": dict(rank=1, fit_scale=True),
    "scale+r2": dict(rank=2, fit_scale=True),
}
LAM_GRID = (0.03, 0.3, 3.0, 30.0, 300.0, 3000.0)
SPACINGS = [
    ("snumat_hse", "snumat_gga"),
    ("jarvis_mbj", "jarvis_optb88vdw"),
    ("expt_zhuo", "mp_gga"),
]
FAMILIES = ["oxide", "halide", "chalcogenide", "nitride"]


def build_events(universe, n_events, min_class_size, max_bridges, seed):
    per_channel, formulas, n_seen, n_used = trial.load_universe(universe)
    ch_edges = trial.within_channel_edges(per_channel, n_events, min_class_size)
    kept = sorted(ch_edges)
    within = []
    for ch in kept:
        for e in ch_edges[ch]:
            within.append(Hyperedge(ranking=tuple(
                PooledPlackettLuce.pooled_key(ch, c) for c in e.ranking),
                metadata={"type": "within"}))
    bridges = trial.bridge_edges(formulas, set(kept), 10 ** 9)
    rng = np.random.default_rng(seed)
    order = rng.permutation(len(bridges))[:max_bridges]
    bridges = [Hyperedge(ranking=bridges[j].ranking,
                         metadata={"type": "bridge",
                                   "formula": bridges[j].metadata.get("formula")})
               for j in sorted(order)]
    # keep only bridge events whose classes also appear in within events,
    # so every held-out entity has a fitted theta
    classes = {k.split(SEP, 1)[1] for e in within for k in e.ranking}
    bridges = [e for e in bridges if e.ranking[0].split(SEP, 1)[1] in classes]
    return within, bridges, kept, n_seen, n_used


def split(events, rng, fr=(0.6, 0.2)):
    idx = rng.permutation(len(events))
    a = int(fr[0] * len(events)); b = int((fr[0] + fr[1]) * len(events))
    return ([events[j] for j in idx[:a]], [events[j] for j in idx[a:b]],
            [events[j] for j in idx[b:]])


def fit_model(name, lam, edges, seed):
    spec = MODELS[name]
    return InteractionPooledPlackettLuce(
        rank=spec["rank"], fit_scale=spec["fit_scale"],
        lam_scale=lam, lam_int=lam, seed=seed).fit(edges)


def boot_ci(x, rng, B=2000):
    x = np.asarray(x, float)
    means = np.array([x[rng.integers(0, len(x), len(x))].mean() for _ in range(B)])
    return [float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))]


def implied_spacings(model):
    M = model.source_class_matrix()
    out = {}
    for hi, lo in SPACINGS:
        if hi in M and lo in M:
            out[f"{hi}-{lo}"] = {f: float(M[hi][f] - M[lo][f])
                                 for f in FAMILIES if f in M[hi]}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--universe",
                    default="output/bandgap_archive_build/bandgap_universe.jsonl")
    ap.add_argument("--n-events", type=int, default=300)
    ap.add_argument("--min-class-size", type=int, default=30)
    ap.add_argument("--max-bridges", type=int, default=12000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--n-scale-boot", type=int, default=100)
    ap.add_argument("--output-dir", default="output/nonadditive_pooled")
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    if args.smoke:
        args.n_events, args.max_bridges, args.n_scale_boot = 60, 1500, 3
        args.output_dir = args.output_dir + "_smoke"

    root = Path(__file__).resolve().parent.parent
    universe = (root / args.universe).resolve()
    sha = hashlib.sha256(universe.read_bytes()).hexdigest()[:16]
    out = (root / args.output_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    within, bridges, kept, n_seen, n_used = build_events(
        universe, args.n_events, args.min_class_size, args.max_bridges, args.seed)
    print(f"[load] sha {sha}; {n_used}/{n_seen} records; sources {kept}; "
          f"{len(within)} within events, {len(bridges)} bridge events")
    rng = np.random.default_rng(args.seed)
    w_tr, w_va, w_te = split(within, rng)
    b_tr, b_va, b_te = split(bridges, rng)
    train, val, test = w_tr + b_tr, w_va + b_va, w_te + b_te

    selection, final = {}, {}
    for name in MODELS:
        grid = [None] if name == "additive" else LAM_GRID
        best = None
        for lam in grid:
            m = fit_model(name, lam if lam is not None else 1.0, train, args.seed)
            v = float(np.mean(m.event_logliks(val)))
            selection.setdefault(name, []).append({"lam": lam, "val_ll": v})
            print(f"[val] {name:9s} lam={lam}: {v:.5f} (conv {m.converged})")
            if best is None or v > best[1]:
                best = (lam, v)
        lam = best[0]
        m = fit_model(name, lam if lam is not None else 1.0, train + val, args.seed)
        final[name] = {"lam": lam, "model": m}

    base_w = final["additive"]["model"].event_logliks(w_te)
    base_b = final["additive"]["model"].event_logliks(b_te)
    brng = np.random.default_rng(args.seed + 7)
    results_models = {}
    for name, d in final.items():
        m = d["model"]
        lw, lb = m.event_logliks(w_te), m.event_logliks(b_te)
        diff_all = np.concatenate([lw - base_w, lb - base_b])
        results_models[name] = {
            "lam": d["lam"], "n_parameters": m.n_parameters(),
            "converged": m.converged,
            "test_ll_per_event": float(np.mean(np.concatenate([lw, lb]))),
            "gain_vs_additive": float(diff_all.mean()),
            "gain_ci": boot_ci(diff_all, brng) if name != "additive" else [0.0, 0.0],
            "gain_within": float(np.mean(lw - base_w)),
            "gain_within_ci": boot_ci(lw - base_w, brng) if name != "additive" else [0.0, 0.0],
            "gain_bridge": float(np.mean(lb - base_b)),
            "gain_bridge_ci": boot_ci(lb - base_b, brng) if name != "additive" else [0.0, 0.0],
            "margin_scales": m.margin_scales,
            "source_class_log_strength": m.source_class_matrix(),
            "implied_spacings": implied_spacings(m),
        }
        print(f"[test] {name:9s} gain {diff_all.mean():+.4f} "
              f"(within {np.mean(lw - base_w):+.4f}, bridge {np.mean(lb - base_b):+.4f})")

    # bootstrap the learned margin scales (event resampling, chosen lambda)
    lam_s = final["scale"]["lam"]
    allev = train + val
    reps = []
    for r in range(args.n_scale_boot):
        rr = np.random.default_rng([args.seed, 55, r])
        samp = [allev[j] for j in rr.integers(0, len(allev), len(allev))]
        reps.append(fit_model("scale", lam_s, samp, args.seed).margin_scales)
    scale_ci = {c: [float(np.percentile([x[c] for x in reps], 2.5)),
                    float(np.percentile([x[c] for x in reps], 97.5))]
                for c in reps[0]} if reps else {}

    fam_path = root / "output/family_offsets/family_offsets.json"
    fam_ref = {}
    if fam_path.exists():
        fj = json.load(open(fam_path))
        for f in FAMILIES:
            if f in fj.get("families", {}):
                fam_ref[f] = fj["families"][f].get("spacings")

    best_name = max((n for n in MODELS if n != "additive"),
                    key=lambda n: results_models[n]["gain_vs_additive"])
    res = {
        "config": vars(args), "universe_sha256_16": sha,
        "sources": kept, "n_within": len(within), "n_bridge": len(bridges),
        "n_test_within": len(w_te), "n_test_bridge": len(b_te),
        "lam_grid": LAM_GRID, "selection": selection,
        "models": results_models, "best_nonadditive": best_name,
        "scale_bootstrap_ci": scale_ci,
        "family_offsets_reference_spacings": fam_ref,
        "runtime_s": round(time.time() - t0, 1),
    }
    (out / "nonadditive_pooled.json").write_text(json.dumps(res, indent=2))

    L = ["# Non-additive pooled model: learned departure from the additive cold start", "",
         f"Archive sha256[:16] = `{sha}`; sources: {', '.join(kept)}.",
         f"Events: {len(within)} within-source class rankings, {len(bridges)} "
         f"cross-source bridge rankings; split 60/20/20 (train/validation/test).",
         f"Penalty grid (pull toward additive): {LAM_GRID}; chosen on validation, "
         "refitted on train+validation, scored on test.", "",
         "| model | chosen lambda | parameters | held-out gain per event vs additive [95% CI] | within-source | bridge |",
         "|---|---|---|---|---|---|"]
    for n, r in results_models.items():
        L.append(f"| {n} | {r['lam']} | {r['n_parameters']} | {r['gain_vs_additive']:+.4f} "
                 f"[{r['gain_ci'][0]:+.4f}, {r['gain_ci'][1]:+.4f}] | {r['gain_within']:+.4f} "
                 f"[{r['gain_within_ci'][0]:+.4f}, {r['gain_within_ci'][1]:+.4f}] | "
                 f"{r['gain_bridge']:+.4f} [{r['gain_bridge_ci'][0]:+.4f}, {r['gain_bridge_ci'][1]:+.4f}] |")
    L += ["", "## Learned margin scales a_c (scale model; 1 = additive)", "",
          "| source | a_c | bootstrap 95% CI |", "|---|---|---|"]
    ms = results_models["scale"]["margin_scales"]
    for c in sorted(ms, key=ms.get, reverse=True):
        ci = scale_ci.get(c, [float("nan")] * 2)
        L.append(f"| {c} | {ms[c]:.3f} | [{ci[0]:.3f}, {ci[1]:.3f}] |")
    L += ["", f"## Class-dependent ladder spacings implied by `{best_name}`", "",
          "| spacing | " + " | ".join(FAMILIES) + " |", "|---|" + "---|" * len(FAMILIES)]
    for sp, v in results_models[best_name]["implied_spacings"].items():
        L.append(f"| {sp} | " + " | ".join(
            f"{v[f]:+.3f}" if f in v else "NA" for f in FAMILIES) + " |")
    if fam_ref:
        L += ["", "Family-resolved spacings fitted separately (scripts/family_offsets.py):", "",
              "| spacing | " + " | ".join(FAMILIES) + " |", "|---|" + "---|" * len(FAMILIES)]
        keymap = {"snumat_hse-snumat_gga": "hse_minus_gga_snumat",
                  "jarvis_mbj-jarvis_optb88vdw": "mbj_minus_optb88_jarvis",
                  "expt_zhuo-mp_gga": "expt_minus_gga_mp"}
        for sp, k in keymap.items():
            L.append(f"| {sp} | " + " | ".join(
                f"{fam_ref[f][k]:+.3f}" if f in fam_ref and fam_ref[f].get(k) is not None
                else "NA" for f in FAMILIES) + " |")
    L += ["", f"Runtime {res['runtime_s']} s."]
    (out / "nonadditive_pooled.md").write_text("\n".join(L) + "\n")
    print(f"[done] {out} ({res['runtime_s']} s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
