"""
Hypergraph rating trial on the frozen band-gap archive band-gap universe.

Materials counterpart of scripts/hypergraph_trial.py. Entities are
(channel, material-class) pairs; two kinds of hyperedges are built from the
frozen o=(x,s,m,y,u,p) inventory, with NO synthetic jitter -- every event
compares real archived records:

  1. Within-channel class events: for event t, each class contributes the
     value of one real record (deterministic stride indexing into its sorted
     value list); classes are ranked by band gap (wide-gap-first guidance).
     These identify the class strengths theta.
  2. Cross-channel bridge events: a formula measured/computed in >= 2
     channels yields one hyperedge ranking its (channel, class) entities by
     the channel-specific gap of the SAME formula. Class cancels, so these
     identify the channel offsets delta -- which should recover the known
     fidelity ordering (hybrid/meta-GGA gaps wider than GGA).

The pooled Plackett-Luce factorisation w = exp(theta_class + delta_channel)
then gives ONE combined class ranking across all 8 band-gap channels plus a
per-channel systematic offset: the "combine all ratings together" view that
pairwise per-channel graphs cannot express.

Only records with gap > 0.01 eV enter (metals are all exactly 0 and carry no
ordering signal; exact ties are dropped by the hyperedge builder anyway).
Class taxonomy mirrors classify_by_anion in aflow_mp_jarvis_phase2.py
(halide > oxide > nitride > chalcogenide > carbide > boride > pnictide >
other).

Usage (from shared/):
    uv run python scripts/hypergraph_trial_bandgap.py
Deterministic: no RNG anywhere. ASCII only.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import argparse
import json
import os
import re
import tempfile

import numpy as np

from src.rating.hypergraph_rating import (
    Hyperedge,
    PlackettLuceModel,
    PooledPlackettLuce,
    clique_expansion_hodge,
    hyperedge_from_scores,
    hyperedge_residuals,
    _kendall_tau,
)

# Mirrors classify_by_anion in scripts/aflow_mp_jarvis_phase2.py
ANION_PRIORITY = [
    ("halide", {"F", "Cl", "Br", "I"}),
    ("oxide", {"O"}),
    ("nitride", {"N"}),
    ("chalcogenide", {"S", "Se", "Te"}),
    ("carbide", {"C"}),
    ("boride", {"B"}),
    ("pnictide", {"P", "As", "Sb", "Bi"}),
]

STRIDE = 104729  # deterministic pseudo-shuffle stride (prime)
MIN_GAP_EV = 0.01


def classify_formula(formula: str) -> str:
    species = set(re.findall(r"[A-Z][a-z]?", formula))
    for name, elems in ANION_PRIORITY:
        if species & elems:
            return name
    return "other"


def load_universe(path: Path):
    """-> per_channel[channel][class] = sorted values; formulas[formula][channel] = value."""
    per_channel = {}
    formulas = {}
    n_used = n_seen = 0
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            n_seen += 1
            formula = (rec.get("x") or {}).get("formula")
            gap = (rec.get("y") or {}).get("bandgap_eV")
            channel = rec.get("s")
            if not formula or gap is None or channel is None:
                continue
            gap = float(gap)
            if gap <= MIN_GAP_EV:
                continue
            n_used += 1
            cls = classify_formula(formula)
            per_channel.setdefault(channel, {}).setdefault(cls, []).append(gap)
            # keep the first record per (formula, channel): deterministic
            formulas.setdefault(formula, {}).setdefault(channel, gap)
    for channel in per_channel:
        for cls in per_channel[channel]:
            per_channel[channel][cls].sort()
    return per_channel, formulas, n_seen, n_used


def within_channel_edges(per_channel, n_events, min_class_size):
    edges = {}
    for channel in sorted(per_channel):
        classes = {c: v for c, v in per_channel[channel].items()
                   if len(v) >= min_class_size}
        if len(classes) < 3:
            continue
        ch_edges = []
        for t in range(n_events):
            scores = {}
            for cls in sorted(classes):
                vals = classes[cls]
                scores[cls] = vals[(t * STRIDE + 12345) % len(vals)]
            e = hyperedge_from_scores(scores, higher_is_better=True,
                                      channel=channel, metadata={"event": t})
            if e is not None:
                ch_edges.append(e)
        if len(ch_edges) >= 10:
            edges[channel] = ch_edges
    return edges


def bridge_edges(formulas, channels_kept, max_bridges):
    out = []
    for formula in sorted(formulas):
        chans = {c: v for c, v in formulas[formula].items() if c in channels_kept}
        if len(chans) < 2:
            continue
        cls = classify_formula(formula)
        # One entity per channel; the class factor is shared by construction,
        # so these events carry pure channel-offset signal.
        scores = {PooledPlackettLuce.pooled_key(c, cls): v
                  for c, v in chans.items()}
        e = hyperedge_from_scores(scores, higher_is_better=True,
                                  metadata={"formula": formula})
        if e is not None:
            out.append(e)
        if len(out) >= max_bridges:
            break
    return out


def atomic_write(path: Path, text: str):
    fd, tmp = tempfile.mkstemp(dir=str(path.parent), suffix=".tmp")
    with os.fdopen(fd, "w", encoding="utf-8") as fh:
        fh.write(text)
    os.replace(tmp, str(path))


def main():
    ap = argparse.ArgumentParser(description="Hypergraph trial on band-gap archive universe")
    ap.add_argument("--universe",
                    default="output/bandgap_archive_build/bandgap_universe.jsonl")
    ap.add_argument("--n-events", type=int, default=300)
    ap.add_argument("--min-class-size", type=int, default=30)
    ap.add_argument("--max-bridges", type=int, default=4000)
    ap.add_argument("--output-dir", default="output/hypergraph_trial_bandgap")
    args = ap.parse_args()

    root = Path(__file__).parent.parent
    universe = (root / args.universe).resolve()
    if not universe.exists():
        print(f"[error] universe not found: {universe}; "
              "run compute/bandgap_archive/build_inventory.py first")
        return 1
    out_dir = (root / args.output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    per_channel, formulas, n_seen, n_used = load_universe(universe)
    print(f"[load] {n_seen} records, {n_used} with gap > {MIN_GAP_EV} eV, "
          f"{len(per_channel)} channels")

    ch_edges = within_channel_edges(per_channel, args.n_events,
                                    args.min_class_size)
    results = {"config": vars(args), "channels": {}, "pooled": None}

    pooled_edges = []
    for channel, edges in sorted(ch_edges.items()):
        model = PlackettLuceModel().fit(edges)
        diag = hyperedge_residuals(model, edges)
        hodge = clique_expansion_hodge(edges)
        mean_order = sorted(
            {c for e in edges for c in e.ranking},
            key=lambda c: (-float(np.mean(per_channel[channel][c])), c))
        results["channels"][channel] = {
            "n_hyperedges": len(edges),
            "classes": sorted({c for e in edges for c in e.ranking}),
            "pl_ranking": model.ranking(),
            "pl_log_strengths": model.log_strengths,
            "mean_edge_residual": diag["mean_edge_residual"],
            "pairwise_cyclic_energy_fraction": hodge["cyclic_energy_fraction"],
            "tau_pl_vs_mean_value_order": _kendall_tau(model.ranking(),
                                                       mean_order),
        }
        print(f"[ok] {channel}: {len(edges)} events, ranking = "
              f"{model.ranking()}, tau vs mean-order = "
              f"{results['channels'][channel]['tau_pl_vs_mean_value_order']:.3f}")
        for e in edges:
            pooled_edges.append(Hyperedge(
                ranking=tuple(PooledPlackettLuce.pooled_key(channel, c)
                              for c in e.ranking)))

    bridges = bridge_edges(formulas, set(ch_edges), args.max_bridges)
    print(f"[ok] {len(bridges)} cross-channel bridge events "
          f"(same formula, >= 2 channels)")
    pooled_edges.extend(bridges)

    if len(ch_edges) >= 2 and bridges:
        pooled = PooledPlackettLuce().fit(pooled_edges)
        offsets = dict(sorted(pooled.channel_offsets.items(),
                              key=lambda kv: -kv[1]))
        results["pooled"] = {
            "n_within_events": len(pooled_edges) - len(bridges),
            "n_bridge_events": len(bridges),
            "converged": pooled.converged,
            "combined_class_ranking": pooled.ranking(),
            "class_log_strengths": pooled.policy_log_strengths,
            "channel_offsets_sorted": offsets,
        }
        print("[ok] pooled combined class ranking: "
              f"{pooled.ranking()}")
        print("[ok] channel offsets (log scale, wide-gap-propensity first):")
        for c, v in offsets.items():
            print(f"       {c:>18s}  {v:+.3f}")

    atomic_write(out_dir / "hypergraph_trial_bandgap.json",
                 json.dumps(results, indent=2, sort_keys=True))

    lines = ["# Hypergraph trial on the frozen band-gap archive universe", "",
             f"Records used (gap > {MIN_GAP_EV} eV): {n_used} of {n_seen}.",
             "One within-channel event ranks material classes by real archived",
             "gaps (deterministic stride draw); one bridge event ranks the",
             "same formula across channels (identifies fidelity offsets).", ""]
    lines.append("| channel | events | PL class ranking (best first) | "
                 "edge residual | pairwise cyclic energy | tau vs mean order |")
    lines.append("|---|---|---|---|---|---|")
    for channel, e in sorted(results["channels"].items()):
        lines.append(f"| {channel} | {e['n_hyperedges']} | "
                     f"{' > '.join(e['pl_ranking'])} | "
                     f"{e['mean_edge_residual']:.3f} | "
                     f"{e['pairwise_cyclic_energy_fraction']:.3f} | "
                     f"{e['tau_pl_vs_mean_value_order']:.3f} |")
    lines.append("")
    if results["pooled"]:
        p = results["pooled"]
        lines += ["## Pooled cross-channel fit", "",
                  f"{p['n_within_events']} within-channel events + "
                  f"{p['n_bridge_events']} bridge events.", "",
                  "Combined class ranking: " +
                  " > ".join(p["combined_class_ranking"]), "",
                  "Channel offsets (log wide-gap propensity, descending):", ""]
        for c, v in p["channel_offsets_sorted"].items():
            lines.append(f"- `{c}`: {v:+.3f}")
        lines.append("")
    atomic_write(out_dir / "hypergraph_trial_bandgap.md", "\n".join(lines) + "\n")
    print(f"[done] wrote {out_dir / 'hypergraph_trial_bandgap.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
