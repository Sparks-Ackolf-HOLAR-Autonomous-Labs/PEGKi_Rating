"""
Regret-based omega from the measured paired-campaign logs (band-gap archive).

The recorded PEGKi limit says the omega term measures CLASS
INFORMATIVENESS (dispersion of class-as-policy MAE-equivalents), not
true campaign regret (writings/OMEGA_DECISION.md, "semantics caveat").
This script computes a REGRET-BASED omega from the committed exp 26
paired-campaign logs and compares it with the interim dispersion omega
on matched directed channel pairs. It is an ANALYSIS of committed logs
-- it runs no new campaigns and queries no new observations.

Definition
----------
omega_regret(A -> B) = mean over matched seeds of the average over
EARLY-WINDOW checkpoints of max(0, tau_cold(b) - tau_full(b)), where
"full" is the fixed-trust full-reuse condition of exp 26 (lambda =
config lam) -- the expected early harm of following donor A's advice
on receiver B. The pairwise arm is primary; the k-way arm is reported.
omega_regret_final uses the final checkpoint only. Percentile-bootstrap
CIs over seeds.

Comparison with the interim dispersion omega
--------------------------------------------
Phase-3 engine names map to universe channels ONLY for:
  mbj_band_gap_theoretical        -> jarvis_mbj
  optb88vdw_band_gap_theoretical  -> jarvis_optb88vdw
  snumat_gga_band_gap_theoretical -> snumat_gga
  snumat_hse_band_gap_theoretical -> snumat_hse
(the gga_band_gap_difference/intersection composites are ignored).
On the 12 matched directed pairs: Spearman and Kendall correlation
between omega_regret and omega_total, a scatter, a per-pair table.
The full 56-cell omega_regret matrix is reported regardless -- the
regret omega exists for every cell the campaigns measured (including
expt_zhuo), which the interim omega does not.

Outputs (default ../output/omega_regret/ = shared/output/omega_regret/)
-----------------------------------------------------------------------
  omega_regret.json / .md / .png

Run (from shared/):
  uv run python scripts/omega_regret.py
Self-test (synthetic curves with known regret incl. zero-regret and
all-harm cases; name-mapping test): --self-test. Exits nonzero on
failure. ASCII only; seeded and reproducible (--seed).

Note: the small checkpoint loader is deliberately duplicated with
campaign_protocols/code/experiments/29_cost_calibration.py so both scripts stay
independently runnable.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.abspath(__file__))

ENGINE_TO_CHANNEL = {
    "mbj_band_gap_theoretical": "jarvis_mbj",
    "optb88vdw_band_gap_theoretical": "jarvis_optb88vdw",
    "snumat_gga_band_gap_theoretical": "snumat_gga",
    "snumat_hse_band_gap_theoretical": "snumat_hse",
}
ARMS = ("pair", "kway")


def _default_input_dir() -> Path:
    return (Path(_HERE) / ".." / "output" / "bandgap_cluster_run2" /
            "paired_campaign_k4").resolve()


def _default_phase3_omega() -> Path:
    return (Path(_HERE) / ".." / "output" / "bandgap_cluster_run1" /
            "jarvis_phase3" / "omega_decomposition.json").resolve()


def _default_output_dir() -> Path:
    return (Path(_HERE) / ".." / "output" / "omega_regret").resolve()


# ---------------------------------------------------------------------------
# Loading (duplicated tiny loader; see module docstring)
# ---------------------------------------------------------------------------
def load_aggregate(input_dir: Path) -> dict:
    with open(input_dir / "paired_campaign.json", "r",
              encoding="utf-8") as f:
        return json.load(f)


def load_checkpoint_rows(input_dir: Path) -> list:
    rows = []
    for path in sorted((input_dir / "checkpoints").glob("*.json")):
        with open(path, "r", encoding="utf-8") as f:
            rows.append(json.load(f))
    return rows


# ---------------------------------------------------------------------------
# Regret computation
# ---------------------------------------------------------------------------
def bootstrap_ci(values: np.ndarray, rng: np.random.Generator,
                 n_boot: int = 2000):
    if len(values) < 2:
        return [float(values.mean()), float(values.mean())]
    means = [float(values[rng.integers(0, len(values), len(values))].mean())
             for _ in range(n_boot)]
    return [float(np.percentile(means, 2.5)),
            float(np.percentile(means, 97.5))]


def seed_regret(cold, full, n_early):
    """Early regret for one seed: mean over the first n_early checkpoints
    of max(0, tau_cold - tau_full)."""
    c = np.asarray(cold[:n_early], dtype=float)
    f = np.asarray(full[:n_early], dtype=float)
    return float(np.mean(np.maximum(0.0, c - f)))


def compute_omega_regret(rows, n_early, seed):
    """-> {arm: {donor: {receiver: {omega_regret, ci, omega_regret_final,
    final_ci, n_seeds}}}} from the fixed-trust full-reuse curves."""
    per_cell: dict = {}
    for row in rows:
        rec = row["receiver"]
        curves = row["curves"]
        for key, taus in curves.items():
            parts = key.split("|")
            if parts[0] != "full":
                continue
            _, arm, donor = parts
            cold = curves.get(f"cold|{arm}")
            if cold is None:
                continue
            early = seed_regret(cold, taus, n_early)
            final = float(max(0.0, cold[-1] - taus[-1]))
            per_cell.setdefault((arm, donor, rec),
                                []).append((early, final))
    rng = np.random.default_rng(seed)
    out: dict = {a: {} for a in ARMS}
    for (arm, donor, rec) in sorted(per_cell):
        vals = per_cell[(arm, donor, rec)]
        early = np.asarray([v[0] for v in vals], dtype=float)
        final = np.asarray([v[1] for v in vals], dtype=float)
        out.setdefault(arm, {}).setdefault(donor, {})[rec] = {
            "omega_regret": float(early.mean()),
            "omega_regret_ci": bootstrap_ci(early, rng),
            "omega_regret_final": float(final.mean()),
            "omega_regret_final_ci": bootstrap_ci(final, rng),
            "n_seeds": len(vals),
        }
    return out


# ---------------------------------------------------------------------------
# Interim-omega matching and rank correlations (no scipy dependency)
# ---------------------------------------------------------------------------
def match_pairs(phase3_omega: dict, regret_pair_arm: dict):
    """-> list of {donor, receiver, engine_pair, omega_total,
    omega_regret} on the mapped directed pairs only."""
    matched = []
    for key, vals in sorted(phase3_omega.items()):
        if "->" not in key:
            continue
        a, b = key.split("->")
        if a not in ENGINE_TO_CHANNEL or b not in ENGINE_TO_CHANNEL:
            continue
        donor, rec = ENGINE_TO_CHANNEL[a], ENGINE_TO_CHANNEL[b]
        cell = regret_pair_arm.get(donor, {}).get(rec)
        if cell is None:
            continue
        matched.append({
            "engine_pair": key,
            "donor": donor,
            "receiver": rec,
            "omega_total_interim": float(vals["omega_total"]),
            "omega_reversal_interim": float(vals.get("omega_reversal",
                                                     float("nan"))),
            "omega_regret": cell["omega_regret"],
            "omega_regret_ci": cell["omega_regret_ci"],
        })
    return matched


def _rankdata(x):
    x = np.asarray(x, dtype=float)
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(len(x), dtype=float)
    i = 0
    while i < len(x):
        j = i
        while j + 1 < len(x) and x[order[j + 1]] == x[order[i]]:
            j += 1
        ranks[order[i:j + 1]] = 0.5 * (i + j) + 1.0
        i = j + 1
    return ranks


def spearman(x, y):
    rx, ry = _rankdata(x), _rankdata(y)
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    denom = float(np.sqrt((rx ** 2).sum() * (ry ** 2).sum()))
    if denom == 0.0:
        return float("nan")
    return float((rx * ry).sum() / denom)


def kendall(x, y):
    """Kendall tau-b (handles ties)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = len(x)
    conc = disc = tx = ty = 0
    for i in range(n):
        for j in range(i + 1, n):
            dx = x[i] - x[j]
            dy = y[i] - y[j]
            if dx == 0 and dy == 0:
                tx += 1
                ty += 1
            elif dx == 0:
                tx += 1
            elif dy == 0:
                ty += 1
            elif dx * dy > 0:
                conc += 1
            else:
                disc += 1
    n0 = n * (n - 1) / 2
    denom = float(np.sqrt((n0 - tx) * (n0 - ty)))
    if denom == 0.0:
        return float("nan")
    return float((conc - disc) / denom)


# ---------------------------------------------------------------------------
# Outputs
# ---------------------------------------------------------------------------
def make_plot(matched, regret_pair, receivers, out_path):
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.4))

    ax = axes[0]
    if matched:
        xs = [m["omega_total_interim"] for m in matched]
        ys = [m["omega_regret"] for m in matched]
        ax.scatter(xs, ys, color="tab:blue", zorder=3)
        for m in matched:
            ax.annotate(
                f"{m['donor'].replace('jarvis_', 'j.').replace('snumat_', 's.')}"
                f"->{m['receiver'].replace('jarvis_', 'j.').replace('snumat_', 's.')}",
                (m["omega_total_interim"], m["omega_regret"]),
                fontsize=6.5, xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel("interim dispersion omega (omega_total, phase 3)")
    ax.set_ylabel("omega_regret (early, pair arm)")
    ax.set_title("Interim dispersion omega vs regret omega\n"
                 "(12 matched directed pairs)")
    ax.grid(alpha=0.3)

    ax = axes[1]
    channels = sorted({d for d in regret_pair} |
                      {r for d in regret_pair.values() for r in d})
    M = np.full((len(channels), len(receivers)), np.nan)
    for i, d in enumerate(channels):
        for j, r in enumerate(receivers):
            cell = regret_pair.get(d, {}).get(r)
            if cell:
                M[i, j] = cell["omega_regret"]
    vmax = float(np.nanmax(M)) if np.isfinite(M).any() else 1.0
    im = ax.imshow(M, cmap="Reds", vmin=0.0, vmax=max(vmax, 1e-9))
    ax.set_xticks(range(len(receivers)),
                  [r.replace("_", "\n") for r in receivers], fontsize=7)
    ax.set_yticks(range(len(channels)),
                  [c.replace("_", " ") for c in channels], fontsize=7)
    ax.set_xlabel("receiver B")
    ax.set_ylabel("donor A")
    ax.set_title("omega_regret(A->B), early window, full reuse, pair arm\n"
                 "(exists for ALL campaign cells, incl. expt_zhuo)")
    fig.colorbar(im, ax=ax, shrink=0.85, label="expected early harm "
                 "(tau units)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def write_markdown(path, note, matched, corr, regret, top_cells):
    lines = [
        "# Regret-based omega from the paired-campaign logs",
        "",
        note,
        "",
        "omega_regret(A->B) = mean over matched seeds of the average "
        "over early-window checkpoints of max(0, tau_cold - tau_full) "
        "under fixed-trust full reuse -- the expected early harm of "
        "following donor A's advice on receiver B. Pairwise arm "
        "primary; k-way reported in omega_regret.json.",
        "",
        "## Correlation with the interim dispersion omega "
        "(12 matched directed pairs)",
        "",
        f"Spearman rho = {corr['spearman']:+.4f}; "
        f"Kendall tau = {corr['kendall']:+.4f} "
        f"(n = {corr['n_pairs']}).",
        "",
        "| donor | receiver | interim omega_total | omega_regret "
        "[95% CI] | omega_regret_final |",
        "|---|---|---|---|---|",
    ]
    reg_pair = regret["pair"]
    for m in sorted(matched, key=lambda m: -m["omega_regret"]):
        cell = reg_pair[m["donor"]][m["receiver"]]
        lines.append(
            f"| {m['donor']} | {m['receiver']} | "
            f"{m['omega_total_interim']:.4f} | "
            f"{m['omega_regret']:.4f} [{m['omega_regret_ci'][0]:.4f}, "
            f"{m['omega_regret_ci'][1]:.4f}] | "
            f"{cell['omega_regret_final']:.4f} |")
    lines += [
        "",
        "## Largest-regret cells (all 56 campaign cells, pair arm)",
        "",
        "| rank | donor | receiver | omega_regret [95% CI] | "
        "omega_regret_final |",
        "|---|---|---|---|---|",
    ]
    for i, (donor, rec, cell) in enumerate(top_cells, start=1):
        lines.append(
            f"| {i} | {donor} | {rec} | {cell['omega_regret']:.4f} "
            f"[{cell['omega_regret_ci'][0]:.4f}, "
            f"{cell['omega_regret_ci'][1]:.4f}] | "
            f"{cell['omega_regret_final']:.4f} |")
    lines += [
        "",
        "## Full pair-arm matrix (omega_regret, early window)",
        "",
    ]
    receivers = sorted({r for d in reg_pair.values() for r in d})
    header = "| donor \\ receiver | " + " | ".join(receivers) + " |"
    lines.append(header)
    lines.append("|" + "---|" * (len(receivers) + 1))
    for donor in sorted(reg_pair):
        cells = []
        for r in receivers:
            c = reg_pair[donor].get(r)
            cells.append(f"{c['omega_regret']:.4f}" if c else "-")
        lines.append(f"| {donor} | " + " | ".join(cells) + " |")
    lines += [
        "",
        "The regret omega exists for every cell the campaigns measured "
        "(including the experimental channel expt_zhuo); the interim "
        "dispersion omega exists only where phase 3 ran. K-way arm, "
        "CIs, and the final-checkpoint variant: omega_regret.json. "
        "Computed entirely from committed campaign logs; no new "
        "observations.",
    ]
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------
def self_test() -> int:
    failures = 0
    n_early = 2

    # zero-regret: full dominates cold everywhere
    r = seed_regret([0.1, 0.2, 0.3, 0.4], [0.5, 0.5, 0.5, 0.5], n_early)
    ok = r == 0.0
    print(f"[{'PASS' if ok else 'FAIL'}] zero-regret case: {r:g} "
          "(expected 0)")
    failures += 0 if ok else 1

    # all-harm: cold beats full by 0.2 at every early checkpoint
    r = seed_regret([0.5, 0.5, 0.1, 0.1], [0.3, 0.3, 0.9, 0.9], n_early)
    ok = abs(r - 0.2) < 1e-12
    print(f"[{'PASS' if ok else 'FAIL'}] all-harm case: {r:g} "
          "(expected 0.2; late checkpoints ignored)")
    failures += 0 if ok else 1

    # mixed: harm 0.2 at checkpoint 1, gain at checkpoint 2 -> mean 0.1
    r = seed_regret([0.4, 0.2], [0.2, 0.4], n_early)
    ok = abs(r - 0.1) < 1e-12
    print(f"[{'PASS' if ok else 'FAIL'}] one-sided clipping: {r:g} "
          "(expected 0.1)")
    failures += 0 if ok else 1

    # end-to-end on synthetic rows (2 seeds), incl. final variant
    rows = [
        {"receiver": "r1", "seed": 0, "curves": {
            "cold|pair": [0.5, 0.5, 0.5, 0.9],
            "full|pair|d1": [0.3, 0.3, 0.9, 0.5]}},
        {"receiver": "r1", "seed": 1, "curves": {
            "cold|pair": [0.5, 0.5, 0.5, 0.5],
            "full|pair|d1": [0.5, 0.5, 0.5, 0.5]}},
    ]
    reg = compute_omega_regret(rows, n_early, seed=0)
    cell = reg["pair"]["d1"]["r1"]
    ok = (abs(cell["omega_regret"] - 0.1) < 1e-12 and
          abs(cell["omega_regret_final"] - 0.2) < 1e-12 and
          cell["n_seeds"] == 2)
    print(f"[{'PASS' if ok else 'FAIL'}] end-to-end: omega_regret "
          f"{cell['omega_regret']:g} (expected 0.1: seeds 0.2 and 0), "
          f"final {cell['omega_regret_final']:g} (expected 0.2)")
    failures += 0 if ok else 1

    # name mapping: composites ignored, only the 4 engines map
    fake = {}
    engines = list(ENGINE_TO_CHANNEL) + ["gga_band_gap_difference_theoretical"]
    for a in engines:
        for b in engines:
            if a != b:
                fake[f"{a}->{b}"] = {"omega_total": 0.1,
                                     "omega_reversal": 0.0}
    reg_fake = {"pair": {}}
    for a in ENGINE_TO_CHANNEL.values():
        for b in ENGINE_TO_CHANNEL.values():
            if a != b:
                reg_fake["pair"].setdefault(a, {})[b] = {
                    "omega_regret": 0.05, "omega_regret_ci": [0.0, 0.1],
                    "omega_regret_final": 0.0,
                    "omega_regret_final_ci": [0.0, 0.0], "n_seeds": 2}
    matched = match_pairs(fake, reg_fake["pair"])
    chans = {(m["donor"], m["receiver"]) for m in matched}
    ok = (len(matched) == 12 and
          all(d in ENGINE_TO_CHANNEL.values() and
              r in ENGINE_TO_CHANNEL.values() for d, r in chans) and
          not any("difference" in m["engine_pair"] or
                  "intersection" in m["engine_pair"] for m in matched))
    print(f"[{'PASS' if ok else 'FAIL'}] name mapping: {len(matched)} "
          "matched pairs (expected 12), composites excluded")
    failures += 0 if ok else 1

    # correlation sanity
    x = [1.0, 2.0, 3.0, 4.0, 5.0]
    y = [2.0, 4.0, 6.0, 8.0, 10.0]
    ok = (abs(spearman(x, y) - 1.0) < 1e-12 and
          abs(kendall(x, y) - 1.0) < 1e-12 and
          abs(spearman(x, y[::-1]) + 1.0) < 1e-12 and
          abs(kendall(x, y[::-1]) + 1.0) < 1e-12)
    print(f"[{'PASS' if ok else 'FAIL'}] rank correlations: monotone -> "
          "+1, reversed -> -1")
    failures += 0 if ok else 1

    print(f"self-test: {'PASS' if failures == 0 else 'FAIL'}")
    return failures


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Regret-based omega from the exp 26 paired-campaign "
                    "logs, compared with the interim dispersion omega")
    ap.add_argument("--input-dir", default=None,
                    help="exp 26 output dir with paired_campaign.json "
                         "and checkpoints/ (default: 50-seed k4 CHPC run)")
    ap.add_argument("--phase3-omega", default=None,
                    help="omega_decomposition.json from phase 3 "
                         "(default: bandgap_cluster_run1/jarvis_phase3)")
    ap.add_argument("--output-dir", default=None)
    ap.add_argument("--seed", type=int, default=42, help="bootstrap seed")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()

    if args.self_test:
        return self_test()

    input_dir = (Path(args.input_dir).resolve() if args.input_dir
                 else _default_input_dir())
    phase3_path = (Path(args.phase3_omega).resolve() if args.phase3_omega
                   else _default_phase3_omega())
    out_dir = (Path(args.output_dir).resolve() if args.output_dir
               else _default_output_dir())
    if not (input_dir / "paired_campaign.json").exists():
        print(f"[error] aggregate not found under {input_dir}")
        return 1
    if not phase3_path.exists():
        print(f"[error] phase-3 omega not found: {phase3_path}")
        return 1
    out_dir.mkdir(parents=True, exist_ok=True)

    agg = load_aggregate(input_dir)
    early = agg["early_window_checkpoints"]
    receivers = agg["receivers"]
    rows = load_checkpoint_rows(input_dir)
    print(f"[load] {len(rows)} per-(receiver, seed) checkpoint files; "
          f"early window = first {len(early)} of "
          f"{len(agg['checkpoints'])} checkpoints")
    with open(phase3_path, "r", encoding="utf-8") as f:
        phase3 = json.load(f)

    regret = compute_omega_regret(rows, len(early), args.seed)
    matched = match_pairs(phase3, regret["pair"])
    corr = {
        "spearman": spearman([m["omega_total_interim"] for m in matched],
                             [m["omega_regret"] for m in matched]),
        "kendall": kendall([m["omega_total_interim"] for m in matched],
                           [m["omega_regret"] for m in matched]),
        "n_pairs": len(matched),
    }

    flat = [(d, r, c) for d, recs in regret["pair"].items()
            for r, c in recs.items()]
    flat.sort(key=lambda t: -t[2]["omega_regret"])
    top_cells = flat[:8]

    results = {
        "inputs": {
            "input_dir": str(input_dir),
            "phase3_omega": str(phase3_path),
            "n_checkpoint_files": len(rows),
            "early_window_checkpoints": early,
            "receivers": receivers,
            "bootstrap_seed": args.seed,
        },
        "definition": (
            "omega_regret(A->B) = mean over matched seeds of the "
            "average over early-window checkpoints of "
            "max(0, tau_cold(b) - tau_full(b)) under fixed-trust full "
            "reuse; omega_regret_final uses the final checkpoint only"),
        "engine_to_channel": ENGINE_TO_CHANNEL,
        "omega_regret": regret,
        "matched_pairs": matched,
        "correlation_interim_vs_regret": corr,
        "largest_regret_cells": [
            {"donor": d, "receiver": r, **c} for d, r, c in top_cells],
    }
    with open(out_dir / "omega_regret.json", "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)

    note = (f"Inputs: {len(rows)} (receiver, seed) checkpoint files from "
            f"{input_dir.name} (matched seeds, fixed-trust full reuse); "
            f"early window = first {len(early)} of "
            f"{len(agg['checkpoints'])} checkpoints; interim omega from "
            f"{phase3_path.name}. Bootstrap CIs over seeds.")
    write_markdown(out_dir / "omega_regret.md", note, matched, corr,
                   regret, top_cells)
    make_plot(matched, regret["pair"], receivers,
              out_dir / "omega_regret.png")

    print(f"[save] {out_dir}")
    print(f"  interim vs regret (12 matched pairs): Spearman "
          f"{corr['spearman']:+.4f}, Kendall {corr['kendall']:+.4f}")
    for d, r, c in top_cells[:5]:
        print(f"  top regret: {d} -> {r}: {c['omega_regret']:.4f} "
              f"[{c['omega_regret_ci'][0]:.4f}, "
              f"{c['omega_regret_ci'][1]:.4f}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
