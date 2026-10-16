"""
Minimum bridge overlap for a stable cross-channel offset
(the certificate-thinness curve).

Motivation (PEGKi limit -> prescription):
  The recorded PEGKi limit says thin cross-channel overlap means no transfer
  certificate -- exactly why Al-Mg-Si failed the retrospective gate
  (compute/almgsi_archive/results/TRACKA_GATE.md): only 3 shared ternary compositions
  between aflow_dHf and mp_dHf, so the ternary-only Kendall tau of 1.0 at n=3
  certifies nothing. This experiment turns that limit into a PRESCRIPTION.
  Using channel pairs of the frozen band-gap universe where the shared-
  composition overlap is rich, we subsample the overlap down to
  Al-Mg-Si-like thinness and measure how the estimated cross-channel offset
  (its sign, its magnitude) and the shared-set rank agreement destabilise.
  The output is a bridge-budget curve: how many shared (bridge) compositions
  a prospective campaign must buy before a transfer certificate is possible.

Data:
  shared/output/bandgap_archive_build/bandgap_universe.jsonl (o-schema JSONL:
  x.formula, s = channel, y.bandgap_eV). Records with gap <= 0.01 eV are
  dropped, mirroring MIN_GAP_EV in scripts/hypergraph_trial_bandgap.py
  (metals are all exactly 0 and carry no ordering signal).

Method per channel pair (A, B):
  1. Formula matching: hypergraph_trial_bandgap.py matches RAW formula
     strings across channels; here we normalise each formula to a canonical
     reduced form instead (parse element counts incl. parentheses and
     fractional stoichiometries via exact Fractions, divide by the GCD,
     emit elements alphabetically), so e.g. "Ba(GaO2)2" and "BaGa2O4"
     match. Unparseable formulas are skipped and counted.
  2. Per (channel, formula) with several records, the MEDIAN gap is used
     (hypergraph_trial_bandgap.py keeps the first record per pair; the
     median is the more defensible single value and is documented here).
  3. GROUND TRUTH: every shared formula contributes one 2-entity bridge
     event ranking channel A's median value against channel B's (exact ties
     dropped by hyperedge_from_scores, the package draw-drop rule). The two
     CHANNELS are the entities of a 2-entity Plackett-Luce fit
     (PlackettLuceModel; for 2-entity edges this is exactly Bradley-Terry).
     The log-strength difference IS the cross-channel offset:
         delta_AB = theta_A - theta_B
     delta_AB > 0 means channel A's gaps run systematically HIGHER than
     channel B's over shared formulas. This 2-entity reduction is chosen
     over reusing PooledPlackettLuce with class-level pooling because the
     offset between exactly two channels needs no policy factor -- the
     policy (the formula) is shared within every event and cancels by
     construction; the reduction is simpler and exactly equivalent for
     this quantity. Full-overlap Kendall tau between the channels' values
     over shared formulas is computed alongside.
  4. THINNING: for each n in the grid (default 3,5,8,12,20,35,60,100,200,
     capped at the available overlap; n=3 is exactly the Al-Mg-Si
     situation), draw R random subsamples of n shared formulas (seeded,
     without replacement) and refit delta_n and tau_n per subsample.
     Reported vs n:
       (a) sign stability      P(sign(delta_n) == sign(delta_full))
       (b) median |delta_n - delta_full|
       (c) tau power           fraction of subsamples with two-sided
                               kendalltau p < alpha (scipy.stats.kendalltau,
                               method="auto": exact permutation null for
                               small n without ties, normal approximation
                               otherwise) -- "tau CI excludes 0"
       (d) tau CI half-width   analytic null normal approximation
                               1.96 * sqrt(2*(2n+5) / (9*n*(n-1)))
                               (a function of n only), plus the empirical
                               2.5-97.5 percentile spread of tau_n.
  5. CERTIFICATE THRESHOLD: n* = smallest tested n with
     sign stability >= 0.95 AND tau power >= 0.8. These thresholds are
     PLANNING DEFAULTS (overridable via --sign-stability / --tau-power),
     not laws.

Honest caveat: band-gap channel pairs are a PROXY for Al-Mg-Si channel
pairs (formation-energy DFT databases, empirical property channels). The
transferable content is the SHAPE of the curve and the ORDER OF MAGNITUDE
of n*, not the exact number.

Known estimator artifact, stated up front: when one channel wins (almost)
every bridge event -- complete separation, e.g. HSE vs GGA of the same
database -- the Bradley-Terry MLE diverges and the fitted |delta| is
limited by the MAP damping of PlackettLuceModel (roughly
|delta| <= log((wins + prior)/prior)), so the median-absolute-offset-error
column can be NON-monotone in n for such pairs. The certificate criteria
(sign stability and tau power) are unaffected by this artifact.

Outputs (to --output-dir, default shared/output/bridge_budget/):
  bridge_budget.json, bridge_budget.md, bridge_budget.png

Usage (from shared/):
  uv run python scripts/bridge_budget.py --self-test
  uv run python scripts/bridge_budget.py --smoke --workers 8
  uv run python scripts/bridge_budget.py --workers 8

Seeded and deterministic given --seed (per-(pair, n) SeedSequence spawning,
independent of worker scheduling). ASCII only.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import argparse
import json
import multiprocessing
import os
import re
import tempfile
from fractions import Fraction
from functools import reduce
from math import gcd

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

from src.rating.hypergraph_rating import (
    PlackettLuceModel,
    hyperedge_from_scores,
)

MIN_GAP_EV = 0.01  # mirrors MIN_GAP_EV in scripts/hypergraph_trial_bandgap.py
DEFAULT_N_GRID = (3, 5, 8, 12, 20, 35, 60, 100, 200)
SMOKE_N_GRID = (3, 20)
DEFAULT_PAIRS = (
    ("jarvis_optb88vdw", "jarvis_mbj"),   # same-database fidelity pair
    ("snumat_gga", "snumat_hse"),         # same-database GGA vs hybrid
    ("mp_gga", "expt_zhuo"),              # DFT vs experiment (prospective-like)
)
ALPHA = 0.05

_ELEM_RE = re.compile(r"([A-Z][a-z]?)(\d*(?:\.\d+)?)")


# ------------------------------------------------------------------
# Formula normalisation
# ------------------------------------------------------------------

def normalise_formula(formula):
    """Canonical reduced formula: element counts reduced, alphabetical.

    Handles parentheses with integer multipliers and fractional
    stoichiometries (exact via Fraction). Returns None if unparseable.
    Examples: "Ba(GaO2)2" -> "BaGa2O4"; "Sn0.05Pb0.95Se1" -> "PbSe20Sn"
    is WRONG -- fractions are scaled exactly: "Sn0.05Pb0.95Se1" ->
    counts {Sn:1, Pb:19, Se:20} -> "Pb19Se20Sn".
    """
    stack = [{}]
    i = 0
    n = len(formula)
    while i < n:
        ch = formula[i]
        if ch == "(":
            stack.append({})
            i += 1
        elif ch == ")":
            if len(stack) < 2:
                return None
            i += 1
            m = re.match(r"\d+(?:\.\d+)?", formula[i:])
            mult = Fraction(m.group()) if m else Fraction(1)
            if m:
                i += m.end()
            top = stack.pop()
            for el, c in top.items():
                stack[-1][el] = stack[-1].get(el, Fraction(0)) + c * mult
        else:
            m = _ELEM_RE.match(formula, i)
            if not m or m.end() == i:
                return None
            el = m.group(1)
            c = Fraction(m.group(2)) if m.group(2) else Fraction(1)
            stack[-1][el] = stack[-1].get(el, Fraction(0)) + c
            i = m.end()
    if len(stack) != 1 or not stack[0]:
        return None
    counts = stack[0]
    if any(c <= 0 for c in counts.values()):
        return None
    # Scale fractions to exact integers, then reduce by the GCD.
    denom_lcm = 1
    for c in counts.values():
        denom_lcm = denom_lcm * c.denominator // gcd(denom_lcm, c.denominator)
    ints = {el: int(c * denom_lcm) for el, c in counts.items()}
    g = reduce(gcd, ints.values())
    parts = []
    for el in sorted(ints):
        k = ints[el] // g
        parts.append(el if k == 1 else f"{el}{k}")
    return "".join(parts)


# ------------------------------------------------------------------
# Universe loading
# ------------------------------------------------------------------

def load_channel_values(universe_path, channels):
    """-> values[channel][canonical_formula] = median gap; plus counters.

    Drops records with gap <= MIN_GAP_EV (hypergraph_trial_bandgap.py
    convention). Median over multiple records per (channel, formula).
    """
    raw = {c: {} for c in channels}
    n_seen = n_used = n_unparsed = 0
    with open(universe_path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            n_seen += 1
            channel = rec.get("s")
            if channel not in raw:
                continue
            formula = (rec.get("x") or {}).get("formula")
            gap = (rec.get("y") or {}).get("bandgap_eV")
            if not formula or gap is None:
                continue
            gap = float(gap)
            if gap <= MIN_GAP_EV:
                continue
            key = normalise_formula(formula)
            if key is None:
                n_unparsed += 1
                continue
            raw[channel].setdefault(key, []).append(gap)
            n_used += 1
    values = {
        c: {f: float(np.median(v)) for f, v in raw[c].items()}
        for c in channels
    }
    return values, n_seen, n_used, n_unparsed


# ------------------------------------------------------------------
# Offset fit (2-entity Bradley-Terry via PlackettLuceModel)
# ------------------------------------------------------------------

def fit_offset(vals_a, vals_b, name_a, name_b):
    """delta_AB = theta_A - theta_B from paired values (one event each).

    delta_AB > 0: channel A's values run systematically higher. Exact ties
    are dropped (draw-drop rule). Returns None if no untied event remains.
    """
    edges = []
    for va, vb in zip(vals_a, vals_b):
        e = hyperedge_from_scores({name_a: float(va), name_b: float(vb)},
                                  higher_is_better=True)
        if e is not None:
            edges.append(e)
    if not edges:
        return None
    model = PlackettLuceModel().fit(edges)
    return float(model.log_strengths[name_a] - model.log_strengths[name_b])


def tau_ci_half_width_null(n):
    """1.96 * null-hypothesis normal-approx sd of Kendall tau at size n."""
    if n < 2:
        return float("nan")
    return float(1.96 * np.sqrt(2.0 * (2 * n + 5) / (9.0 * n * (n - 1))))


# ------------------------------------------------------------------
# Thinning worker (parallelised over (pair, n))
# ------------------------------------------------------------------

def thin_task(task):
    """Evaluate one (pair, n) cell: R seeded subsamples of n shared formulas."""
    a = np.asarray(task["a"], dtype=float)
    b = np.asarray(task["b"], dtype=float)
    n = task["n"]
    R = task["resamples"]
    alpha = task["alpha"]
    sign_full = task["sign_full"]
    delta_full = task["delta_full"]
    name_a, name_b = task["name_a"], task["name_b"]
    rng = np.random.default_rng(
        np.random.SeedSequence(entropy=(task["seed"], task["pair_index"], n)))

    n_shared = len(a)
    sign_hits = 0
    abs_errors = []
    taus = []
    power_hits = 0
    n_degenerate = 0
    deltas = []
    for _ in range(R):
        idx = rng.choice(n_shared, size=n, replace=False)
        sa, sb = a[idx], b[idx]
        delta_n = fit_offset(sa, sb, name_a, name_b)
        if delta_n is None:
            n_degenerate += 1  # all events tied: no offset estimate
        else:
            deltas.append(delta_n)
            abs_errors.append(abs(delta_n - delta_full))
            if np.sign(delta_n) == sign_full:
                sign_hits += 1
        res = stats.kendalltau(sa, sb)
        tau_n, p_n = float(res.statistic), float(res.pvalue)
        if np.isfinite(tau_n):
            taus.append(tau_n)
        if np.isfinite(p_n) and p_n < alpha:
            power_hits += 1

    taus = np.array(taus) if taus else np.array([np.nan])
    deltas = np.array(deltas) if deltas else np.array([np.nan])
    return {
        "pair": task["pair"],
        "n": n,
        "resamples": R,
        "sign_stability": sign_hits / R,
        "median_abs_offset_error": (float(np.median(abs_errors))
                                    if abs_errors else float("nan")),
        "tau_power": power_hits / R,
        "tau_ci_half_width_null": tau_ci_half_width_null(n),
        "tau_median": float(np.nanmedian(taus)),
        "tau_p2_5": float(np.nanpercentile(taus, 2.5)),
        "tau_p97_5": float(np.nanpercentile(taus, 97.5)),
        "delta_median": float(np.nanmedian(deltas)),
        "delta_p2_5": float(np.nanpercentile(deltas, 2.5)),
        "delta_p97_5": float(np.nanpercentile(deltas, 97.5)),
        "n_degenerate": n_degenerate,
    }


# ------------------------------------------------------------------
# Per-pair pipeline
# ------------------------------------------------------------------

def build_pair_ground_truth(values, chan_a, chan_b):
    """Shared formulas (sorted), aligned median values, full-overlap fit."""
    shared = sorted(set(values[chan_a]) & set(values[chan_b]))
    a = np.array([values[chan_a][f] for f in shared], dtype=float)
    b = np.array([values[chan_b][f] for f in shared], dtype=float)
    delta_full = fit_offset(a, b, chan_a, chan_b) if len(shared) >= 2 else None
    if len(shared) >= 2:
        res = stats.kendalltau(a, b)
        tau_full, tau_full_p = float(res.statistic), float(res.pvalue)
    else:
        tau_full, tau_full_p = float("nan"), float("nan")
    return shared, a, b, delta_full, tau_full, tau_full_p


def find_n_star(curve, sign_thresh, power_thresh):
    for row in sorted(curve, key=lambda r: r["n"]):
        if (row["sign_stability"] >= sign_thresh
                and row["tau_power"] >= power_thresh):
            return row["n"]
    return None


# ------------------------------------------------------------------
# Outputs
# ------------------------------------------------------------------

def atomic_write(path, text):
    fd, tmp = tempfile.mkstemp(dir=str(path.parent), suffix=".tmp")
    with os.fdopen(fd, "w", encoding="utf-8") as fh:
        fh.write(text)
    os.replace(tmp, str(path))


def write_plot(pairs_out, out_path, sign_thresh, power_thresh):
    n_panels = max(len(pairs_out), 1)
    fig, axes = plt.subplots(1, n_panels, figsize=(4.6 * n_panels, 4.0),
                             squeeze=False)
    for ax, po in zip(axes[0], pairs_out):
        curve = sorted(po["curve"], key=lambda r: r["n"])
        ns = [r["n"] for r in curve]
        ax.plot(ns, [r["sign_stability"] for r in curve], "o-",
                color="#1f6f8b", label="sign stability")
        ax.plot(ns, [r["tau_power"] for r in curve], "s-",
                color="#c86b3c", label="tau power (p < 0.05)")
        ax.axhline(sign_thresh, color="#1f6f8b", ls=":", lw=1)
        ax.axhline(power_thresh, color="#c86b3c", ls=":", lw=1)
        ax.axvline(3, color="#a02c2c", ls="--", lw=1.2)
        ax.text(3, 0.02, " n=3 (Al-Mg-Si gate)", color="#a02c2c",
                fontsize=8, rotation=90, va="bottom")
        if po["n_star"] is not None:
            ax.axvline(po["n_star"], color="#2c7a2c", ls="--", lw=1.2)
            ax.text(po["n_star"], 0.02, f" n*={po['n_star']}",
                    color="#2c7a2c", fontsize=8, rotation=90, va="bottom")
        ax.set_xscale("log")
        ax.set_ylim(-0.03, 1.05)
        ax.set_xlabel("bridge compositions n (log)")
        ax.set_title(f"{po['channel_a']} vs {po['channel_b']}\n"
                     f"overlap={po['n_shared']}, "
                     f"delta_full={po['delta_full']:+.3f}, "
                     f"tau_full={po['tau_full']:.3f}", fontsize=9)
        ax.legend(fontsize=8, loc="center right")
    axes[0][0].set_ylabel("fraction of subsamples")
    fig.suptitle("Bridge budget: offset sign stability and tau power vs "
                 "shared-composition count", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def write_markdown(results, out_path):
    cfg = results["config"]
    lines = [
        "# Bridge budget: minimum shared-composition overlap for a "
        "stable cross-channel offset",
        "",
        f"Universe: `{cfg['universe']}` "
        f"({results['n_seen']} records seen, {results['n_used']} used with "
        f"gap > {MIN_GAP_EV} eV in the requested channels, "
        f"{results['n_unparsed']} unparseable formulas skipped).",
        "",
        "Per pair, ALL shared (canonical reduced) formulas define the "
        "ground-truth cross-channel offset delta_full (2-entity "
        "Bradley-Terry log-strength difference; delta_AB > 0 means channel "
        "A's gaps run higher) and full-overlap Kendall tau. The overlap is "
        f"then thinned to n compositions ({cfg['resamples']} seeded "
        "subsamples per n) and the offset and tau are refitted.",
        "",
        "Tau significance per subsample: scipy.stats.kendalltau two-sided "
        f"p < {cfg['alpha']} (exact permutation null at small n without "
        "ties, normal approximation otherwise). The tau CI half-width "
        "column is the analytic null normal approximation "
        "1.96*sqrt(2(2n+5)/(9n(n-1))).",
        "",
        f"Certificate threshold (planning defaults): n* = smallest n with "
        f"sign stability >= {cfg['sign_stability']} AND tau power >= "
        f"{cfg['tau_power']}.",
        "",
        "Estimator note: under (near-)complete separation (one channel "
        "wins almost every bridge event, e.g. hybrid vs GGA of the same "
        "database) the Bradley-Terry MLE diverges and the fitted |delta| "
        "is limited by the MAP damping of PlackettLuceModel, so the "
        "median-absolute-offset-error column can be non-monotone in n for "
        "such pairs. Sign stability and tau power are unaffected.",
        "",
    ]
    for po in results["pairs"]:
        lines += [
            f"## {po['channel_a']} vs {po['channel_b']}",
            "",
            f"- records used: {po['n_records_a']} ({po['channel_a']}), "
            f"{po['n_records_b']} ({po['channel_b']}); "
            f"shared formulas: **{po['n_shared']}**",
            f"- delta_full = {po['delta_full']:+.4f} "
            f"(sign: {'A higher' if po['delta_full'] > 0 else 'B higher'}), "
            f"tau_full = {po['tau_full']:.4f} (p = {po['tau_full_p']:.2e})",
            f"- n* = {po['n_star'] if po['n_star'] is not None else 'not reached in tested grid'}",
            "",
            "| n | sign stability | median abs offset error | tau power | "
            "tau CI half-width (null) | tau median [2.5, 97.5 pct] |",
            "|---|---|---|---|---|---|",
        ]
        for row in sorted(po["curve"], key=lambda r: r["n"]):
            lines.append(
                f"| {row['n']} | {row['sign_stability']:.3f} | "
                f"{row['median_abs_offset_error']:.4f} | "
                f"{row['tau_power']:.3f} | "
                f"{row['tau_ci_half_width_null']:.3f} | "
                f"{row['tau_median']:.3f} "
                f"[{row['tau_p2_5']:.3f}, {row['tau_p97_5']:.3f}] |")
        lines.append("")

    lines += ["## Bridge-budget table (n* per pair)", "",
              "| pair | overlap | delta_full | tau_full | n* | n=3 status |",
              "|---|---|---|---|---|---|"]
    for po in results["pairs"]:
        row3 = next((r for r in po["curve"] if r["n"] == 3), None)
        if row3 is None:
            status3 = "n=3 not tested"
        else:
            status3 = (f"sign stab {row3['sign_stability']:.2f}, "
                       f"tau power {row3['tau_power']:.2f}")
        nstar = po["n_star"] if po["n_star"] is not None else "> max tested"
        lines.append(f"| {po['channel_a']} vs {po['channel_b']} | "
                     f"{po['n_shared']} | {po['delta_full']:+.3f} | "
                     f"{po['tau_full']:.3f} | {nstar} | {status3} |")
    lines += [
        "",
        "## What n=3 means for the Al-Mg-Si archive gate",
        "",
        "The Al-Mg-Si retrospective gate (compute/almgsi_archive/results/"
        "TRACKA_GATE.md) found exactly 3 shared ternary compositions "
        "between aflow_dHf and mp_dHf, with ternary-only Kendall tau = 1.0 "
        "at n=3. The thinning curves above show what n=3 buys in channel "
        "pairs where the truth is KNOWN from a rich overlap: the exact "
        "permutation test at n=3 can never reach p < 0.05 (its smallest "
        "two-sided p-value is 1/3), so the tau power at n=3 is identically "
        "zero -- a perfect tau over 3 compositions certifies nothing -- and "
        "the offset sign estimated from 3 bridge events flips across "
        "subsamples unless the channel separation is very large. The n* "
        "column is the bridge budget a prospective Al-Mg-Si campaign must "
        "purchase (shared compositions measured/computed in BOTH channels) "
        "before a transfer certificate is possible.",
        "",
        "Caveat, stated honestly: these are band-gap channel pairs, a proxy "
        "for the Al-Mg-Si formation-energy/empirical channel pairs. The "
        "transferable content is the SHAPE of the destabilisation curve "
        "and the ORDER OF MAGNITUDE of n*, not the exact number.",
        "",
    ]
    atomic_write(out_path, "\n".join(lines) + "\n")


# ------------------------------------------------------------------
# Self-test (no universe needed)
# ------------------------------------------------------------------

def _synthetic_pair(rng, n_formulas, shift, noise, rho_break=0.0):
    """Two channels over shared synthetic formulas with known offset.

    channel B = latent + shift + noise; channel A = latent + noise.
    shift > 0 makes B systematically higher, so delta_AB < 0.
    rho_break > 0 replaces B's latent with an independent draw with that
    probability scale (rho_break=1: channels are unrelated -> tau near 0).
    """
    latent = rng.normal(0.0, 1.0, size=n_formulas)
    latent_b = latent.copy()
    if rho_break > 0:
        latent_b = (1.0 - rho_break) * latent + rho_break * rng.normal(
            0.0, 1.0, size=n_formulas)
    a = latent + rng.normal(0.0, noise, size=n_formulas)
    b = latent_b + shift + rng.normal(0.0, noise, size=n_formulas)
    return a, b


def _run_thin(a, b, n, R, seed, alpha=ALPHA):
    delta_full = fit_offset(a, b, "chan_a", "chan_b")
    task = {
        "pair": "selftest", "pair_index": 0, "n": n, "resamples": R,
        "alpha": alpha, "seed": seed, "a": a, "b": b,
        "delta_full": delta_full, "sign_full": float(np.sign(delta_full)),
        "name_a": "chan_a", "name_b": "chan_b",
    }
    return thin_task(task), delta_full


def _self_test():
    print("bridge_budget self-test")
    rng = np.random.default_rng(12345)
    failures = 0

    # 0. Formula normalisation
    cases = [("Ba(GaO2)2", "BaGa2O4"), ("BaGa2O4", "BaGa2O4"),
             ("Sn2O4", "O2Sn"), ("Sn0.05Pb0.95Se1", "Pb19Se20Sn"),
             ("Sr2Cu(ClO)2", "Cl2CuO2Sr2")]
    ok0 = True
    for raw, want in cases:
        got = normalise_formula(raw)
        if got != want:
            ok0 = False
            print(f"  [0] FAIL: {raw} -> {got}, expected {want}")
    print(f"  [0] formula normalisation: {'ok' if ok0 else 'FAIL'}")
    if not ok0:
        failures += 1

    # 1. Offset sign recovery: B higher by construction -> delta_AB < 0
    a, b = _synthetic_pair(rng, 500, shift=0.6, noise=0.5)
    d_ab = fit_offset(a, b, "chan_a", "chan_b")
    print(f"  [1] delta_AB with B higher: {d_ab:+.3f} (expect < 0)")
    if not (d_ab is not None and d_ab < 0):
        failures += 1
        print("  [1] FAIL")

    # 2. Orientation both ways (sign rule): swapping channels flips the
    #    sign exactly; reversing the shift flips it too.
    d_ba = fit_offset(b, a, "chan_b", "chan_a")
    a2, b2 = _synthetic_pair(rng, 500, shift=-0.6, noise=0.5)
    d_neg = fit_offset(a2, b2, "chan_a", "chan_b")
    ok2 = (abs(d_ab + d_ba) < 1e-9) and (d_neg > 0)
    print(f"  [2] orientation: delta_BA = {d_ba:+.3f} "
          f"(= -delta_AB: {abs(d_ab + d_ba) < 1e-9}), "
          f"reversed-shift delta_AB = {d_neg:+.3f} (expect > 0)")
    if not ok2:
        failures += 1
        print("  [2] FAIL")

    # 3. Magnitude direction: larger true shift -> larger |delta|
    a3, b3 = _synthetic_pair(rng, 800, shift=0.3, noise=0.5)
    a4, b4 = _synthetic_pair(rng, 800, shift=1.2, noise=0.5)
    d_small = fit_offset(a3, b3, "chan_a", "chan_b")
    d_large = fit_offset(a4, b4, "chan_a", "chan_b")
    print(f"  [3] |delta|(shift=0.3) = {abs(d_small):.3f}, "
          f"|delta|(shift=1.2) = {abs(d_large):.3f}")
    if not abs(d_large) > abs(d_small):
        failures += 1
        print("  [3] FAIL: magnitude not monotone in true shift")

    # 4. Thinning converges to truth at large n and is unstable at n=3
    a5, b5 = _synthetic_pair(rng, 500, shift=0.25, noise=0.6)
    row3, d_full5 = _run_thin(a5, b5, 3, 300, seed=99)
    row200, _ = _run_thin(a5, b5, 200, 300, seed=99)
    ok4 = (row200["sign_stability"] >= 0.95
           and row200["sign_stability"] >= row3["sign_stability"]
           and row200["median_abs_offset_error"]
           < row3["median_abs_offset_error"])
    print(f"  [4] thinning (delta_full={d_full5:+.3f}): "
          f"sign stab n=3 {row3['sign_stability']:.2f} vs "
          f"n=200 {row200['sign_stability']:.2f}; "
          f"med abs err n=3 {row3['median_abs_offset_error']:.3f} vs "
          f"n=200 {row200['median_abs_offset_error']:.3f}")
    if not ok4:
        failures += 1
        print("  [4] FAIL: large-n thinning does not converge to truth")

    # 5. Tau power: correlated channels gain power with n; n=3 is exactly
    #    powerless (exact test cannot reach p < 0.05); unrelated channels
    #    stay near the alpha level.
    row60, _ = _run_thin(a5, b5, 60, 300, seed=7)
    au, bu = _synthetic_pair(rng, 500, shift=0.0, noise=0.5, rho_break=1.0)
    rowu, _ = _run_thin(au, bu, 60, 300, seed=7)
    tau_full5 = float(stats.kendalltau(a5, b5).statistic)
    ok5 = (row3["tau_power"] == 0.0 and row60["tau_power"] >= 0.8
           and rowu["tau_power"] <= 0.25 and tau_full5 > 0.3)
    print(f"  [5] tau power: n=3 {row3['tau_power']:.2f} (expect 0), "
          f"n=60 correlated {row60['tau_power']:.2f} (expect >= 0.8), "
          f"n=60 unrelated {rowu['tau_power']:.2f} (expect <= 0.25); "
          f"tau_full = {tau_full5:.3f}")
    if not ok5:
        failures += 1
        print("  [5] FAIL")

    # 6. Known perfect tau: zero noise, zero shift-noise -> tau_full = 1
    lat = np.sort(rng.normal(0.0, 1.0, size=50))
    tau_perfect = float(stats.kendalltau(lat, lat + 0.5).statistic)
    print(f"  [6] noiseless tau_full = {tau_perfect:.3f} (expect 1.0)")
    if abs(tau_perfect - 1.0) > 1e-12:
        failures += 1
        print("  [6] FAIL")

    # 7. Determinism: same seed -> identical thinning row
    row_a, _ = _run_thin(a5, b5, 12, 50, seed=4242)
    row_b, _ = _run_thin(a5, b5, 12, 50, seed=4242)
    ok7 = row_a == row_b
    print(f"  [7] deterministic reruns identical: {ok7}")
    if not ok7:
        failures += 1
        print("  [7] FAIL")

    print(f"self-test {'PASS' if failures == 0 else 'FAIL'} "
          f"({failures} failures)")
    return failures


# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------

def parse_pairs(spec):
    pairs = []
    for item in spec.split(","):
        item = item.strip()
        if not item:
            continue
        if ":" not in item:
            raise ValueError(f"pair must be 'chanA:chanB', got {item!r}")
        a, b = item.split(":", 1)
        pairs.append((a.strip(), b.strip()))
    if not pairs:
        raise ValueError("no pairs parsed")
    return pairs


def main():
    ap = argparse.ArgumentParser(
        description="Minimum bridge overlap for a stable cross-channel "
                    "offset (certificate-thinness curve)")
    ap.add_argument("--universe",
                    default="output/bandgap_archive_build/bandgap_universe.jsonl",
                    help="JSONL universe, relative to shared/ (or absolute)")
    ap.add_argument("--pairs",
                    default=",".join(f"{a}:{b}" for a, b in DEFAULT_PAIRS),
                    help="comma-separated chanA:chanB pairs")
    ap.add_argument("--resamples", type=int, default=300)
    ap.add_argument("--seed", type=int, default=20260926)
    ap.add_argument("--output-dir", default="output/bridge_budget")
    ap.add_argument("--workers", type=int, default=1,
                    help="parallelise over (pair, n) cells")
    ap.add_argument("--alpha", type=float, default=ALPHA,
                    help="two-sided significance level for tau power")
    ap.add_argument("--sign-stability", type=float, default=0.95,
                    help="n* threshold on sign stability (planning default)")
    ap.add_argument("--tau-power", type=float, default=0.8,
                    help="n* threshold on tau power (planning default)")
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--smoke", action="store_true",
                    help="R=30, two n values (3 and 20)")
    args = ap.parse_args()

    if args.self_test:
        return _self_test()

    root = Path(__file__).parent.parent
    universe = Path(args.universe)
    if not universe.is_absolute():
        universe = (root / universe).resolve()
    if not universe.exists():
        print(f"[error] universe not found: {universe}; "
              "run compute/bandgap_archive/build_inventory.py first")
        return 1
    out_dir = Path(args.output_dir)
    if not out_dir.is_absolute():
        out_dir = (root / out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    pairs = parse_pairs(args.pairs)
    resamples = 30 if args.smoke else args.resamples
    n_grid = SMOKE_N_GRID if args.smoke else DEFAULT_N_GRID
    channels = sorted({c for p in pairs for c in p})

    values, n_seen, n_used, n_unparsed = load_channel_values(
        universe, channels)
    print(f"[load] {n_seen} records seen; {n_used} used "
          f"(gap > {MIN_GAP_EV} eV, requested channels); "
          f"{n_unparsed} unparseable formulas skipped")
    for c in channels:
        print(f"[load] {c}: {len(values[c])} distinct reduced formulas")

    pairs_out = []
    tasks = []
    for pair_index, (chan_a, chan_b) in enumerate(pairs):
        if not values.get(chan_a) or not values.get(chan_b):
            print(f"[skip] {chan_a} vs {chan_b}: channel missing from "
                  "universe")
            continue
        shared, a, b, delta_full, tau_full, tau_full_p = \
            build_pair_ground_truth(values, chan_a, chan_b)
        if len(shared) < 3 or delta_full is None:
            print(f"[skip] {chan_a} vs {chan_b}: overlap {len(shared)} "
                  "too thin to define a ground truth (this IS the "
                  "Al-Mg-Si situation)")
            continue
        grid = sorted({n for n in n_grid if n <= len(shared)})
        print(f"[pair] {chan_a} vs {chan_b}: overlap = {len(shared)}, "
              f"delta_full = {delta_full:+.4f}, tau_full = {tau_full:.4f} "
              f"(p = {tau_full_p:.2e}), n grid = {grid}")
        po = {
            "channel_a": chan_a, "channel_b": chan_b,
            "n_records_a": len(values[chan_a]),
            "n_records_b": len(values[chan_b]),
            "n_shared": len(shared),
            "delta_full": delta_full,
            "delta_sign_convention":
                "delta_AB = theta_A - theta_B > 0 means channel A's gaps "
                "run systematically higher over shared formulas",
            "tau_full": tau_full, "tau_full_p": tau_full_p,
            "curve": [], "n_star": None,
        }
        pairs_out.append(po)
        for n in grid:
            tasks.append({
                "pair": f"{chan_a}:{chan_b}", "pair_index": pair_index,
                "n": n, "resamples": resamples, "alpha": args.alpha,
                "seed": args.seed, "a": a.tolist(), "b": b.tolist(),
                "delta_full": delta_full,
                "sign_full": float(np.sign(delta_full)),
                "name_a": chan_a, "name_b": chan_b,
            })

    if not pairs_out:
        print("[error] no pair with usable overlap")
        return 1

    if args.workers > 1:
        with multiprocessing.Pool(processes=args.workers) as pool:
            rows = pool.map(thin_task, tasks, chunksize=1)
    else:
        rows = [thin_task(t) for t in tasks]

    by_pair = {f"{po['channel_a']}:{po['channel_b']}": po
               for po in pairs_out}
    for row in rows:
        by_pair[row["pair"]]["curve"].append(row)
    for po in pairs_out:
        po["curve"].sort(key=lambda r: r["n"])
        po["n_star"] = find_n_star(po["curve"], args.sign_stability,
                                   args.tau_power)
        nstar = po["n_star"] if po["n_star"] is not None else "> max tested"
        print(f"[n*] {po['channel_a']} vs {po['channel_b']}: n* = {nstar} "
              f"(sign stab >= {args.sign_stability}, "
              f"tau power >= {args.tau_power})")
        for r in po["curve"]:
            print(f"     n={r['n']:>4d}  sign_stab={r['sign_stability']:.3f}"
                  f"  med|d_err|={r['median_abs_offset_error']:.4f}"
                  f"  tau_power={r['tau_power']:.3f}"
                  f"  tau_ci_hw={r['tau_ci_half_width_null']:.3f}")

    results = {
        "config": {
            "universe": str(universe), "pairs": [list(p) for p in pairs],
            "resamples": resamples, "seed": args.seed,
            "n_grid": list(n_grid), "alpha": args.alpha,
            "sign_stability": args.sign_stability,
            "tau_power": args.tau_power, "smoke": bool(args.smoke),
            "min_gap_eV": MIN_GAP_EV,
        },
        "n_seen": n_seen, "n_used": n_used, "n_unparsed": n_unparsed,
        "pairs": pairs_out,
    }
    atomic_write(out_dir / "bridge_budget.json",
                 json.dumps(results, indent=2, sort_keys=True))
    write_markdown(results, out_dir / "bridge_budget.md")
    write_plot(pairs_out, out_dir / "bridge_budget.png",
               args.sign_stability, args.tau_power)
    print(f"[done] wrote {out_dir / 'bridge_budget.json'}, .md, .png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
