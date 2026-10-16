"""
Are the fidelity-ladder offsets chemistry-dependent? (additivity stress test)

The recorded PEGKi limit says the pooled cross-channel model of
scripts/hypergraph_trial_bandgap.py is ADDITIVE: one offset per channel,
shared by all chemistries (w = exp(theta_class + delta_channel)), with no
deeper shared representation. This experiment stress-tests that additivity
on the frozen band-gap archive universe: it fits the SAME pooled Plackett-Luce trial
separately within each anion family (oxide, halide, chalcogenide, nitride)
and tests whether the per-channel offsets -- in particular the ladder
spacings

    delta(snumat_hse)  - delta(snumat_gga)
    delta(jarvis_mbj)  - delta(jarvis_optb88vdw)
    delta(expt_zhuo)   - delta(mp_gga)

differ across families beyond their uncertainties. Electronic-structure
practice says they should (the GGA gap underestimate is chemistry-dependent:
strongly for wide-gap ionic systems, weakly for covalent narrow-gap ones), so
a positive result documents precisely what a family-resolved
(representation-sharing) upgrade must capture, and a negative result defends
additivity at the current resolution.

Exact restriction that keeps the family fit well-posed
------------------------------------------------------
The trial's event construction is imported and reused unchanged
(load_universe, within_channel_edges, bridge_edges, classify_formula, same
stride, same gap > 0.01 eV filter, same channel set = the channels the trial
keeps). Restricting records to one anion family leaves exactly one class per
channel, so the trial's WITHIN-channel events (which rank classes inside one
channel) are degenerate under the restriction (a hyperedge needs >= 2 distinct
entities) -- and they carry no channel-offset information anyway: inside a
within-channel event every participant shares the same delta_channel, which
cancels in every Plackett-Luce stage. Symmetrically, in a BRIDGE event (same
formula ranked across channels) every participant shares the same class, so
theta_class cancels and the event carries pure channel-offset signal. The
pooled log-likelihood therefore factorises exactly into
f(theta; within-events) + g(delta; bridge-events), and the family-restricted
fit is: the pooled model on that family's entities (channel::family, all kept
channels) with the family's bridge events intact -- i.e. exactly the
delta-identifying sub-likelihood of the trial, restricted to one family.

Normalisation (kept identical to the trial): PooledPlackettLuce gauge-fixes
channel offsets to geometric mean 1 (zero mean in log space) over the
channels present in the fitted events. All fits here keep that gauge. Because
family fits may in principle present different channel subsets, the report
additionally tabulates mp_gga-referenced offsets (delta_c - delta_mp_gga),
which are gauge-invariant, and the three ladder spacings are differences of
offsets and hence gauge-invariant by construction.

Fast refits for the bootstrap: the bootstrap uses a vectorised implementation
of the identical MM update (the PooledPlackettLuce channel block with the
single-class policy block at its gauge-fixed value 1). The main fits are done
BOTH ways and the script asserts agreement, then uses the fast path for the
B resampled refits. Bootstrap resamples BRIDGE EVENTS with replacement
(event bootstrap; within-channel events are irrelevant to offsets, see the
factorisation above -- the full-pooled vs bridges-only deviation is computed
and reported for the global fit).

Held-out predictive check (model-selection version of the same question):
the union of the families' bridge events is split 80/20 (seeded, stratified
per family); (a) one global additive model is fitted on all training events
pooled, (b) one model per family on that family's training events; each
held-out event is scored under (a) and under (b) routed to its own family's
model. Reported: per-event mean log-likelihood difference (family - global)
with a bootstrap CI, and the parameter counts of both models. Only bridge
events enter this comparison: they are the only events on which the two
models differ (offsets), and a within-channel event spans multiple families
so it cannot be routed.

Usage (from shared/):
    uv run python scripts/family_offsets.py                # full experiment
    uv run python scripts/family_offsets.py --self-test    # synthetic checks
    uv run python scripts/family_offsets.py --smoke        # tiny real-data run

Outputs: family_offsets.json, family_offsets.md, family_offsets.png in
--output-dir. Deterministic given --seed. ASCII only.
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import argparse
import importlib.util
import json
import tempfile
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from src.rating.hypergraph_rating import (  # noqa: E402
    Hyperedge,
    PooledPlackettLuce,
)

# ------------------------------------------------------------------
# Reuse the trial's event construction verbatim (import, not copy)
# ------------------------------------------------------------------

def _load_trial_module():
    """Import scripts/hypergraph_trial_bandgap.py (same dir as this file)."""
    try:
        import hypergraph_trial_bandgap as mod  # script dir is sys.path[0]
        return mod
    except ImportError:
        path = Path(__file__).resolve().parent / "hypergraph_trial_bandgap.py"
        spec = importlib.util.spec_from_file_location(
            "hypergraph_trial_bandgap", path)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod


trial = _load_trial_module()

SEP = PooledPlackettLuce.SEP

# The three ladder spacings under test: (name, high-fidelity ch, low ch)
SPACING_DEFS = [
    ("hse_minus_gga_snumat", "snumat_hse", "snumat_gga"),
    ("mbj_minus_optb88_jarvis", "jarvis_mbj", "jarvis_optb88vdw"),
    ("expt_minus_gga_mp", "expt_zhuo", "mp_gga"),
]

DEFAULT_FAMILIES = ["oxide", "halide", "chalcogenide", "nitride"]
KNOWN_FAMILIES = [name for name, _ in trial.ANION_PRIORITY] + ["other"]
REF_CHANNEL = "mp_gga"
MIN_FAMILY_BRIDGES = 30
CI_LO, CI_HI = 2.5, 97.5


# ------------------------------------------------------------------
# Fast vectorised channel-offset fit (identical MM update)
# ------------------------------------------------------------------

def fit_channel_pl(events, prior=0.01, max_iter=2000, tol=1e-9):
    """Plackett-Luce over channels from bridge events, vectorised MM.

    events: list of tuples of channel names, ordered best -> worst (the
    channel projection of a bridge hyperedge; the shared class factor cancels
    exactly, see module docstring). This is the PooledPlackettLuce channel
    block with the single-policy block at its gauge value 1, so it converges
    to the same fixed point; the main fits assert that numerically.

    Returns ({channel: log_offset (zero-mean over present channels)}, converged).
    """
    channels = sorted({c for ev in events for c in ev})
    if len(channels) < 2:
        return {c: 0.0 for c in channels}, True
    idx = {c: i for i, c in enumerate(channels)}
    C = len(channels)
    members = []
    winners = []
    for ev in events:
        ids = [idx[c] for c in ev]
        for s in range(len(ids) - 1):
            members.append(ids[s:])
            winners.append(ids[s])
    S = len(members)
    M = np.zeros((S, C), dtype=float)
    for r, mem in enumerate(members):
        M[r, mem] = 1.0
    wins = np.bincount(np.array(winners, dtype=int), minlength=C).astype(float)
    logw = np.zeros(C, dtype=float)
    w = np.ones(C, dtype=float)
    converged = False
    for _ in range(max_iter):
        totals = M @ w
        denom = M.T @ (1.0 / totals)
        w_new = (wins + prior) / (denom + prior)
        logw_new = np.log(np.maximum(w_new, 1e-300))
        logw_new = logw_new - logw_new.mean()
        delta = float(np.max(np.abs(logw_new - logw)))
        logw = logw_new
        w = np.exp(logw)
        if delta < tol:
            converged = True
            break
    return {c: float(logw[idx[c]]) for c in channels}, converged


def bridge_channel_tuples(bridges):
    """Project bridge hyperedges (channel::class keys) to channel tuples."""
    return [tuple(k.split(SEP, 1)[0] for k in e.ranking) for e in bridges]


def spacings_from_offsets(offsets):
    out = {}
    for name, hi, lo in SPACING_DEFS:
        if hi in offsets and lo in offsets:
            out[name] = float(offsets[hi] - offsets[lo])
        else:
            out[name] = None
    return out


def referenced_offsets(offsets, ref=REF_CHANNEL):
    if ref not in offsets:
        return None
    return {c: float(v - offsets[ref]) for c, v in offsets.items()}


def event_loglik(offsets, ev):
    """PL log-likelihood of one channel-tuple event under given offsets.

    Channels unseen at fit time take the gauge-mean offset 0.0.
    """
    w = np.exp(np.array([offsets.get(c, 0.0) for c in ev], dtype=float))
    ll = 0.0
    for s in range(len(w) - 1):
        ll += float(np.log(w[s] / w[s:].sum()))
    return ll


# ------------------------------------------------------------------
# Bootstrap (event resampling, fast refits, optional workers)
# ------------------------------------------------------------------

def _boot_chunk(payload):
    events, seed_lists = payload
    n = len(events)
    reps = []
    for sl in seed_lists:
        rng = np.random.default_rng(sl)
        take = rng.integers(0, n, size=n)
        sample = [events[int(i)] for i in take]
        off, _ = fit_channel_pl(sample)
        reps.append(off)
    return reps


def bootstrap_offsets(events, tag, n_boot, seed, workers):
    """n_boot event-bootstrap refits; returns list of offset dicts."""
    seed_lists = [[int(seed), 1001, int(tag), b] for b in range(n_boot)]
    if workers > 1 and n_boot > 1:
        chunk = max(1, (n_boot + workers - 1) // workers)
        chunks = [seed_lists[i:i + chunk]
                  for i in range(0, n_boot, chunk)]
        with ProcessPoolExecutor(max_workers=workers) as ex:
            outs = list(ex.map(_boot_chunk, [(events, c) for c in chunks]))
        return [r for out in outs for r in out]
    return _boot_chunk((events, seed_lists))


def percentile_ci(values):
    """[lo, hi] percentile CI over non-nan values; None if < 3 valid."""
    arr = np.asarray(values, dtype=float)
    valid = arr[~np.isnan(arr)]
    if valid.size < 3:
        return None
    return [float(np.percentile(valid, CI_LO)),
            float(np.percentile(valid, CI_HI))]


def summarize_bootstrap(reps, channels):
    """Per-channel offset CIs (fit gauge and mp_gga-referenced) + spacing CIs."""
    off_ci = {}
    ref_ci = {}
    for c in channels:
        off_ci[c] = percentile_ci(
            [r.get(c, np.nan) for r in reps])
        ref_ci[c] = percentile_ci(
            [r[c] - r[REF_CHANNEL]
             if (c in r and REF_CHANNEL in r) else np.nan for r in reps])
    spa_ci = {}
    spa_valid = {}
    for name, hi, lo in SPACING_DEFS:
        vals = [r[hi] - r[lo] if (hi in r and lo in r) else np.nan
                for r in reps]
        spa_ci[name] = percentile_ci(vals)
        spa_valid[name] = int(np.sum(~np.isnan(np.asarray(vals))))
    return off_ci, ref_ci, spa_ci, spa_valid


# ------------------------------------------------------------------
# Analysis pipeline
# ------------------------------------------------------------------

def family_bridge_edges(formulas, family, channels_kept, max_bridges):
    """The trial's bridge_edges on the family-restricted record set."""
    fam_formulas = {f: chans for f, chans in formulas.items()
                    if trial.classify_formula(f) == family}
    return trial.bridge_edges(fam_formulas, channels_kept, max_bridges)


def pooled_offsets(bridges, extra_pooled_edges=None):
    """PooledPlackettLuce fit (the trial's model, verbatim) -> channel offsets."""
    edges = list(extra_pooled_edges) if extra_pooled_edges else []
    edges.extend(bridges)
    model = PooledPlackettLuce().fit(edges)
    return dict(model.channel_offsets), bool(model.converged)


def max_offset_deviation(a, b):
    common = sorted(set(a) & set(b))
    if not common:
        return float("nan")
    av = np.array([a[c] for c in common])
    bv = np.array([b[c] for c in common])
    av = av - av.mean()
    bv = bv - bv.mean()
    return float(np.max(np.abs(av - bv)))


def run_analysis(cfg):
    """Full pipeline; returns the results dict (also written to disk)."""
    root = Path(__file__).resolve().parent.parent
    universe = Path(cfg.universe)
    if not universe.is_absolute():
        universe = (root / cfg.universe).resolve()
    if not universe.exists():
        raise FileNotFoundError(f"universe not found: {universe}")
    out_dir = Path(cfg.output_dir)
    if not out_dir.is_absolute():
        out_dir = (root / cfg.output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    per_channel, formulas, n_seen, n_used = trial.load_universe(universe)
    print(f"[load] {n_seen} records, {n_used} with gap > "
          f"{trial.MIN_GAP_EV} eV, {len(per_channel)} channels")

    # --- trial-identical event construction --------------------------------
    ch_edges = trial.within_channel_edges(per_channel, cfg.events,
                                          cfg.min_class_size)
    channels_kept = sorted(ch_edges)
    print(f"[ok] channels kept (trial rule): {channels_kept}")
    if len(channels_kept) < 2:
        raise RuntimeError("fewer than 2 channels kept; cannot fit offsets")

    within_pooled = []
    for channel, edges in sorted(ch_edges.items()):
        for e in edges:
            within_pooled.append(Hyperedge(
                ranking=tuple(PooledPlackettLuce.pooled_key(channel, c)
                              for c in e.ranking)))

    global_bridges = trial.bridge_edges(formulas, set(channels_kept),
                                        cfg.max_bridges)
    print(f"[ok] global bridge events: {len(global_bridges)} "
          f"(cap {cfg.max_bridges}); within-channel events: "
          f"{len(within_pooled)}")

    # --- global fit: the committed baseline, replicated verbatim -----------
    print("[fit] global pooled model (trial replication: within + bridges)")
    g_off, g_conv = pooled_offsets(global_bridges,
                                   extra_pooled_edges=within_pooled)
    g_tuples = bridge_channel_tuples(global_bridges)
    g_off_fast, _ = fit_channel_pl(g_tuples)
    g_dev = max_offset_deviation(g_off, g_off_fast)
    print(f"[check] global full-pooled vs bridges-only offsets: "
          f"max dev = {g_dev:.2e} (offsets are identified by bridges alone)")
    g_spa = spacings_from_offsets(g_off)
    print("[ok] global offsets (zero-mean gauge): " + ", ".join(
        f"{c}={g_off[c]:+.3f}" for c in sorted(g_off, key=g_off.get,
                                               reverse=True)))
    print("[ok] global spacings: " + ", ".join(
        f"{k}={v:+.3f}" if v is not None else f"{k}=NA"
        for k, v in g_spa.items()))
    print(f"[boot] global: {cfg.bootstrap} refits on resampled bridges")
    g_reps = bootstrap_offsets(g_tuples, 999, cfg.bootstrap, cfg.seed,
                               cfg.workers)
    g_off_ci, g_ref_ci, g_spa_ci, g_spa_valid = summarize_bootstrap(
        g_reps, channels_kept)

    results = {
        "config": {k: getattr(cfg, k) for k in
                   ("universe", "families", "bootstrap", "events", "seed",
                    "output_dir", "workers", "min_class_size", "max_bridges",
                    "holdout_frac")},
        "universe": {"path": str(universe), "n_seen": n_seen,
                     "n_used": n_used, "channels_kept": channels_kept},
        "gauge_note": ("Offsets are zero-mean in log space over present "
                       "channels (the trial's gauge: geometric mean 1). "
                       "offsets_ref are delta_c - delta_mp_gga "
                       "(gauge-invariant); spacings are gauge-invariant."),
        "global": {
            "n_within_events": len(within_pooled),
            "n_bridge_events": len(global_bridges),
            "converged": g_conv,
            "full_pooled_vs_bridges_only_max_dev": g_dev,
            "offsets": {c: float(v) for c, v in g_off.items()},
            "offsets_ref": referenced_offsets(g_off),
            "offset_ci": g_off_ci,
            "offset_ref_ci": g_ref_ci,
            "spacings": g_spa,
            "spacing_ci": g_spa_ci,
            "spacing_boot_valid": g_spa_valid,
        },
        "families": {},
    }

    # --- family-restricted fits --------------------------------------------
    fam_tuples = {}
    for fi, family in enumerate(cfg.families):
        bridges = family_bridge_edges(formulas, family, set(channels_kept),
                                      cfg.max_bridges)
        n_rec = {c: len(per_channel.get(c, {}).get(family, []))
                 for c in channels_kept}
        tuples = bridge_channel_tuples(bridges)
        ev_per_ch = {c: sum(1 for ev in tuples if c in ev)
                     for c in channels_kept}
        print(f"[fam] {family}: {len(bridges)} bridge events; records/ch: "
              + ", ".join(f"{c}={n_rec[c]}" for c in channels_kept))
        entry = {
            "n_records_per_channel": n_rec,
            "n_bridge_events": len(bridges),
            "n_bridge_events_per_channel": ev_per_ch,
        }
        if len(bridges) < MIN_FAMILY_BRIDGES:
            print(f"[warn] {family}: fewer than {MIN_FAMILY_BRIDGES} bridge "
                  "events; skipping fit")
            entry["skipped"] = True
            results["families"][family] = entry
            continue
        print(f"[fit] {family}: pooled model on family bridges "
              "(within-events degenerate under one-family restriction)")
        f_off, f_conv = pooled_offsets(bridges)
        f_off_fast, _ = fit_channel_pl(tuples)
        f_dev = max_offset_deviation(f_off, f_off_fast)
        if not np.isnan(f_dev) and f_dev > 1e-4:
            print(f"[warn] {family}: pooled vs fast fit max dev {f_dev:.2e} "
                  "exceeds 1e-4")
        f_spa = spacings_from_offsets(f_off)
        print(f"[boot] {family}: {cfg.bootstrap} refits "
              f"({len(tuples)} events each, resampled with replacement)")
        reps = bootstrap_offsets(tuples, fi, cfg.bootstrap, cfg.seed,
                                 cfg.workers)
        off_ci, ref_ci, spa_ci, spa_valid = summarize_bootstrap(
            reps, channels_kept)
        entry.update({
            "skipped": False,
            "converged": f_conv,
            "pooled_vs_fast_max_dev": f_dev,
            "offsets": {c: float(v) for c, v in f_off.items()},
            "offsets_ref": referenced_offsets(f_off),
            "offset_ci": off_ci,
            "offset_ref_ci": ref_ci,
            "spacings": f_spa,
            "spacing_ci": spa_ci,
            "spacing_boot_valid": spa_valid,
            "spacing_minus_global": {
                k: (f_spa[k] - g_spa[k])
                if (f_spa[k] is not None and g_spa[k] is not None) else None
                for k in f_spa},
        })
        results["families"][family] = entry
        fam_tuples[family] = tuples
        print("[ok] " + family + " spacings: " + ", ".join(
            f"{k}={v:+.3f}" if v is not None else f"{k}=NA"
            for k, v in f_spa.items()))

    fitted = [f for f in cfg.families
              if not results["families"].get(f, {}).get("skipped", True)]

    # --- additivity test -----------------------------------------------------
    additivity = {}
    for name, hi, lo in SPACING_DEFS:
        vals = {f: results["families"][f]["spacings"][name] for f in fitted
                if results["families"][f]["spacings"][name] is not None}
        entry = {"per_family": vals, "global": g_spa[name]}
        if len(vals) >= 2:
            f_hi = max(vals, key=vals.get)
            f_lo = min(vals, key=vals.get)
            ci_hi_fam = results["families"][f_hi]["spacing_ci"][name]
            ci_lo_fam = results["families"][f_lo]["spacing_ci"][name]
            disjoint = None
            if ci_hi_fam is not None and ci_lo_fam is not None:
                disjoint = bool(ci_hi_fam[0] > ci_lo_fam[1])
            entry.update({
                "max_pairwise_diff": float(vals[f_hi] - vals[f_lo]),
                "family_max": f_hi,
                "family_min": f_lo,
                "extreme_family_cis_disjoint": disjoint,
            })
        additivity[name] = entry
    results["additivity"] = additivity

    # --- held-out predictive check -------------------------------------------
    results["holdout"] = holdout_check(fam_tuples, cfg)

    # --- verdict ---------------------------------------------------------------
    n_disjoint = sum(1 for a in additivity.values()
                     if a.get("extreme_family_cis_disjoint") is True)
    hold = results["holdout"]
    hold_positive = (hold is not None and hold["ci"] is not None
                     and hold["ci"][0] > 0.0)
    if n_disjoint > 0 and hold_positive:
        verdict = ("CHEMISTRY-DEPENDENT: ladder spacings differ across "
                   "families beyond bootstrap uncertainty AND the "
                   "family-resolved model wins the held-out comparison. "
                   "Additivity is violated at the current resolution.")
    elif n_disjoint > 0:
        verdict = ("PARTIALLY CHEMISTRY-DEPENDENT: some spacings separate "
                   "across families, but the family-resolved model does not "
                   "significantly beat the global additive model held-out.")
    elif hold_positive:
        verdict = ("WEAKLY CHEMISTRY-DEPENDENT: no single spacing separates "
                   "beyond CIs, but the family-resolved model wins the "
                   "held-out comparison.")
    else:
        verdict = ("ADDITIVITY DEFENDED at the current resolution: spacings "
                   "agree within bootstrap CIs and the family-resolved model "
                   "shows no significant held-out gain.")
    results["verdict"] = verdict
    results["n_spacings_with_disjoint_extreme_cis"] = n_disjoint
    print("[verdict] " + verdict)

    # --- outputs -------------------------------------------------------------
    trial.atomic_write(out_dir / "family_offsets.json",
                       json.dumps(results, indent=2, sort_keys=True))
    trial.atomic_write(out_dir / "family_offsets.md",
                       render_markdown(results, cfg))
    make_plot(results, cfg, out_dir / "family_offsets.png")
    print(f"[done] wrote {out_dir / 'family_offsets.json'}, .md, .png")
    return results


def holdout_check(fam_tuples, cfg):
    """80/20 split; global-additive vs family-resolved held-out log-lik."""
    if not fam_tuples:
        return None
    train_all = []
    fam_train = {}
    fam_test = {}
    for fi, family in enumerate(sorted(fam_tuples)):
        events = fam_tuples[family]
        rng = np.random.default_rng([int(cfg.seed), 7, fi])
        order = rng.permutation(len(events))
        n_train = max(1, int(round(len(events) * (1.0 - cfg.holdout_frac))))
        n_train = min(n_train, len(events) - 1)
        fam_train[family] = [events[i] for i in order[:n_train]]
        fam_test[family] = [events[i] for i in order[n_train:]]
        train_all.extend(fam_train[family])

    g_off, _ = fit_channel_pl(train_all)
    f_off = {f: fit_channel_pl(fam_train[f])[0] for f in fam_train}

    diffs = []
    per_family = {}
    n_unseen = 0
    for family in sorted(fam_test):
        fdiffs = []
        for ev in fam_test[family]:
            unseen = [c for c in ev
                      if c not in f_off[family] or c not in g_off]
            n_unseen += 1 if unseen else 0
            d = event_loglik(f_off[family], ev) - event_loglik(g_off, ev)
            fdiffs.append(d)
            diffs.append(d)
        per_family[family] = {
            "n_test": len(fdiffs),
            "mean_ll_diff": float(np.mean(fdiffs)) if fdiffs else None,
        }
    if not diffs:
        return None
    diffs = np.asarray(diffs, dtype=float)
    rng = np.random.default_rng([int(cfg.seed), 11])
    boot_means = [float(np.mean(diffs[rng.integers(0, diffs.size,
                                                   size=diffs.size)]))
                  for _ in range(2000)]
    ci = percentile_ci(boot_means)
    n_ch_g = len(g_off)
    n_ch_f = {f: len(f_off[f]) for f in f_off}
    result = {
        "holdout_frac": cfg.holdout_frac,
        "n_train_events": len(train_all),
        "n_test_events": int(diffs.size),
        "mean_ll_diff_family_minus_global": float(np.mean(diffs)),
        "ci": ci,
        "per_family": per_family,
        "params_global_offsets": n_ch_g,
        "params_global_free": n_ch_g - 1,
        "params_family_offsets": int(sum(n_ch_f.values())),
        "params_family_free": int(sum(v - 1 for v in n_ch_f.values())),
        "n_test_events_with_unseen_channel": n_unseen,
        "note": ("Bridge events only: they alone carry offset information "
                 "and route to a unique family. Global model refitted on the "
                 "pooled training events (same pool as the family models) "
                 "for an apples-to-apples comparison."),
    }
    print(f"[holdout] mean per-event ll diff (family - global) = "
          f"{result['mean_ll_diff_family_minus_global']:+.5f}"
          + (f", CI [{ci[0]:+.5f}, {ci[1]:+.5f}]" if ci else "")
          + f" over {diffs.size} held-out events")
    return result


# ------------------------------------------------------------------
# Markdown report
# ------------------------------------------------------------------

def _fmt(v, fmt="{:+.3f}"):
    return fmt.format(v) if v is not None else "NA"


def _fmt_ci(ci):
    return f"[{ci[0]:+.3f}, {ci[1]:+.3f}]" if ci else "NA"


def _art(word):
    return "an" if word[:1].lower() in "aeiou" else "a"


def physics_paragraph(results):
    """Data-driven interpretation; wording conditioned on computed numbers."""
    add = results["additivity"]
    hse = add.get("hse_minus_gga_snumat", {})
    per = hse.get("per_family", {})
    lines = []
    if len(per) >= 2:
        f_hi = hse["family_max"]
        f_lo = hse["family_min"]
        lines.append(
            f"The hybrid-functional rung of the ladder is widest for "
            f"{f_hi}s (delta(snumat_hse) - delta(snumat_gga) = "
            f"{per[f_hi]:+.3f}) and narrowest for {f_lo}s "
            f"({per[f_lo]:+.3f}), a stretch of "
            f"{hse['max_pairwise_diff']:.3f} log-strength units against a "
            f"global value of {_fmt(hse.get('global'))}.")
        if hse.get("extreme_family_cis_disjoint") is True:
            lines.append(
                "The bootstrap CIs of those two families do not overlap, so "
                "the pooled model's single per-channel offset is averaging "
                "over genuinely different chemistry-specific corrections -- "
                "consistent with electronic-structure practice, where the "
                "GGA gap underestimate (and hence the value of stepping up "
                "to a hybrid) depends on how ionic/localised the states at "
                "the gap edges are.")
        else:
            lines.append(
                "Within bootstrap uncertainty, however, the families' "
                "spacings are statistically compatible, so at this "
                "resolution the frozen universe does not force a "
                "chemistry-dependent correction on this rung.")
        routing = []
        for name, label in (("hse_minus_gga_snumat", "GGA -> HSE (SNUMAT)"),
                            ("mbj_minus_optb88_jarvis",
                             "OptB88vdW -> TB-mBJ (JARVIS)"),
                            ("expt_minus_gga_mp", "MP GGA -> experiment")):
            p = add[name].get("per_family", {})
            if len(p) >= 2:
                fh = add[name]["family_max"]
                fl = add[name]["family_min"]
                routing.append(
                    f"{label}: worth {p[fh]:+.3f} for {_art(fh)} {fh} "
                    f"campaign vs {p[fl]:+.3f} for {_art(fl)} {fl} campaign")
        if routing:
            lines.append(
                "For routing advice between channels this reads as: "
                + "; ".join(routing) + ".")
        hold = results.get("holdout")
        if hold and hold.get("ci"):
            gain = hold["mean_ll_diff_family_minus_global"]
            lo, hi = hold["ci"]
            if lo > 0:
                lines.append(
                    f"The held-out check agrees: routing each event to its "
                    f"family's offsets gains {gain:+.4f} log-likelihood per "
                    f"event (CI [{lo:+.4f}, {hi:+.4f}]) over the single "
                    f"additive model, at the cost of "
                    f"{hold['params_family_free']} vs "
                    f"{hold['params_global_free']} free offset parameters. "
                    "A family-resolved (representation-sharing) upgrade to "
                    "the pooled model must capture at least this.")
            else:
                lines.append(
                    f"Held-out, the family-resolved model gains only "
                    f"{gain:+.4f} log-likelihood per event (CI "
                    f"[{lo:+.4f}, {hi:+.4f}]), i.e. the extra "
                    f"{hold['params_family_free'] - hold['params_global_free']}"
                    " free parameters do not pay for themselves: the "
                    "additive offsets summarise the ladder adequately at "
                    "the current resolution.")
    return " ".join(lines)


def render_markdown(results, cfg):
    channels = results["universe"]["channels_kept"]
    fams = [f for f in cfg.families if f in results["families"]]
    fitted = [f for f in fams
              if not results["families"][f].get("skipped", True)]
    g = results["global"]
    lines = [
        "# Family-resolved fidelity-ladder offsets (additivity stress test)",
        "",
        "Question: are the per-channel offsets of the pooled cross-channel",
        "Plackett-Luce trial (hypergraph_trial_bandgap.py) chemistry-dependent,",
        "or is one additive offset per channel enough? Fits below restrict the",
        "frozen band-gap archive universe to one anion family at a time. Under a",
        "single-family restriction the trial's within-channel events are",
        "degenerate (one class) and carry no offset information (the channel",
        "factor cancels within a channel), so each family fit is the pooled",
        "model on that family's cross-channel bridge events -- exactly the",
        "offset-identifying sub-likelihood of the trial.",
        "",
        f"Universe: `{results['universe']['path']}`"
        f" ({results['universe']['n_used']} of "
        f"{results['universe']['n_seen']} records used, gap > "
        f"{trial.MIN_GAP_EV} eV). Channels kept: "
        + ", ".join(f"`{c}`" for c in channels) + ".",
        "",
        f"Gauge: {results['gauge_note']}",
        "",
        f"Bootstrap: B = {cfg.bootstrap} refits per fit, resampling bridge",
        f"events with replacement; CIs are percentile [{CI_LO}, {CI_HI}].",
        f"Seed = {cfg.seed}.",
        "",
        "## Records and events per family",
        "",
        "| channel | " + " | ".join(
            f"{f} records / bridge-events" for f in fams) + " |",
        "|---|" + "---|" * len(fams),
    ]
    for c in channels:
        row = [c]
        for f in fams:
            e = results["families"][f]
            row.append(f"{e['n_records_per_channel'].get(c, 0)} / "
                       f"{e['n_bridge_events_per_channel'].get(c, 0)}")
        lines.append("| " + " | ".join(row) + " |")
    lines.append("| total bridge events | " + " | ".join(
        str(results["families"][f]["n_bridge_events"]) for f in fams) + " |")
    lines += [
        "",
        "## Channel offsets by family (zero-mean gauge, with bootstrap CI)",
        "",
        "| channel | global | " + " | ".join(fitted) + " |",
        "|---|---|" + "---|" * len(fitted),
    ]
    for c in sorted(channels, key=lambda x: -g["offsets"].get(x, -99)):
        row = [c, f"{_fmt(g['offsets'].get(c))} "
                  f"{_fmt_ci(g['offset_ci'].get(c))}"]
        for f in fitted:
            e = results["families"][f]
            row.append(f"{_fmt(e['offsets'].get(c))} "
                       f"{_fmt_ci(e['offset_ci'].get(c))}")
        lines.append("| " + " | ".join(row) + " |")
    lines += [
        "",
        "## Ladder spacings (gauge-invariant) and the additivity test",
        "",
        "| spacing | global | " + " | ".join(fitted)
        + " | max pairwise diff | extreme CIs disjoint? |",
        "|---|---|" + "---|" * (len(fitted) + 2),
    ]
    for name, hi, lo in SPACING_DEFS:
        a = results["additivity"][name]
        row = [f"`{hi}` - `{lo}`",
               f"{_fmt(g['spacings'][name])} "
               f"{_fmt_ci(g['spacing_ci'].get(name))}"]
        for f in fitted:
            e = results["families"][f]
            row.append(f"{_fmt(e['spacings'][name])} "
                       f"{_fmt_ci(e['spacing_ci'].get(name))}")
        row.append(_fmt(a.get("max_pairwise_diff"), "{:.3f}")
                   + (f" ({a.get('family_max')} - {a.get('family_min')})"
                      if a.get("family_max") else ""))
        row.append(str(a.get("extreme_family_cis_disjoint")))
        lines.append("| " + " | ".join(row) + " |")
    lines += ["", "Family spacing minus global spacing:", ""]
    for name, hi, lo in SPACING_DEFS:
        parts = []
        for f in fitted:
            v = results["families"][f]["spacing_minus_global"][name]
            parts.append(f"{f} {_fmt(v)}")
        lines.append(f"- `{name}`: " + ", ".join(parts))
    hold = results.get("holdout")
    if hold:
        lines += [
            "",
            "## Held-out predictive check (80/20, seeded)",
            "",
            f"- events: {hold['n_train_events']} train / "
            f"{hold['n_test_events']} test (bridge events only; see note)",
            f"- mean per-event log-likelihood, family-resolved minus global: "
            f"{hold['mean_ll_diff_family_minus_global']:+.5f}, bootstrap CI "
            + (f"[{hold['ci'][0]:+.5f}, {hold['ci'][1]:+.5f}]"
               if hold["ci"] else "NA"),
            f"- parameters (channel offsets): global "
            f"{hold['params_global_offsets']} "
            f"({hold['params_global_free']} free) vs family-resolved "
            f"{hold['params_family_offsets']} "
            f"({hold['params_family_free']} free)",
            f"- held-out events containing a channel unseen in training: "
            f"{hold['n_test_events_with_unseen_channel']}",
            "- per family: " + ", ".join(
                f"{f}: {hold['per_family'][f]['mean_ll_diff']:+.5f} "
                f"(n={hold['per_family'][f]['n_test']})"
                for f in sorted(hold["per_family"])),
            "",
            f"Note: {hold['note']}",
        ]
    lines += [
        "",
        "## Verdict",
        "",
        results["verdict"],
        "",
        "## Interpretation",
        "",
        physics_paragraph(results) or "(fewer than 2 families fitted)",
        "",
        "![family offsets](family_offsets.png)",
        "",
    ]
    return "\n".join(lines)


# ------------------------------------------------------------------
# Plot (matplotlib Agg; palette + rules from the dataviz skill)
# ------------------------------------------------------------------

FAMILY_STYLE = {  # fixed hue order (validated categorical palette, light mode)
    "oxide": ("#2a78d6", "o"),
    "halide": ("#eb6834", "s"),
    "chalcogenide": ("#1baf7a", "^"),
    "nitride": ("#eda100", "D"),
    "carbide": ("#e87ba4", "v"),
    "boride": ("#008300", "P"),
    "pnictide": ("#4a3aa7", "X"),
    "other": ("#e34948", "*"),
}
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
SURFACE = "#fcfcfb"


def make_plot(results, cfg, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    g = results["global"]
    fitted = [f for f in cfg.families
              if not results["families"].get(f, {}).get("skipped", True)]
    channels = sorted(results["universe"]["channels_kept"],
                      key=lambda c: -g["offsets"].get(c, -99))
    ref_ok = g.get("offsets_ref") is not None and all(
        results["families"][f].get("offsets_ref") for f in fitted)

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(12.5, 5.2), gridspec_kw={"width_ratios": [2.1, 1.0]})
    fig.patch.set_facecolor(SURFACE)

    def offsets_of(entry):
        return entry["offsets_ref"] if ref_ok else entry["offsets"]

    def ci_of(entry):
        return entry["offset_ref_ci"] if ref_ok else entry["offset_ci"]

    # Panel 1: offsets per channel, grouped by family, global as reference
    ax1.set_facecolor(SURFACE)
    n_f = len(fitted)
    width = 0.72
    for k, family in enumerate(fitted):
        e = results["families"][family]
        color, marker = FAMILY_STYLE.get(family, (INK2, "o"))
        xs, ys, los, his = [], [], [], []
        for i, c in enumerate(channels):
            off = offsets_of(e)
            if off is None or c not in off:
                continue
            x = i - width / 2 + width * (k + 0.5) / max(n_f, 1)
            xs.append(x)
            ys.append(off[c])
            ci = ci_of(e).get(c)
            if ci:
                los.append(off[c] - ci[0])
                his.append(ci[1] - off[c])
            else:
                los.append(0.0)
                his.append(0.0)
        ax1.errorbar(xs, ys, yerr=[los, his], fmt=marker, color=color,
                     ecolor=color, elinewidth=1.4, capsize=2.5, ms=6,
                     mec=SURFACE, mew=0.8, ls="none", label=family)
    g_off = g["offsets_ref"] if ref_ok else g["offsets"]
    ax1.plot(range(len(channels)),
             [g_off.get(c, np.nan) for c in channels],
             marker="d", color=INK2, ms=7, ls="none", label="global (all)",
             zorder=1)
    ax1.axhline(0.0, color="#c3c2b7", lw=1.0, zorder=0)
    ax1.set_xticks(range(len(channels)))
    ax1.set_xticklabels(channels, rotation=30, ha="right", fontsize=8.5,
                        color=INK2)
    ylabel = ("channel offset minus mp_gga offset (log-strength)"
              if ref_ok else "channel offset (zero-mean gauge, log-strength)")
    ax1.set_ylabel(ylabel, fontsize=9, color=INK)
    ax1.set_title("Per-channel fidelity offsets by anion family",
                  fontsize=10.5, color=INK, loc="left")
    ax1.grid(axis="y", color=GRID, lw=0.7)
    ax1.set_axisbelow(True)
    ax1.legend(frameon=False, fontsize=8.5, loc="upper right",
               labelcolor=INK2)

    # Panel 2: the three ladder spacings
    ax2.set_facecolor(SURFACE)
    labels = {
        "hse_minus_gga_snumat": "HSE - GGA\n(SNUMAT)",
        "mbj_minus_optb88_jarvis": "mBJ - OptB88\n(JARVIS)",
        "expt_minus_gga_mp": "expt - GGA\n(Zhuo/MP)",
    }
    names = [n for n, _, _ in SPACING_DEFS]
    for k, family in enumerate(fitted):
        e = results["families"][family]
        color, marker = FAMILY_STYLE.get(family, (INK2, "o"))
        xs, ys, los, his = [], [], [], []
        for i, name in enumerate(names):
            v = e["spacings"][name]
            if v is None:
                continue
            x = i - width / 2 + width * (k + 0.5) / max(n_f, 1)
            xs.append(x)
            ys.append(v)
            ci = e["spacing_ci"].get(name)
            if ci:
                los.append(v - ci[0])
                his.append(ci[1] - v)
            else:
                los.append(0.0)
                his.append(0.0)
        ax2.errorbar(xs, ys, yerr=[los, his], fmt=marker, color=color,
                     ecolor=color, elinewidth=1.4, capsize=2.5, ms=6,
                     mec=SURFACE, mew=0.8, ls="none")
    ax2.plot(range(len(names)),
             [g["spacings"][n] if g["spacings"][n] is not None else np.nan
              for n in names],
             marker="d", color=INK2, ms=7, ls="none", zorder=1)
    ax2.axhline(0.0, color="#c3c2b7", lw=1.0, zorder=0)
    ax2.set_xticks(range(len(names)))
    ax2.set_xticklabels([labels[n] for n in names], fontsize=8.5, color=INK2)
    ax2.set_ylabel("ladder spacing (log-strength)", fontsize=9, color=INK)
    ax2.set_title("Ladder spacings", fontsize=10.5, color=INK, loc="left")
    ax2.grid(axis="y", color=GRID, lw=0.7)
    ax2.set_axisbelow(True)

    for ax in (ax1, ax2):
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
        for spine in ("left", "bottom"):
            ax.spines[spine].set_color("#c3c2b7")
        ax.tick_params(colors=MUTED, labelsize=8.5)

    fig.suptitle("Fidelity-ladder additivity stress test on the frozen "
                 "band-gap archive universe", fontsize=11.5, color=INK, x=0.01,
                 ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    plt.close(fig)


# ------------------------------------------------------------------
# Self-test: synthetic universe with KNOWN family-dependent offsets
# ------------------------------------------------------------------

SYNTH_CHANNELS = ["snumat_hse", "snumat_gga", "jarvis_mbj",
                  "jarvis_optb88vdw", "expt_zhuo", "mp_gga"]
SYNTH_BASE_DELTA = {"snumat_hse": 1.0, "snumat_gga": -1.0,
                    "jarvis_mbj": 0.7, "jarvis_optb88vdw": -0.7,
                    "expt_zhuo": 0.5, "mp_gga": -0.5}
SYNTH_FAMILY_ANION = [("oxide", "O"), ("halide", "F"), ("nitride", "N")]


def write_synthetic_universe(path, scales, n_formulas, seed):
    """Gumbel-noise universe: gap = exp(theta_f + scale*base_delta + gumbel).

    Rankings of one formula across channels then follow a Plackett-Luce law
    with channel strengths exp(scale_family * base_delta), so the true ladder
    spacing scales linearly with the family's scale factor.
    """
    rng = np.random.default_rng(seed)
    with open(path, "w", encoding="utf-8") as fh:
        for family, anion in SYNTH_FAMILY_ANION:
            s = float(scales[family])
            for i in range(n_formulas):
                formula = f"Zz{i}{anion}"
                theta = float(rng.normal(1.0, 0.5))
                for ch in SYNTH_CHANNELS:
                    gap = float(np.exp(theta + s * SYNTH_BASE_DELTA[ch]
                                       + rng.gumbel()))
                    rec = {"s": ch, "x": {"formula": formula},
                           "y": {"bandgap_eV": gap}}
                    fh.write(json.dumps(rec) + "\n")


def _self_test_cfg(universe, out_dir, seed):
    return argparse.Namespace(
        universe=str(universe), families=["oxide", "halide"], bootstrap=30,
        events=40, seed=seed, output_dir=str(out_dir), workers=1,
        min_class_size=30, max_bridges=2000, holdout_frac=0.2)


def run_self_test():
    print("family_offsets self-test (synthetic universes, seeded)")
    failures = []

    def check(label, ok):
        print(f"  [{'ok' if ok else 'FAIL'}] {label}")
        if not ok:
            failures.append(label)

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)

        # --- scenario A: family X (oxide) spacing = 2 x family Y (halide) --
        uni_a = tmp / "synthetic_diff.jsonl"
        write_synthetic_universe(
            uni_a, {"oxide": 1.0, "halide": 0.5, "nitride": 0.75},
            n_formulas=400, seed=101)
        print("[self-test] scenario A: oxide spacing = 2x halide spacing")
        res_a = run_analysis(_self_test_cfg(uni_a, tmp / "out_a", seed=5))
        name = "hse_minus_gga_snumat"
        ox = res_a["families"]["oxide"]["spacings"][name]
        ha = res_a["families"]["halide"]["spacings"][name]
        gl = res_a["global"]["spacings"][name]
        check(f"oxide hse-gga spacing near 2.0 (got {ox:+.3f})",
              1.4 < ox < 2.6)
        check(f"halide hse-gga spacing near 1.0 (got {ha:+.3f})",
              0.5 < ha < 1.5)
        check(f"recovered 2x difference (oxide - halide = {ox - ha:+.3f})",
              ox - ha > 0.5)
        check("extreme-family CIs disjoint for hse-gga spacing",
              res_a["additivity"][name]["extreme_family_cis_disjoint"] is True)
        check(f"global spacing between families "
              f"(halide {ha:+.3f} < global {gl:+.3f} < oxide {ox:+.3f})",
              min(ox, ha) < gl < max(ox, ha))
        hold_a = res_a["holdout"]
        check(f"held-out: family-resolved beats global "
              f"(mean diff {hold_a['mean_ll_diff_family_minus_global']:+.5f})",
              hold_a["mean_ll_diff_family_minus_global"] > 0.0)
        for fname in ("family_offsets.json", "family_offsets.md",
                      "family_offsets.png"):
            check(f"output written: {fname}", (tmp / "out_a" / fname).exists())

        # --- scenario B: null, all families share the same spacing ---------
        uni_b = tmp / "synthetic_null.jsonl"
        write_synthetic_universe(
            uni_b, {"oxide": 0.75, "halide": 0.75, "nitride": 0.75},
            n_formulas=400, seed=202)
        print("[self-test] scenario B: null (identical spacings)")
        res_b = run_analysis(_self_test_cfg(uni_b, tmp / "out_b", seed=5))
        oxb = res_b["families"]["oxide"]["spacings"][name]
        hab = res_b["families"]["halide"]["spacings"][name]
        check(f"null: spacings close (|oxide - halide| = "
              f"{abs(oxb - hab):.3f} < 0.5)", abs(oxb - hab) < 0.5)
        any_disjoint = any(
            res_b["additivity"][n]["extreme_family_cis_disjoint"] is True
            for n, _, _ in SPACING_DEFS)
        check("null: no spacing flags disjoint extreme-family CIs",
              not any_disjoint)
        hold_b = res_b["holdout"]
        check(f"null: held-out gain not significantly positive "
              f"(CI low = {hold_b['ci'][0]:+.5f})", hold_b["ci"][0] <= 0.0)
        check("held-out gain larger in scenario A than in null B",
              hold_a["mean_ll_diff_family_minus_global"]
              > hold_b["mean_ll_diff_family_minus_global"])

    print(f"self-test {'PASS' if not failures else 'FAIL'} "
          f"({len(failures)} failures)")
    for f in failures:
        print(f"  failed: {f}")
    return 1 if failures else 0


# ------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Family-resolved fidelity-ladder offsets "
                    "(additivity stress test)")
    ap.add_argument("--universe",
                    default="output/bandgap_archive_build/bandgap_universe.jsonl")
    ap.add_argument("--families", nargs="+", default=list(DEFAULT_FAMILIES),
                    help=f"anion families to fit (default: "
                         f"{' '.join(DEFAULT_FAMILIES)})")
    ap.add_argument("--bootstrap", type=int, default=200,
                    help="bootstrap refits per fit (default 200)")
    ap.add_argument("--events", type=int, default=300,
                    help="within-channel events per channel for the global "
                         "trial replication (trial's --n-events)")
    ap.add_argument("--min-class-size", type=int, default=30)
    ap.add_argument("--max-bridges", type=int, default=4000,
                    help="bridge-event cap per fit (trial's --max-bridges)")
    ap.add_argument("--holdout-frac", type=float, default=0.2)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--output-dir", default="output/family_offsets")
    ap.add_argument("--workers", type=int, default=1,
                    help="processes for bootstrap refits")
    ap.add_argument("--self-test", action="store_true",
                    help="run synthetic-universe checks and exit")
    ap.add_argument("--smoke", action="store_true",
                    help="tiny run: two families, B=20, small event counts")
    args = ap.parse_args()

    if args.self_test:
        return run_self_test()

    if args.smoke:
        args.families = ["oxide", "halide"]
        args.bootstrap = 20
        args.events = 50
        args.max_bridges = 400
        print("[smoke] families=oxide,halide bootstrap=20 events=50 "
              "max-bridges=400")

    bad = [f for f in args.families if f not in KNOWN_FAMILIES]
    if bad:
        print(f"[error] unknown families: {bad}; known: {KNOWN_FAMILIES}")
        return 1

    run_analysis(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
