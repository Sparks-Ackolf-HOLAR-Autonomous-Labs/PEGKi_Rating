"""
decision_body_cycles.py
-----------------------
Graded cyclic structure among DECISION BODIES, with material classes as
the comparison case.

Design considerations carried over from the earlier triangle analysis
(local_cycles.py), three points:
  1. Cyclic behaviour is often loose: Blade-Chest vectors that overlap or
     are sparse, loops that are not one regular rock-paper-scissors
     shape. "No cycle anywhere" is a much stronger claim than "some loose
     cyclic behaviour", and the need for exploration versus exploitation
     itself presupposes that no single choice wins everywhere.
  2. A yes/no "is this triangle a cycle" test is too binary.
  3. Cycles are expected among the decision bodies (learners,
     acquisition rules, reuse strategies), each with its own topology,
     not only among material classes, which were one example.

So this script (a) replaces the binary test with graded quantities and
(b) runs them on every committed decision-body tournament on materials
data, with the material-class tournaments recomputed the same way.

Graded quantities (all computed on the observed comparison graph)
-----------------------------------------------------------------
  loop share        Hodge split of the pairwise log-odds flow into a
                    consistent ranking (gradient) and loops (local curl
                    plus global harmonic). Loop share = curl + harmonic
                    energy fraction, in [0, 1]. Bootstrap interval.
  expected cycles   For each fully observed triangle, the cyclic
                    strength c = max over the two orientations of
                    min(P_ab, P_bc, P_ca) - 1/2 (c > 0: a majority cycle;
                    c slightly below 0: a loose, nearly closed loop).
                    The expected number of cyclic triangles is the sum
                    over triangles of the bootstrap probability that
                    c > 0. No threshold on any single triangle.
  loop regularity   The loop part of the flow is a skew-symmetric matrix;
                    its singular values come in pairs, one pair per
                    Blade-Chest plane. Regularity = energy in the top
                    plane / total loop energy. Near 1: one regular
                    rock-paper-scissors cycle; low: loops spread over
                    many overlapping planes (the "loose" case).
  context reversal  For bodies i, j and contexts (environments or
                    campaign stages), reversal strength
                    r_ij = min(max_ctx (P_ij - 1/2), max_ctx (P_ji - 1/2)).
                    r > 0: i wins in one context and loses in another;
                    its size is the smaller of the two margins.
  commitment regret For contexts A and B: commit to A's best body and
                    use it in B; regret = best score in B minus the
                    committed body's score in B. The value of continuing
                    to explore which decision body to use.

Strata, and what each one means
-------------------------------
  Replicates (seeds, bootstrap resamples) are NOT strata. They measure
  noise and only feed the intervals.
  Environment: the data source (receiver), the property, or the colour
    engine. Does a decision body's standing survive a change of the
    knowledge space it acts in? This is the transfer question itself.
  Campaign stage: early budget (checkpoints within the pre-registered
    early window) versus late budget. Does the body that wins while
    little is known also win once much is known? This is the
    exploration-versus-exploitation question.
  Query format: pairwise versus k-way events, kept as separate
    tournaments because the format is a physical batch choice.
  Material classes: the contexts are the eight data sources (the same
    class pair judged under different gap definitions). Crystal-system
    and element-count splits (local_cycles.py) are scope restrictions a
    campaign can adopt ("search only ternaries"); they are not rerun
    here. Splitting by gap window is rejected: conditioning on the
    outcome being ranked truncates the distributions and manufactures
    near-ties.

Decision-body tournaments (all from committed artifacts)
-------------------------------------------------------
  ML models   JARVIS-Leaderboard test MAEs: machine-learning models are
              the bodies, each property is one k-way round (lower MAE
              wins). Contexts: property family (electronic, energetic,
              mechanical). Eight properties are too few for a
              bootstrap, so intervals leave one property out at a time.
              DFT reference entries (es_dft_*) are excluded; they are
              not learners.
  Reuse       Exp 26 paired campaigns (CHPC, 50 matched seeds): bodies
  strategies  are campaign strategies (no reuse; pooled reuse; fixed-
              trust reuse of donor d; learned-trust reuse of donor d);
              one round = one (receiver, seed, budget checkpoint), the
              score is Kendall tau to the receiver's reference.
              Separate tournaments for pairwise queries and for 4-way
              and 6-way events. Contexts: receiver, stage. Clusters:
              (receiver, seed).
  Acquisition Exp 27: pairwise random, 4-way random, minimum-spread,
  rules       Fisher-greedy, trained policy. One round = one (held-out
              source, checkpoint); scores are seed means (per-seed
              curves were not stored), so intervals leave one held-out
              source out at a time. Contexts: held-out source, stage.
  More rounds (learning campaigns run for far more than eight rounds):
  ML models,  every JARVIS-Leaderboard AI single-property task scored by
  all tasks   MAE with at least three models (43 tasks; scraped 28 Sept
              2026 into data/jarvis/leaderboard_scores_ai_all.json).
              Cluster bootstrap over tasks. Contexts: property family.
  Acquisition Exp 27b: per-seed curves, 100 seeds per held-out source,
  rules, long budget to 1920 evaluations. One round = one (source, seed,
  horizon     checkpoint); cluster bootstrap over (source, seed). Stage:
              early = up to 240 evaluations, late = beyond.
  Colour      The nine acquisition algorithms of the colour testbed
  testbed     (regression layer only; not materials evidence). One round
              = one (engine, experiment); each experiment stores only its
              final round, so the only context is the engine.

ASCII only. Usage (from shared/):
    uv run python scripts/decision_body_cycles.py [--universe PATH]
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import argparse
import hashlib
import itertools
import json

import numpy as np

import importlib.util


def _load(name):
    path = Path(__file__).resolve().parent / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


N_BOOT = 200
MIN_BOOT_CLUSTERS = 20
EARLY_MAX = 120          # exp 26 pre-registered early window ends at 120
LB_FAMILY = {
    "dft_3d_optb88vdw_bandgap": "electronic", "dft_3d_mbj_bandgap": "electronic",
    "dft_3d_slme": "electronic",
    "dft_3d_formation_energy_peratom": "energetic", "dft_3d_ehull": "energetic",
    "dft_3d_exfoliation_energy": "energetic",
    "dft_3d_bulk_modulus_kv": "mechanical", "dft_3d_shear_modulus_gv": "mechanical",
}
LB_MIN_PROPS = 4
LB_ALL_MIN_MODELS = 3


# ------------------------------------------------------------ round tensors
class Rounds:
    """Scores per round (higher is better; NaN = body absent), with the
    cluster id used by the bootstrap and context labels per round."""

    def __init__(self, bodies, S, cluster, contexts):
        self.bodies = list(bodies)
        self.S = np.asarray(S, float)
        self.cluster = np.asarray(cluster)
        self.contexts = {k: np.asarray(v) for k, v in contexts.items()}
        n = len(self.bodies)
        pres = ~np.isnan(self.S)
        A = np.where(pres, self.S, 0.0)
        both = pres[:, :, None] & pres[:, None, :]
        gt = (A[:, :, None] > A[:, None, :]) & both
        eq = (A[:, :, None] == A[:, None, :]) & both
        idx = np.arange(n)
        both[:, idx, idx] = False
        eq[:, idx, idx] = False
        self.T = (gt + 0.5 * eq).reshape(len(self.S), -1).astype(float)
        self.M = both.reshape(len(self.S), -1).astype(float)
        self.ucl, self.cinv = np.unique(self.cluster, return_inverse=True)

    def wins(self, w=None):
        n = len(self.bodies)
        if w is None:
            w = np.ones(len(self.S))
        return (w @ self.T).reshape(n, n), (w @ self.M).reshape(n, n)

    def resample_weights(self, rng, n_boot):
        """Cluster bootstrap when there are at least MIN_BOOT_CLUSTERS
        clusters; otherwise leave-one-cluster-out (jackknife), because a
        bootstrap over 8 properties or 8 held-out sources mostly
        duplicates clusters and says little about uncertainty."""
        k = len(self.ucl)
        if k >= MIN_BOOT_CLUSTERS:
            for _ in range(n_boot):
                draw = np.bincount(rng.integers(0, k, k), minlength=k)
                yield draw[self.cinv].astype(float)
        else:
            for c in range(k):
                yield (self.cinv != c).astype(float)

    def scheme(self):
        return "cluster bootstrap" if len(self.ucl) >= MIN_BOOT_CLUSTERS else "leave-one-cluster-out"


# ------------------------------------------------------------ graded measures
def hodge(W, N):
    """Hodge split of y_ij = log((W_ij + 1/2) / (W_ji + 1/2)) on the
    observed graph. Returns fractions and the loop residual matrix."""
    n = W.shape[0]
    edges = [(i, j) for i in range(n) for j in range(i + 1, n) if N[i, j] > 0]
    if not edges:
        return None
    eidx = {e: k for k, e in enumerate(edges)}
    y = np.array([np.log((W[i, j] + 0.5) / (W[j, i] + 0.5)) for i, j in edges])
    E = len(edges)
    B1 = np.zeros((E, n))
    for k, (i, j) in enumerate(edges):
        B1[k, i], B1[k, j] = 1.0, -1.0
    tris = [t for t in itertools.combinations(range(n), 3)
            if (t[0], t[1]) in eidx and (t[1], t[2]) in eidx and (t[0], t[2]) in eidx]
    tot = float(y @ y)
    g = B1 @ np.linalg.lstsq(B1, y, rcond=None)[0]
    r = y - g
    if tris:
        B2 = np.zeros((E, len(tris)))
        for c, (i, j, k) in enumerate(tris):
            B2[eidx[(i, j)], c] = 1.0
            B2[eidx[(j, k)], c] = 1.0
            B2[eidx[(i, k)], c] = -1.0
        cp = B2 @ np.linalg.lstsq(B2, r, rcond=None)[0]
    else:
        cp = np.zeros(E)
    h = r - cp
    R = np.zeros((n, n))
    for k, (i, j) in enumerate(edges):
        R[i, j], R[j, i] = r[k], -r[k]
    sv = np.linalg.svd(R, compute_uv=False)
    le = float((sv ** 2).sum())
    reg = float((sv[:2] ** 2).sum() / le) if le > 1e-12 else float("nan")
    if tot <= 0:
        return {"grad": 1.0, "curl": 0.0, "harm": 0.0, "loop": 0.0,
                "regularity": float("nan"), "n_edges": E, "n_tri": len(tris)}
    return {"grad": float(g @ g) / tot, "curl": float(cp @ cp) / tot,
            "harm": float(h @ h) / tot, "loop": float(r @ r) / tot,
            "regularity": reg, "n_edges": E, "n_tri": len(tris)}


def regularity_null(N, rng, n_draw=200):
    """Regularity expected when the loop part is structureless: i.i.d.
    Gaussian flows on the same observed graph, gradient removed."""
    n = N.shape[0]
    vals = []
    for _ in range(n_draw):
        Z = rng.normal(size=(n, n))
        Z = np.triu(Z, 1)
        Z = Z - Z.T
        h = hodge(np.exp(Z), N)   # y = log((e^z + 1/2) / (e^-z + 1/2)), odd in z
        if h and not np.isnan(h["regularity"]):
            vals.append(h["regularity"])
    return float(np.mean(vals)) if vals else float("nan")


def top_plane(W, N, names, k=4):
    """Bodies with the largest weight in the top Blade-Chest plane of the
    loop residual (the plane holding the most loop energy)."""
    n = W.shape[0]
    Y = np.zeros((n, n))
    obs = N > 0
    Y[obs] = np.log((W[obs] + 0.5) / (W.T[obs] + 0.5))
    edges = [(i, j) for i in range(n) for j in range(i + 1, n) if obs[i, j]]
    B1 = np.zeros((len(edges), n))
    y = np.array([Y[i, j] for i, j in edges])
    for e, (i, j) in enumerate(edges):
        B1[e, i], B1[e, j] = 1.0, -1.0
    s = np.linalg.lstsq(B1, y, rcond=None)[0]
    R = np.where(obs, Y - (s[:, None] - s[None, :]), 0.0)
    U, sv, Vt = np.linalg.svd(R)
    if sv[0] < 1e-9:
        return []
    load = U[:, 0] ** 2 + U[:, 1] ** 2
    order = np.argsort(-load)[:k]
    return [{"body": names[i], "weight": float(load[i] / load.sum())} for i in order]


def triangle_strength(P, obs):
    """Cyclic strength c per fully observed triangle (see module doc)."""
    n = P.shape[0]
    out = {}
    for a, b, c in itertools.combinations(range(n), 3):
        if not (obs[a, b] and obs[b, c] and obs[a, c]):
            continue
        s1 = min(P[a, b], P[b, c], P[c, a])
        s2 = min(P[b, a], P[c, b], P[a, c])
        out[(a, b, c)] = max(s1, s2) - 0.5
    return out


def pmat(W, N):
    with np.errstate(invalid="ignore", divide="ignore"):
        P = np.where(N > 0, W / np.maximum(N, 1e-12), 0.5)
    return P, N > 0


def reversal(Pc, Oc, n):
    """Pc: list of per-context P matrices. Returns {(i, j): r_ij}."""
    out = {}
    for i, j in itertools.combinations(range(n), 2):
        m = [P[i, j] - 0.5 for P, O in zip(Pc, Oc) if O[i, j]]
        if len(m) < 2:
            continue
        out[(i, j)] = float(min(max(m), max(-x for x in m)))
    return out


def commitment_regret(ctx_scores):
    """ctx_scores: {ctx: {body: mean score}} -> mean and max regret of
    committing to one context's winner in another context."""
    ctxs = sorted(ctx_scores)
    regs = []
    for a, b in itertools.permutations(ctxs, 2):
        common = set(ctx_scores[a]) & set(ctx_scores[b])
        if len(common) < 2:
            continue
        win = max(common, key=lambda k: ctx_scores[a][k])
        best = max(ctx_scores[b][k] for k in common)
        regs.append(best - ctx_scores[b][win])
    if not regs:
        return None
    return {"mean": float(np.mean(regs)), "max": float(np.max(regs)),
            "share_nonzero": float(np.mean([x > 1e-12 for x in regs])),
            "n_pairs": len(regs)}


def analyse(rounds, rng, n_boot, regret_note):
    n = len(rounds.bodies)
    W, N = rounds.wins()
    P, O = pmat(W, N)
    base = hodge(W, N)
    tri = triangle_strength(P, O)
    loops, regs = [], []
    hit = {t: 0 for t in tri}
    nres = 0
    for w in rounds.resample_weights(rng, n_boot):
        nres += 1
        Wb, Nb = rounds.wins(w)
        hb = hodge(Wb, Nb)
        if hb:
            loops.append(hb["loop"])
            if not np.isnan(hb["regularity"]):
                regs.append(hb["regularity"])
        Pb, Ob = pmat(Wb, Nb)
        for t, c in triangle_strength(Pb, Ob).items():
            if t in hit and c > 0:
                hit[t] += 1
    freq = {t: h / nres for t, h in hit.items()}
    name = rounds.bodies
    top = sorted(tri, key=lambda t: -tri[t])[:5]
    res = {
        "n_bodies": n, "n_rounds": int(len(rounds.S)),
        "n_clusters": int(len(rounds.ucl)),
        "resampling": rounds.scheme(), "n_resamples": nres,
        "hodge": base,
        "top_plane_bodies": top_plane(W, N, name),
        "regularity_null": regularity_null(N, rng),
        "loop_share_ci90": [float(np.percentile(loops, 5)), float(np.percentile(loops, 95))] if loops else None,
        "regularity_ci90": [float(np.percentile(regs, 5)), float(np.percentile(regs, 95))] if regs else None,
        "n_triangles": len(tri),
        "point_cycles": int(sum(c > 0 for c in tri.values())),
        "expected_cycles": float(sum(freq.values())),
        "triangles_undecided": int(sum(0.05 < f < 0.95 for f in freq.values())),
        "strength_quantiles": ([float(q) for q in np.percentile(list(tri.values()), [5, 25, 50, 75, 95])]
                               if tri else None),
        "strengths": [float(c) for c in tri.values()],
        "strongest_triangles": [
            {"triangle": " > ".join(name[i] for i in t), "strength": float(tri[t]),
             "p_cycle": float(freq[t])} for t in top],
        "contexts": {},
    }
    for cname, lab in rounds.contexts.items():
        levels = sorted(set(lab.tolist()))
        Pc, Oc, scores, within = [], [], {}, {}
        for lv in levels:
            w = (lab == lv).astype(float)
            Wc, Nc = rounds.wins(w)
            Pl, Ol = pmat(Wc, Nc)
            hc = hodge(Wc, Nc)
            within[str(lv)] = hc["loop"] if hc else float("nan")
            Pc.append(Pl)
            Oc.append(Ol)
            sub = rounds.S[lab == lv]
            sc = {}
            for k, b in enumerate(name):
                col = sub[:, k]
                col = col[~np.isnan(col)]
                if len(col):
                    sc[b] = float(col.mean())
            scores[str(lv)] = sc
        rv = reversal(Pc, Oc, n)
        vals = np.array(list(rv.values())) if rv else np.array([0.0])
        best = {lv: max(sc, key=sc.get) for lv, sc in scores.items() if sc}
        strongest = sorted(rv, key=lambda p: -rv[p])[:3]
        res["contexts"][cname] = {
            "levels": [str(x) for x in levels],
            "best_body": best,
            "n_distinct_best": len(set(best.values())),
            "within_context_loop_share": within,
            "reversal_pairs": int((vals > 0).sum()),
            "n_pairs": int(len(rv)),
            "reversal_mean_positive": float(vals[vals > 0].mean()) if (vals > 0).any() else 0.0,
            "reversal_max": float(vals.max()),
            "strongest_reversals": [
                {"pair": f"{name[i]} vs {name[j]}", "strength": rv[(i, j)]} for i, j in strongest
                if rv[(i, j)] > 0],
            "commitment_regret": commitment_regret(scores),
            "regret_unit": regret_note,
        }
    return res


# ------------------------------------------------------------ builders
def build_leaderboard(root):
    d = json.load(open(root / "data/jarvis/leaderboard_scores_combined.json"))
    count = {}
    for prop, sc in d.items():
        for m in sc:
            if not m.startswith("es_dft"):
                count[m] = count.get(m, 0) + 1
    bodies = sorted(m for m, c in count.items() if c >= LB_MIN_PROPS)
    S, cl, fam = [], [], []
    for prop in sorted(d):
        row = [(-np.log(d[prop][m]) if m in d[prop] else np.nan) for m in bodies]
        S.append(row)
        cl.append(prop)
        fam.append(LB_FAMILY.get(prop, "other"))
    return Rounds(bodies, S, cl, {"property family": fam})


def _arm_body(key):
    parts = key.split("|")
    if parts[0] in ("cold", "combined"):
        return parts[0], parts[1]
    return f"{parts[0]}:{parts[2]}", parts[1]


def build_campaign(ckdir, query, checkpoints):
    rows, bodies = [], set()
    for f in sorted(Path(ckdir).glob("*.json")):
        d = json.load(open(f))
        for ci, b in enumerate(checkpoints):
            sc = {}
            for key, curve in d["curves"].items():
                body, q = _arm_body(key)
                if q == query:
                    sc[body] = float(curve[ci])
            bodies |= set(sc)
            rows.append((d["receiver"], d["seed"], b, sc))
    bodies = sorted(bodies, key=lambda s: (s.split(":")[0] not in ("cold", "combined"), s))
    S = [[r[3].get(b, np.nan) for b in bodies] for r in rows]
    cl = [f"{r[0]}|{r[1]}" for r in rows]
    return Rounds(bodies, S, cl, {
        "receiver": [r[0] for r in rows],
        "stage": ["early" if r[2] <= EARLY_MAX else "late" for r in rows]})


def build_acquisition(root):
    base = root.parent / "campaign_protocols/results/learned_acquisition"
    agg = json.load(open(base / "learned_acquisition.json"))
    cps = agg["checkpoints"]
    half = cps[len(cps) // 2]
    bodies = ["pair_random", "kway_random", "kway_active", "kway_fisher", "kway_learned"]
    S, cl, src, st = [], [], [], []
    for f in sorted(base.glob("fold_*.json")):
        d = json.load(open(f))
        for ci, b in enumerate(cps):
            S.append([float(d["curves"][k][ci]) for k in bodies])
            cl.append(d["held_out"])
            src.append(d["held_out"])
            st.append("early" if b <= half else "late")
    return Rounds(bodies, S, cl, {"held-out source": src, "stage": st})


def lb_family(task):
    """Property family of a JARVIS-Leaderboard AI task, from its name."""
    t = task.lower()
    if t.startswith("hmof"):
        return "porous (MOF)"
    if any(k in t for k in ("encut", "kpoint")):
        return "calculation settings"
    if any(k in t for k in ("modulus", "_kv", "_gv")):
        return "mechanical"
    if any(k in t for k in ("gap", "mass", "slme", "eps", "spillage", "efg",
                            "powerfact", "piezo", "dielectric", "magmom",
                            "magnetic")):
        return "electronic"
    return "energetic"


def build_leaderboard_all(root):
    """All JARVIS-Leaderboard AI single-property tasks scored by MAE with
    at least LB_ALL_MIN_MODELS models (scraped 28 Sept 2026; see the
    _meta block of the data file). One task = one round."""
    d = json.load(open(root / "data/jarvis/leaderboard_scores_ai_all.json"))
    metric, scores = d["_meta"]["metric_column"], d["scores"]
    tasks = sorted(t for t in scores if metric.get(t) == "MAE"
                   and len(scores[t]) >= LB_ALL_MIN_MODELS)
    count = {}
    for t in tasks:
        for m in scores[t]:
            count[m] = count.get(m, 0) + 1
    bodies = sorted(m for m, c in count.items() if c >= LB_MIN_PROPS)
    S = [[(-np.log(scores[t][m]) if m in scores[t] and scores[t][m] > 0 else np.nan)
          for m in bodies] for t in tasks]
    return Rounds(bodies, S, tasks, {"property family": [lb_family(t) for t in tasks]})


def build_acquisition_long(root):
    """Exp 27b: per-seed curves, 100 seeds per held-out source, budget to
    1920. One round = one (held-out source, seed, checkpoint); cluster =
    (held-out source, seed). Stage: up to the exp 27 horizon of 240
    evaluations, 240-960, and 960-1920."""
    base = root.parent / "campaign_protocols/results/learned_acquisition_long"
    bodies = ["pair_random", "kway_random", "kway_active", "kway_fisher", "kway_learned"]
    S, cl, src, st = [], [], [], []
    for f in sorted(base.glob("fold_*.json")):
        d = json.load(open(f))
        cps = d["checkpoints"]
        for si, seed in enumerate(d["seeds"]):
            for ci, b in enumerate(cps):
                S.append([d["curves"][k][si][ci] for k in bodies])
                cl.append(f"{d['held_out']}|{seed}")
                src.append(d["held_out"])
                st.append("1: to 240" if b <= 240 else
                          ("2: 240-960" if b <= 960 else "3: 960-1920"))
    if not S:
        return None
    return Rounds(bodies, S, cl, {"held-out source": src, "stage": st})


def build_colour(root):
    try:
        ht = _load("hypergraph_trial")
    except Exception as exc:  # pragma: no cover
        print("[skip] colour testbed:", exc)
        return None
    per = {}
    for eng in ("db_spectral", "db_mixbox", "db_ryb", "db_km"):
        db = root / "output" / eng
        if not (db / "policies").exists():
            continue
        reader = ht.PolicyDatabaseReader(db)
        for pol in sorted(reader.list_policies()):
            for exp_id in sorted(reader.list_experiments(pol)):
                exp = reader.load_experiment(pol, exp_id)
                for rnd in exp.rounds:
                    per.setdefault((eng, exp_id, rnd.round_id), {})[pol] = -float(rnd.best_color_distance)
    if not per:
        return None
    bodies = sorted({p for sc in per.values() for p in sc})
    keys = sorted(per)
    S = [[per[k].get(b, np.nan) for b in bodies] for k in keys]
    # each colour experiment stores one final round, so there is no stage axis
    return Rounds(bodies, S, [f"{k[0]}|{k[1]}" for k in keys], {
        "engine": [k[0].replace("db_", "") for k in keys]})


def class_tournaments(universe, rng, n_boot):
    """Material classes judged per data source, same graded measures.
    Bootstrap resamples records within each class."""
    lc = _load("local_cycles")
    recs = lc.load(universe)
    out, perctx = {}, {}
    names_all = None
    for src in lc.LABEL:
        cl = {}
        for s, f, c, g, cs, es in recs:
            if s == src:
                cl.setdefault(c, []).append(g)
        cl = {c: np.asarray(v) for c, v in cl.items() if len(v) >= lc.MIN_CLASS}
        names, Pd = lc.pc.dominance(cl)
        n = len(names)
        P = np.full((n, n), 0.5)
        for (a, b), p in Pd.items():
            P[names.index(a), names.index(b)] = p
        W, N = P * 1000.0, np.ones((n, n)) * 1000.0 - np.eye(n) * 1000.0
        # log-odds on the exact dominance probabilities (large pseudo-count)
        base = hodge(W, N)
        tri = triangle_strength(P, N > 0)
        hit = {t: 0 for t in tri}
        loops = []
        for _ in range(n_boot):
            rs = {c: rng.choice(v, size=len(v), replace=True) for c, v in cl.items()}
            nb, Pb = lc.pc.dominance(rs)
            Pm = np.full((n, n), 0.5)
            for (a, b), p in Pb.items():
                Pm[nb.index(a), nb.index(b)] = p
            hb = hodge(Pm * 1000.0, N)
            loops.append(hb["loop"])
            for t, c in triangle_strength(Pm, N > 0).items():
                hit[t] += c > 0
        out[src] = {"n_bodies": n, "hodge": base, "regularity_null": regularity_null(N, rng),
                    "loop_share_ci90": [float(np.percentile(loops, 5)), float(np.percentile(loops, 95))],
                    "n_triangles": len(tri),
                    "point_cycles": int(sum(c > 0 for c in tri.values())),
                    "expected_cycles": float(sum(h / n_boot for h in hit.values())),
                    "strength_quantiles": [float(q) for q in np.percentile(list(tri.values()), [5, 25, 50, 75, 95])],
                    "strengths": [float(c) for c in tri.values()]}
        perctx[src] = (names, P)
        print("classes", src, round(base["loop"], 4), out[src]["expected_cycles"], flush=True)
    common = sorted(set.intersection(*[set(v[0]) for v in perctx.values()]))
    Pc, Oc = [], []
    for src, (names, P) in perctx.items():
        ix = [names.index(c) for c in common]
        Pc.append(P[np.ix_(ix, ix)])
        Oc.append(np.ones((len(ix), len(ix)), bool))
    rv = reversal(Pc, Oc, len(common))
    vals = np.array(list(rv.values()))
    top = sorted(rv, key=lambda p: -rv[p])[:3]
    out["_across_sources"] = {
        "classes": common, "reversal_pairs": int((vals > 0).sum()), "n_pairs": len(rv),
        "reversal_max": float(vals.max()),
        "strongest_reversals": [{"pair": f"{common[i]} vs {common[j]}", "strength": rv[(i, j)]}
                                for i, j in top if rv[(i, j)] > 0]}
    return out


# ------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--universe", default="output/bandgap_archive_build/bandgap_universe.jsonl")
    ap.add_argument("--output-dir", default="output/decision_body_cycles")
    ap.add_argument("--n-boot", type=int, default=N_BOOT)
    ap.add_argument("--skip-classes", action="store_true")
    args = ap.parse_args()
    root = Path(__file__).resolve().parent.parent
    out = (root / args.output_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(0)
    res = {"config": {"n_boot": args.n_boot, "early_max": EARLY_MAX,
                      "leaderboard_min_properties": LB_MIN_PROPS}, "decision_bodies": {}}

    tours = {}
    tours["ML models (JARVIS-Leaderboard)"] = (build_leaderboard(root),
                                               "log MAE ratio (0.10 = 10.5% more error)")
    chpc = root / "output/bandgap_cluster_run2"
    cps = json.load(open(chpc / "paired_campaign_k4/paired_campaign.json"))["checkpoints"]
    tours["Reuse strategies, pairwise queries"] = (
        build_campaign(chpc / "paired_campaign_k4/checkpoints", "pair", cps), "Kendall tau")
    tours["Reuse strategies, 4-way events"] = (
        build_campaign(chpc / "paired_campaign_k4/checkpoints", "kway", cps), "Kendall tau")
    tours["Reuse strategies, 6-way events"] = (
        build_campaign(chpc / "paired_campaign_k6/checkpoints", "kway", cps), "Kendall tau")
    tours["Acquisition rules (exp 27)"] = (build_acquisition(root), "Kendall tau")
    tours["ML models, all leaderboard tasks"] = (build_leaderboard_all(root),
                                                 "log MAE ratio (0.10 = 10.5% more error)")
    acq_long = build_acquisition_long(root)
    if acq_long is not None:
        tours["Acquisition rules, long horizon (exp 27b)"] = (acq_long, "Kendall tau")
    col = build_colour(root)
    if col is not None:
        tours["Acquisition algorithms (colour testbed)"] = (col, "colour distance")

    for name, (r, unit) in tours.items():
        res["decision_bodies"][name] = analyse(r, rng, args.n_boot, unit)
        h = res["decision_bodies"][name]
        print(name, h["n_bodies"], h["n_rounds"], round(h["hodge"]["loop"], 4),
              round(h["expected_cycles"], 2), flush=True)

    if not args.skip_classes:
        universe = Path(args.universe)
        if not universe.is_absolute():
            universe = (root / universe).resolve()
        res["universe_sha256_16"] = hashlib.sha256(universe.read_bytes()).hexdigest()[:16]
        res["material_classes"] = class_tournaments(universe, rng, args.n_boot)

    (out / "decision_body_cycles.json").write_text(json.dumps(res, indent=1, sort_keys=True))
    write_md(res, out / "decision_body_cycles.md")
    print("[done]", out)


def _f(x, d=3):
    return "n/a" if x is None or (isinstance(x, float) and np.isnan(x)) else f"{x:.{d}f}"


def write_md(res, path):
    L = ["# Graded cycles: decision bodies versus material classes", ""]
    if "universe_sha256_16" in res:
        L.append(f"Archive sha256[:16] = `{res['universe_sha256_16']}`; bootstrap {res['config']['n_boot']} resamples.")
        L.append("")
    L += ["## Whole tournaments", "",
          "| tournament | bodies | rounds | loop share (90% CI) | local curl | global harmonic | regularity (null) | triangles | point cycles | expected cycles | undecided triangles |",
          "|---|---|---|---|---|---|---|---|---|---|---|"]
    for name, h in res["decision_bodies"].items():
        ci = h["loop_share_ci90"]
        L.append(f"| {name} | {h['n_bodies']} | {h['n_rounds']} | {_f(h['hodge']['loop'])} "
                 f"[{_f(ci[0])}, {_f(ci[1])}] | {_f(h['hodge']['curl'])} | {_f(h['hodge']['harm'])} | "
                 f"{_f(h['hodge']['regularity'], 2)} ({_f(h['regularity_null'], 2)}) | {h['n_triangles']} | {h['point_cycles']} | "
                 f"{_f(h['expected_cycles'], 1)} | {h['triangles_undecided']} |")
    for src, h in res.get("material_classes", {}).items():
        if src.startswith("_"):
            continue
        ci = h["loop_share_ci90"]
        L.append(f"| Material classes, {src} | {h['n_bodies']} | records | {_f(h['hodge']['loop'])} "
                 f"[{_f(ci[0])}, {_f(ci[1])}] | {_f(h['hodge']['curl'])} | {_f(h['hodge']['harm'])} | "
                 f"{_f(h['hodge']['regularity'], 2)} ({_f(h['regularity_null'], 2)}) | {h['n_triangles']} | {h['point_cycles']} | "
                 f"{_f(h['expected_cycles'], 1)} | n/a |")
    L += ["", "Loop share = curl + harmonic energy fraction of the pairwise log-odds flow. "
          "Regularity = share of the loop energy in the top Blade-Chest plane (1 = one regular cycle); "
          "null = the same share for structureless random flows on the same comparison graph. "
          "Expected cycles = sum over triangles of the bootstrap probability that the triangle is cyclic. "
          "Undecided = triangles cyclic in between 5% and 95% of resamples.", ""]
    L += ["## Strongest near-cycles per tournament", ""]
    for name, h in res["decision_bodies"].items():
        L.append(f"**{name}**: " + "; ".join(
            f"{t['triangle']} (strength {t['strength']:+.3f}, P(cycle) {t['p_cycle']:.2f})"
            for t in h["strongest_triangles"][:3]))
        L.append("")
    L += ["## Contexts: reversal and commitment regret", "",
          "| tournament | context | levels | loop share within one level (mean, max) | distinct winners | reversed pairs | max reversal | mean regret | max regret | regret unit |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for name, h in res["decision_bodies"].items():
        for c, x in h["contexts"].items():
            cr = x["commitment_regret"] or {}
            wl = [v for v in x["within_context_loop_share"].values() if not np.isnan(v)]
            L.append(f"| {name} | {c} | {len(x['levels'])} | {_f(float(np.mean(wl)))}, {_f(float(np.max(wl)))} | {x['n_distinct_best']} | "
                     f"{x['reversal_pairs']} / {x['n_pairs']} | {_f(x['reversal_max'])} | "
                     f"{_f(cr.get('mean'))} | {_f(cr.get('max'))} | {x['regret_unit']} |")
    L.append("")
    for name, h in res["decision_bodies"].items():
        for c, x in h["contexts"].items():
            L.append(f"- {name}, by {c}: winners {x['best_body']}; strongest reversals: " + (
                "; ".join(f"{r['pair']} ({r['strength']:+.3f})" for r in x["strongest_reversals"]) or "none"))
    if "material_classes" in res:
        a = res["material_classes"]["_across_sources"]
        L += ["", "## Material classes across data sources", "",
              f"{a['reversal_pairs']} of {a['n_pairs']} class pairs reverse between at least two sources; "
              f"max reversal strength {a['reversal_max']:+.3f}. Strongest: " +
              ("; ".join(f"{r['pair']} ({r['strength']:+.3f})" for r in a["strongest_reversals"]) or "none")]
    path.write_text("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
