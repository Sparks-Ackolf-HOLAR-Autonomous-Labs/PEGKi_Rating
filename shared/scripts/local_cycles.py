"""
Is the absence of rock-paper-scissors cycles global only, or local as well?

pair_combinatorics.py showed that the class-level majority tournament of
every whole data source (and of every set-operation region) has zero directed
3-cycles. That is a statement about one aggregate table per source. This
script tests the same property at four local levels, from the version-locked
archive only:

  (L1) Sub-populations. The class tournament is rebuilt inside each crystal
       system (from the space-group number) and inside each element-count
       stratum (binary, ternary, quaternary and higher) of every source.
       Classes need at least MIN_STRATUM_CLASS records in the stratum.
  (L2) Robustness of every triangle. Records are bootstrap-resampled within
       each class; a triangle is a robust cycle if it is cyclic in at least
       half of the resamples, and fragile if it is cyclic in some but fewer.
       Applied to the whole sources and to every stratum of (L1).
  (L3) Every step of a campaign. A random k-way campaign (one evaluation =
       one archived record, k distinct classes per round, a round is the
       sorted order of its draws) is replayed; after each budget checkpoint
       the empirical win table of the observed pairs is decomposed by the
       Hodge theorem on the observed comparison graph into a gradient part,
       a local curl part (triangles), and a global harmonic part (loops of
       length four or more that no triangle fills). Directed 3-cycles among
       the fully observed triangles are counted. The population tournament
       has no cycle, so every cycle here is a sampling cycle; the question
       is how often they occur and when they vanish.
  (L4) Sources as players. For every pair of sources, P_ab is the fraction
       of shared compounds on which source a reports the wider gap. The
       tournament over the eight sources is tested for cycles globally and
       within each anion family (the local version of the fidelity ladder).

ASCII only. Usage (from shared/):
    uv run python scripts/local_cycles.py [--universe PATH] [--seeds 50]
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import argparse
import hashlib
import itertools
import json
import re

import numpy as np

import importlib.util


def _load(name):
    path = Path(__file__).resolve().parent / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


pc = _load("pair_combinatorics")
sos = pc.sos
trial = sos.trial
LABEL = sos.LABEL

MIN_CLASS = 30           # whole-source class size (set_operation_study default)
MIN_STRATUM_CLASS = 20   # class size inside a stratum
MIN_STRATUM_CLASSES = 3  # a tournament needs at least one triangle
N_BOOT = 200
BUDGETS = (8, 16, 24, 32, 48, 64, 96, 128, 176, 240, 480, 960, 2400)
CAMPAIGN_BUDGET = 240  # the paper's campaign budget; later checkpoints extend it

CRYSTAL = [(2, "triclinic"), (15, "monoclinic"), (74, "orthorhombic"),
           (142, "tetragonal"), (167, "trigonal"), (194, "hexagonal"),
           (230, "cubic")]


def crystal_system(spg):
    if spg is None:
        return None
    for hi, name in CRYSTAL:
        if int(spg) <= hi:
            return name
    return None


def n_elements(formula):
    return len(set(re.findall(r"[A-Z][a-z]?", formula)))


def element_stratum(formula):
    n = n_elements(formula)
    return {1: None, 2: "binary", 3: "ternary"}.get(n, "quaternary+")


def load(path):
    """-> list of (source, formula, class, gap, crystal system, element stratum)."""
    out = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            r = json.loads(line)
            x = r.get("x") or {}
            f = x.get("formula")
            g = (r.get("y") or {}).get("bandgap_eV")
            s = r.get("s")
            if not f or g is None or s is None or float(g) <= trial.MIN_GAP_EV:
                continue
            out.append((s, f, trial.classify_formula(f), float(g),
                        crystal_system(x.get("spg")), element_stratum(f)))
    return out


# ---------------------------------------------------------------- tournaments
def cyclic_triangles(names, P):
    """Triangles whose majority edges form a directed 3-cycle."""
    cyc = []
    for a, b, c in itertools.combinations(names, 3):
        ab, bc, ca = P[(a, b)] > 0.5, P[(b, c)] > 0.5, P[(c, a)] > 0.5
        if P[(a, b)] == 0.5 or P[(b, c)] == 0.5 or P[(c, a)] == 0.5:
            continue
        if ab == bc == ca:
            cyc.append((a, b, c))
    return cyc


def curl_stats(names, P):
    """Hodge split of Y = logit(P) on the complete graph: cyclic energy
    fraction and the largest triangle curl |Y_ab + Y_bc + Y_ca|."""
    n = len(names)
    idx = {c: i for i, c in enumerate(names)}
    Y = np.zeros((n, n))
    for (a, b), p in P.items():
        p = min(max(p, 1e-3), 1 - 1e-3)
        Y[idx[a], idx[b]] = np.log(p / (1 - p))
    s = Y.mean(axis=1)                        # least-squares potential on K_n
    G = s[:, None] - s[None, :]
    iu = np.triu_indices(n, 1)
    tot = float((Y[iu] ** 2).sum())
    res = float(((Y - G)[iu] ** 2).sum())
    curls = [abs(Y[i, j] + Y[j, k] + Y[k, i])
             for i, j, k in itertools.combinations(range(n), 3)]
    return (res / tot if tot > 0 else 0.0), float(max(curls) if curls else 0.0)


def bootstrap_cycles(classes, rng, n_boot):
    names = sorted(classes)
    arrs = {c: np.asarray(classes[c], float) for c in names}
    tri = list(itertools.combinations(names, 3))
    hits = {t: 0 for t in tri}
    any_cycle = 0
    for _ in range(n_boot):
        rs = {c: rng.choice(v, size=len(v), replace=True) for c, v in arrs.items()}
        _, P = pc.dominance(rs)
        cyc = cyclic_triangles(names, P)
        any_cycle += bool(cyc)
        for t in cyc:
            hits[t] += 1
    freq = {t: h / n_boot for t, h in hits.items()}
    return {
        "n_triangles": len(tri),
        "p_any_cycle": any_cycle / n_boot,
        "robust_cycles": [ "->".join(t) for t, f in freq.items() if f >= 0.5],
        "fragile_triangles": sum(1 for f in freq.values() if 0 < f < 0.5),
        "max_triangle_cycle_freq": float(max(freq.values()) if freq else 0.0),
        "worst_triangle": "->".join(max(freq, key=freq.get)) if freq and max(freq.values()) > 0 else None,
    }


def tournament_report(classes, rng, n_boot):
    names, P = pc.dominance(classes)
    cyc = cyclic_triangles(names, P)
    cf, mc = curl_stats(names, P)
    m = [abs(2 * P[p] - 1) for p in itertools.combinations(names, 2)]
    out = {"n_classes": len(names), "n_records": int(sum(len(v) for v in classes.values())),
           "majority_3cycles": len(cyc), "cycles": ["->".join(t) for t in cyc],
           "cyclic_energy_fraction": cf, "max_triangle_curl": mc,
           "resolvability": float(np.mean(m)),
           "near_tie_pairs": int(sum(v < 0.1 for v in m))}
    out["bootstrap"] = bootstrap_cycles(classes, rng, n_boot)
    return out


# ------------------------------------------------------------- campaign steps
def hodge_observed(W, n, near=None):
    """Hodge split of the empirical flow on the observed comparison graph.
    W[i, j] = wins of i over j. Returns energy fractions (gradient, curl,
    harmonic), the number of fully observed triangles, the number of those
    that are directed 3-cycles, and the observed edge count."""
    edges = [(i, j) for i in range(n) for j in range(i + 1, n) if W[i, j] + W[j, i] > 0]
    E = len(edges)
    if E == 0:
        return None
    eidx = {e: k for k, e in enumerate(edges)}
    y = np.array([np.log((W[i, j] + 0.5) / (W[j, i] + 0.5)) for i, j in edges])
    B1 = np.zeros((E, n))
    for k, (i, j) in enumerate(edges):
        B1[k, i], B1[k, j] = 1.0, -1.0
    tris = [t for t in itertools.combinations(range(n), 3)
            if (t[0], t[1]) in eidx and (t[1], t[2]) in eidx and (t[0], t[2]) in eidx]
    B2 = np.zeros((E, len(tris)))
    for c, (i, j, k) in enumerate(tris):
        B2[eidx[(i, j)], c] = 1.0
        B2[eidx[(j, k)], c] = 1.0
        B2[eidx[(i, k)], c] = -1.0
    tot = float(y @ y)
    if tot == 0:
        return {"grad": 0.0, "curl": 0.0, "harm": 0.0, "n_tri": len(tris),
                "n_cyc": 0, "n_cyc_near": 0, "n_edges": E}
    g = B1 @ np.linalg.lstsq(B1, y, rcond=None)[0]
    r = y - g
    cpart = B2 @ np.linalg.lstsq(B2, r, rcond=None)[0] if tris else np.zeros(E)
    h = r - cpart
    ncyc, ncyc_near = 0, 0
    for i, j, k in tris:
        d = [np.sign(W[i, j] - W[j, i]), np.sign(W[j, k] - W[k, j]), np.sign(W[k, i] - W[i, k])]
        if 0 not in d and d[0] == d[1] == d[2]:
            ncyc += 1
            if near is not None and ({(i, j), (j, k), (i, k)} & near):
                ncyc_near += 1
    return {"grad": float(g @ g) / tot, "curl": float(cpart @ cpart) / tot,
            "harm": float(h @ h) / tot, "n_tri": len(tris), "n_cyc": ncyc,
            "n_cyc_near": ncyc_near, "n_edges": E}


def near_tie_pairs(classes, thr=0.1):
    names, P = pc.dominance(classes)
    return {(i, j) for i, j in itertools.combinations(range(len(names)), 2)
            if abs(2 * P[(names[i], names[j])] - 1) < thr}


def campaign(classes, k, budget_max, rng, near=None):
    names = sorted(classes)
    n = len(names)
    arrs = [np.asarray(classes[c], float) for c in names]
    W = np.zeros((n, n))
    spent, out, bi = 0, {}, 0
    while spent < budget_max:
        kk = min(k, n)
        pick = rng.choice(n, size=kk, replace=False)
        draws = [(arrs[c][rng.integers(len(arrs[c]))] + 1e-9 * rng.random(), c) for c in pick]
        draws.sort(reverse=True)
        for a in range(kk):
            for b in range(a + 1, kk):
                W[draws[a][1], draws[b][1]] += 1
        spent += kk
        while bi < len(BUDGETS) and spent >= BUDGETS[bi]:
            out[BUDGETS[bi]] = hodge_observed(W, n, near)
            bi += 1
    return out


# --------------------------------------------------------- sources as players
def source_tournament(recs, sources, classes_filter=None, min_shared=30, rng=None, n_boot=0):
    fv = {}
    for s, f, c, g, _cs, _es in recs:
        if classes_filter is not None and c != classes_filter:
            continue
        fv.setdefault(f, {}).setdefault(s, g)
    P, nshared, outc = {}, {}, {}
    for a, b in itertools.combinations(sources, 2):
        sh = [d for d in fv.values() if a in d and b in d]
        if len(sh) < min_shared:
            continue
        o = np.array([1.0 if d[a] > d[b] else 0.5 if d[a] == d[b] else 0.0 for d in sh])
        w = float(o.mean())
        P[(a, b)], P[(b, a)] = w, 1 - w
        nshared[f"{a}|{b}"] = len(sh)
        outc[(a, b)] = o
    players = sorted({a for a, _ in P})
    cyc = [t for t in itertools.combinations(players, 3)
           if all(e in P for e in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])))
           and t in cyclic_triangles_sub(t, P)]
    closed = sum(1 for t in itertools.combinations(players, 3)
                 if all(e in P for e in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0]))))
    boot = {}
    if rng is not None and n_boot and cyc:
        # resample the shared compounds of every pair; how often is each
        # point-estimate cycle still a cycle?
        for t in cyc:
            hit = 0
            for _ in range(n_boot):
                Pb = {}
                for e in ((t[0], t[1]), (t[1], t[2]), (t[0], t[2])):
                    key = e if e in outc else (e[1], e[0])
                    o = outc[key]
                    w = float(rng.choice(o, size=len(o), replace=True).mean())
                    Pb[key], Pb[(key[1], key[0])] = w, 1 - w
                hit += bool(cyclic_triangles(list(t), Pb))
            boot["->".join(t)] = hit / n_boot
    # the same triangle on the compounds that all three sources report: each
    # compound then orders the three sources transitively, so a surviving
    # cycle would be a Condorcet (majority-vote) cycle; a vanishing one is a
    # coverage cycle created by comparing each pair on different compounds
    common = {}
    for t in cyc:
        cm = [d for d in fv.values() if all(x in d for x in t)]
        Pc = {}
        for a, b in itertools.permutations(t, 2):
            Pc[(a, b)] = float(np.mean([1.0 if d[a] > d[b] else 0.5 if d[a] == d[b] else 0.0
                                        for d in cm])) if cm else 0.5
        common["->".join(t)] = {"n_common": len(cm),
                                "cyclic_on_common_set": bool(cm) and bool(cyclic_triangles(list(t), Pc)),
                                "P_common": {f"{a}>{b}": Pc[(a, b)] for a, b in
                                             ((t[0], t[1]), (t[1], t[2]), (t[2], t[0]))}}
    edges = {}
    for t in cyc:
        for e in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
            edges[f"{e[0]}>{e[1]}"] = {"P": P[e], "n_shared": len(outc.get(e, outc.get((e[1], e[0]))))}
    return {"players": players, "n_pairs": len(nshared), "closed_triangles": closed,
            "majority_3cycles": len(cyc), "cycles": ["->".join(t) for t in cyc],
            "cycle_bootstrap_freq": boot, "cycle_edges": edges,
            "cycle_on_common_compounds": common,
            "wins": {a: sum(1 for (x, y), p in P.items() if x == a and p > 0.5) for a in players}}


def cyclic_triangles_sub(t, P):
    return cyclic_triangles(list(t), P)


# ------------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--universe", default="output/bandgap_archive_build/bandgap_universe.jsonl")
    ap.add_argument("--output-dir", default="output/local_cycles")
    ap.add_argument("--seeds", type=int, default=50)
    ap.add_argument("--n-boot", type=int, default=N_BOOT)
    args = ap.parse_args()
    root = Path(__file__).resolve().parent.parent
    universe = Path(args.universe)
    if not universe.is_absolute():
        universe = (root / universe).resolve()
    sha = hashlib.sha256(universe.read_bytes()).hexdigest()[:16]
    out = (root / args.output_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    recs = load(universe)
    sources = list(LABEL)
    rng = np.random.default_rng(0)
    res = {"universe_sha256_16": sha,
           "config": {"min_class": MIN_CLASS, "min_stratum_class": MIN_STRATUM_CLASS,
                      "n_boot": args.n_boot, "seeds": args.seeds, "budgets": BUDGETS}}

    def classes_where(src, pred, min_size):
        d = {}
        for s, f, c, g, cs, es in recs:
            if s == src and pred(cs, es):
                d.setdefault(c, []).append(g)
        return {c: sorted(v) for c, v in d.items() if len(v) >= min_size}

    # whole sources (global, with bootstrap) -------------------------------
    whole = {}
    for s in sources:
        cl = classes_where(s, lambda cs, es: True, MIN_CLASS)
        whole[s] = cl
        res.setdefault("global", {})[s] = tournament_report(cl, rng, args.n_boot)
        print("global", s, res["global"][s]["majority_3cycles"],
              res["global"][s]["bootstrap"]["p_any_cycle"], flush=True)

    # (L1)+(L2) strata ---------------------------------------------------------
    strata = {"crystal_system": [n for _, n in CRYSTAL],
              "element_count": ["binary", "ternary", "quaternary+"]}
    res["strata"] = {}
    for kind, values in strata.items():
        for v in values:
            for s in sources:
                if kind == "crystal_system":
                    pred = (lambda cs, es, v=v: cs == v)
                else:
                    pred = (lambda cs, es, v=v: es == v)
                cl = classes_where(s, pred, MIN_STRATUM_CLASS)
                if len(cl) < MIN_STRATUM_CLASSES:
                    continue
                rep = tournament_report(cl, rng, args.n_boot)
                res["strata"][f"{kind}|{v}|{s}"] = rep
                print(kind, v, s, rep["n_classes"], rep["majority_3cycles"],
                      rep["bootstrap"]["p_any_cycle"], flush=True)

    # (L3) campaign steps ----------------------------------------------------
    res["campaign"] = {}
    for k in (2, 4):
        agg = {b: {"grad": [], "curl": [], "harm": [], "cyc_frac": [], "any_cyc": [],
                   "edges": [], "n_cyc": [], "n_cyc_near": []} for b in BUDGETS}
        per_source = {}
        for s in sources:
            ps = {b: {"any_cyc": [], "cyc_frac": []} for b in BUDGETS}
            near = near_tie_pairs(whole[s])
            for seed in range(args.seeds):
                r = np.random.default_rng(1000 * k + seed + 17 * sources.index(s))
                steps = campaign(whole[s], k, BUDGETS[-1], r, near)
                for b, h in steps.items():
                    if h is None:
                        continue
                    a = agg[b]
                    a["grad"].append(h["grad"]); a["curl"].append(h["curl"])
                    a["harm"].append(h["harm"]); a["edges"].append(h["n_edges"])
                    a["n_cyc"].append(h["n_cyc"]); a["n_cyc_near"].append(h["n_cyc_near"])
                    frac = h["n_cyc"] / h["n_tri"] if h["n_tri"] else 0.0
                    a["cyc_frac"].append(frac); a["any_cyc"].append(float(h["n_cyc"] > 0))
                    ps[b]["any_cyc"].append(float(h["n_cyc"] > 0))
                    ps[b]["cyc_frac"].append(frac)
            per_source[s] = {str(b): {"p_any_cycle": float(np.mean(v["any_cyc"])),
                                      "mean_cyclic_triangle_fraction": float(np.mean(v["cyc_frac"]))}
                             for b, v in ps.items() if v["any_cyc"]}
        res["campaign"][f"k{k}"] = {
            "mean": {str(b): {**{key: float(np.mean(v)) for key, v in a.items() if v},
                              "share_cycles_touching_near_tie":
                                  float(np.sum(a["n_cyc_near"]) / max(1, np.sum(a["n_cyc"])))}
                     for b, a in agg.items()},
            "per_source": per_source}
        print("campaign k", k, {b: round(np.mean(a["any_cyc"]), 3) for b, a in agg.items()},
              flush=True)

    # (L4) sources as players ---------------------------------------------------
    fams = sorted({c for _, _, c, _, _, _ in recs})
    res["sources_as_players"] = {"all": source_tournament(recs, sources, rng=rng, n_boot=args.n_boot)}
    for c in fams:
        rep = source_tournament(recs, sources, classes_filter=c, rng=rng, n_boot=args.n_boot)
        if rep["closed_triangles"]:
            res["sources_as_players"][c] = rep

    (out / "local_cycles.json").write_text(json.dumps(res, indent=2))
    write_md(res, out)
    plot(res, out)
    return 0


def write_md(res, out):
    L = ["# Rock-paper-scissors cycles: global versus local", "",
         f"Archive sha256[:16] = `{res['universe_sha256_16']}`; bootstrap "
         f"{res['config']['n_boot']} resamples; campaign replays {res['config']['seeds']} seeds.", "",
         "## Whole data sources (global)", "",
         "| source | classes | 3-cycles | cyclic energy | max triangle curl | P(any cycle, bootstrap) | fragile triangles |",
         "|---|---|---|---|---|---|---|"]
    for s, r in res["global"].items():
        b = r["bootstrap"]
        L.append(f"| {LABEL[s]} | {r['n_classes']} | {r['majority_3cycles']} | "
                 f"{r['cyclic_energy_fraction']:.3f} | {r['max_triangle_curl']:.2f} | "
                 f"{b['p_any_cycle']:.3f} | {b['fragile_triangles']} / {b['n_triangles']} |")
    L += ["", "## Sub-populations (local, L1 and L2)", "",
          "| stratum | source | classes | records | 3-cycles | P(any cycle, bootstrap) | robust cycles | worst triangle (freq) |",
          "|---|---|---|---|---|---|---|---|"]
    tot = {"tables": 0, "with_cycle": 0, "robust": 0, "tri": 0}
    for key, r in res["strata"].items():
        kind, v, s = key.split("|")
        b = r["bootstrap"]
        tot["tables"] += 1; tot["tri"] += b["n_triangles"]
        tot["with_cycle"] += r["majority_3cycles"] > 0
        tot["robust"] += len(b["robust_cycles"])
        wt = f"{b['worst_triangle']} ({b['max_triangle_cycle_freq']:.2f})" if b["worst_triangle"] else "none"
        L.append(f"| {v} | {LABEL[s]} | {r['n_classes']} | {r['n_records']} | "
                 f"{r['majority_3cycles']} {'(' + ', '.join(r['cycles']) + ')' if r['cycles'] else ''} | "
                 f"{b['p_any_cycle']:.3f} | {', '.join(b['robust_cycles']) or 'none'} | {wt} |")
    L += ["", f"Summary: {tot['tables']} stratum tournaments, {tot['tri']} triangles; "
          f"{tot['with_cycle']} tournaments contain a point-estimate 3-cycle; "
          f"{tot['robust']} triangles are cyclic in at least half of the bootstrap resamples."]
    res["strata_summary"] = tot
    L += ["", "## Every step of a random campaign (L3)", "",
          "Hodge energy fractions of the empirical flow on the observed comparison graph; "
          "cycles are directed 3-cycles among fully observed triangles.", "",
          "| k | budget | observed pairs | gradient | local curl | global harmonic | triangles cyclic | P(any cycle) | cycles touching a near-tie pair |",
          "|---|---|---|---|---|---|---|---|---|"]
    for k, d in res["campaign"].items():
        for b, m in d["mean"].items():
            L.append(f"| {k[1:]} | {b} | {m['edges']:.1f} | {m['grad']:.3f} | {m['curl']:.3f} | "
                     f"{m['harm']:.3f} | {m['cyc_frac']:.3f} | {m['any_cyc']:.3f} | "
                     f"{m['share_cycles_touching_near_tie']:.2f} |")
    L += ["", "## Data sources as players (L4)", "",
          "| scope | sources | pairs with >= 30 shared compounds | closed triangles | 3-cycles (bootstrap persistence; cyclic on the common compounds?) |",
          "|---|---|---|---|---|"]
    for c, r in res["sources_as_players"].items():
        cy = ", ".join(f"{x} ({r['cycle_bootstrap_freq'].get(x, float('nan')):.2f}; "
                       f"{'yes' if r['cycle_on_common_compounds'][x]['cyclic_on_common_set'] else 'no'}, "
                       f"n = {r['cycle_on_common_compounds'][x]['n_common']})" for x in r["cycles"])
        L.append(f"| {c} | {len(r['players'])} | {r['n_pairs']} | {r['closed_triangles']} | "
                 f"{r['majority_3cycles']} {'(' + cy + ')' if cy else ''} |")
    L += ["", "Edges of the source-level cycles (P = share of shared compounds where the first source reports the wider gap):", ""]
    for c, r in res["sources_as_players"].items():
        for e, v in r["cycle_edges"].items():
            L.append(f"- {c}: {e} P = {v['P']:.3f} over {v['n_shared']} compounds")
    (out / "local_cycles.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


def plot(res, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "writings" / "papers_2026"))
    try:
        import paper_fig_style as st
        st.apply_style()
    except Exception:
        st = None
    fig, axes = plt.subplots(1, 2, figsize=(6.73, 2.5))
    colors = {"k2": "#0072B2", "k4": "#D55E00"}
    ax = axes[0]
    for k, d in res["campaign"].items():
        bs = [int(b) for b in d["mean"]]
        ax.plot(bs, [d["mean"][str(b)]["any_cyc"] for b in bs], "-o", ms=3,
                color=colors[k], label=f"{'pairs' if k == 'k2' else '4-sample batches'}")
    ax.set_xscale("log")
    ax.axvline(CAMPAIGN_BUDGET, color="#555555", lw=0.8, ls="--")
    ax.text(CAMPAIGN_BUDGET * 1.08, 0.04, "campaign\nbudget", fontsize=7, color="#555555")
    ax.set_xlabel("Evaluations spent (log scale)")
    ax.set_ylabel("Share of campaigns with a cycle")
    ax.set_ylim(0, 1.02)
    ax.set_title("(a) Cycles seen during a campaign", loc="left")
    ax.legend(frameon=False)
    ax = axes[1]
    d = res["campaign"]["k2"]["mean"]
    bs = [int(b) for b in d]
    g = np.array([d[str(b)]["grad"] for b in bs])
    c = np.array([d[str(b)]["curl"] for b in bs])
    h = np.array([d[str(b)]["harm"] for b in bs])
    ax.stackplot(bs, g, c, h, colors=["#9E9E9E", "#E69F00", "#CC79A7"],
                 labels=["consistent ranking", "local loops (triangles)",
                         "global loops (4+ classes)"])
    ax.set_xscale("log")
    ax.axvline(CAMPAIGN_BUDGET, color="#333333", lw=0.8, ls="--")
    ax.set_xlabel("Evaluations spent, pairs (log scale)")
    ax.set_ylabel("Share of comparison signal")
    ax.set_ylim(0, 1.0)
    ax.set_title("(b) Where the signal sits", loc="left")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(facecolor=cc, edgecolor="#333333", linewidth=0.5, label=ll) for cc, ll in
                       (("#9E9E9E", "consistent ranking"), ("#E69F00", "local loops (triangles)"),
                        ("#CC79A7", "global loops (4+ classes)"))],
              loc="lower right", frameon=True, facecolor="white", framealpha=1.0)
    fig.tight_layout()
    fig.savefig(out / "local_cycles.png", dpi=300)
    fig.savefig(out / "local_cycles.pdf")


if __name__ == "__main__":
    raise SystemExit(main())
