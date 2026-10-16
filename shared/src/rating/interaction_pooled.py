"""
Interaction-pooled Plackett-Luce: the additive pooled model as a cold start

The pooled model of hypergraph_rating.PooledPlackettLuce is additive by
construction:

    log w[source, class] = theta[class] + delta[source]

so every source is forced to share one class ordering AND one set of class
margins, and a source can only shift everything up or down. Nothing in the
data requires that. Here the additive form is kept as the starting point
(cold start) and the data are allowed to move away from it:

    log w[c, i] = a_c * theta_i + delta_c + sum_r U[c, r] * V[i, r]
    a_c = exp(b_c)

  * a_c  : per-source margin scale. a_c > 1 means source c separates the
           classes more widely than the pooled consensus; a_c < 1 means it
           compresses them. The additive model is a_c = 1.
  * U V' : rank-R source x class interaction. It lets a source reorder
           classes, and it lets the cross-source offset depend on the class
           (the family-resolved offsets of scripts/family_offsets.py are the
           full-rank version of this). The additive model is U = V = 0.

The departure from the additive model is penalised,

    pen = lam_scale * sum_c b_c^2 + lam_int * (||U||^2 + ||V||^2)
          + ridge * (||theta||^2 + ||delta||^2),

so lam -> infinity recovers the additive fit exactly (the cold start), and
lam is chosen on held-out events, i.e. learned from the data rather than
assumed. The fit maximises the Plackett-Luce log-likelihood minus pen with
L-BFGS and an analytic gradient.

Entities use the same "source::class" keys as PooledPlackettLuce, so the same
within-source events and cross-source bridge events feed both models.

Usage:
    uv run python -m src.rating.interaction_pooled     # self-test
"""

from __future__ import annotations

from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import minimize

from .hypergraph_rating import Hyperedge, PooledPlackettLuce

SEP = PooledPlackettLuce.SEP


class _EventTensor:
    """Padded (E x K) index arrays for vectorised Plackett-Luce."""

    def __init__(self, edges: Sequence[Hyperedge], cidx: Dict[str, int],
                 iidx: Dict[str, int]):
        K = max(len(e.ranking) for e in edges)
        E = len(edges)
        self.c = np.zeros((E, K), dtype=int)
        self.i = np.zeros((E, K), dtype=int)
        self.mask = np.zeros((E, K), dtype=bool)
        self.w = np.array([e.weight for e in edges], dtype=float)
        self.k = np.array([len(e.ranking) for e in edges], dtype=int)
        for r, e in enumerate(edges):
            for j, key in enumerate(e.ranking):
                ch, cl = key.split(SEP, 1)
                self.c[r, j] = cidx[ch]
                self.i[r, j] = iidx[cl]
                self.mask[r, j] = True
        # stage j is a scored stage iff j < k - 1
        self.stage = np.arange(K)[None, :] < (self.k[:, None] - 1)


def _pl_loglik_and_grad_s(S: np.ndarray, T: _EventTensor):
    """Per-event PL log-likelihood and d ll / d S (both E-shaped / E x K).

    S: E x K log-strengths (padding entries ignored).
    """
    Sm = np.where(T.mask, S, -np.inf)
    # reverse cumulative logsumexp: L[:, j] = logsumexp(S[:, j:])
    L = np.logaddexp.accumulate(Sm[:, ::-1], axis=1)[:, ::-1]
    stage = T.stage
    with np.errstate(invalid="ignore"):
        ll = np.where(stage, Sm - L, 0.0).sum(axis=1)
    # d ll / d S_m = 1[m is a stage winner] - sum_{j <= m, stage j} exp(S_m - L_j)
    invL = np.where(stage, np.exp(-np.where(stage, L, 0.0)), 0.0)
    cum = np.cumsum(invL, axis=1)
    with np.errstate(over="ignore", invalid="ignore"):
        prob_sum = np.where(T.mask, np.exp(np.where(T.mask, S, 0.0)) * cum, 0.0)
    g = stage.astype(float) - prob_sum
    return ll, g


class InteractionPooledPlackettLuce:
    """Pooled Plackett-Luce with learnable departures from additivity.

    rank: number of interaction factors R (0 = no interaction term).
    fit_scale: learn per-source margin scales a_c (False = a_c fixed at 1).
    lam_scale, lam_int: penalties pulling the fit back to the additive
        cold start. Large values reproduce the additive model.
    """

    def __init__(self, rank: int = 0, fit_scale: bool = True,
                 lam_scale: float = 1.0, lam_int: float = 1.0,
                 ridge: float = 1e-3, seed: int = 0):
        self.rank = int(rank)
        self.fit_scale = bool(fit_scale)
        self.lam_scale = float(lam_scale)
        self.lam_int = float(lam_int)
        self.ridge = float(ridge)
        self.seed = int(seed)
        self.sources: List[str] = []
        self.classes: List[str] = []
        self.theta = None
        self.delta = None
        self.b = None
        self.U = None
        self.V = None
        self.converged = False
        self.n_iterations = 0

    # -- parameter packing -------------------------------------------------
    def _unpack(self, x):
        C, I, R = len(self.sources), len(self.classes), self.rank
        o = 0
        theta = x[o:o + I]; o += I
        delta = x[o:o + C]; o += C
        if self.fit_scale:
            b = x[o:o + C]; o += C
        else:
            b = np.zeros(C)
        U = x[o:o + C * R].reshape(C, R); o += C * R
        V = x[o:o + I * R].reshape(I, R); o += I * R
        return theta, delta, b, U, V

    def _pack(self, theta, delta, b, U, V):
        parts = [theta, delta]
        if self.fit_scale:
            parts.append(b)
        parts += [U.ravel(), V.ravel()]
        return np.concatenate(parts)

    def log_strength_matrix(self, theta=None, delta=None, b=None, U=None,
                            V=None) -> np.ndarray:
        theta = self.theta if theta is None else theta
        delta = self.delta if delta is None else delta
        b = self.b if b is None else b
        U = self.U if U is None else U
        V = self.V if V is None else V
        return (np.exp(b)[:, None] * theta[None, :] + delta[:, None]
                + U @ V.T)

    # -- objective ---------------------------------------------------------
    def _objective(self, x, T: _EventTensor):
        theta, delta, b, U, V = self._unpack(x)
        M = self.log_strength_matrix(theta, delta, b, U, V)
        S = M[T.c, T.i]
        ll, g = _pl_loglik_and_grad_s(S, T)
        g = g * T.w[:, None]
        f = -(T.w * ll).sum()
        # scatter dL/dS to the C x I matrix
        G = np.zeros_like(M)
        np.add.at(G, (T.c[T.mask], T.i[T.mask]), g[T.mask])
        a = np.exp(b)
        g_theta = (G * a[:, None]).sum(axis=0)
        g_delta = G.sum(axis=1)
        g_b = (G * theta[None, :]).sum(axis=1) * a
        g_U = G @ V
        g_V = G.T @ U
        # penalties (on the negative log-likelihood)
        f += self.ridge * (theta @ theta + delta @ delta)
        f += self.lam_scale * (b @ b)
        f += self.lam_int * ((U * U).sum() + (V * V).sum())
        grad = -self._pack(g_theta, g_delta, g_b, g_U, g_V)
        grad += self._pack(2 * self.ridge * theta, 2 * self.ridge * delta,
                           2 * self.lam_scale * b, 2 * self.lam_int * U,
                           2 * self.lam_int * V)
        return f, grad

    # -- fitting -----------------------------------------------------------
    def fit(self, hyperedges: Iterable[Hyperedge],
            init: Optional["InteractionPooledPlackettLuce"] = None,
            max_iter: int = 2000) -> "InteractionPooledPlackettLuce":
        edges = list(hyperedges)
        if not edges:
            raise ValueError("No hyperedges to fit")
        keys = [k for e in edges for k in e.ranking]
        for k in keys:
            if SEP not in k:
                raise ValueError(f"entities must be 'source{SEP}class', got {k!r}")
        self.sources = sorted({k.split(SEP, 1)[0] for k in keys})
        self.classes = sorted({k.split(SEP, 1)[1] for k in keys})
        cidx = {c: j for j, c in enumerate(self.sources)}
        iidx = {c: j for j, c in enumerate(self.classes)}
        T = _EventTensor(edges, cidx, iidx)
        C, I, R = len(self.sources), len(self.classes), self.rank
        rng = np.random.default_rng(self.seed)
        if init is not None and init.sources == self.sources \
                and init.classes == self.classes:
            theta, delta = init.theta.copy(), init.delta.copy()
            b = init.b.copy() if self.fit_scale else np.zeros(C)
            U0 = np.zeros((C, R)); V0 = np.zeros((I, R))
            r0 = min(R, init.rank)
            if r0:
                U0[:, :r0] = init.U[:, :r0]; V0[:, :r0] = init.V[:, :r0]
            U = U0 + 0.01 * rng.standard_normal((C, R))
            V = V0 + 0.01 * rng.standard_normal((I, R))
        else:
            theta, delta, b = np.zeros(I), np.zeros(C), np.zeros(C)
            U = 0.01 * rng.standard_normal((C, R))
            V = 0.01 * rng.standard_normal((I, R))
        x0 = self._pack(theta, delta, b, U, V)
        res = minimize(self._objective, x0, args=(T,), jac=True,
                       method="L-BFGS-B",
                       options={"maxiter": max_iter, "gtol": 1e-7})
        self.converged = bool(res.success)
        self.n_iterations = int(res.nit)
        self.theta, self.delta, self.b, self.U, self.V = self._unpack(res.x)
        # gauge: delta zero-mean (theta shift is absorbed by delta via a_c)
        return self

    # -- evaluation --------------------------------------------------------
    def event_logliks(self, hyperedges: Iterable[Hyperedge]) -> np.ndarray:
        """Per-event log-likelihood; unseen sources/classes raise KeyError."""
        edges = list(hyperedges)
        cidx = {c: j for j, c in enumerate(self.sources)}
        iidx = {c: j for j, c in enumerate(self.classes)}
        T = _EventTensor(edges, cidx, iidx)
        S = self.log_strength_matrix()[T.c, T.i]
        ll, _ = _pl_loglik_and_grad_s(S, T)
        return ll

    @property
    def margin_scales(self) -> Dict[str, float]:
        return {c: float(np.exp(self.b[j])) for j, c in enumerate(self.sources)}

    def source_class_matrix(self) -> Dict[str, Dict[str, float]]:
        M = self.log_strength_matrix()
        return {c: {k: float(M[j, m]) for m, k in enumerate(self.classes)}
                for j, c in enumerate(self.sources)}

    def n_parameters(self) -> int:
        C, I, R = len(self.sources), len(self.classes), self.rank
        return I + C + (C if self.fit_scale else 0) + R * (C + I)


def _self_test() -> int:
    failures = 0
    rng = np.random.default_rng(3)
    sources = ["s0", "s1", "s2", "s3"]
    classes = ["A", "B", "C", "D", "E"]
    theta = np.array([1.5, 0.8, 0.0, -0.7, -1.6])
    delta = np.array([0.6, -0.2, 0.3, -0.7])

    def sample(M, n_within=150, n_bridge=300):
        edges = []
        for j, s in enumerate(sources):
            for _ in range(n_within):
                g = M[j] + rng.gumbel(size=len(classes))
                order = np.argsort(-g)
                edges.append(Hyperedge(ranking=tuple(
                    f"{s}{SEP}{classes[m]}" for m in order)))
        for _ in range(n_bridge):
            m = int(rng.integers(len(classes)))
            g = M[:, m] + rng.gumbel(size=len(sources))
            order = np.argsort(-g)
            edges.append(Hyperedge(ranking=tuple(
                f"{sources[j]}{SEP}{classes[m]}" for j in order)))
        return edges

    # 1. gradient check
    M_add = theta[None, :] + delta[:, None]
    edges = sample(M_add, 20, 40)
    mdl = InteractionPooledPlackettLuce(rank=2, lam_scale=0.3, lam_int=0.3)
    mdl.sources, mdl.classes = sources, classes
    T = _EventTensor(edges, {c: j for j, c in enumerate(sources)},
                     {c: j for j, c in enumerate(classes)})
    x = rng.standard_normal(mdl._pack(np.zeros(5), np.zeros(4), np.zeros(4),
                                      np.zeros((4, 2)), np.zeros((5, 2))).size) * 0.3
    f0, g0 = mdl._objective(x, T)
    num = np.zeros_like(x)
    for q in range(x.size):
        e = np.zeros_like(x); e[q] = 1e-6
        num[q] = (mdl._objective(x + e, T)[0] - mdl._objective(x - e, T)[0]) / 2e-6
    err = float(np.max(np.abs(num - g0)) / max(1.0, np.max(np.abs(num))))
    ok = err < 1e-4
    print(f"  [1] analytic vs numeric gradient: rel err = {err:.2e}")
    failures += 0 if ok else 1

    # 2. cold start: heavy penalty reproduces the additive MM fit
    edges = sample(M_add)
    add = PooledPlackettLuce(prior_weight=1e-6).fit(edges, max_iter=5000)
    cold = InteractionPooledPlackettLuce(rank=2, lam_scale=1e6, lam_int=1e6,
                                         ridge=1e-8).fit(edges)
    d_add = np.array([add.channel_offsets[s] for s in sources])
    d_cold = cold.delta - cold.delta.mean()
    t_add = np.array([add.policy_log_strengths[c] for c in classes])
    t_add = t_add - t_add.mean()
    t_cold = cold.theta - cold.theta.mean()
    dev = float(max(np.max(np.abs(d_add - d_cold)), np.max(np.abs(t_add - t_cold))))
    ok = dev < 2e-2
    print(f"  [2] heavy-penalty fit = additive MM fit: max dev = {dev:.2e}")
    failures += 0 if ok else 1

    # 3. recovers a margin scale when one source stretches margins
    a_true = np.array([1.0, 2.0, 1.0, 0.5])
    M_sc = a_true[:, None] * theta[None, :] + delta[:, None]
    edges = sample(M_sc, 400, 400)
    sc = InteractionPooledPlackettLuce(rank=0, lam_scale=0.01).fit(edges)
    a_hat = np.array([sc.margin_scales[s] for s in sources])
    ratio = a_hat / np.exp(np.mean(np.log(a_hat)))
    ratio_true = a_true / np.exp(np.mean(np.log(a_true)))
    err = float(np.max(np.abs(np.log(ratio) - np.log(ratio_true))))
    ok = err < 0.2
    print(f"  [3] margin scales recovered: true {np.round(ratio_true, 2)} "
          f"fit {np.round(ratio, 2)} (max log err {err:.3f})")
    failures += 0 if ok else 1

    # 4. held-out: scale model beats additive on stretched data, and does
    #    not lose on additive data
    def heldout_gain(M):
        tr, te = sample(M, 300, 300), sample(M, 100, 100)
        a0 = InteractionPooledPlackettLuce(rank=0, fit_scale=False).fit(tr)
        a1 = InteractionPooledPlackettLuce(rank=0, lam_scale=0.1).fit(tr)
        return float(np.mean(a1.event_logliks(te) - a0.event_logliks(te)))
    g_st, g_add = heldout_gain(M_sc), heldout_gain(M_add)
    ok = g_st > 0.02 and g_add > -0.01
    print(f"  [4] held-out gain of scale model: stretched data {g_st:+.4f}, "
          f"additive data {g_add:+.4f} per event")
    failures += 0 if ok else 1

    print(f"self-test {'PASS' if failures == 0 else 'FAIL'} ({failures} failures)")
    return failures


if __name__ == "__main__":
    raise SystemExit(_self_test())
