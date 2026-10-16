# Learned transfer-score weights and routing (fixed choices as cold start)

Archive sha256[:16] = `f2158ecdc78e0229`. Target: measured early transfer gain T_early of the 50-seed paired campaign (TENKi exp 26, CHPC, fixed-trust full reuse). All numbers are leave-one-receiver-out.

## Query arm: pair

| quantity | fixed weights (1,1,1) | learned weights |
|---|---|---|
| mean within-receiver Spearman (donor ordering) | +0.342 | +0.349 |
| pooled Spearman over 56 cells | +0.317 | +0.361 |
| sign accuracy (transfer helps or hurts) | 0.875 raw, 0.875 calibrated | 0.875 |
| mean squared error of T_early | 0.0247 (calibrated) | 0.0243 |
| receivers whose best donor is ranked first | 2/8 | 3/8 |

Routing value (mean measured T_early obtained by the routing rule on held-out receivers):

| rule | value |
|---|---|
| lp_fixed | +0.2306 |
| lp_learned | +0.2093 |
| learned_from_lp_fixed_score | +0.3397 |
| learned_from_lp_learned_score | +0.3475 |
| uniform | +0.2311 |
| oracle | +0.4243 |

Weights fitted on all 56 cells (lam = 0.01 by leave-one-receiver-out): normalised (alpha, gamma, delta) = (+0.429, -0.073, +0.498); the cold start is (+0.333, +0.333, +0.333).

## Query arm: kway

| quantity | fixed weights (1,1,1) | learned weights |
|---|---|---|
| mean within-receiver Spearman (donor ordering) | +0.360 | +0.306 |
| pooled Spearman over 56 cells | +0.254 | +0.388 |
| sign accuracy (transfer helps or hurts) | 0.786 raw, 0.786 calibrated | 0.786 |
| mean squared error of T_early | 0.0199 (calibrated) | 0.0182 |
| receivers whose best donor is ranked first | 2/8 | 2/8 |

Routing value (mean measured T_early obtained by the routing rule on held-out receivers):

| rule | value |
|---|---|
| lp_fixed | +0.1467 |
| lp_learned | +0.0991 |
| learned_from_lp_fixed_score | +0.2340 |
| learned_from_lp_learned_score | +0.2432 |
| uniform | +0.1467 |
| oracle | +0.3073 |

Weights fitted on all 56 cells (lam = 0.01 by leave-one-receiver-out): normalised (alpha, gamma, delta) = (+0.371, -0.147, +0.482); the cold start is (+0.333, +0.333, +0.333).

