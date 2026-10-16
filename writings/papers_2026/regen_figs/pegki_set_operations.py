#!/usr/bin/env python3
"""PEGKi figure: set-operation studies on the band-gap archive.

Panels: (a) for each fidelity step, the transfer score of the donor into the
region where one rigid scissor shift reproduces the receiver gap within
0.3 eV minus its score into the region where it does not; (b) for each of
the six regions, the score of the best of the eight donors minus that of
the donor that defines the split, labelled with the best donor. Bars: 95%
intervals over 50 record-bootstrap replicates, as stored in
shared/output/set_operation_study/set_operation_study.json.

Size and fonts: RSC column (3.27 in), paper_fig_style fonts.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from paper_fig_style import (  # noqa: E402
    RSC_COL, PEGKI_SOURCE_NAMES, SOURCE_COLOURS, GREY, INK, apply_style, save,
)
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
ART = ROOT / "shared/output/set_operation_study/set_operation_study.json"
OUT = ROOT / "writings/papers_2026/pegki_digital_discovery/figures/set_operations.png"

TITLES = {("mp_gga", "expt_zhuo"): "MP GGA to experiment",
          ("snumat_gga", "snumat_hse"): "SNUMAT GGA to HSE06",
          ("jarvis_optb88vdw", "jarvis_mbj"): "OptB88vdW to TB-mBJ"}


def main():
    apply_style()
    d = json.load(open(ART, encoding="utf-8"))
    h1 = d["h1_ownership"]
    h2 = {h["receiver_study"]: h for h in d["h2_rerouting"]}
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(RSC_COL, 3.9),
                                 gridspec_kw={"height_ratios": [1.0, 1.7]})
    # (a) agreement minus disagreement region, per step
    for i, h in enumerate(h1):
        a, x = h["donor"], h["receiver"]
        lo, hi = h["ci95"]
        a1.errorbar([h["difference"]], [i], xerr=[[h["difference"] - lo], [hi - h["difference"]]],
                    fmt="o", ms=4.5, color=SOURCE_COLOURS[a], ecolor=SOURCE_COLOURS[a],
                    elinewidth=1.1, capsize=2)
    a1.set_yticks(range(len(h1)))
    a1.set_yticklabels([f"{TITLES[(h['donor'], h['receiver'])]}\n"
                        f"({h['records_agree']:,} / {h['records_disagree']:,} records)" for h in h1])
    a1.invert_yaxis()
    a1.axvline(0, color=INK, lw=0.7)
    a1.set_xlim(-0.75, 0.25)
    a1.set_xlabel(r"$\Delta$(agree region) $-$ $\Delta$(disagree region)")
    a1.set_title("(a) Same donor, two regions", loc="left", fontsize=8)
    a1.grid(axis="y", visible=False)
    # (b) best donor minus defining donor, per region
    rows = []
    for h in h1:
        for kind in ("agree", "disagree"):
            rows.append((h, kind, h2[h[kind]]))
    for i, (h, kind, s) in enumerate(rows):
        lo, hi = s["ci95"]
        v = s["best_minus_defining"]
        a2.errorbar([v], [i], xerr=[[v - lo], [hi - v]], fmt="D", ms=4,
                    color=SOURCE_COLOURS[s["best_donor"]], ecolor=SOURCE_COLOURS[s["best_donor"]],
                    elinewidth=1.1, capsize=2)
        a2.annotate(PEGKI_SOURCE_NAMES[s["best_donor"]], (hi, i), xytext=(4, 0),
                    textcoords="offset points", va="center", fontsize=7, color=INK)
    a2.set_yticks(range(len(rows)))
    a2.set_yticklabels([f"{TITLES[(h['donor'], h['receiver'])].split(' to ')[1]}, {kind}"
                        for h, kind, _ in rows])
    a2.invert_yaxis()
    a2.axvline(0, color=INK, lw=0.7)
    a2.set_xlim(-0.1, 0.72)
    a2.set_xticks([0, 0.2, 0.4])
    a2.set_xlabel(r"$\Delta$(best donor) $-$ $\Delta$(splitting donor)")
    a2.set_title("(b) Best of eight donors per region", loc="left", fontsize=8)
    a2.grid(axis="y", visible=False)
    save(fig, OUT)


if __name__ == "__main__":
    main()
