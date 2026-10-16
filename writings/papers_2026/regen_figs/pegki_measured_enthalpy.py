#!/usr/bin/env python3
"""PEGKi annex figure: the measured formation-enthalpy receiver over all
combinations of the beta, gamma and Mg2Si measurements.

Data: compute/almgsi_campaign/results/lab_receivers/measured_enthalpy_combinations.json.
(a) One row per compound; the grey band spans the two measurements, the
    measurement labels mark its ends; markers are the computed values of
    five sources (beta from each source's Al-Mg hull at x_Mg = 0.389).
(b) Mean absolute error of each source against the three compounds: the bar
    spans the minimum to maximum over the 8 measurement combinations.
Takeaway: the measurements disagree on the beta-gamma order, yet the order of
the sources by error is the same in all 8 combinations.
Size and fonts: paper_fig_style, RSC full width at print size.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from paper_fig_style import RSC_FULL, MP_BLUE, JARVIS_GREEN, SNUMAT_PURPLE, EXPT_GOLD, ACCENT_RED, INK, GREY, apply_style, save  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "compute" / "almgsi_campaign" / "results" / "lab_receivers" / "measured_enthalpy_combinations.json"
OUT = ROOT / "writings" / "papers_2026" / "pegki_digital_discovery" / "figures" / "measured_enthalpy.png"
KJ = 96.485332

SHOW = [("qe_pbe", "QE PBE", MP_BLUE, "o"), ("qe_hse06", "QE HSE06", SNUMAT_PURPLE, "s"),
        ("archive_mp", "MP archive", GREY, "D"), ("mace_relax", "MACE-MP-0", JARVIS_GREEN, "^"),
        ("gpumd_nep89_md", "NEP89 MD", ACCENT_RED, "v")]
COMPOUNDS = [("Mg2Si", "Mg$_2$Si", "Mg2Si"), ("gamma", "$\\gamma$ Al$_{12}$Mg$_{17}$", "gamma"),
             ("beta", "$\\beta$ Al$_3$Mg$_2$", "beta_hull")]
MLABEL = {"BP1970": "B&P", "PH1978": "P&H", "JANAF": "JANAF", "Feufel1997": "Feufel"}


def main():
    apply_style()
    d = json.loads(DATA.read_text())
    meas, comp = d["measured_kJ_per_mol_atoms"], d["computed_eV_per_atom"]
    fig, (a, b) = plt.subplots(1, 2, figsize=(RSC_FULL, 2.5), gridspec_kw={"width_ratios": [1.35, 1]})
    for y, (key, lab, ckey) in enumerate(COMPOUNDS):
        vals = {k: v / KJ for k, v in meas[key].items()}
        lo, hi = min(vals.values()), max(vals.values())
        a.fill_betweenx([y - 0.28, y + 0.28], lo, hi, color=EXPT_GOLD, alpha=0.35, lw=0)
        for i, (s, sl, col, mk) in enumerate(SHOW):
            a.scatter(comp[s][ckey], y + (i - 2) * 0.09, color=col, marker=mk, s=16, zorder=3,
                      label=sl if y == 0 else None, lw=0)
    a.axvline(0, color=GREY, lw=0.6, ls=":")
    a.set_yticks(range(len(COMPOUNDS)))
    a.set_yticklabels([c[1] for c in COMPOUNDS])
    a.set_ylim(-0.6, len(COMPOUNDS) - 0.25)
    a.set_xlabel("formation enthalpy (eV/atom)")
    a.set_title("(a) measured band and computed values", loc="left")
    a.legend(loc="upper left", fontsize=6, handletextpad=0.2, borderaxespad=0.2)
    a.grid(axis="y", visible=False)

    summ = d["summary"]["beta_hull"]
    order = sorted(summ, key=lambda s: summ[s]["mae_min"])
    for y, s in enumerate(order):
        r = summ[s]
        b.plot([r["mae_min"], r["mae_max"]], [y, y], color=MP_BLUE, lw=4, solid_capstyle="butt")
    b.set_yticks(range(len(order)))
    b.set_yticklabels([summ[s]["label"] for s in order], fontsize=6.5)
    b.invert_yaxis()
    b.set_xlim(0, None)
    b.set_xlabel("mean absolute error (eV/atom)")
    b.set_title("(b) range over 8 combinations", loc="left")
    b.grid(axis="y", visible=False)
    save(fig, OUT)


if __name__ == "__main__":
    main()
