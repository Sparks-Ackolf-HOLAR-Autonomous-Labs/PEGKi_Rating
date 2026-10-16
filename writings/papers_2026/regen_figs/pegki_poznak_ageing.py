#!/usr/bin/env python3
"""PEGKi annex figure: measured age hardening of six Al-Mg-Si alloys
(Poznak et al. 2018) as a receiver for the computed sources.

Data: compute/almgsi_campaign/results/lab_receivers/poznak2018_ageing.json
and poznak_receiver.json.
(a) Peak hardness at 175 C against Mg/Si (weight ratio) for the two solute
    levels; filled = 0.5 h natural ageing before ageing, open = 4 h.
(b) Agreement (Kendall tau) of each source's precipitation driving force with
    the measured peak hardness, 90% interval from resampling alloys; the label
    gives the precipitate the source finds most favourable; dashed line =
    ordering by the beta'' (Mg5Si6) solute available, with no source energy.
Takeaway: the alloys harden in the order of beta'' solute available; sources
whose lowest-energy precipitate is Mg2Si rank the Mg-rich alloys too high.
Size and fonts: paper_fig_style, RSC full width at print size.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from paper_fig_style import RSC_FULL, MP_BLUE, JARVIS_GREEN, SNUMAT_PURPLE, ACCENT_RED, INK, GREY, apply_style, save  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
LAB = ROOT / "compute" / "almgsi_campaign" / "results" / "lab_receivers"
OUT = ROOT / "writings" / "papers_2026" / "pegki_digital_discovery" / "figures" / "poznak_ageing.png"
SHOW = [("qe_pbe", "QE PBE"), ("qe_hse06", "QE HSE06"), ("mace_relax", "MACE-MP-0 0 K"),
        ("mace_md", "MACE-MP-0 MD"), ("gpumd_nep89_md", "NEP89 MD"), ("lammps_meam_md", "MEAM MD")]
PHASE = {"mp_1367_Mg2Si": "Mg$_2$Si", "mp_568306_Mg5Si6": "$\\beta''$", "mp_542854_Mg9Si5": "$\\beta'$"}


def main():
    apply_style()
    data = json.loads((LAB / "poznak2018_ageing.json").read_text())
    res = json.loads((LAB / "poznak_receiver.json").read_text())
    fig, (a, b) = plt.subplots(1, 2, figsize=(RSC_FULL, 2.4), gridspec_kw={"width_ratios": [1, 1.15]})
    for lvl, col, lab in (("HS", MP_BLUE, "Mg+Si 1.4 wt%"), ("LS", ACCENT_RED, "Mg+Si 0.7 wt%")):
        al = sorted([x for x in data["alloys"] if x["name"].startswith(lvl)], key=lambda x: x["Mg_wt"] / x["Si_wt"])
        r = [x["Mg_wt"] / x["Si_wt"] for x in al]
        a.plot(r, [x["peak_HV_NA05"] for x in al], "-o", color=col, ms=4, lw=1, label=lab)
        a.plot(r, [x["peak_HV_NA4"] for x in al], "o", mfc="white", mec=col, ms=4, lw=0)
    a.axvline(5 / 6 * 24.305 / 28.085, color=GREY, ls=":", lw=0.8)
    a.text(5 / 6 * 24.305 / 28.085 * 1.06, 70, "$\\beta''$\nMg$_5$Si$_6$", fontsize=6.5, color=INK)
    a.set_xscale("log")
    a.set_xticks([0.25, 0.5, 1, 2])
    a.set_xticklabels(["0.25", "0.5", "1", "2"])
    a.set_xlabel("Mg/Si (weight ratio)")
    a.set_ylabel("peak hardness at 175 $^\\circ$C (HV)")
    a.set_title("(a) measured (open: 4 h natural ageing)", loc="left")
    a.legend(loc="center left", fontsize=6.5, bbox_to_anchor=(0.0, 0.62))

    rows = [(lab, res["donors"][s]) for s, lab in SHOW if s in res["donors"]]
    for y, (lab, d) in enumerate(rows):
        t, (lo, hi) = d["agreement_peak_HV"], d["agreement_peak_HV_ci90"]
        col = JARVIS_GREEN if t > 0.5 else (SNUMAT_PURPLE if t > 0 else ACCENT_RED)
        b.plot([lo, hi], [y, y], color=col, lw=2.2, alpha=0.45, solid_capstyle="butt")
        b.plot(t, y, "o", color=col, ms=4.5)
        ph = sorted({PHASE[p] for p in d["phase_used"].values()})
        b.text(1.08, y, "/".join(ph), va="center", fontsize=6.5, color=INK)
    base = res["composition_only_beta_dprime"]["agreement_peak_HV"]
    b.axvline(base, color=INK, ls="--", lw=0.8)
    b.text(base - 0.03, -0.75, "$\\beta''$ solute only", ha="right", va="center", fontsize=6.5)
    b.axvline(0, color=GREY, lw=0.6, ls=":")
    b.set_yticks(range(len(rows)))
    b.set_yticklabels([r[0] for r in rows])
    b.set_ylim(-1.1, len(rows) - 0.4)
    b.invert_yaxis()
    b.set_xlim(-1.05, 1.35)
    b.set_xticks([-1, -0.5, 0, 0.5, 1])
    b.set_xlabel("agreement with peak hardness (Kendall $\\tau$)")
    b.set_title("(b) driving-force ordering, 6 alloys", loc="left")
    b.grid(axis="y", visible=False)
    save(fig, OUT)


if __name__ == "__main__":
    main()
