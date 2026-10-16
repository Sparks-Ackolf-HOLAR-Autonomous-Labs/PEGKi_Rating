"""PEGKi figure: the added Al-Mg-Si receivers, Kendall tau of every donor.

Columns are receivers closer to the alloy than the DFT formation energies:
MACE-MP-0 MD at 300 K, quasi-harmonic free energies at 453 and 803 K,
CALPHAD (COST 507) hull energies at 453 K, the MP archive, vacancy binding,
and the measured enthalpies when more than three compounds match. Rows are
donors. MD donors are taken at the receiver temperature where one exists.
Cells show tau with the number of shared items in the column header; blank
where the donor is the receiver or shares fewer than five items.

Reads only committed artifacts under compute/almgsi_campaign/results.
Output: pegki_digital_discovery/figures/receiver_matrix.png (full width).
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import kendalltau

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
LIVE = HERE.parents[2] / "compute" / "almgsi_campaign"
sys.path.insert(0, str(LIVE))

import paper_fig_style as st  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import receiver_extension_study as rx  # noqa: E402

ROWS = [("mace_relax", "MACE-MP-0, 0 K"), ("qe_pbe", "QE PBE"), ("qe_hse06", "QE HSE06"),
        ("archive_mp", "MP archive"), ("mace_md", "MACE-MP-0 MD"), ("meam_md", "MEAM MD"),
        ("nep_md", "NEP89 MD")]


def md_name(row, T):
    eng = {"mace_md": "mace_md", "meam_md": "lammps_meam_md", "nep_md": "gpumd_nep89_md"}[row]
    if T == 300:
        return "mace_md_300" if eng == "mace_md" else eng
    return f"{eng}_{T}"


def tau_on(a, b):
    common = sorted(set(a) & set(b))
    if len(common) < 5:
        return None, len(common)
    return float(kendalltau([a[c] for c in common], [b[c] for c in common]).statistic), len(common)


def main():
    src, idx = rx.load_src()
    import re
    # binaries and ternaries only, as in the campaigns (elemental structures are references)
    multi = {c for c, r in idx.items() if len(set(re.findall(r"[A-Z][a-z]?", r["formula"]))) >= 2}
    val = {s: {c: v["value"] for c, v in t.items() if c in multi} for s, t in src.items()}
    hull = json.loads((LIVE / "results" / "analysis" / "hull_receiver.json").read_text(encoding="utf-8"))
    cols = [("MACE-MP-0 MD\n300 K", "mace_md_300", 300), ("Quasi-harmonic\n453 K", "qha_453", 453),
            ("Quasi-harmonic\n803 K", "qha_803", 803), ("CALPHAD hull\n453 K", "calphad_453", 453),
            ("MP archive", "archive_mp", 0)]
    M = np.full((len(ROWS), len(cols)), np.nan)
    ns = []
    for j, (_, rec, T) in enumerate(cols):
        nmax = 0
        for i, (row, _) in enumerate(ROWS):
            if rec.startswith("calphad"):
                key = {"mace_relax": "mace_relax", "qe_pbe": "qe_pbe", "qe_hse06": "qe_hse06",
                       "archive_mp": "archive_mp"}.get(row) or md_name(row, T)
                h = hull.get(str(T), {}).get(key)
                if h:
                    M[i, j] = h["tau"]
                    nmax = max(nmax, h["n"])
                continue
            d = md_name(row, T if T else 300) if row.endswith("_md") else row
            if d == rec or d not in val or rec not in val:
                continue
            t, n = tau_on(val[d], val[rec])
            if t is not None:
                M[i, j] = t
                nmax = max(nmax, n)
        ns.append(nmax)
    vac = json.loads((LIVE / "results" / "analysis" / "vacancy_binding_n3.json").read_text(encoding="utf-8"))
    vrec = "qe_pbe" if len(vac["sources"]["qe_pbe"]["binding_eV"]) >= 5 else "mace_relax"
    vcol = np.full(len(ROWS), np.nan)
    vmap = {"mace_relax": "mace_relax", "qe_pbe": "qe_pbe", "meam_md": "lammps_meam_relax", "nep_md": "gpumd_nep89_relax"}
    nv = 0
    for i, (row, _) in enumerate(ROWS):
        k = f"{vmap.get(row)}->{vrec}"
        if k in vac["tau"]:
            vcol[i] = vac["tau"][k]["tau"]
            nv = vac["tau"][k]["n"]
    M = np.column_stack([M, vcol])
    cols.append((f"Vacancy binding\n({'PBE' if vrec == 'qe_pbe' else 'MACE-MP-0'})", "vac", 0))
    ns.append(nv)

    st.apply() if hasattr(st, "apply") else None
    fig, ax = plt.subplots(figsize=(st.RSC_FULL, 2.55))
    im = ax.imshow(M, cmap="RdBu", vmin=-1, vmax=1, aspect="auto")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if np.isfinite(M[i, j]):
                ax.text(j, i, f"{M[i, j]:+.2f}", ha="center", va="center", fontsize=7,
                        color="white" if abs(M[i, j]) > 0.6 else st.INK)
    ax.set_xticks(range(len(cols)), [f"{c[0]}\n(n = {n})" for c, n in zip(cols, ns)], fontsize=7)
    ax.set_yticks(range(len(ROWS)), [r[1] for r in ROWS], fontsize=7)
    ax.tick_params(length=0)
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
    cb.set_label("Kendall $\\tau$ with the receiver", fontsize=7)
    cb.ax.tick_params(labelsize=7)
    fig.tight_layout()
    out = HERE.parent / "pegki_digital_discovery" / "figures" / "receiver_matrix.png"
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=600)
    json.dump({"columns": [c[0].replace("\n", " ") for c in cols], "n": ns, "rows": [r[1] for r in ROWS],
               "tau": [[None if not np.isfinite(x) else round(float(x), 3) for x in r] for r in M]},
              open(out.with_suffix(".json"), "w", encoding="utf-8"), indent=1)
    print(out)


if __name__ == "__main__":
    main()
