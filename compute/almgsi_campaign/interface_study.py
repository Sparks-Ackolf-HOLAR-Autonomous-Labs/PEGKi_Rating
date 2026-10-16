"""Al-Mg-Si contact on a Si device: does Si leave the substrate?

Uses records already computed for the live Al-Mg-Si campaign (108-site
Al solute cells, elemental references, Mg2Si) to ask two device
questions for an Al-Mg-Si metallisation, bond pad, or heat sink placed
on a Si wafer and annealed:

1. Si dissolution (junction spiking). The enthalpy to move one Si atom
   from the diamond-Si substrate into fcc Al as a dilute solute,
       dH_Si = E(Al107Si) - E(Al108) + e(Al) - e(Si),
   sets the equilibrium Si content of the contact, x_Si ~ exp(-dH_Si/kT)
   (configurational entropy only; vibrational entropy left out, so x_Si
   is an ordering device across sources, not a solubility prediction).
2. Interfacial Mg2Si. With Mg dissolved in the Al, the reaction
       2 Mg(in Al) + Si(substrate) -> Mg2Si
   has energy per Si atom
       dE_rxn = Ef_fu(Mg2Si) - 2 dH_Mg,
   dH_Mg = E(Al107Mg) - E(Al108) + e(Al) - e(Mg).
   dE_rxn < 0 means the contact draws Si out of the substrate and grows
   an interlayer of Mg2Si, itself a narrow-gap semiconductor.
   A third number, dH_Si - B(Mg-Si 1nn), is the dissolution enthalpy of
   a Si atom that pairs with a Mg solute (B > 0 attractive).

Every source that computed the needed cells enters: QE PBE, GPAW r2SCAN
(single points), MACE-MP-0 relaxation, and the three MD potentials at
300, 453 and 803 K. The band-gap side (Mg2Si gap per source) is read
from the computed semiconductor campaign.

Writes results/analysis/interface_study.{json,md}.
"""
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
RES = HERE / "results"
SEMI = HERE.parent / "semiconductor_campaign" / "results"
KB = 8.617333262e-5
T_ANNEAL = 723.0  # 450 C, a typical contact sinter anneal

STATIC = ["qe_pbe", "gpaw_r2scan", "mace_relax"]
MD = ["mace_md", "lammps_meam_md", "gpumd_nep89_md"]
LABEL = {"qe_pbe": "QE PBE", "gpaw_r2scan": "GPAW r2SCAN", "mace_relax": "MACE-MP-0, 0 K",
         "mace_md": "MACE-MP-0 MD", "lammps_meam_md": "MEAM MD", "gpumd_nep89_md": "NEP89 MD"}
REFS = {"Al": "mp_134_Al", "Si": "mp_149_Si", "Mg": "mp_153_Mg"}


def load(d, cid):
    p = d / f"{cid}.json"
    if not p.exists():
        return None
    r = json.loads(p.read_text())
    return r if r.get("status") == "ok" and r.get("energy_per_atom") is not None else None


def tot(r):
    s = (r.get("energy_per_atom_sem") or 0.0)
    return r["energy_per_atom"] * r["natoms"], s * r["natoms"]


def one(d):
    need = {k: load(d, c) for k, c in REFS.items()}
    cells = {c: load(d, f"clu_{c}") for c in ("Al108", "Si1", "Mg1", "MgSi_1nn")}
    mg2si = load(d, "mp_1367_Mg2Si")
    if any(v is None for v in need.values()) or cells["Al108"] is None:
        return None
    e = {k: v["energy_per_atom"] for k, v in need.items()}
    se = {k: (v.get("energy_per_atom_sem") or 0.0) for k, v in need.items()}
    h, sh = tot(cells["Al108"])
    out = {}
    for el, c in (("Si", "Si1"), ("Mg", "Mg1")):
        if cells[c] is None:
            continue
        t, st = tot(cells[c])
        out[f"dH_{el}"] = t - h + e["Al"] - e[el]
        out[f"dH_{el}_sem"] = math.sqrt(st ** 2 + sh ** 2 + se["Al"] ** 2 + se[el] ** 2)
    if mg2si is not None:
        ef = 3 * mg2si["energy_per_atom"] - 2 * e["Mg"] - e["Si"]
        out["Ef_Mg2Si_fu"] = ef
        if "dH_Mg" in out:
            out["dE_rxn"] = ef - 2 * out["dH_Mg"]
    if cells["MgSi_1nn"] is not None and cells["Si1"] is not None and cells["Mg1"] is not None:
        b = -(tot(cells["MgSi_1nn"])[0] - h) + (tot(cells["Si1"])[0] - h) + (tot(cells["Mg1"])[0] - h)
        out["B_MgSi_1nn"] = b
        out["dH_Si_paired"] = out["dH_Si"] - b
    if "dH_Si" in out:
        out["x_Si_723K"] = math.exp(-out["dH_Si"] / (KB * T_ANNEAL))
    return out


def gap_side():
    cid = "arc_Mg2Si_mp-1367"
    g = {}
    for s in ("qe_pbe", "qe_hse06"):
        p = SEMI / s / f"{cid}.json"
        if p.exists():
            g[s] = json.loads(p.read_text()).get("gap_eV")
    arch = json.loads((SEMI / "archive_donors.json").read_text())
    for k, v in arch.items():
        if cid in v:
            g[k] = v[cid]["value"]
    return g


def main():
    rows = {}
    for s in STATIC:
        r = one(RES / s)
        if r:
            rows[f"{s}|0K"] = r
    for T in ("T300", "T453", "T803"):
        for s in MD:
            r = one(RES.parent / "results_T" / T / s) if T != "T300" or (RES.parent / "results_T" / T / s).exists() else None
            if r:
                rows[f"{s}|{T[1:]}K"] = r
    gaps = gap_side()
    out = {"T_anneal_K": T_ANNEAL, "rows": rows, "mg2si_gap_eV": gaps}
    (RES / "analysis").mkdir(exist_ok=True)
    (RES / "analysis" / "interface_study.json").write_text(json.dumps(out, indent=1))
    L = ["# Al-Mg-Si contact on Si: interface energetics per source", "",
         "dH in eV per atom moved; dE_rxn in eV per Si taken from the substrate; "
         "x_Si = exp(-dH_Si/kT) at 723 K (configurational entropy only).", "",
         "| source | T | dH_Si | dH_Mg | Ef(Mg2Si)/f.u. | dE_rxn | B(Mg-Si) | dH_Si paired | x_Si(723 K) |",
         "|---|---|---|---|---|---|---|---|---|"]
    f = lambda v: "" if v is None else f"{v:+.3f}"
    for k, r in rows.items():
        s, T = k.split("|")
        L.append(f"| {LABEL[s]} | {T} | {f(r.get('dH_Si'))} | {f(r.get('dH_Mg'))} | {f(r.get('Ef_Mg2Si_fu'))} | "
                 f"{f(r.get('dE_rxn'))} | {f(r.get('B_MgSi_1nn'))} | {f(r.get('dH_Si_paired'))} | "
                 f"{'' if r.get('x_Si_723K') is None else f'{r['x_Si_723K']:.1e}'} |")
    L += ["", "Mg2Si band gap per source (eV): " + ", ".join(f"{k} {v:.2f}" for k, v in gaps.items())]
    (RES / "analysis" / "interface_study.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
