#!/usr/bin/env python
"""Archived r2SCAN formation energies from the Materials Project for the pool.

For every pool compound that has a Materials Project id (candidates_index.json),
fetch the formation energy per atom from the MP thermo endpoint for the
thermo types r2SCAN (pure r2SCAN, VASP PAW) and GGA_GGA+U. MP now serves ids
in a base-26 letter form (mp-134 -> mp-aaaaaafe); both forms are recorded.
Needs MP_API in the environment. Output: results/mp_r2scan/mp_r2scan_formation.json
"""
import datetime
import json
import os
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent


def new_id(mpid: str) -> str:
    n = int(mpid.split("-")[1])
    s = ""
    for _ in range(8):
        n, r = divmod(n, 26)
        s = chr(97 + r) + s
    return "mp-" + s


def main():
    idx = json.loads((HERE / "candidates_index.json").read_text())
    pool = [d for d in idx if d.get("mpid")]
    ids = {new_id(d["mpid"]): d for d in pool}
    q = urllib.parse.urlencode({"material_ids": ",".join(ids), "_limit": 1000,
                                "_fields": "material_id,thermo_type,formation_energy_per_atom,energy_above_hull"})
    req = urllib.request.Request("https://api.materialsproject.org/materials/thermo/?" + q,
                                 headers={"X-API-KEY": os.environ["MP_API"], "User-Agent": "singular_obj"})
    data = json.load(urllib.request.urlopen(req, timeout=120))["data"]
    recs = {d["cid"]: {"mpid": d["mpid"], "mpid_new": new_id(d["mpid"]), "formula": d["formula"]} for d in pool}
    for row in data:
        d = ids.get(row["material_id"])
        if d is None:
            continue
        tag = {"r2SCAN": "r2scan", "GGA_GGA+U": "gga"}.get(row["thermo_type"])
        if tag:
            recs[d["cid"]][f"{tag}_ef"] = row["formation_energy_per_atom"]
            recs[d["cid"]][f"{tag}_ehull"] = row["energy_above_hull"]
    out = {"source": "Materials Project thermo endpoint, thermo_type r2SCAN (VASP PAW r2SCAN) and GGA_GGA+U",
           "accessed": datetime.date.today().isoformat(), "records": recs}
    (HERE / "results" / "mp_r2scan").mkdir(parents=True, exist_ok=True)
    (HERE / "results" / "mp_r2scan" / "mp_r2scan_formation.json").write_text(json.dumps(out, indent=1))
    n = sum("r2scan_ef" in r for r in recs.values())
    print(f"{len(recs)} pool compounds with an MP id; {n} have an r2SCAN formation energy")


if __name__ == "__main__":
    main()
