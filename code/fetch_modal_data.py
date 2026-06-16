import csv
import io
import json
import os
import re
import sys
import time
import zipfile

import requests

CODE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(CODE_DIR, "data")
os.makedirs(DATA, exist_ok=True)
MODEL_CACHE = os.path.join(CODE_DIR, "iML1515.json")

UA = {"User-Agent": "Mozilla/5.0 (X11; Linux x86_64) Chrome/120 Safari/537.36"}
SESS = requests.Session()
SESS.headers.update(UA)


def get(url, timeout=90, tries=3, **kw):
    last = None
    for i in range(tries):
        try:
            return SESS.get(url, timeout=timeout, **kw)
        except Exception as e:
            last = e
            print(f"    retry {i+1}/{tries} {type(e).__name__}", flush=True)
            time.sleep(2)
    print(f"    FAILED {url}: {last}", flush=True)
    return None


# MediaDB component name -> BiGG exchange id(s) (standard M9).
# One salt supplies multiple ions (e.g. MgSO4 -> Mg2+ & SO4; KH2PO4 -> K+ & Pi).
# Substring match: ions from all matching keys are unioned.
COMPOUND_TO_EX = {
    "d-glucose": ["EX_glc__D_e"], "glucose": ["EX_glc__D_e"],
    "ammonium chloride": ["EX_nh4_e", "EX_cl_e"], "ammonium": ["EX_nh4_e"],
    "phosphate": ["EX_pi_e"], "phosphoric acid": ["EX_pi_e"],
    "potassium phosphate": ["EX_k_e", "EX_pi_e"],
    "monopotassium phosphate": ["EX_k_e", "EX_pi_e"],
    "dipotassium phosphate": ["EX_k_e", "EX_pi_e"],
    "sodium phosphate": ["EX_na1_e", "EX_pi_e"],
    "disodium phosphate": ["EX_na1_e", "EX_pi_e"],
    "magnesium sulfate": ["EX_mg2_e", "EX_so4_e"],
    "magnesium sulphate": ["EX_mg2_e", "EX_so4_e"],
    "magnesium chloride": ["EX_mg2_e", "EX_cl_e"], "magnesium": ["EX_mg2_e"],
    "sulfate": ["EX_so4_e"], "sulphate": ["EX_so4_e"],
    "calcium chloride": ["EX_ca2_e", "EX_cl_e"], "calcium": ["EX_ca2_e"],
    "sodium chloride": ["EX_na1_e", "EX_cl_e"],
    "potassium chloride": ["EX_k_e", "EX_cl_e"], "potassium": ["EX_k_e"],
    "cobalt chloride": ["EX_cobalt2_e", "EX_cl_e"], "cobalt": ["EX_cobalt2_e"],
    "cupric chloride": ["EX_cu2_e", "EX_cl_e"], "cupric": ["EX_cu2_e"],
    "copper": ["EX_cu2_e"],
    "zinc sulfate": ["EX_zn2_e", "EX_so4_e"], "zinc chloride": ["EX_zn2_e", "EX_cl_e"],
    "zinc": ["EX_zn2_e"],
    "ferrous sulfate": ["EX_fe2_e", "EX_so4_e"], "ferrous": ["EX_fe2_e"],
    "ferric chloride": ["EX_fe3_e", "EX_cl_e"], "ferric": ["EX_fe3_e"],
    "iron(ii)": ["EX_fe2_e"], "iron(iii)": ["EX_fe3_e"], "iron": ["EX_fe3_e"],
    "manganese chloride": ["EX_mn2_e", "EX_cl_e"], "manganese": ["EX_mn2_e"],
    "nickel chloride": ["EX_ni2_e", "EX_cl_e"], "nickel": ["EX_ni2_e"],
    "sodium molybdate": ["EX_na1_e", "EX_mobd_e"], "molybdate": ["EX_mobd_e"],
    "sodium selenite": ["EX_na1_e", "EX_sel_e"], "selenite": ["EX_sel_e"],
    "sodium": ["EX_na1_e"], "chloride": ["EX_cl_e"],
    "thiamine": ["EX_thm_e"], "thiamine hydrochloride": ["EX_thm_e"],
    "boric acid": [], "borate": [], "edta": [],  # no iML1515 exchange
}


def map_compound(name):
    """Return list of BiGG EX ids for a compound (union over matching substrings)."""
    n = name.strip().lower()
    if n in COMPOUND_TO_EX:
        return list(COMPOUND_TO_EX[n])
    hits = []
    for key, exs in COMPOUND_TO_EX.items():
        if key in n:
            hits.extend(exs)
    return list(dict.fromkeys(hits)) if hits else "UNMAPPED"


def fetch_mediadb(media_id=147):
    print(f"[MediaDB] media/{media_id} (M9 MG1655)...", flush=True)
    from bs4 import BeautifulSoup
    r = get(f"https://mediadb.systemsbiology.net/defined_media/media/{media_id}/")
    if not r:
        return None
    s = BeautifulSoup(r.text, "lxml")
    comps = []
    for t in s.find_all("table"):
        hdr = [x.get_text(strip=True) for x in t.find_all("th")]
        if not any("compound" in h.lower() for h in hdr):
            continue
        for row in t.find_all("tr")[1:]:
            cells = [c.get_text(strip=True) for c in row.find_all("td")]
            if len(cells) >= 2:
                comps.append((cells[0], cells[1]))
    print(f"    {len(comps)} components", flush=True)

    medium, provenance = {}, {}
    unmapped = []
    for name, amount in comps:
        exs = map_compound(name)
        try:
            mM = float(re.sub(r"[^\d.eE+-]", "", amount))
        except ValueError:
            mM = None
        if exs == "UNMAPPED":
            unmapped.append(name)
            continue
        for ex in exs:  # one salt -> multiple ions (e.g. MgSO4 -> Mg2+ & SO4)
            uptake = 10.0 if ex == "EX_glc__D_e" else 1000.0
            medium[ex] = uptake
            provenance.setdefault(ex, []).append(
                {"compound": name, "amount_mM": mM, "uptake_set": uptake})
    # implicit aerobic / proton / water
    for ex in ["EX_o2_e", "EX_h2o_e", "EX_h_e", "EX_co2_e"]:
        if ex not in medium:
            medium[ex] = 1000.0
            provenance.setdefault(ex, []).append({"compound": "implicit aerobic", "uptake_set": 1000.0})
    # Backfill essential ions: ensure every ion of the iML1515-default M9 is open
    # (MediaDB sometimes omits trace/cation entries; without K/Mg the model is
    # infeasible). Backfilled ions are flagged in provenance.
    from modal_context import M9_GLUCOSE_AEROBIC
    backfilled = []
    for ex, up in M9_GLUCOSE_AEROBIC.items():
        if ex not in medium:
            medium[ex] = up
            backfilled.append(ex)
            provenance.setdefault(ex, []).append(
                {"compound": "BACKFILL iML1515-default M9 (tak ada di MediaDB list)", "uptake_set": up})

    out = {"_meta": {"source": "MediaDB media/147 (M9; Fischer et al), MG1655",
                     "doi": "10.1371/journal.pone.0103548",
                     "note": "MediaDB composition (multi-ion per salt); canonical uptake "
                             "bounds (glc=10, others=1000); mM amount in provenance; "
                             "essential ions absent from MediaDB backfilled from iML1515-default M9",
                     "unmapped_compounds": unmapped,
                     "backfilled_ions": backfilled},
           "medium": medium, "provenance": provenance}
    path = os.path.join(DATA, "medium_m9.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"    -> {path}  ({len(medium)} exchange; unmapped={unmapped})", flush=True)
    return out


BNIDS = {
    "110422": "ATP non-growth maintenance (aerobic, glucose)",
    "111285": "NGAM",
    "110421": "GAM",
    "111006": "ATP concentration intracellular",
    "104673": "ATP concentration, glucose-fed exponential",
    "101983": "ATP requirement growth on glucose",
}


def fetch_bionumbers():
    print("[BioNumbers] BNID NGAM/GAM/[ATP]...", flush=True)
    import html as htmlmod
    rows = []
    for bnid, desc in BNIDS.items():
        r = get(f"https://bionumbers.hms.harvard.edu/bionumber.aspx?id={bnid}", tries=3, timeout=60)
        val, unit = "", ""
        if r:
            txt = htmlmod.unescape(re.sub(r"<[^>]+>", " ", r.text))
            txt = re.sub(r"\s+", " ", txt)
            m = re.search(r"Value\s+([\d.,]+(?:[eE][+-]?\d+)?)\s*([A-Za-z][\w/%.°]*)", txt)
            if m:
                val, unit = m.group(1), m.group(2)
            else:
                m2 = re.search(r"([\d.]+)\s*(mmol\s*ATP[/\w]*)", txt)
                if m2:
                    val, unit = m2.group(1), m2.group(2)
        rows.append({"bnid": bnid, "quantity": desc, "value": val, "unit": unit})
        print(f"    BNID {bnid}: {val} {unit}  ({desc})", flush=True)
    path = os.path.join(DATA, "bionumbers_atpm.csv")
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["bnid", "quantity", "value", "unit"])
        w.writeheader()
        w.writerows(rows)
    print(f"    -> {path}", flush=True)
    return rows


def build_kegg_to_bigg():
    """Map KEGG compound id -> BiGG metabolite base id from iML1515 annotation."""
    import cobra
    m = cobra.io.load_json_model(MODEL_CACHE)
    k2b = {}
    for met in m.metabolites:
        base = met.id.rsplit("_", 1)[0]  # drop compartment suffix
        kegg = met.annotation.get("kegg.compound")
        if kegg:
            for kid in (kegg if isinstance(kegg, list) else [kegg]):
                k2b.setdefault(kid, base)
    return k2b


def fetch_ecmdb(max_pages=60):
    print("[ECMDB] concentrations (best-effort)...", flush=True)
    from bs4 import BeautifulSoup
    # 1) map met_id (ECMDB) -> kegg from bulk zip
    zpath = os.path.join(DATA, "ecmdb.json.zip")
    metid_to_kegg = {}
    try:
        if not os.path.exists(zpath):
            r = get("https://ecmdb.ca/download/ecmdb.json.zip", timeout=120)
            if r:
                open(zpath, "wb").write(r.content)
        if os.path.exists(zpath):
            zf = zipfile.ZipFile(zpath)
            data = json.loads(zf.read(zf.namelist()[0]).decode("utf-8", "replace"))
            for rec in (data if isinstance(data, list) else data.values()):
                mid = rec.get("met_id") or rec.get("m2m_id") or rec.get("id")
                if mid and rec.get("kegg_id"):
                    metid_to_kegg[str(mid)] = rec["kegg_id"]
            print(f"    bulk map: {len(metid_to_kegg)} met_id->kegg", flush=True)
    except Exception as e:
        print(f"    bulk zip failed: {type(e).__name__}: {e}", flush=True)

    k2b = build_kegg_to_bigg()
    print(f"    iML1515 kegg->bigg: {len(k2b)}", flush=True)

    # 2) harvest concentration table via HTML pagination
    rows, seen = [], set()
    for page in range(1, max_pages + 1):
        r = get(f"https://ecmdb.ca/concentrations?page={page}", timeout=90, tries=2)
        if not r or r.status_code >= 500:
            print(f"    page {page}: {getattr(r,'status_code','ERR')} -> stop", flush=True)
            break
        s = BeautifulSoup(r.text, "lxml")
        t = s.find("table")
        if not t:
            break
        trs = t.find_all("tr")[1:]
        new = 0
        for tr in trs:
            c = [x.get_text(strip=True) for x in tr.find_all("td")]
            if len(c) < 3:
                continue
            ecmdb_id, name, conc = c[0], c[1], c[2]
            cond = c[3] if len(c) > 3 else ""
            cite = c[-1] if len(c) > 4 else ""
            key = (ecmdb_id, conc, cond)
            if key in seen:
                continue
            seen.add(key)
            new += 1
            # convert conc "1± 0 uM" -> mM
            mu = re.search(r"([\d.]+).*?(uM|mM|M|nM)", conc)
            mM = None
            if mu:
                v = float(mu.group(1)); u = mu.group(2)
                mM = v / 1000 if u == "uM" else v if u == "mM" else v * 1000 if u == "M" else v / 1e6
            kegg = metid_to_kegg.get(ecmdb_id)
            bigg = k2b.get(kegg) if kegg else None
            rows.append({"ecmdb_id": ecmdb_id, "name": name, "kegg_id": kegg or "",
                         "bigg_id": bigg or "", "conc_mM": mM if mM is not None else "",
                         "raw_conc": conc, "condition": cond, "ref": cite})
        print(f"    page {page}: +{new} rows (total {len(rows)})", flush=True)
        if new == 0:
            break
        time.sleep(1)

    path = os.path.join(DATA, "ecmdb_conc.csv")
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["ecmdb_id", "name", "kegg_id", "bigg_id",
                                          "conc_mM", "raw_conc", "condition", "ref"])
        w.writeheader()
        w.writerows(rows)
    mapped = sum(1 for r in rows if r["bigg_id"])
    print(f"    -> {path}  ({len(rows)} concentrations, {mapped} mapped to BiGG id)", flush=True)
    return rows


def main():
    which = sys.argv[1:] or ["mediadb", "bionumbers", "ecmdb"]
    summary = {}
    if "mediadb" in which:
        try:
            r = fetch_mediadb(); summary["mediadb"] = len(r["medium"]) if r else "FAIL"
        except Exception as e:
            summary["mediadb"] = f"ERR {type(e).__name__}: {e}"
    if "bionumbers" in which:
        try:
            r = fetch_bionumbers(); summary["bionumbers"] = len(r)
        except Exception as e:
            summary["bionumbers"] = f"ERR {type(e).__name__}: {e}"
    if "ecmdb" in which:
        try:
            r = fetch_ecmdb(); summary["ecmdb"] = len(r)
        except Exception as e:
            summary["ecmdb"] = f"ERR {type(e).__name__}: {e}"
    print("\n==== SUMMARY ====", flush=True)
    for k, v in summary.items():
        print(f"  {k}: {v}", flush=True)


if __name__ == "__main__":
    main()
