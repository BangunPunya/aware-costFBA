from __future__ import annotations

import json
import os
import math

CODE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(CODE_DIR, "data")
CACHE_PATH = os.path.join(DATA_DIR, "dG_cache.json")

# Ambang screen kelayakan (kJ/mol).
#
# Reaksi pusat-karbon (PGK +19.5, FBA +23, MDH +26.5 kJ/mol) punya ΔG'° positif pada
# keadaan standar (1 M) maupun fisiologis (1 mM), namun berjalan maju in vivo karena
# rasio konsentrasi nyata (mis. ATP/ADP tinggi, produk rendah) menggeser ΔG' aktual
# < 0. Screen tanda-ΔG naif akan keliru memblokir reaksi ini -> glikolisis mati ->
# model infeasible (false-positive klasik).
#
# Screen memakai ΔG' = ΔG'° + RT·ln(Q), Q = Π [produk]^koef / Π [substrat]^|koef|
# (mM -> M; H2O & H+ aktivitas = 1, dikecualikan). Konsentrasi dari
# ModalContext.metabolite_conc (ECMDB), default DEFAULT_CONC_M (1 mM) bila tak
# terukur. Arah-diizinkan di-flag infeasible hanya bila ΔG'(arah itu) jelas > 0
# melampaui ketidakpastian ΔG'° (sd) + marjin:
#     fwd infeasible <=>  (ΔG'°_fwd + RT·ln Q) - margin_sd*sd > MIN_MARGIN_KJ
# Bukan TFA LP penuh (tak ada variabel ΔG di-kopel ke fluks), tapi memakai
# konsentrasi dan mampu mem-flag arah infeasible.
FEASIBILITY_MARGIN_SD = 1.0     # berapa sd ΔG'° di atas ΔG' sebelum dibilang "jelas" infeasible
# MIN_MARGIN_KJ menyerap ketidakpastian konsentrasi metabolomik ECMDB (titik-tunggal,
# ~1 orde besar per metabolit). Faktor-10 per spesies = RT·ln(10) ≈ 5.9 kJ/mol.
# Reaksi pusat-karbon (PGK/GAPD/MDH/TPI) tampil ΔG' maju ~+10-19 kJ/mol di ECMDB
# titik-tunggal namun JELAS berjalan maju in vivo (false-positive klasik). Maka hanya
# flag bila ΔG' maju > MIN_MARGIN_KJ (≈ 2× RT·ln10) di atas sd - cukup konservatif agar
# iML1515 kurasi tetap tumbuh, tapi tetap mem-flag arah yg ΔG'-nya benar-benar
# besar-positif (propagasi ketidakpastian ΔG'° -> ΔG_margin). Bisa di-override per-panggilan.
MIN_MARGIN_KJ = 12.0            # marjin mutlak (kJ/mol) ≈ ketidakpastian konsentrasi metabolomik
RT_KJ = 2.5776                  # R*T pada 310.15 K (37 °C) (kJ/mol) - fisiologis E. coli
DEFAULT_CONC_M = 1e-3           # default konsentrasi metabolit tak-terukur = 1 mM (Molar)

# Sufiks kompartemen iML1515 yang umum (untuk strip id metabolit cobra -> base bigg).
_COMPARTMENTS = ("_c", "_p", "_e")

# Singleton ComponentContribution (init pertama mengunduh cache ~1.3 GB; berikutnya beberapa detik).
_CC = None


def _get_cc():
    """Lazy-load ComponentContribution. Hanya dipanggil saat butuh hitung (bukan saat baca cache)."""
    global _CC
    if _CC is None:
        from equilibrator_api import ComponentContribution
        _CC = ComponentContribution()
    return _CC


def _load_cache() -> dict:
    if os.path.exists(CACHE_PATH):
        try:
            with open(CACHE_PATH) as f:
                return json.load(f)
        except (json.JSONDecodeError, OSError):
            return {}
    return {}


def _save_cache(cache: dict) -> None:
    os.makedirs(DATA_DIR, exist_ok=True)
    with open(CACHE_PATH, "w") as f:
        json.dump(cache, f, indent=2, sort_keys=True)


def _bigg_accession(met) -> str | None:
    """Map metabolit cobra -> accession eQuilibrator `bigg.metabolite:<id>` (dari anotasi bigg.metabolite, atau strip sufiks kompartemen dari met.id)."""
    ann = met.annotation.get("bigg.metabolite") if met.annotation else None
    if isinstance(ann, (list, tuple)) and ann:
        ann = ann[0]
    if isinstance(ann, str) and ann:
        return f"bigg.metabolite:{ann}"
    # fallback: strip sufiks kompartemen dari id cobra
    mid = met.id
    for suf in _COMPARTMENTS:
        if mid.endswith(suf):
            return f"bigg.metabolite:{mid[: -len(suf)]}"
    return f"bigg.metabolite:{mid}"


def _build_formula(rxn) -> str | None:
    """Bangun formula reaksi eQuilibrator (`<koef> <acc> + ... = ...`) dari metabolit cobra; None bila ada metabolit tak-terpetakan atau reaksi tanpa metabolit."""
    left, right = [], []
    for met, coef in rxn.metabolites.items():
        acc = _bigg_accession(met)
        if acc is None:
            return None
        c = abs(coef)
        term = f"{c} {acc}" if c != 1 else acc
        (left if coef < 0 else right).append(term)
    if not left or not right:
        return None  # exchange/demand/sink/biomass - bukan reaksi internal seimbang
    return " + ".join(left) + " = " + " + ".join(right)


def _is_skippable(rxn) -> tuple[bool, str]:
    """Reaksi yang eQuilibrator tak bisa/tak relevan resolve: exchange, demand, sink,
    biomass, dan transport antar-kompartemen (senyawa sama dua sisi -> ΔG'° transport
    tak ditangani screen sederhana ini)."""
    rid = rxn.id
    if rid.startswith(("EX_", "DM_", "SK_")):
        return True, "exchange/demand/sink"
    if "biomass" in rid.lower() or "BIOMASS" in rid:
        return True, "biomass"
    if len(rxn.metabolites) < 2:
        return True, "boundary/degenerate"
    # Transmembrane proton translocation (ATPS4rpp, CYTBO3_4pp, NADH16pp): H+ di >1
    # kompartemen -> kelayakannya diatur PROTON-MOTIVE FORCE (Δp = ΔpH + ΔΨ membran),
    # BUKAN ΔG'° kimiawi. eQuilibrator ΔG'° tak memuat potensial-membran -> screen
    # kimiawi akan KELIRU mem-flag arah fisiologis (ATP synthase iML1515 maju ADP->ATP
    # didorong Δp). Maka di-skip dari screen ini (Δp di luar lingkup).
    h_comps = {met.compartment for met in rxn.metabolites if _met_base(met.id) == "h"}
    if len(h_comps) > 1:
        return True, "translokasi-proton transmembran (Δp, di luar ΔG'° kimiawi)"
    # transport: himpunan base-id (tanpa kompartemen) identik kiri==kanan -> reaksi transport
    bases = {}
    for met, coef in rxn.metabolites.items():
        base = met.id
        for suf in _COMPARTMENTS:
            if base.endswith(suf):
                base = base[: -len(suf)]
                break
        bases.setdefault(base, 0.0)
        bases[base] += coef
    # bila semua base bersih (net) ~0 -> murni transport (pindah kompartemen)
    if all(abs(v) < 1e-9 for v in bases.values()):
        return True, "transport (antar-kompartemen)"
    return False, ""


def _compute_one(rxn) -> tuple[tuple[float, float] | None, str]:
    """Hitung ΔG'° (mean, sd) kJ/mol untuk satu reaksi. Return ((mean, sd), note) atau (None, alasan)."""
    skip, why = _is_skippable(rxn)
    if skip:
        return None, f"skip: {why}"
    formula = _build_formula(rxn)
    if formula is None:
        return None, "skip: metabolit tak-terpetakan / tak seimbang"
    cc = _get_cc()
    try:
        eq_rxn = cc.parse_reaction_formula(formula)
    except Exception as e:  # parse gagal (id tak dikenal eQuilibrator)
        return None, f"parse gagal: {type(e).__name__}"
    if not eq_rxn.is_balanced():
        # reaksi tak seimbang (mis. tanpa H/charge) -> ΔG'° tak andal; flag tapi tetap coba? -> skip
        return None, "skip: reaksi tak seimbang (atom/charge)"
    try:
        dg = cc.standard_dg_prime(eq_rxn)
    except Exception as e:
        return None, f"standard_dg_prime gagal: {type(e).__name__}"
    mean = float(dg.value.m_as("kJ/mol"))
    sd = float(dg.error.m_as("kJ/mol"))
    return (mean, sd), "ok"


def get_dG0(model, rxn_ids, use_cache: bool = True, verbose: bool = False
            ) -> dict[str, tuple[float, float]]:
    """ΔG'° (mean, sd) kJ/mol per reaksi untuk daftar rxn_ids; di-cache ke code/data/dG_cache.json. Reaksi tak-resolvable (transport/exchange/biomass/parse-gagal) di-skip dengan catatan agar tak dihitung ulang."""
    cache = _load_cache() if use_cache else {}
    result: dict[str, tuple[float, float]] = {}
    dirty = False

    for rid in rxn_ids:
        if use_cache and rid in cache:
            entry = cache[rid]
            val = entry.get("dG0") if isinstance(entry, dict) else entry
            if val is not None:
                result[rid] = (float(val[0]), float(val[1]))
            continue
        try:
            rxn = model.reactions.get_by_id(rid)
        except KeyError:
            cache[rid] = {"dG0": None, "note": "rxn tak ada di model"}
            dirty = True
            continue
        val, note = _compute_one(rxn)
        cache[rid] = {"dG0": list(val) if val else None, "note": note}
        dirty = True
        if val is not None:
            result[rid] = val
        if verbose:
            print(f"  {rid}: {note}" + (f"  ΔG'°={val[0]:.1f}±{val[1]:.1f}" if val else ""))

    if dirty and use_cache:
        _save_cache(cache)
    return result


def _met_base(met_id: str) -> str:
    """Strip sufiks kompartemen -> base id (mis. `atp_c` -> `atp`). Untuk deteksi H+/H2O."""
    for suf in _COMPARTMENTS:
        if met_id.endswith(suf):
            return met_id[: -len(suf)]
    return met_id


def _ln_Q(rxn, conc: dict[str, float] | None) -> float:
    """ln(Q) arah maju, Q = Π a_i^ν_i (ν>0 produk, ν<0 substrat; a_i = konsentrasi MOLAR; H2O & H+ aktivitas=1 dikecualikan; metabolit tak-terukur default DEFAULT_CONC_M = 1 mM)."""
    conc = conc or {}
    ln_q = 0.0
    for met, coef in rxn.metabolites.items():
        base = _met_base(met.id)
        if base in ("h", "h2o"):
            continue  # aktivitas = 1 (konvensi standar) - jangan masukkan ke Q
        mM = conc.get(met.id)
        a = (mM * 1e-3) if (mM is not None and mM > 0) else DEFAULT_CONC_M  # Molar
        ln_q += coef * math.log(a)
    return ln_q


def flag_infeasible(model, dG0: dict[str, tuple[float, float]],
                    conc: dict[str, float] | None = None,
                    margin_sd: float = FEASIBILITY_MARGIN_SD,
                    min_margin_kj: float = MIN_MARGIN_KJ) -> dict[str, str]:
    """Screen kelayakan sadar-konsentrasi: ΔG' = ΔG'° + RT·ln(Q). Flag arah-diizinkan
    bila ΔG'(arah itu) jelas > 0 melampaui sd + marjin mutlak:
        ΔG'_fwd = mean + RT·ln(Q);  ΔG'_rev = -ΔG'_fwd
        fwd infeasible <=> ΔG'_fwd - margin_sd*sd > min_margin_kj
    conc = {met_id_c: mM} dari ECMDB (default 1 mM; H2O & H+ aktivitas=1). Memakai
    konsentrasi nyata (mass-action) lewat RT·ln(Q) sehingga PGK/FBA/MDH (ΔG'° positif)
    yg berjalan maju in vivo tidak keliru di-flag. Bukan TFA LP penuh.
    Return {rxn_id: arah}: "forward" -> upper_bound=0, "reverse" -> lower_bound=0.
    Hanya reaksi yg arah-infeasible-nya MASIH diizinkan bound.
    """
    flagged: dict[str, str] = {}
    for rid, (mean, sd) in dG0.items():
        try:
            rxn = model.reactions.get_by_id(rid)
        except KeyError:
            continue
        # Lewati reaksi yg screen kimiawi tak boleh nilai (translokasi-proton transmembran,
        # transport, dll) - feasibilitasnya diatur Δp/transport, bukan ΔG'° kimiawi.
        skip, _why = _is_skippable(rxn)
        if skip:
            continue
        sd = abs(sd) if not math.isnan(sd) else 0.0
        dG_fwd = mean + RT_KJ * _ln_Q(rxn, conc)  # ΔG' arah maju pada konsentrasi diberikan
        fwd_allowed = rxn.upper_bound > 1e-9
        rev_allowed = rxn.lower_bound < -1e-9

        fwd_infeasible = (dG_fwd - margin_sd * sd) > min_margin_kj
        rev_infeasible = (-dG_fwd - margin_sd * sd) > min_margin_kj

        if fwd_infeasible and fwd_allowed:
            flagged[rid] = "forward"
        elif rev_infeasible and rev_allowed:
            flagged[rid] = "reverse"
    return flagged


if __name__ == "__main__":
    # smoke test ringan: butuh cache eQuilibrator sudah terbangun.
    from cobra.io import load_json_model
    m = load_json_model(os.path.join(CODE_DIR, "iML1515.json"))
    test_ids = ["PGK", "PGI", "FBA", "TPI", "PGM", "ENO", "PYK", "ATPM"]
    print("Menghitung ΔG'° untuk:", test_ids)
    dg = get_dG0(m, test_ids, verbose=True)
    print("\nHasil:")
    for rid, (mean, sd) in dg.items():
        print(f"  {rid}: ΔG'° = {mean:+.1f} ± {sd:.1f} kJ/mol")
    print("\nflag_infeasible:", flag_infeasible(m, dg))
