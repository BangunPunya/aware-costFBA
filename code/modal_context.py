from __future__ import annotations

import csv
import json
import os
from dataclasses import dataclass, field

import cobra
from cobra.io import load_json_model

CODE_DIR = os.path.dirname(os.path.abspath(__file__))
MODEL_CACHE = os.path.join(CODE_DIR, "iML1515.json")
DATA_DIR = os.path.join(CODE_DIR, "data")
NGAM_BNID = "110422"  # ATP non-growth maintenance (aerobic, glucose) -> ATPM

# M9 aerobic glucose-minimal = iML1515 default medium.
# Uptake in mmol gDW^-1 h^-1 (lower_bound = -value). Glucose limiting (-10);
# salts/minerals & O2 unlimited (-1000) per the BiGG curation.
M9_GLUCOSE_AEROBIC = {
    "EX_glc__D_e": 10.0, "EX_o2_e": 1000.0, "EX_nh4_e": 1000.0,
    "EX_pi_e": 1000.0, "EX_so4_e": 1000.0, "EX_h2o_e": 1000.0,
    "EX_h_e": 1000.0, "EX_co2_e": 1000.0, "EX_k_e": 1000.0,
    "EX_na1_e": 1000.0, "EX_cl_e": 1000.0, "EX_ca2_e": 1000.0,
    "EX_mg2_e": 1000.0, "EX_mn2_e": 1000.0, "EX_fe2_e": 1000.0,
    "EX_fe3_e": 1000.0, "EX_zn2_e": 1000.0, "EX_cu2_e": 1000.0,
    "EX_cobalt2_e": 1000.0, "EX_ni2_e": 1000.0, "EX_mobd_e": 1000.0,
    "EX_sel_e": 1000.0, "EX_slnt_e": 1000.0, "EX_tungs_e": 1000.0,
}

# ATPM non-growth maintenance (mmol ATP gDW^-1 h^-1).
# iML1515 default = 6.86 (curated). Override from BioNumbers NGAM
# (BNID 110422/111285) when available.
ATPM_DEFAULT = 6.86


@dataclass
class ModalContext:
    base: cobra.Model
    medium: dict[str, float] = field(default_factory=lambda: dict(M9_GLUCOSE_AEROBIC))
    atpm: float = ATPM_DEFAULT
    # Pending-data layers - attached when data arrives (web fetch / equilibrator)
    metabolite_conc: dict[str, float] | None = None   # ECMDB
    dG0: dict[str, tuple[float, float]] | None = None  # eQuilibrator
    kcat: dict | None = None                            # BRENDA (optional)
    toggles: dict[str, bool] = field(default_factory=lambda: {
        "medium": True, "atpm": True,
        "ecmdb_pool": False, "tfa_dG": False, "ec_kcat": False,
    })
    provenance: dict = field(default_factory=dict)

    @classmethod
    def load(cls, model_path: str = MODEL_CACHE) -> "ModalContext":
        if not os.path.exists(model_path):
            raise FileNotFoundError(
                f"{model_path} not found - run m0_baseline.py first (download+cache iML1515)."
            )
        base = load_json_model(model_path)
        mc = cls(base=base)
        mc.provenance = {
            "model": "iML1515 (BiGG, local cache)",
            "medium": "M9 aerobic glucose = iML1515 default",
            "atpm": f"curated default {ATPM_DEFAULT}; BioNumbers NGAM override = pending",
            "ecmdb_pool": "pending data (ECMDB dump)",
            "tfa_dG": "pending data (equilibrator-api cache 1.3GB)",
        }
        mc._attach_data_files()   # auto-load data/ produced by fetch_modal_data.py if present
        return mc

    def _attach_data_files(self, data_dir: str = DATA_DIR) -> None:
        """Auto-load modal data from data/ (produced by fetch_modal_data.py); absent files fall back to defaults (constant M9 + curated ATPM). Pending toggles not switched on automatically - only data prepared so the layer can be enabled for ablation."""
        # medium from MediaDB M9
        p = os.path.join(data_dir, "medium_m9.json")
        if os.path.exists(p):
            d = json.load(open(p))
            if d.get("medium"):
                self.medium = {k: float(v) for k, v in d["medium"].items()}
                self.provenance["medium"] = d.get("_meta", {}).get("source", p)
        # ATPM from BioNumbers NGAM
        p = os.path.join(data_dir, "bionumbers_atpm.csv")
        if os.path.exists(p):
            for row in csv.DictReader(open(p)):
                if row.get("bnid") == NGAM_BNID and row.get("value"):
                    try:
                        self.atpm = float(row["value"]); self.provenance["atpm"] = \
                            f"BioNumbers BNID {NGAM_BNID} NGAM = {self.atpm}"
                    except ValueError:
                        pass
        # metabolite concentrations from ECMDB -> metabolite_conc {met_id_c: mM}
        p = os.path.join(data_dir, "ecmdb_conc.csv")
        if os.path.exists(p):
            conc = {}
            for row in csv.DictReader(open(p)):
                bigg, mM = row.get("bigg_id"), row.get("conc_mM")
                if bigg and mM not in (None, ""):
                    try:
                        conc[f"{bigg}_c"] = float(mM)
                    except ValueError:
                        pass
            if conc:
                self.metabolite_conc = conc
                self.provenance["ecmdb_pool"] = f"ECMDB {len(conc)} concentrations mapped to BiGG"
        # ΔG'° cache from equilibrator -> dG0 {rxn_id: (mean, sd)}
        p = os.path.join(data_dir, "dG_cache.json")
        if os.path.exists(p):
            d = json.load(open(p))
            dg = {}
            for rid, v in d.items():
                if isinstance(v, dict):       # {rxn:{"dG0":[mean,sd],"note":..}}
                    v = v.get("dG0")
                if isinstance(v, (list, tuple)) and len(v) >= 2:
                    try:
                        dg[rid] = (float(v[0]), float(v[1]))
                    except (ValueError, TypeError):
                        pass
            if dg:
                self.dG0 = dg
                self.provenance["tfa_dG"] = f"equilibrator ΔG'° cache: {len(dg)} reactions"
        return None

    def as_polos(self) -> "ModalContext":
        """Plain-FBA baseline: rich medium + ATPM=0, no TFA. Deliberately permissive to contrast with modal layers."""
        rich = {ex.id: 1000.0 for ex in self.base.exchanges}
        polos = ModalContext(
            base=self.base, medium=rich, atpm=0.0,
            toggles={"medium": True, "atpm": True,
                     "ecmdb_pool": False, "tfa_dG": False, "ec_kcat": False},
        )
        polos.provenance = {"mode": "plain FBA: rich medium (all EX uptake -1000), ATPM=0"}
        return polos

    def apply(self, model: cobra.Model) -> cobra.Model:
        """Apply active modal layers to model (order: medium->ATPM->pool->TFA->ec)."""
        if self.toggles.get("medium"):
            self._apply_medium(model)
        if self.toggles.get("atpm"):
            self._apply_atpm(model)
        if self.toggles.get("ecmdb_pool"):
            self._apply_ecmdb_pool(model)
        if self.toggles.get("tfa_dG"):
            self._apply_tfa(model)
        if self.toggles.get("ec_kcat"):
            self._apply_kcat(model)
        return model

    # live layers
    def _apply_medium(self, model: cobra.Model) -> None:
        # Close all exchange uptake first, then open per the medium.
        for ex in model.exchanges:
            ex.lower_bound = 0.0
        for ex_id, uptake in self.medium.items():
            try:
                model.reactions.get_by_id(ex_id).lower_bound = -abs(uptake)
            except KeyError:
                pass  # exchange not in model - ignore

    def _apply_atpm(self, model: cobra.Model) -> None:
        try:
            model.reactions.ATPM.lower_bound = self.atpm
        except (KeyError, AttributeError):
            pass

    # pending-data layers (toggle present, need external data)
    def _apply_ecmdb_pool(self, model: cobra.Model) -> None:
        """ECMDB metabolite-pool layer: supplies measured metabolite concentrations (self.metabolite_conc, {met_id_c: mM}) for the concentration-aware thermo screen in _apply_tfa. Validates concentrations available (no bound change). Effect appears when tfa_dG is on: _apply_tfa passes self.metabolite_conc to dG.flag_infeasible so ΔG' is at measured concentrations. Toggle: tfa_dG only -> ΔG' at 1 mM default (Q≈1); tfa_dG + ecmdb_pool -> ΔG' at measured ECMDB concentrations."""
        if not self.metabolite_conc:
            raise NotImplementedError(
                "ecmdb_pool on but metabolite_conc empty - fetch ECMDB first."
            )
        # No standalone bound. Record concentration availability for audit;
        # _apply_tfa consumes self.metabolite_conc when ecmdb_pool is on.
        self.provenance["ecmdb_pool_applied"] = {
            "n_conc": len(self.metabolite_conc),
            "role": "feed measured concentrations to TFA, not a standalone pool constraint",
        }

    def _apply_tfa(self, model: cobra.Model) -> None:
        """Concentration-aware thermodynamic screen: ΔG' = ΔG'° + RT·ln(Q). For each reaction with ΔG'°, dG.flag_infeasible computes forward ΔG' at metabolite concentrations and flags directions where ΔG' > 0 (beyond sd + margin); flagged directions tightened to bound 0. Toggles: tfa_dG only -> conc=None -> all mets 1 mM (Q≈1, ~sign of ΔG'°); tfa_dG + ecmdb_pool -> conc=self.metabolite_conc -> ΔG' at measured ECMDB (unmeasured default 1 mM, H2O/H+ activity=1). RT = 2.5776 kJ/mol @ 310.15 K (dG.RT_KJ). Uses concentrations + flags infeasible directions but is NOT a full TFA LP (no ΔG_r variables coupled to flux direction via big-M/indicator). A standard TFA LP (out of scope) would add per-reaction ΔG_r + sign(ΔG_r) <-> flux-direction coupling + ΔG_r = ΔG'° + RT·Σν·ln(c) with bounded LP c, giving simultaneous flux-thermo feasibility."""
        if not self.dG0:
            raise NotImplementedError(
                "tfa_dG on but dG0 empty - pip install equilibrator-api + cache."
            )
        import dG as _dG  # local import: avoid equilibrator dependency when tfa_dG off
        # ecmdb_pool on -> measured concentrations; else None -> 1 mM default in flag_infeasible.
        use_conc = self.metabolite_conc if self.toggles.get("ecmdb_pool") else None
        flagged = _dG.flag_infeasible(model, self.dG0, conc=use_conc)
        n_constrained = 0
        for rid, direction in flagged.items():
            try:
                rxn = model.reactions.get_by_id(rid)
            except KeyError:
                continue
            if direction == "forward" and rxn.upper_bound > 0.0:
                rxn.upper_bound = 0.0  # forbid forward direction (forward ΔG' > 0)
                n_constrained += 1
            elif direction == "reverse" and rxn.lower_bound < 0.0:
                rxn.lower_bound = 0.0  # forbid reverse direction
                n_constrained += 1
        # record to provenance for audit (reactions screened / tightened)
        self.provenance["tfa_dG_screen"] = {
            "n_dG0": len(self.dG0),
            "n_flagged": len(flagged),
            "n_constrained": n_constrained,
            "conc_source": ("ECMDB (ecmdb_pool on)" if use_conc else "default 1 mM (ecmdb_pool off)"),
            "method": "ΔG'+concentration screen (RT·ln Q), not a full TFA LP",
        }

    def _apply_kcat(self, model: cobra.Model) -> None:
        if not self.kcat:
            raise NotImplementedError("ec_kcat on but kcat empty - BRENDA/SABIO-RK.")
        # enzyme-capacity constraint (ec-FBA/GECKO)
