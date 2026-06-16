#!/usr/bin/env python3
from __future__ import annotations

import json
import os
from dataclasses import dataclass, field, asdict
from typing import Optional

import requests

UNIPROT_REST = "https://rest.uniprot.org/uniprotkb/{acc}.json"
_DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


# Fetch / cache
def fetch_uniprot(accession: str = "P08660", cache: bool = True) -> dict:
    """Fetch a UniProtKB entry as JSON, caching to data/<acc>.json."""
    os.makedirs(_DATA_DIR, exist_ok=True)
    path = os.path.join(_DATA_DIR, f"{accession}.json")
    if cache and os.path.exists(path):
        with open(path) as fh:
            return json.load(fh)
    resp = requests.get(UNIPROT_REST.format(acc=accession), timeout=60)
    resp.raise_for_status()
    j = resp.json()
    if cache:
        with open(path, "w") as fh:
            json.dump(j, fh, indent=2)
    return j


# Feature parsing
@dataclass
class MutationProposal:
    position: int                       # 1-based residue position
    wt_aa: str                          # wild-type residue
    candidate_mutations: list[str]      # proposed substitution residues
    rationale: str                      # documented reason for targeting
    region_tag: str                     # which WHERE-gate bucket matched
    priority: float                     # ranking score (higher = better)

    def as_tuple(self):
        """(position, wt_aa, candidate_mutations, rationale) per spec."""
        return (self.position, self.wt_aa, self.candidate_mutations, self.rationale)


# Substitution heuristic tables
_CHARGE = {
    "D": "-", "E": "-",
    "K": "+", "R": "+", "H": "+",
}
_SMALL = set("GASCP")
_LARGE = set("FWYRKHQE")

# Disruptive substitution preferences per WT residue class.
_DISRUPT = {
    # charged -> charge reversal / neutralization
    "D": ["A", "N", "K"],
    "E": ["A", "Q", "K"],
    "K": ["E", "A", "Q"],
    "R": ["E", "A", "Q"],
    "H": ["A", "D", "F"],
    # polar -> hydrophobic / size change
    "S": ["F", "A", "L"],
    "T": ["I", "A", "F"],
    "N": ["D", "A", "L"],
    "Q": ["E", "A", "L"],
    "Y": ["A", "F", "S"],
    # hydrophobic -> charge introduction / size change
    "M": ["I", "K", "T"],
    "L": ["R", "A", "F"],
    "I": ["T", "R", "A"],
    "V": ["D", "A", "F"],
    "F": ["S", "A", "D"],
    "W": ["A", "S", "R"],
    "A": ["D", "V", "F"],
    "G": ["D", "A", "R"],
    "P": ["A", "S", "L"],
    "C": ["S", "A", "D"],
}


def _disrupt_candidates(wt: str) -> list[str]:
    return _DISRUPT.get(wt, ["A", "D", "F"])


def _iter_features(uniprot_json: dict):
    for f in uniprot_json.get("features", []):
        loc = f.get("location", {})
        try:
            s = loc["start"]["value"]
            e = loc["end"]["value"]
        except (KeyError, TypeError):
            continue
        if s is None or e is None:
            continue
        yield f, int(s), int(e)


def parse_regulatory_features(uniprot_json: dict) -> dict:
    """Extract regulatory/allosteric annotation defining the WHERE-gate. Returns dict: lysine_sites (set[int], Binding-site ligand=L-lysine), act_domain ((start,end)|None), interface (list of (start,end) lower-priority), mutagenesis (list of experimentally-tested dicts)."""
    lysine_sites: set[int] = set()
    act_domain = None
    interface: list[tuple[int, int]] = []
    mutagenesis: list[dict] = []

    for f, s, e in _iter_features(uniprot_json):
        ftype = f.get("type")
        desc = (f.get("description") or "")
        if ftype == "Binding site":
            lig = (f.get("ligand") or {}).get("name", "") or ""
            if "lysine" in lig.lower():
                lysine_sites.update(range(s, e + 1))
        elif ftype == "Domain":
            if "ACT" in desc.upper():
                act_domain = (s, e)
        elif ftype == "Region":
            d = desc.lower()
            if "interface" in d or "homodimer" in d or "dimer" in d:
                interface.append((s, e))
        elif ftype == "Mutagenesis":
            alt = f.get("alternativeSequence", {}) or {}
            mutagenesis.append({
                "position": s,
                "end": e,
                "wt": alt.get("originalSequence", ""),
                "mut": alt.get("alternativeSequences", []),
                "description": desc,
            })

    return {
        "lysine_sites": lysine_sites,
        "act_domain": act_domain,
        "interface": interface,
        "mutagenesis": mutagenesis,
    }


# WHERE-gate + proposal generation
# Region priority weights (higher = more likely to host an fbr mutation).
_W_LYS_SITE = 3.0      # annotated L-lysine contact residue
_W_LYS_FLANK = 2.0     # +/- flank window around a lysine contact
_W_ACT = 1.0           # anywhere in the ACT regulatory domain
_W_INTERFACE = 0.3     # dimer interface (weakly regulatory)
_LYS_FLANK = 2         # residues either side of a lysine contact to include


def propose_mutations(
    uniprot_json: dict,
    strategy: str = "lysine_site",
) -> list[MutationProposal]:
    """Propose ranked candidate fbr mutations. strategy: "lysine_site" -> lysine contacts + flank + ACT domain (default); "act_domain" -> whole ACT domain; "all_regulatory" -> ACT + interface/dimer. Returns list[MutationProposal] ranked by priority (desc), each with position, wt residue, disruptive substitutions, rationale, region tag, numeric priority."""
    seq = uniprot_json["sequence"]["value"]
    n = len(seq)
    reg = parse_regulatory_features(uniprot_json)
    lys = reg["lysine_sites"]
    act = reg["act_domain"]
    interface = reg["interface"]

    # Build per-position priority (max over matching buckets) + tag.
    prio: dict[int, float] = {}
    tag: dict[int, str] = {}

    def bump(pos: int, weight: float, label: str):
        if 1 <= pos <= n and weight > prio.get(pos, -1.0):
            prio[pos] = weight
            tag[pos] = label

    if strategy in ("lysine_site", "all_regulatory", "act_domain"):
        # ACT domain baseline
        if act:
            for p in range(act[0], act[1] + 1):
                bump(p, _W_ACT, "ACT_domain")
    if strategy in ("lysine_site", "all_regulatory"):
        for p in lys:
            bump(p, _W_LYS_SITE, "lysine_binding_site")
            for d in range(1, _LYS_FLANK + 1):
                bump(p - d, _W_LYS_FLANK, "lysine_site_flank")
                bump(p + d, _W_LYS_FLANK, "lysine_site_flank")
    if strategy == "all_regulatory":
        for (s, e) in interface:
            for p in range(s, e + 1):
                bump(p, _W_INTERFACE, "interface")

    proposals: list[MutationProposal] = []
    for pos, weight in prio.items():
        wt = seq[pos - 1]
        label = tag[pos]
        cands = _disrupt_candidates(wt)
        rationale = _rationale(label, wt, lys, act)
        proposals.append(MutationProposal(
            position=pos, wt_aa=wt, candidate_mutations=cands,
            rationale=rationale, region_tag=label, priority=weight,
        ))

    # Rank: priority desc, then proximity to nearest lysine contact (closer
    # first), then position for determinism.
    def _near_lys(p: int) -> int:
        return min((abs(p - l) for l in lys), default=10_000)

    proposals.sort(key=lambda m: (-m.priority, _near_lys(m.position), m.position))
    return proposals


def _rationale(label: str, wt: str, lys: set[int], act) -> str:
    if label == "lysine_binding_site":
        return (f"Annotated L-lysine allosteric contact ({wt}); disrupting "
                "effector binding relieves feedback inhibition (fbr).")
    if label == "lysine_site_flank":
        return (f"Flanks an L-lysine allosteric contact ({wt}); local "
                "perturbation can weaken effector binding.")
    if label == "ACT_domain":
        return (f"ACT regulatory domain residue ({wt}); fbr mutations cluster "
                "in the ACT effector-binding fold.")
    if label == "interface":
        return (f"Dimer/interface residue ({wt}); allostery is transmitted "
                "across the interface (weak regulatory signal).")
    return f"Regulatory-region residue ({wt})."


# Optional ESM2 fold-safety gate (NOT the intent-selector)
def esm2_fold_scores(
    seq: str,
    proposals: list[MutationProposal],
    model_name: str = "facebook/esm2_t6_8M_UR50D",
    top_n: Optional[int] = None,
) -> Optional[dict]:
    """Fold-safety score: log P(mt) - log P(wt) per position (more positive = less likely to break fold). Fold gate, not intent gate. Delegated to esm2_gate.score_proposals (dual-env: torch CPU wheels exist for Python 3.10-3.12 only, but this module runs under 3.14; esm2_gate re-execs a conda env (3.11 + torch CPU) as worker, returns None if unreachable). Returns {"{wt}{pos}{mt}": delta_logp} or None (deferred). CPU only. Default wt-marginal scoring."""
    try:
        import esm2_gate  # local module; self-contained dual-env dispatch
    except Exception:  # noqa: BLE001
        return None

    sel = proposals if top_n is None else proposals[:top_n]
    # Flatten MutationProposal -> per-candidate items; esm2_gate accepts the
    # MutationProposal objects directly (it reads .position/.wt_aa/.candidate_mutations).
    return esm2_gate.score_proposals(seq, sel, model_name=model_name)


def fold_gate_proposals(
    seq: str,
    proposals: list["MutationProposal"],
    threshold: Optional[float] = None,
    top_n: Optional[int] = None,
) -> Optional[list[dict]]:
    """Post-filter propose_mutations() through ESM2 fold gate. Returns list of dicts (per wt->mt candidate) with masked/wt-marginal score + fold_pass verdict, or None if DEFERRED. Passes when score >= threshold (default esm2_gate.THRESHOLD_T, permissive guardrail calibrated to reject backbone/motif breakers while passing documented fbr mutants)."""
    try:
        import esm2_gate
    except Exception:  # noqa: BLE001
        return None
    sel = proposals if top_n is None else proposals[:top_n]
    t = esm2_gate.THRESHOLD_T if threshold is None else threshold
    return esm2_gate.fold_gate(seq, sel, threshold=t)


# CLI demo
if __name__ == "__main__":
    j = fetch_uniprot("P08660")
    seq = j["sequence"]["value"]
    reg = parse_regulatory_features(j)
    print(f"P08660  len={len(seq)}  ACT={reg['act_domain']}  "
          f"lysine_sites={sorted(reg['lysine_sites'])}")
    props = propose_mutations(j, strategy="lysine_site")
    print(f"\n{len(props)} candidate positions (top 20):")
    for m in props[:20]:
        print(f"  {m.wt_aa}{m.position:<4} -> {','.join(m.candidate_mutations):<8} "
              f"prio={m.priority}  [{m.region_tag}]")
    scores = esm2_fold_scores(seq, props, top_n=10)
    if scores is None:
        print("\nESM2 fold-gate: DEFERRED (transformers/torch unavailable).")
    else:
        print("\nESM2 fold-gate (delta logP, top positions):")
        for k, v in list(scores.items())[:15]:
            print(f"  {k}: {v:+.3f}")
