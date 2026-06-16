#!/usr/bin/env python3
from __future__ import annotations

import json
import os
from dataclasses import dataclass, field, asdict
from typing import Optional

import cobra  # noqa: F401  (version asserted in __main__)

# COMPOSE (import, never reimplement).
import operator_mutate as om
from operator_mutate import (
    fetch_uniprot,
    propose_mutations,
    parse_regulatory_features,
    MutationProposal,
)
import reverse_translate as rt
import loop as loop_mod
from loop import run_loop

CODE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(CODE_DIR, "data")
PROTEOME_FASTA = os.path.join(DATA_DIR, "ecoli_k12_proteome.fasta")
OUT_JSON = os.path.join(DATA_DIR, "proteome_gen_results.json")
OUT_FASTA = os.path.join(DATA_DIR, "proteome_gen_ranked.fasta")

# GC band the reverse-translator warns on (from reverse_translate defaults).
GC_LOW, GC_HIGH = 0.45, 0.60


# Pathway -> recruited enzyme map (L-lysine demo). dapA/DHDPS (P0A6L2) is the
# feedback-committed step -> fbr-relaxation router (modal Δ>0, plain Δ≈0).
# lysC/ASPK (P08660) = aspartate kinase III; 3 isozymes (lysC/thrA/metL) ->
# overexpression -> Δ=0 (FBA cannot reward a deregulated isozyme).
PATHWAY_ENZYMES: dict[str, list[dict]] = {
    "L-lysine": [
        {
            "gene": "dapA",
            "uniprot": "P0A6L2",
            "rxn_hint": "DHDPS",
            "role": "first committed, feedback-inhibited step of the DAP/lysine "
                    "branch; L-lysine allosteric inhibitor binds a Site cluster "
                    "(not an ACT domain); fbr-relaxation router (modal Δ>0).",
        },
        {
            "gene": "lysC",
            "uniprot": "P08660",
            "rxn_hint": "ASPK",
            "role": "aspartate kinase III; lysine-sensitive ACT domain. 3 isozymes "
                    "(lysC/thrA/metL); overexpression router -> Δ=0 (FBA cannot "
                    "reward a deregulated isozyme).",
        },
    ],
}


# Provenance container for a concrete variant protein.
@dataclass
class Variant:
    variant_id: str
    gene: str
    uniprot: str
    wt_aa: str
    position: int          # 1-based
    mut_aa: str
    mutation: str          # e.g. "M318I"
    region_tag: str        # WHERE-gate bucket that matched
    rationale: str
    where_rank: int        # rank from operator_mutate's ranked proposal list
    where_priority: float  # operator_mutate priority weight
    protein: str           # the concrete mutated protein string

    def header_prov(self) -> dict:
        return {
            "variant_id": self.variant_id,
            "gene": self.gene,
            "uniprot": self.uniprot,
            "mutation": self.mutation,
            "region_tag": self.region_tag,
            "where_rank": self.where_rank,
            "where_priority": self.where_priority,
            "rationale": self.rationale,
        }


# proteome FASTA reader (recruit source).
def _read_proteome(path: str = PROTEOME_FASTA) -> dict[str, dict]:
    """Parse the cached E. coli proteome FASTA -> {acc: {gene, seq, desc}}."""
    out: dict[str, dict] = {}
    acc = gene = desc = None
    buf: list[str] = []

    def _flush():
        if acc is not None:
            out[acc] = {"gene": gene, "seq": "".join(buf), "desc": desc}

    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                _flush()
                buf = []
                # >sp|P0A6L2|DAPA_ECOLI ... GN=dapA ...
                parts = line[1:].split("|")
                acc = parts[1] if len(parts) >= 2 else None
                desc = line[1:].rstrip("\n")
                gene = None
                for tok in line.split():
                    if tok.startswith("GN="):
                        gene = tok[3:]
                        break
            else:
                buf.append(line.strip())
    _flush()
    return out


# OPERATOR 1 - recruit
def recruit(target_pathway: str = "L-lysine") -> list[tuple[str, str, str]]:
    """Pull the pathway's enzymes from the cached proteome -> [(gene, uniprot_acc, protein_seq), ...]."""
    spec = PATHWAY_ENZYMES.get(target_pathway)
    if not spec:
        raise ValueError(f"No recruited enzyme set for pathway {target_pathway!r}")
    proteome = _read_proteome() if os.path.exists(PROTEOME_FASTA) else {}
    recruited: list[tuple[str, str, str]] = []
    for e in spec:
        acc = e["uniprot"]
        seq = (proteome.get(acc) or {}).get("seq", "")
        if not seq:
            # offline-safe fallback: cached UniProt json carries the sequence too.
            seq = fetch_uniprot(acc)["sequence"]["value"]
        recruited.append((e["gene"], acc, seq))
    return recruited


# WHERE-gate adapter. operator_mutate.propose_mutations handles the lysC-style
# gate (Binding-site with ligand=L-lysine, ACT Domain); DHDPS (dapA) annotates its
# lysine-inhibitor contacts as `Site` features instead. This adapter reuses
# operator_mutate's feature iterator + substitution heuristic + dataclass to
# surface those Site residues - a new feature selector, not new substitution logic.
_W_LYS_INHIBITOR_SITE = 3.0  # same weight class as a lysine binding-site contact


def _propose_lysine_inhibitor_sites(uniprot_json: dict) -> list[MutationProposal]:
    """Propose substitutions at `Site` residues marked L-lysine inhibitor (allosteric) contacts; works for dapA P0A6L2."""
    seq = uniprot_json["sequence"]["value"]
    props: list[MutationProposal] = []
    seen: set[int] = set()
    for f, s, e in om._iter_features(uniprot_json):  # reused iterator
        if f.get("type") != "Site":
            continue
        desc = (f.get("description") or "").lower()
        if "lysine" not in desc or "inhibitor" not in desc:
            continue
        for pos in range(s, e + 1):
            if pos in seen or not (1 <= pos <= len(seq)):
                continue
            seen.add(pos)
            wt = seq[pos - 1]
            cands = om._disrupt_candidates(wt)  # reused heuristic table
            props.append(MutationProposal(
                position=pos, wt_aa=wt, candidate_mutations=cands,
                rationale=(f"Annotated L-lysine allosteric INHIBITOR contact ({wt}); "
                           "disrupting effector binding relieves feedback inhibition "
                           "(fbr). [Site-feature WHERE-gate]"),
                region_tag="lysine_inhibitor_site",
                priority=_W_LYS_INHIBITOR_SITE,
            ))
    props.sort(key=lambda m: (-m.priority, m.position))
    return props


def _where_gate(uniprot_json: dict) -> list[MutationProposal]:
    """Run the regulatory WHERE-gate: operator_mutate native proposals (lysC/ACT), else the Site-feature lysine-inhibitor selector (dapA)."""
    props = propose_mutations(uniprot_json, strategy="lysine_site")
    if props:
        return props
    return _propose_lysine_inhibitor_sites(uniprot_json)


# OPERATOR 2 - mutate
def mutate(protein_seq: str,
           uniprot_json: dict,
           n_variants: int = 3,
           gene: str = "?",
           uniprot: str = "?") -> list[Variant]:
    """Apply WHERE-gate substitutions onto the WT protein -> concrete variant proteins (one deregulating point substitution each)."""
    wt = uniprot_json["sequence"]["value"]
    if protein_seq and protein_seq != wt:
        # Position indices come from the UniProt sequence; align to it. The cached
        # proteome seq and UniProt seq are the same entry; this is explicit to
        # avoid an off-by-one.
        wt = wt  # authoritative
    proposals = _where_gate(uniprot_json)
    variants: list[Variant] = []
    for where_rank, m in enumerate(proposals[:n_variants], start=1):
        mut_aa = m.candidate_mutations[0]  # top disruptive substitution
        pos = m.position
        if not (1 <= pos <= len(wt)) or wt[pos - 1] != m.wt_aa:
            continue
        var_seq = wt[:pos - 1] + mut_aa + wt[pos:]
        mutation = f"{m.wt_aa}{pos}{mut_aa}"
        variants.append(Variant(
            variant_id=f"{gene}_{mutation}",
            gene=gene, uniprot=uniprot,
            wt_aa=m.wt_aa, position=pos, mut_aa=mut_aa, mutation=mutation,
            region_tag=m.region_tag, rationale=m.rationale,
            where_rank=where_rank, where_priority=m.priority,
            protein=var_seq,
        ))
    return variants


# OPERATOR 3 - recombine  (MVP STUB / future work - not implemented)
def recombine(variants: list[Variant], *args, **kwargs):
    """OPTIONAL domain-aware recombination (stub / future work). Returns input variants unchanged with stub status."""
    return {
        "status": "NotImplemented (future-work stub)",
        "reason": "domain-aware recombination needs domain boundaries + an "
                  "epistasis/fold filter (ESM2 fold-gate, currently deferred). "
                  "Single-substitution variants suffice to close the loop.",
        "variants": variants,
    }


# variant protein -> lab-ready DNA (reverse_translate), GC-band aware.
def _to_lab_ready_dna(variant: Variant, seed: int = 1234) -> dict:
    """reverse_translate the variant protein; retry 'weighted' if GC outside 0.45-0.60."""
    res = rt.reverse_translate(variant.protein, strategy="max", seed=seed)
    used = "max"
    if not res["gc_in_range"]:
        alt = rt.reverse_translate(variant.protein, strategy="weighted", seed=seed)
        if alt["gc_in_range"] or abs(alt["gc"] - 0.525) < abs(res["gc"] - 0.525):
            res, used = alt, "weighted"
    res["strategy_used"] = used
    return res


# annotatability check - does the MUTATED DNA still map to the intended reaction?
# Reuses loop's annotate-cascade machinery (annotate.annotate).
def _annotatability(dna: str) -> dict:
    """Annotate variant DNA and report reaction mapping + identity (mutation should not break annotatability)."""
    import annotate
    genemap = annotate.load_genemap()
    mc = loop_mod.ModalContext.load()
    res = annotate.annotate(dna, mc.base, genemap=genemap)
    top_rxn = res.reactions[0].reaction_id if res.reactions else None
    ident = res.assignments[0].identity if res.assignments else 0.0
    label = res.assignments[0].confidence_label if res.assignments else "none"
    return {
        "recovered_reaction": top_rxn,
        "identity_pct": round(float(ident), 2),
        "cascade_confidence": round(float(res.confidence), 4),
        "confidence_label": label,
        "low_confidence": bool(res.low_confidence),
    }


# generate - the slot-2 entrypoint
def generate(target_id: str = "EX_lys__L_e",
             target_pathway: str = "L-lysine",
             n_variants: int = 5,
             per_enzyme_variants: int = 3,
             growth_fraction: float = loop_mod.GROWTH_FRACTION_DEFAULT,
             throttle_frac: float = loop_mod.THROTTLE_FRAC_DEFAULT,
             write_outputs: bool = True) -> dict:
    """Full proteome-anchored generative loop: recruit -> mutate -> reverse_translate -> run_loop -> ranked."""
    os.makedirs(DATA_DIR, exist_ok=True)

    # 1. recruit
    recruited = recruit(target_pathway)

    # 2. mutate each recruited enzyme via the WHERE-gate
    all_variants: list[Variant] = []
    recruit_report: list[dict] = []
    for gene, acc, seq in recruited:
        j = fetch_uniprot(acc)
        vs = mutate(seq, j, n_variants=per_enzyme_variants, gene=gene, uniprot=acc)
        all_variants.extend(vs)
        recruit_report.append({
            "gene": gene, "uniprot": acc, "protein_len": len(seq),
            "n_where_gate_variants": len(vs),
            "variants": [v.mutation for v in vs],
        })
    # cap to n_variants total; the per-enzyme cap already balances, so truncate
    # deterministically.
    all_variants = all_variants[:n_variants] if n_variants else all_variants

    # 3. recombine (MVP stub - recorded, not applied)
    recombine_status = recombine(all_variants)["status"]

    # 4. reverse_translate each variant -> lab-ready DNA
    dna_inputs: list[tuple[str, str]] = []
    variant_index: dict[str, dict] = {}
    for v in all_variants:
        rtres = _to_lab_ready_dna(v)
        dna = rtres["dna"]
        dna_inputs.append((v.variant_id, dna))
        variant_index[v.variant_id] = {
            "provenance": v.header_prov(),
            "dna_len_bp": len(dna),
            "gc": round(rtres["gc"], 4),
            "gc_in_range": rtres["gc_in_range"],
            "rt_strategy": rtres["strategy_used"],
            "banned_sites": rtres["banned_sites"],
            "annotatability": _annotatability(dna),
            "protein": v.protein,
            "dna": dna,
        }

    # 5. score + rank through the closed loop (plain + modal FBA oracle).
    # NOTE: run_loop writes loop_results.json/loop_ranked.fasta as a side-effect;
    # we keep our own proteome_gen_* outputs and let run_loop manage its own.
    ranked = run_loop(target_id, dna_inputs=dna_inputs,
                      growth_fraction=growth_fraction,
                      throttle_frac=throttle_frac,
                      write_outputs=True)

    # 6. weave provenance + annotatability back into the ranked rows
    for r in ranked:
        vid = r["id"]
        vinfo = variant_index.get(vid, {})
        r["variant_provenance"] = vinfo.get("provenance")
        r["annotatability"] = vinfo.get("annotatability")
        r["gc"] = vinfo.get("gc")
        r["rt_strategy"] = vinfo.get("rt_strategy")

    payload = {
        "title": "proteome-anchored generator (local, no API, no GPU): "
                 "recruit -> mutate(WHERE-gate) -> "
                 "reverse_translate -> plain+modal FBA oracle -> ranked candidates",
        "target_id": target_id,
        "target_pathway": target_pathway,
        "recruited_enzymes": recruit_report,
        "recombine_status": recombine_status,
        "n_variants_scored": len(dna_inputs),
        "caveats": [
            "(H1) Modal Δ reflects the intervention TYPE the recovered reaction "
            "routes to (DHDPS->fbr-relaxation; ASPK->overexpression->Δ=0), not the "
            "specific point mutation. Stoichiometric FBA is blind to sequence/kcat, "
            "so this does not claim FBA validated the mutation. The mutation's "
            "deregulating-plausibility is a separate evidence stream from the "
            "operator_mutate WHERE-gate (benchmarked there: recall=1.0 vs documented "
            "fbr ground-truth). Two independent evidence streams.",
            "(H2) In-silico re-ranking + lab-ready candidate emission, not a "
            "viability/titre/growth prediction.",
            "(H3) ESM2 fold-gate deferred (transformers/torch not installed); "
            "variants are not yet fold-filtered. The WHERE-gate supplies intent; "
            "ESM2 would supply fold safety. Known gap.",
            "(H4) Near-natural bias: proteome-anchored variants are 1 substitution "
            "from a real folded E. coli enzyme -> low novelty, low fold-risk, high "
            "orderability. Opposite axis-end from Evo2's de-novo explore mode.",
            "(H5) fbr throttle magnitude is illustrative (needs a measured Ki); the "
            "qualitative claim: plain Δ≈0, modal Δ>0 when the throttle binds.",
        ],
        "results": [
            {k: v for k, v in r.items() if k != "dna"} for r in ranked
        ],
        "cobra_version": cobra.__version__,
    }

    if write_outputs:
        with open(OUT_JSON, "w") as f:
            json.dump(payload, f, indent=2)
        _write_fasta(ranked, variant_index, target_id)

    return {
        "payload": payload,
        "ranked": ranked,
        "variant_index": variant_index,
        "recruited": recruited,
        "variants": all_variants,
    }


# ranked lab-ready FASTA (provenance-rich header).
def _write_fasta(ranked: list[dict], variant_index: dict, target_id: str) -> None:
    """Write ranked lab-ready candidate DNA with provenance-rich FASTA headers (scored rows only, in rank order)."""
    scored = [r for r in ranked if r.get("scored")]
    scored.sort(key=lambda r: (r.get("rank") is None, r.get("rank") or 1e9))
    lines: list[str] = []
    for r in scored:
        vid = r["id"]
        vinfo = variant_index.get(vid, {})
        prov = vinfo.get("provenance", {}) or {}
        annot = vinfo.get("annotatability", {}) or {}
        hdr = (
            f">rank{r.get('rank')}|{vid}"
            f"|gene={prov.get('gene')}|uniprot={prov.get('uniprot')}"
            f"|mut={prov.get('mutation')}|where_rank={prov.get('where_rank')}"
            f"|rxn={r.get('recovered_reaction')}"
            f"|intervention={r.get('intervention_type')}"
            f"|plainD={r.get('plain_delta')}|modalD={r.get('modal_delta')}"
            f"|conf={r.get('confidence')}"
            f"|annot_id={annot.get('identity_pct')}%"
            f"|gc={vinfo.get('gc')}"
            f" lab-ready | NOTE: modalΔ=intervention-type evidence, "
            f"mutation-plausibility=WHERE-gate evidence (separate streams)"
        )
        lines.append(hdr)
        seq = vinfo.get("dna", "")
        for i in range(0, len(seq), 70):
            lines.append(seq[i:i + 70])
    with open(OUT_FASTA, "w") as f:
        f.write("\n".join(lines) + ("\n" if lines else ""))


# pretty report
def print_report(out: dict, target_id: str) -> None:
    ranked = out["ranked"]
    line = "-" * 118
    print("\n" + line)
    print("PROTEOME-ANCHORED GENERATOR (slot 2) - ranked candidates for %s" % target_id)
    print(line)
    hdr = "%-5s | %-14s | %-9s | %-9s | %-16s | %10s | %10s | %5s | %7s"
    print(hdr % ("rank", "variant", "mut", "rxn", "intervention",
                 "plain Δ", "modal Δ", "conf", "annot%"))
    print(line)
    row = "%-5s | %-14s | %-9s | %-9s | %-16s | %10s | %10s | %5s | %7s"
    for r in ranked:
        prov = r.get("variant_provenance") or {}
        annot = r.get("annotatability") or {}
        rank = str(r.get("rank")) if r.get("rank") is not None else "-"
        pd = "%+.5f" % r["plain_delta"] if r.get("plain_delta") is not None else "n/a"
        md = "%+.5f" % r["modal_delta"] if r.get("modal_delta") is not None else "n/a"
        print(row % (rank, r["id"], prov.get("mutation", "?"),
                     str(r.get("recovered_reaction")),
                     r.get("intervention_type") or "-", pd, md,
                     "%.2f" % r.get("confidence", 0.0),
                     "%.1f" % (annot.get("identity_pct") or 0.0)))
    print(line)
    print("JSON  ->", OUT_JSON)
    print("FASTA ->", OUT_FASTA)
    print(line + "\n")


if __name__ == "__main__":
    assert cobra.__version__ == "0.31.1", f"cobra {cobra.__version__} != 0.31.1"
    out = generate(target_id="EX_lys__L_e", n_variants=5)
    print_report(out, "EX_lys__L_e")
    print("cobra version:", cobra.__version__)
