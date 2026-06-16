from __future__ import annotations

import json
import os

import cobra  # noqa: F401  (version asserted in __main__)

import evo2_nim
import annotate
from loop import run_loop, print_report

CODE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(CODE_DIR, "data")
CDS_CACHE = os.path.join(DATA_DIR, "cds_cache.json")
OUT_JSON = os.path.join(DATA_DIR, "loop_evo2_results.json")
OUT_FASTA = os.path.join(DATA_DIR, "loop_evo2_ranked.fasta")

TARGET = "EX_lys__L_e"

# Prompt genes -> their known iML1515 reaction (ground truth for comparison).
# Restricted to the lysine pathway (dapA primary; lysC secondary).
PROMPT_GENES = [
    ("dapA", {"DHDPS"}),   # first committed, single-gene lysine-specific step
    ("lysC", {"ASPK"}),    # aspartate kinase III (isozyme-bearing)
]

PROMPT_BP = 400        # length of real-gene prefix used as the Evo2 prompt
NUM_TOKENS = 300       # continuation length requested from Evo2
TEMPERATURE = 0.7
TOP_K = 3
RANDOM_SEED = 42       # reproducibility

# Degeneracy / realism gates for the GENERATED continuation.
GC_MIN, GC_MAX = 0.35, 0.65    # realistic E. coli coding GC band
MAX_HOMOPOLYMER = 12           # reject if a single base repeats > this many times
MIN_DISTINCT_FRAC = 0.0        # (kept explicit; we require all 4 bases present below)


# sequence helpers (pure; no external deps)
def gc_content(seq: str) -> float:
    seq = seq.upper()
    if not seq:
        return 0.0
    gc = sum(1 for b in seq if b in "GC")
    return gc / len(seq)


def max_homopolymer(seq: str) -> int:
    seq = seq.upper()
    best = run = 0
    prev = ""
    for b in seq:
        run = run + 1 if b == prev else 1
        prev = b
        best = max(best, run)
    return best


def nt_identity(a: str, b: str) -> float:
    """% identical bases over overlapping length (position-wise, no alignment). Drift metric: continuation vs real gene next bases. High -> Evo2 reproduced gene; low -> drifted to novel sequence."""
    a, b = a.upper(), b.upper()
    n = min(len(a), len(b))
    if n == 0:
        return 0.0
    same = sum(1 for i in range(n) if a[i] == b[i])
    return 100.0 * same / n


def is_degenerate(seq: str) -> tuple[bool, str]:
    """Return (degenerate?, reason). Realistic coding continuation needs GC in band, all four bases, no absurd homopolymer run."""
    seq = seq.upper()
    if len(seq) < 30:
        return True, f"too short ({len(seq)} nt)"
    gc = gc_content(seq)
    if not (GC_MIN <= gc <= GC_MAX):
        return True, f"GC={gc*100:.1f}% outside realistic band {GC_MIN*100:.0f}-{GC_MAX*100:.0f}%"
    present = set(b for b in seq if b in "ACGT")
    if present != set("ACGT"):
        return True, f"not all 4 bases present ({''.join(sorted(present))})"
    hp = max_homopolymer(seq)
    if hp > MAX_HOMOPOLYMER:
        return True, f"homopolymer run {hp} > {MAX_HOMOPOLYMER}"
    return False, "ok"


# Evo2 generation (real API) - prompt-seeded continuation
def generate_for_gene(gene: str, full_cds: str) -> dict:
    """Prompt Evo2 with first PROMPT_BP of real gene; return record with prompt, full output, continuation-only segment, GC%, degeneracy, and drift metric (continuation vs real-gene next bases)."""
    prompt = full_cds[:PROMPT_BP].upper()
    print(f"[evo2] {gene}: prompt={len(prompt)}bp num_tokens={NUM_TOKENS} seed={RANDOM_SEED}",
          flush=True)
    out = evo2_nim.generate_dna(prompt, num_tokens=NUM_TOKENS,
                                temperature=TEMPERATURE, top_k=TOP_K,
                                random_seed=RANDOM_SEED)
    full_gen = (out.get("sequence", "") or "").upper()
    # NIM returns prompt+continuation; isolate the continuation.
    if full_gen.startswith(prompt):
        cont = full_gen[len(prompt):]
    else:
        cont = full_gen
    # Real gene's next bases after the prompt (for the drift / identity metric).
    real_next = full_cds[PROMPT_BP:PROMPT_BP + len(cont)].upper()

    degen, reason = is_degenerate(cont)
    rec = {
        "gene": gene,
        "prompt_bp": len(prompt),
        "full_gen_bp": len(full_gen),
        "continuation_bp": len(cont),
        "gc_full": round(gc_content(full_gen) * 100, 1),
        "gc_continuation": round(gc_content(cont) * 100, 1),
        "max_homopolymer_cont": max_homopolymer(cont),
        "degenerate": degen,
        "degenerate_reason": reason,
        "cont_vs_realgene_identity_pct": round(nt_identity(cont, real_next), 1),
        "elapsed_ms": out.get("elapsed_ms"),
        "prompt": prompt,
        "full_gen": full_gen,
        "continuation": cont,
        "real_next": real_next,
    }
    print(f"   -> full={len(full_gen)}bp cont={len(cont)}bp "
          f"GC(cont)={rec['gc_continuation']}% degenerate={degen} ({reason}) "
          f"drift_identity={rec['cont_vs_realgene_identity_pct']}% vs real gene")
    return rec


# annotate comparison: prompt-gene KNOWN reaction vs Evo2 output reaction
def annotate_segment(dna: str, model, genemap) -> dict:
    res = annotate.annotate(dna, model, genemap=genemap)
    top_rxn = res.reactions[0] if res.reactions else None
    top_asn = res.assignments[0] if res.assignments else None
    return {
        "n_orfs": len(res.orfs),
        "recovered_reactions": [r.reaction_id for r in res.reactions],
        "top_reaction": top_rxn.reaction_id if top_rxn else None,
        "cascade_confidence": res.confidence,
        "low_confidence": res.low_confidence,
        "top_assigned_gene": top_asn.gene if top_asn else None,
        "top_assigned_gene_name": top_asn.gene_name if top_asn else None,
        "top_identity_pct": round(top_asn.identity, 1) if top_asn else None,
        "top_method": top_asn.method if top_asn else None,
        "top_confidence_label": top_asn.confidence_label if top_asn else None,
        "method": res.method,
    }


# main
def run() -> int:
    os.makedirs(DATA_DIR, exist_ok=True)
    cds = json.load(open(CDS_CACHE))

    # shared annotate resources (cached model + genemap)
    print("DIAMOND available:", annotate._diamond_bin() or "NO -> python-pairwise fallback")
    annot_model = cobra.io.load_json_model(annotate.MODEL_CACHE)
    genemap = annotate.load_genemap()

    print("=" * 80)
    print("STEP 1-2: Evo2-40b generation (prompt-seeded) + annotate")
    print("=" * 80)

    records: list[dict] = []
    evo2_inputs: list[tuple[str, str]] = []   # (id, dna) fed into run_loop
    for gene, known_rxns in PROMPT_GENES:
        if gene not in cds:
            print(f"[skip] {gene}: not in cds_cache.json")
            continue
        gen = generate_for_gene(gene, cds[gene])

        # annotate BOTH the full output and the continuation-only segment
        ann_full = annotate_segment(gen["full_gen"], annot_model, genemap)
        ann_cont = annotate_segment(gen["continuation"], annot_model, genemap)

        evo2_id = f"evo2_{gene}_seed{RANDOM_SEED}"
        preserved = bool(known_rxns & set(ann_full["recovered_reactions"]))

        rec = {
            "evo2_id": evo2_id,
            "prompt_gene": gene,
            "prompt_gene_known_reaction": sorted(known_rxns),
            **gen,
            "annotate_full": ann_full,
            "annotate_continuation_only": ann_cont,
            "preserves_prompt_gene_reaction": preserved,
        }
        records.append(rec)

        print(f"   annotate(full): top_rxn={ann_full['top_reaction']} "
              f"gene={ann_full['top_assigned_gene']} id={ann_full['top_identity_pct']}% "
              f"conf={ann_full['cascade_confidence']} lowC={ann_full['low_confidence']} "
              f"=> preserves {sorted(known_rxns)}? {preserved}")
        print(f"   annotate(cont): top_rxn={ann_cont['top_reaction']} "
              f"recovered={ann_cont['recovered_reactions']} "
              f"conf={ann_cont['cascade_confidence']} lowC={ann_cont['low_confidence']}")

        # Feed the full Evo2 output into the loop unless degenerate.
        if gen["degenerate"]:
            print(f"   [GATE] {evo2_id} degenerate ({gen['degenerate_reason']}) "
                  f"-> EXCLUDED from loop scoring (reported only).")
        else:
            evo2_inputs.append((evo2_id, gen["full_gen"]))

    # STEP 3: run the closed loop on the Evo2 DNA
    print("\n" + "=" * 80)
    print(f"STEP 3: run_loop(target={TARGET}) on {len(evo2_inputs)} non-degenerate "
          f"Evo2 sequence(s)")
    print("=" * 80)

    ranked = []
    if evo2_inputs:
        # run_loop writes loop_results.json/loop_ranked.fasta; we want our OWN files,
        # so disable its writes and serialize here.
        ranked = run_loop(TARGET, evo2_inputs, write_outputs=False)
        print_report(ranked, TARGET)
    else:
        print("No non-degenerate Evo2 sequences to score.")

    # index loop results by id for joining back into per-sequence records
    loop_by_id = {r["id"]: r for r in ranked}
    closed_arc = False
    for rec in records:
        lr = loop_by_id.get(rec["evo2_id"])
        if lr is None:
            rec["loop"] = {"scored": False,
                           "skip_reason": "excluded (degenerate) - not scored"}
            continue
        rec["loop"] = {
            "rank": lr.get("rank"),
            "recovered_reaction": lr.get("recovered_reaction"),
            "intervention_type": lr.get("intervention_type"),
            "confidence": lr.get("confidence"),
            "low_confidence": lr.get("low_confidence"),
            "plain_delta": lr.get("plain_delta"),
            "modal_delta": lr.get("modal_delta"),
            "growth_retained": lr.get("growth_retained"),
            "scored": lr.get("scored"),
            "skip_reason": lr.get("skip_reason"),
        }
        if lr.get("scored") and lr.get("modal_delta") is not None and lr["modal_delta"] > 1e-6:
            closed_arc = True

    # write artifacts
    payload = {
        "title": "Evo2-40b DNA -> annotate -> Candidate -> "
                 "plain+modal FBA oracle -> ranking",
        "target_id": TARGET,
        "generation": {
            "model": "nvidia evo2-40b (NIM hosted)",
            "prompt_bp": PROMPT_BP,
            "num_tokens": NUM_TOKENS,
            "temperature": TEMPERATURE,
            "top_k": TOP_K,
            "random_seed": RANDOM_SEED,
            "prompt_genes": [g for g, _ in PROMPT_GENES],
        },
        "degeneracy_gates": {
            "gc_band": [GC_MIN, GC_MAX],
            "max_homopolymer": MAX_HOMOPOLYMER,
            "require_all_four_bases": True,
        },
        "caveats": [
            "Evo2 generation is prompt-seeded on E. coli CDS; the continuation is "
            "the generated novel segment.",
            "Degenerate continuations (GC out of band / homopolymer / missing bases) "
            "are excluded from scoring but reported.",
            "cont_vs_realgene_identity is a drift metric: high => Evo2 reproduced the "
            "gene (annotatable), low => it drifted (may be unmappable).",
            "Cascade confidence propagates to ranking; low_confidence is flagged.",
            "This is in-silico re-ranking through the generative path, not a "
            "viability / titre prediction.",
        ],
        "closed_arc": closed_arc,
        "n_sequences": len(records),
        "n_scored": sum(1 for r in records if r["loop"].get("scored")),
        # strip bulky DNA strings from the top-level results array; keep in FASTA
        "results": [_strip_seqs(r) for r in records],
        "cobra_version": cobra.__version__,
    }
    with open(OUT_JSON, "w") as f:
        json.dump(payload, f, indent=2)

    _write_fasta(records, ranked)

    # final console summary
    print("\n" + "=" * 80)
    print("GENERATIVE ARC SUMMARY")
    print("=" * 80)
    hdr = (f"{'prompt':6s} {'known':8s} {'cont_bp':>7s} {'GC%':>5s} {'drift_id%':>9s} "
           f"{'evo2_rxn':10s} {'conf':>5s} {'lowC':>4s} {'modalΔ':>9s} {'preserved':>9s}")
    print(hdr)
    print("-" * len(hdr))
    for r in records:
        known = ",".join(r["prompt_gene_known_reaction"])
        af = r["annotate_full"]
        lp = r["loop"]
        md = lp.get("modal_delta")
        md_s = f"{md:+.5f}" if isinstance(md, (int, float)) else "n/a"
        deg = " DEGEN" if r["degenerate"] else ""
        print(f"{r['prompt_gene']:6s} {known:8s} {r['continuation_bp']:7d} "
              f"{r['gc_continuation']:5.1f} {r['cont_vs_realgene_identity_pct']:9.1f} "
              f"{str(af['top_reaction']):10s} {af['cascade_confidence']:5.2f} "
              f"{'Y' if af['low_confidence'] else 'n':>4s} {md_s:>9s} "
              f"{str(r['preserves_prompt_gene_reaction']):>9s}{deg}")
    print("-" * len(hdr))
    print(f"closed_arc (any modal Δ>0 from real Evo2 DNA): {closed_arc}")
    print(f"JSON  -> {OUT_JSON}")
    print(f"FASTA -> {OUT_FASTA}")
    print(f"cobra version: {cobra.__version__}")
    return 0


def _strip_seqs(rec: dict) -> dict:
    """Drop bulky raw-DNA strings from the JSON results array (kept in FASTA)."""
    drop = {"prompt", "full_gen", "continuation", "real_next"}
    out = {k: v for k, v in rec.items() if k not in drop}
    # keep short hashes/lengths already present; nothing else to do
    return out


def _write_fasta(records: list[dict], ranked: list[dict]) -> None:
    """All Evo2 outputs as FASTA, headers carry generation + annotate + loop provenance. Scored sequences ordered by loop rank first."""
    rank_by_id = {r["id"]: r.get("rank") for r in ranked}
    # order: scored-by-rank, then degenerate/unscored
    def sort_key(r):
        rk = rank_by_id.get(r["evo2_id"])
        return (0, rk) if isinstance(rk, int) else (1, 9999)
    ordered = sorted(records, key=sort_key)

    lines: list[str] = []
    for r in ordered:
        lp = r["loop"]
        af = r["annotate_full"]
        rk = lp.get("rank")
        rk_s = f"rank{rk}" if isinstance(rk, int) else ("DEGENERATE" if r["degenerate"]
                                                         else "unscored")
        md = lp.get("modal_delta")
        md_s = f"{md:+.5f}" if isinstance(md, (int, float)) else "n/a"
        hdr = (f">{r['evo2_id']}|{rk_s}|prompt={r['prompt_gene']}"
               f"|known={','.join(r['prompt_gene_known_reaction'])}"
               f"|evo2_rxn={af['top_reaction']}|conf={af['cascade_confidence']}"
               f"|drift_id={r['cont_vs_realgene_identity_pct']}%"
               f"|GCcont={r['gc_continuation']}%|modalD={md_s}"
               f"|preserved={r['preserves_prompt_gene_reaction']}")
        lines.append(hdr)
        seq = r["full_gen"]
        for i in range(0, len(seq), 70):
            lines.append(seq[i:i + 70])
    with open(OUT_FASTA, "w") as f:
        f.write("\n".join(lines) + ("\n" if lines else ""))


if __name__ == "__main__":
    assert cobra.__version__ == "0.31.1", f"cobra must be 0.31.1, got {cobra.__version__}"
    raise SystemExit(run())
