from __future__ import annotations
import json, os, sys

import cobra
import evo2_nim
import annotate
from loop_evo2 import (nt_identity, gc_content, max_homopolymer, is_degenerate,
                       annotate_segment, PROMPT_BP, NUM_TOKENS, TOP_K, RANDOM_SEED,
                       DATA_DIR, CDS_CACHE)

TEMPERATURES = [0.7, 1.0, 1.3]
GENES = [("dapA", "DHDPS"), ("lysC", "ASPK")]
OUT_JSON = os.path.join(DATA_DIR, "evo2_highT_results.json")


def one(gene: str, expected_rxn: str, full_cds: str, temperature: float,
        model, genemap) -> dict:
    prompt = full_cds[:PROMPT_BP].upper()
    out = evo2_nim.generate_dna(prompt, num_tokens=NUM_TOKENS,
                                temperature=temperature, top_k=TOP_K,
                                random_seed=RANDOM_SEED)
    full_gen = (out.get("sequence", "") or "").upper()
    cont = full_gen[len(prompt):] if full_gen.startswith(prompt) else full_gen
    real_next = full_cds[PROMPT_BP:PROMPT_BP + len(cont)].upper()
    degen, reason = is_degenerate(cont)
    ann = annotate_segment(full_gen, model, genemap)
    rec = {
        "gene": gene, "expected_reaction": expected_rxn, "temperature": temperature,
        "continuation_bp": len(cont),
        "gc_continuation_pct": round(gc_content(cont) * 100, 1),
        "max_homopolymer": max_homopolymer(cont),
        "degenerate": degen, "degenerate_reason": reason,
        "drift_identity_pct": round(nt_identity(cont, real_next), 1),  # low = novel
        "top_reaction": ann["top_reaction"],
        "reaction_correct": ann["top_reaction"] == expected_rxn,
        "top_identity_pct": ann["top_identity_pct"],
        "cascade_confidence": ann["cascade_confidence"],
        "confidence_label": ann["top_confidence_label"],
        "still_annotates": (ann["top_reaction"] == expected_rxn) and not ann["low_confidence"],
    }
    print(f"[{gene} T={temperature}] drift={rec['drift_identity_pct']}% "
          f"rxn={rec['top_reaction']}({'ok' if rec['reaction_correct'] else 'X'}) "
          f"id={rec['top_identity_pct']}% conf={rec['cascade_confidence']} "
          f"novel_and_annotates={rec['still_annotates'] and rec['drift_identity_pct'] < 80}",
          flush=True)
    return rec


def run() -> int:
    assert cobra.__version__ == "0.31.1", f"cobra {cobra.__version__} != 0.31.1"
    cds = json.load(open(CDS_CACHE))
    model = cobra.io.load_json_model(annotate.MODEL_CACHE)
    genemap = annotate.load_genemap()
    results = []
    for gene, rxn in GENES:
        if gene not in cds:
            print(f"[skip] {gene} not in cds_cache.json"); continue
        for T in TEMPERATURES:
            try:
                results.append(one(gene, rxn, cds[gene], T, model, genemap))
            except Exception as e:
                print(f"[err] {gene} T={T}: {e}", flush=True)
                results.append({"gene": gene, "temperature": T, "error": str(e)})
    json.dump({"prompt_bp": PROMPT_BP, "num_tokens": NUM_TOKENS, "top_k": TOP_K,
               "random_seed": RANDOM_SEED, "temperatures": TEMPERATURES,
               "results": results}, open(OUT_JSON, "w"), indent=2)
    print(f"\nsaved -> {OUT_JSON}")
    # summary
    print("\n=== SUMMARY (novel = drift<80% AND still annotates correctly) ===")
    for r in results:
        if "error" in r: continue
        novel = r["drift_identity_pct"] < 80 and r["still_annotates"]
        print(f"  {r['gene']:5s} T={r['temperature']}: drift={r['drift_identity_pct']:5.1f}%  "
              f"correct_rxn={r['reaction_correct']}  conf={r['cascade_confidence']}  "
              f"NOVEL+ANNOTATES={novel}")
    return 0


if __name__ == "__main__":
    sys.exit(run())
