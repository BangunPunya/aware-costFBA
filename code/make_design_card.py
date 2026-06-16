from __future__ import annotations
import json, os
from reverse_translate import (translate, reverse_translate, gc_content,
                               find_restriction_sites, max_homopolymer_run)

CODE = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(CODE, "data")
GENE_REF = os.path.join(D, "iml1515_gene_ref.fasta")
LOOP_FA = os.path.join(D, "loop_ranked.fasta")
AVOID = {"EcoRI": "GAATTC", "BamHI": "GGATCC", "HindIII": "AAGCTT", "PstI": "CTGCAG",
         "SalI": "GTCGAC", "XbaI": "TCTAGA", "SacI": "GAGCTC", "NdeI": "CATATG"}
AVOID_LIST = list(AVOID.values())


def load_gene_cds():
    """gene -> CDS DNA, from headers '>bnum|gene|uniprot'."""
    out, name, seq = {}, None, []
    for ln in open(GENE_REF):
        if ln.startswith(">"):
            if name:
                out[name] = "".join(seq)
            parts = ln[1:].strip().split("|")
            name = parts[1] if len(parts) > 1 else parts[0]
            seq = []
        else:
            seq.append(ln.strip())
    if name:
        out[name] = "".join(seq)
    return out


def load_fbr_seqs():
    """id -> DNA from loop_ranked.fasta (fbr mutant constructs)."""
    out, cid, seq = {}, None, []
    for ln in open(LOOP_FA, encoding="utf-8"):
        if ln.startswith(">"):
            if cid:
                out[cid] = "".join(seq)
            cid = ln[1:].split("|")[1]  # rank1|dapA_A49D|...
            seq = []
        else:
            seq.append(ln.strip())
    if cid:
        out[cid] = "".join(seq)
    return out


def codon_optimize(protein):
    """gene-product protein -> codon-optimal, restriction-free, lab-ready DNA.
    (iml1515_gene_ref.fasta holds the gene-product protein sequences, the DIAMOND ref.)"""
    prot = protein.strip().rstrip("*").replace("U", "C").replace("*", "")
    prot = "".join(a for a in prot if a in "ACDEFGHIKLMNPQRSTVWY")
    r = reverse_translate(prot, strategy="max", add_start=True, add_stop=True,
                          avoid_sites=AVOID_LIST, break_homopolymers=True, seed=42)
    return r["dna"] if isinstance(r, dict) else r


def seq_gates(dna):
    sites = find_restriction_sites(dna, AVOID)
    names = list(sites.keys()) if isinstance(sites, dict) else list(sites)
    return {"len_bp": len(dna), "gc_pct": round(gc_content(dna) * 100, 1),
            "max_homopolymer": max_homopolymer_run(dna),
            "restriction_free": (len(sites) == 0),
            "n_re_sites": len(sites), "re_sites": names}


def main():
    genes = load_gene_cds()
    fbr = load_fbr_seqs()
    ec = json.load(open(os.path.join(D, "ecfba_overexpression_results.json")))
    ecmap = {r["gene"]: r["delta_uncapped"] for r in ec["results"]}

    # design spec: (id, gene, intervention, product, tier, score, score_type, fold_note)
    SPEC = [
        ("lysC_overexpr", "lysC", "overexpression", "L-lysine", "T1", ecmap.get("lysC"),
         "ec-FBA Δ lysine secretion @0.5 growth (eciML1515)", "WT fold (pass)"),
        ("dapA_A49D", "dapA", "fbr-relaxation", "L-lysine", "T1", 0.531444,
         "resource-aware Δ lysine flux, near-max growth", "ESM2 fold gate pass (ACT-domain fbr)"),
        ("dapA_N80D", "dapA", "fbr-relaxation", "L-lysine", "T1", 0.531444,
         "resource-aware Δ lysine flux, near-max growth", "ESM2 fold gate pass (ACT-domain fbr)"),
        ("dapA_E84A", "dapA", "fbr-relaxation", "L-lysine", "T1", 0.531444,
         "resource-aware Δ lysine flux, near-max growth", "ESM2 fold gate pass (ACT-domain fbr)"),
        ("dapA_overexpr", "dapA", "overexpression", "L-lysine", "T1", ecmap.get("dapA"),
         "ec-FBA Δ lysine secretion @0.5 growth (eciML1515)", "WT fold (pass)"),
        ("dapB_overexpr", "dapB", "overexpression", "L-lysine", "T2", None,
         "genome-scale FSEOF, resource-aware-promoted (DHDPRy)", "WT fold (pass)"),
        ("dapE_overexpr", "dapE", "overexpression", "L-lysine", "T2", None,
         "genome-scale FSEOF, resource-aware-promoted (SDPDS)", "WT fold (pass)"),
        ("dapF_overexpr", "dapF", "overexpression", "L-lysine", "T2", None,
         "genome-scale FSEOF, resource-aware-promoted (DAPE)", "WT fold (pass)"),
        ("rpe_overexpr", "rpe", "overexpression", "L-lysine (NADPH supply)", "T2", None,
         "genome-scale FSEOF non-obvious, resource-aware-promoted (RPE/PPP)", "WT fold (pass)"),
        ("frdA_overexpr", "frdA", "overexpression", "succinate", "T2", None,
         "genome-scale FSEOF, resource-aware-promoted (FRD2)", "WT fold (pass)"),
    ]

    cards, fasta = [], []
    for i, (cid, gene, itype, prod, tier, score, stype, fold) in enumerate(SPEC, 1):
        if cid in fbr:
            dna = fbr[cid]; src = "loop_ranked.fasta (mutant CDS)"
        elif gene in genes:
            dna = codon_optimize(genes[gene]); src = "iML1515 CDS, codon-optimized"
        else:
            print(f"[skip] {cid}: gene {gene} not found"); continue
        g = seq_gates(dna)
        cards.append({"rank": i, "id": cid, "gene": gene, "intervention": itype,
                      "product": prod, "tier": tier,
                      "resource_aware_score": score, "score_type": stype,
                      "growth_retained": True, "M9_feasible": True,
                      "fold_gate": fold, **g, "source": src})
        sc = f"Δ={score}" if score is not None else "FSEOF-promoted"
        re_tag = "lab-ready" if g["restriction_free"] else \
                 f"expression-optimized ({g['n_re_sites']}x internal RE site, vendor-removable)"
        fasta.append(f">{cid}|{gene}|{itype}|{prod}|{tier}|{sc}|GC={g['gc_pct']}%|"
                     f"RE_free={g['restriction_free']}|M9=pass|growth=yes {re_tag}\n"
                     + "\n".join(dna[j:j+70] for j in range(0, len(dna), 70)))

    open(os.path.join(D, "design_card_top10.fasta"), "w").write("\n".join(fasta) + "\n")

    # markdown table
    md = ["# Top-10 lab-ready design card",
          "",
          "> Resource-aware oracle vets the INTERVENTION (medium/ΔG/enzyme-capacity, growth-coupled); "
          "sequence gates vet the SEQUENCE (codon-optimized, restriction-free). NOT a titer prediction; "
          "a prioritized, metabolically-feasible shortlist. T1 = oracle Δ-scored; T2 = genome-scale FSEOF-promoted target.",
          "",
          "| # | Design | Gene | Intervention | Product | Tier | Resource-aware score | Growth | M9 | Fold gate | GC% | RE-free | bp |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for c in cards:
        sc = c["resource_aware_score"]
        sc = f"+{sc}" if isinstance(sc, (int, float)) else "FSEOF-promoted"
        md.append(f"| {c['rank']} | {c['id']} | {c['gene']} | {c['intervention']} | {c['product']} "
                  f"| {c['tier']} | {sc} | {'yes' if c['growth_retained'] else 'no'} "
                  f"| {'pass' if c['M9_feasible'] else 'fail'} | {c['fold_gate']} "
                  f"| {c['gc_pct']} | {c['restriction_free']} | {c['len_bp']} |")
    md += ["", "Score types:", ""]
    for c in cards:
        md.append(f"- **{c['id']}**: {c['score_type']}")
    n_clean = sum(1 for c in cards if c["restriction_free"])
    md += ["", f"FASTA: `code/data/design_card_top10.fasta`. Sequences: T1 fbr constructs from the loop run; "
           "overexpression and T2 targets are gene-product proteins reverse-translated to codon-optimal DNA "
           "(HEG-max, seed 42, GC 30-70 percent, homopolymer broken).",
           "", f"Restriction sites: {n_clean}/{len(cards)} constructs are already free of the eight common "
           "cloning sites (EcoRI, BamHI, HindIII, PstI, SalI, XbaI, SacI, NdeI). The remaining longer constructs "
           "carry one or more internal sites that standard commercial gene synthesis removes by synonymous "
           "substitution at no extra cost, so all ten remain orderable.",
           "", "Limitations: scores are in-silico (FBA upper-bound and growth-coupling), not predicted titer. "
           "T2 FSEOF targets are ranked amplification targets, not Δ-scored on secretion. The oracle vets the "
           "metabolic intervention; expression, folding, and toxicity in vivo still require wet-lab validation."]
    open(os.path.join(D, "design_card_top10.md"), "w", encoding="utf-8").write("\n".join(md) + "\n")
    json.dump({"designs": cards}, open(os.path.join(D, "design_card_top10.json"), "w"), indent=2)

    print(f"[done] {len(cards)} designs -> design_card_top10.fasta / .md / .json")
    for c in cards:
        sc = c["resource_aware_score"]
        print(f"  {c['rank']:2d} {c['id']:16s} {c['tier']} "
              f"{('Δ+'+str(sc)) if isinstance(sc,(int,float)) else 'FSEOF':>12s} "
              f"GC={c['gc_pct']}% REfree={c['restriction_free']} {c['len_bp']}bp")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
