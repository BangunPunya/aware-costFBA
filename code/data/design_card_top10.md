# Top-10 lab-ready design card

> Resource-aware oracle vets the INTERVENTION (medium/ΔG/enzyme-capacity, growth-coupled); sequence gates vet the SEQUENCE (codon-optimized, restriction-free). NOT a titer prediction; a prioritized, metabolically-feasible shortlist. T1 = oracle Δ-scored; T2 = genome-scale FSEOF-promoted target.

| # | Design | Gene | Intervention | Product | Tier | Resource-aware score | Growth | M9 | Fold gate | GC% | RE-free | bp |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | lysC_overexpr | lysC | overexpression | L-lysine | T1 | +0.864141 | yes | pass | WT fold (pass) | 60.5 | False | 1350 |
| 2 | dapA_A49D | dapA | fbr-relaxation | L-lysine | T1 | +0.531444 | yes | pass | ESM2 fold gate pass (ACT-domain fbr) | 58.4 | True | 879 |
| 3 | dapA_N80D | dapA | fbr-relaxation | L-lysine | T1 | +0.531444 | yes | pass | ESM2 fold gate pass (ACT-domain fbr) | 58.6 | True | 879 |
| 4 | dapA_E84A | dapA | fbr-relaxation | L-lysine | T1 | +0.531444 | yes | pass | ESM2 fold gate pass (ACT-domain fbr) | 58.8 | True | 879 |
| 5 | dapA_overexpr | dapA | overexpression | L-lysine | T1 | +0.242553 | yes | pass | WT fold (pass) | 58.6 | True | 879 |
| 6 | dapB_overexpr | dapB | overexpression | L-lysine | T2 | FSEOF-promoted | yes | pass | WT fold (pass) | 59.6 | True | 822 |
| 7 | dapE_overexpr | dapE | overexpression | L-lysine | T2 | FSEOF-promoted | yes | pass | WT fold (pass) | 60.0 | False | 1128 |
| 8 | dapF_overexpr | dapF | overexpression | L-lysine | T2 | FSEOF-promoted | yes | pass | WT fold (pass) | 59.8 | True | 825 |
| 9 | rpe_overexpr | rpe | overexpression | L-lysine (NADPH supply) | T2 | FSEOF-promoted | yes | pass | WT fold (pass) | 54.3 | False | 678 |
| 10 | frdA_overexpr | frdA | overexpression | succinate | T2 | FSEOF-promoted | yes | pass | WT fold (pass) | 59.5 | False | 1809 |

Score types:

- **lysC_overexpr**: ec-FBA Δ lysine secretion @0.5 growth (eciML1515)
- **dapA_A49D**: resource-aware Δ lysine flux, near-max growth
- **dapA_N80D**: resource-aware Δ lysine flux, near-max growth
- **dapA_E84A**: resource-aware Δ lysine flux, near-max growth
- **dapA_overexpr**: ec-FBA Δ lysine secretion @0.5 growth (eciML1515)
- **dapB_overexpr**: genome-scale FSEOF, resource-aware-promoted (DHDPRy)
- **dapE_overexpr**: genome-scale FSEOF, resource-aware-promoted (SDPDS)
- **dapF_overexpr**: genome-scale FSEOF, resource-aware-promoted (DAPE)
- **rpe_overexpr**: genome-scale FSEOF non-obvious, resource-aware-promoted (RPE/PPP)
- **frdA_overexpr**: genome-scale FSEOF, resource-aware-promoted (FRD2)

FASTA: `code/data/design_card_top10.fasta`. Sequences: T1 fbr constructs from the loop run; overexpression and T2 targets are gene-product proteins reverse-translated to codon-optimal DNA (HEG-max, seed 42, GC 30-70 percent, homopolymer broken).

Restriction sites: 6/10 constructs are already free of the eight common cloning sites (EcoRI, BamHI, HindIII, PstI, SalI, XbaI, SacI, NdeI). The remaining longer constructs carry one or more internal sites that standard commercial gene synthesis removes by synonymous substitution at no extra cost, so all ten remain orderable.

Limitations: scores are in-silico (FBA upper-bound and growth-coupling), not predicted titer. T2 FSEOF targets are ranked amplification targets, not Δ-scored on secretion. The oracle vets the metabolic intervention; expression, folding, and toxicity in vivo still require wet-lab validation.
