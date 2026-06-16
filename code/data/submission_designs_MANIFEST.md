# submission_designs.fasta - MANIFEST

Generated: 2026-06-16

**Total unique records: 14**  
(Section 1 primary curated: 10 | Section 2 supplementary generator: 4)  
**Duplicates removed (exact DNA sequence match): 8**  
Total input design records read: 22

**Dedup rule:** exact DNA sequence string. Section 1 takes precedence over Section 2; within Section 2, first-seen source order (loop > proteome > evo2) wins for identical sequences. Reference genome/proteome inputs were not touched.

**File layout:** Section 1 (records 1-10) then Section 2 (records 11-14), in order. No comment lines - file is a valid plain FASTA. Section 2 headers carry an appended `|provenance=` tag.

| # | ID | Gene | Intervention | Product / Rxn | Tier / Provenance | Length (bp) | Score | Section |
|---|----|------|--------------|---------------|-------------------|-------------|-------|---------|
| 1 | lysC_overexpr | lysC | overexpression | L-lysine | T1 | 1350 | Δ=0.864141 | 1-primary |
| 2 | dapA_A49D | dapA | fbr-relaxation | L-lysine | T1 | 879 | Δ=0.531444 | 1-primary |
| 3 | dapA_N80D | dapA | fbr-relaxation | L-lysine | T1 | 879 | Δ=0.531444 | 1-primary |
| 4 | dapA_E84A | dapA | fbr-relaxation | L-lysine | T1 | 879 | Δ=0.531444 | 1-primary |
| 5 | dapA_overexpr | dapA | overexpression | L-lysine | T1 | 879 | Δ=0.242553 | 1-primary |
| 6 | dapB_overexpr | dapB | overexpression | L-lysine | T2 | 822 | - | 1-primary |
| 7 | dapE_overexpr | dapE | overexpression | L-lysine | T2 | 1128 | - | 1-primary |
| 8 | dapF_overexpr | dapF | overexpression | L-lysine | T2 | 825 | - | 1-primary |
| 9 | rpe_overexpr | rpe | overexpression | L-lysine (NADPH supply) | T2 | 678 | - | 1-primary |
| 10 | frdA_overexpr | frdA | overexpression | succinate | T2 | 1809 | - | 1-primary |
| 11 | rank4 (lysC_M318I) | lysC | overexpression | rxn:ASPK | loop | 1350 | modalΔ=0.0 | 2-supplementary |
| 12 | rank5 (lysC_S321F) | lysC | overexpression | rxn:ASPK | loop | 1350 | modalΔ=0.0 | 2-supplementary |
| 13 | evo2_dapA_seed42 | dapA | generator-denovo | - | evo2 | 300 | modalD=+0.53144 | 2-supplementary |
| 14 | evo2_lysC_seed42 | lysC | generator-denovo | - | evo2 | 300 | modalD=+0.00000 | 2-supplementary |

## Duplicates removed (sequences already present, not re-included)
| Source | Record | Reason |
|--------|--------|--------|
| loop_ranked | rank1 dapA_A49D | identical to Section 1 dapA_A49D |
| loop_ranked | rank2 dapA_N80D | identical to Section 1 dapA_N80D |
| loop_ranked | rank3 dapA_E84A | identical to Section 1 dapA_E84A |
| proteome_gen_ranked | rank1 dapA_A49D | identical to Section 1 dapA_A49D |
| proteome_gen_ranked | rank2 dapA_N80D | identical to Section 1 dapA_N80D |
| proteome_gen_ranked | rank3 dapA_E84A | identical to Section 1 dapA_E84A |
| proteome_gen_ranked | rank4 lysC_M318I | identical to loop_ranked rank4 (kept as provenance=loop) |
| proteome_gen_ranked | rank5 lysC_S321F | identical to loop_ranked rank5 (kept as provenance=loop) |

Total duplicates removed: 8
