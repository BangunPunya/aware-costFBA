#!/usr/bin/env python3
from __future__ import annotations

import random
from collections import OrderedDict

# Codon usage table: E. coli K-12 (Kazusa, taxid 83333), frequency per 1000.
# Source: https://www.kazusa.or.jp/codon/ ("Escherichia coli K12").
# Each amino acid maps to its synonymous codons sorted (descending) by
# frequency. The first codon is the optimal ('max') choice.
# HEG = highly-expressed-gene preferred codon (Sharp & Li 1987) noted inline.
CODON_USAGE = {
    'A': [('GCG', 33.7), ('GCC', 25.5), ('GCA', 20.2), ('GCT', 15.3)],          # Ala; HEG GCT/GCA also strong
    'R': [('CGC', 22.0), ('CGT', 20.9), ('CGG', 5.4), ('CGA', 3.6),
          ('AGA', 2.1), ('AGG', 1.2)],                                          # Arg; HEG CGT
    'N': [('AAC', 24.4), ('AAT', 17.7)],                                        # Asn; HEG AAC
    'D': [('GAT', 32.1), ('GAC', 19.1)],                                        # Asp; HEG GAT
    'C': [('TGC', 6.5), ('TGT', 5.2)],                                          # Cys; HEG TGC
    'Q': [('CAG', 28.8), ('CAA', 15.3)],                                        # Gln; HEG CAG
    'E': [('GAA', 39.4), ('GAG', 17.8)],                                        # Glu; HEG GAA
    'G': [('GGC', 29.6), ('GGT', 24.7), ('GGG', 11.1), ('GGA', 8.0)],           # Gly; HEG GGT/GGC
    'H': [('CAT', 12.9), ('CAC', 9.7)],                                         # His; HEG CAT/CAC
    'I': [('ATT', 30.3), ('ATC', 24.4), ('ATA', 4.4)],                          # Ile; HEG ATC
    'L': [('CTG', 52.6), ('TTA', 13.9), ('TTG', 13.7), ('CTC', 11.1),
          ('CTT', 11.0), ('CTA', 3.9)],                                         # Leu; HEG CTG
    'K': [('AAA', 33.6), ('AAG', 10.3)],                                        # Lys; HEG AAA
    'M': [('ATG', 27.9)],                                                       # Met
    'F': [('TTT', 22.1), ('TTC', 16.6)],                                        # Phe; HEG TTC
    'P': [('CCG', 23.2), ('CCA', 8.4), ('CCT', 7.0), ('CCC', 5.5)],             # Pro; HEG CCG
    'S': [('AGC', 16.1), ('TCT', 8.5), ('TCC', 8.6), ('AGT', 8.8),
          ('TCG', 8.9), ('TCA', 7.2)],                                          # Ser; HEG TCT/AGC
    'T': [('ACC', 23.4), ('ACG', 14.4), ('ACT', 9.0), ('ACA', 7.1)],           # Thr; HEG ACC
    'W': [('TGG', 15.2)],                                                       # Trp
    'Y': [('TAT', 16.2), ('TAC', 12.2)],                                        # Tyr; HEG TAC
    'V': [('GTG', 26.4), ('GTT', 18.3), ('GTC', 15.3), ('GTA', 10.9)],          # Val; HEG GTG/GTT
    '*': [('TAA', 2.0), ('TGA', 0.9), ('TAG', 0.2)],                            # Stop; TAA preferred in E. coli
}

# Sort each amino acid's codon list by descending frequency (defensive).
for _aa in CODON_USAGE:
    CODON_USAGE[_aa] = sorted(CODON_USAGE[_aa], key=lambda c: c[1], reverse=True)

# Standard genetic code for translation (sense codons identical to NCBI table 11)
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

# Banned restriction sites (palindromic / common cloning enzymes).
# These are scanned on the forward strand; all listed are palindromes except
# NotI which is also a palindrome, so forward-strand scan is sufficient.
BANNED_SITES = OrderedDict([
    ('EcoRI', 'GAATTC'),
    ('BamHI', 'GGATCC'),
    ('HindIII', 'AAGCTT'),
    ('XhoI', 'CTCGAG'),
    ('NdeI', 'CATATG'),
    ('NotI', 'GCGGCCGC'),
])

# Strong consensus Shine-Dalgarno / RBS (optional). Canonical strong E. coli RBS:
# AGGAGG SD core with 8 nt spacer to ATG (Anderson/iGEM BBa_B0034 lineage:
# AAAGAGGAGAAA + spacer). Documented so users know what gets prepended.
DEFAULT_RBS = 'AAAGAGGAGAAATACTAG'  # SD core AGGAGG, with spacer before ATG
# NOTE: when add_rbs=True the RBS is prepended *before* the start ATG; the
# spacer length here places the SD ~7-8 nt upstream of the ATG. Users targeting
# a specific vector should supply their own rbs_seq.

_DNA_BASES = set('ACGT')


# Translation (round-trip verification)
def translate(dna: str) -> str:
    """Translate DNA coding sequence to protein (NCBI table 11 sense codons). Stops render as '*'. Trailing partial codons ignored. Raises on unknown codons."""
    dna = dna.upper().replace('U', 'T')
    protein = []
    for i in range(0, len(dna) - len(dna) % 3, 3):
        codon = dna[i:i + 3]
        if codon not in CODON_TABLE:
            raise ValueError(f"Unknown codon '{codon}' at position {i}")
        protein.append(CODON_TABLE[codon])
    return ''.join(protein)


# Sequence-hygiene helpers
def gc_content(dna: str) -> float:
    """Fraction of G/C bases (0..1). Empty string -> 0.0."""
    if not dna:
        return 0.0
    dna = dna.upper()
    gc = sum(1 for b in dna if b in 'GC')
    return gc / len(dna)


def find_restriction_sites(dna: str, sites=BANNED_SITES) -> list:
    """Return list of (enzyme, motif, position) for every banned-site occurrence on the forward strand."""
    dna = dna.upper()
    found = []
    for enzyme, motif in sites.items():
        start = 0
        while True:
            idx = dna.find(motif, start)
            if idx == -1:
                break
            found.append((enzyme, motif, idx))
            start = idx + 1
    return found


def max_homopolymer_run(dna: str):
    """Return (base, run_length) of the longest single-base run. ('', 0) if empty."""
    if not dna:
        return ('', 0)
    dna = dna.upper()
    best_base, best_len = dna[0], 1
    cur_base, cur_len = dna[0], 1
    for b in dna[1:]:
        if b == cur_base:
            cur_len += 1
        else:
            cur_base, cur_len = b, 1
        if cur_len > best_len:
            best_base, best_len = cur_base, cur_len
    return (best_base, best_len)


# Core codon selection
def _ordered_codons(aa: str, strategy: str, rng: random.Random):
    """Return synonymous codons for aa ordered by preference. 'max' -> descending frequency (deterministic); 'weighted' -> frequency-weighted random shuffle (top = sampled choice, rest = fallback)."""
    table = CODON_USAGE[aa]
    if strategy == 'max':
        return [c for c, _ in table]
    elif strategy == 'weighted':
        # Frequency-weighted sampling without replacement to get an ordering.
        remaining = list(table)
        ordered = []
        while remaining:
            total = sum(f for _, f in remaining)
            r = rng.uniform(0, total)
            acc = 0.0
            for i, (c, f) in enumerate(remaining):
                acc += f
                if r <= acc:
                    ordered.append(c)
                    remaining.pop(i)
                    break
            else:
                ordered.append(remaining.pop()[0])
        return ordered
    else:
        raise ValueError(f"Unknown strategy '{strategy}' (use 'max' or 'weighted')")


def _creates_problem(dna_so_far: str, candidate_codon: str, avoid_sites: bool,
                     break_homopolymers: bool, homopolymer_max: int):
    """Check whether appending candidate_codon introduces a banned restriction site or over-long homopolymer run (only the junction needs checking). Returns True if problematic."""
    # Window: longest motif is 8 (NotI). Check tail+codon region.
    window = dna_so_far[-8:] + candidate_codon
    if avoid_sites:
        for motif in BANNED_SITES.values():
            if motif in window:
                return True
    if break_homopolymers:
        # Check homopolymer run within the junction window (tail + codon).
        tail = dna_so_far[-(homopolymer_max):] + candidate_codon
        _, run = max_homopolymer_run(tail)
        if run >= homopolymer_max:
            return True
    return False


# Main entry point
def reverse_translate(protein_seq: str,
                      strategy: str = 'max',
                      add_start: bool = True,
                      add_stop: bool = True,
                      add_rbs: bool = False,
                      rbs_seq: str = DEFAULT_RBS,
                      avoid_sites: bool = True,
                      break_homopolymers: bool = True,
                      homopolymer_max: int = 8,
                      gc_low: float = 0.45,
                      gc_high: float = 0.60,
                      seed=None) -> dict:
    """Reverse-translate a protein into codon-optimized, lab-orderable E. coli K-12 DNA.

    strategy {'max','weighted'}: 'max' uses most frequent codon; 'weighted' samples
    by usage freq. add_start: prepend ATG if protein lacks M. add_stop: append TAA
    (E. coli preferred) unless trailing '*'. add_rbs/rbs_seq: prepend Shine-Dalgarno.
    avoid_sites/break_homopolymers/homopolymer_max: swap synonyms to avoid BANNED_SITES
    and keep runs < max (default 8). gc_low/gc_high: GC warning window (0.45-0.60).
    seed: RNG for 'weighted'. Trailing '*' tolerated/stripped; case-insensitive.

    Returns dict: dna, cds, protein_in, gc, gc_in_range, banned_sites,
    unavoidable_sites, homopolymer_max_run, strategy, added_start/stop/rbs.
    """
    rng = random.Random(seed)

    prot = protein_seq.strip().upper().replace(' ', '').replace('\n', '')
    ended_in_stop = prot.endswith('*')
    if ended_in_stop:
        prot = prot[:-1]

    # Validate residues.
    for aa in prot:
        if aa not in CODON_USAGE or aa == '*':
            raise ValueError(f"Unsupported residue '{aa}' in protein sequence "
                             "(ambiguous codes B/Z/X/J/U/O and internal '*' "
                             "are not allowed)")

    added_start = False
    encode = prot
    if add_start and (not prot or prot[0] != 'M'):
        encode = 'M' + prot
        added_start = True

    # Build CDS codon-by-codon with constraint-aware swapping.
    #
    # Many banned sites straddle a codon boundary (e.g. NdeI CAT|ATG = His-Met,
    # where Met has only the ATG codon). A purely forward greedy pass cannot fix
    # those by swapping the *current* codon, so we allow single-codon backtracking:
    # if no synonym of the current residue avoids the problem, we re-pick the
    # *previous* codon with its next-best synonym and retry. Only if the previous
    # residue is also exhausted do we declare the site genuinely unavoidable.
    cds = []
    aa_list = list(encode)
    # candidate ordering per position (computed once; deterministic for 'max',
    # frequency-sampled once for 'weighted')
    cand_lists = [_ordered_codons(aa, strategy, rng) for aa in aa_list]
    cand_idx = [0] * len(aa_list)           # which synonym we're currently on
    unavoidable = []

    def prefix(upto):
        return ''.join(cds[:upto])

    i = 0
    while i < len(aa_list):
        dna_str = ''.join(cds[:i])
        candidates = cand_lists[i]
        chosen = None
        # try synonyms starting at the current index for this position
        for j in range(cand_idx[i], len(candidates)):
            codon = candidates[j]
            if not _creates_problem(dna_str, codon, avoid_sites,
                                    break_homopolymers, homopolymer_max):
                chosen = codon
                cand_idx[i] = j
                break
        if chosen is not None:
            if i < len(cds):
                cds[i] = chosen
            else:
                cds.append(chosen)
            i += 1
            continue

        # No synonym at this position works. Try to backtrack to position i-1
        # and advance it to its next synonym (single-codon lookback).
        backtracked = False
        if i > 0 and cand_idx[i - 1] + 1 < len(cand_lists[i - 1]):
            cand_idx[i - 1] += 1
            cand_idx[i] = 0          # reset current position to retry from top
            cds.pop()                # remove previous codon; it will be re-chosen
            i -= 1
            backtracked = True
        if not backtracked:
            # Genuinely unavoidable: take the top synonym, record the site(s).
            chosen = candidates[cand_idx[i]] if cand_idx[i] < len(candidates) else candidates[0]
            dna_str = ''.join(cds[:i])
            for enzyme, motif in BANNED_SITES.items():
                if motif in (dna_str[-8:] + chosen) and motif not in dna_str[-8:]:
                    unavoidable.append((enzyme, motif, len(dna_str)))
            if i < len(cds):
                cds[i] = chosen
            else:
                cds.append(chosen)
            i += 1

    dna_str = ''.join(cds)

    # Append stop codon.
    added_stop = False
    if add_stop and not ended_in_stop:
        stop_candidates = _ordered_codons('*', 'max', rng)  # TAA first
        stop_chosen = None
        for codon in stop_candidates:
            if not _creates_problem(dna_str, codon, avoid_sites,
                                    break_homopolymers, homopolymer_max):
                stop_chosen = codon
                break
        if stop_chosen is None:
            stop_chosen = stop_candidates[0]
        cds.append(stop_chosen)
        dna_str += stop_chosen
        added_stop = True

    cds_seq = ''.join(cds)

    # Prepend RBS.
    added_rbs = False
    final_dna = cds_seq
    if add_rbs:
        final_dna = rbs_seq.upper() + cds_seq
        added_rbs = True

    banned = find_restriction_sites(final_dna)
    gc = gc_content(final_dna)
    homo = max_homopolymer_run(final_dna)

    return {
        'dna': final_dna,
        'cds': cds_seq,
        'protein_in': prot,
        'gc': gc,
        'gc_in_range': gc_low <= gc <= gc_high,
        'banned_sites': banned,
        'unavoidable_sites': unavoidable,
        'homopolymer_max_run': homo,
        'strategy': strategy,
        'added_start': added_start,
        'added_stop': added_stop,
        'added_rbs': added_rbs,
    }


if __name__ == '__main__':
    demo = "MSEIVK"
    res = reverse_translate(demo)
    print("protein:", demo)
    print("dna    :", res['dna'])
    print("back   :", translate(res['cds']))
    print("gc     :", round(res['gc'], 3), "in_range:", res['gc_in_range'])
    print("banned :", res['banned_sites'])
    print("homo   :", res['homopolymer_max_run'])
