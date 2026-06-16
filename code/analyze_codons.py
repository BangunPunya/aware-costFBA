#!/usr/bin/env python3
import os
import sys
import json
import time
import statistics
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
GENOME_FASTA = os.path.join(DATA, "ecoli_k12_genome.fasta")
REPORT_JSON = os.path.join(DATA, "codon_analysis.json")
ACC = "NC_000913.3"

START_CODONS = ["ATG", "GTG", "TTG"]
STOP_CODONS = ["TAA", "TAG", "TGA"]

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]


def fetch_genome():
    """Download NC_000913.3 nucleotide FASTA via NCBI EUtils efetch (requests)."""
    if os.path.exists(GENOME_FASTA) and os.path.getsize(GENOME_FASTA) > 4_000_000:
        with open(GENOME_FASTA) as fh:
            return fh.read()
    import requests
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = dict(db="nuccore", id=ACC, rettype="fasta", retmode="text")
    for attempt in range(4):
        r = requests.get(url, params=params, timeout=120)
        if r.status_code == 200 and r.text.startswith(">"):
            with open(GENOME_FASTA, "w") as fh:
                fh.write(r.text)
            return r.text
        time.sleep(3 * (attempt + 1))
    raise RuntimeError(f"genome fetch failed: HTTP {r.status_code}")


def parse_fasta_single(txt):
    lines = txt.splitlines()
    header = lines[0]
    seq = "".join(l.strip() for l in lines[1:] if not l.startswith(">"))
    return header, seq.upper()


def fetch_official_cds_count():
    """Annotated protein-coding CDS count for NC_000913.3.

    Parses GenBank feature table (rettype=gb); counts CDS features with a
    /translation qualifier and gathers protein aa-lengths.
    """
    import requests
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = dict(db="nuccore", id=ACC, rettype="gbwithparts", retmode="text")
    gb_path = os.path.join(DATA, "ecoli_k12_genome.gb")
    txt = None
    if os.path.exists(gb_path) and os.path.getsize(gb_path) > 5_000_000:
        with open(gb_path) as fh:
            txt = fh.read()
    if txt is None:
        for attempt in range(4):
            r = requests.get(url, params=params, timeout=300)
            if r.status_code == 200 and "LOCUS" in r.text[:200]:
                txt = r.text
                with open(gb_path, "w") as fh:
                    fh.write(txt)
                break
            time.sleep(4 * (attempt + 1))
        if txt is None:
            return dict(error=f"gb fetch failed HTTP {r.status_code}")

    # Parse the FEATURES table. Count CDS features and gather /translation lengths.
    lines = txt.splitlines()
    in_features = False
    cds_count = 0
    pseudo_in_current = False
    cur_is_cds = False
    translations = []
    collecting_trans = False
    trans_buf = []
    prot_lengths = []

    def flush_translation():
        nonlocal trans_buf
        if trans_buf:
            seq = "".join(trans_buf).replace('"', "").replace(" ", "")
            if seq:
                prot_lengths.append(len(seq))
        trans_buf = []

    for line in lines:
        if line.startswith("FEATURES"):
            in_features = True
            continue
        if in_features and line.startswith("ORIGIN"):
            flush_translation()
            break
        if not in_features:
            continue
        # Feature key lines start at col 5 (5 spaces) with key at col 6
        if len(line) > 5 and line[5] != " " and not line.startswith(" " * 21):
            # new feature
            if collecting_trans:
                flush_translation()
                collecting_trans = False
            key = line[5:21].strip()
            if key == "CDS":
                cds_count += 1
                cur_is_cds = True
            else:
                cur_is_cds = False
        else:
            # qualifier / continuation line (indented ~21 spaces)
            stripped = line.strip()
            if collecting_trans:
                if stripped.startswith("/"):
                    flush_translation()
                    collecting_trans = False
                else:
                    trans_buf.append(stripped.replace('"', ""))
                    continue
            if cur_is_cds and stripped.startswith("/translation="):
                collecting_trans = True
                trans_buf = [stripped.split("=", 1)[1].replace('"', "")]
    flush_translation()

    return dict(
        cds_feature_count=cds_count,
        cds_with_translation=len(prot_lengths),
        protein_lengths_aa=prot_lengths,
    )


def count_codons_anywhere(seq, codons):
    """Count overlapping-free, non-frame-restricted occurrences of each codon."""
    out = {}
    for c in codons:
        # count non-overlapping occurrences across the whole string
        out[c] = seq.count(c)
    return out


def count_codons_framed(seq, codons):
    """Count codons restricted to each of 3 frames (frame 0,1,2)."""
    res = {}
    n = len(seq)
    for frame in range(3):
        cnt = Counter()
        for i in range(frame, n - 2, 3):
            cod = seq[i:i + 3]
            if cod in codons:
                cnt[cod] += 1
        res[f"frame{frame}"] = dict(cnt)
    return res


def enumerate_orfs(seq, min_aa_thresholds):
    """Enumerate naive start->in-frame-stop ORFs across all 3 frames of `seq`.

    ORF = start codon followed (in frame) by first stop codon. Length in aa =
    (stop_pos - start_pos)/3. Returns (all_orf_aa, first_start_orf_aa): the first
    counts every start->stop (over-counts nested starts); the second keeps one ORF
    per stop-segment (first start after previous stop).
    """
    start_set = set(START_CODONS)
    stop_set = set(STOP_CODONS)
    n = len(seq)
    all_orf_aa = []          # every start->stop (naive, over-counts nested)
    first_start_orf_aa = []  # one per stop region (first start after prev stop)

    for frame in range(3):
        last_stop = frame - 3  # position just before frame start
        seen_start_in_segment = None
        i = frame
        while i <= n - 3:
            cod = seq[i:i + 3]
            if cod in stop_set:
                last_stop = i
                seen_start_in_segment = None
                i += 3
                continue
            if cod in start_set:
                # find next in-frame stop
                j = i + 3
                while j <= n - 3:
                    c2 = seq[j:j + 3]
                    if c2 in stop_set:
                        aa = (j - i) // 3  # codons from start up to (not incl) stop
                        all_orf_aa.append(aa)
                        if seen_start_in_segment is None:
                            first_start_orf_aa.append(aa)
                            seen_start_in_segment = i
                        break
                    j += 3
            i += 3
    return all_orf_aa, first_start_orf_aa


def threshold_counts(aa_list, thresholds):
    return {f">={t}aa": sum(1 for a in aa_list if a >= t) for t in thresholds}


def run_pyrodigal(seq):
    import pyrodigal
    # single/normal mode, trained on this genome (not meta)
    orf_finder = pyrodigal.GeneFinder(meta=False)
    orf_finder.train(seq)
    genes = orf_finder.find_genes(seq)
    lengths_aa = []
    for g in genes:
        prot = g.translate()
        # translate() returns protein incl trailing '*'
        L = len(prot.rstrip("*"))
        lengths_aa.append(L)
    return len(genes), lengths_aa


def prot_stats(lengths):
    if not lengths:
        return {}
    return dict(
        n=len(lengths),
        min_aa=min(lengths),
        median_aa=statistics.median(lengths),
        mean_aa=round(statistics.mean(lengths), 1),
        max_aa=max(lengths),
        count_lt_100aa=sum(1 for x in lengths if x < 100),
        count_lt_50aa=sum(1 for x in lengths if x < 50),
    )


def main():
    print("Fetching genome...", flush=True)
    txt = fetch_genome()
    header, seq = parse_fasta_single(txt)
    n = len(seq)
    gc = (seq.count("G") + seq.count("C")) / n * 100.0
    rc = revcomp(seq)
    print(f"  genome {ACC}: {n:,} bp, GC {gc:.2f}%", flush=True)

    thresholds = [30, 60, 100, 150, 200, 250]

    # 2. Raw codon counts.
    start_fwd_any = count_codons_anywhere(seq, START_CODONS)
    stop_fwd_any = count_codons_anywhere(seq, STOP_CODONS)
    start_rc_any = count_codons_anywhere(rc, START_CODONS)
    stop_rc_any = count_codons_anywhere(rc, STOP_CODONS)
    start_both_any = {c: start_fwd_any[c] + start_rc_any[c] for c in START_CODONS}
    stop_both_any = {c: stop_fwd_any[c] + stop_rc_any[c] for c in STOP_CODONS}

    start_fwd_framed = count_codons_framed(seq, set(START_CODONS))
    stop_fwd_framed = count_codons_framed(seq, set(STOP_CODONS))
    start_rc_framed = count_codons_framed(rc, set(START_CODONS))
    stop_rc_framed = count_codons_framed(rc, set(STOP_CODONS))

    print("Enumerating naive ORFs (6 frames)...", flush=True)
    fwd_all, fwd_first = enumerate_orfs(seq, thresholds)
    rc_all, rc_first = enumerate_orfs(rc, thresholds)
    naive_all_aa = fwd_all + rc_all          # every start->stop, all 6 frames
    naive_first_aa = fwd_first + rc_first    # one per stop-segment, all 6 frames

    naive_all_thr = threshold_counts(naive_all_aa, thresholds)
    naive_first_thr = threshold_counts(naive_first_aa, thresholds)

    print("Running pyrodigal (single/normal mode, trained)...", flush=True)
    pyro_n, pyro_lengths = run_pyrodigal(seq)
    print(f"  pyrodigal predicted {pyro_n} genes", flush=True)

    print("Fetching official NCBI CDS annotation...", flush=True)
    official = fetch_official_cds_count()
    off_cds = official.get("cds_with_translation") or official.get("cds_feature_count")
    off_lengths = official.get("protein_lengths_aa", [])
    print(f"  official CDS (with translation): {off_cds}", flush=True)

    # protein length stats.
    pyro_stats = prot_stats(pyro_lengths)
    off_stats = prot_stats(off_lengths) if off_lengths else {}

    # reference real gene count for ratios
    ref = off_cds if off_cds else pyro_n
    median_aa = off_stats.get("median_aa") or pyro_stats.get("median_aa")
    median_bp = int(median_aa * 3) if median_aa else None

    # over-counting ratios
    overcount = {}
    for t in thresholds:
        k = f">={t}aa"
        overcount[k] = dict(
            naive_all_orfs=naive_all_thr[k],
            naive_first_start_orfs=naive_first_thr[k],
            ratio_naiveAll_over_official=round(naive_all_thr[k] / ref, 2),
            ratio_naiveFirst_over_official=round(naive_first_thr[k] / ref, 2),
        )

    report = dict(
        genome=dict(accession=ACC, header=header, length_bp=n, gc_percent=round(gc, 3)),
        raw_codon_counts=dict(
            start_codons_forward_anywhere=start_fwd_any,
            stop_codons_forward_anywhere=stop_fwd_any,
            start_codons_both_strands_anywhere=start_both_any,
            stop_codons_both_strands_anywhere=stop_both_any,
            start_codons_forward_framed=start_fwd_framed,
            stop_codons_forward_framed=stop_fwd_framed,
            start_codons_revcomp_framed=start_rc_framed,
            stop_codons_revcomp_framed=stop_rc_framed,
        ),
        naive_orf_enumeration=dict(
            thresholds_aa=thresholds,
            note="naive_all_orfs = every in-frame start->stop in all 6 frames (counts nested starts). naive_first_start_orfs = one ORF per stop-segment (first start after prev stop), 6 frames.",
            counts_naive_all=naive_all_thr,
            counts_naive_first_start=naive_first_thr,
            total_naive_all_orfs_unfiltered=len(naive_all_aa),
            total_naive_first_orfs_unfiltered=len(naive_first_aa),
        ),
        pyrodigal=dict(
            mode="single/normal (meta=False), self-trained on NC_000913.3",
            gene_count=pyro_n,
            protein_length_stats=pyro_stats,
        ),
        official_ncbi=dict(
            accession=ACC,
            cds_feature_count=official.get("cds_feature_count"),
            cds_with_translation=official.get("cds_with_translation"),
            protein_length_stats=off_stats,
            error=official.get("error"),
        ),
        comparison=dict(
            reference_real_gene_count=ref,
            reference_source="official NCBI CDS-with-translation" if off_cds else "pyrodigal",
            median_real_protein_aa=median_aa,
            median_real_protein_bp=median_bp,
            overcounting_by_threshold=overcount,
        ),
    )

    with open(REPORT_JSON, "w") as fh:
        json.dump(report, fh, indent=2)

    # pretty print.
    print("\n" + "=" * 78)
    print("E. coli K-12 MG1655 (NC_000913.3) - naive codon/ORF vs real genes")
    print("=" * 78)
    print(f"Genome length : {n:,} bp   GC: {gc:.2f}%")
    print(f"\nRaw START codons (forward strand, anywhere): {start_fwd_any}")
    print(f"Raw STOP  codons (forward strand, anywhere): {stop_fwd_any}")
    print(f"Raw START codons (both strands, anywhere)  : {start_both_any}")
    print(f"Raw STOP  codons (both strands, anywhere)  : {stop_both_any}")
    tot_start_both = sum(start_both_any.values())
    tot_stop_both = sum(stop_both_any.values())
    print(f"  total start-codon sites (both strands): {tot_start_both:,}")
    print(f"  total stop-codon  sites (both strands): {tot_stop_both:,}")

    print("\n--- COMPARISON TABLE: naive ORF count vs real gene count ---")
    print(f"{'min len':>8} | {'naive ALL':>11} | {'naive 1st':>10} | "
          f"{'pyrodigal':>9} | {'official':>8} | {'ratio all/off':>13}")
    print("-" * 78)
    for t in thresholds:
        k = f">={t}aa"
        print(f"{t:>6}aa | {naive_all_thr[k]:>11,} | {naive_first_thr[k]:>10,} | "
              f"{pyro_n:>9,} | {ref:>8,} | "
              f"{overcount[k]['ratio_naiveAll_over_official']:>13}x")
    print("-" * 78)
    print(f"pyrodigal genes: {pyro_n:,}   |   official annotated CDS: {ref:,}")

    print("\n--- PROTEIN LENGTH DISTRIBUTION ---")
    if off_stats:
        print(f"Official NCBI CDS  : n={off_stats['n']}, min={off_stats['min_aa']}, "
              f"median={off_stats['median_aa']}, mean={off_stats['mean_aa']}, "
              f"max={off_stats['max_aa']} aa | <100aa: {off_stats['count_lt_100aa']}, "
              f"<50aa: {off_stats['count_lt_50aa']}")
    print(f"Pyrodigal genes    : n={pyro_stats['n']}, min={pyro_stats['min_aa']}, "
          f"median={pyro_stats['median_aa']}, mean={pyro_stats['mean_aa']}, "
          f"max={pyro_stats['max_aa']} aa | <100aa: {pyro_stats['count_lt_100aa']}, "
          f"<50aa: {pyro_stats['count_lt_50aa']}")
    print(f"\nMedian real protein length: {median_aa} aa  (~{median_bp} bp)")
    print(f"\nReport written: {REPORT_JSON}")
    return report


if __name__ == "__main__":
    main()
