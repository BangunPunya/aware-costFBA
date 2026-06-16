from __future__ import annotations

import json
import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field

import cobra

CODE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(CODE_DIR, "data")
MODEL_CACHE = os.path.join(CODE_DIR, "iML1515.json")

PROTEOME_FASTA = os.path.join(DATA_DIR, "ecoli_k12_proteome.fasta")
REF_FASTA = os.path.join(DATA_DIR, "iml1515_gene_ref.fasta")
REF_DMND = os.path.join(DATA_DIR, "iml1515_gene_ref.dmnd")
GENEMAP_JSON = os.path.join(DATA_DIR, "iml1515_genemap.json")

# Lokasi binary DIAMOND (cek beberapa lokasi umum instalasi)
_DIAMOND_CANDIDATES = [
    shutil.which("diamond"),
    os.path.expanduser("~/.local/bin/diamond"),
    "/usr/bin/diamond",
    "/usr/local/bin/diamond",
]

# Genetic code 11 (bacterial) - for fallback translate if needed
_CODON = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L",
    "CTA": "L", "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S",
    "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A",
    "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R",
    "CGA": "R", "CGG": "R", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# Ambang kepercayaan (DIAMOND %id), kalibrasi konservatif.
ID_HIGH = 90.0     # >= -> high confidence (homolog dekat / native)
ID_MED = 50.0      # >= -> medium
ID_LOW = 30.0      # >= -> low (twilight zone homologi); < ID_LOW -> di-flag/uncertain
EVALUE_MAX = 1e-5  # hit dengan e-value > ini di-flag tidak-signifikan


# dataclasses
@dataclass
class ORF:
    index: int
    protein: str
    dna: str
    begin: int
    end: int
    strand: int
    partial: bool
    confidence: float           # 0..1 (kelengkapan ORF: lengkap=1.0, parsial<1)
    note: str = ""


@dataclass
class Assignment:
    gene: str | None = None         # b-number gen target (subject DIAMOND)
    gene_name: str | None = None
    ec: list[str] = field(default_factory=list)
    identity: float = 0.0           # % identitas DIAMOND
    evalue: float = 1.0
    bitscore: float = 0.0
    aln_cov: float = 0.0            # cakupan query (qcovhsp/100)
    method: str = ""               # "diamond-blastp" | "fallback-python-pairwise"
    confidence: float = 0.0        # 0..1 turunan identity/evalue/cov
    confidence_label: str = "none"  # high|medium|low|none
    flag: str = ""                 # peringatan bila rendah/tidak-signifikan

    def as_dict(self) -> dict:
        return dict(self.__dict__)


@dataclass
class ReactionHit:
    reaction_id: str
    via: str                       # "gene:bNNNN" | "ec:x.x.x.x"
    confidence: float
    note: str = ""


@dataclass
class AnnotationResult:
    orfs: list[ORF] = field(default_factory=list)
    assignments: list[Assignment] = field(default_factory=list)
    reactions: list[ReactionHit] = field(default_factory=list)
    confidence: float = 0.0        # confidence kaskade gabungan (min rantai terbaik)
    method: str = ""
    low_confidence: bool = False
    notes: list[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        return {
            "orfs": [vars(o) for o in self.orfs],
            "assignments": [a.as_dict() for a in self.assignments],
            "reactions": [vars(r) for r in self.reactions],
            "confidence": self.confidence,
            "method": self.method,
            "low_confidence": self.low_confidence,
            "notes": self.notes,
        }


# util
def _diamond_bin() -> str | None:
    for c in _DIAMOND_CANDIDATES:
        if c and os.path.exists(c) and os.access(c, os.X_OK):
            return c
    return None


def _translate(dna: str) -> str:
    dna = dna.upper().replace("U", "T")
    aa = []
    for i in range(0, len(dna) - 2, 3):
        aa.append(_CODON.get(dna[i:i + 3], "X"))
    return "".join(aa)


def load_genemap() -> dict:
    """b-number -> {name, uniprot, reactions, ec}. Bangun bila belum ada."""
    if not os.path.exists(GENEMAP_JSON):
        build_reference_db()
    with open(GENEMAP_JSON) as fh:
        return json.load(fh)


def _confidence_from_hit(identity: float, evalue: float, cov: float) -> tuple[float, str, str]:
    """Map (identity%, evalue, coverage) -> (confidence 0..1, label, flag); nonsignificant e-value or identity < ID_LOW -> low/none + flag."""
    flag = ""
    if evalue > EVALUE_MAX:
        return 0.05, "none", f"e-value {evalue:.1e} > {EVALUE_MAX:.0e} (tidak signifikan)"
    if identity >= ID_HIGH:
        label = "high"
    elif identity >= ID_MED:
        label = "medium"
    elif identity >= ID_LOW:
        label = "low"
        flag = f"identitas {identity:.1f}% di zona-rendah (>= {ID_LOW}% tapi < {ID_MED}%)"
    else:
        return 0.10, "none", f"identitas {identity:.1f}% < {ID_LOW}% (twilight zone, tak dipercaya)"
    # confidence kontinu: identitas dominan, dipotong cakupan
    base = min(1.0, identity / 100.0)
    conf = round(base * (0.5 + 0.5 * min(1.0, cov)), 4)
    return conf, label, flag


# reference DB builder (idempotent; cache in data/)
def build_reference_db(force: bool = False) -> dict:
    """Build reference FASTA + genemap + (if diamond present) DIAMOND db; ref = iML1515 gene products (b-number), sequences via UniProt UP000000625 accessions so each hit maps to model gene -> reaction (GPR)."""
    import requests  # lokal: hanya saat membangun DB

    os.makedirs(DATA_DIR, exist_ok=True)
    # 1. proteome (cache)
    if force or not os.path.exists(PROTEOME_FASTA) or os.path.getsize(PROTEOME_FASTA) < 100000:
        r = requests.get(
            "https://rest.uniprot.org/uniprotkb/stream",
            params={"query": "proteome:UP000000625", "format": "fasta", "compressed": "false"},
            timeout=180,
        )
        r.raise_for_status()
        with open(PROTEOME_FASTA, "w") as fh:
            fh.write(r.text)
    # parse acc -> seq
    seqs: dict[str, str] = {}
    acc, buf = None, []
    with open(PROTEOME_FASTA) as fh:
        for line in fh:
            if line.startswith(">"):
                if acc:
                    seqs[acc] = "".join(buf)
                parts = line[1:].split("|")
                acc = parts[1] if len(parts) >= 2 else None
                buf = []
            else:
                buf.append(line.strip())
    if acc:
        seqs[acc] = "".join(buf)

    # 2. model genes -> ref fasta + genemap
    model = cobra.io.load_json_model(MODEL_CACHE)
    genemap: dict[str, dict] = {}
    with open(REF_FASTA, "w") as out:
        for g in model.genes:
            uni = g.annotation.get("uniprot")
            ecs = set()
            for r in g.reactions:
                ec = r.annotation.get("ec-code")
                if isinstance(ec, list):
                    ecs.update(ec)
                elif isinstance(ec, str):
                    ecs.add(ec)
            genemap[g.id] = {
                "name": g.name,
                "uniprot": uni,
                "reactions": [r.id for r in g.reactions],
                "ec": sorted(e for e in ecs if e),
            }
            if uni and uni in seqs:
                out.write(f">{g.id}|{g.name}|{uni}\n{seqs[uni]}\n")
    with open(GENEMAP_JSON, "w") as fh:
        json.dump(genemap, fh)

    # 3. DIAMOND db (jika binary ada)
    db = _diamond_bin()
    if db:
        subprocess.run(
            [db, "makedb", "--in", REF_FASTA, "--db", os.path.join(DATA_DIR, "iml1515_gene_ref")],
            check=True, capture_output=True,
        )
    return genemap


# 1. ORF calling - pyrodigal (meta mode utk seq pendek/tunggal)
def orf_call(dna: str, meta: bool = True, min_aa: int = 30) -> list[ORF]:
    """Call ORFs from raw DNA via pyrodigal (meta mode = no training); each ORF -> protein + completeness confidence (complete=1.0, partial lower); ORFs < min_aa dropped."""
    import pyrodigal

    dna = dna.strip().upper().replace("U", "T")
    orf_finder = pyrodigal.GeneFinder(meta=meta)
    genes = orf_finder.find_genes(dna.encode())
    out: list[ORF] = []
    idx = 0
    starts = ("ATG", "GTG", "TTG")
    stops = ("TAA", "TAG", "TGA")
    for gene in genes:
        prot = gene.translate().rstrip("*")
        if len(prot) < min_aa:
            continue
        # Codon ORF (untai +/- sudah ditangani: ambil substring lalu rev-comp bila -)
        seq = dna[gene.begin - 1:gene.end]
        if gene.strand == -1:
            comp = str.maketrans("ACGTN", "TGCAN")
            seq = seq.translate(comp)[::-1]
        has_start = seq[:3] in starts
        has_stop = seq[-3:] in stops if len(seq) >= 3 else False
        # pyrodigal menandai partial bila ORF menyentuh tepi input (artefak utk CDS
        # terisolasi); kelengkapan ditentukan dari codon start+stop nyata.
        codon_complete = has_start and has_stop
        edge_partial = bool(gene.partial_begin or gene.partial_end)
        if codon_complete:
            conf, note = 1.0, "ORF lengkap (start+stop codon ada)"
        elif edge_partial:
            conf, note = 0.6, "ORF parsial (menyentuh tepi sekuens; tanpa start/stop penuh)"
        else:
            conf, note = 0.85, "ORF tengah-genom tanpa codon batas lengkap"
        out.append(ORF(
            index=idx,
            protein=prot,
            dna=seq,
            begin=gene.begin,
            end=gene.end,
            strand=gene.strand,
            partial=not codon_complete,
            confidence=conf,
            note=note,
        ))
        idx += 1
    return out


# 2. assign_function - DIAMOND blastp vs ref db (fallback: pure-python pairwise)
def assign_function(protein: str, genemap: dict | None = None,
                    max_target_seqs: int = 5) -> Assignment:
    """Assign {gene, ec, identity, evalue} for one protein via DIAMOND blastp vs data/iml1515_gene_ref.dmnd; falls back to pairwise identity vs iML1515 ref seqs (method="fallback-python-pairwise") if DIAMOND/db absent."""
    if genemap is None:
        genemap = load_genemap()
    db = _diamond_bin()
    if db and os.path.exists(REF_DMND):
        return _assign_diamond(protein, genemap, db, max_target_seqs)
    return _assign_fallback(protein, genemap)


def _assign_diamond(protein: str, genemap: dict, db_bin: str,
                    max_target_seqs: int) -> Assignment:
    with tempfile.TemporaryDirectory() as td:
        qf = os.path.join(td, "q.fasta")
        with open(qf, "w") as fh:
            fh.write(f">query\n{protein}\n")
        of = os.path.join(td, "out.tsv")
        cmd = [
            db_bin, "blastp", "--db", os.path.join(DATA_DIR, "iml1515_gene_ref"),
            "--query", qf, "--out", of, "--quiet",
            "--max-target-seqs", str(max_target_seqs),
            "--outfmt", "6", "sseqid", "pident", "evalue", "bitscore", "qcovhsp",
        ]
        try:
            subprocess.run(cmd, check=True, capture_output=True, timeout=120)
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
            a = _assign_fallback(protein, genemap)
            a.flag = (a.flag + " | DIAMOND gagal, jatuh ke fallback").strip(" |")
            return a
        rows = []
        if os.path.exists(of):
            with open(of) as fh:
                for line in fh:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) >= 5:
                        rows.append(parts)
    if not rows:
        return Assignment(method="diamond-blastp", confidence=0.0,
                          confidence_label="none", flag="tak ada hit DIAMOND signifikan")
    # top hit (rows sudah terurut bitscore oleh diamond)
    sseqid, pident, evalue, bitscore, qcov = rows[0][:5]
    bnum = sseqid.split("|")[0]
    info = genemap.get(bnum, {})
    ident = float(pident)
    ev = float(evalue)
    cov = float(qcov) / 100.0
    conf, label, flag = _confidence_from_hit(ident, ev, cov)
    return Assignment(
        gene=bnum, gene_name=info.get("name"), ec=info.get("ec", []),
        identity=ident, evalue=ev, bitscore=float(bitscore), aln_cov=cov,
        method="diamond-blastp", confidence=conf, confidence_label=label, flag=flag,
    )


# fallback: ungapped pairwise identity (k-mer prefilter)
def _load_ref_seqs() -> dict[str, str]:
    if not os.path.exists(REF_FASTA):
        build_reference_db()
    seqs: dict[str, str] = {}
    bnum, buf = None, []
    with open(REF_FASTA) as fh:
        for line in fh:
            if line.startswith(">"):
                if bnum:
                    seqs[bnum] = "".join(buf)
                bnum = line[1:].split("|")[0]
                buf = []
            else:
                buf.append(line.strip())
    if bnum:
        seqs[bnum] = "".join(buf)
    return seqs


def _best_ungapped_identity(q: str, s: str) -> tuple[float, float]:
    """Best ungapped identity over all offsets (q shifted on s); returns (identity%, coverage_query). Enough for round-trip (native gene -> ~identical protein); not Smith-Waterman."""
    if not q or not s:
        return 0.0, 0.0
    best_id, best_cov = 0.0, 0.0
    lq, ls = len(q), len(s)
    # geser q relatif s: offset dari -(lq-1)..(ls-1), batasi rentang utk kecepatan
    for off in range(-(lq - 1), ls):
        matches = aligned = 0
        for i in range(lq):
            j = off + i
            if 0 <= j < ls:
                aligned += 1
                if q[i] == s[j]:
                    matches += 1
        if aligned >= 0.5 * lq and aligned > 0:
            ident = 100.0 * matches / aligned
            if ident > best_id:
                best_id = ident
                best_cov = aligned / lq
    return best_id, best_cov


def _assign_fallback(protein: str, genemap: dict) -> Assignment:
    refs = _load_ref_seqs()
    # prefilter k-mer: kandidat yang berbagi >=1 6-mer dgn query (hemat O(N))
    k = 6
    qk = {protein[i:i + k] for i in range(0, max(0, len(protein) - k + 1))}
    cand = []
    for bnum, s in refs.items():
        sk = {s[i:i + k] for i in range(0, max(0, len(s) - k + 1), 3)}  # subsample
        if qk & sk:
            cand.append((bnum, s))
    if not cand:  # tak ada k-mer cocok -> skor semua (jarang; query sangat divergen)
        cand = list(refs.items())[:200]
    best = None
    for bnum, s in cand:
        ident, cov = _best_ungapped_identity(protein, s)
        if best is None or ident > best[1]:
            best = (bnum, ident, cov)
    bnum, ident, cov = best
    info = genemap.get(bnum, {})
    # e-value tak terdefinisi di fallback -> pakai proxy konservatif dari identitas/cov
    pseudo_ev = 1e-30 if (ident >= ID_MED and cov >= 0.8) else (1e-3 if ident >= ID_LOW else 1.0)
    conf, label, flag = _confidence_from_hit(ident, pseudo_ev, cov)
    flag = (flag + " | FALLBACK pure-python (DIAMOND tak tersedia); e-value adalah proxy").strip(" |")
    return Assignment(
        gene=bnum, gene_name=info.get("name"), ec=info.get("ec", []),
        identity=round(ident, 2), evalue=pseudo_ev, bitscore=0.0, aln_cov=round(cov, 3),
        method="fallback-python-pairwise", confidence=conf,
        confidence_label=label, flag=flag,
    )


# 3. map_to_reactions - GPR iML1515 (b-number / EC -> reaction ids)
def map_to_reactions(gene_or_ec: str, model: cobra.Model,
                     genemap: dict | None = None,
                     base_confidence: float = 1.0) -> list[ReactionHit]:
    """Map gene (b-number/name) or EC to iML1515 reactions via GPR: b-number (e.g. b4024) -> model.genes; gene name (e.g. lysC) -> .name match; EC (e.g. 2.7.2.4) -> ec-code annotation match. Reaction confidence = base_confidence (from assign) times GPR-complexity penalty."""
    if genemap is None:
        genemap = load_genemap()
    hits: list[ReactionHit] = []
    tok = gene_or_ec.strip()

    # EC code?
    if tok.count(".") == 3 and all(p.isdigit() or p == "-" for p in tok.split(".")):
        for r in model.reactions:
            ec = r.annotation.get("ec-code")
            ec_list = ec if isinstance(ec, list) else ([ec] if isinstance(ec, str) else [])
            if tok in ec_list:
                # EC bisa banyak->banyak: penalti
                hits.append(ReactionHit(
                    reaction_id=r.id, via=f"ec:{tok}",
                    confidence=round(base_confidence * 0.7, 4),
                    note="dipetakan via EC (bisa many-to-many; confidence dipotong)",
                ))
        return hits

    # b-number langsung
    gene = None
    if tok in model.genes:
        gene = model.genes.get_by_id(tok)
    else:
        # nama gen
        for g in model.genes:
            if g.name == tok:
                gene = g
                break
    if gene is None:
        return hits
    rxns = list(gene.reactions)
    for r in rxns:
        # gen yang mengkatalisis >1 reaksi: ambiguitas -> sedikit potong
        penalty = 1.0 if len(rxns) == 1 else 0.85
        hits.append(ReactionHit(
            reaction_id=r.id, via=f"gene:{gene.id}",
            confidence=round(base_confidence * penalty, 4),
            note=("GPR langsung (gen 1 reaksi)" if len(rxns) == 1
                  else f"gen mengkatalisis {len(rxns)} reaksi (ambigu)"),
        ))
    return hits


# 4. annotate - full cascade
def annotate(dna: str, model: cobra.Model | None = None,
             genemap: dict | None = None, prefer_ec: bool = False) -> AnnotationResult:
    """Full cascade: DNA -> ORF -> gene/EC (homology) -> iML1515 reaction; cascade confidence = max over ORF of (orf_conf * assign_conf * rxn_conf); low_confidence set if < ID_LOW-equivalent or no reaction mapped."""
    if model is None:
        model = cobra.io.load_json_model(MODEL_CACHE)
    if genemap is None:
        genemap = load_genemap()

    res = AnnotationResult()
    orfs = orf_call(dna)
    res.orfs = orfs
    if not orfs:
        res.notes.append("Tak ada ORF terpanggil (DNA terlalu pendek / non-coding).")
        res.low_confidence = True
        return res

    method_seen = set()
    best_chain = 0.0
    for orf in orfs:
        asn = assign_function(orf.protein, genemap=genemap)
        res.assignments.append(asn)
        method_seen.add(asn.method)
        # pilih jalur pemetaan: gen (default) atau EC
        rxn_hits: list[ReactionHit] = []
        if asn.gene and not prefer_ec:
            rxn_hits = map_to_reactions(asn.gene, model, genemap=genemap,
                                        base_confidence=asn.confidence)
        if not rxn_hits and asn.ec:
            for ec in asn.ec:
                rxn_hits += map_to_reactions(ec, model, genemap=genemap,
                                             base_confidence=asn.confidence)
        # chain confidence
        for rh in rxn_hits:
            chain = round(orf.confidence * rh.confidence, 4)
            rh.confidence = chain
            best_chain = max(best_chain, chain)
        res.reactions.extend(rxn_hits)

    # dedup reaksi (ambil confidence tertinggi per reaction_id)
    dedup: dict[str, ReactionHit] = {}
    for rh in res.reactions:
        if rh.reaction_id not in dedup or rh.confidence > dedup[rh.reaction_id].confidence:
            dedup[rh.reaction_id] = rh
    res.reactions = sorted(dedup.values(), key=lambda x: x.confidence, reverse=True)

    res.confidence = round(best_chain, 4)
    res.method = "+".join(sorted(method_seen)) if method_seen else "none"
    res.low_confidence = (
        best_chain < (ID_LOW / 100.0) or not res.reactions
        or "fallback-python-pairwise" in method_seen
    )
    if "fallback-python-pairwise" in method_seen:
        res.notes.append("MODE FALLBACK aktif (DIAMOND tak tersedia): "
                         "homologi via pairwise pure-python; e-value adalah proxy.")
    if not res.reactions:
        res.notes.append("Tak ada reaksi iML1515 terpeta dari ORF manapun.")
    return res


# CLI ringan utk inspeksi manual
if __name__ == "__main__":
    import sys
    dna = sys.argv[1] if len(sys.argv) > 1 else ""
    if not dna:
        print("usage: python3 annotate.py <DNA>  (atau import sbg modul)")
        print("DIAMOND:", _diamond_bin() or "TIDAK ADA -> mode fallback")
        sys.exit(0)
    m = cobra.io.load_json_model(MODEL_CACHE)
    r = annotate(dna, m)
    print(json.dumps(r.as_dict(), indent=2, ensure_ascii=False))
