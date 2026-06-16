ModalFBA-GenLoop is an in-silico pipeline that turns raw generative-model DNA
into ranked, lab-ready metabolic-engineering designs.

## Installation

Python 3.10+ recommended.

```bash
pip install -r requirements.txt
```

Core dependencies (pinned where it matters):

- `cobra==0.31.1`  **(PINNED - required; later cobra changes GECKO handling)**
- `numpy`, `scipy`, `matplotlib`
- `pyrodigal`  (gene calling for the annotation cascade)
- `requests`  (NVIDIA NIM / data fetch)

External / optional tools (not installed by pip):

- **DIAMOND** - homology step of the annotation cascade. Install separately
  (`conda install -c bioconda diamond` or download a release binary) and make
  sure `diamond` is on `$PATH`.
- **Evo2 via NVIDIA NIM** - optional generative front-end. Requires a free
  NVIDIA API key from <https://build.nvidia.com/arc/evo2-40b> ("Get API Key").
  **The key must be user-supplied and is never committed.** Provide it via
  either:
  ```bash
  export NVIDIA_API_KEY=nvapi-xxxx          # in your shell / ~/.bashrc
  # or write it to a file readable by all shells:
  echo 'nvapi-xxxx' > ~/.nvidia_evo2_key
  ```
  `code/evo2_nim.py` reads `$NVIDIA_API_KEY` first, then falls back to
  `~/.nvidia_evo2_key`. Neither location is inside the repo.

The FBA / ec-FBA stack is developed and run under **WSL** (Linux). On Windows,
set `PYTHONUTF8=1` before running any script (the figure/design-card writers
emit UTF-8).

## Data

Large reference data are **not committed**. Download them into `code/data/`
(and `code/` for the iML1515 cache) before running:

| File | Source |
|------|--------|
| `iML1515.json` / `iML1515.xml` | BiGG: <http://bigg.ucsd.edu/static/models/iML1515.xml> |
| `eciML1515_batch.*` (GECKO ec-model) | SysBioChalmers ecModels: <https://github.com/SysBioChalmers/ecModels> (`eciML1515/model/`) |
| `bigg_universal_model.json` | BiGG universal model: <http://bigg.ucsd.edu/static/namespace/universal_model.json> |
| `ecoli_k12_genome.fasta` / `.gb` | NCBI / Ensembl Bacteria (E. coli K-12 MG1655) |
| `ecoli_k12_proteome.fasta` | UniProt proteome UP000000625 (E. coli K-12) |
| `iml1515_gene_ref.fasta` / `.dmnd` | Regenerate from iML1515 gene sequences; build the DIAMOND DB with `diamond makedb --in iml1515_gene_ref.fasta -d iml1515_gene_ref` |
| `ecmdb.json.zip` | ECMDB: <https://ecmdb.ca> |
| `ymdb.json.zip` | YMDB: <https://ymdb.ca> |

Small result artifacts (JSON reports, design FASTAs, figures, the M9 medium
definition, etc.) are committed so reviewers can inspect results without any
downloads. See `REPO-STRUCTURE.md` for the full include/exclude list.

## Usage

Run from `code/` (prefix with `PYTHONUTF8=1` on Windows):

```bash
# Resource-aware benchmark (M3 growth-coupling)
python benchmark_m3.py

# AUROC with DeLong confidence intervals
python auroc_ci_delong.py

# Genome-scale screens
python screen_genomewide.py coupling
python screen_genomewide.py fseof
python screen_genomewide.py heterologous

# ec-FBA overexpression demo (eciML1515) - shows plain-FBA Δ=0 vs modal Δ>0
python ecfba_overexpression_demo.py

# Blind held-out benchmark
python benchmark_blind_heldout.py

# Lab-ready output (top-10 design card FASTA + scorecard)
python make_design_card.py
```

Note: this repository ships the minimal essential set (core pipeline + scripts that
reproduce the reported results and the lab-ready output). Unit tests, figure-rendering
scripts, and intermediate milestone demos are omitted for a clean package; the rendered
figures and all result artifacts are included under `code/data/`.

## Results

- **Annotation cascade:** 4/4 round-trip recovery on the test set.
- **Feedback-resistant lysine (modal vs plain FBA):** plain stoichiometric FBA
  Δ = 0; resource-aware modal oracle Δ ≈ **+0.531**.
- **ec-FBA overexpression (eciML1515):** overexpression designs become
  rankable (e.g. dapA +0.24, lysC +0.86) where plain FBA cannot distinguish them.
- **AUROC (reconciled):** **0.90–0.99** across configurations, with DeLong
  significance.
- **Genome-scale screens:** coupling, FSEOF, and heterologous-audit screens run
  across the whole secretable exchange set.
- **Blind held-out benchmark:** AUROC **0.96**.
- **Design card:** ranked **top-10** lab-ready FASTA with honest in-silico
  limits stated.

## Citation

> *Placeholder - update on acceptance.*
> Bangun & Kusumawaty. "A resource-aware flux balance for ranking
> generative-model DNA designs in *Escherichia coli* metabolic engineering."
> *ACS Synthetic Biology* (in preparation), 2026.

## License

MIT (suggested). See `LICENSE`.
