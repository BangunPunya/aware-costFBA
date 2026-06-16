# data/ - download instructions for large external files

Small result artifacts (benchmark/screen JSON, design FASTA, figures, CSV caches)
are committed here. The large reference models/genomes below are **not committed**
(see `.gitignore`); download them into this folder before running the pipeline.

| File | Source |
|------|--------|
| `iML1515.json` (and `code/iML1515.json`) | BiGG: http://bigg.ucsd.edu/static/models/iML1515.xml (convert to JSON, or load XML) |
| `eciML1515_batch.xml`, `eciML1515_batch.mat` | GECKO ecModels: https://github.com/SysBioChalmers/ecModels (path `eciML1515/model/`) |
| `eciML1515_batch_fixed.mat` | Regenerate from `eciML1515_batch.mat` by sanitising the one whitespace reaction id (`protein pseudoreaction` -> `protein_pseudoreaction`); see `ecfba_overexpression_demo.py` notes |
| `bigg_universal_model.json` | BiGG universal model: http://bigg.ucsd.edu/static/namespace/universal_model.json |
| `ecoli_k12_genome.fasta`, `ecoli_k12_genome.gb` | NCBI/Ensembl E. coli K-12 MG1655 (NC_000913.3) |
| `ecoli_k12_proteome.fasta` | UniProt proteome UP000000625 |
| `iml1515_gene_ref.fasta` | Regenerate from iML1515 gene products (DIAMOND reference) |
| `iml1515_gene_ref.dmnd` | Rebuild: `diamond makedb --in iml1515_gene_ref.fasta -d iml1515_gene_ref` |
| `ecmdb.json.zip`, `ymdb.json.zip` | ECMDB (https://ecmdb.ca) / YMDB (https://ymdb.ca) metabolite dumps |

After downloading, run scripts with `PYTHONUTF8=1` on Windows. See top-level `README.md` for usage.
