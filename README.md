# PanSynteny: de novo discovery of recurrent syntenic gene neighborhoods across bacterial pangenomes

> PanSynteny is a pangenome-based framework for de novo discovery of recurrent syntenic gene neighborhoods. Starting from user-defined focal genes, PanSynteny extracts local genomic neighborhoods, canonicalizes orientation, collapses recurring gene-cluster paths across genomes, and reports conserved co-localized neighborhood blocks with genome-level traceability.

**Status:** Testing release v0.2.0. Prototype — the YAML schema and output column set may change between 0.x.y versions. Not yet for production.

Given (1) a pan-genome (centroid_80 clusters from MIDAS) and (2) a user-curated focal-gene table (which centroid_80 clusters to investigate), this pipeline reconstructs **de novo recurring operons** the focal genes live in and harmonizes them across the strains in the species-level pangenome. Outputs are emitted at three granularity levels: coarse canonical operon, length-variant isoform, and per-genome instance - so you can ask either "what operon is this" or "which strains carry which variant." Optionally, it also extracts the contiguous focal-gene blocks (grouped by `focal_label` direction) within each operon. 

---

## Installation

The full install (R + Python + conda env) is documented in **[SETUP.md](docs/SETUP.md)**. TL;DR:

```bash
# One-shot conda env (pulls R + pyyaml + gffutils)
conda env create -f environment.yml
conda activate strain-aware-operon

# One R package not in conda
R -e "install.packages('randomcoloR', repos='https://cloud.r-project.org/')"

# Sanity check
python -c "import yaml, gffutils; print(yaml.__version__, gffutils.__version__)"
```

---

## Quickstart

1. **Edit `example.yaml` for your data:**
   ```bash
   # in example.yaml, set: job.species_id, job.proj_dir, job.input_dir, data.midasdb_dir,
   # data.data_dir, data.focal_meta, data.clusters_80_updated, and the sources: list
   ```
   See [USER_GUIDE.md §Configuration](docs/USER_GUIDE.md#configuration-yaml) for every key. (Tip: `cp example.yaml my_config.yaml` first if you want to keep the template clean.)

2. **Run the four ordered commands** - all read the same YAML:

   ```bash
   # Step 0a — build the unified genome catalog from sources: in the YAML
   python build_genome_catalog.py example.yaml

   # Step 0 — snapshot config, cache focal_meta, enumerate missing neighbor TSVs
   Rscript prepare.R example.yaml

   # Step 0 — materialise any missing per-focal neighbor TSVs
   bash run_species.sh example.yaml

   # Steps 1–6 — the analytical pipeline
   Rscript pipeline.R example.yaml
   ```

   **Tip — if `conda activate` doesn't take** (IDE terminals, scripts that skip shell init), invoke the env's binaries directly:
   ```bash
   ENV=$(conda info --base)/envs/strain-aware-operon
   PY=$ENV/bin/python
   RSC=$ENV/bin/Rscript
   export LD_LIBRARY_PATH=$ENV/lib:$LD_LIBRARY_PATH

   $PY build_genome_catalog.py example.yaml
   $RSC prepare.R example.yaml
   bash run_species.sh example.yaml
   $RSC pipeline.R example.yaml
   ```

3. **Read your outputs** under `{proj_dir}/`:

   | Folder / file | What's in it |
   | --- | --- |
   | `step1_setup/` | Run config snapshot, focal_meta cache, genome catalog |
   | `step2_neighbors/` | Per-genome operon graphs (cache for Step 1) |
   | `step3_path/canonical_paths*.tsv` | Operons at three granularity levels (coarse → fine → per-genome) |
   | `step4_parse/` | Operon summaries, fine-isoform selection, BLAST gene-id lists |
   | `step5_figures/*.pdf` | gggenes operon visualizations |
   | `step6_blocks/rep.tsv` | Trait-associated blocks + per-strain attribution (gated by `blocks.skip_block`) |

   Column-level schemas live in [SCHEMA.md](docs/SCHEMA.md).

---

## Documentation

| Doc | Read when... |
| --- | --- |
| **[USER_GUIDE.md](docs/USER_GUIDE.md)** | You're configuring a run — YAML schema, tunables, output reference. |
| **[SETUP.md](docs/SETUP.md)** | You're installing dependencies. |
| **[STEPS.md](docs/STEPS.md)** | You need to know exactly what one step does (inputs, outputs, logic, caveats). |
| **[SCHEMA.md](docs/SCHEMA.md)** | You need column-level schemas for any file the pipeline reads or writes. |
| **[PIPELINE.md](docs/PIPELINE.md)** | You need the c80-column glossary or the truncation/fragmentation flag semantics. |
| **[diagram.md](docs/diagram.md)** | You want the YAML-key → consumer-step data flow. |
| **[CLAUDE.md](CLAUDE.md)** | You're a Claude / AI agent working on this repo (keystones, conventions, gotchas). |

---

## Citation

If you use this software, please cite:

> [Manuscript in preparation. Contact the author for citation guidance.]

## License

MIT — see [LICENSE](LICENSE).

## Contact

Chunyu Zhao — <chunyu.zhao@gladstone.ucsf.edu>
