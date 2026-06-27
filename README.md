# PanSynteny: de novo discovery of recurrent syntenic gene neighborhoods across bacterial pangenomes

> Find the recurrent, co-localized syntenic gene neighborhoods your genes of interest live in, and how they vary across the strains of a bacterial species.

**Status:** Testing release v0.5.0. Prototype - the YAML schema and output column set may change between 0.x.y versions. Not yet for production.

**What it does.** Point PanSynteny at a pangenome (centroid_80 clusters from MIDAS) and a short list of focal genes you care about, and it finds the **recurrent, co-localized syntenic gene neighborhoods** those genes sit in, and shows how each one varies across the strains of the species. You get plain TSV tables and gene-arrow (gggenes) figures at three levels of detail: the consensus neighborhood, its length variants, and the exact genes in each genome. Optionally, it also pulls out the trait-associated gene blocks within each neighborhood.

> **A note on terminology.** Elsewhere in these docs (and in the code), a recurrent, co-localized syntenic gene neighborhood is called an **operon-like structure** or **operon** as shorthand. PanSynteny makes **no claim** of shared transcription, strand co-orientation, or operonic regulation - the structures are defined purely by **recurrent physical co-localization** across genomes.

---

## Installation

The full install (R + Python + conda env) is documented in **[SETUP.md](docs/SETUP.md)**. TL;DR:

```bash
# One-shot conda env (pulls R + pyyaml + gffutils)
conda env create -f environment.yml
conda activate pansynteny

# Sanity check
python -c "import yaml, gffutils; print(yaml.__version__, gffutils.__version__)"
```

---

## Quickstart

> **Prerequisites:** the `pansynteny` conda env (see [SETUP.md](docs/SETUP.md)) and a focal-gene table (TSV).

```bash
conda activate pansynteny
```

**1. Make your run config** - copy the template and fill in your paths:

```bash
cp example.yaml my_run.yaml
# edit my_run.yaml: job.{species_id, proj_dir, input_dir},
# data.{midasdb_dir, data_dir, focal_meta, clusters_80_updated}, and the sources: list
```

Every key is documented in [USER_GUIDE.md section Configuration](docs/USER_GUIDE.md#configuration-yaml).

**2. Run the four ordered commands** - all read the same config:

```bash
python build_genome_catalog.py my_run.yaml   # Step 0a - build the genome catalog from sources:
Rscript prepare.R my_run.yaml                # Step 0b - cache focal_meta, list missing neighbor TSVs
bash build_neighbor_lists.sh my_run.yaml     # Step 0c - materialise the neighbor TSVs
Rscript pipeline.R my_run.yaml               # Steps 1-6 - neighborhoods -> figures -> trait blocks
```

**3. Read your outputs** under `{proj_dir}/`:

| Folder / file | What's in it |
| --- | --- |
| `step1_setup/` | Run config snapshot, focal_meta cache, genome catalog |
| `step2_neighbors/` | Per-genome neighborhood graphs (cache for Step 1) |
| `step3_path/canonical_paths*.tsv` | Gene neighborhoods at three granularity levels (coarse -> fine -> per-genome) |
| `step4_parse/` | Neighborhood summaries, fine-isoform selection, BLAST gene-id lists |
| `step5_figures/*.pdf` | gggenes neighborhood visualizations |
| `step6_blocks/rep.tsv` | Trait-associated blocks + per-strain attribution (gated by `blocks.skip_block`) |

Column-level schemas live in [SCHEMA.md](docs/SCHEMA.md).

> **If `conda activate` doesn't stick** (IDE terminals, non-login shells), call the env's interpreters directly - set the lib path and prepend `$ENV/bin/`:
> ```bash
> ENV=$(conda info --base)/envs/pansynteny
> export LD_LIBRARY_PATH=$ENV/lib:$LD_LIBRARY_PATH
> # then: $ENV/bin/python build_genome_catalog.py ..., $ENV/bin/Rscript prepare.R ...  (bash scripts run as-is)
> ```

---

## Documentation

| Doc | Read when... |
| --- | --- |
| **[USER_GUIDE.md](docs/USER_GUIDE.md)** | You're configuring a run - YAML schema, tunables, output reference. |
| **[SETUP.md](docs/SETUP.md)** | You're installing dependencies. |
| **[STEPS.md](docs/STEPS.md)** | You need to know exactly what one step does (inputs, outputs, logic, caveats). |
| **[SCHEMA.md](docs/SCHEMA.md)** | You need column-level schemas for any file the pipeline reads or writes. |
| **[PIPELINE.md](docs/PIPELINE.md)** | You need the c80-column glossary or the truncation/fragmentation flag semantics. |
| **[diagram.md](docs/diagram.md)** | You want the YAML-key -> consumer-step data flow. |
| **[CLAUDE.md](CLAUDE.md)** | You're a Claude / AI agent working on this repo (keystones, conventions, gotchas). |

---

## Citation

If you use this software, please cite:

> [Manuscript in preparation. Contact the author for citation guidance.]

## License

MIT - see [LICENSE](LICENSE).

## Contact

Chunyu Zhao - <chunyu.zhao@gladstone.ucsf.edu>
