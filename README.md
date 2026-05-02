# pangenome-operons

> Strain-resolved operon catalogs from microbial pan-genome data.

**Status:** Testing release v0.1.0. Prototype — the YAML schema and output column set may change between 0.x.y versions. Not yet for production.

## What it does

Given (1) a pan-genome (centroid_80 clusters from MIDAS) and (2) per-gene trait statistics, this pipeline reconstructs the **operons** that trait-associated focal genes live in, harmonizes them across strains, and surfaces the contiguous trait-associated blocks within each operon. Outputs are emitted at three granularity levels — coarse canonical operon, length-variant isoform, and per-genome instance — so you can ask either "what operon is this" or "which strains carry which variant."

The distinguishing feature versus existing pan-genome tools (Roary, ppanggolin, PIRATE) is **operon-level rather than gene-level analysis**, with **per-strain attribution preserved** all the way through.

## Quickstart

```bash
# 1. Install dependencies (R packages + conda env). See README_SETUP.md.
Rscript install_packages.R

# 2. Configure: edit example.yaml to point at your data.
cp example.yaml my_run.yaml

# 3. Run.
Rscript prepare.R my_run.yaml      # Step 0: focal selection + missing-TSV check
Rscript pipeline.R my_run.yaml     # Steps 1-6: full analysis
```

Per-focal neighbor TSVs are produced by an external preprocessing step;
`prepare.R` reports any missing ones.

## What you get

Outputs land under `proj_dir` (with a per-species subdirectory; see [USER_GUIDE.md](USER_GUIDE.md) for the YAML schema):

| Artifact | Step | Purpose |
| --- | --- | --- |
| `step3_path/canonical_paths*.tsv` | 3 | Operons at three granularity levels |
| `step4_block/representative_path.tsv` | 4 | Trait-associated blocks, ranked per locus |
| `step4_block/rep.tsv` | 4 | Per-strain attribution for each block |
| `step5_parse/fine_*.tsv` | 5 | BLAST gene-id lists for downstream sequence work |
| `step6_figures/*.pdf` | 6 | gggenes operon visualizations |

## Documentation

- **[USER_GUIDE.md](USER_GUIDE.md)** — YAML schema, tunables, output reference, when to use which level.
- **[STEPS.md](STEPS.md)** — per-step input / output / logic / caveats. The deep reference.
- **[README_SETUP.md](README_SETUP.md)** — install, conda env, troubleshooting.
- **[diagram.md](diagram.md)** — YAML key → consumer-step data flow.

## Citation

If you use this software, please cite:

> [Manuscript in preparation. Contact the author for citation guidance.]

## License

MIT — see [LICENSE](LICENSE).

## Contact

Chunyu Zhao — <chunyu.zhao@gladstone.ucsf.edu>
