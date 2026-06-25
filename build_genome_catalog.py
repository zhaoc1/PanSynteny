#!/usr/bin/env python3
"""
build_genome_catalog.py <config.yaml>

Step 0a - build the unified genome catalog the neighbor-extraction chain and
the R pipeline both consume. Reads the `sources:` list from the YAML and
normalises every source into two source-agnostic artefacts (mirrors model.R
get_target):

  {proj_dir}/step1_setup/catalog_genes_info.tsv
      header:  gene_id <TAB> centroid_80 <TAB> gene_length
      body:    union of (gene, c80, length) rows across every source
  {proj_dir}/step1_setup/catalog_genome_toc.tsv
      header:  genome_id <TAB> genes_file_path
      body:    union of (genome_id, path) rows across every source

Catalog lives under proj_dir (not data_dir) so each run carries its own copy -
different `sources:` lists across proj_dirs no longer clobber a shared catalog.
The expensive per-genome `.genes` files (prokka source conversions) stay under
`<genomes_dir>` and remain shared across runs via the `-s` idempotency guard.

The 3-column genes_info schema is exactly what `load_c80_tables` (midas.R)
projects out of the legacy midasdb `genes_annotated.tsv`, so this file is a
drop-in for pipeline.R via `get_target("catalog_genes_info")`.

Each `sources:` entry:
  name         free label (for logging)
  type         midas | prokka
  genes_info   membership file: gene_id in col 1, centroid_80 in col `c80_col`,
               gene_length in col `length_col`
  c80_col      1-based column of centroid_80 in `genes_info`
  length_col   1-based column of gene_length in `genes_info`. Defaults to 8
               for midas (UHGG schema). **Required for prokka** - the prokka
               membership schema is user-controlled, so it must be declared.
  genomes_dir  directory of <genome_id>/<genome_id>.genes  (uniform contract)

Both source types read gene_length the same way: straight from the column.
For prokka, `<genome_id>/<genome_id>.gff` is still converted IN PLACE to
`<genome_id>.genes` via gff_to_genes.parse_gff_to_tsv when missing (idempotent
on `-s`) - that file is needed by downstream get_neighbor.sh, but the catalog
build no longer scans it for gene length.

Global placeholders expanded in every source string field: `{species_id}`,
`{midasdb_dir}`, `{proj_dir}`, `{input_dir}` (mirrors config.R's
`load_job_config`). Each
source entry may also declare its own *local* placeholders by adding extra
string/int fields (any key other than the reserved `genes_info` / `genomes_dir`
/ `c80_col` / `length_col` / `name` / `type` becomes a `{key}` placeholder
usable inside `genes_info` and `genomes_dir` for that same source). Example:

  - name:        ecor
    type:        prokka
    ecor_dir:    "/data/.../ecor72"
    genes_info:  "{ecor_dir}/ecor_gene_centroid80.tsv"
    genomes_dir: "{ecor_dir}/ecor_prokka"

A genome_id must come from exactly one source: the genome->path map has to be
unambiguous, so a duplicate genome_id across sources is reported and the build
STOPS (the partially-written genes_info.tsv is removed; genome_toc.tsv is only
written after the dup check passes, so it's never left in a bad state).

`sources:` is consumed only here - config.R ignores it (it is a list, not a
scalar section), so the R pipeline is unaffected.

Run order: build_genome_catalog.py -> prepare.R -> build_neighbor_lists.sh -> pipeline.R

Requires the pansynteny conda env (pyyaml + gffutils).
"""

import argparse
import os
import sys
import time
from collections import defaultdict
from pathlib import Path


# --- Hard-fail early with a clear message if the runtime is wrong -------------
try:
    import yaml
except ImportError:
    sys.exit("ERROR: pyyaml not importable - activate the pansynteny "
             "conda env (or set its python on PATH).")

# Direct import beats subprocess: one python process, no env-detection dance.
# gff_to_genes.py lives alongside the other helpers under scripts/.
sys.path.insert(0, str(Path(__file__).resolve().parent / "scripts"))
try:
    from gff_to_genes import parse_gff_to_tsv
except ImportError as exc:
    sys.exit(f"ERROR: cannot import gff_to_genes ({exc}) - activate the "
             "pansynteny conda env (it needs gffutils).")


GENES_INFO_NAME    = "catalog_genes_info.tsv"
GENOME_TOC_NAME    = "catalog_genome_toc.tsv"
GENES_INFO_HEADER  = "gene_id\tcentroid_80\tgene_length\n"
GENOME_TOC_HEADER  = "genome_id\tgenes_file_path\n"
PROGRESS_EVERY     = 5_000_000   # log a heartbeat every N membership rows scanned
DEFAULT_LENGTH_COL = 8           # UHGG `genes_info.tsv` puts gene_length in col 8


def genome_id_from_gene_id(gene_id: str) -> str:
    """Strip the trailing _NNNNN field from a gene_id.

    Shared contract with focal_neighbor_list.sh's awk join. Works for both
    `GUT_GENOME000040_00388` -> `GUT_GENOME000040` and
    `GCF_900448275.1_00001` -> `GCF_900448275.1`.
    """
    return gene_id.rsplit("_", 1)[0]


def load_config(config_path: Path):
    """Parse the YAML, return (species_id, proj_dir, [sources]).

    `proj_dir` is used as-is (no implicit species_id suffix). To isolate
    multi-species runs, set `proj_dir: "/path/.../<species_id>"` in the YAML.
    """
    with open(config_path) as fh:
        cfg = yaml.safe_load(fh)

    species_id  = str(cfg["job"]["species_id"])
    proj_dir    = os.path.expanduser(str(cfg["job"]["proj_dir"]))
    # input_dir is required (no backward-compat fallback)
    if not cfg["job"].get("input_dir"):
        sys.exit("ERROR: job.input_dir is required in <config.yaml>. Add the absolute "
                 "path to your user-provided inputs under job:.")
    input_dir   = os.path.expanduser(str(cfg["job"]["input_dir"]))
    midasdb_dir = os.path.expanduser(str(cfg["data"]["midasdb_dir"]))

    raw_sources = cfg.get("sources") or []
    if not raw_sources:
        sys.exit(f"ERROR: no `sources:` list in {config_path}")

    global_placeholders = {
        "{species_id}":  species_id,
        "{midasdb_dir}": midasdb_dir,
        "{proj_dir}":    proj_dir,
        "{input_dir}":   input_dir,
    }

    def expand_globals(value: str) -> str:
        out = str(value)
        for ph, repl in global_placeholders.items():
            out = out.replace(ph, repl)
        return out

    # Keys whose values describe paths/columns the catalog reads directly.
    # Anything else (e.g. `ecor_dir: "..."`) is treated as a per-source local
    # placeholder usable inside `genes_info` / `genomes_dir` of the same entry.
    RESERVED = {"name", "type", "genes_info", "c80_col", "length_col", "genomes_dir"}

    sources = []
    required = ("name", "type", "genes_info", "c80_col", "genomes_dir")
    for s in raw_sources:
        missing = [k for k in required if k not in s]
        if missing:
            sys.exit(f"ERROR: source '{s.get('name', '?')}' missing keys: {missing}")
        if s["type"] not in ("midas", "prokka"):
            sys.exit(f"ERROR: source '{s['name']}' has unknown type {s['type']!r} "
                     "(expected 'midas' or 'prokka')")
        # length_col: midas defaults to UHGG's col 8; prokka must declare it,
        # since its membership schema is user-controlled.
        if s["type"] == "prokka" and "length_col" not in s:
            sys.exit(f"ERROR: source '{s['name']}' (prokka) requires `length_col` - "
                     "prokka membership must include gene_length as a named column.")
        length_col = int(s.get("length_col", DEFAULT_LENGTH_COL))

        # Local placeholders: any non-reserved scalar field. Expand globals
        # inside the local values too so e.g. `ecor_dir: "{proj_dir}/ecor"` works.
        local_placeholders = {
            "{" + k + "}": expand_globals(v)
            for k, v in s.items()
            if k not in RESERVED and isinstance(v, (str, int, float))
        }

        def expand_all(value: str) -> str:
            out = expand_globals(value)
            for ph, repl in local_placeholders.items():
                out = out.replace(ph, str(repl))
            return out

        sources.append({
            "name":        str(s["name"]),
            "type":        str(s["type"]),
            "genes_info":  Path(os.path.expanduser(expand_all(s["genes_info"]))),
            "c80_col":     int(s["c80_col"]),
            "length_col":  length_col,
            "genomes_dir": Path(os.path.expanduser(expand_all(s["genomes_dir"]))),
        })
    return species_id, Path(proj_dir), sources


def _scan_membership(src: dict, genes_info_out, must_have_idx: int) -> set:
    """Single pass over a source's membership file: write rows + collect genomes.

    Each kept row writes `gene_id <TAB> centroid_80 <TAB> gene_length` to
    `genes_info_out` (length from column `src["length_col"]`). Returns the set
    of unique genome_ids contributed by this source.
    """
    t0 = time.monotonic()
    c80_idx    = src["c80_col"]    - 1
    length_idx = src["length_col"] - 1
    src_genomes = set()
    n_rows, n_skipped = 0, 0
    with open(src["genes_info"]) as fh:
        next(fh, None)  # header
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) <= must_have_idx:
                n_skipped += 1
                continue
            gene_id, c80, gl = cols[0], cols[c80_idx], cols[length_idx]
            if not (gene_id and c80 and gl):
                n_skipped += 1
                continue
            genes_info_out.write(f"{gene_id}\t{c80}\t{gl}\n")
            src_genomes.add(genome_id_from_gene_id(gene_id))
            n_rows += 1
            if n_rows % PROGRESS_EVERY == 0:
                print(f"    ... {n_rows:>12,} rows scanned, "
                      f"{len(src_genomes):>6,} genomes  "
                      f"({time.monotonic() - t0:5.1f}s)", flush=True)
    print(f"    {n_rows:,} membership rows, {len(src_genomes):,} unique genomes "
          f"({time.monotonic() - t0:.1f}s)" +
          (f"  [{n_skipped:,} malformed rows skipped]" if n_skipped else ""),
          flush=True)
    return src_genomes


def process_midas_source(src: dict, genes_info_out, toc_entries: list):
    """Midas source - single pass: gene_length comes straight from a column
    in the membership file (UHGG's `genes_info.tsv` puts it at col 8). The
    .genes files are part of the midasdb and assumed present (not validated).
    """
    must_have = max(src["c80_col"], src["length_col"]) - 1
    src_genomes = _scan_membership(src, genes_info_out, must_have)
    # Emit TOC rows. Sorted for stable, deterministic output.
    for g in sorted(src_genomes):
        genes_fp = src["genomes_dir"] / g / f"{g}.genes"
        toc_entries.append((g, str(genes_fp), src["name"]))


def process_prokka_source(src: dict, genes_info_out, toc_entries: list):
    """Prokka source - gene_length must already be in the membership file
    (column `length_col`), like midas. The user's upstream BLAST + annotation
    pipeline is responsible for populating it (typically `end - start + 1`
    from the prokka GFF).

    Two stages (parallel to midas):
      1. Single pass over membership: write rows + collect genome set.
      2. For each genome: ensure <g>/<g>.genes exists (convert .gff -> .genes
         via gff_to_genes if missing) so downstream get_neighbor.sh can read
         it. Emit one TOC row per genome.
    """
    must_have = max(src["c80_col"], src["length_col"]) - 1
    src_genomes = _scan_membership(src, genes_info_out, must_have)
    n_genomes = len(src_genomes)

    n_converted = 0
    for g in sorted(src_genomes):
        genes_fp = src["genomes_dir"] / g / f"{g}.genes"
        if not (genes_fp.exists() and genes_fp.stat().st_size > 0):
            gff_fp = src["genomes_dir"] / g / f"{g}.gff"
            if not (gff_fp.exists() and gff_fp.stat().st_size > 0):
                sys.exit(f"ERROR: no .gff for {g}: {gff_fp}")
            parse_gff_to_tsv(str(gff_fp), str(genes_fp), keep_db=False)
            n_converted += 1
            if n_converted % 10 == 0:
                print(f"    ... derived {n_converted}/{n_genomes} .genes from .gff",
                      flush=True)
        toc_entries.append((g, str(genes_fp), src["name"]))
    print(f"    {n_converted} genes derived from gff "
          f"({n_genomes - n_converted} already present)", flush=True)


def process_source(src: dict, genes_info_out, toc_entries: list):
    """Dispatch to the per-type handler."""
    print(f">>> source '{src['name']}' (type={src['type']})", flush=True)
    if not src["genes_info"].is_file():
        sys.exit(f"ERROR: genes_info not found: {src['genes_info']}")
    if not src["genomes_dir"].is_dir():
        sys.exit(f"ERROR: genomes_dir not found: {src['genomes_dir']}")
    if src["type"] == "midas":
        process_midas_source(src, genes_info_out, toc_entries)
    else:  # prokka (validated in load_config)
        process_prokka_source(src, genes_info_out, toc_entries)


def check_duplicates(toc_entries):
    """Return {genome_id: [source_names]} for any genome_id that came from
       more than one source. Empty dict means clean."""
    by_genome = defaultdict(list)
    for genome_id, _path, src_name in toc_entries:
        by_genome[genome_id].append(src_name)
    return {g: srcs for g, srcs in by_genome.items() if len(srcs) > 1}


def main():
    ap = argparse.ArgumentParser(description="Build the unified genome catalog.")
    ap.add_argument("config", type=Path, help="Path to <config.yaml>.")
    args = ap.parse_args()

    if not args.config.is_file():
        sys.exit(f"ERROR: config not found: {args.config}")

    species_id, proj_dir, sources = load_config(args.config)
    # proj_dir is used as-is (no implicit species_id append).
    catalog_dir    = proj_dir / "step1_setup"
    catalog_dir.mkdir(parents=True, exist_ok=True)
    genes_info_out = catalog_dir / GENES_INFO_NAME
    genome_toc_out = catalog_dir / GENOME_TOC_NAME

    print(f">>> building genome catalog for species {species_id}")
    print(f"    catalog_dir = {catalog_dir}")

    # Stream genes_info to disk incrementally (it can be huge); hold the toc
    # in memory (small, ~thousands of rows) so the dup check can run before
    # any toc file is on disk.
    toc_entries: list[tuple[str, str, str]] = []
    with open(genes_info_out, "w") as gi_w:
        gi_w.write(GENES_INFO_HEADER)
        for src in sources:
            process_source(src, gi_w, toc_entries)

    # Duplicate-genome_id check: the genome->path map has to be unambiguous.
    dups = check_duplicates(toc_entries)
    if dups:
        print("ERROR: duplicate genome_id(s) across sources - catalog would be "
              "ambiguous:", file=sys.stderr)
        for genome_id, src_names in sorted(dups.items()):
            print(f"  {genome_id}: {', '.join(src_names)}", file=sys.stderr)
        genes_info_out.unlink(missing_ok=True)   # leave nothing half-written
        sys.exit(1)

    # Write the TOC now that we know it's clean.
    with open(genome_toc_out, "w") as toc_w:
        toc_w.write(GENOME_TOC_HEADER)
        for genome_id, path, _src_name in toc_entries:
            toc_w.write(f"{genome_id}\t{path}\n")

    # -1 for the header row we wrote.
    n_gi_rows = sum(1 for _ in open(genes_info_out)) - 1
    print(">>> catalog built:")
    print(f"    genes_info : {genes_info_out}  ({n_gi_rows:,} rows + header)")
    print(f"    genome_toc : {genome_toc_out}  ({len(toc_entries):,} genomes + header)")


if __name__ == "__main__":
    main()
