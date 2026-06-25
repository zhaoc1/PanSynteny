#!/usr/bin/env python3
"""
Standalone Stage-2 converter: Prokka GFF3  ->  per-genome `.genes` TSV.

This reproduces MIDAS `parse_gff_to_tsv()` (midas/common/utilities.py) as a
self-contained script. The only external dependency is `gffutils`.

Output `.genes` columns (TAB-delimited, with header):
    gene_id  contig_id  start  end  strand  gene_type

Usage:
    ./gff_to_genes.py GENOME.gff [-o GENOME.genes] [--keep-db]
    ./gff_to_genes.py GENOME.gff           # writes GENOME.genes next to input
"""
import argparse
import gzip
import os
import sys

import gffutils


# Column order of the `.genes` file (genes_feature_schema in midas/params/schemas.py)
GENES_COLUMNS = ["gene_id", "contig_id", "start", "end", "strand", "gene_type"]


def _open_out(path):
    """Open output for writing, transparently gzip-ing if path ends in .gz."""
    if path.endswith(".gz"):
        return gzip.open(path, "wt")
    return open(path, "w")


def parse_gff_to_tsv(gff3_file, genes_file, keep_db=False):
    """Convert a Prokka GFF3 file into MIDAS `.genes` feature TSV.

    Mirrors midas/common/utilities.py:parse_gff_to_tsv
    """
    db_path = f"{gff3_file}.db"

    # Clean any stale outputs / sqlite db from a previous run.
    for stale in (genes_file, db_path):
        if os.path.exists(stale):
            os.remove(stale)

    # gffutils builds an on-disk sqlite db from the GFF3.
    db = gffutils.create_db(gff3_file, db_path)

    n_written = 0
    with _open_out(genes_file) as stream:
        stream.write("\t".join(GENES_COLUMNS) + "\n")
        for feature in db.all_features():
            # Prokka emits a `prokka`-sourced record per contig; skip those.
            if feature.source == "prokka":
                continue
            # Features without an ID (e.g. CRISPR repeats) are not genes.
            if "ID" not in feature.attributes:
                continue

            seqid = feature.seqid
            start = feature.start
            stop = feature.stop
            strand = feature.strand
            gene_id = feature.attributes["ID"][0]
            locus_tag = feature.attributes["locus_tag"][0]
            assert gene_id == locus_tag, \
                f"ID ({gene_id}) != locus_tag ({locus_tag}) for a feature in {gff3_file}"
            gene_type = feature.featuretype

            stream.write("\t".join([
                gene_id, seqid, str(start), str(stop), strand, gene_type,
            ]) + "\n")
            n_written += 1

    if not keep_db and os.path.exists(db_path):
        os.remove(db_path)

    return n_written


def main():
    parser = argparse.ArgumentParser(
        description="Convert a Prokka GFF3 file into a MIDAS per-genome .genes TSV.")
    parser.add_argument("gff3", help="Input Prokka GFF3 file ({genome_id}.gff)")
    parser.add_argument(
        "-o", "--output",
        help="Output .genes path (default: input with .gff replaced by .genes)")
    parser.add_argument(
        "--keep-db", action="store_true",
        help="Keep the intermediate gffutils sqlite db ({gff3}.db)")
    args = parser.parse_args()

    if not os.path.isfile(args.gff3):
        sys.exit(f"ERROR: input GFF3 not found: {args.gff3}")

    if args.output:
        genes_file = args.output
    elif args.gff3.endswith(".gff"):
        genes_file = args.gff3[:-len(".gff")] + ".genes"
    else:
        genes_file = args.gff3 + ".genes"

    n = parse_gff_to_tsv(args.gff3, genes_file, keep_db=args.keep_db)
    print(f"Wrote {n} gene features -> {genes_file}")


if __name__ == "__main__":
    main()
