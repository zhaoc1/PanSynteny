#!/bin/bash
# ------------------------------------------------------------------------------
# get_neighbor.sh  <gene_id> <genes_file> [n_genes]
#
# For one gene member in one genome: emit that gene plus up to n_genes flanking
# genes on either side along the same contig.
#   - double-passes the .genes file: pass 1 records the contig(s) the gene is
#     on, pass 2 emits every gene on that contig
#   - `grep -C n_genes` keeps the member +/- n_genes neighbours
#   - prepends the gene member id as a new first column
#
# n_genes defaults to 20 if not supplied; generate_neighbor_list.sh passes the
# value resolved from neighbor.n_genes in the YAML config.
#
# Innermost script of the neighbor-extraction chain
#   run_species.sh -> generate_neighbor_list.sh -> get_neighbor.sh
# Copied verbatim from
#   mwas-neighbor-pangraph/pipeline_v1/step1_gene_neighbors/get_neighbor.sh
# (no hardcoded paths — nothing to integrate).
# ------------------------------------------------------------------------------

set -e

if [ $# -lt 2 ] || [ $# -gt 3 ]; then
    echo "Usage: $0 <gene_id> <genes_file> [n_genes]"
    exit 1
fi

gene_id="$1"
gene_file="$2"
n_genes="${3:-20}"  # default to 20 if not given

awk -v pat="$gene_id" 'FNR==NR { if ($1 ~ pat) a[$2]; next } ($2 in a)' $gene_file $gene_file | \
    grep -C $n_genes $gene_id | \
    awk -v OFS="\t" -v var="$gene_id" '{print var, $0}'
