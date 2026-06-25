#!/bin/bash
# ------------------------------------------------------------------------------
# focal_neighbor_list.sh  <config.yaml> <query>
#
# For one focal centroid_80 (`query`): pull every gene member from the unified
# genome-catalog genes_info, resolve each member's per-genome .genes file via
# the genome_toc, and extract up to n_genes same-contig neighbours per member
# (via get_neighbor.sh). Writes one <out_dir>/<query>.tsv (7 cols, no header):
#   gene_member  neighbor_gene_id  contig_id  start  end  strand  gene_type
#
# Config-driven — reads the same YAML the R pipeline loads:
#   job.species_id  +  job.proj_dir   -> step1_setup/ (catalog lookup)
#   job.species_id  +  data.data_dir  -> list_of_neighbors/ (output)
#   data.n_genes                      -> flank size for get_neighbor.sh (def 20)
#
# Inputs (built by build_genome_catalog.py, mirror model.R get_target):
#   {proj_dir}/step1_setup/catalog_genes_info.tsv  gene_id <TAB> centroid_80 <TAB> gene_length
#   {proj_dir}/step1_setup/catalog_genome_toc.tsv  genome_id <TAB> .genes path
# Output:
#   {data_dir}/{species_id}/list_of_neighbors/<query>.tsv
#
# Source-agnostic: it never sees "UHGG vs ECOR" — build_genome_catalog.py has
# already normalised every source into the two catalog files above. The
# genome_id is derived from the gene_id (strip trailing _NNNNN); that derivation
# is the shared contract with build_genome_catalog.py's genome_id_from_gene_id().
#
# Copied from
#   mwas-neighbor-pangraph/pipeline_v1/step1_gene_neighbors/focal_neighbor_list.sh
# and adapted: hardcoded midasdb_dir + path-convention .genes lookup replaced by
# the config-driven genome catalog.
# ------------------------------------------------------------------------------

set -e

if [ $# -ne 2 ]; then
    echo "Usage: $0 <config.yaml> <query>"
    exit 1
fi

config="$1"
query="$2"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --- read config values from the YAML (same file the R pipeline loads) --------
yaml_get() {
    # yaml_get <config.yaml> <dotted.key> [default]
    python3 - "$@" <<'PYEOF'
import yaml, sys, os
cfg, key = sys.argv[1], sys.argv[2]
default = sys.argv[3] if len(sys.argv) > 3 else None
d = yaml.safe_load(open(cfg))
for part in key.split('.'):
    if isinstance(d, dict) and part in d:
        d = d[part]
    else:
        if default is not None:
            print(default); sys.exit(0)
        sys.exit("Missing key '%s' in %s" % (key, cfg))
print(os.path.expanduser(str(d)))
PYEOF
}

species_id=$(yaml_get "$config" job.species_id)
proj_dir=$(yaml_get "$config" job.proj_dir)
data_dir=$(yaml_get "$config" data.data_dir)
n_genes=$(yaml_get "$config" data.n_genes 20)

# inputs from the genome catalog; output dir (mirror model.R get_target)
catalog_dir="${proj_dir}/step1_setup"
genes_info_fp="${catalog_dir}/catalog_genes_info.tsv"
genome_toc_fp="${catalog_dir}/catalog_genome_toc.tsv"
out_dir="${data_dir}/${species_id}/list_of_neighbors"

for f in "$genes_info_fp" "$genome_toc_fp"; do
    [ -s "$f" ] || { echo "ERROR: catalog file missing/empty: $f" >&2
                     echo "       run build_genome_catalog.py $config first." >&2
                     exit 1; }
done

mkdir -p "$out_dir"
outfile="${out_dir}/${query}.tsv"

# idempotency guard: skip if already materialised (-s = exists and non-empty,
# so a truncated TSV from a failed run still gets regenerated)
[[ -s "$outfile" ]] && { echo "skip $query (exists)"; exit 0; }

# every gene member of this centroid_80 (exact match, catalog col 2), paired
# with its genome's .genes file resolved through the genome_toc. Both catalog
# files carry a header row (FNR==1) which we skip on each file.
list_of_genes_with_file=$(awk -F'\t' -v c80="$query" '
    FNR == 1   { next }                            # skip header in both files
    NR == FNR  { toc[$1] = $2; next }              # genome_toc: genome_id -> path
    $2 == c80 {
        # derive genome_id from gene_id: strip the trailing _NNNNN field
        n = split($1, a, "_"); g = a[1];
        for (i = 2; i < n; i++) g = g "_" a[i];
        if (g in toc) print $1 "\t" toc[g];
        else print "WARN: no genome_toc entry for genome " g " (gene " $1 ")" > "/dev/stderr";
    }' "$genome_toc_fp" "$genes_info_fp")

if [ -z "$list_of_genes_with_file" ]; then
    echo "WARN: no gene members for $query — writing empty $outfile" >&2
    : > "$outfile"
    exit 0
fi

# one get_neighbor.sh per (gene member, .genes file), n_genes from the config
echo "$list_of_genes_with_file" | xargs -l -P 4 bash -c \
    "bash \"$script_dir/get_neighbor.sh\" \"\$0\" \"\$1\" $n_genes" > "${outfile}"
