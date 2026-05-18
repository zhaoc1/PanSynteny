#!/bin/bash
# ------------------------------------------------------------------------------
# run_species.sh  <config.yaml>
#
# Step 0 entry point — materialise the per-focal neighbor TSVs that pipeline.R consumes. 
# Fans generate_neighbor_list.sh across every focal centroid listed in gene_list.tsv.
#
# Run order:
#   python  build_genome_catalog.py <config.yaml>  # builds the genome catalog
#   Rscript prepare.R               <config.yaml>  # enumerates missing neighbor TSVs
#   bash    run_species.sh          <config.yaml>  # <-- materialises them
#   Rscript pipeline.R              <config.yaml>  # consumes them
#
# Config-driven: reads the same YAML the R pipeline loads:
#   job.proj_dir                     -> gene_list.tsv path (under step1_setup/; proj_dir is used as-is, no species_id append)
#   data.data_dir + job.species_id   -> {data_dir}/{species_id}/list_of_neighbors (shared species-level cache)
#   job.parallel_jobs                -> xargs -P (required)
# ------------------------------------------------------------------------------

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <config.yaml>"
    exit 1
fi

config="$1"
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
parallel_jobs=$(yaml_get "$config" job.parallel_jobs)   # -P for the xargs fan-out across focals (required; yaml_get errors if absent)

gene_list="${proj_dir}/step1_setup/gene_list.tsv"
out_dir="${data_dir}/${species_id}/list_of_neighbors"

if [[ ! -f "$gene_list" ]]; then
    echo "No gene_list.tsv at $gene_list. Nothing to materialise."
    exit 0
fi

mkdir -p "$out_dir"
total=$(wc -l < "$gene_list")
echo ">>> species_id=$species_id"
echo "    gene_list=$gene_list ($total focal centroids)"
echo "    out_dir=$out_dir"
echo "    parallel_jobs=$parallel_jobs"

# one generate_neighbor_list.sh per focal centroid, parallel_jobs in parallel.
# Capture generate_neighbor_list.sh's per-focal "skip" messages so we can report
# materialised-vs-skipped truthfully at the end.
log_tmp=$(mktemp)
trap 'rm -f "$log_tmp"' EXIT
cat "$gene_list" | xargs -I{} -P "$parallel_jobs" bash -c \
    "bash \"$script_dir/generate_neighbor_list.sh\" \"$config\" \"{}\"" \
    | tee "$log_tmp"

skipped=$(grep -c '^skip ' "$log_tmp" 2>/dev/null || true)
skipped=${skipped:-0}
materialised=$((total - skipped))

if [[ "$materialised" -gt 0 ]]; then
    echo ">>> Done. Materialised $materialised/$total neighbor TSV(s); $skipped already present. Re-run prepare.R / pipeline.R to pick them up."
else
    echo ">>> Done. All $total neighbor TSV(s) were already present: nothing new to materialise."
fi
