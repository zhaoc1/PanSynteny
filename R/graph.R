# ------------------------------------------------------------------------------
# graph.R
#
# Graph-algorithmic support for the strain-aware operon pipeline.
# 
# Author:  Chunyu Zhao <chunyu.zhao@gladstone.ucsf.edu>
# Created: 2025-07-14
# Updated: 2026-04-28
# ------------------------------------------------------------------------------

library(dplyr)
library(igraph)
library(ggraph)
library(tibble)
library(purrr)

longest_common_subseq_len <- function(a, b) {
  m <- length(a)
  n <- length(b)
  dp <- matrix(0, m+1, n+1)
  for (i in 1:m) {
    for (j in 1:n) {
      if (a[i] == b[j]) {
        dp[i+1, j+1] <- dp[i, j] + 1
      } else {
        dp[i+1, j+1] <- max(dp[i+1, j], dp[i, j+1])
      }
    }
  }
  return(dp[m+1, n+1])
}


lcs_similarity <- function(a, b) {
  lcs_len <- longest_common_subseq_len(a, b)
  lcs_len / max(length(a), length(b))
}


generate_path_order <- function(pdf) {
  path_genes <- pdf %>% select(path_label, gene, order) %>%
    mutate(across(where(is.character), trimws)) %>%
    arrange(path_label, order) %>% group_by(path_label) %>% summarise(gene_list = list(gene), .groups = "drop")
  labels <- path_genes$path_label
  n <- length(labels)
  if (n > 1) {
    sim_matrix <- matrix(0, n, n, dimnames = list(labels, labels))
    for (i in 1:n) {
      for (j in 1:n) {
        a <- path_genes$gene_list[[i]]
        b <- path_genes$gene_list[[j]]
        sim_matrix[i, j] <- lcs_similarity(a, b)
      }
    }
    dist_matrix <- as.dist(1 - sim_matrix)
    hc <- hclust(dist_matrix, method = "average")  # or "complete", "ward.D"
    # Step 5: Get ordered path_label list
    ordered_path_labels <- labels[hc$order]
  } else {
    ordered_path_labels <- labels
  }
  return(ordered_path_labels)
}


#' Build directed adjacency edges from a per-path positional neighbor table
#'
#' Convert each focal-gene-centered neighborhood into directed edges between
#' positionally adjacent neighbors. Within each `gene_member` path, neighbors
#' are sorted by the input `position` column (small -> large) and each neighbor
#' is linked to the one immediately following it. The result is a per-path,
#' per-edge table suitable for aggregation into a support-weighted graph and
#' for extracting maximal paths via [get_maximal_paths_by_type()].
#'
#' The edges encode **adjacency in the observed neighborhood**, not
#' chromosomal adjacency: if an intermediate position was `NA` and got
#' filtered out upstream, the resulting edge skips that gap and connects the
#' next two present neighbors.
#'
#' **Edge direction is whatever the caller put in `position`.** The Step 1
#' helper `orient_focal_gene_neighbors()` populates `position` with focal-
#' relative coordinates, so edges read upstream -> downstream relative to the
#' focal gene. The Step 2 driver `stitch_paths_across_focal_genes()`
#' overwrites `position` with chromosomal rank derived from
#' `neighbor_gene_start`, so edges read in chromosomal order centered on the
#' focal. This function does not know which convention is in use; consumers
#' interpreting direction must check the upstream caller.
#'
#' @details
#' **Self-loops are removed.** If the same `neighbor_gene_id` appears at two
#' positions in the same path, the resulting `from == to` edge is dropped
#' before returning. Callers can assume the output is free of self-loops.
#'
#' **Per-path edge counts.** A path with `n` surviving neighbors contributes
#' at most `n − 1` edges (fewer if any are self-loops, which are dropped).
#' Different paths can contribute the same `(from, to)` pair; those rows are
#' kept distinct in the output (differentiated by `gene_member`) so that
#' downstream aggregation can count per-edge support.
#'
#' @export
make_positional_edges <- function(df, assoc = "pos") {
  df %>%
    filter(path_label == assoc) %>%
    arrange(gene_member, position) %>%
    group_by(gene_member) %>%
    mutate(to = lead(neighbor_gene_id)) %>%
    filter(!is.na(to)) %>%
    ungroup() %>% 
    select(focal_c80, gene_member, from = neighbor_gene_id, to, type = path_label) %>%
    filter(from != to) %>%
    distinct()
}


get_edges_from_paths <- function(df, edge_type) {
  # Expand paths into edges
  df %>%
    filter(path_type == edge_type) %>%
    rowwise() %>%
    mutate(c80_genes = strsplit(c80_path_coarse_canonical, " → ")) %>%
    ungroup() %>%
    select(collapsed_path_id, c80_genes) %>% #neighbor_genome
    unnest_wider(c80_genes, names_sep = "_") %>%
    pivot_longer(cols = starts_with("c80_genes_"), values_to = "gene", names_prefix = "c80_genes_") %>%
    group_by(collapsed_path_id) %>% #neighbor_genome
    mutate(next_gene = lead(gene)) %>%
    filter(!is.na(next_gene)) %>%
    ungroup() %>%
    select(from = gene, to = next_gene, collapsed_path_id) #neighbor_genome
}


#' Enumerate maximal source-to-sink paths for one edge type
#'
#' Build a directed graph from `edges`, restrict to edges of a single
#' `edge_type`, and return every maximal path - i.e., every path from an
#' in-degree-0 node (source) to an out-degree-0 node (sink). Within a DAG
#' these are exactly the paths that cannot be extended at either end, so the
#' output captures all complete operon chains observed at this edge type.
#'
#' The graph is assumed to be a DAG; a cycle triggers an error rather than
#' silent truncation, because the depth-first enumeration below has no
#' principled way to choose an edge to break. The function then splits the
#' graph into weakly connected components and enumerates paths per component
#' via depth-first search from each source node.
#'
#' @details
#' **DAG assumption.** The function errors with `"Graph must be a DAG."` if
#' the filtered graph contains any directed cycle. In principle an unusual
#' tandem-duplication pattern could produce a two-cycle (`A → B → A`); in
#' practice positional edges derived from focal-gene-centered neighborhoods
#' are acyclic because `position` is monotonic within each path.
#'
#' **Enumeration is exhaustive.** Every source-to-sink path is returned, not
#' just the longest. In a densely branched DAG the number of paths can grow
#' exponentially in the number of nodes; for the roughly chain-shaped
#' neighborhoods produced upstream this is not a practical concern, but be
#' aware if you ever feed in a lattice-like graph.
#'
#' **Support is not attached here.** Edge-support information (e.g.,
#' `n_support_operon` or `support_genes` from the `edge_support` aggregation)
#' is ignored by this function; only the graph structure is used. To recover
#' per-edge support for a returned path, re-join by `(from, to)` against the
#' original edge table.
#'
#' **Node ID preservation.** Nodes are identified by their `from` / `to`
#' string values in the input `edges`; the function performs no translation.
#' Callers that want to label paths by `neighbor_c80_coarse` (rather than
#' `neighbor_gene_id`) do so after this function returns.
#' @export
get_maximal_paths_by_type <- function(edges, edge_type = "pos") {
  # Step 1: Filter by edge type
  edges <- subset(edges, type == edge_type)
  
  g <- graph_from_data_frame(edges, directed = TRUE)

  # By construction the per-genome graph must be a DAG: edges are emitted by
  # make_positional_edges() in monotone chromosomal order (lead() within each
  # gene_member), self-loops are dropped, and there is no mechanism that could
  # introduce a back-edge from a later position to an earlier one. A cycle here
  # therefore signals a real upstream problem (e.g., a tandem A → B → A in
  # chromosomal order, or a corrupted neighbor TSV) and we hard-stop rather
  # than silently producing nonsense paths - DFS path enumeration has no
  # principled way to break cycles.
  if (!is_dag(g)) stop("Graph must be a DAG.")
  
  # Step 2: Identify weakly connected components
  comps <- components(g, mode = "weak")
  
  path_id <- 1
  path_records <- list()
  
  for (comp in unique(comps$membership)) {
    sub_nodes <- names(comps$membership[comps$membership == comp])
    sub_g <- induced_subgraph(g, vids = sub_nodes)
    
    # Step 3: Get source nodes (in-degree == 0)
    source_nodes <- V(sub_g)[igraph::degree(sub_g, mode = "in") == 0]$name
    
    for (src in source_nodes) {
      # Step 4: DFS from source, collect maximal paths
      dfs_collect <- function(node, visited) {
        visited <- c(visited, node)
        out_nbrs <- neighbors(sub_g, node, mode = "out")$name
        if (length(out_nbrs) == 0) {
          return(list(visited))
        } else {
          paths <- list()
          for (nbr in out_nbrs) {
            if (!(nbr %in% visited)) {
              subpaths <- dfs_collect(nbr, visited)
              paths <- c(paths, subpaths)
            }
          }
          return(paths)
        }
      }
      paths_from_src <- dfs_collect(src, character(0))
      for (p in paths_from_src) {
        path_records[[length(path_records) + 1]] <- data.frame(
          path_id = path_id,
          path_component_id = comp,
          path_length = length(p) - 1,
          path_start = p[1],
          path_end = p[length(p)],
          path_string = paste(p, collapse = " → "),
          stringsAsFactors = FALSE
        )
        path_id <- path_id + 1
      }
    }
  }
  df_paths <- do.call(rbind, path_records)
  return(df_paths)
}


#' Assemble operon paths and edge support across genomes
#'
#' Iterate over each genome present in `gene_neighbors`, build a directed
#' adjacency graph from the positional neighbor data, aggregate per-edge
#' support across all operons in that genome, and enumerate every maximal
#' source-to-sink path through the graph. Results from all genomes are
#' concatenated into two tables - one of paths, one of per-edge support -
#' with node IDs translated into cluster labels at both coarse (cluster-level,
#' `neighbor_c80_coarse`) and fine (length-variant, `neighbor_c80_fine`) resolution.
#'
#' This function is the orchestration layer on top of [make_positional_edges()]
#' and [get_maximal_paths_by_type()]: it wires together the per-genome loop,
#' edge-support aggregation, graph traversal, and c80-label enrichment that
#' would otherwise be duplicated at each caller.
#'
#' @details
#' **Why both resolutions?** Downstream analyses vary in whether they want to
#' treat length variants of the same cluster as distinct (e.g., truncated vs.
#' full-length copies) or collapse them. Keeping both `c80_path_coarse`
#' (coarse) and `c80_path_fine` (fine) in parallel columns means no
#' downstream code has to re-derive either mapping - each consumer picks the
#' resolution that fits the question.
#'
#' **Per-genome scope.** Edges and paths are computed within each genome
#' independently; no cross-genome edges are created. Cross-genome
#' aggregation (e.g., collapsing equivalent paths across strains) happens
#' downstream on the returned tables, not here.
#'
#' **DAG assumption.** Each genome's subgraph is expected to be acyclic;
#' [get_maximal_paths_by_type()] errors out on cycles. Because `position` is
#' monotonic within a focal path, cycles shouldn't arise from well-behaved
#' input.
#'
#' **Per-genome path key for Step 3.** [collapse_paths_across_genomes()]
#' and downstream Step 3 helpers key each per-genome path by
#' `path_genome_comp`, which is not added here. The caller is expected to
#' derive it from the columns in `path_df` before passing downstream:
#' \preformatted{
#' path_df <- path_df \%>\%
#'   mutate(path_genome_comp = paste(path_genome, path_type,
#'                                   path_component_id, sep = "||"))
#' }
#' See [pipeline.R] for the canonical wiring.
#'
#' @export
stitch_paths_across_focal_genes <- function(gene_neighbors) {
  list_of_genomes <- unique(gene_neighbors$gene_member_genome)
  print(paste("Total number of genomes with at least one operon:", length(list_of_genomes)))
  
  list_of_path <- list()
  list_of_esupport <- list()
  for (ge in list_of_genomes) {
    df <- gene_neighbors %>% filter(gene_member_genome == ge)
    
    # all gene_members from the same genome regardless of their Step 1 orientation, 
    # produce edges in the same chromosomal direction
    df <- df %>%
      group_by(gene_member) %>%
      arrange(neighbor_gene_start, .by_group = TRUE) %>%
      mutate(
        anchor_index = match(gene_member, neighbor_gene_id),
        position = row_number() - anchor_index
      ) %>%
      ungroup() %>%
      select(-anchor_index)
    
    # enumerate directed edges between positionally adjacent neighbors,
    # separately for each association class (pos/neg), then stack.
    asso <- unique(df$path_label)
    edge_df <- do.call(rbind, lapply(asso, function(a) {make_positional_edges(df, a)}))
    
    # aggregate: per (from, to, type), count how many distinct focal paths
    # support the edge and record the supporting gene_members.
    edge_support <- edge_df %>%
      group_by(from, to, type) %>%
      summarise(
        n_support_operon = n_distinct(gene_member),
        support_genes = paste(sort(unique(gene_member)), collapse = ";"),
        .groups = "drop"
      )
    
    # walk the support-weighted graph: for each edge type, extract every
    # source-to-sink path (DAG assumption enforced inside the helper).
    paths <- do.call(rbind, lapply(asso, function(a) get_maximal_paths_by_type(edge_support, a) %>% mutate(path_type = a)))
    
    # Keep both resolutions: coarse (neighbor_c80_coarse) is the primary
    # string used downstream; fine (neighbor_c80_fine) is preserved so
    # length-sensitive analyses can be run if needed. Length-variant
    # duplicates would otherwise inflate apparent diversity without signal.
    id_to_c80 <- df %>% select(neighbor_gene_id, neighbor_c80_coarse) %>% distinct() %>% deframe()
    id_to_c80_label <- df %>% select(neighbor_gene_id, neighbor_c80_fine) %>% distinct() %>% deframe()
    
    # build both path-string views by mapping each node in the path
    paths$c80_path_coarse <- sapply(paths$path_string, function(path) {
      ids <- strsplit(path, " → ")[[1]]
      paste(id_to_c80[ids], collapse = " → ")
    })
    paths$c80_path_fine <- sapply(paths$path_string, function(path) {
      ids <- strsplit(path, " → ")[[1]]
      paste(id_to_c80_label[ids], collapse = " → ")
    })
    
    # endpoint lookups (coarse only; label versions are the first/last elements of c80_path_fine if needed)
    paths$c80_start <- id_to_c80[paths$path_start]
    paths$c80_end <- id_to_c80[paths$path_end]
    
    # edge-level support expressed as c80-label strings (both resolutions)
    edge_support <- edge_support %>%
      mutate(
        support_c80s = sapply(support_genes, function(gene_str) {
          gene_ids <- unlist(strsplit(gene_str, ";"))
          paste(id_to_c80[gene_ids], collapse = ";")
        }),
        support_c80_labels = sapply(support_genes, function(gene_str) {
          gene_ids <- unlist(strsplit(gene_str, ";"))
          paste(id_to_c80_label[gene_ids], collapse = ";")
        })
      )
    
    # per-genome provenance columns (reading genome from the outer loop var)
    paths$path_genome <- ge
    edge_support$support_genome <- ge
    
    # put the most useful columns first for inspection
    paths <- paths %>%
      select(path_genome, path_id:path_length, path_type, c80_path_coarse, c80_path_fine, c80_start, c80_end, everything())
    
    list_of_path[[ge]] <- paths
    list_of_esupport[[ge]] <- edge_support
  }
  
  # concatenate across genomes; .id column holds the genome identifier
  path_df <- bind_rows(list_of_path, .id = "neighbor_genome")
  esupport_df <- bind_rows(list_of_esupport, .id = "neighbor_genome")
 
  list(path_df = path_df, esupport_df = esupport_df) 
}


# Function to compute longest contiguous overlap
max_overlap <- function(a, b) {
  max_len <- 0
  for (i in seq_along(a)) {
    for (j in seq_along(b)) {
      k <- 0
      while ((i + k) <= length(a) && (j + k) <= length(b) && a[i + k] == b[j + k]) {
        k <- k + 1
      }
      if (k > max_len) max_len <- k
    }
  }
  return(max_len)
}


#' Collapse per-genome maximal paths into a non-redundant cross-genome set
#'
#' Reduce the per-genome path table produced by
#' [stitch_paths_across_focal_genes()] to a non-redundant set where rows that
#' share an identical coarse-label path string (and edge `type`) collapse into
#' one, annotated with the genomes that contributed and with pointers back to
#' the per-genome path instances. This is the first step of cross-genome
#' consolidation; direction canonicalization (forward vs. reverse) is
#' deferred to [generate_canonical_path()].
#'
#' @details
#' **Collapse key is the raw coarse string.** Grouping is on
#' `c80_path_coarse` verbatim - no direction canonicalization is applied
#' here. Because Step 2 sorts by chromosomal `neighbor_gene_start` and
#' discards Step 1's per-focal orientation, two genomes carrying the same
#' operon on opposite strands emit mirror-image strings
#' (`A → B → C` vs. `C → B → A`) and remain as two separate
#' `collapsed_path_id` rows. Forward/reverse unification happens one step
#' later in [generate_canonical_path()], which applies [normalize_path()]
#' per row.
#'
#' **Resolution.** Collapse is on `c80_path_coarse` (coarse cluster-level
#' labels). The length-variant `c80_path_fine` column, if present,
#' is not used and is not propagated downstream of this function - length
#' variants of the same coarse cluster string therefore collapse here.
#'
#' **`path_length` in the grouping is redundant** (fully implied by
#' `c80_path_coarse`), kept in the key for schema clarity.
#'
#' @export
collapse_paths_across_genomes <- function(path_df) {
  # non-redundant set of paths, annotated with which genomes they appear in, how many genomes, and their genome-specific instances.
  collapsed_paths <- path_df %>%
    select(neighbor_genome, path_id, path_length, c80_path_coarse, path_type, path_genome_comp) %>%
    group_by(c80_path_coarse, path_length, path_type) %>%
    summarise(
      n_genomes = n_distinct(neighbor_genome),
      neighbor_genomes = paste(sort(unique(neighbor_genome)), collapse = ";"),
      per_genome_path_w_ids = paste(sort(unique(path_genome_comp)), collapse = ";"), 
      .groups = "drop"
    ) %>%
    mutate(collapsed_path_id = paste0("Path_", row_number())) %>% 
    select(collapsed_path_id, everything())
  return(collapsed_paths)
}


# Build the cleaned token vector used for orientation decisions:
# strip synthetic small-ORF tokens (`_`-prefixed; produced by
# compute_short_gene_prevalence in midas.R), then collapse consecutive
# duplicates (a b b c -> a b c). Used by normalize_path and should_flip.
# The chosen direction is applied to the original full vector, so synthetic
# ORFs and adjacent duplicates remain in the stored canonical.
clean_for_orientation <- function(tokens) {
  tokens <- tokens[!startsWith(tokens, "_")]
  if (length(tokens) <= 1) return(tokens)
  tokens[c(TRUE, tokens[-1] != tokens[-length(tokens)])]
}


#' Normalize a path string to a canonical direction
#'
#' Decide a deterministic direction for `path_string` (forward vs. reversed)
#' by lex-comparing the two orientations of a *cleaned* token vector
#' (synthetic small-ORF tokens stripped, adjacent duplicates collapsed),
#' and return the chosen direction applied to the **full** original tokens.
#'
#' **Orientation-only contract.** This function only ever returns either
#' `forward_full` or `reverse_full` of the input - it never reorders tokens
#' within a direction, never drops tokens from the stored output, and never
#' changes token content. Downstream callers
#' ([explode_canonical_into_collapsed_paths()] and the per-isoform / per-
#' genome expansion functions in path.R) rely on this invariant so that one
#' coarse-string flip decision drives flips at fine resolution and
#' gene-id resolution via a simple `rev()`. **Do not introduce
#' token-reordering operations here without re-deriving fine-level
#' direction logic in the callers.**
normalize_path <- function(path_string) {
  # Split into vector
  genes <- strsplit(path_string, " → ", fixed = TRUE)[[1]]
  forward_full <- paste(genes, collapse = " → ")
  reverse_full <- paste(rev(genes), collapse = " → ")

  cleaned <- clean_for_orientation(genes)
  # Degenerate (no real-gene backbone after cleaning, or everything collapsed
  # into one token e.g. palindromic-cleaned). Fall back to full-string lex-min
  # so output stays deterministic. Rare for operons >= 5 genes.
  if (length(cleaned) <= 1) {
    return(if (forward_full < reverse_full) forward_full else reverse_full)
  }
  forward_cleaned <- paste(cleaned, collapse = " → ")
  reverse_cleaned <- paste(rev(cleaned), collapse = " → ")
  if (forward_cleaned < reverse_cleaned) forward_full else reverse_full
}


split_path_string <- function(path_str) {
  # Split on the " → " arrow separator. Returns character(0) for NA/empty input
  # so downstream `rev()` / indexing stays safe.
  if (is.na(path_str) || length(path_str) == 0) character(0)
  else strsplit(path_str, " → ", fixed = TRUE)[[1]]
}


#' Unify forward and reverse path renderings into canonical operon identities
#'
#' Close the direction-ambiguity loop left open by
#' [collapse_paths_across_genomes()]: apply [normalize_path()] to every row so
#' that a path and its exact reverse share one string, assign a surrogate
#' `canonical_path_id` keyed on that normalized string, aggregate genome
#' support across the direction-mirror rows, and apply the final minimum-
#' support gate. This is where Step 2's orientation discard is finally
#' canonicalized.
#' 
#' @details
#' **Mixed-type sanity check.** The function warns if any `canonical_path_id`
#' spans more than one `type` value (e.g., a `pos` path whose reverse was
#' recorded as `neg`). The diagnostic data frame inside the `if` block is
#' built but not assigned or printed - only the `warning()` surfaces.
#' Execution continues; if the invariant is actually violated, the
#' subsequent `summarise(type = unique(type))` aborts with a length-mismatch
#' error rather than with the warned groups.
#'
#' **`n_genomes` double-count.** Summing `n_genomes` across forward + reverse
#' rows double-counts any genome that happens to carry *both* directions of
#' the same operon. In practice operons are strand-bound within a genome, so
#' this is rare; the column is not guaranteed to equal
#' `n_distinct(neighbor_genome)` across the underlying `path_df`.
#'
#' **`neighbor_genomes` dedup granularity.** `unique()` is applied to the
#' pre-`;`-joined strings inherited from upstream, not to individual genome
#' IDs. If forward and reverse rows contribute `"A;B"` and `"B;C"` the output
#' is `"A;B;B;C"`. Re-split and dedupe downstream if a gene-level list is
#' needed.
#'
#' **Mixed separators, inherited.** `collapsed_path_id` is comma-separated;
#' `neighbor_genomes` is semicolon-separated. Not harmonized here.
#'
#' **Scope vs. [orient_paths_within_component()].** This function provides
#' *identity* canonicalization only (one surrogate per path+reverse pair).
#' Consistent orientation *within a joint component* - so that paths sharing
#' a subpath align left-to-right - is the job of
#' [orient_paths_within_component()] and is applied later.
#'
#' @export
generate_canonical_path <- function(collapsed_paths, path_min_genomes) {
  # Add normalized path column
  collapsed_paths$c80_path_coarse_canonical <- vapply(collapsed_paths$c80_path_coarse, normalize_path, character(1))
  
  # Assign new collapsed IDs (e.g., cp_1, cp_2, ...)
  collapsed_paths$canonical_path_id <- paste0("cp_", match(collapsed_paths$c80_path_coarse_canonical, unique(collapsed_paths$c80_path_coarse_canonical)))
  
  # Make sure there is no mix
  if (nrow(collapsed_paths %>% group_by(canonical_path_id) %>% filter(n_distinct(path_type) > 1) %>% ungroup()) > 0) {
    warning("some component has mixed path")
    collapsed_paths %>%
      group_by(canonical_path_id) %>% filter(n_distinct(path_type) > 1) %>% ungroup() %>%
      select(collapsed_path_id, canonical_path_id, n_genomes)
  }

  canonical_paths <- collapsed_paths %>%
    group_by(canonical_path_id) %>%
    summarise(
      path_type = unique(path_type),
      collapsed_path_id = paste(unique(collapsed_path_id), collapse = ","),
      c80_path_coarse_canonical = unique(c80_path_coarse_canonical),
      n_genomes = sum(n_genomes, na.rm = TRUE),
      neighbor_genomes = paste(unique(neighbor_genomes), collapse = ";")
    ) %>%
    ungroup() %>% 
    select(canonical_path_id, everything())
  
  canonical_paths <- canonical_paths %>% filter(n_genomes >= path_min_genomes) #<--- bare minimal
}


#' Group canonical paths into cross-type connected components at gene level
#'
#' Build the gene-level adjacency graph spanning all canonical paths (across
#' every edge type, with direction and type information discarded) and
#' compute its connected components. Each component is a maximal set of
#' `centroid_80` gene nodes that are transitively linked by observed
#' adjacency anywhere in the corpus - a candidate "genomic locus identity"
#' that unifies all canonical-path variants touching that gene neighborhood.
#'
#' This scaffold is the grouping scope for downstream operations:
#' [decorate_paths_with_components()] projects component membership back onto
#' each canonical path, [orient_paths_within_component()] uses the component
#' as the reference frame for within-component direction alignment, and
#' component-level aggregation is the unit for focal-block analysis.
#'
#' @details
#' **Edge type is intentionally discarded.** Per-type edges are collected
#' then collapsed to distinct `(from, to)` pairs before graph construction.
#' A gene that participates in both a `pos` and a `neg` path therefore
#' bridges those subgraphs into one component. The justification (see
#' inline comment at the top of the function) is that bidirectional cycles
#' in merged operon paths are most likely artifacts rather than real
#' biological signal, so treating the gene-adjacency relation as undirected
#' and type-agnostic is the correct identity for locus grouping.
#'
#' **Undirected graph.** `graph_from_data_frame(..., directed = FALSE)`
#' means forward and reverse renderings of the same adjacency unify
#' automatically, and `components(mode = "weak")` returns the same result
#' as `"strong"` on an undirected graph.
#'
#' **Hard stop on empty input.** If no edges are produced (empty
#' `canonical_paths` or a filter that rejects every row), the function
#' `stop()`s rather than returning an empty map. Callers should guard
#' upstream.
#'
#' **Resolution is coarse.** Nodes are the `centroid_80` cluster labels
#' embedded in `c80_path_coarse_canonical`. Length-variant `c80_path_fine`
#' information is not considered here - length isoforms always co-compose
#' with their parent cluster.
#'
#' **Hub-gene merging risk.** A c80 gene that participates in many
#' unrelated operons (e.g., a promiscuous regulator or a ribosomal gene)
#' bridges those operons into one megacomponent. The undirected +
#' type-collapsed choice amplifies this. No current code flags such
#' components; a manual sanity check on component size is advised if
#' results look unexpectedly merged.
#'
#' @export
compute_joint_components <- function(canonical_paths, edge_types = c("pos", "neg")) {
  # 1) Collect edges
  # # Bidirectional edges forming cycles in your merged operon paths are most 
  # likely artifacts, not biologically valid operons.
  etypes <- as.character(edge_types)
  combined_edges <- do.call(rbind, lapply(etypes, function(a) get_edges_from_paths(canonical_paths, a)))
  
  if (is.null(combined_edges) || nrow(combined_edges) == 0) {
    stop("No edges were produced. Check edge_types and get_edges_from_paths")
  }
  
  #### Build joint graph from all edges for connected components
  # Use all edges (both pos and neg), drop type info for now
  joint_edges <- combined_edges %>% select(from, to) %>% distinct()
  
  # build un-directed graph
  joint_graph <- igraph::graph_from_data_frame(joint_edges, directed = FALSE)
  
  # Compute connected components
  joint_comps <- components(joint_graph, mode = "weak")$membership
  
  # Map each node to its joint_component_id
  joint_component_map <- data.frame(node = names(joint_comps), 
                                    joint_component_id = joint_comps, stringsAsFactors = FALSE)
  return(joint_component_map)
}


#' Project gene-level joint-component membership onto canonical paths
#'
#' For each canonical path, determine which joint component(s) its genes
#' belong to and attach the result as a single string column. This is the
#' path-level projection of the node-level membership map produced by
#' [compute_joint_components()], and it supplies the grouping key that
#' [orient_paths_within_component()] uses next.
#'
#' @details
#' **Column is a string, not an integer.** `"3"` and `"3;7"` are both valid
#' values; numeric comparison does not apply.
#'
#' **Unmapped genes become literal `"NA"`.** A gene absent from the
#' component map (e.g., a gene that appeared only in a single-node path so
#' produced no edges in [compute_joint_components()]) joins as `NA_integer_`
#' and is rendered as the string `"NA"` by `paste(...)`. Downstream
#' `is.na()` does not catch these; string equality `== "NA"` does.
#'
#' **Multi-component paths are rare in practice.** Because every adjacent
#' pair in a canonical path contributes an edge in
#' [compute_joint_components()], all genes of a path are guaranteed
#' co-component by graph connectivity. Multi-value strings arise only when
#' a path's genes were lost from the map (unusual) or when the input map
#' was computed from a different edge set.
#'
#' **Cross-type shared nodes do not split the path.** A c80 that appears in
#' both `pos` and `neg` operons cannot push `joint_component_ids` to a
#' multi-value string, because [compute_joint_components()] drops type info
#' before building the undirected graph. The opposite happens: a shared
#' node fuses every operon it touches into one component, so all paths
#' through it inherit the same single id. The string stays single-valued,
#' but the value can label a much larger blob than any one operon - see
#' the megacomponent caveat in [compute_joint_components()].
#'
#' @export
decorate_paths_with_components <- function(canonical_paths, joint_component_map) {
  # Split gene_path into individual genes
  canonical_paths_genes <- canonical_paths %>%
    select(canonical_path_id, c80_path_coarse_canonical) %>%
    mutate(gene = strsplit(c80_path_coarse_canonical, " → ")) %>%
    unnest(gene)
  
  # Map each gene to its graph/component ID.
  canonical_paths_genes <- canonical_paths_genes %>%
    left_join(joint_component_map, by = c("gene" = "node"))
  
  # Summarise, per path, which components it covers.
  path_component_map <- canonical_paths_genes %>%
    group_by(canonical_path_id) %>%
    summarise(joint_component_ids = paste(sort(unique(joint_component_id)), collapse = ";"))
  
  canonical_paths <- canonical_paths %>% 
    left_join(path_component_map, by=c("canonical_path_id")) %>%
    select(joint_component_ids, canonical_path_id, collapsed_path_id, path_type, everything())
}


#' Decide whether a path should be flipped to align with a reference
#'
#' Compute the longest contiguous-substring overlap (via
#' [max_overlap()]) between the reference and the target in both
#' orientations of the target, and return `TRUE` when the reversed
#' target scores strictly higher than the forward target.
#'
#' Both inputs pass through [clean_for_orientation()] before scoring,
#' so synthetic small-ORF tokens and adjacent duplicates do not
#' contribute to the decision. [max_overlap()] is symmetric in its
#' two arguments, so the result is invariant to argument order.
#'
#' **Forward bias on ties.** The comparison is strict `>`, so when
#' `rev_score == fwd_score` the function returns `FALSE` and the
#' target is left as-is. The return value alone does not distinguish
#' a confident "no flip" from an arbitrary tied "no flip".
#'
#' **Tie cases (return `FALSE`).**
#' \itemize{
#'   \item Reference contains a palindromic substring that matches
#'     the target equally in both directions.
#'   \item Target shares zero tokens with the reference (peripheral
#'     path; reaches the reference only transitively through other
#'     paths in the same component).
#'   \item Target shares only one token with the reference - a
#'     single token matches itself in either direction.
#' }
#' Empty cleaned vectors short-circuit to `FALSE`.
#'
#' @return `TRUE` if `rev(target)` aligns more with `reference` than
#'   `target` does; `FALSE` otherwise (including all ties).
should_flip <- function(reference, target) {
  ref_c <- clean_for_orientation(reference)
  tgt_c <- clean_for_orientation(target)
  if (length(ref_c) == 0L || length(tgt_c) == 0L) return(FALSE)
  fwd_score <- max_overlap(tgt_c, ref_c)
  rev_score <- max_overlap(rev(tgt_c), ref_c)
  rev_score > fwd_score
}


#' Return the best-oriented version of a query relative to a reference
#'
#' Decide-and-apply wrapper around [should_flip()]: when the query
#' should be flipped to align with the reference, return the reversed
#' query; otherwise return the query unchanged. The flip is applied
#' to the **full** query token vector - synthetic ORFs and adjacent
#' duplicates survive in the returned data even though they did not
#' drive the decision (cleaning happens inside `should_flip` for the
#' decision only).
best_orient <- function(ref, qry) {
  if (should_flip(qry, ref)) {
    return(rev(qry))
  } else {
    return(qry)
  }
}

collapse_path <- function(gene_list) paste(gene_list, collapse = " → ")


#' Align path directions within each joint component
#'
#' Within each joint component, choose one path as the direction reference
#' (the longest path in the component) and flip every other path to the
#' orientation that maximizes its contiguous-substring overlap against the
#' reference. The result is a set of canonical paths that, within a
#' component, all read left-to-right in a mutually consistent chromosomal
#' direction - which is what downstream visualization, sub-pattern
#' detection, and focal-block analysis assume.
#'
#' This function does **not** aim for a globally meaningful direction: the
#' reference's own orientation is inherited from [normalize_path()]'s
#' lexicographic rule, which is biology-agnostic. Only relative consistency
#' within a component is guaranteed.
#'
#' **Orientation-only contract.** Like [normalize_path()], this function only
#' ever picks between `forward` and `rev(forward)` of each path's token
#' vector - never reorders tokens within a direction or changes token
#' content. Downstream callers
#' ([explode_canonical_into_collapsed_paths()] and the per-isoform / per-
#' genome expansion functions in path.R) rely on this so the
#' `needs_flip` boolean computed on coarse strings can drive fine-level
#' and gene-id-level flips via a simple `rev()`. Future modifiers must
#' preserve this invariant.
#'
#' @details
#' **Algorithm.** For each distinct value of `joint_component_ids`:
#' \enumerate{
#'   \item Pick the longest path in the group as reference (`which.max`,
#'     first-max tie-break).
#'   \item For every path in the group, call
#'     [best_orient()]\code{(ref_path, p)}, which uses [should_flip()] to
#'     compare forward vs. reverse [max_overlap()] against the reference and
#'     reverses the path if the reverse scores strictly higher.
#'   \item Write the oriented token vectors back into the group, then into
#'     `dd`.
#' }
#' After all groups are processed, `c80_path_coarse_canonical` is rebuilt from the
#' oriented token vectors via [collapse_path()].
#'
#' **Edge cases with practical impact.**
#' \describe{
#'   \item{Longest-path tie}{`which.max` returns the first index. If two
#'     paths share the maximum length and point in opposite directions, the
#'     chosen reference depends on input row order. A different reference
#'     produces a whole-component mirror flip. Correctness within the
#'     component is unaffected - only the absolute direction differs across
#'     runs.}
#'   \item{Zero overlap between query and reference}{A peripheral path that
#'     connects to the reference only transitively (through an intermediate
#'     path) can have `fwd_score == rev_score == 0`. The strict `>` in
#'     [should_flip()] keeps it forward by default, so its orientation is
#'     chosen arbitrarily rather than by its actual relationship to its
#'     non-reference neighbors.}
#'   \item{Palindromic overlap}{If the reference contains a substring whose
#'     token sequence is its own reverse (e.g., `[A, B, C, B, A]`), a query
#'     matching that region scores equally forward and reverse. The tie
#'     keeps it forward.}
#'   \item{Single-gene bridge}{A query that shares only one gene with the
#'     reference (or with any other oriented path) is fundamentally
#'     direction-blind under a contiguous-substring metric: one token
#'     matches itself in either direction. No fix inside this function can
#'     recover direction from a single-gene bridge.}
#' }
#'
#' **Benign cases.** A singleton component, the reference compared to
#' itself, and a path whose `joint_component_ids` is a multi-component
#' string all reduce to no-ops (no flip) and produce correct output.
#'
#' @export
orient_paths_within_component <- function(dd) {
  dd$splited_path_string <- lapply(dd$c80_path_coarse_canonical, split_path_string)

  # Align all paths within each component
  dd$oriented_path <- vector("list", nrow(dd))
  grouped <- split(dd, dd$joint_component_ids)

  for (gid in names(grouped)) {
    group <- grouped[[gid]]
    # Pick reference by *cleaned* length so a path padded with synthetic
    # small-ORF tokens or fragmentation duplicates does not dominate
    # orientation. Tie behavior (first-index) is unchanged. ref_path stays
    # the FULL token vector at the chosen index.
    lengths <- sapply(group$splited_path_string, function(p) length(clean_for_orientation(p)))
    ref_idx <- which.max(lengths)
    ref_path <- group$splited_path_string[[ref_idx]]
    group$oriented_path <- lapply(group$splited_path_string, function(p) best_orient(ref_path, p))
    dd$oriented_path[dd$joint_component_ids == gid] <- group$oriented_path
  }
  dd$c80_path_coarse_canonical <- sapply(dd$oriented_path, collapse_path)
  dd <- dd %>% select(-one_of(c("splited_path_string", "oriented_path")))

  return(dd)
}


#' Run Step 2 - per-genome path stitching
#'
#' Orchestrator for Step 2. Calls [stitch_paths_across_focal_genes()] to
#' build per-genome maximal paths from `gene_neighbors`, derives the
#' composite per-genome path key `path_genome_comp` (the join column Step
#' 3's cross-genome aggregation needs), drops the now-redundant
#' `path_genome` and `path_component_id` columns, and persists the two
#' Step 2 caches: `path_df.rds` (per-genome paths) and `esupport_df.rds`
#' (per-edge focal support).
#'
#' Idempotent: if the `path_df` cache already exists on disk, the
#' build/save block is skipped and the cached frame is returned. To force
#' a re-run, delete the file at `get_target("path_df")`.
#'
#' @param gene_neighbors Output of Step 1 - one row per (focal, genome,
#'   neighbor position), with `neighbor_c80_coarse` and
#'   `neighbor_c80_fine` populated.
#'
#' @return The `path_df` frame (one row per per-genome maximal path).
#'   On cache hit, returned from `readRDS`; on cache miss, returned in
#'   memory after `saveRDS`. The companion `esupport_df` is persisted
#'   but not returned (Step 3 does not consume it; ad-hoc inspection
#'   can `readRDS(get_target("esupport_df"))`).
#'
#' @export
run_step2_path_stitching <- function(gene_neighbors) {
  path_df_rds <- get_target("path_df")
  if (file.exists(path_df_rds)) {
    return(readRDS(path_df_rds))
  }

  res <- stitch_paths_across_focal_genes(gene_neighbors)

  path_df <- res$path_df %>%
    mutate(path_genome_comp = paste(path_genome, path_type, path_component_id, sep = "||")) %>%
    select(-one_of(c("path_genome", "path_component_id")))

  saveRDS(path_df, get_target("path_df"))
  saveRDS(res$esupport_df, get_target("esupport_df"))

  path_df
}
