suppressPackageStartupMessages({
  library(ggplot2)
  library(pheatmap)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop(
    paste0(
      "Usage: Rscript 20_transcriptome_descriptive_plots.R ",
      "gene_TPM.tsv immune_panel.tsv gene_map.tsv outdir"
    )
  )
}

tpm_file   <- args[1]
panel_file <- args[2]
map_file   <- args[3]
outdir     <- args[4]

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

########################################
# 1. Read TPM table
########################################

tpm <- read.delim(
  tpm_file,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  row.names = 1
)

if (nrow(tpm) == 0) {
  stop("No genes found in TPM table.")
}

X <- log2(as.matrix(tpm) + 1)

########################################
# 2. PCA
########################################

v <- apply(X, 1, var)

keep <- is.finite(v) & v > 0
Xnz <- X[keep, , drop = FALSE]

if (nrow(Xnz) < 2) {
  stop("Too few variable genes for PCA.")
}

pc <- prcomp(
  t(Xnz),
  center = TRUE,
  scale. = FALSE
)

var_exp <- pc$sdev^2 / sum(pc$sdev^2)

pca_df <- data.frame(
  sample = rownames(pc$x),
  PC1 = pc$x[, 1],
  PC2 = pc$x[, 2],
  stringsAsFactors = FALSE
)

write.table(
  pca_df,
  file.path(outdir, "PCA_scores.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  data.frame(
    PC = paste0("PC", seq_along(var_exp)),
    variance_explained = var_exp
  ),
  file.path(outdir, "PCA_variance.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

p <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    label = sample
  )
) +
  geom_point(size = 3) +
  geom_text(vjust = -0.7) +
  theme_classic() +
  xlab(
    sprintf(
      "PC1 (%.1f%%)",
      100 * var_exp[1]
    )
  ) +
  ylab(
    sprintf(
      "PC2 (%.1f%%)",
      100 * var_exp[2]
    )
  )

ggsave(
  file.path(
    outdir,
    "Figure5A_PCA.pdf"
  ),
  p,
  width = 6,
  height = 5
)

########################################
# 3. Top 100 variable genes
########################################

valid_v <- v[is.finite(v)]

n_top <- min(
  100,
  length(valid_v)
)

top_ids <- names(
  sort(
    valid_v,
    decreasing = TRUE
  )
)[seq_len(n_top)]

top_raw <- X[
  top_ids,
  ,
  drop = FALSE
]

top_scaled <- t(
  scale(
    t(top_raw)
  )
)

top_scaled[
  !is.finite(top_scaled)
] <- 0

write.table(
  top_raw,
  file.path(
    outdir,
    "top100_log2TPM.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

write.table(
  top_scaled,
  file.path(
    outdir,
    "top100_row_scaled.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

pdf(
  file.path(
    outdir,
    "FigureS2_top100_heatmap.pdf"
  ),
  width = 8,
  height = 10
)

pheatmap(
  top_scaled,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "none"
)

dev.off()

########################################
# 4. Read gene ID -> gene symbol map
########################################

gene_map <- read.delim(
  map_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_map <- c(
  "gene_id",
  "gene_name"
)

if (
  !all(
    required_map %in% colnames(gene_map)
  )
) {
  stop(
    "gene_map.tsv must contain gene_id and gene_name columns."
  )
}

gene_map <- gene_map[
  !duplicated(gene_map$gene_id),
  ,
  drop = FALSE
]

rownames(gene_map) <- gene_map$gene_id

########################################
# 5. Read predefined immune / viral-entry panel
########################################

panel <- read.delim(
  panel_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_panel <- c(
  "gene",
  "category"
)

if (
  !all(
    required_panel %in% colnames(panel)
  )
) {
  stop(
    "immune_panel.tsv must contain at least gene and category columns."
  )
}

########################################
# 6. Match panel gene symbols to gene loci
########################################

available_ids <- rownames(X)

map_sub <- gene_map[
  intersect(
    rownames(gene_map),
    available_ids
  ),
  ,
  drop = FALSE
]

hits <- merge(
  data.frame(
    gene_id = rownames(map_sub),
    gene_name = map_sub$gene_name,
    stringsAsFactors = FALSE
  ),
  panel,
  by.x = "gene_name",
  by.y = "gene",
  all = FALSE
)

if (nrow(hits) == 0) {
  warning(
    "No predefined immune or viral-entry genes matched detected loci."
  )
}

########################################
# 7. Add locus suffixes for duplicated symbols
########################################

hits <- hits[
  order(
    hits$gene_name,
    hits$gene_id
  ),
  ,
  drop = FALSE
]

hits$display_name <- hits$gene_name

dup_symbols <- names(
  which(
    table(hits$gene_name) > 1
  )
)

for (g in dup_symbols) {

  idx <- which(
    hits$gene_name == g
  )

  hits$display_name[idx] <- paste0(
    g,
    "-",
    seq_along(idx)
  )
}

write.table(
  hits,
  file.path(
    outdir,
    "immune_panel_detected.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

########################################
# 8. Helper function for heatmap
########################################

make_heatmap <- function(
  hit_df,
  outfile,
  output_table,
  cluster_rows = TRUE,
  min_height = 6,
  row_height = 0.18
) {

  if (nrow(hit_df) == 0) {
    return(invisible(NULL))
  }

  gene_ids <- hit_df$gene_id

  mat <- X[
    gene_ids,
    ,
    drop = FALSE
  ]

  rownames(mat) <- hit_df$display_name

  z <- t(
    scale(
      t(mat)
    )
  )

  z[
    !is.finite(z)
  ] <- 0

  write.table(
    z,
    file.path(
      outdir,
      output_table
    ),
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )

  h <- max(
    min_height,
    row_height * nrow(z)
  )

  pdf(
    file.path(
      outdir,
      outfile
    ),
    width = 8,
    height = h
  )

  pheatmap(
    z,
    cluster_rows = cluster_rows,
    cluster_cols = TRUE,
    scale = "none"
  )

  dev.off()
}

########################################
# 9. Immune-related loci
########################################

immune_hits <- hits[
  hits$category != "Viral-entry",
  ,
  drop = FALSE
]

make_heatmap(
  immune_hits,
  "Figure5B_immune_heatmap.pdf",
  "immune_panel_row_scaled.tsv",
  cluster_rows = TRUE,
  min_height = 6,
  row_height = 0.18
)

########################################
# 10. Viral-entry-associated loci
########################################

entry_hits <- hits[
  hits$category == "Viral-entry",
  ,
  drop = FALSE
]

make_heatmap(
  entry_hits,
  "Figure5C_viral_entry_heatmap.pdf",
  "viral_entry_row_scaled.tsv",
  cluster_rows = FALSE,
  min_height = 4,
  row_height = 0.45
)

########################################
# 11. Summary
########################################

cat(
  "Transcriptomic descriptive analysis completed.\n"
)

cat(
  "Genes in TPM matrix:",
  nrow(tpm),
  "\n"
)

cat(
  "Genes used for PCA:",
  nrow(Xnz),
  "\n"
)

cat(
  "Predefined loci detected:",
  nrow(hits),
  "\n"
)

cat(
  "Immune-related loci:",
  nrow(immune_hits),
  "\n"
)

cat(
  "Viral-entry-associated loci:",
  nrow(entry_hits),
  "\n"
)
