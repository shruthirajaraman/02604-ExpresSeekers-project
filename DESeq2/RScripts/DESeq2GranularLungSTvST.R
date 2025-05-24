# Load required libraries
library(DESeq2)
library(tidyverse)
library(ggrepel)
library(biomaRt)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Read count data
#count_data <- read.csv("DESeq2Input/lungFilteredGenesCounts.csv", header = TRUE, row.names = 1, check.names = FALSE)
count_data <- read.csv("DESeq2Input/subsetted_lung_counts_ml_model.csv", header = TRUE, row.names = 1, check.names = FALSE)

# Read metadata file
#metadata <- read.csv("DESeq2Input/lung_DiseaseType_Metadata.csv", header=TRUE, row.names=1)
metadata <- read.csv("DESeq2Input/lung_subset_DiseaseType_Metadata.csv", header=TRUE, row.names=1)

# Ensure column names in count data match row names in metadata
if (!all(colnames(count_data) %in% rownames(metadata))) {
  stop("Mismatch between count data columns and metadata row names!")
}

# Define subtype comparisons
#comparison_pairs <- list(
 # c("squamous_cell_neoplasms", "adenomas_and_adenocarcinomas"),
  #c("squamous_cell_neoplasms", "acinar_cell_neoplasms"),
  #c("adenomas_and_adenocarcinomas", "acinar_cell_neoplasms")
#)
comparison_pairs <- list(
  c("squamous_cell_neoplasms", "adenomas_and_adenocarcinomas")
)

label_map <- c(
  squamous_cell_neoplasms = "Squamous Cell Neoplasms",
  adenomas_and_adenocarcinomas = "Adenomas and Adenocarcinomas",
  acinar_cell_neoplasms = "Acinar Neoplasms"
)

table(metadata$condition)

summary_list <- list()

# Loop
for (pair in comparison_pairs) {
  group1 <- pair[1]
  group2 <- pair[2]
  
  comp_label <- paste(group1, "vs", group2, sep = "_")
  plot_label_1 <- label_map[[group1]]
  plot_label_2 <- label_map[[group2]]
  plot_title <- paste(plot_label_1, "vs", plot_label_2)
  
  message("Running comparison: ", comp_label)
  print(table(metadata$condition)[c(group1, group2)])
  
  subset_meta <- metadata[metadata$condition %in% c(group1, group2), , drop = FALSE]
  subset_counts <- count_data[, rownames(subset_meta)]
  subset_meta$condition <- factor(subset_meta$condition, levels = c(group1, group2))
  
  # Check sample count
  if (any(table(subset_meta$condition) < 4)) {
    message("Skipping ", comp_label, " — not enough samples.")
    next
  }
  
  # Create output directory
  out_dir <- paste0("results/", comp_label, "/")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # DESeq2 workflow
  dds <- DESeqDataSetFromMatrix(countData = subset_counts,
                                colData = subset_meta,
                                design = ~ condition)
  dds <- DESeq(dds)
  res <- results(dds, alpha = 0.05)
  coef_name <- resultsNames(dds)[grepl(group2, resultsNames(dds))]
  if (length(coef_name) == 0) {
    message("Skipping ", comp_label, " — no matching coefficient found.")
    next
  }
  
  res_shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm")
  sig_genes <- subset(res_shrunk, padj < 0.05 & abs(log2FoldChange) > 1)
  message("Significant DEGs: ", nrow(sig_genes))
  
  up_genes <- sum(sig_genes$log2FoldChange > 0, na.rm = TRUE)
  down_genes <- sum(sig_genes$log2FoldChange < 0, na.rm = TRUE)
  top_gene <- rownames(sig_genes)[which.min(sig_genes$padj)[1]]
  
  summary_list[[comp_label]] <- data.frame(
    comparison = comp_label,
    num_sig_genes = nrow(sig_genes),
    upregulated = up_genes,
    downregulated = down_genes,
    top_gene = top_gene
  )
  
  # Save results
  write.csv(sig_genes[order(sig_genes$padj, na.last = NA), ],
            file = paste0(out_dir, "sig_genes.csv"))
  write.csv(as.data.frame(res[order(res$padj, na.last = NA), ]),
            file = paste0(out_dir, "all_genes_results.csv"))
  
  # MA plot
  png(paste0(out_dir, "MA_plot.png"), width = 1000, height = 800, res = 150)
  plotMA(res, main = paste("MA Plot:", plot_title), ylim = c(-5, 5))
  abline(h = c(-1, 1), col = "blue", lty = 2)
  dev.off()
  
  # Gene IDs to gene names
  ensembl_ids <- gsub("\\..*", "", rownames(res_df))
  
  # Get mapping
  mapping <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = ensembl_ids,
    mart = ensembl
  )
  res_df <- as.data.frame(res_shrunk)
  
  res_df$ensembl_gene_id <- gsub("\\..*", "", rownames(res_df))
  
  # Merge
  res_annotated <- merge(res_df, mapping, by = "ensembl_gene_id", all.x = TRUE)
  
  # Use gene name (HGNC symbol) where available; fallback to Ensembl ID
  res_annotated$gene_label <- ifelse(res_annotated$hgnc_symbol != "",
                                     res_annotated$hgnc_symbol,
                                     res_annotated$ensembl_gene_id)
  # Custom MA plot
  res_df <- as.data.frame(res_annotated)
  res_df$mean <- res_df$baseMean
  res_df$sig <- with(res_df, ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "Significant", "NS"))
  
  res_annotated$mean <- as.numeric(res_annotated$baseMean)
  res_annotated$sig <- ifelse(res_annotated$padj < 0.05 & abs(res_annotated$log2FoldChange) > 1,
                              "Significant", "NS")
  
  gg_ma <- ggplot(res_annotated, aes(x = log10(mean + 1), y = log2FoldChange)) +
    geom_point(aes(color = sig), alpha = 0.5) +
    scale_color_manual(values = c("Significant" = "blue", "NS" = "gray")) +
    geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "blue") +
    ylim(-5, 5) +
    labs(title = paste("MA Plot:", plot_title),
         x = "log10 Mean Normalized Count",
         y = "Log2 Fold Change") +
    theme_minimal()

  top_genes_ma <- res_annotated %>%
    filter(padj < 0.05 & !is.na(log2FoldChange) & abs(log2FoldChange) <= 5) %>%
    arrange(padj) %>%
    slice(1:10)
  
  gg_ma <- gg_ma +
    geom_text_repel(data = top_genes_ma,
                    aes(label = gene_label),
                    max.overlaps = 10,
                    size = 3,
                    box.padding = 0.3,
                    point.padding = 0.2,
                    segment.color = "grey50")

  
  ggsave(paste0(out_dir, "MA_plot_custom.png"), plot = gg_ma, width = 7.5, height = 6)
  
  # Unshrunk MA plot
  
  res_df <- as.data.frame(res)
  
  res_df$ensembl_gene_id <- gsub("\\..*", "", rownames(res_df))
  
  # Merge
  res_annotated <- merge(res_df, mapping, by = "ensembl_gene_id", all.x = TRUE)
  
  # Use gene name (HGNC symbol) where available; fallback to Ensembl ID
  res_annotated$gene_label <- ifelse(res_annotated$hgnc_symbol != "",
                                     res_annotated$hgnc_symbol,
                                     res_annotated$ensembl_gene_id)
  # Custom MA plot
  res_df <- as.data.frame(res_annotated)
  res_df$mean <- res_df$baseMean
  res_df$sig <- with(res_df, ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "Significant", "NS"))
  
  res_annotated$mean <- as.numeric(res_annotated$baseMean)
  res_annotated$sig <- ifelse(res_annotated$padj < 0.05 & abs(res_annotated$log2FoldChange) > 1,
                              "Significant", "NS")
  
  gg_ma <- ggplot(res_annotated, aes(x = log10(mean + 1), y = log2FoldChange)) +
    geom_point(aes(color = sig), alpha = 0.5) +
    scale_color_manual(values = c("Significant" = "blue", "NS" = "gray")) +
    geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "blue") +
    ylim(-5, 5) +
    labs(title = paste("MA Plot:", plot_title),
         x = "log10 Mean Normalized Count",
         y = "Log2 Fold Change") +
    theme_minimal()
  
  top_genes_ma <- res_annotated %>%
    filter(padj < 0.05 & !is.na(log2FoldChange) & abs(log2FoldChange) <= 5) %>%
    arrange(padj) %>%
    slice(1:10)
  
  gg_ma <- gg_ma +
    geom_text_repel(data = top_genes_ma,
                    aes(label = gene_label),
                    max.overlaps = 10,
                    size = 3,
                    box.padding = 0.3,
                    point.padding = 0.2,
                    segment.color = "grey50")
  
  
  ggsave(paste0(out_dir, "MA_plot_unshrunk.png"), plot = gg_ma, width = 7.5, height = 6)
  
  
  
  # Summary file
  sink(paste0(out_dir, "summary.txt"))
  summary(res)
  sink()
  
  # Volcano plot
  res_df <- as.data.frame(res_shrunk)
  res_df$gene <- rownames(res_df)
  res_df$volcano_color <- "NS"
  res_df$volcano_color[res_df$padj < 0.05 & res_df$log2FoldChange > 0] <- "Up"
  res_df$volcano_color[res_df$padj < 0.05 & res_df$log2FoldChange < 0] <- "Down"
  
  top_genes <- res_df[order(res_df$padj), ][1:10, ]
  
  volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = volcano_color), alpha = 0.7) +
    geom_text_repel(data = top_genes,
                    aes(label = gene),
                    max.overlaps = 10,
                    size = 3,
                    box.padding = 0.3,
                    point.padding = 0.2,
                    segment.color = "grey50") +
    scale_color_manual(values = c("NS" = "black", "Up" = "red", "Down" = "blue")) +
    theme_minimal() +
    ggtitle(paste("Volcano Plot:", plot_title)) +
    xlab("Log2 Fold Change") +
    ylab("-log10 Adjusted p-value")
  
  ggsave(paste0(out_dir, "volcano_plot.png"),
         plot = volcano, width = 7.47, height = 6.22)
}

# Save summary table
summary_df <- bind_rows(summary_list)
write.csv(summary_df, file = "results/subtype_DEG_summary_table.csv", row.names = FALSE)

message("Subtype vs Subtype DESeq2 analysis complete!")
