# Load required libraries
library(DESeq2)
library(tidyverse)
library(ggrepel)

# Read count data
count_data <- read.csv("DESeq2Input/breastFilteredGenesCounts.csv", header = TRUE, row.names = 1, check.names = FALSE)

# Read metadata file
#metadata <- read.csv("DESeq2Input/breast_diagnosis_grouped_Metadata.csv", header=TRUE, row.names=1)
metadata <- read.csv("DESeq2Input/breastMetadata.csv", header=TRUE, row.names=1)
# Ensure column names in count data match row names in metadata
if (!all(colnames(count_data) %in% rownames(metadata))) {
  stop("Mismatch between count data columns and metadata row names!")
}

# Set up
normal_condition <- "normal"
unique_conditions <- unique(metadata$condition)
cancer_conditions <- unique_conditions[unique_conditions != normal_condition]
summary_list <- list()

table(metadata$condition)
print(sum(metadata$condition == normal_condition))

label_map <- c(
  infiltrating_duct_and_lobular_carcinoma = "Infiltrating Duct and Lobular Carcinoma",
  infiltrating_lobular_carcinoma = "Infiltraing Lobular Carcinoma",
  other_adenocarcinoma  = "Other Adenocarcinoma",
  infiltrating_ductal_carcinoma = "Infiltrating Ductal Carcinoma",
  other_carcinoma = "Other Carcinoma",
  not_reported = "Not Reported",
  cancer = "Tumor"
)

# Loop
for (cancer_type in cancer_conditions) {
  print(sum(metadata$condition == cancer_type))
  
  subset_meta <- metadata[metadata$condition %in% c(normal_condition, cancer_type), , drop = FALSE]
  subset_counts <- count_data[, rownames(subset_meta)]
  subset_meta$condition <- factor(subset_meta$condition, levels = c(normal_condition, cancer_type))
  
  if (sum(subset_meta$condition == cancer_type) < 6) {
    message("Skipping ", cancer_type, " — too few samples for comparison.")
    next
  }
  
  
  plot_label <- label_map[[cancer_type]]
  if (is.null(plot_label)) plot_label <- cancer_type
  message("Running comparison: ", cancer_type)
  # Create output folder
  out_dir <- paste0("results/DiagnosisGrouped_Breast/", cancer_type, "_vs_normal/")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  dds <- DESeqDataSetFromMatrix(countData = subset_counts,
                                colData = subset_meta,
                                design = ~ condition)
  dds <- DESeq(dds)
  res <- results(dds, alpha = 0.05)
  coef_name <- resultsNames(dds)[grepl(cancer_type, resultsNames(dds))]
  
  if (length(coef_name) == 0) {
    message("Skipping ", cancer_type, " — no matching coefficient found.")
    next
  }
  
  res_shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm")
  sig_genes <- subset(res_shrunk, padj < 0.05 & abs(log2FoldChange) > 1)
  message("Significant DEGs: ", nrow(sig_genes))
  
  # Summary stats
  up_genes <- sum(sig_genes$log2FoldChange > 1, na.rm = TRUE)
  down_genes <- sum(sig_genes$log2FoldChange < -1, na.rm = TRUE)
  top_gene <- rownames(sig_genes)[which.min(sig_genes$padj)[1]]
  
  summary_list[[cancer_type]] <- data.frame(
    comparison = paste0(cancer_type, "_vs_normal"),
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
  
  # MA Plot
  png(paste0(out_dir, "MA_plot.png"), width = 1000, height = 800, res = 150)
  plotMA(res, main = paste("MA Plot:", plot_label, "vs Normal"), ylim = c(-5, 5))
  abline(h = c(-1, 1), col = "blue", lty = 2)
  dev.off()
  
  # Save summary
  sink(paste0(out_dir, "summary.txt"))
  summary(res)
  sink()
  
  # Volcano Plot
  res_df <- as.data.frame(res_shrunk)
  res_df$gene <- rownames(res_df)
  
  # Categorize for color
  res_df$volcano_color <- "NS"
  res_df$volcano_color[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Up"
  res_df$volcano_color[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Down"
  
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
    ggtitle(paste("Volcano Plot:", plot_label, "vs Normal")) +
    xlab("Log2 Fold Change") +
    ylab("-log10 Adjusted p-value")
  
  ggsave(paste0(out_dir, "volcano_plot.png"),
         plot = volcano, width = 7.47, height = 6.22)
}


# Save summary CSV
summary_df <- bind_rows(summary_list)
write.csv(summary_df, file = "results/DEG_breast_new_summary_table.csv", row.names = FALSE)

message("DESeq2 analysis complete. Results saved.")
