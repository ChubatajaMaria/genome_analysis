#!/bin/bash
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 02:00:00
#SBATCH -J deseq2_analysis
#SBATCH --mail-type=ALL
#SBATCH --mail-user mariia.chubataia.8934@student.uu.se
#SBATCH --output=%x.%j.out

# Load required modules
module load bioinfo-tools
module load R_packages/4.1.1
module load R/4.1.1

# Set working directory
cd /home/mach8934/genome_analysis/data/09_read_counting/

# Create R script within SLURM job
cat > deseq2_analysis.R << 'EOF'
# Set library path to home directory
.libPaths(c(paste0(Sys.getenv("HOME"), "/R/library"), .libPaths()))

# Check and install required packages if needed
packages_needed <- c("DESeq2", "apeglm", "pheatmap")

# Function to check and install packages
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cran.rstudio.com/")
    }
    library(BiocManager)
    if (pkg %in% c("DESeq2", "apeglm")) {
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    } else {
      install.packages(pkg, repos = "https://cran.rstudio.com/")
    }
  }
}

# Install required packages
sapply(packages_needed, install_if_missing)

# Load essential libraries
library("DESeq2")
library("pheatmap")
library("ggplot2")

# Import gene count data
gene_counts <- read.delim("counts.txt", comment.char = "#", stringsAsFactors = FALSE)

# Data structure: Geneid, Chr, Start, End, Strand, Length, Control_1, Control_2, Control_3, Treatment_1, Treatment_2, Treatment_3
# Keep only Geneid and count columns (remove genomic coordinates)
count_matrix <- gene_counts[, c(1, 7:ncol(gene_counts))]

# Set gene IDs as row names and remove the Geneid column
rownames(count_matrix) <- count_matrix$Geneid
count_matrix <- count_matrix[-1]

# Create sample metadata
sample_info <- data.frame(
  row.names = colnames(count_matrix),
  group = factor(rep(c("control", "treatment"), each = 3)),
  batch = factor(c(1, 2, 3, 1, 2, 3))  # Added batch information
)

# Create DESeq2 dataset object
dds_object <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_info,
  design = ~ group
)

# Display dataset information
print(dds_object)

# Quality control: remove genes with low counts
# Keep genes with at least 5 counts across samples
expressed_genes <- rowSums(counts(dds_object)) >= 5
dds_object <- dds_object[expressed_genes, ]

# Set reference level for comparison
dds_object$group <- relevel(dds_object$group, ref = "control")

# Perform differential expression analysis
dds_object <- DESeq(dds_object)
differential_results <- results(dds_object)

# Show results summary
print(differential_results)
print(resultsNames(dds_object))

# Apply log fold change shrinkage
shrunken_results <- lfcShrink(dds_object, coef = "group_treatment_vs_control", type = "apeglm")

# Display summaries
summary(differential_results)
summary(shrunken_results)

# Generate MA plots
png("ma_plots.png", width = 1200, height = 600)
par(mfrow = c(1, 2))
plotMA(differential_results, main = "Unshrunken MA Plot")
plotMA(shrunken_results, ylim = c(-5, 5), main = "Shrunken MA Plot")
par(mfrow = c(1, 1))
dev.off()

# Plot top differentially expressed gene
top_gene <- rownames(differential_results)[which.min(differential_results$padj)]
png("top_gene_counts.png", width = 800, height = 600)
plotCounts(dds_object, gene = top_gene, intgroup = "group", 
          main = paste("Expression of", top_gene))
dev.off()

# Export results to CSV
results_df <- as.data.frame(shrunken_results)
results_df$gene_id <- rownames(results_df)
write.csv(results_df, "differential_expression_results.csv", row.names = FALSE)

# Create expression categories based on fold change
results_df$expression_category <- cut(
  results_df$log2FoldChange,
  breaks = c(-Inf, -1.5, -0.5, 0.5, 1.5, Inf),
  labels = c("Strongly Down", "Moderately Down", "No Change", 
             "Moderately Up", "Strongly Up")
)

# Enhanced scatter plot
png("expression_scatter.png", width = 1000, height = 800)
expression_plot <- ggplot(results_df, aes(x = baseMean, y = log2FoldChange, 
                                         color = expression_category)) +
  geom_point(size = 0.8, alpha = 0.7) +
  scale_x_log10() +
  scale_color_brewer(type = "div", palette = "RdYlBu", 
                     name = "Expression Change") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(title = "Gene Expression vs. Mean Expression Level",
       x = "Average Normalized Expression",
       y = "Log2 Fold Change (Shrunken)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

print(expression_plot)
dev.off()

# Simple volcano plot
results_df$significance <- ifelse(results_df$padj < 0.05 & abs(results_df$log2FoldChange) > 1, 
                                  ifelse(results_df$log2FoldChange > 1, "Upregulated", "Downregulated"), 
                                  "Not significant")

png("volcano_plot.png", width = 1000, height = 800)
vol_plot <- ggplot(results_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not significant" = "gray50")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(title = "Volcano Plot: Treatment vs Control",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

print(vol_plot)
dev.off()

# Hierarchical clustering heatmap
# Apply variance stabilizing transformation
vst_data <- vst(dds_object, blind = FALSE)

# Get top 40 significant genes
significant_genes <- shrunken_results[!is.na(shrunken_results$padj) & 
                                     shrunken_results$padj < 0.05, ]
top_genes_list <- head(rownames(significant_genes[order(significant_genes$padj), ]), 40)

# Extract transformed counts for heatmap
heatmap_matrix <- assay(vst_data)[top_genes_list, ]
heatmap_scaled <- t(scale(t(heatmap_matrix)))

# Create annotation for heatmap columns
col_annotation <- as.data.frame(colData(dds_object)[, "group", drop = FALSE])
colnames(col_annotation) <- "Treatment"

# Generate heatmap
png("heatmap.png", width = 1000, height = 1200)
pheatmap(
  heatmap_scaled,
  annotation_col = col_annotation,
  show_rownames = TRUE,
  fontsize_row = 8,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "none",
  color = colorRampPalette(c("#2166AC", "white", "#D6604D"))(100),
  main = "Top 40 Differentially Expressed Genes"
)
dev.off()

# Principal Component Analysis
pca_results <- plotPCA(vst_data, intgroup = "group", returnData = TRUE)
percent_variance <- round(100 * attr(pca_results, "percentVar"))

png("pca_plot.png", width = 1000, height = 800)
pca_plot <- ggplot(pca_results, aes(PC1, PC2, color = group, shape = group)) +
  geom_point(size = 5) +
  labs(title = "Principal Component Analysis",
       x = paste0("PC1: ", percent_variance[1], "% of variance"),
       y = paste0("PC2: ", percent_variance[2], "% of variance")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.title = element_text(size = 12))

print(pca_plot)
dev.off()

# Summary statistics
cat("\n=== Analysis Summary ===\n")
cat("Total genes analyzed:", nrow(dds_object), "\n")
cat("Significantly upregulated genes (padj < 0.05, log2FC > 1):", 
    sum(results_df$padj < 0.05 & results_df$log2FoldChange > 1, na.rm = TRUE), "\n")
cat("Significantly downregulated genes (padj < 0.05, log2FC < -1):", 
    sum(results_df$padj < 0.05 & results_df$log2FoldChange < -1, na.rm = TRUE), "\n")
cat("Analysis completed successfully!\n")

# Save session info for reproducibility
sessionInfo()
EOF

# Run the R script
echo "Starting DESeq2 analysis..."
Rscript deseq2_analysis.R

echo "Analysis completed at $(date)"
echo "Check the following output files:"
echo "- differential_expression_results.csv"
echo "- ma_plots.png"
echo "- top_gene_counts.png"
echo "- expression_scatter.png"
echo "- volcano_plot.png"
echo "- heatmap.png"
echo "- pca_plot.png"
