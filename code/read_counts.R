# Load libraries
library(ggplot2)
library(reshape2)

# Set working directory
setwd('/home/mach8934/genome_analysis/data/09_read_counting/')

# Check if file exists
if (!file.exists("counts.txt")) {
    cat("ERROR: counts.txt file not found in current directory!\n")
    cat("Current working directory:", getwd(), "\n")
    cat("Files in directory:\n")
    print(list.files())
    stop("File not found")
}

cat("Loading data from:", file.path(getwd(), "counts.txt"), "\n")

# Load the data
counts <- read.table("counts.txt", header = TRUE, sep = "\t", comment.char = "")

# Extract the count columns
# Columns 7 onwards contain the count data
count_data <- counts[, 7:ncol(counts)]

# Compute total count per gene
counts$TotalCount <- rowSums(count_data)

# Remove zeros to avoid log10(0)
counts_nonzero <- counts[counts$TotalCount > 0, ]

# Apply log10 to total counts
counts_nonzero$Log10Count <- log10(counts_nonzero$TotalCount)

# Create output directory for plots
dir.create("plots", showWarnings = FALSE)

# Plot and save as PDF
pdf("plots/gene_counts_histogram.pdf", width = 8, height = 6)
hist(counts_nonzero$Log10Count,
     breaks = 50,
     col = "steelblue",
     border = "black",
     main = "Histogram of Gene Counts (log10)",
     xlab = "log10 Total Counts per Gene",
     ylab = "Frequency")
dev.off()

# Also save as PNG for easy viewing
png("plots/gene_counts_histogram.png", width = 800, height = 600)
hist(counts_nonzero$Log10Count,
     breaks = 50,
     col = "steelblue",
     border = "black",
     main = "Histogram of Gene Counts (log10)",
     xlab = "log10 Total Counts per Gene",
     ylab = "Frequency")
dev.off()

# Print summary statistics
cat("Summary of gene count analysis:\n")
cat("Total number of genes:", nrow(counts), "\n")
cat("Genes with non-zero counts:", nrow(counts_nonzero), "\n")
cat("Percentage of genes with counts:", round(100 * nrow(counts_nonzero) / nrow(counts), 2), "%\n")
cat("\nSummary statistics of log10 counts:\n")
print(summary(counts_nonzero$Log10Count))

# Save the processed data
write.table(counts_nonzero, "processed_counts.txt", sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nAnalysis completed successfully!\n")
cat("Plots saved in: plots/gene_counts_histogram.pdf and plots/gene_counts_histogram.png\n")
cat("Processed data saved as: processed_counts.txt\n")
