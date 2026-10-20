## Load in the library
library(tidyverse)
library(tximport)
library(edgeR)

# Path to input files
input_dir <- "/absolute/path/to/salmon/output"
ref_fasta <- "/absolute/path/to/ensembl/Homo_sapiens.GRCh37.cdna.all.fa"
sample_sheet <- "sample_sheet.csv"
counts_file <- "hervk_samples.csv"

# Input Sample Sheet
# "sample" column is the output folder name for that sample, patient is numbered 1-14, tissues and N "normal" or T "tumour"
samples <- read_csv(sample_sheet, col_types = cols(sample = col_character(), patient = col_integer(), tissue = col_character())) %>%
	mutate(tissue = factor(tissue, levels = c("N", "T")), label = paste0(patient, tissue)) %>%
	# order the tissues as 1N, 1T, 2N, 2T, ...14N, 14T
	arrange(patient, tissue)

# Build the count matrix
tx2gene <- tibble(header = str_subset(read_lines(ref_fasta), "^>\\S+")) %>%
	transmutate(tx_id = str_remove(str_sxact(header, "^>\\S+"), "^>"),
		gene_id = str_extract(header, "(?<+\\sgene:)\\S+"))

files <- file.path(input_dir, samples$sample, "quant.sf")
names(files) <- samples$label

# Get counts that have been corrected for transcript length
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")

# Write out to .csv file
txi$counts %>%
	as.data.frame() %>%
	rownames_to_column("gene_id") %>%
	write_csv(count_file)

# Read in the count matrix
Samples <- read_csv(counts_file, show_col_types = FALSE) %>%
	column_to_rownames("gene_id") %>%
	as.matrix()

Samples <- Samples[, samples$label]
tissue <- samples$tissue
patient <- factor(samples$patient)

# Filter low expression genes
# Get counts per million
cpm_log <- cpm(Samples, log = TRUE)
median_log2_cpm <- apply(cpm_log, 1, median)
# Visual so you know where to set the cut off
hist(median_log2_cpm)
# Manual input of where cut off line should go
expr_cutoff <- -1
# Add cut off line
abline(v = expr_cutoff, col = "red", lwd = 3)

# Indicate what data will be plotted
# where median log2 of the CPM is greater than the cut off
sum(median_log2_cpm > expr_cutoff)

# Remove the cut off data
data_clean <- Samples[median_log2_cpm > expr_cutoff, ]

# Get CPM of new figures sans cut off figures
cpm_log <- cpm(data_clean, log = TRUE)
hist(cpm_log)
heatmap(cor(cpm_log))

# Principal Component Analysis, black for normal and red for tumour cells
pca <- prcomp(t(cpm_log), scale. = TRUE)
plot(pca$x[, 1], pca$x[, 2], type = "n", xlab = "PC1", ylab = "PC2")
text(pca$x[, 1], pca$x[, 2], labels = colnames(cpm_log),col = as.integer(tissue))
summary(pca)

# Comparison between normal and cancer cells
y <- DGEList(counts = data_clean, group = tissue)
y
y <- calcNormFactors(y)
y$samples

# DE between normal and tumour pair match
design <- model.matrix(~ patient + tissue)
y <- estimateDisp(y, design)
sqrt(y$common.dispersion) # biological coefficient of variation
plotBCV(y)

# DE with EdgeR quasi-likelihood F-test (recommended by edgeR)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = "tissueT")
results_edgeR <- topTags(qlf, n = Inf, sort.by = "none")$table
head(results_edgeR)

# Estimate False Discovery Rate at 10\%
sum(results_edgeR$FDR < .1)
plotMD(qlf, status = decideTests(qlf, p.value = 0.1))
abline(h = c(-2, 2), col = "blue")

# Select data based on P-value
results_edgeR <- results_edgeR %>%
	rownames_to_column("gene_id") %>%
	mutate(Bonf = p.adjust(PValue, method = "bonferroni"))

Sig_0.05 <- results_edgeR %>% filter(PValue <- 0.05) %>% arrange(Bonf)
Sig_0.01 <- results_edgeR %>% filter(PValue <- 0.01) %>% arrange(Bonf)

# Write out to .csv files
write_csv(results_edgeR, "all_edgeR_results.csv")
write_csv(Sig_0.05, "p05_edgeR_results.csv")
write_csv(Sig_0.01, "p01_edgeR_results.csv")
