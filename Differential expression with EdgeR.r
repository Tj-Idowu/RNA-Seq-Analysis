```{r}
## Load in the library
library(edgeR)
# read in count dataset
Samples <- as.matrix(read.csv("hervK_samples.csv"))
## Change Column names
colnames(Samples)<- c("1N","1T","2N","2T","3N","3T","4N",
"4T","5N","5T","6N","6T","7N","7T","8N","8T","9N","9T","10N","10T",
"11N","11T","12N","12T","13N","13T","14N","14T")
group<-c(1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2)

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

# Pricipal Component Analysis
pca <- prcomp(t(cpm_log), scale. = TRUE)
plot(pca$x[, 1], pca$x[, 2], pch = ".", xlab = "PC1", ylab = "PC2")
text(pca$x[, 1], pca$x[, 2], labels = colnames(cpm_log),col = group)
summary(pca)

# Comparison between normal and cancer cells
y <- DGEList(counts = data_clean, group = group)
y
y <- calcNormFactors(y)
y$samples

# DE between normal and tumour pair match
y <- estimateDisp(y)
sqrt(y$common.dispersion) # biological coefficient of variation
plotBCV(y)

# Fisher's Exact Test
Fisher <- exactTest(y)
results_edgeR <- topTags(Fisher, n = nrow(data_clean), sort.by = "none")
head(results_edgeR$table)

# Estimate False Discovery Rate at 10\%
sum(results_edgeR$table$FDR < .1)
plotSmear(Fisher, de.tags = rownames(results_edgeR)[results_edgeR$table$FDR < .1])
abline(h = c(-2, 2), col = "blue")

# Select data based on P-value
Sig_0.05<-subset(results_edgeR$table, results_edgeR$table$PValue<=0.05)
Bonf<-p.adjust(Sig_0.05$PValue,method="bonferroni")
Sig_0.05<-cbind(Sig_0.05,Bonf)
Sig_0.05 <- Sig_0.05[order(-Sig_0.05$Bonf),]

Sig_0.01<-subset(results_edgeR$table, results_edgeR$table$PValue<=0.01)
# Get adjusted P-value
Bonf<-p.adjust(Sig_0.01$PValue,method="bonferroni")
Sig_0.01<-cbind(Sig_0.01,Bonf)
Sig_0.01 <- Sig_0.01[order(Sig_0.01$Bonf),]
```
