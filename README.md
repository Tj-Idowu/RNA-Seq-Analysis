# HERV-K RNA-Seq Differential Expression Pipeline

The pipeline uses two paired-end RNA sequencing reads from matched normal and tumour prostate samples. It does quality checks, trims and quantifies transcript expression using Salmon. It then tests the differential expression of the normal versus tumour samples with EdgeR.

It is a two-step process.

1. **`quality_check_and_analysis.sh`**, a Bash script run on an HPC cluster, does QC with FastQC, adapter and quality trimming with Trimmomatic, and transcript quantification with Salmon.
2. **`Differential expression with EdgeR.r`**, an R script, does filtering, normalisation, data exploration and differential expression tests on the counts produced in the first step.

---

## Requiremnets

### Step 1
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [MultiQC](https://multiqc.info/)
- Java for TRimmomatic use: `$ module load java`
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic). The script expects Trimmomatic 0.38; if a different version is used, then change it in the script.
- The `TruSeq3-PE.fa` adapter file is in the working directory in the script. If this has changed, then edit the script.
- [Salmon](https://combine-lab.github.io/salmon/)
- Ensembl human cDNA reference (GRCh37)

  ### Step 2
  R with the Bioconductor packages;
  - [edgeR](https://bioconductor.org/packages/edgeR/)
  - [tximport](https://bioconductor.org/packages/tximport/)

  and `tidyverse` 
  installed

```r
install.packages("tidyverse")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("edgeR", "tximport"))
```

---

## Input data

- The bash script takes gzipped paired-end fastq files
- The R script takes 14 matched **normal (N)** and **tumour (T)** samples

---

## How it works

### Step 1: Quality Control and Read Quantification

The bash script runs on a SLURM HPC cluster (make sure to edit the conditions before running or change if using a different environment).

#### 1. Quality Check on Raw Reads
First, FastQC is run on each raw `*.fastq.gz` file, and a summary report is made with MultiQC

#### Trimming with Trimmomatic
Trimmomatic is run on each pair with:
| Step | Setting | Purpose |
|---|---|---|
| `ILLUMINACLIP` | `TruSeq3-PE.fa:2:30:10` | Remove Illumina TruSeq adapters |
| `LEADING` | `3` | Trim low-quality bases from the start |
| `TRAILING` | `3` | Trim low-quality bases from the end |
| `SLIDINGWINDOW` | `4:15` | Cut when average quality in a 4-base window drops below 15 |
| `MINLEN` | `36` | Drop reads shorter than 36 bases |

Outputs:
Outputs:
- `trimmed/`: surviving paired reads
- `unpaired/`: reads whose mate was discarded
- `trimlog/`: per-sample trimming logs

#### 3. Quality check on trimmed reads
Repeat FastQC and MultiQC on trimmed reads. Go through the report to check read quality

#### 4. Quantification with Salmon
- Build an index file for Salmon from the GRCh37 Ensembl reference
- Quantify each trimmed read against the index with library type detection, GC and sequence bias correction, and equivalence class output.

### Step 2: Differential Expression
Salmon creates a `quant.sf` file for each sample. For the differential expression analysis, all the sample outputs need to be in a single file. Save the single concatenated count matrix into a .csv file so that it does not need to be made again for subsequent reanalysis.

#### 1. Load and relabel the data
Read the count matrix, renaming the columns "1N" to "14N" and "1T" to "14T" where N = normal and T = tumour.

#### 2. Filter low-expression features
- Calculate the log2 counts per million (CPM) and the median log2 CPM for each feature.
- PLot a histogram of the medians to choose a cutoff; the default is `expr_cutoff <- -1`, which will be shown by a red line.
- Keep only features with a median log2 CPM above the cutoff.

#### 3. Exploratory Analysis
- Create a histogram of the log2 CPM values.
- Create a heatmap of the sample-to-sample correlations.
- Do principal component analysis (PCA), with samples labelled and coloured by group.

#### 4. Normalisation and Dispersion
- Build a `DGEList` and apply trimmed mean of M-values (TMM) normalisation with `calcNormFactors`.
- Estimate dispersion with `estimateDisp` to print the biological coefficient of variation and plot this with `plotBCV`.

#### 5. Differential Expression
- Test the normal versus tumour samples with edgeR's `exactTest`.
- Report how many features pass (false discovery rate) FDR <0.1 and show them on a smear plot with blue lines at log2 fold change ±2.

#### 6. Significant Results Tables
- Set the table in descending order of the Bonferroni-adjusted p-value column and isolate rows with p ≤ 0.05 into a new table.
- Set the table in descending order of the Bonferroni-adjusted p-value column and isolate rows with p ≤ 0.01 into a new table.
- Save these tables in .csv files
