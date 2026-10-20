#!/bin/bash
#SBATCH --job-name=NAME
#SBATCH -N 4
#SBATCH --ntasks-per-node 5
#SBATCH -t 2:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=MAIL.ME@email.com

# Path to Directories

# Path to Raw Reads
raw_dir="path/to/raw_reads"

# Path to Project
proj_dir="path/to/project/dir"

# Path to Trimmomatic
trimmomatic_dir="$project_dir/package/Trimmomatic-0.38"
trimmomatic_jar="$trimmomatic_dir/trimmomatic-0.38.jar"
adapters="$trimmomatic_dir/adapters/TruSeq3-PE.fa"

# Path to Reference Directories
ref_dir="$proj_dir/ensembl"
ref_fasta="$ref_dir/Homo_sapiens.GRCh37.cdna.fa"
salmon_index="$ref_dir/refcdna_index"

# Other Analysis Directories
trimmed_dir="$raw_dir/trimmed"
unpaired_dir="$raw_dir/unpaired"
trimlog_dir="$raw_dir/trimlog"
qc_raw_dir="$raw_dir/qc_raw"
qc_trim_dir="$raw_dir/qc_trim"
salmon_out="$proj_dir/salmon"

threads=8

mkdir -p "$trimmed_dir" "$unpaired_dir" "$trimlog_dir" "$qc_raw_dir" "$qc_trim_dir" "$salmon_out"

# Quality Check on the Raw Files with FastQC and report with multiQC

for f in "$raw_dir"/*.fastq.gz; 
do
        fastqc -o "$qc_raw_dir" "$f" 
done

multiqc --force -o "$qc_raw_dir" "$qc_raw_dir"

# Trim the Raw Reads with Trimmomatic

for f1 in "$raw_dir"/*_1.fastq.gz;
do
        f2=${f1/_1.fastq.gz/_2.fastq.gz}
        sample="$(basename "$f1" _1.fastq.gz)"
        r1="$(basename "$f1")"
        r2="$(basename "$f2")"
        
        java -jar "$trimmomatic_jar" PE -threads "$threads" -phred33 \
        -trimlog "$trimlog_dir/${sample}_trimlog.txt" \
	"$f1" "$f2" \
	"$trimmed_dir"/"$r1" "$unpaired_dir"/"$r1" \
	"$trimmed_dir"/"$r2" "$unpaired_dir"/"$r2" \
	ILLUMINACLIP:"$adapters":2:30:10 LEADING:3 TRAILING:3 \
	SLIDINGWINDOW:4:15 MINLEN:36
done

# Quality Check on the Trimmed Files with FastQC and report with multiQC

for f in "$trimmed_dir"/*.fastq.gz; 
do 
        fastqc -o "$qc_trim_dir" "$f" 
done

multiqc --force -o "$qc_trim_dir" "$qc_trim_dir"

# Salmon Index Build

salmon index -t "$ref_fasta" -i "$salmon_index" -p "$threads"

# Read Quantification with Salmon
for fn in "$trimmed_dir"/*_1.fastq.gz;
do
	fn2="${fn/_1.fastq.gz/_2.fastq.gz}"
	sample="$(basename "$fn" _1.fastq.gz)"
	
	salmon quant -i "$salomn_index" -l A --dumpEq --gcBias --seqBias \
         -1 ${fn} \
         -2 ${fn2} \
         -p "$threads" -o "$salmon_out/$sample"
done
