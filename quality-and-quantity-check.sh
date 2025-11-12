# Shell script to trim reads
#!/bin/bash
# Trimmomatic
mkdir trimlog
mkdir trimmed
mkdir unpaired
for f1 in *_1.fastq.gz
do
        f2=${f1/_1.fastq.gz/_2.fastq.gz}
        java -jar ../packages/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 8 -phred33 \
        -trimlog trimlog/${f1/_1.fastq.gz}_trimlog.txt \
"${f1}" "${f2}" trimmed/"$f1" unpaired/"$f1" trimmed/"$f2" \
unpaired/"$f2" \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 \
SLIDINGWINDOW:4:15 MINLEN:36
done

# Quality Check on the Trimmed Files
#!/bin/bash
# Quality check

for f in *.fastq.gz; 
do 
        fastqc $f; 
done

multiqc .;

# Read Quantification with Salmon
#!/bin/bash
# indexing the reference
salmon index -t ensembl/Homo_sapiens.GRCh37.cdna.all.fa \
-i ensembl/refcDNA_index

# Analyse Samples
for fn in *_1.fastq.gz;
do
fn2=${fn/_1.fastq.gz/_2.fastq.gz}
salmon quant -i ../ensembl/refcDNA_index -l A --dumpEq --gcBias --seqBias \
         -1 ${fn} \
         -2 ${fn2} \
         -p 8 -o ../salmon/${fn/_1.fastq.gz}
done
