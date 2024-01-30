#!/bin/bash

set -eux

# Script to call germline variants in a human WGS paired end reads 2 X 100bp
# Following GATK4 best practices workflow - https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-
if false;
then
#To download the data 
wget -P /mnt/a/ani_linux/pipeline/DNASeq_pipeline/reads https://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_1.filt.fastq.gz
wget -P /mnt/a/ani_linux/pipeline/DNASeq_pipeline/reads https://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_2.filt.fastq.gz

echo "Run prep files..."

#download reference file
wget -P /mnt/a/ani_linux/pipeline/DNASeq_pipeline/data https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip /mnt/a/ani_linux/pipeline/DNASeq_pipeline/data/hg38.fa.gz


# index ref - .fai file before running haplotype caller
samtools faidx /mnt/a/ani_linux/pipeline/DNASeq_pipeline/data/hg38.fa


# ref dict - .dict file before running haplotype caller
gatk CreateSequenceDictionary R=/mnt/a/ani_linux/pipeline/DNASeq_pipeline/data/hg38.fa o=/mnt/a/ani_linux/pipeline/DNASeq_pipeline/data/hg38.dict 

# download known sites files for BQSR from GATK resource bundle
wget -P /mnt/a/ani_linux/pipeline/DNASeq_pipeline/data/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget -P /mnt/a/ani_linux/pipeline/DNASeq_pipeline/data/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx
fi

#directories
ref="/mnt/a/ani_linux/pipeline/DNASeq_pipeline/data/hg38.fa"
known_sites="/mnt/a/ani_linux/pipeline/DNASeq_pipeline/data/Homo_sapiens_assembly38.dbsnp138.vcf"
aligned_reads="/mnt/a/ani_linux/pipeline/DNASeq_pipeline/aligned_reads"
reads="/mnt/a/ani_linux/pipeline/DNASeq_pipeline/reads"
results="/mnt/a/ani_linux/pipeline/DNASeq_pipeline/results"
data="/mnt/a/ani_linux/pipeline/DNASeq_pipeline/data"

fastqc ${reads}/SRR062634_1.filt.fastq.gz -o= ${reads}/
fastqc ${reads}/SRR062634_2.filt.fastq.gz -o= ${reads}/

#map using bwa-mem

#bwa index reference

bwa index ${ref}

bwa alignment

bwa mem -t 4 -R "@RG\tID:SRR062634\tPL:ILLUMINA\tSM:SRR062634" ${ref} ${reads}/small_SRR062634_1.filt.fastq.gz ${reads}/small_SRR062634_2.filt.fastq.gz > ${aligned_reads}/small_SRR062634.paired.sam

# STEP 3: Mark Duplicates and Sort - GATK4
gatk SortSam -I ${aligned_reads}/small_SRR062634.paired.sam -O ${aligned_reads}/small_SRR062634.sorted.sam -SO coordinate

#echo "STEP 3: Mark Duplicates and Sort - GATK4"
gatk MarkDuplicates -I ${aligned_reads}/small_SRR062634.sorted.sam -O ${aligned_reads}/small_SRR062634_sorted_dedup_reads.bam -M dup_entriers -AS true

sudo chmod +rwx ${aligned_reads}/SRR062634_sorted_dedup_reads.bam

#STEP 4: Base quality recalibration

echo "STEP 4: Base quality recalibration"

# 1. build the model
gatk BaseRecalibrator -I ${aligned_reads}/small_SRR062634_sorted_dedup_reads.bam -R ${ref} --known-sites ${known_sites} -O ${data}/small_recal_data.table


# 2. Apply the model to adjust the base quality scores

gatk ApplyBQSR -I ${aligned_reads}/small_SRR062634_sorted_dedup_reads.bam -R ${ref} --bqsr-recal-file ${data}/recal_data.table -O ${aligned_reads}/small_SRR062634_sorted_dedup_bqsr_reads.bam 

#-----------------------------------------------
#STEP 5: Collect Alignment & Insert Size Metrics
#-----------------------------------------------


echo "STEP 5: Collect Alignment & Insert Size Metrics"

gatk CollectAlignmentSummaryMetrics R=${ref} I=${aligned_reads}/small_SRR062634_sorted_dedup_bqsr_reads.bam O=${aligned_reads}/small_alignment_metrics.txt
gatk CollectInsertSizeMetrics INPUT=${aligned_reads}/small_SRR062634_sorted_dedup_bqsr_reads.bam OUTPUT=${aligned_reads}/small_insert_size_metrics.txt HISTOGRAM_FILE=${aligned_reads}/small_insert_size_histogram.pdf



# ----------------------------------------------
# STEP 6: Call Variants - gatk haplotype caller
# ----------------------------------------------

echo "STEP 6: Call Variants - gatk haplotype caller"

gatk HaplotypeCaller -R ${ref} -I ${aligned_reads}/small_SRR062634_sorted_dedup_bqsr_reads.bam -O ${results}/small_raw_variants.vcf



# extract SNPs & INDELS

gatk SelectVariants -R ${ref} -V ${results}/small_raw_variants.vcf --select-type SNP -O ${results}/small_raw_snps.vcf
gatk SelectVariants -R ${ref} -V ${results}/small_raw_variants.vcf --select-type INDEL -O ${results}/small_raw_indels.vcf
