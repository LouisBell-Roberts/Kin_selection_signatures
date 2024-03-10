#!/bin/bash
#SBATCH --job-name=Camp_flo_SNP
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --clusters=htc
#SBATCH --partition=medium
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=16G
#SBATCH --mail-user=louis.bell-roberts@sjc.ox.ac.uk
module purge
module load Anaconda3//2020.11
source activate $DATA/myenv

#Louis Bell-Roberts
#SNP calling pipeline - place script within working directory containing the WGS reads, reference sequence and Truseq file. Set job name.
##Truseq file - specify whether single or paired reads. This will also alter the SNP calling script.
###Must specify ploidy

#Create a file with all sample names derived from the prefix of the fastq files in the working directory
##Names used for parallel computation
ls *.fastq.gz | cut -d '_' -f 1 | sort | uniq > names.txt

##### Trimmomatic #############################################
##Run Trimmomatic in parallel for each sample using parallel command. Parallel command will iterate over the lines of the names.txt file at the {} placeholder
##Trim adapters and low-quality bases from sequence data
cat names.txt | parallel -t "java -jar /data/zool-cooperation/newc4920/Trimmomatic-0.39/trimmomatic-0.39.jar \
PE -phred33 \
{}_1.fastq.gz {}_2.fastq.gz \
{}_1_out_paired.fq.gz {}_1_out_unpaired.fq.gz \
{}_2_out_paired.fq.gz {}_2_out_unpaired.fq.gz \
ILLUMINACLIP:TruSeq2-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 CROP:150"

###########BWA###################################
##Align paired-end reads to the reference sequence using parallel computation
##Use wildcard to identify the reference sequence and index it
refseq=*.fna
bwa index $refseq
cat names.txt | parallel -t "bwa aln $refseq {}_1_out_paired.fq.gz > {}_1.sai"
cat names.txt | parallel -t "bwa aln $refseq {}_2_out_paired.fq.gz > {}_2.sai"
cat names.txt | parallel -t "bwa sampe $refseq {}_1.sai {}_2.sai {}_1_out_paired.fq.gz {}_2_out_paired.fq.gz > {}_1_aln.sam"

###########SAMTOOLS########################
##Indexes the reference genome
##Convert SAM to BAM format, sort BAM files, and create index files for sorted BAM files
samtools faidx $refseq
cat names.txt | parallel -t "samtools view -S -b {}_1_aln.sam > {}_1_aln.bam"
cat names.txt | parallel -t "samtools sort -o sorted{}_1.bam {}_1_aln.bam"
cat names.txt | parallel -t "samtools index sorted{}_1.bam"

###########BCFTOOLS########################
##Call variants from the aligned sequence data
##Create a pileup file using the sample names stored in names.txt. Length determined by the number of names
bcftools mpileup -f $refseq $(cat names.txt | xargs -I{} echo "sorted{}_1.bam") > my-raw.bcf

##SNP calling
bcftools call -vc -Ov my-raw.bcf > my-raw.vcf

##Filter SNPs
###Quality score above 20
bcftools filter -i 'TYPE="snp" && QUAL>20' my-raw.vcf > step1.vcf
vcfutils.pl varFilter step1.vcf > OUTPUT.vcf
