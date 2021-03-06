#!/bin/bash
#SBATCH -J FilterMetagenomicReads
#SBATCH -t 100:00:00
#SBATCH -n 24
#SBATCH -N 1
#SBATCH --mem=50g
#SBATCH -o FilterMetagenomicReads.out
#SBATCH -e FilterMetagenomicReads.err

module load samtools
module load seqtk
module load bowtie2/2.3.5.1
module load bbmap

for file in `cat filelist.txt`
do
java -jar /gpfs/data/rbeinart/bin/trinityrnaseq-v2.11.0/trinity-plugins/Trimmomatic/trimmomatic-0.36.jar PE -threads 24 -phred33 /gpfs/data/rbeinart/Raw_reads/Snail_metagenomics_UCD/${file}_*R1_001.fastq.gz /gpfs/data/rbeinart/Raw_reads/Snail_metagenomics_UCD/${file}_*R2_001.fastq.gz ${file}_R1_paired.fq ${file}_R1_unpaired.fq ${file}_R2_paired.fq ${file}_R2_unpaired.fq ILLUMINACLIP:/gpfs/data/rbeinart/cbreusing/Adapters/Tn5.fasta:2:30:10 SLIDINGWINDOW:4:20 LEADING:5 TRAILING:5 MINLEN:75
bowtie2 -p 24 -x /gpfs/data/rbeinart/Databases/contaminants -1 ${file}_R1_paired.fq -2 ${file}_R2_paired.fq | samtools view -bS -h -@ 24 - | samtools sort -@ 24 - > ${file}.bowtie2.cont.sorted.bam
samtools view -@ 24 -f12 ${file}.bowtie2.cont.sorted.bam > ${file}.cont.unmapped.sam
cut -f1 ${file}.cont.unmapped.sam | sort | uniq > ${file}.cont.unmapped_ids.lst
seqtk subseq ${file}_R1_paired.fq ${file}.cont.unmapped_ids.lst > ${file}_R1_clean.fastq
seqtk subseq ${file}_R2_paired.fq ${file}.cont.unmapped_ids.lst > ${file}_R2_clean.fastq
done

# Separate host reads from symbiont reads for downstream analyses
bbsplit.sh -Xmx50g ref=/gpfs/data/rbeinart/cbreusing/Genomes/Ifr_SOX.fasta,/gpfs/data/rbeinart/cbreusing/Genomes/Gamma1.fasta,/gpfs/data/rbeinart/cbreusing/Genomes/GammaLau.fna,/gpfs/data/rbeinart/cbreusing/Genomes/Epsilon.fasta,/gpfs/data/rbeinart/cbreusing/Genomes/Ifr_MOX.fasta build=1 path=/gpfs/data/rbeinart/cbreusing/Genomes/ in=${file}_R1_clean.fastq in2=${file}_R2_clean.fastq ambiguous=best ambiguous2=toss outu=${file}_unmapped_#.fq

