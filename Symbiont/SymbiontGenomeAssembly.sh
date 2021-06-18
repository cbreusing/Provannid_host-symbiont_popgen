#!/bin/bash
#SBATCH -J SymbiontGenomeAssembly
#SBATCH -t 100:00:00
#SBATCH -n 24
#SBATCH -N 1
#SBATCH --mem=100g
#SBATCH -o SymbiontGenomeAssembly.out
#SBATCH -e SymbiontGenomeAssembly.err

module load python/3.6.6
module load gcc/5.4
module load samtools
module load seqtk
module load bowtie2/2.3.5.1
module load bbmap

file=("Alv337")

# Filtering of Illumina reads for symbiont genome assembly
java -jar /gpfs/data/rbeinart/bin/trinityrnaseq-v2.11.0/trinity-plugins/Trimmomatic/trimmomatic-0.36.jar PE -threads 24 -phred33 ${file}_R1.fastq.gz ${file}_R2.fastq.gz ${file}_R1_paired.fq ${file}_R1_unpaired.fq ${file}_R2_paired.fq ${file}_R2_unpaired.fq ILLUMINACLIP:/gpfs/data/rbeinart/cbreusing/Adapters/Illumina.fa:2:30:10 SLIDINGWINDOW:4:20 LEADING:5 TRAILING:5 MINLEN:50
bowtie2 -p 24 -x /gpfs/data/rbeinart/Databases/contaminants -1 ${file}_R1_paired.fq -2 ${file}_R2_paired.fq | samtools view -bS -h -@ 24 - | samtools sort -@ 24 - > ${file}.bowtie2.cont.sorted.bam
samtools view -@ 24 -f12 ${file}.bowtie2.cont.sorted.bam > ${file}.cont.unmapped.sam
cut -f1 ${file}.cont.unmapped.sam | sort | uniq > ${file}.cont.unmapped_ids.lst
seqtk subseq ${file}_R1_paired.fq ${file}.cont.unmapped_ids.lst > ${file}_R1_clean.fastq
seqtk subseq ${file}_R2_paired.fq ${file}.cont.unmapped_ids.lst > ${file}_R2_clean.fastq
# Separate reads based on draft genomes assembled in Beinart et al. (2019)
bbsplit.sh -Xmx50g ref=/gpfs/data/rbeinart/cbreusing/Genomes/Gamma1.fna,/gpfs/data/rbeinart/cbreusing/Genomes/GammaLau.fna,/gpfs/data/rbeinart/cbreusing/Genomes/Epsilon.fna,/gpfs/data/rbeinart/cbreusing/Genomes/Ifr_SOX.fna build=1 path=/gpfs/data/rbeinart/cbreusing/Genomes/ in=${file}_R1_clean.fastq in2=${file}_R2_clean.fastq ambiguous=best ambiguous2=toss basename=${file}_%_#.fq outu=${file}_unmapped_#.fq

# Filtering of Nanopore reads for symbiont genome assembly
/gpfs/data/rbeinart/Software/Porechop/porechop-runner.py -i ${file}_nanopore.fastq -o ${file}_nanopore_trim.fastq --threads 24 -v 1
minimap2 -t 24 -x map-ont -a Epsilon.fna ${file}_nanopore_trim.fastq -t 24 | samtools view -Sb -h -@ 24 - | samtools sort -@ 24 - > ${file}.nanopore.sorted.bam
samtools view -@ 24 -F4 ${file}.nanopore.sorted.bam > ${file}.nanopore.mapped.sam
cut -f1 ${file}.nanopore.mapped.sam | sort | uniq > ${file}.mapped_ont_ids.lst
seqtk subseq ${file}_nanopore_trim.fastq ${file}.mapped_ont_ids.lst > ${file}_symbiont.nanopore.fastq

# Genome assembly
metaspades.py -1 ${file}_Epsilon_1.fq -2 ${file}_Epsilon_2.fq --nanopore ${file}_symbiont.nanopore.fastq -k 21,31,41,51,61,71,81,91 -t 24 -o metaSpades

# Detailed metagenome binning instructions are provided at: https://github.com/kbseah/genome-bin-tools/wiki
