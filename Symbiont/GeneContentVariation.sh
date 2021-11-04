#!/bin/bash
#SBATCH -J GeneContentVariation
#SBATCH -t 10:00:00
#SBATCH -n 24
#SBATCH -N 1
#SBATCH --mem=100g
#SBATCH -o GeneContentVariation.out
#SBATCH -e GeneContentVariation.err

module load samtools
module load bowtie2/2.3.5.1
module load python/3.7.4
module load bbmap
module load seqtk
module load perl
module load R

for file in `cat subset.list`
do
# Extract reads from downsampled bam file for analysis
samtools view -h -@ 24 -o ${file}.Epsilon.subsampled.sam ${file}.Epsilon.subsampled.bam
cut -f1 ${file}.Epsilon.subsampled.sam | sort | uniq > ${file}_ids.lst
seqtk subseq ${file}_R1_clean.fastq ${file}_ids.lst > ${file}_R1.fastq
seqtk subseq ${file}_R2_clean.fastq ${file}_ids.lst > ${file}_R2.fastq
reformat.sh in=${file}_R1.fastq in2=${file}_R2.fastq out=${file}.fastq
# Map interleaved reads to (pan)genome reference for presence/absence analysis (PanPhlan)
panphlan_map.py -i ${file}.fastq -o map_results/${file}_mapping.csv --nproc 4 --indexes panphlan/Epsilon -p panphlan/Epsilon_pangenome.tsv
# Map reads against protein-coding gene set for gene abundance analysis (Salmon)
bbmap.sh ref=Epsilon_pangenome.fasta nodisk in=${file}_R1.fastq in2=${file}_R2.fastq out=${file}.sam sam=1.4 local=t requirecorrectstrand=t samestrandpairs=f xstag=fs
salmon quant -t Epsilon_pangenome.fasta -l IU --meta --rangeFactorizationBins 4 --numBootstraps 1000 --seqBias --gcBias -s -u -a ${file}.sam -o ${file}
done

# Calculate gene presence/absence
panphlan_profiling.py -p panphlan/Epsilon_pangenome.tsv -i map_results --o_matrix gene_presence_absence.csv --min_coverage 1 --left_max 1.70 --right_min 0.30 --add_ref

# Normalize gene coverage across samples
perl /gpfs/data/rbeinart/bin/trinityrnaseq-v2.11.0/util/abundance_estimates_to_matrix.pl --est_method salmon --gene_trans_map none --cross_sample_norm TMM --name_sample_by_basedir --quant_files SalmonQuant.txt

