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
samtools view -h -@ 24 -o ${file}_subsampled.sam ${file}_subsampled.bam
cut -f1 ${file}_subsampled.sam | sort | uniq > ${file}_ids.lst
seqtk subseq ${file}_R1_clean.fastq ${file}_ids.lst > ${file}_R1.fastq
seqtk subseq ${file}_R2_clean.fastq ${file}_ids.lst > ${file}_R2.fastq
reformat.sh in=${file}_R1.fastq in2=${file}_R2.fastq out=${file}.fastq
# Map interleaved reads to (pan)genome reference for presence/absence analysis (PanPhlan)
panphlan_map.py -i ${file}.fastq -o map_results/${file}_mapping.csv --nproc 24 --indexes pangenome/Epsilon -p pangenome/Epsilon_pangenome_RAST.tsv
# Map reads against protein-coding gene set for gene abundance analysis (Salmon)
bbmap.sh ref=Epsilon.ffn nodisk in=${file}_R1.fastq in2=${file}_R2.fastq out=${file}.sam sam=1.4 local=t requirecorrectstrand=t samestrandpairs=f xstag=fs
salmon quant -t Epsilon.ffn -l IU --meta --rangeFactorizationBins 4 --numBootstraps 1000 --seqBias --gcBias -s -u -a ${file}.sam -o ${file}
done

# Calculate gene presence/absence
panphlan_profiling.py -p pangenome/Epsilon_pangenome_RAST.tsv -i map_results --o_matrix gene_presence_absence.csv --min_coverage 1 --left_max 1.70 --right_min 0.30 --o_covmat coverage_matrix.csv --add_ref --o_covplot_normed Epsilon_covplot

# Normalize gene coverage across samples
perl /gpfs/data/rbeinart/bin/trinityrnaseq-v2.11.0/util/abundance_estimates_to_matrix.pl --est_method salmon --gene_trans_map none --cross_sample_norm TMM --name_sample_by_basedir --quant_files SalmonQuant.txt

