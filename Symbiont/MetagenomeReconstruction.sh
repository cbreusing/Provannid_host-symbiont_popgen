#!/bin/bash
#SBATCH -J MAGReconstruction
#SBATCH -t 100:00:00
#SBATCH -n 24
#SBATCH -N 1
#SBATCH --mem=100g
#SBATCH -o MAGReconstruction.out
#SBATCH -e MAGReconstruction.err

module load samtools/1.12
module load bowtie2/2.3.5.1
module load hmmer
module load idba
module load prodigal
module load R
module load ruby

for file in `cat filelist.txt`
do
metaspades.py -1 ${file}_R1_clean.fastq -2 ${file}_R2_clean.fastq -k 21,31,41,51,61,71,81,91,101,111,121 -t 24 -o metaSpades_${file}
cd metaSpades_${file}
bowtie2-build scaffolds.fasta scaffolds.fasta
bowtie2 -p 24 -x scaffolds.fasta -1 ../${file}_R1_clean.fastq -2 ../${file}_R2_clean.fastq | samtools view -bS -h -@ 24 - | samtools sort -@ 24 - > ${file}.bowtie2.sorted.bam
jgi_summarize_bam_contig_depths --outputDepth depth.txt ${file}.bowtie2.sorted.bam
metabat2 -i scaffolds.fasta -a depth.txt -m 1500 -o metabat/${file} --unbinned -t 24
python /gpfs/data/rbeinart/cbreusing/miniconda3/envs/graphbin/bin/prepresult.py --binned metabat --output metabat --prefix ${file}
graphbin --graph assembly_graph_with_scaffolds.gfa --binned metabat/${file}_initial_contig_bins.csv --output GraphBin --prefix ${file} --assembler SPAdes --paths scaffolds.paths --contigs scaffolds.fasta
mkdir maxbin
run_MaxBin.pl -contig scaffolds.fasta -reads ../${file}_R1_clean.fastq -reads2 ../${file}_R2_clean.fastq -out maxbin/${file} -thread 24 -min_contig_length 500
/gpfs/data/rbeinart/Software/DAS_Tool/src/Fasta_to_Scaffolds2Bin.sh -i maxbin -e fasta > maxbin.scaffolds2bin.tsv
/gpfs/data/rbeinart/Software/DAS_Tool/src/Fasta_to_Scaffolds2Bin.sh -i metabat -e fa > metabat.scaffolds2bin.tsv
/gpfs/data/rbeinart/Software/DAS_Tool/src/Fasta_to_Scaffolds2Bin.sh -i GraphBin/${file}_bins -e fasta > graphbin.scaffolds2bin.tsv
DAS_Tool -i maxbin.scaffolds2bin.tsv,metabat.scaffolds2bin.tsv,graphbin.scaffolds2bin.tsv -l MaxBin2,MetaBAT2,GraphBin --score_threshold 0 -c scaffolds.fasta -o DAS_Tool/${file} --write_bins 1 -t 24
cd ..
done
