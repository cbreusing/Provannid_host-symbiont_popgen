#!/bin/bash
#SBATCH -J PangenomeReconstruction
#SBATCH -t 100:00:00
#SBATCH -n 24
#SBATCH -N 1
#SBATCH --mem=50g
#SBATCH -o PangenomeReconstruction.out
#SBATCH -e PangenomeReconstruction.err

for file in `cat filelist.txt`
do
prokka --locustag ${file} --force --outdir prokka --prefix ${file} --gcode 11 --cpus 24 ${file}.fasta
done

mkdir panaroo-qc
panaroo-qc -t 24 --graph_type all -i prokka/*.gff --ref_db /gpfs/data/rbeinart/Databases/refseq.genomes.k21s1000.msh -o panaroo-qc

mkdir panaroo
panaroo -i prokka/*.gff -o panaroo --mode strict -a core --aligner mafft -c 0.98 --core_threshold 0.95 -t 24 --refind_prop_match 0.5 --search_radius 1000 -f 0.7
