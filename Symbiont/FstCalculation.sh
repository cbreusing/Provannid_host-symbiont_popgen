#!/bin/bash
#SBATCH -J FstCalculation
#SBATCH -t 24:00:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=50g
#SBATCH -o FstCalculation.out
#SBATCH -e FstCalculation.err

python /gpfs/data/rbeinart/cbreusing/Scripts/Fst_per_snp_and_pop.py Epsilon.Freebayes.FINAL.recode.vcf Epsilon.meta.txt ABE KM > Fst_per_snp_ABE-KM.txt
python /gpfs/data/rbeinart/cbreusing/Scripts/Fst_per_snp_and_pop.py Epsilon.Freebayes.FINAL.recode.vcf Epsilon.meta.txt ABE TC > Fst_per_snp_ABE-TC.txt
python /gpfs/data/rbeinart/cbreusing/Scripts/Fst_per_snp_and_pop.py Epsilon.Freebayes.FINAL.recode.vcf Epsilon.meta.txt KM TC > Fst_per_snp_KM-TC.txt

python /gpfs/data/rbeinart/cbreusing/Scripts/Fst_per_gene.py Epsilon.Freebayes.FINAL.genes.sorted.vcf Epsilon.meta.txt ABE KM > Fst_per_gene_ABE-KM.txt
python /gpfs/data/rbeinart/cbreusing/Scripts/Fst_per_gene.py Epsilon.Freebayes.FINAL.genes.sorted.vcf Epsilon.meta.txt ABE TC > Fst_per_gene_ABE-TC.txt
python /gpfs/data/rbeinart/cbreusing/Scripts/Fst_per_gene.py Epsilon.Freebayes.FINAL.genes.sorted.vcf Epsilon.meta.txt KM TC > Fst_per_gene_KM-TC.txt
