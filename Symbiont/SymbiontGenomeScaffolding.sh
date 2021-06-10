#!/bin/bash
#SBATCH -J SymbiontGenomeScaffolding
#SBATCH -t 48:00:00
#SBATCH -n 24
#SBATCH -N 1
#SBATCH --mem=200g
#SBATCH -o SymbiontGenomeScaffolding.out
#SBATCH -e SymbiontGenomeScaffolding.err

module load samtools
module load bowtie
module load perl
module load bbmap
module load bwa
module load idba
module load bamtools
module load velvet
module load bowtie2/2.3.5.1
module load python/3.7.4

perl /gpfs/data/rbeinart/Software/GapFiller/GapFiller.pl -l libAlv337.txt -s metaSpades.bin.fasta -T 24 -b Alv337.polished
/gpfs/data/rbeinart/Software/LR_Gapcloser/src/LR_Gapcloser.sh -i Alv337.polished/Alv337.polished.gapfilled.final.fa -l Alv337_nanopore_trim.fasta -s n -t 24 -r 10 -o Alv337.gapclosed 
ra2.py -i Alv337.gapclosed/iteration-10/gapclosed.fasta -1 Alv337_Epsilon_1.fq -2 Alv337_Epsilon_2.fq -t 24
fix_fasta.py gapclosed.curated/re_assembled.fa | nr_fasta.py rename - > Alv337.curated.fa
perl /gpfs/data/rbeinart/Software/SSPACE/SSPACE_Standard_v3.0.pl -l libAlv337.txt -s Alv337.curated.fa -x 1 -T 24 -b Alv337.curated
bwa index Alv337.curated/Alv337.curated.final.scaffolds.fasta
bwa mem -a Alv337.curated/Alv337.curated.final.scaffolds.fasta Alv337.curated/Alv337.curated.final.scaffolds.fasta > align-self.sam
samtools view -Sb align-self.sam > align-self.bam
bwa mem -t24 -k11 -W20 -r10 -A1 -B1 -O1 -E1 -L0 -a -Y Alv337.curated/Alv337.curated.final.scaffolds.fasta Alv337_symbiont.nanopore.fasta > align.sam
samtools view -Sb align.sam > align.bam
SLR -c Alv337.curated/Alv337.curated.final.scaffolds.fasta -r align.bam -d align-self.bam -p SLR -n 2000 -x 50
perl /gpfs/data/rbeinart/Software/SSPACE-LongRead/SSPACE-LongRead.pl -c SLR/unique-contig-set.fa -p Alv337_symbiont.nanopore.fastq -b Alv337.longread
bwa index Alv337.longread/scaffolds.fasta
bwa mem Alv337.longread/scaffolds.fasta SLR/unique-contig-set.fa > unique-contig.sam
samtools view -Sb unique-contig.sam > unique-contig.bam
SLR-unique-ambiguous -c Alv337.curated/Alv337.curated.final.scaffolds.fasta -r align.bam -d align-self.bam -u SLR/unique-contig-set.fa -s Alv337.longread/scaffolds.fasta -b unique-contig.bam -p SLR_unique_ambiguous

# Iterate above steps up to ten times and proceed with best scaffolding result from either SLR or SLR-unique-ambiguous  

mv SLR_unique_ambiguous/scaffold_set.fa Epsilon_consensus.fasta

bwa index Epsilon_consensus.fasta
bowtie2-build Epsilon_consensus.fasta Epsilon_consensus.fasta
bowtie2 -p 24 -x Epsilon_consensus.fasta -1 Alv337_Epsilon_1.fq -2 Alv337_Epsilon_2.fq | samtools view -bS -h -@ 24 - | samtools sort -@ 24 - > Epsilon.bowtie2.sorted.bam
samtools index Epsilon.bowtie2.sorted.bam
bwa mem -t24 -k11 -W20 -r10 -A1 -B1 -O1 -E1 -L0 -a -Y Epsilon_consensus.fasta Alv337_nanopore_trim.fasta > Epsilon.bwa.sam
samtools view -Sb -h -@ 24 Epsilon.bwa.sam | samtools sort -@ 24 - > Epsilon.bwa.sorted.bam
samtools index Epsilon.bwa.sorted.bam
pilon -Xmx200G --genome Epsilon_consensus.fasta --frags Epsilon.bowtie2.sorted.bam --unpaired Epsilon.bwa.sorted.bam --fix snps,indels,amb --output Epsilon --outdir pilon --threads 24

