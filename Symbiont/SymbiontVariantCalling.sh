#!/bin/bash
#SBATCH -J Freebayes
#SBATCH -t 48:00:00
#SBATCH -n 24
#SBATCH -N 1
#SBATCH --mem=50g
#SBATCH -o Freebayes.out
#SBATCH -e Freebayes.err

module load bowtie2/2.3.5.1
module load samtools
module load gatk/4.1.6.0
module load picard-tools/2.17.11
module load R
module load vcftools
module load bcftools

#Build reference indices and dictionary.
bowtie2-build Epsilon.fasta Epsilon.fasta
java -Xmx2g -jar /gpfs/runtime/opt/picard-tools/2.17.11/picard.jar CreateSequenceDictionary R=Epsilon.fasta O=Epsilon.dict
samtools faidx Epsilon.fasta

# Provide a list with filtered symbiont samples
for file in `cat subset.list`
do
#Map reads to reference
bowtie2 --very-sensitive -p 24 -x Epsilon.fasta -1 ${file}_Epsilon_1.fq -2 ${file}_Epsilon_2.fq | samtools view -bS -h -@ 24 - | samtools sort -@ 24 - > ${file}.Epsilon.sorted.bam
#Mark and remove duplicates
java -Xmx2g -jar /gpfs/runtime/opt/picard-tools/2.17.11/picard.jar MarkDuplicates I=${file}.Epsilon.sorted.bam O=${file}.Epsilon.dedup.bam M=${file}.Epsilon.metrics.txt REMOVE_DUPLICATES=true
samtools index ${file}.Epsilon.dedup.bam
#Indel realignment
lofreq viterbi -f Epsilon.fasta -k ${file}.Epsilon.dedup.bam | samtools sort -@ 24 - > ${file}.realigned.bam
#Base recalibration
lofreq indelqual -f Epsilon.fasta --dindel -o ${file}.indelqual.bam ${file}.realigned.bam
samtools index ${file}.indelqual.bam
samtools view -bS -h -F4 ${file}.indelqual.bam > ${file}.filt.dedup.bam
#Check number of alignments for each BAM file, select the lowest number of alignments for downsampling, in this case it was 199000
samtools view -c ${file}.filt.dedup.bam >> alignment_counts.txt
done

num=199000

for file in `cat subset.list`
do
count=`samtools view -c ${file}.filt.dedup.bam`
if [ $num -le $count ]
    then
    frac=`bc -l <<< $num/$count`
    samtools view -h -bs $frac ${file}.filt.dedup.bam > ${file}_subsampled.bam
fi
samtools index ${file}_subsampled.bam
java -Xmx2g -jar /gpfs/runtime/opt/picard-tools/2.17.11/picard.jar AddOrReplaceReadGroups I=${file}_subsampled.bam O=${file}_subsampled.RG.bam RGID=${file} RGLB=LIB_${file} RGPL=ILLUMINA RGPU=FLOWCELL1 RGSM=${file} VALIDATION_STRINGENCY=SILENT
done

cat *_subsampled.RG.bam > bam.list

# Run variant calling analysis using parameters for metagenomic data
freebayes -f Epsilon.fasta -L bam.list -v Epsilon.Freebayes.vcf -F 0.01 -C 1 -p 1 --pooled-continuous -g 1000 -P 0.05 -m 30 -q 20 --min-coverage 10 --haplotype-length 0 --report-monomorphic

cat Epsilon.Freebayes.vcf | vcf-sort -c > Epsilon.Freebayes.sorted.vcf
bgzip -c Epsilon.Freebayes.sorted.vcf > Epsilon.Freebayes.vcf.gz
tabix -p vcf Epsilon.Freebayes.vcf.gz
# Filter VCF file based on strand bias, position around indels, base quality and read depth and exclude sites/individuals with more than 25% missing data
bcftools filter -g 5 -i 'REF!="N" && SRP > 5 && SAP > 5 && EPP > 5 && QUAL > 20 && INFO/DP > 10' -o Epsilon.Freebayes.filtered.vcf -O v Epsilon.Freebayes.vcf.gz
vcftools --vcf Epsilon.Freebayes.filtered.vcf --max-missing 0.75 --remove-indv B16_897 --remove-indv B16_812 --remove-indv B16_808 --recode --recode-INFO-all --out Epsilon.Freebayes
vcftools --vcf Epsilon.Freebayes.recode.vcf --max-meanDP 25 --recode --recode-INFO-all --out Epsilon.Freebayes.FINAL
# Export haplotype and allele count information
vcftools --vcf Epsilon.Freebayes.FINAL.recode.vcf --extract-FORMAT-info GT --out Epsilon.Freebayes.FINAL
gatk VariantsToTable -V Epsilon.Freebayes.FINAL.recode.vcf -O Epsilon.Freebayes.FINAL.RD.FORMAT -F CHROM -F POS -F REF -ASGF RO
gatk VariantsToTable -V Epsilon.Freebayes.FINAL.recode.vcf -O Epsilon.Freebayes.FINAL.AD.FORMAT -F CHROM -F POS -F ALT -ASGF AO --split-multi-allelic


