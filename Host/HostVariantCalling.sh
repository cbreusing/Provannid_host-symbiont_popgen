#!/bin/bash
#SBATCH -J HostVariantCalling
#SBATCH -t 48:00:00
#SBATCH -n 24
#SBATCH -N 1
#SBATCH --mem=50g
#SBATCH -o HostVariantCalling.out
#SBATCH -e HostVariantCalling.err

module load bowtie2/2.3.5.1
module load samtools
module load gatk/4.1.6.0
module load picard-tools/2.17.11
module load R

# Index the host transcriptome reference (renamed from the original file for brevity)
bowtie2-build A_boucheti_reference.fa A_boucheti_reference.fa
java -Xmx2g -jar /gpfs/runtime/opt/picard-tools/2.17.11/picard.jar CreateSequenceDictionary R=A_boucheti_reference.fa O=A_boucheti_reference.dict
samtools faidx A_boucheti_reference.fa

# Map host genomic reads against reference, mark duplicates, realign around indel regions and recalibrate base scores
for file in `cat filelist.txt`
do
bowtie2 --very-sensitive -p 24 -x A_boucheti_reference.fa -1 ../${file}_unmapped_1.fq -2 ../${file}_unmapped_2.fq | samtools view -bS -h -@ 24 - | samtools sort -@ 24 - > ${file}.sorted.bam
java -Xmx2g -jar /gpfs/runtime/opt/picard-tools/2.17.11/picard.jar MarkDuplicates I=${file}.sorted.bam O=${file}.dedup.bam M=${file}.metrics.txt
samtools index ${file}.dedup.bam
lofreq viterbi -f A_boucheti_reference.fa -k ${file}.dedup.bam | samtools sort -@ 24 - > ${file}.realigned.bam
lofreq indelqual -f A_boucheti_reference.fa --dindel -o ${file}.indelqual.bam ${file}.realigned.bam
samtools index ${file}.indelqual.bam
done

ls *indelqual.bam > bam.list

# Get an overview of the read depth distribution (choose minimum and maximum values accordingly in downstream filtering)
angsd -P 24 -bam bam.list -ref A_boucheti_reference.fa -out ANGSD/A_boucheti.qc -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -C 50 -baq 1 -minMapQ 30 -doQsDist 1 -doDepth 1 -doCounts 1 -maxDepth 500
Rscript /gpfs/data/rbeinart/Software/ngsTools/Scripts/plotQC.R ANGSD/A_boucheti.qc

# Calculate inbreeding coefficients and site frequency spectra for each host population
for POP in TC KM ABE NS;
do
IND=`cat ${POP}.list | wc -l`
angsd -P 24 -bam ${POP}.list -ref A_boucheti_reference.fa -gl 1 -baq 1 -C 50 -minInd ${IND} -minMapQ 30 -minQ 20 -setMinDepth 2 -setMaxDepth 40 -doCounts 1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -doMaf 1 -doMajorMinor 1 -SNP_pval 1e-6 -doGlf 3 -minMaf 0.01 -skipTriallelic 1 -out ANGSD/${POP}
NSITES=`zcat ANGSD/${POP}.mafs.gz | tail -n+2 | wc -l`
zcat ANGSD/${POP}.glf.gz > ANGSD/${POP}.glf
/gpfs/data/rbeinart/Software/ngsTools/ngsF/ngsF.sh --n_ind ${IND} --n_sites ${NSITES} --glf ANGSD/${POP}.glf --out ANGSD/${POP}.indF
angsd -P 24 -bam ${POP}.list -ref A_boucheti_reference.fa -anc A_boucheti_reference.fa -gl 1 -baq 1 -C 50 -minInd ${IND} -minMapQ 30 -minQ 20 -setMinDepth 2 -setMaxDepth 40 -doCounts 1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -doMaf 1 -doMajorMinor 4 -doSaf 2 -indF ANGSD/${POP}.indF -out ANGSD/${POP}
/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/realSFS ANGSD/${POP}.saf.idx -tole 1e-6 -maxIter 5000 -P 24 -fold 1 > ANGSD/$POP.sfs
angsd -P 24 -bam ${POP}.list -ref A_boucheti_reference.fa -anc A_boucheti_reference.fa -gl 1 -baq 1 -C 50 -minInd ${IND} -minMapQ 30 -minQ 20 -setMinDepth 2 -setMaxDepth 40 -doCounts 1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -doMaf 1 -doMajorMinor 4 -doSaf 2 -indF ANGSD/${POP}.indF -doThetas 1 -pest ANGSD/${POP}.sfs -out ANGSD/${POP}
done 

# Calculate pairwise FSTs between populations
/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/realSFS ANGSD/TC.saf.idx ANGSD/KM.saf.idx -tole 1e-6 -maxIter 5000 -fold 1 -P 24 > ANGSD/TC-KM.folded.sfs
/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/realSFS ANGSD/TC.saf.idx ANGSD/ABE.saf.idx -tole 1e-6 -maxIter 5000 -fold 1 -P 24 > ANGSD/TC-ABE.folded.sfs
/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/realSFS ANGSD/TC.saf.idx ANGSD/NS.saf.idx -tole 1e-6 -maxIter 5000 -fold 1 -P 24 > ANGSD/TC-NS.folded.sfs
/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/realSFS ANGSD/KM.saf.idx ANGSD/ABE.saf.idx -tole 1e-6 -maxIter 5000 -fold 1 -P 24 > ANGSD/KM-ABE.folded.sfs
/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/realSFS ANGSD/KM.saf.idx ANGSD/NS.saf.idx -tole 1e-6 -maxIter 5000 -fold 1 -P 24 > ANGSD/KM-NS.folded.sfs
/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/realSFS ANGSD/ABE.saf.idx ANGSD/NS.saf.idx -tole 1e-6 -maxIter 5000 -fold 1 -P 24 > ANGSD/ABE-NS.folded.sfs

/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/realSFS fst index ANGSD/TC.saf.idx ANGSD/KM.saf.idx -sfs ANGSD/TC-KM.folded.sfs -fold 1 -fstout ANGSD/TC-KM -whichFST 1 -maxIter 5000
/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/realSFS fst stats ANGSD/TC-KM.fst.idx -whichFST 1 -maxIter 5000 > ANGSD/TC-KM.fst.txt
/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/realSFS fst index ANGSD/TC.saf.idx ANGSD/ABE.saf.idx -sfs ANGSD/TC-ABE.folded.sfs -fold 1 -fstout ANGSD/TC-ABE -whichFST 1 -maxIter 5000
/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/realSFS fst stats ANGSD/TC-ABE.fst.idx -whichFST 1 -maxIter 5000 > ANGSD/TC-ABE.fst.txt
/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/realSFS fst index ANGSD/TC.saf.idx ANGSD/NS.saf.idx -sfs ANGSD/TC-NS.folded.sfs -fold 1 -fstout ANGSD/TC-NS -whichFST 1 -maxIter 5000
/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/realSFS fst stats ANGSD/TC-NS.fst.idx -whichFST 1 -maxIter 5000 > ANGSD/TC-NS.fst.txt
/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/realSFS fst index ANGSD/KM.saf.idx ANGSD/ABE.saf.idx -sfs ANGSD/KM-ABE.folded.sfs -fold 1 -fstout ANGSD/KM-ABE -whichFST 1 -maxIter 5000
/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/realSFS fst stats ANGSD/KM-ABE.fst.idx -whichFST 1 -maxIter 5000 > ANGSD/KM-ABE.fst.txt
/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/realSFS fst index ANGSD/ABE.saf.idx ANGSD/NS.saf.idx -sfs ANGSD/ABE-NS.folded.sfs -fold 1 -fstout ANGSD/ABE-NS -whichFST 1 -maxIter 5000
/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/realSFS fst stats ANGSD/ABE-NS.fst.idx -whichFST 1 -maxIter 5000 > ANGSD/ABE-NS.fst.txt 
/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/realSFS fst index ANGSD/KM.saf.idx ANGSD/NS.saf.idx -sfs ANGSD/KM-NS.folded.sfs -fold 1 -fstout ANGSD/KM-NS -whichFST 1 -maxIter 5000
/gpfs/data/rbeinart/Software/ngsTools/angsd/misc/realSFS fst stats ANGSD/KM-NS.fst.idx -whichFST 1 -maxIter 5000 > ANGSD/KM-NS.fst.txt

# Create a list of individual inbreeding coefficients to correct for deviations from HWE and run the final genotype likelihood estimation
cat ANGSD/TC.indF ANGSD/KM.indF ANGSD/ABE.indF ANGSD/NS.indF > ANGSD/A_boucheti.indF
angsd -P 24 -bam bam.list -ref A_boucheti_reference.fa -gl 1 -baq 1 -C 50 -minInd 12 -minMapQ 30 -minQ 20 -setMinDepth 2 -setMaxDepth 40 -doCounts 1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -doMaf 1 -doMajorMinor 1 -SNP_pval 1e-6 -dosnpstat 1 -doHWE 1 -sb_pval 0.05 -hetbias_pval 0.05 -makeMatrix 1 -doIBS 1 -doCov 1 -doPost 1 -doGeno 2 -minMaf 0.01 -indF ANGSD/A_boucheti.indF -skipTriallelic 1 -out ANGSD/A_boucheti

# Repeat analysis using only sites that are highly differentiated between populations (based on OutFLANK results)
sort -k1 sites.txt > sites_sorted.txt
angsd sites index sites_sorted.txt 
angsd -P 24 -bam bam.list -ref A_boucheti_reference.fa -sites sites_sorted.txt -gl 1 -baq 1 -C 50 -minInd 12 -minMapQ 30 -minQ 20 -setMinDepth 2 -setMaxDepth 40 -doCounts 1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -doMaf 1 -doMajorMinor 1 -SNP_pval 1e-6 -dosnpstat 1 -doHWE 1 -sb_pval 0.05 -hetbias_pval 0.05 -makeMatrix 1 -doIBS 1 -doCov 1 -doPost 1 -doGeno 2 -minMaf 0.01 -indF ANGSD/A_boucheti.indF -skipTriallelic 1 -out ANGSD/A_boucheti.highFST
