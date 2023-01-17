#!/bin/bash
#SBATCH -J TranscriptomeAssembly
#SBATCH -t 100:00:00
#SBATCH -n 24
#SBATCH -N 1
#SBATCH --mem=150g
#SBATCH -o TranscriptomeAssembly.out
#SBATCH -e TranscriptomeAssembly.err

module load fastqc
module load perl
module load samtools
module load bowtie2/2.3.5.1
module load bbmap
module load jellyfish/2.2.10
module load python/3.7.4
module load cdhit
module load transdecoder/5.4.0
module load blast/2.6.0+ hmmer/3.1b2
module load seqtk

# This is an initial quality check.
fastqc -t 24 *.fastq.gz

# Provide a sample list of raw sequencing reads 
for file in `cat allreads.txt`
do

# Quality and adaptor trimming.
java -jar /gpfs/data/rbeinart/bin/trinityrnaseq-v2.11.0/trinity-plugins/Trimmomatic/trimmomatic-0.36.jar PE -threads 24 -phred33 ${file}_R1.fastq.gz ${file}_R2.fastq.gz ${file}_R1_paired.fastq ${file}_R1_unpaired.fastq ${file}_R2_paired.fastq ${file}_R2_unpaired.fastq ILLUMINACLIP:/gpfs/data/rbeinart/cbreusing/Adapters/Illumina.fa:2:30:10 SLIDINGWINDOW:4:20 LEADING:5 TRAILING:5 MINLEN:25

# These commands remove contaminants from the data, such as human DNA and the PhiX standard.
bowtie2 -p 24 -x /gpfs/data/rbeinart/Databases/contaminants -1 ${file}_R1_paired.fastq -2 ${file}_R2_paired.fastq | samtools view -bS -h -@ 24 - | samtools sort -@ 24 - > ${file}.bowtie2.sorted.bam
samtools view -@ 24 -f12 ${file}.bowtie2.sorted.bam > ${file}.unmapped.sam
cut -f1 ${file}.unmapped.sam | sort | uniq > ${file}.unmapped_ids.lst
seqtk subseq ${file}_R1_paired.fastq ${file}.unmapped_ids.lst > ${file}_R1_clean.fastq
seqtk subseq ${file}_R2_paired.fastq ${file}.unmapped_ids.lst > ${file}_R2_clean.fastq

# The contaminant free samples are then merged and used in sortmerna to remove rRNA
/gpfs/data/rbeinart/Software/sortmerna-2.1b/scripts/merge-paired-reads.sh ${file}_R1_clean.fastq ${file}_R2_clean.fastq ${file}_merged.fastq

sortmerna --ref /gpfs/data/rbeinart/Software/sortmerna-2.1b/rRNA_databases/silva-bac-16s-id90.fasta,/gpfs/data/rbeinart/Software/sortmerna-2.1b/index/silva-bac-16s-db:/gpfs/data/rbeinart/Software/sortmerna-2.1b/rRNA_databases/silva-bac-23s-id98.fasta,/gpfs/data/rbeinart/Software/sortmerna-2.1b/index/silva-bac-23s-db:/gpfs/data/rbeinart/Software/sortmerna-2.1b/rRNA_databases/silva-arc-16s-id95.fasta,/gpfs/data/rbeinart/Software/sortmerna-2.1b/index/silva-arc-16s-db:/gpfs/data/rbeinart/Software/sortmerna-2.1b/rRNA_databases/silva-arc-23s-id98.fasta,/gpfs/data/rbeinart/Software/sortmerna-2.1b/index/silva-arc-23s-db:/gpfs/data/rbeinart/Software/sortmerna-2.1b/rRNA_databases/silva-euk-18s-id95.fasta,/gpfs/data/rbeinart/Software/sortmerna-2.1b/index/silva-euk-18s-db:/gpfs/data/rbeinart/Software/sortmerna-2.1b/rRNA_databases/silva-euk-28s-id98.fasta,/gpfs/data/rbeinart/Software/sortmerna-2.1b/index/silva-euk-28s:/gpfs/data/rbeinart/Software/sortmerna-2.1b/rRNA_databases/rfam-5s-database-id98.fasta,/gpfs/data/rbeinart/Software/sortmerna-2.1b/index/rfam-5s-db:/gpfs/data/rbeinart/Software/sortmerna-2.1b/rRNA_databases/rfam-5.8s-database-id98.fasta,/gpfs/data/rbeinart/Software/sortmerna-2.1b/index/rfam-5.8s-db --reads ${file}_merged.fastq --sam --best 1 --min_lis 2 --fastx --paired_in --aligned ${file}_rRNA --other ${file}_non_rRNA --log -a 24 -v 

# Separate host and symbiont reads with BBSplit
bbsplit.sh -Xmx50g ref=/gpfs/data/rbeinart/cbreusing/Genomes/Gamma1.fasta,/gpfs/data/rbeinart/cbreusing/Genomes/GammaLau.fna,/gpfs/data/rbeinart/cbreusing/Genomes/Epsilon.fasta,/gpfs/data/rbeinart/cbreusing/Genomes/Ifr_SOX.fasta,/gpfs/data/rbeinart/cbreusing/Genomes/Ifr_MOX.fasta build=1 path=/gpfs/data/rbeinart/cbreusing/Genomes/ in=${file}_non_rRNA.fastq ambiguous=best ambiguous2=toss basename=${file}_%_#.fq outu=${file}_unmapped_#.fq out=${file}_bbmap.sam

done

# Provide a sample list of separated host reads
for file in `cat hostreads.txt`
do

# These commands correct sequencing errors and remove flagged overrepresented sequences
perl /gpfs/data/rbeinart/Software/rcorrector/run_rcorrector.pl -1 ${file}_1.fq -2 ${file}_2.fq -t 24
python /gpfs/data/rbeinart/Software/TranscriptomeAssemblyTools/FilterUncorrectabledPEfastq.py -1 ${file}_1.cor.fq -2 ${file}_2.cor.fq -s ${file}
fastqc -t 24 unfixrm_${file}_1.cor.fq unfixrm_${file}_2.cor.fq
unzip unfixrm_${file}_1.cor_fastqc.zip
unzip unfixrm_${file}_2.cor_fastqc.zip
python /gpfs/data/rbeinart/Software/TranscriptomeAssemblyTools/RemoveFastqcOverrepSequenceReads.py -1 unfixrm_${file}_1.cor.fq -2 unfixrm_${file}_2.cor.fq -fql unfixrm_${file}_1.cor_fastqc/fastqc_data.txt -fqr unfixrm_${file}_2.cor_fastqc/fastqc_data.txt

done

# Assembly
/gpfs/data/rbeinart/bin/trinityrnaseq-v2.11.0/Trinity --seqType fq --max_memory 150G --samples_file A_boucheti.txt --SS_lib_type RF --output A_boucheti.PasaFly.Trinity --bfly_algorithm PASAFLY --CPU 24 --full_cleanup
/gpfs/data/rbeinart/bin/trinityrnaseq-v2.11.0/Trinity --seqType fa --max_memory 100G --single A_boucheti_KM433.fasta,A_boucheti_TC432.fasta --output A_boucheti.454.Trinity --bfly_algorithm PASAFLY --no_normalize_reads --CPU 24 --full_cleanup

cat A_boucheti*Trinity*fasta > A_boucheti.Trinity.combined.fasta

# Transcript clustering to remove redundancies
cd-hit-est -i A_boucheti.Trinity.combined.fasta -o A_boucheti.Trinity.merged95.fasta -c 0.95 -n 10 -r 1 -T 24 -d 0 -M 50000 -g 1 -aS 0.95

# This command is just for quality control to see how well the assembly reflects the original input reads
bowtie2 -p 24 -q --no-unal -x A_boucheti.Trinity.merged95.fasta -1 rmoverrep_unfixrm_109_unmapped_1.cor.fq,rmoverrep_unfixrm_119_unmapped_1.cor.fq,rmoverrep_unfixrm_145_unmapped_1.cor.fq,rmoverrep_unfixrm_159_unmapped_1.cor.fq,rmoverrep_unfixrm_2_unmapped_1.cor.fq,rmoverrep_unfixrm_10_unmapped_1.cor.fq,rmoverrep_unfixrm_121_unmapped_1.cor.fq,rmoverrep_unfixrm_147_unmapped_1.cor.fq,rmoverrep_unfixrm_241_unmapped_1.cor.fq,rmoverrep_unfixrm_4_unmapped_1.cor.fq,rmoverrep_unfixrm_115_unmapped_1.cor.fq,rmoverrep_unfixrm_123_unmapped_1.cor.fq,rmoverrep_unfixrm_149_unmapped_1.cor.fq,rmoverrep_unfixrm_243_unmapped_1.cor.fq,rmoverrep_unfixrm_117_unmapped_1.cor.fq,rmoverrep_unfixrm_12_unmapped_1.cor.fq,rmoverrep_unfixrm_14_unmapped_1.cor.fq,rmoverrep_unfixrm_24_unmapped_1.cor.fq -2 rmoverrep_unfixrm_109_unmapped_2.cor.fq,rmoverrep_unfixrm_119_unmapped_2.cor.fq,rmoverrep_unfixrm_145_unmapped_2.cor.fq,rmoverrep_unfixrm_159_unmapped_2.cor.fq,rmoverrep_unfixrm_2_unmapped_2.cor.fq,rmoverrep_unfixrm_10_unmapped_2.cor.fq,rmoverrep_unfixrm_121_unmapped_2.cor.fq,rmoverrep_unfixrm_147_unmapped_2.cor.fq,rmoverrep_unfixrm_241_unmapped_2.cor.fq,rmoverrep_unfixrm_4_unmapped_2.cor.fq,rmoverrep_unfixrm_115_unmapped_2.cor.fq,rmoverrep_unfixrm_123_unmapped_2.cor.fq,rmoverrep_unfixrm_149_unmapped_2.cor.fq,rmoverrep_unfixrm_243_unmapped_2.cor.fq,rmoverrep_unfixrm_117_unmapped_2.cor.fq,rmoverrep_unfixrm_12_unmapped_2.cor.fq,rmoverrep_unfixrm_14_unmapped_2.cor.fq,rmoverrep_unfixrm_24_unmapped_2.cor.fq 2>align_stats.txt

# ORF prediction
TransDecoder.LongOrfs -t A_boucheti.Trinity.merged95.fasta

blastp -query A_boucheti.Trinity.merged95.fasta.transdecoder_dir/longest_orfs.pep -db /gpfs/data/rbeinart/Databases/uniref90.fasta -max_target_seqs 1 -outfmt "6 std stitle" -evalue 1e-5 -num_threads 24 > blastp.outfmt6
/gpfs/data/rbeinart/bin/trinityrnaseq-v2.11.0/util/misc/blast_outfmt6_group_segments.pl blastp.outfmt6 A_boucheti.Trinity.merged95.fasta.transdecoder_dir/longest_orfs.pep /gpfs/data/rbeinart/Databases/uniref90.fasta > blastp.outfmt6.grouped
hmmscan --cpu 24 --domtblout pfam.domtblout /gpfs/data/rbeinart/Databases/Pfam-A.hmm A_boucheti.Trinity.merged95.fasta.transdecoder_dir/longest_orfs.pep

TransDecoder.Predict -t A_boucheti.Trinity.merged95.fasta --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6

# Removal of remaining bacterial contaminants with BlobTools
cat rmoverrep_unfixrm*_unmapped.1.cor.fq > A_boucheti_1.fq
cat rmoverrep_unfixrm*_unmapped.2.cor.fq > A_boucheti_2.fq

bowtie2-build A_boucheti.Trinity.merged95.fasta.transdecoder.cds A_boucheti.Trinity.merged95.fasta.transdecoder.cds
bowtie2 --very-sensitive -p 24 -x A_boucheti.Trinity.merged95.fasta.transdecoder.cds -1 A_boucheti_1.fq -2 A_boucheti_2.fq | samtools view -bS -h -@ 24 - | samtools sort -@ 24 - > A_boucheti.sorted.bam
samtools view -h -@ 24 -o A_boucheti.sam A_boucheti.sorted.bam 

blobtools create -i A_boucheti.Trinity.merged95.fasta.transdecoder.cds -o A_boucheti -s A_boucheti.sam -t A_boucheti_hits.txt
blobtools view -i A_boucheti.blobDB.json -r all
blobtools plot -i A_boucheti.blobDB.json --notitle -r superkingdom --format pdf --colours colors.txt
grep "Eukaryota" A_boucheti.blobDB.table.txt > seqs_to_keep.txt
grep "no-hit" A_boucheti.blobDB.table.txt >> seqs_to_keep.txt
perl -anle 'print $F[0]' seqs_to_keep.txt > euk_seqs.txt
seqtk subseq A_boucheti.Trinity.merged95.fasta.transdecoder.cds euk_seqs.txt > A_boucheti.Trinity.merged95.filtered.cds

# Assembly statistics
perl /gpfs/data/rbeinart/bin/trinityrnaseq-v2.11.0/util/TrinityStats.pl A_boucheti.Trinity.merged95.filtered.cds > TrinityStats.txt

# Assembly completeness
source activate busco

busco -i A_boucheti.Trinity.merged95.filtered.cds -l mollusca_odb10 -o busco_mollusca_merged -m tran -c 24
