#!/usr/bin/env bash

THREADS=16

GENOME="genome/Nothobranchius_guentheri.contigs.fsa"
PROTEINS="proteins/tldc_proteins.faa"

RNASEQ_DIR="rnaseq"

mkdir -p \
alignments \
transcripts \
orfs \
analysis \
tmp

#########################################################
# STEP 1
# MINIPROT INDEX
#########################################################

echo "STEP 1: miniprot index"

miniprot -d genome/genome.mpi $GENOME

#########################################################
# STEP 2
# MAP PROTEINS TO GENOME
#########################################################

echo "STEP 2: protein mapping"

miniprot \
--gff \
genome/genome.mpi \
$PROTEINS \
> alignments/proteins.gff

#########################################################
# STEP 3
# HISAT2 INDEX
#########################################################

echo "STEP 3: hisat2 index"

hisat2-build \
$GENOME \
genome/genome_idx

#########################################################
# STEP 4
# RNA-SEQ ALIGNMENT
#########################################################

echo "STEP 4: RNA-seq alignment"

hisat2 \
-p $THREADS \
-x genome/genome_idx \
-1 $READ1 \
-2 $READ2 \
| samtools sort -@ $THREADS \
-o alignments/rnaseq.bam

samtools index alignments/rnaseq.bam

#########################################################
# STEP 5
# BAM QC
#########################################################

echo "STEP 5: BAM statistics"

samtools flagstat \
alignments/rnaseq.bam \
> analysis/bam_stats.txt

#########################################################
# STEP 6
# TRANSCRIPT ASSEMBLY
#########################################################

echo "STEP 6: transcript assembly"

stringtie \
alignments/rnaseq.bam \
-p $THREADS \
-o transcripts/transcripts.gtf

#########################################################
# STEP 7
# EXTRACT TRANSCRIPTS FASTA
#########################################################

echo "STEP 7: transcript fasta"

gffread \
transcripts/transcripts.gtf \
-g $GENOME \
-w transcripts/transcripts.fa

#########################################################
# STEP 8
# ORF PREDICTION
#########################################################

echo "STEP 8: ORF prediction"

TransDecoder.LongOrfs \
-t transcripts/transcripts.fa

TransDecoder.Predict \
-t transcripts/transcripts.fa

#########################################################
# STEP 9
# MMSEQS DATABASES
#########################################################

echo "STEP 9: mmseqs databases"

mmseqs createdb \
$PROTEINS \
analysis/queryDB

mmseqs createdb \
transcripts/transcripts.fa.transdecoder.pep \
analysis/targetDB

#########################################################
# STEP 10
# PROTEIN SEARCH
#########################################################

echo "STEP 10: mmseqs search"

mmseqs search \
analysis/queryDB \
analysis/targetDB \
analysis/resultDB \
tmp \
--threads $THREADS \
--min-seq-id 0.7 \
-e 1e-5

#########################################################
# STEP 11
# EXPORT TSV
#########################################################

echo "STEP 11: export TSV"

mmseqs convertalis \
analysis/queryDB \
analysis/targetDB \
analysis/resultDB \
analysis/result.tsv

#########################################################
# STEP 12
# CONVERT GFF -> GTF
#########################################################

echo "STEP 12: GFF to GTF"

gffread \
alignments/proteins.gff \
-T \
-o alignments/proteins.gtf

#########################################################
# STEP 13
# CREATE BED FILES
#########################################################

echo "STEP 13: BED files"

awk '$3=="CDS"' alignments/proteins.gtf \
| awk 'BEGIN{OFS="\t"}{
print $1,$4-1,$5,$10,".",$7
}' \
| sed 's/"//g' \
| sed 's/;//g' \
> analysis/proteins.bed
