#!/bin/bash
#SBATCH --job-name=hmmsearch
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=24:00:00
#SBATCH --mem=128G
HMM_PROFILE=/mnt/hot/makulinin/nothobranchius2/PF07534.hmm

DANIO_PROT=/mnt/hot/makulinin/nothobranchius2/danio/protein.faa
DANIO_CDS=/mnt/hot/makulinin/nothobranchius2/danio/cds_from_genomic.fna

NOTHO_PROT=/mnt/hot/makulinin/nothobranchius2/notho/protein.faa
NOTHO_CDS=/mnt/hot/makulinin/nothobranchius2/notho/cds_from_genomic.fna


# searching for domains 
hmmsearch --cpu 32 --domtblout tldc_danio.domtblout $HMM_PROFILE $DANIO_PROT

# filtering by e-value
awk '$7 <= 1e-5' tldc_danio.domtblout > tldc_danio_hits_filtered.tsv

# getting ids
awk '!/^#/ {print $1}' tldc_danio_hits_filtered.tsv | sort | uniq > tldc_danio_ids.txt

# getting protein sequences
seqtk subseq $DANIO_PROT tldc_danio_ids.txt > tldc_danio_proteins.faa

# getting nucleotide sequences
awk '
BEGIN {
  while ((getline < "tldc_danio_ids.txt") > 0) ids[$1]=1
}
/^>/ {
  keep=0
  for (id in ids)
    if ($0 ~ "protein_id="id) keep=1
}
keep { print }
' $DANIO_CDS > cds_danio_selected.fna



# searching for domains 
hmmsearch --cpu 32 --domtblout tldc_notho.domtblout $HMM_PROFILE $NOTHO_PROT

# filtering by e-value
awk '$7 <= 1e-5' tldc_notho.domtblout > tldc_notho_hits_filtered.tsv

# getting ids
awk '!/^#/ {print $1}' tldc_notho_hits_filtered.tsv | sort | uniq > tldc_notho_ids.txt

# getting protein sequences
seqtk subseq $NOTHO_PROT tldc_notho_ids.txt > tldc_notho_proteins.faa

# getting nucleotide sequences
awk '
BEGIN {
  while ((getline < "tldc_notho_ids.txt") > 0) ids[$1]=1
}
/^>/ {
  keep=0
  for (id in ids)
    if ($0 ~ "protein_id="id) keep=1
}
keep { print }
' $NOTHO_CDS > cds_notho_selected.fna