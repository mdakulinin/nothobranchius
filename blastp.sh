#!/bin/bash
#SBATCH --job-name=blastp
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=24:00:00
#SBATCH --mem=128G

BLAST_DB=/mnt/hot/makulinin/nothobranchius2/danio_base/tldc_danio
NOTHO_PROT=/mnt/hot/makulinin/nothobranchius2/notho/protein.faa

blastp \
  -query $NOTHO_PROT \
  -db $BLAST_DB \
  -out notho_vs_danio_tldc.blastp.tsv \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle" \
  -evalue 1e-5 \
  -num_threads 32

cut -f1 notho_vs_danio_tldc.blastp.tsv | sort -u > ids.txt

awk 'BEGIN{while((getline<"ids.txt")>0) ids[$1]=1}
     /^>/{
        header=$0;
        id=substr($1,2);
        keep = (id in ids)
     }
     keep{print}' \
$NOTHO_PROT \
> notho_tldc_hits.faa


