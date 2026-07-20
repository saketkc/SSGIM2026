#!/bin/bash
# Fetch everything 01_BulkDE.qmd needs. Run from bulk/. Resumable; re-run to continue.
set -eo pipefail
cd "$(dirname "$0")"

if command -v wget >/dev/null 2>&1; then
  DL() { wget -c -O "$2" "$1"; }
else
  DL() { curl -L -C - --fail -o "$2" "$1"; }
fi

# 1. Per-sample counts. seqout fetches GSE135251_RAW.tar but does not extract it.
echo ">>> [1/6] per-sample counts"
mkdir -p GSE135251
curl -sS "https://seqout.org/api/project/GSE135251/supplementary/download" | bash
if [ -f GSE135251/GSE135251_RAW.tar ] && [ ! -f GSE135251/GSE135251_RAW.tar.extracted ]; then
  tar -xf GSE135251/GSE135251_RAW.tar -C GSE135251/
  touch GSE135251/GSE135251_RAW.tar.extracted
fi

# 2. Sample metadata
echo ">>> [2/6] metadata"
DL "https://www.dropbox.com/scl/fi/uk6o1s2cegl77gcma7r81/GSE135251_samples.csv?rlkey=2nzbodzezbpx5qm79t35stwr8&dl=1" GSE135251_samples.csv

# 3. cdna transcriptome; t2g.csv built from its fasta headers
echo ">>> [3/6] cdna + t2g"
mkdir -p cdna
DL "https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz" cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
if [ ! -s cdna/t2g.csv ]; then
  zcat cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz | grep '>' | \
    awk 'BEGIN{FS=" "; print "tx_id,gene_id,gene_name"}
         {print substr($1,2) "," substr($4,6) "," substr($7,13)}' > cdna/t2g.csv
fi

# Shortcut: replace steps 4-6 (~3 GB + several min) with the prebuilt kallisto output.
# echo ">>> prebuilt kallisto output"
# DL "https://www.dropbox.com/scl/fi/gbhh42dp5qrtf4n3p5xvc/kallisto_output.zip?rlkey=jvnz6j40uxn2926s1z51ulex4&dl=1" kallisto_output.zip
# mkdir -p rnaseq_data && unzip -o kallisto_output.zip -d rnaseq_data/

# 4. Reads for the kallisto walk-through (~2.8 GB)
echo ">>> [4/6] fastq"
mkdir -p rnaseq_data/fastq
DL "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR988/001/SRR9883171/SRR9883171_1.fastq.gz" rnaseq_data/fastq/SRR9883171_1.fastq.gz
DL "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR988/001/SRR9883171/SRR9883171_2.fastq.gz" rnaseq_data/fastq/SRR9883171_2.fastq.gz

# 5. kallisto (Linux build; on macOS `mamba install kallisto` and skip this)
echo ">>> [5/6] kallisto"
if [ ! -x kallisto/kallisto ]; then
  DL "https://github.com/pachterlab/kallisto/releases/download/v0.52.0/kallisto_linux-v0.52.0.tar.gz" kallisto_linux-v0.52.0.tar.gz
  tar -zxf kallisto_linux-v0.52.0.tar.gz
fi
KALLISTO="./kallisto/kallisto"
command -v kallisto >/dev/null 2>&1 && KALLISTO="kallisto"

# 6. index + quant -> rnaseq_data/kallisto_output/abundance.tsv
echo ">>> [6/6] kallisto index + quant"
mkdir -p rnaseq_data/kallisto_output
[ -s human_cdna.idx ] || "$KALLISTO" index -t 8 -i human_cdna.idx cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
[ -s rnaseq_data/kallisto_output/abundance.tsv ] || "$KALLISTO" quant -i human_cdna.idx -o rnaseq_data/kallisto_output/ -t 8 \
  rnaseq_data/fastq/SRR9883171_1.fastq.gz rnaseq_data/fastq/SRR9883171_2.fastq.gz

echo ""
echo ">>> DONE. Checking notebook inputs:"
for f in GSE135251_samples.csv cdna/t2g.csv rnaseq_data/kallisto_output/abundance.tsv; do
  [ -s "$f" ] && echo "  OK  $f" || echo "  MISSING  $f"
done
echo "  GSE135251/*.counts.txt.gz : $(ls GSE135251/*.counts.txt.gz 2>/dev/null | wc -l | tr -d ' ') files"
