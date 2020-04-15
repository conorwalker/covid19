sample=$1;

outdir=fastq_files/"$sample"

lsout=$(ls -A fastq_files/"$sample")

if [ -z "$lsout" ]; then
    mkdir -p "$outdir"
    fasterq-dump -e 2 -t sra_tmpdir --outdir "$outdir" "$sample"
fi
