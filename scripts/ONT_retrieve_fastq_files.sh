sample=$1;

outdir=ONT_fastq_files/"$sample"

lsout=$(ls -A ONT_fastq_files/"$sample")

if [ -z "$lsout" ]; then
    mkdir -p "$outdir"
    fasterq-dump -e 2 -t sra_tmpdir --outdir "$outdir" "$sample"
fi
