sample=$1;
technology=$2;

outdir=fastq_files/"$technology"/"$sample"

lsout=$(ls -A fastq_files/"$technology"/"$sample")

if [ -z "$lsout" ]; then
    mkdir -p "$outdir"
    fasterq-dump -e 2 -t sra_tmpdir --outdir "$outdir" "$sample"
fi
