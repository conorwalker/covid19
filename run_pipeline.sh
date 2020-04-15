# make some directories necessary for pipeline
if [ ! -d backup_accessions ]; then
    mkdir backup_accessions;
fi
if [ ! -d sra_tmpdir ]; then
    mkdir sra_tmpdir;
fi


# retrieve accession IDs from illumina runs
mv current_illumina.txt backup_accessions/
mv current_metadata.txt backup_accessions/
wget https://raw.githubusercontent.com/galaxyproject/SARS-CoV-2/master/genomics/4-Variation/current_illumina.txt
wget https://raw.githubusercontent.com/galaxyproject/SARS-CoV-2/master/genomics/4-Variation/current_metadata.txt


# fetch reference genome and index if it doesn't exist locally
ref_genome_path=ref_genome/SARS-CoV-2.fa
if [ ! -f "$ref_genome_path" ]; then
    mkdir ref_genome
    efetch -db nucleotide -id MN908947.3 -format fasta > "$ref_genome_path"
    # index the genome for downstream analysis
    bwa index "$ref_genome_path"
    samtools faidx "$ref_genome_path"
    echo
fi

# check if all fastq files up-to-date
while read sample; do
    outdir=fastq_files/"$sample"
    if [ ! -d "$outdir" ]; then
        mkdir "$outdir"
    fi
done < current_illumina.txt

# run snakemake pipeline
#snakemake all --nocolor --cores 2
snakemake all --latency-wait 120 --jobs 500 --cluster 'bsub -e /dev/null -o /dev/null -M 10000 -R "rusage[mem=10000]" -n 2' -k --nocolor

# clean up
rm fastp.html
rm fastp.json

# produce count files for PoMo
if [ ! -d count_files ]; then
    mkdir count_files
fi
python scripts/vcf_to_cf.py variant_calls/standard count_files/standard.cf 
python scripts/vcf_to_cf.py variant_calls/strict count_files/strict.cf

tar cvf - count_files/standard.cf | gzip -9 > standard.cf.tar.gz
tar cvf - count_files/strict.cf | gzip -9 > strict.cf.tar.gz
