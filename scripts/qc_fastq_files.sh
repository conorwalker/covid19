sample=$1;
unpaired_in=fastq_files/"$sample"/"$sample".fastq
read1_in=fastq_files/"$sample"/"$sample"_1.fastq
read2_in=fastq_files/"$sample"/"$sample"_2.fastq

unpaired_out=qc_fastq_files/"$sample"/"$sample".qc.fastq
read1_out=qc_fastq_files/"$sample"/"$sample"_1.qc.fastq
read2_out=qc_fastq_files/"$sample"/"$sample"_2.qc.fastq

strict_unpaired_out=qc_fastq_files/"$sample"/"$sample".qc.strict.fastq
strict_read1_out=qc_fastq_files/"$sample"/"$sample"_1.qc.strict.fastq
strict_read2_out=qc_fastq_files/"$sample"/"$sample"_2.qc.strict.fastq


# check percentage of N bases in reads
sample_n_bases=$(python /nfs/research1/goldman/conor/covid19/within_host_variation/process_sra_runs/scripts/calculate_n_in_fastq.py "$sample" 2)


if [ "$sample_n_bases" == "acceptable" ]; then
    # paired and unpaired reads
    if [[ -f fastq_files/"$sample"/"$sample".fastq && -f fastq_files/"$sample"/"$sample"_1.fastq ]]; then
        # qc
        fastp -q 30 -u 10 -l 50 -x --cut_tail --cut_tail_mean_quality 30 -i "$unpaired_in" -o "$unpaired_out"
        fastp -q 30 -u 10 -l 50 -x --cut_tail --cut_tail_mean_quality 30 -i "$read1_in" -I "$read2_in" -o "$read1_out" -O "$read2_out"
        # strict qc
        fastp -q 35 -u 10 -l 50 -x --cut_tail --cut_tail_mean_quality 30 --trim_front1 15 --trim_tail1 15 -i "$unpaired_in" -o "$strict_unpaired_out"
        fastp -q 35 -u 10 -l 50 -x --cut_tail --cut_tail_mean_quality 30 --trim_front1 15 --trim_tail1 15 -i "$read1_in" -I "$read2_in" -o "$strict_read1_out" -O "$strict_read2_out"
    fi

    # paired reads
    if [[ ! -f fastq_files/"$sample"/"$sample".fastq && -f fastq_files/"$sample"/"$sample"_1.fastq ]]; then
        # qc
        fastp -q 30 -u 10 -l 50 -x --cut_tail --cut_tail_mean_quality 30 -i "$read1_in" -I "$read2_in" -o "$read1_out" -O "$read2_out"
        # strict qc
        fastp -q 35 -u 10 -l 50 -x --cut_tail --cut_tail_mean_quality 30 --trim_front1 15 --trim_tail1 15 -i "$read1_in" -I "$read2_in" -o "$strict_read1_out" -O "$strict_read2_out"
    fi

    # unpaired reads
    if [[ -f fastq_files/"$sample"/"$sample".fastq && ! -f fastq_files/"$sample"/"$sample"_1.fastq ]]; then
        # qc
        fastp -q 30 -u 10 -l 50 -x --cut_tail --cut_tail_mean_quality 30 -i "$unpaired_in" -o "$unpaired_out"
        # strict qc
        fastp -q 35 -u 10 -l 50 -x --cut_tail --cut_tail_mean_quality 30 --trim_front1 15 --trim_tail1 15 -i "$unpaired_in" -o "$strict_unpaired_out"
    fi

    echo "fin" > qc_fastq_files/"$sample"/"$sample".log
else
    echo "excluded" > qc_fastq_files/"$sample"/"$sample".log
fi
