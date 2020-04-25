sample=$1;

ref_genome_path=ref_genome/SARS-CoV-2.fa
picard_metrics=/tmp/"$sample"_picard_metrics.txt

tmp_bam_out=bam_alignments/"$sample"/"$sample".qc.tmp.bam
tmp_bam_out_paired=bam_alignments/"$sample"/"$sample".qc.tmp.paired.bam
final_bam_out=bam_alignments/"$sample"/"$sample".bam

tmp_strict_bam_out=bam_alignments/"$sample"/"$sample".qc.tmp.strict.bam
tmp_strict_bam_out_paired=bam_alignments/"$sample"/"$sample".qc.tmp.strict.paired.bam
final_strict_bam_out=bam_alignments/"$sample"/"$sample".strict.bam

qc_log_fi=qc_fastq_files/"$sample"/"$sample".log

if [[ $(cat $qc_log_fi) == "fin" ]]; then
    # paired and unpaired reads
    if [[ -f qc_fastq_files/"$sample"/"$sample".qc.fastq && -f qc_fastq_files/"$sample"/"$sample"_1.qc.fastq ]]; then
        unpaired_fastq_qc=qc_fastq_files/"$sample"/"$sample".qc.fastq
        read1_fastq_qc=qc_fastq_files/"$sample"/"$sample"_1.qc.fastq
        read2_fastq_qc=qc_fastq_files/"$sample"/"$sample"_2.qc.fastq
        bwa mem -t 2 "$ref_genome_path" "$unpaired_fastq_qc" | samtools sort -o "$tmp_bam_out"
        bwa mem -t 2 "$ref_genome_path" "$read1_fastq_qc" "$read2_fastq_qc" | samtools sort -o "$tmp_bam_out_paired"
        picard MarkDuplicates INPUT="$tmp_bam_out" INPUT="$tmp_bam_out_paired" output="$final_bam_out" METRICS_FILE="$picard_metrics" REMOVE_DUPLICATES=true TMP_DIR=/tmp
        rm $tmp_bam_out
        rm $tmp_bam_out_paired

        strict_unpaired_fastq_qc=qc_fastq_files/"$sample"/"$sample".qc.strict.fastq
        strict_read1_fastq_qc=qc_fastq_files/"$sample"/"$sample"_1.qc.strict.fastq
        strict_read2_fastq_qc=qc_fastq_files/"$sample"/"$sample"_2.qc.strict.fastq
        bwa mem -t 2 "$ref_genome_path" "$strict_unpaired_fastq_qc" | awk '$6 !~ /H|S/' | samtools sort -o "$tmp_strict_bam_out"
        bwa mem -t 2 "$ref_genome_path" "$strict_read1_fastq_qc" "$strict_read2_fastq_qc" | awk '$6 !~ /H|S/' | samtools sort -o "$tmp_strict_bam_out_paired"
        picard MarkDuplicates INPUT="$tmp_strict_bam_out" INPUT="$tmp_strict_bam_out_paired" output="$final_strict_bam_out" METRICS_FILE="$picard_metrics" REMOVE_DUPLICATES=true TMP_DIR=/tmp
        rm $tmp_strict_bam_out
        rm $tmp_strict_bam_out_paired
    fi


    # paired reads
    if [[ ! -f qc_fastq_files/"$sample"/"$sample".qc.fastq && -f qc_fastq_files/"$sample"/"$sample"_1.qc.fastq ]]; then
        read1_fastq_qc=qc_fastq_files/"$sample"/"$sample"_1.qc.fastq
        read2_fastq_qc=qc_fastq_files/"$sample"/"$sample"_2.qc.fastq
        bwa mem -t 2 "$ref_genome_path" "$read1_fastq_qc" "$read2_fastq_qc" | samtools sort -o "$tmp_bam_out_paired"
        picard MarkDuplicates INPUT="$tmp_bam_out_paired" output="$final_bam_out" METRICS_FILE="$picard_metrics" REMOVE_DUPLICATES=true TMP_DIR=/tmp
        rm $tmp_bam_out_paired
        
        strict_read1_fastq_qc=qc_fastq_files/"$sample"/"$sample"_1.qc.strict.fastq
        strict_read2_fastq_qc=qc_fastq_files/"$sample"/"$sample"_2.qc.strict.fastq
        bwa mem -t 2 "$ref_genome_path" "$strict_read1_fastq_qc" "$strict_read2_fastq_qc" | awk '$6 !~ /H|S/' | samtools sort -o "$tmp_strict_bam_out_paired"
        picard MarkDuplicates INPUT="$tmp_strict_bam_out_paired" output="$final_strict_bam_out" METRICS_FILE="$picard_metrics" REMOVE_DUPLICATES=true TMP_DIR=/tmp
        rm $tmp_strict_bam_out_paired
    fi

    # unpaired reads
    if [[ -f qc_fastq_files/"$sample"/"$sample".qc.fastq && ! -f qc_fastq_files/"$sample"/"$sample"_1.qc.fastq ]]; then
        unpaired_fastq_qc=qc_fastq_files/"$sample"/"$sample".qc.fastq
        bwa mem -t 2 "$ref_genome_path" "$unpaired_fastq_qc" | samtools sort -o "$tmp_bam_out"
        picard MarkDuplicates INPUT="$tmp_bam_out" output="$final_bam_out" METRICS_FILE="$picard_metrics" REMOVE_DUPLICATES=true TMP_DIR=/tmp
        rm $tmp_bam_out

        strict_unpaired_fastq_qc=qc_fastq_files/"$sample"/"$sample".qc.strict.fastq
        bwa mem -t 2 "$ref_genome_path" "$strict_unpaired_fastq_qc" | awk '$6 !~ /H|S/' | samtools sort -o "$tmp_strict_bam_out"
        picard MarkDuplicates INPUT="$tmp_strict_bam_out" output="$final_strict_bam_out" METRICS_FILE="$picard_metrics" REMOVE_DUPLICATES=true TMP_DIR=/tmp
        rm $tmp_strict_bam_out
    fi
else
    echo "" > "$final_bam_out"
    echo "" > "$final_strict_bam_out"
fi
