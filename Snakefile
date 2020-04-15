samples = []
for line in shell("cat current_illumina.txt", iterable=True):
    samples.append(line)


rule all:
    input:
        ancient(expand("variant_calls/strict/{sample}.strict.vcf", sample=samples))


rule retrieve_and_qc_fastq_files:
    """
    Filter reads using fastp.
    Two sets of filtered reads are produced, one stricter.
    """
    input:
        lambda wildcards: ancient(expand("fastq_files/{sample}", sample=wildcards.sample))
    output:
        protected("qc_fastq_files/{sample}/{sample}.log")
    log:
        "logs/retrieve_and_qc_fastq_files/{sample}.log"
    threads:
        2
    shell:
        """
        sh scripts/retrieve_fastq_files.sh {wildcards.sample}
        sh scripts/qc_fastq_files.sh {wildcards.sample}
        """

rule align_reads:
    """
    Align reads using bwa mem.
    """
    input:
        lambda wildcards: ancient(expand("qc_fastq_files/{sample}/{sample}.log", sample=wildcards.sample)),
        lambda wildcards: ancient(expand("qc_fastq_files/{sample}/{sample}.log", sample=wildcards.sample))
    output:
        protected("bam_alignments/{sample}/{sample}.bam"),
        protected("bam_alignments/{sample}/{sample}.strict.bam")
    log:
        "logs/align_reads/{sample}.log"
    shell:
        "sh scripts/align_reads.sh {wildcards.sample}"


rule call_variants:
    """
    Call variants using bcftools mpileup | bcftools call | bcftools view.
    """
    input:
        lambda wildcards: ancient(expand("bam_alignments/{sample}/{sample}.bam", sample=wildcards.sample)),
        lambda wildcards: ancient(expand("bam_alignments/{sample}/{sample}.strict.bam", sample=wildcards.sample)),
    output:
        protected("variant_calls/standard/{sample}.vcf"),
        protected("variant_calls/strict/{sample}.strict.vcf")
    log:
        "logs/call_variants/{sample}.log"
    shell:
        "sh scripts/call_variants.sh {wildcards.sample}"
