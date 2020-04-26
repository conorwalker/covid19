samples = []
for line in shell("cat current_illumina.txt", iterable=True):
    samples.append(line)


nanopore_samples = []
for line in shell("cat current_ONT.txt", iterable=True):
    samples.append(line)


variant_callsets = ["alt_only_variant_calls", "lofreq_variant_calls", "variant_calls"]
filtering = ["standard", "strict"]


rule all:
    input:
        ancient(expand("{varcalls}/illumina/{filterset}/{sample}.strict.vcf.gz", sample=samples, filterset=filtering, varcalls=variant_callsets),
        ancient(expand("{varcalls}/ONT/{sample}.strict.vcf.gz", sample=samples, varcalls=variant_callsets)
        )


rule retrieve_and_qc_fastq_files:
    """
    Filter reads using fastp.
    Two sets of filtered reads are produced, one stricter.
    """
    input:
        lambda wildcards: ancient(expand("fastq_files/illumina/{sample}", sample=wildcards.sample))
    output:
        "qc_fastq_files/illumina/{sample}/{sample}.log"
    log:
        "logs/retrieve_and_qc_fastq_files/illumina/{sample}.log"
    threads:
        2
    shell:
        """
        sh scripts/retrieve_fastq_files.sh {wildcards.sample} illumina
        sh scripts/qc_fastq_files.sh {wildcards.sample}
        """


rule ONT_retrieve_and_qc_fastq_files:
    """
    Filter reads using fastp.
    Two sets of filtered reads are produced, one stricter.
    """
    input:
        lambda wildcards: ancient(expand("fastq_files/ONT/{sample}", sample=wildcards.sample))
    output:
        "qc_fastq_files/ONT/{sample}/{sample}.log"
    log:
        "logs/retrieve_and_qc_fastq_files/ONT/{sample}.log"
    threads:
        2
    shell:
        """
        sh scripts/ONT_retrieve_fastq_files.sh {wildcards.sample} ONT
        sh scripts/ONT_qc_fastq_files.sh {wildcards.sample}
        """


rule align_reads:
    """
    Align reads using bwa mem.
    """
    input:
        lambda wildcards: ancient(expand("qc_fastq_files/{sample}/{sample}.log", sample=wildcards.sample)),
        lambda wildcards: ancient(expand("qc_fastq_files/{sample}/{sample}.log", sample=wildcards.sample))
    output:
        "bam_alignments/{sample}/{sample}.bam",
        "bam_alignments/{sample}/{sample}.strict.bam"
    log:
        "logs/align_reads/{sample}.log"
    threads:
        2
    shell:
        "sh scripts/align_reads.sh {wildcards.sample}"


rule ONT_align_reads:
    """
    Align reads using bwa mem.
    """
    input:
        lambda wildcards: ancient(expand("qc_fastq_files/{sample}/{sample}.log", sample=wildcards.sample)),
        lambda wildcards: ancient(expand("qc_fastq_files/{sample}/{sample}.log", sample=wildcards.sample))
    output:
        "bam_alignments/{sample}/{sample}.bam",
        "bam_alignments/{sample}/{sample}.strict.bam"
    log:
        "logs/align_reads/{sample}.log"
    threads:
        2
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
        "variant_calls/standard/{sample}.vcf.gz",
        "variant_calls/strict/{sample}.strict.vcf.gz"
    log:
        "logs/call_variants/{sample}.log"
    threads:
        2
    shell:
        "sh scripts/call_variants.sh {wildcards.sample}"


rule subset_alt_only_variants:
    """
    Subset variant call files to output vcfs containing only alternate genotype positions.
    """
    input:
        lambda wildcards: ancient(expand("variant_calls/standard/{sample}.vcf.gz", sample=wildcards.sample)),
        lambda wildcards: ancient(expand("variant_calls/strict/{sample}.strict.vcf.gz", sample=wildcards.sample))
    output:
        "alt_only_variant_calls/standard/{sample}.vcf.gz",
        "alt_only_variant_calls/strict/{sample}.strict.vcf.gz"
    log:
        "logs/alt_only_variant_calls/{sample}.log"
    shell:
        "sh scripts/subset_variant_sites.sh {wildcards.sample}"    


rule call_lofreq_variants:
    """
    Call variants using bcftools mpileup | bcftools call | bcftools view.
    """
    input:
        lambda wildcards: ancient(expand("bam_alignments/{sample}/{sample}.bam", sample=wildcards.sample)),
        lambda wildcards: ancient(expand("bam_alignments/{sample}/{sample}.strict.bam", sample=wildcards.sample))
    output:
        "lofreq_variant_calls/standard/{sample}.vcf.gz",
        "lofreq_variant_calls/strict/{sample}.strict.vcf.gz"
    log:
        "logs/lofreq_call_variants/{sample}.log"
    shell:
        "sh scripts/lofreq_call_variants.sh {wildcards.sample}"
