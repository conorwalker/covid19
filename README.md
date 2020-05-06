To run the pipeline:
sh run_pipeline.sh

Two sets of QC'd FASTQ, bam, vcf and cf files are produced, according to two sets of filtering criteria:

Standard filtering:
- Base quality >=30 at >=90% of the read
- Base quality >=30 at 3' end (mean quality in window, default window size in fastp is 4)
- Read length >=50
- Minimum base quality of 20 at site to report variant
- Minimum coverage of 5 mapped reads per site to report any variants

Strict filtering:
- Base quality >=35 at >=90% of the read
- Base quality >=30 at 3' end (mean quality in window, default window size in fastp is 4)
- 15 bases are trimmed from the start and end of each read
- Read length >= 50 (after trimming)
- Soft or hard clipped reads removed from the alignment
- Minimum base quality of 30 at site to report variant
- Minimum coverage of 5 mapped reads per site to report any variants
