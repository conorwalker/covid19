To run the pipeline:
sh run_pipeline.sh

This will retrieve the latest "current_illumina.txt" and "current_metadata.txt" files from https://raw.githubusercontent.com/galaxyproject/SARS-CoV-2/master/genomics/4-Variation/

Two sets of QC'd FASTQ, bam, vcf and cf files are produced, according to two sets of filtering criteria:

Standard filtering:
- Base quality >=30 at >=90% of the read
- Base quality >=30 at 3' end (mean quality in window, default window size in fastp is 4)
- Read length >=50
- Minimum coverage of 5 mapped reads per site to report any variants

Strict filtering:
- Base quality >=35 at >=90% of the read
- Base quality >=30 at 3' end (mean quality in window, default window size in fastp is 4)
- 15 bases are trimmed from the start and end of each read
- Read length >= 50 (after trimming)
- Soft or hard clipped reads removed from the alignment
- Minimum coverage of 5 mapped reads per site to report any variants
