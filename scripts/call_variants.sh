sample=$1;

ref_genome_path=ref_genome/SARS-CoV-2.fa

final_bam_out=bam_alignments/"$sample"/"$sample".bam
final_strict_bam_out=bam_alignments/"$sample"/"$sample".strict.bam


bcftools mpileup --min-BQ 20 --max-depth 50000 --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR -f "$ref_genome_path" "$final_bam_out" | bcftools call --ploidy 1 -Ov -m -A -M | bcftools view --include 'INFO/DP>=5' -Oz -o variant_calls/standard/"$sample".vcf.gz

bcftools mpileup --min-BQ 30 --max-depth 50000 --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR -f "$ref_genome_path" "$final_strict_bam_out" | bcftools call --ploidy 1 -Ov -m -A -M | bcftools view --include 'INFO/DP>=5' -Oz -o variant_calls/strict/"$sample".strict.vcf.gz
