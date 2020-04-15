#!/usr/local/bin/python

from glob import glob
from sys import argv


def get_reference_seq_len(ref_path):
    with open(ref_path, "r") as f:
        lines = [i.strip() for i in f.readlines()]
    seqlen = 0
    for l in lines:
        if l.startswith(">"):
            pass
        else:
            seqlen += len(l)
    return seqlen


def main():
    # get number of positions in reference genome
    ref_genome_path = "ref_genome/SARS-CoV-2.fa"
    ref_len = get_reference_seq_len(ref_genome_path)

    # get list of vcf files
    vcf_list = glob("{0}/*vcf".format(argv[1]))
    
    sample_names = []
    sample_positions = {}

    for vcf in vcf_list:
        sample = vcf.split("/")[-1].split(".")[0]
        sample_names.append(sample)
        sample_positions[sample] = ["0,0,0,0"] * ref_len
        
        with open(vcf, "r") as vcf_fi:
            vcf_lines = [i.strip() if i[0] != "#" else "" for i in vcf_fi.readlines()]
        
        # get allele counts for each line of each VCF file, store in the
        # sample_positions dictionary
        for line in vcf_lines:
            # if line started with a comment
            if not line:
                continue
            line = line.split()
            seq_position = int(line[1])
            ref_allele = line[3]
            alt_allele = line[4]
            info_field = line[7]
            info_dic = {}
            split_info = info_field.split(";")
            # do not process line if the INDEL flag is present,
            # the surrounding non-indel rows should give base counts for this
            # region
            if split_info[0] == "INDEL":
                continue

            for val in split_info:
                split_val = val.split("=")
                info_dic[split_val[0]] = split_val[1]

            allelic_depths = line[-1].split(":")[-1].split(",")
            counts = {"A":0, "C":0, "G":0, "T":0}
            # count ref allele
            counts[ref_allele] += int(allelic_depths[0])
            total_alt_alleles = 0
            split_alt_alleles = alt_allele.split(",")

            if alt_allele != ".":
                total_alt_alleles = len(split_alt_alleles)
            # if no alt allele present
            if total_alt_alleles == 0:
                out_val = ",".join([str(i) for i in [counts["A"],counts["C"],counts["G"], counts["T"]]])
                sample_positions[sample][seq_position-1] = out_val
                continue
            # if alt allele(s) present, process allele counts
            for i in range(total_alt_alleles):
                try:
                    counts[split_alt_alleles[i]] = allelic_depths[i+1]
                except IndexError:
                    counts[split_alt_alleles[i]] = 0
            out_val = ",".join([str(i) for i in [counts["A"],counts["C"],counts["G"], counts["T"]]])
            sample_positions[sample][seq_position-1] = out_val
    
    # check number of positions in each sample matches reference length
    # should always be true
    for val in sample_positions.values():
        assert len(val) == ref_len

    
    # write output file for PoMo
    header1 = "COUNTSFILE\tNPOP {0}\tNSITES {1}".format(len(vcf_list), ref_len)
    header2 = "CHROM\tPOS\t" + "\t".join(sample_names)
    with open("{0}".format(argv[2]), "w") as outfi:
        outfi.write(header1 + "\n")
        outfi.write(header2 + "\n")
        for i in range(ref_len):
            line_out_vals = []
            for sv in sample_names:
                line_out_vals.append(sample_positions[sv][i])
            # add POS
            line_out_vals.insert(0, str(i+1))
            # add CHROM placeholder
            line_out_vals.insert(0, "42")
            outfi.write("\t".join(line_out_vals) + "\n")


if __name__ == "__main__":
    main()
