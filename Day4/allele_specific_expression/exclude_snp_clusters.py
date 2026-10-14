"""
Exclude SNPs that occur in dense clusters from a VCF.

A SNP is excluded if it belongs to a window of length `read_length`
containing more than `mismatches` SNPs. The output VCF retains all
non-clustered variants.
"""

import argparse
import os
import sys
from collections import defaultdict

import vcfpy

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Remove SNPs found in clusters from a VCF file.")
    parser.add_argument("input",help="Input VCF file.",)
    parser.add_argument("-l", "--read_length",type=int,help="RNA-seq read length in bases.",)
    parser.add_argument("-m", "--mismatches",type=int,help="Maximum number of SNPs allowed within one read-length window.",)
    return parser.parse_args()

def get_snp_positions(vcf_reader):
    """Collect SNP genomic positions, grouped by chromosome."""
    positions_by_chrom = defaultdict(list)
    for record in vcf_reader:
        positions_by_chrom[record.CHROM].append(record.POS)
    return positions_by_chrom

def find_clustered_snps(snp_positions, read_length, max_mismatches):
    """Identify SNPs that occur in clusters."""
    excluded = set()
    for i, start_pos in enumerate(snp_positions):
        cluster = [start_pos]

        for pos in snp_positions[i + 1:]:
            if pos > start_pos + read_length:
                break
            cluster.append(pos)

        # More SNPs than allowed mismatches means this is a problematic cluster.
        if len(cluster) > max_mismatches:
            excluded.update(cluster)
    return excluded

def get_clustered_snps(snp_positions_by_chrom, read_length, max_mismatches):
    """Find clustered SNPs independently for each chromosome."""
    return {
        chrom: find_clustered_snps(positions, read_length, max_mismatches)
        for chrom, positions in snp_positions_by_chrom.items()
    }

def filter_clustered_snps(infile, outfile, read_length, max_mismatches):
    """Read a VCF, identify clustered SNPs, and write a filtered VCF."""
    # First pass: collect SNP positions and identify clusters.
    with vcfpy.Reader.from_path(infile) as reader:
        snp_positions = get_snp_positions(reader)

    excluded_snps = get_clustered_snps(
        snp_positions, read_length, max_mismatches
    )

    # Second pass: copy the VCF while omitting clustered SNPs.
    with vcfpy.Reader.from_path(infile) as reader, \
         vcfpy.Writer.from_path(outfile, reader.header) as writer:

        for record in reader:
            if record.POS not in excluded_snps[record.CHROM]:
                writer.write_record(record)

def make_output_path(infile):
    """Create the output VCF path."""
    directory, filename = os.path.split(infile)
    base_name = filename.split(".")[0]
    return os.path.join(directory, f"{base_name}_noclusters.vcf")

def main():
    args = parse_args()
    outfile = make_output_path(args.input)

    print(f"Filtering file: {args.input}")
    print(f"Writing filtered VCF to: {outfile}")

    filter_clustered_snps(
        infile=args.input,
        outfile=outfile,
        read_length=args.read_length,
        max_mismatches=args.mismatches,
    )


if __name__ == "__main__":
    main()
