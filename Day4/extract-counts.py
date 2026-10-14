"""
Merges individual HTSeq-count .txt files into a single tab-delimited 
count matrix with genes as rows and samples as columns.
"""

import argparse
import os
import sys
from collections import defaultdict


def parse_arguments():
    """Parse and validate command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Merge HTSeq-count output files into a single expression matrix."
    )
    parser.add_argument(
        "infolder",
        type=str,
        help="Path to folder containing HTSeq-count .txt files"
    )
    parser.add_argument(
        "outfile",
        type=str,
        help="Path for the output merged tab-delimited file"
    )

    # Display full help message if no arguments are provided
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_args()


def get_count_files(folder_path):
    """Return a sorted list of .txt filenames, excluding hidden files."""
    if not os.path.isdir(folder_path):
        sys.exit(f"Error: Directory not found: {folder_path}")

    files = sorted(
        f for f in os.listdir(folder_path)
        if f.endswith(".txt") and not f.startswith(".")
    )

    if not files:
        sys.exit(f"Error: No .txt count files found in {folder_path}")

    return files


def load_counts(folder_path, filenames):
    counts = defaultdict(dict)
    samples = []

    for fname in filenames:
        # Strip xtensions to create a clean column header
        sample_name = fname.replace("_htseqcount.txt", "").replace(".txt", "")
        samples.append(sample_name)

        file_path = os.path.join(folder_path, fname)
        with open(file_path, "r") as fh:
            for line in fh:
                line = line.strip()
                # Skip blank lines and HTSeq summary statistics (e.g., __no_feature)
                if not line or line.startswith("__"):
                    continue

                parts = line.split("\t")
                if len(parts) >= 2:
                    gene_id, count_val = parts[0], float(parts[1])
                    counts[gene_id][sample_name] = count_val

    return samples, counts


def write_matrix(output_path, samples, counts):
    """Write the nested count dictionary to a tab-delimited matrix file."""
    with open(output_path, "w") as out:
        # Write column header
        out.write("Geneid\t" + "\t".join(samples) + "\n")
        # Write gene rows sorted alphabetically
        for gene in sorted(counts.keys()):
            # Default to 0.0 if a gene was absent in a particular sample
            row_vals = [str(counts[gene].get(s, 0.0)) for s in samples]
            out.write(f"{gene}\t" + "\t".join(row_vals) + "\n")

def main():
    args = parse_arguments()

    print(f"Scanning directory: {args.infolder}")
    filenames = get_count_files(args.infolder)
    print(f"Found {len(filenames)} count file(s). Extracting counts...")

    samples, counts = load_counts(args.infolder, filenames)
    print(f"Extracted {len(counts)} unique genes across {len(samples)} samples.")

    print(f"Writing merged matrix to: {args.outfile}")
    write_matrix(args.outfile, samples, counts)
    print("Done.")


if __name__ == "__main__":
    main()