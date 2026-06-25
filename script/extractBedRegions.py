#!/usr/bin/env python3
"""
extractBedRegions.py
Extract FASTA subsequences specified by a BED file.
Requires: biopython
"""

import argparse
import sys
from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser(description="Extract FASTA regions from BED")
    parser.add_argument("--bed", required=True, help="Input BED file")
    parser.add_argument("--fasta", required=True, help="Input FASTA file")
    parser.add_argument("--out-prefix", required=True, help="Output prefix (dir + basename)")
    args = parser.parse_args()

    # Load all sequences into memory
    seq_dict = {}
    for rec in SeqIO.parse(args.fasta, "fasta"):
        seq_dict[rec.id] = str(rec.seq)

    with open(args.bed) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue

            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            if chrom not in seq_dict:
                print(f"[extractBed] Warning: {chrom} not found in FASTA, skipping", file=sys.stderr)
                continue

            seq = seq_dict[chrom]
            subseq = seq[start:end]
            if not subseq:
                continue

            # Name: chrom_start_end
            region_name = f"{chrom}_{start}_{end}"
            out_path = f"{args.out_prefix}_{region_name}.fa"
            with open(out_path, "w") as out:
                out.write(f">{region_name}\n{subseq}\n")
            print(f"[extractBed] {region_name} ({len(subseq)} bp) -> {out_path}")


if __name__ == "__main__":
    main()
