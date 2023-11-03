#!/usr/bin/env python
import argparse
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description='Rename reads based on barcodes')
    parser.add_argument('-f', '--fastq', help='Input FASTQ file after removing barcodes', required=True)
    parser.add_argument('-b', '--barcode', help='Barcodes file', required=True)
    parser.add_argument('-m', '--mapping', help='Mapping file', required=True)
    parser.add_argument('-o', '--output', help='Output directory', required=True)
    return parser.parse_args()

def read_mapping(mapping_file):
    mapping = {}
    with open(mapping_file, 'r') as f:
        for line in f:
            cols = line.strip().split('\t')
            barcode = cols[1]
            mapping[barcode] = line.strip()
    return mapping

def rename_reads(fastq_file, barcode_file, mapping, output_dir):
    output_mapping_file = output_dir + '/mapping_dealed.txt'
    output_fastq_file = output_dir + '/reads_renamed.fastq'

    renamed_reads = []
    with open(barcode_file, 'r') as b, open(output_mapping_file, 'w') as out_mapping:
        for read, barcode in zip(SeqIO.parse(fastq_file, 'fastq'), b):
            barcode = barcode.strip()
            if barcode in mapping:
                cols = mapping[barcode].split('\t')
                cols.append(barcode)
                out_mapping.write('\t'.join(cols) + '\n')
                read.id = cols[0]  # Replace read ID
                read.description = ''  # Remove description if any
                renamed_reads.append(read)

    SeqIO.write(renamed_reads, output_fastq_file, 'fastq')

def main():
    args = parse_args()
    mapping = read_mapping(args.mapping)
    rename_reads(args.fastq, args.barcode, mapping, args.output)

if __name__ == '__main__':
    main()

