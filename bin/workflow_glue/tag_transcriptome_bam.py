#!/usr/bin/env python
"""Tag transcriptome-aligned BAM file with information from read_summary.tsv"""

import pysam
import csv
import argparse
from .util import get_named_logger, wf_parser  # noqa: ABS101

def main():
    parser = wf_parser("tag_transcriptome_bam")
    parser.add_argument("input_bam", help="Input transcriptome-aligned BAM file")
    parser.add_argument("output_bam", help="Output tagged BAM file")
    parser.add_argument("read_summary", help="read_summary.tsv file")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads to use")
    args = parser.parse_args()

    logger = get_named_logger("TagTranscriptomeBAM")
    
    read_info = load_read_summary(args.read_summary)
    tag_bam(args.input_bam, args.output_bam, read_info, args.threads, logger)

def load_read_summary(file_path):
    read_info = {}
    with open(file_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            read_info[row['read_id']] = row
    return read_info

def tag_bam(input_bam, output_bam, read_info, threads, logger):
    with pysam.AlignmentFile(input_bam, "rb", threads=threads) as in_bam, \
         pysam.AlignmentFile(output_bam, "wb", template=in_bam, threads=threads) as out_bam:
        
        total_reads = 0
        tagged_reads = 0
        
        for read in in_bam:
            total_reads += 1
            read_id = read.query_name
            if read_id in read_info:
                tagged_reads += 1
                info = read_info[read_id]
                read.set_tag("CB", info['corrected_barcode'], "Z")
                read.set_tag("CR", info['uncorrected_barcode'], "Z")
                read.set_tag("CY", info['quality_barcode'], "Z")
                read.set_tag("UB", info['corrected_umi'], "Z")
                read.set_tag("UR", info['uncorrected_umi'], "Z")
                read.set_tag("UY", info['quality_umi'], "Z")
                read.set_tag("GN", info['gene'], "Z")
                read.set_tag("TR", info['transcript'], "Z")
                
                # Add transcriptome-specific tags
                read.set_tag("TS", read.reference_name, "Z")  # Transcript name
                read.set_tag("TL", read.reference_length, "i")  # Transcript length
            
            out_bam.write(read)

    logger.info(f"Total reads processed: {total_reads}")
    logger.info(f"Reads tagged: {tagged_reads} ({tagged_reads/total_reads:.2%})")

if __name__ == "__main__":
    main()
