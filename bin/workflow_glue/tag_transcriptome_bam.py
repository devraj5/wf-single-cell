#!/usr/bin/env python
import argparse
import pysam
from collections import namedtuple
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

TagData = namedtuple('TagData', ['cb', 'cr', 'cy', 'ub', 'ur', 'uy', 'gn', 'tr'])

def read_tag_file(filename):
    read_tags = {}
    with open(filename, "r") as f:
        next(f)  # Skip header
        for line in f:
            cols = line.strip().split("\t")
            read_id = cols[0]
            read_tags[read_id] = TagData(*cols[2:10])
    return read_tags

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_bam', help='Input BAM file')
    parser.add_argument('output_bam', help='Output BAM file')
    parser.add_argument('reads_summary', help='TSV file with read tags')
    parser.add_argument('--threads', type=int, default=1, help='Number of threads to use')
    args = parser.parse_args()

    read_tags = read_tag_file(args.reads_summary)
    logger.info(f"Loaded {len(read_tags)} read tags")

    with pysam.AlignmentFile(args.input_bam, "rb", threads=args.threads) as bamfile, \
         pysam.AlignmentFile(args.output_bam, "wb", template=bamfile, threads=args.threads) as output_bam:
        
        total_reads = 0
        tagged_reads = 0

        for ref in bamfile.references:
            logger.info(f"Processing reference: {ref}")
            for read in bamfile.fetch(ref):
                total_reads += 1
                read_id = read.query_name
                if read_id in read_tags:
                    tagged_reads += 1
                    tags = read_tags[read_id]
                    read.set_tag("CB", tags.cb)
                    read.set_tag("CR", tags.cr)
                    read.set_tag("CY", tags.cy)
                    read.set_tag("UB", tags.ub)
                    read.set_tag("UR", tags.ur)
                    read.set_tag("UY", tags.uy)
                    read.set_tag("GN", tags.gn)
                    read.set_tag("TR", tags.tr)
                output_bam.write(read)

        logger.info(f"Processed {total_reads} reads")
        logger.info(f"Tagged {tagged_reads} reads ({tagged_reads/total_reads:.2%})")

if __name__ == "__main__":
    main()
