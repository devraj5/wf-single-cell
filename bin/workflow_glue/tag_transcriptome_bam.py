"""Tag transcriptome-aligned BAM file with information from read_summary.tsv"""

import pysam
import csv
from .util import get_named_logger, wf_parser

def argparser():
    parser = wf_parser("tag_transcriptome_bam")
    parser.add_argument("input_bam", help="Input transcriptome-aligned BAM file")
    parser.add_argument("output_bam", help="Output tagged BAM file")
    parser.add_argument("read_summary", help="read_summary.tsv file")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads to use")
    return parser

def main(args):
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
        skipped_reads = 0
        
        for read in in_bam:
            total_reads += 1
            read_id = read.query_name
            if read_id in read_info:
                info = read_info[read_id]
                
                # Checking if barcode and UMI tags are present
                if (info['corrected_barcode'] != '-' and info['corrected_umi'] != '-'):
                    tagged_reads += 1
                    
                    # Settingtags only if they are present in the read_info
                    if info['corrected_barcode'] != '-':
                        read.set_tag("CB", info['corrected_barcode'], "Z")
                    if info['uncorrected_barcode'] != '-':
                        read.set_tag("CR", info['uncorrected_barcode'], "Z")
                    if info['quality_barcode'] != '-':
                        read.set_tag("CY", info['quality_barcode'], "Z")
                    if info['corrected_umi'] != '-':
                        read.set_tag("UB", info['corrected_umi'], "Z")
                    if info['uncorrected_umi'] != '-':
                        read.set_tag("UR", info['uncorrected_umi'], "Z")
                    if info['quality_umi'] != '-':
                        read.set_tag("UY", info['quality_umi'], "Z")
                    if info['gene'] != '-':
                        read.set_tag("GN", info['gene'], "Z")
                    if info['transcript'] != '-':
                        read.set_tag("TR", info['transcript'], "Z")
                    
                    # Adding transcriptome-specific tags
                    read.set_tag("TS", read.reference_name, "Z")  # Transcript name
                    read.set_tag("TL", read.reference_length, "i")  # Transcript length
                    
                    out_bam.write(read)
                else:
                    skipped_reads += 1
            else:
                skipped_reads += 1

    logger.info(f"Total reads processed: {total_reads}")
    logger.info(f"Reads tagged and written: {tagged_reads} ({tagged_reads/total_reads:.2%})")
    logger.info(f"Reads skipped: {skipped_reads} ({skipped_reads/total_reads:.2%})")
