"""Tag transcriptome-aligned BAM file with information from read_summary.tsv"""

import pysam
import pandas as pd
import collections
from umi_tools import UMIClusterer
from .util import get_named_logger, wf_parser

def argparser():
    parser = wf_parser("tag_transcriptome_bam")
    parser.add_argument("input_bam", help="Input transcriptome-aligned BAM file")
    parser.add_argument("output_bam", help="Output tagged BAM file")
    parser.add_argument("read_summary", help="read_summary.tsv file")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads to use")
    return parser

def process_read_summary(file_path):
    """Load, filter, and deduplicate reads."""
    # Reading needed columns
    df = pd.read_csv(
        file_path, 
        sep='\t',
        usecols=['read_id', 'corrected_barcode', 'uncorrected_barcode', 'quality_barcode',
                 'corrected_umi', 'uncorrected_umi', 'quality_umi', 'gene', 'transcript']
    )
    
    # Filtering for valid barcodes and UMIs
    mask = (df['corrected_barcode'] != '-') & (df['corrected_umi'] != '-')
    df = df[mask]
    
    # Grouping by barcode and gene for UMI deduplication
    df['gene_cell'] = df['gene'] + ':' + df['corrected_barcode']
    
    # Deduplicating UMIs
    clusterer = UMIClusterer(cluster_method="directional")
    dedup_rows = []
    
    for _, group in df.groupby('gene_cell'):
        if len(group) > 1:
            umis = group['uncorrected_umi'].tolist()
            umi_counts = collections.Counter(umis)
            clusters = clusterer(umi_counts, threshold=2)
            
            # Keeping only one read per UMI cluster
            for cluster in clusters:
                dedup_rows.append(group[group['uncorrected_umi'] == cluster[0]].iloc[0])
        else:
            dedup_rows.append(group.iloc[0])
    
    dedup_df = pd.concat(dedup_rows, axis=1).T
    
    # Converting to dictionary for faster lookup
    read_dict = {}
    for _, row in dedup_df.iterrows():
        read_dict[row['read_id']] = {
            'CB': row['corrected_barcode'],
            'CR': row['uncorrected_barcode'],
            'CY': row['quality_barcode'],
            'UB': row['corrected_umi'],
            'UR': row['uncorrected_umi'],
            'UY': row['quality_umi'],
            'GN': row['gene'],
            'TR': row['transcript']
        }
    
    return read_dict

def tag_bam(input_bam, output_bam, read_info, threads, logger):
    with pysam.AlignmentFile(input_bam, "rb", threads=threads) as in_bam, \
         pysam.AlignmentFile(output_bam, "wb", template=in_bam, threads=threads) as out_bam:
        
        written = 0
        for read in in_bam:
            if read.query_name in read_info:
                tags = read_info[read.query_name]
                for tag, value in tags.items():
                    read.set_tag(tag, value if value != '-' else '-', "Z")
                out_bam.write(read)
                written += 1
        
        logger.info(f"Written {written:,} reads to output BAM")

def main(args):
    logger = get_named_logger("TagTranscriptomeBAM")
    
    logger.info("Processing read summary")
    read_info = process_read_summary(args.read_summary)
    
    logger.info("Tagging BAM file")
    tag_bam(args.input_bam, args.output_bam, read_info, args.threads, logger)