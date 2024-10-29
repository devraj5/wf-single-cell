"""Tag transcriptome-aligned BAM file with information from read_summary.tsv and deduplicate UMIs"""

import collections
import itertools
from pathlib import Path
import pysam
import csv
import pandas as pd
from editdistance import eval as edit_distance
from umi_tools import UMIClusterer
from .util import get_named_logger, wf_parser

def argparser():
    parser = wf_parser("tag_transcriptome_bam")
    parser.add_argument("input_bam", help="Input transcriptome-aligned BAM file")
    parser.add_argument("output_bam", help="Output tagged BAM file")
    parser.add_argument("read_summary", help="read_summary.tsv file")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads to use")
    return parser

def get_adj_list_directional_lev(umis, counts, threshold=2):
    """Use Levenshtein distance for UMI clustering instead of hamming."""
    adj_list = {umi: [] for umi in umis}
    iter_umi_pairs = itertools.combinations(umis, 2)
    for umi1, umi2 in iter_umi_pairs:
        if edit_distance(umi1, umi2) <= threshold:
            if counts[umi1] >= (counts[umi2] * 2) - 1:
                adj_list[umi1].append(umi2)
            if counts[umi2] >= (counts[umi1] * 2) - 1:
                adj_list[umi2].append(umi1)
    return adj_list

def cluster_umis(umis):
    """Cluster UMIs using UMI-tools directional method."""
    if len(umis) == 1:
        return umis[0]
    
    clusterer = UMIClusterer(cluster_method="directional")
    clusterer._get_adj_list_directional = get_adj_list_directional_lev
    
    umi_counts = collections.Counter(umis)
    clusters = clusterer(umi_counts, threshold=2)
    
    # Return representative UMI from largest cluster
    return clusters[0][0]

def load_and_process_read_summary(file_path):
    """Load and process read summary, performing UMI deduplication."""
    df = pd.read_csv(file_path, sep='\t')
    
    # Create gene/cell index for UMI deduplication
    df['gene_cell_umi'] = df.apply(
        lambda x: f"{x['gene']}:{x['corrected_barcode']}:{x['corrected_umi']}", 
        axis=1
    )
    
    # Group by gene and cell barcode
    grouped = df.groupby(['gene', 'corrected_barcode'])
    
    # Perform UMI deduplication within each group
    deduplicated_reads = set()
    for _, group in grouped:
        umis = group['corrected_umi'].tolist()
        if len(umis) > 0:
            representative_umi = cluster_umis(umis)
            # Keep the first read with the representative UMI
            read_to_keep = group[group['corrected_umi'] == representative_umi].iloc[0].name
            deduplicated_reads.add(read_to_keep)
    
    # Convert to dictionary for quick lookup
    read_info = {}
    for _, row in df.iterrows():
        if row.name in deduplicated_reads:
            read_info[row['read_id']] = {
                'CB': row['corrected_barcode'],
                'CR': row['uncorrected_barcode'],
                'CY': row['quality_barcode'],
                'UB': row['corrected_umi'],
                'UR': row['uncorrected_umi'],
                'UY': row['quality_umi'],
                'GN': row['gene'],
                'TR': row['transcript']
            }
    
    return read_info

def tag_bam(input_bam, output_bam, read_info, threads, logger):
    with pysam.AlignmentFile(input_bam, "rb", threads=threads) as in_bam, \
         pysam.AlignmentFile(output_bam, "wb", template=in_bam, threads=threads) as out_bam:
        
        total_reads = 0
        written_reads = 0
        skipped_reads = 0
        
        for read in in_bam:
            total_reads += 1
            read_id = read.query_name
            
            if read_id in read_info:
                written_reads += 1
                info = read_info[read_id]
                
                # Set all required tags, using '-' for missing values
                tags_to_set = {
                    'CB': info.get('CB', '-'),
                    'CR': info.get('CR', '-'),
                    'CY': info.get('CY', '-'),
                    'UB': info.get('UB', '-'),
                    'UR': info.get('UR', '-'),
                    'UY': info.get('UY', '-'),
                    'GN': info.get('GN', '-'),
                    'TR': info.get('TR', '-')
                }
                
                for tag, value in tags_to_set.items():
                    read.set_tag(tag, value, "Z")
                
                out_bam.write(read)
            else:
                skipped_reads += 1

        logger.info(f"Total reads processed: {total_reads}")
        logger.info(f"Unique reads after deduplication: {written_reads} ({written_reads/total_reads:.2%})")
        logger.info(f"Duplicate/skipped reads: {skipped_reads} ({skipped_reads/total_reads:.2%})")

def main(args):
    logger = get_named_logger("TagTranscriptomeBAM")
    
    logger.info("Loading and processing read summary file...")
    read_info = load_and_process_read_summary(args.read_summary)
    
    logger.info("Tagging BAM file with deduplicated reads...")
    tag_bam(args.input_bam, args.output_bam, read_info, args.threads, logger)
