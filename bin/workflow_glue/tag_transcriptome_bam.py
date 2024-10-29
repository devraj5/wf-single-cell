"""Tag transcriptome-aligned BAM file with information from read_summary.tsv and perform UMI deduplication."""

import collections
import itertools
from pathlib import Path
import pysam
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
    """Cluster UMIs using directional method with Levenshtein distance."""
    if len(umis) == 1:
        return umis
    
    clusterer = UMIClusterer(cluster_method="directional")
    clusterer._get_adj_list_directional = get_adj_list_directional_lev
    
    umi_counts = collections.Counter(umis)
    clusters = clusterer(umi_counts, threshold=2)
    
    if len(clusters) == len(umis):
        return umis
    
    # Create mapping of corrected UMIs
    umi_map = {}
    for cluster in clusters:
        if len(cluster) > 1:
            for umi in cluster[1:]:
                umi_map[umi] = cluster[0]
    
    return pd.Series(umis).replace(umi_map) if umi_map else umis

def deduplicate_and_tag_bam(input_bam, output_bam, read_info, threads, logger):
    """Process BAM file with UMI deduplication and tagging."""
    # First pass: group reads by barcode-gene-transcript combination and deduplicate UMIs
    read_groups = collections.defaultdict(list)
    
    with pysam.AlignmentFile(input_bam, "rb", threads=threads) as in_bam:
        for read in in_bam:
            read_id = read.query_name
            if read_id in read_info:
                info = read_info[read_id]
                if info['corrected_barcode'] != '-' and info['corrected_umi'] != '-':
                    key = (
                        info['corrected_barcode'],
                        info.get('gene', '-'),
                        info.get('transcript', '-')
                    )
                    read_groups[key].append((read_id, info['corrected_umi']))

    # Perform UMI deduplication within each group
    kept_reads = set()
    for group_reads in read_groups.values():
        read_ids, umis = zip(*group_reads)
        deduplicated_umis = cluster_umis(umis)
        # Keep only the first occurrence of each UMI after deduplication
        seen_umis = set()
        for read_id, umi in zip(read_ids, deduplicated_umis):
            if umi not in seen_umis:
                kept_reads.add(read_id)
                seen_umis.add(umi)

    # Second pass: write deduplicated and tagged BAM
    with pysam.AlignmentFile(input_bam, "rb", threads=threads) as in_bam, \
         pysam.AlignmentFile(output_bam, "wb", template=in_bam, threads=threads) as out_bam:
        
        total_reads = 0
        written_reads = 0
        
        for read in in_bam:
            total_reads += 1
            read_id = read.query_name
            
            if read_id in kept_reads:
                info = read_info[read_id]
                
                # Set required tags
                read.set_tag("CB", info.get('corrected_barcode', '-'), "Z")
                read.set_tag("CR", info.get('uncorrected_barcode', '-'), "Z")
                read.set_tag("CY", info.get('quality_barcode', '-'), "Z")
                read.set_tag("UB", info.get('corrected_umi', '-'), "Z")
                read.set_tag("UR", info.get('uncorrected_umi', '-'), "Z")
                read.set_tag("UY", info.get('quality_umi', '-'), "Z")
                read.set_tag("GN", info.get('gene', '-'), "Z")
                read.set_tag("TR", info.get('transcript', '-'), "Z")
                
                out_bam.write(read)
                written_reads += 1

    logger.info(f"Total reads processed: {total_reads}")
    logger.info(f"Unique reads after deduplication: {written_reads} ({written_reads/total_reads:.2%})")

def main(args):
    logger = get_named_logger("TagTranscriptomeBAM")
    
    # Load read summary information
    logger.info("Loading read summary information")
    read_info = {}
    with open(args.read_summary, 'r') as f:
        for row in pd.read_csv(f, sep='\t').to_dict('records'):
            read_info[row['read_id']] = row
    
    # Process BAM file
    logger.info("Processing BAM file")
    deduplicate_and_tag_bam(args.input_bam, args.output_bam, read_info, args.threads, logger)