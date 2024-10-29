"""Tag transcriptome-aligned BAM file with information from read_summary.tsv and perform UMI deduplication."""

import collections
import itertools
from pathlib import Path
import pysam
import pandas as pd
import numpy as np
from editdistance import eval as edit_distance
from umi_tools import UMIClusterer
from .util import get_named_logger, wf_parser

def argparser():
    parser = wf_parser("tag_transcriptome_bam")
    parser.add_argument("input_bam", help="Input transcriptome-aligned BAM file")
    parser.add_argument("output_bam", help="Output tagged BAM file")
    parser.add_argument("read_summary", help="read_summary.tsv file")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads to use")
    parser.add_argument(
        "--ref_interval", type=int, default=1000,
        help="Size of genomic window (bp) to assign as gene name if no gene assigned")
    return parser

def get_adj_list_directional_lev(self, umis, counts, threshold=2):
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

def umi_clusterer_call(self, umis, threshold):
    """Modified UMIClusterer.__call__ to handle variable length UMIs."""
    counts = umis
    umis = list(umis.keys())

    self.positions += 1
    number_of_umis = len(umis)
    self.total_umis_per_position += number_of_umis

    if number_of_umis > self.max_umis_per_position:
        self.max_umis_per_position = number_of_umis

    adj_list = self.get_adj_list(umis, counts, threshold)
    clusters = self.get_connected_components(umis, adj_list, counts)
    final_umis = [list(x) for x in self.get_groups(clusters, adj_list, counts)]

    return final_umis

def create_region_name(chr_name, start, end, ref_interval):
    """Create a fake gene name from alignment coordinates."""
    midpoint = int((start + end) / 2)
    interval_start = int(np.floor(midpoint / ref_interval) * ref_interval)
    interval_end = int(np.ceil(midpoint / ref_interval) * ref_interval)
    return f"{chr_name}_{interval_start}_{interval_end}"

def cluster(umis):
    """Cluster UMIs exactly as in create_matrix."""
    if len(umis) == 1:  # early return
        return umis
    
    # Monkey-patch the UMI-tools clusterer
    UMIClusterer._get_adj_list_directional = get_adj_list_directional_lev
    UMIClusterer.__call__ = umi_clusterer_call
    
    clusterer = UMIClusterer(cluster_method="directional")
    umi_counts = collections.Counter(umis)
    clusters = clusterer(umi_counts, threshold=2)

    if len(clusters) == len(umis):  # no corrections
        return umis

    # create list of corrections
    umi_map = dict()
    for clust in clusters:
        if len(clust) > 1:
            for umi in clust[1:]:
                umi_map[umi] = clust[0]
                
    if len(umi_map) > 0:
        umis = pd.Series(umis).replace(umi_map)
    return umis

def process_reads(input_bam, read_info, ref_interval):
    """Process reads to match create_matrix logic."""
    # Create dataframe from read info
    df = pd.DataFrame.from_dict(read_info, orient='index')
    
    # Handle unassigned genes
    df['no_gene'] = False
    df_no_gene = df.loc[df.gene == '-']
    if len(df_no_gene) > 0:
        # Create temporary gene names based on location
        regions = df_no_gene.apply(
            lambda x: create_region_name(x.chr, x.start, x.end, ref_interval), 
            axis=1)
        df.loc[regions.index, 'gene'] = regions
        df.loc[df.index.isin(regions.index), 'no_gene'] = True

    # Create gene/cell index for clustering
    df["gene_cell"] = df["gene"] + ":" + df["corrected_barcode"]
    df['read_id'] = df.index
    df.set_index('gene_cell', inplace=True, drop=True)
    
    # Cluster UMIs
    groups = df.groupby("gene_cell")["uncorrected_umi"]
    df["corrected_umi"] = groups.transform(cluster)
    df.set_index('read_id', drop=True, inplace=True)
    
    # Reset unassigned genes to '-'
    df.loc[df.no_gene, 'gene'] = '-'
    
    return df

def deduplicate_and_tag_bam(input_bam, output_bam, read_info, ref_interval, threads, logger):
    """Process BAM file with UMI deduplication matching create_matrix."""
    # Process reads using create_matrix logic
    processed_df = process_reads(input_bam, read_info, ref_interval)
    
    # Create set of unique combinations after UMI clustering
    unique_combinations = set(
        zip(processed_df['corrected_barcode'], 
            processed_df['gene'],
            processed_df['corrected_umi']))
    
    # Write deduplicated and tagged BAM
    with pysam.AlignmentFile(input_bam, "rb", threads=threads) as in_bam, \
         pysam.AlignmentFile(output_bam, "wb", template=in_bam, threads=threads) as out_bam:
        
        total_reads = 0
        written_reads = 0
        seen_combinations = set()
        
        for read in in_bam:
            total_reads += 1
            read_id = read.query_name
            
            if read_id in processed_df.index:
                info = processed_df.loc[read_id]
                combination = (info['corrected_barcode'], 
                             info['gene'],
                             info['corrected_umi'])
                
                # Only write first occurrence of each unique combination
                if combination in unique_combinations and combination not in seen_combinations:
                    seen_combinations.add(combination)
                    
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
    
    # Load read summary
    read_info = {}
    df = pd.read_csv(args.read_summary, sep='\t')
    read_info = df.set_index('read_id').to_dict('index')
    
    # Process BAM file
    deduplicate_and_tag_bam(
        args.input_bam, 
        args.output_bam, 
        read_info, 
        args.ref_interval,
        args.threads, 
        logger)