"""Correct UMIs in tagged BAM files using clustering."""

import collections
import itertools
import pysam
from editdistance import eval as edit_distance
from umi_tools import UMIClusterer
from .util import get_named_logger, wf_parser  # Import from workflow-glue utils

def argparser():
    """Create argument parser."""
    parser = wf_parser("correct_umis")
    parser.add_argument("-i", "--input", required=True, help="Input BAM file")
    parser.add_argument("-o", "--output", required=True, help="Output BAM file")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads for sorting")
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
    
    adj_list = self.get_adj_list(umis, counts, threshold)
    clusters = self.get_connected_components(umis, adj_list, counts)
    final_umis = [list(x) for x in self.get_groups(clusters, adj_list, counts)]
    
    return final_umis

# Monkey-patch the UMI-tools clusterer
UMIClusterer._get_adj_list_directional = get_adj_list_directional_lev
UMIClusterer.__call__ = umi_clusterer_call

def cluster_umis(umis):
    """Cluster UMIs using directional method."""
    if len(umis) == 1:
        return umis[0]
    
    clusterer = UMIClusterer(cluster_method="directional")
    umi_counts = collections.Counter(umis)
    clusters = clusterer(umi_counts, threshold=2)
    
    # Return representative UMI from largest cluster
    return clusters[0][0]

def process_bam(in_bam, out_bam):
    """Process BAM file to correct UMIs."""
    logger = get_named_logger("UMICorrect")
    logger.info(f"Processing BAM file: {in_bam}")
    
    with pysam.AlignmentFile(in_bam, "rb") as input_bam:
        with pysam.AlignmentFile(out_bam, "wb", header=input_bam.header) as output_bam:
            current_group = []
            current_key = None
            
            for read in input_bam.fetch(until_eof=True):
                cell_barcode = read.get_tag("CB")
                gene = read.get_tag("GN")
                key = f"{cell_barcode}_{gene}"
                
                if key != current_key:
                    if current_group:
                        umis = [read.get_tag("UB") for read in current_group]
                        corrected_umi = cluster_umis(umis)
                        for read in current_group:
                            read.set_tag("UB", corrected_umi)
                            output_bam.write(read)
                    
                    current_group = [read]
                    current_key = key
                else:
                    current_group.append(read)
            
            if current_group:
                umis = [read.get_tag("UB") for read in current_group]
                corrected_umi = cluster_umis(umis)
                for read in current_group:
                    read.set_tag("UB", corrected_umi)
                    output_bam.write(read)

def main(args):
    """Entry point for correct_umis command."""
    process_bam(args.input, args.output)