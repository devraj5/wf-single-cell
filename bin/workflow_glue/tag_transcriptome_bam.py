"""Tag transcriptome-aligned BAM file with information from read_summary.tsv"""

import pysam
import pandas as pd
import collections
import itertools
from umi_tools import UMIClusterer
from .util import get_named_logger, wf_parser

def argparser():
    parser = wf_parser("tag_transcriptome_bam")
    parser.add_argument("input_bam", help="Input transcriptome-aligned BAM file")
    parser.add_argument("output_bam", help="Output tagged BAM file")
    parser.add_argument("read_summary", help="read_summary.tsv file")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads to use")
    return parser

def get_adj_list_directional_lev(self, umis, counts, threshold=2):
    """Use Levenshtein distance for UMIclustering instead of hamming."""
    from editdistance import eval as edit_distance
    adj_list = {umi: [] for umi in umis}
    for umi1, umi2 in itertools.combinations(umis, 2):
        if edit_distance(umi1, umi2) <= threshold:
            if counts[umi1] >= (counts[umi2] * 2) - 1:
                adj_list[umi1].append(umi2)
            if counts[umi2] >= (counts[umi1] * 2) - 1:
                adj_list[umi2].append(umi1)
    return adj_list

def cluster_umis(umis):
    """Cluster UMIs using UMI-tools directional method."""
    if len(umis) == 1:
        return umis
    
    # Monkey-patch the UMI-tools clusterer
    UMIClusterer._get_adj_list_directional = get_adj_list_directional_lev
    
    clusterer = UMIClusterer(cluster_method="directional")
    umi_counts = collections.Counter(umis)
    clusters = clusterer(umi_counts, threshold=2)
    
    if len(clusters) == len(umis):
        return umis
    
    # Create mapping of corrected UMIs
    umi_map = {}
    for clust in clusters:
        if len(clust) > 1:
            for umi in clust[1:]:
                umi_map[umi] = clust[0]
    
    return [umi_map.get(umi, umi) for umi in umis]

def load_read_summary(file_path):
    """Load and process read summary with UMI deduplication."""
    # Read the TSV file
    df = pd.read_csv(file_path, sep='\t')
    
    # Filter for reads with both corrected barcode and UMI
    df = df[
        (df['corrected_barcode'] != '-') & 
        (df['corrected_umi'] != '-')
    ]
    
    # Group by barcode and gene for UMI deduplication
    df['gene_cell'] = df['gene'] + ':' + df['corrected_barcode']
    grouped = df.groupby('gene_cell')
    
    # Perform UMI clustering for each group
    for name, group in grouped:
        umis = group['uncorrected_umi'].tolist()
        corrected_umis = cluster_umis(umis)
        df.loc[group.index, 'corrected_umi'] = corrected_umis
    
    # Convert to dictionary format
    return df.set_index('read_id').to_dict('index')

def tag_bam(input_bam, output_bam, read_info, threads, logger):
    with pysam.AlignmentFile(input_bam, "rb", threads=threads) as in_bam, \
         pysam.AlignmentFile(output_bam, "wb", template=in_bam, threads=threads) as out_bam:
        
        total_reads = 0
        written_reads = 0
        skipped_reads = 0
        
        for read in in_bam:
            total_reads += 1
            read_id = read.query_name
            
            # Only process reads that are in our filtered read_info
            if read_id in read_info:
                info = read_info[read_id]
                
                # Set all required tags
                read.set_tag("CB", info['corrected_barcode'], "Z")
                read.set_tag("CR", info['uncorrected_barcode'], "Z")
                read.set_tag("CY", info.get('quality_barcode', '-'), "Z")
                read.set_tag("UB", info['corrected_umi'], "Z")
                read.set_tag("UR", info['uncorrected_umi'], "Z")
                read.set_tag("UY", info.get('quality_umi', '-'), "Z")
                read.set_tag("GN", info.get('gene', '-'), "Z")
                read.set_tag("TR", info.get('transcript', '-'), "Z")
                
                out_bam.write(read)
                written_reads += 1
            else:
                skipped_reads += 1

    logger.info(f"Total reads processed: {total_reads}")
    logger.info(f"Reads written: {written_reads} ({written_reads/total_reads:.2%})")
    logger.info(f"Reads skipped: {skipped_reads} ({skipped_reads/total_reads:.2%})")

def main(args):
    logger = get_named_logger("TagTranscriptomeBAM")
    
    logger.info("Loading and processing read summary with UMI deduplication")
    read_info = load_read_summary(args.read_summary)
    
    logger.info("Tagging BAM file")
    tag_bam(args.input_bam, args.output_bam, read_info, args.threads, logger)
