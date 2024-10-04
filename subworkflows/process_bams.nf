import java.util.ArrayList;


process split_gtf_by_chroms {
    label "singlecell"
    cpus 1
    memory "1 GB"
    input:
        path("ref.gtf")
    output:
        path("*"), emit: chrom_gtf
    """
    gawk '/^[^#]/ {print>\$1".gtf"}' ref.gtf 
    """
}   


process generate_whitelist{
    label "singlecell"
    cpus 4
    memory "4 GB"
    publishDir "${params.out_dir}/${meta.alias}", mode: 'copy'
    input:
        tuple val(meta),
              path("barcodes/?_barcode.tsv")
    output:
        tuple val(meta),
              path("whitelist.tsv"),
              emit: whitelist
        tuple val(meta),
              path("kneeplot.png"),
              emit: kneeplot
    script:
    def no_thresholding_opt = meta.kit.split(':')[0] == 'visium' ? '--no_cell_filter' : ""
    """
    workflow-glue create_shortlist \
        barcodes whitelist.tsv \
        --counts \
        --method quantile \
        --exp_cells ${meta['expected_cells']} \
        --plot "kneeplot.png" \
        --counts_out "high_qual_bc_counts.tsv" \
        --threads ${task.cpus} \
        ${no_thresholding_opt}
    """
}


process assign_barcodes{
    label "singlecell"
    cpus 1
    memory "2 GB"
    input:
         tuple val(meta),
               path("whitelist.tsv"),
               path("extract_barcodes.tsv")
    output:
        tuple val(meta),
              path("bc_assign_counts.tsv"),
              emit: chrom_assigned_barcode_counts
        tuple val(meta),
              path("extract_barcodes_with_bc.tsv"),
              emit: tags
    """
    workflow-glue assign_barcodes \
        whitelist.tsv extract_barcodes.tsv \
        extract_barcodes_with_bc.tsv bc_assign_counts.tsv \
        --max_ed ${params.barcode_max_ed} \
        --min_ed_diff ${params.barcode_min_ed_diff} \
        --use_kmer_index
    """
}


process merge_bams {
    label "wf_common"
    cpus params.threads
    memory "8 GB"
    input:
        tuple val(meta),
            path('bams/*aln.bam'),
            path('bams/*aln.bam.bai')
    output:
        tuple val(meta),
              path("merged.sorted.bam"),
              path("merged.sorted.bam.bai"),
              emit: merged_bam
    script:
    """
    samtools merge -@ ${task.cpus -1} --write-index -o "merged.sorted.bam##idx##merged.sorted.bam.bai" bams/*.bam
    """
}


process cat_tags_by_chrom {
    label "wf_common"
    cpus params.threads
    memory "8 GB"
    input:
        tuple val(meta),
              path('tags/*tags.tsv')
    output:
         tuple val(meta),
              path("chr_tags/*"),
              emit: merged_tags

    """
    mkdir chr_tags
    files=(tags/*)
    chr_col=\$(awk -v RS='\t' '/chr/{print NR; exit}' "\${files[0]}")

    awk -F'\t' -v chr_col=\$chr_col 'FNR==1{hdr=\$0; next} \
    {if (!seen[\$chr_col]++) \
        print hdr>"chr_tags/"\$chr_col".tsv"; \
        print>"chr_tags/"\$chr_col".tsv"}' tags/*
    """
}


process stringtie {
    label "singlecell"
    cpus params.threads
    memory = { 3.GB * task.attempt }
    maxRetries = 3
    errorStrategy = { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    input:
        path 'ref_genome.fa'
        path 'ref_genome.fa.fai'
        tuple val(meta),
              path("align.bam"),
              path("align.bam.bai"),
              val(chr),
              path("chr.gtf")

    output:
        tuple val(meta),
              val(chr),
              path("transcriptome.fa"),
              path("chr.gtf"),
              path("stringtie.gff"),
              path("reads.fastq.gz"),
              emit: read_tr_map
    script:
    """
    samtools view -h align.bam ${chr}  \
        | tee >(
            stringtie -L ${params.stringtie_opts} -p ${task.cpus} \
                -G chr.gtf -l "${chr}.stringtie" -o "stringtie.gff" - ) \
        | samtools fastq \
        | bgzip --threads 2 -c > reads.fastq.gz
    gffread -g ref_genome.fa -w "transcriptome.fa" "stringtie.gff"
    """
}


process align_to_transcriptome {
    label "singlecell"
    cpus params.threads
    memory = "32 GB"
    input:
        tuple val(meta),
              val(chr),
              path('transcriptome.fa'),
              path('chr.gtf'),
              path('stringtie.gff'),
              path("reads.fq.gz")
    output:
        tuple val(meta),
              val(chr),
              path("chr.gtf"),
              path("tr_align.bam"),
              path("tr_align.bam.bai"),
              path('stringtie.gff'),
              emit: read_tr_map
    script:
    def view_threads = 1
    def sort_threads = 3
    def mm2_threads = Math.max(task.cpus - view_threads - sort_threads, 4)
    """
    minimap2 --eqx -N 100 -ax map-ont \
        --cap-kalloc 100m --cap-sw-mem 50m \
        --end-bonus 10 -p 0.9 -N 3 -t $mm2_threads \
        transcriptome.fa reads.fq.gz \
    | samtools view -h -@ $view_threads -b -F 2052 - \
    | samtools sort -n -@ $sort_threads --no-PG - > tr_align.bam

    samtools index -@ ${task.cpus} tr_align.bam
    """
}


process assign_features {
    label "singlecell"
    cpus 1
    memory { 1.0.GB.toBytes() + (tags.size() * 2 ) }
    input:
        tuple val(meta),
              val(chr),
              path("chr.gtf"),
              path("tr_align.bam"),
              path("stringtie.gff"),
              path(tags, stageAs: "tags.tsv")
    output:
        tuple val(meta),
              val(chr),
              path("feature_assigns.tsv"),
              emit: feature_assigns
        tuple val(meta),
              path("gffcompare.annotated.gtf"),
              emit: annotation
    """
    gffcompare -o gffcompare -r chr.gtf stringtie.gff

    workflow-glue assign_features \
        tr_align.bam \
        gffcompare.stringtie.gff.tmap \
        chr.gtf \
        tags.tsv \
        feature_assigns.tsv \
        --min_mapq ${params.gene_assigns_minqv}
    """
}


process tag_transcriptome_bam {
    label "singlecell"
    cpus 4
    memory "16 GB"
    publishDir "${params.out_dir}/${meta.alias}", mode: 'copy'
    input:
        tuple val(meta),
              val(chr),
              path('tr_align.bam'),
              path('tr_align.bam.bai'),
              path('feature_assigns.tsv')
    output:
        tuple val(meta),
              val(chr),
              path("tagged_tr_align.bam"),
              path("tagged_tr_align.bam.bai")
    script:
    """
    mkdir tags
    cp feature_assigns.tsv tags/
    workflow-glue tag_bam \
        tr_align.bam tagged_tr_align.bam tags \
        --threads ${task.cpus}
    samtools index -@ ${task.cpus} "tagged_tr_align.bam"
    """
}


process merge_tagged_transcriptome_bams {
    label "singlecell"
    cpus params.threads
    memory "8 GB"
    input:
        tuple val(meta),
              path('tagged_tr_align.bam'),
              path('tagged_tr_align.bam.bai')
    output:
        tuple val(meta),
              path("tagged_tr_align_merged.bam"),
              path("tagged_tr_align_merged.bam.bai"),
              emit: tagged_tr_bam
    script:
    """
    samtools merge -@ ${task.cpus -1} --write-index -o "tagged_tr_align_merged.bam##idx##tagged_tr_align_merged.bam.bai" tagged_tr_align.bam
    """
}


process create_matrix {
    label "singlecell"
    cpus 1
    memory {1.0.GB.toBytes()  + (read_tags.size() * 20) }
    input:
        tuple val(meta), val(chr), path("features.tsv"), path(read_tags, stageAs: "barcodes.tsv")
    output:
        tuple val(meta), val(chr), path("summary.tsv"), emit: summary
        tuple val(meta), val(chr), val("gene"), path("expression.gene.hdf"), emit: gene
        tuple val(meta), val(chr), val("transcript"), path("expression.transcript.hdf"), emit: transcript
        tuple val(meta), val(chr), path("stats.json"), emit: stats
    """
    workflow-glue create_matrix \
        ${chr} barcodes.tsv features.tsv \
        --tsv_out summary.tsv \
        --hdf_out expression.hdf \
        --stats stats.json
    """
}


process process_matrix {
    label "singlecell"
    cpus  1
    memory "16 GB"
    publishDir "${params.out_dir}/${meta.alias}", mode: 'copy', pattern: "*{mito,umap,raw,processed}*"
    input:
        tuple val(meta), val(feature), path('inputs/matrix*.hdf')
    output:
        tuple val(meta), val(feature), path("${feature}_raw_feature_bc_matrix"), emit: raw
        tuple val(meta), val(feature), path("${feature}_processed_feature_bc_matrix"), emit: processed
        tuple val(meta), val(feature), path("${feature}.expression.mean-per-cell.tsv"), emit: meancell
        tuple val(meta), val(feature), path("gene.expression.mito-per-cell.tsv"), emit: mitocell, optional: true
        tuple val(meta), val(feature), path("${feature}.expression.umap*.tsv"), emit: umap
    script:
    def mito_prefixes = params.mito_prefix.replaceAll(',', ' ')
    """
    export NUMBA_NUM_THREADS=${task.cpus}
    workflow-glue process_matrix \
        inputs/matrix*.hdf \
        --feature ${feature} \
        --raw ${feature}_raw_feature_bc_matrix \
        --processed ${feature}_processed_feature_bc_matrix \
        --per_cell_mito ${feature}.expression.mito-per-cell.tsv \
        --per_cell_expr ${feature}.expression.mean-per-cell.tsv \
        --umap_tsv ${feature}.expression.umap.tsv \
        --enable_filtering \
        --min_features $params.matrix_min_genes \
        --min_cells $params.matrix_min_cells \
        --max_mito $params.matrix_max_mito \
        --mito_prefixes $mito_prefixes \
        --norm_count $params.matrix_norm_count \
        --enable_umap \
        --replicates 3
    """
}


process merge_transcriptome {
    label "singlecell"
    cpus 2
    memory "2GB"
    publishDir "${params.out_dir}/${meta.alias}", mode: 'copy'
    input:
        tuple val(meta),
            path('fasta/?.fa'),
            path('gffs/?.gff')
    output:
        tuple val(meta),
            path("transcriptome.gff.gz"),
            path("transcriptome.fa.gz"),
            emit: merged_annotation
    """
    find fasta/ -name '*.fa' -exec cat {} + \
        | bgzip --threads ${task.cpus} -c  \
        > "transcriptome.fa.gz"
    find gffs/ -name '*.gff' -exec cat {} + \
        | grep -v '^#' \
        | bgzip --threads ${task.cpus} -c  \
        > "transcriptome.gff.gz"
    """
}


process combine_final_tag_files {
    label "singlecell"
    cpus 1
    memory "1 GB"
    publishDir "${params.out_dir}/${meta.alias}", mode: 'copy'
    input:
        tuple val(meta),
              path("tags*.tsv")
    output:
        tuple val(meta),
              path("read_summary.tsv")
    """
    awk 'FNR>1 || NR==1' *.tsv > "read_summary.tsv"
    """
}


process umi_gene_saturation {
    label "singlecell"
    cpus 4
    memory "32 GB"
    input:
        tuple val(meta),
              path("read_tags.tsv")
    output:
        tuple val(meta),
              path("saturation_curves.png"),
              emit: saturation_curve
    """
    export POLARS_MAX_THREADS=$task.cpus

    workflow-glue calc_saturation \
        --output "saturation_curves.png" \
        --read_tags read_tags.tsv
    """
}


process pack_images {
    label "singlecell"
    cpus 1
    memory "1 GB"
    input:
        tuple val(meta),
              path("images_${meta.alias}/*")
    output:
         tuple val(meta),
              path("images_${meta.alias}")
    """
    echo packing images
    """
}


process tag_bam {
    label "singlecell"
    cpus 4
    memory "16 GB"
    publishDir "${params.out_dir}/${meta.alias}", mode: 'copy'
    input:
        tuple val(meta), path('align.bam'), path('align.bam.bai'), path('tags/tag_*.tsv')
    output:
         tuple val(meta), path("tagged.bam"), path('tagged.bam.bai')
    script:
    """
    workflow-glue tag_bam \
        align.bam tagged.bam tags \
        --threads ${task.cpus}
    samtools index -@ ${task.cpus} "tagged.bam"
    """
}


process tag_tr_bam {
    label "singlecell"
    cpus 4
    memory "16 GB"
    publishDir "${params.out_dir}/${meta.alias}", mode: 'copy'
    input:
        tuple val(meta), val(chr), path('tr_align.bam'), path('tr_align.bam.bai'), path('feature_assigns.tsv')
    output:
        tuple val(meta), val(chr), path("tr_tagged.bam"), path('tr_tagged.bam.bai'), emit: tagged_tr_bam
    script:
    """
    workflow-glue tag_transcriptome_bam \
        tr_align.bam tr_tagged.bam feature_assigns.tsv \
        --threads ${task.cpus}
    samtools index -@ ${task.cpus} "tr_tagged.bam"
    """
}


workflow process_bams {
    take:
        bam
        extracted_barcodes
        high_qual_bc_counts
        gtf
        ref_genome_fasta
        ref_genome_idx
    main:
        chr_gtf = split_gtf_by_chroms(gtf)
            .flatten()
            .map {fname -> tuple(fname.baseName, fname)}  // [chr, gtf]

        generate_whitelist(high_qual_bc_counts)

        assign_barcodes(
            generate_whitelist.out.whitelist
            .cross(extracted_barcodes)
            .map {it ->
                meta = it[0][0]
                whitelist = it[0][1]
                barcodes = it[1][1]
                [meta, whitelist, barcodes]})

        chr_tags = cat_tags_by_chrom(assign_barcodes.out.tags.groupTuple())
            .transpose()
            .map {meta, file -> [meta, file.baseName, file]}

        merge_bams(bam.groupTuple())

        stringtie(
            ref_genome_fasta,
            ref_genome_idx,
            merge_bams.out.merged_bam
                .combine(chr_gtf))

        align_to_transcriptome(stringtie.out.read_tr_map)

        assign_features(
            align_to_transcriptome.out.read_tr_map
                .join(chr_tags, by: [0, 1])
        )

        // Tag the transcriptome-mapped BAM file
        tag_tr_bam(
            align_to_transcriptome.out.read_tr_map
                .join(assign_features.out.feature_assigns, by: [0,1])
        )

        create_matrix(
            assign_features.out.feature_assigns
                .join(chr_tags, by: [0, 1])
        )

        process_matrix(
            create_matrix.out.gene.groupTuple(by: [0, 2])
            .mix(
                create_matrix.out.transcript.groupTuple(by: [0, 2]))
            .map {meta, feature, hdfs -> [meta, feature, hdfs]})

        merge_transcriptome(
            assign_features.out.annotation.groupTuple()
                .join(stringtie.out.read_tr_map.groupTuple())
                .map{
                    meta, ann_tr_gff, chr, tr_fa, ref_gtf, str_gff, fastq ->
                    [meta, tr_fa, ann_tr_gff]})

        tags_by_sample = create_matrix.out.summary
            .groupTuple()
            .map{meta, chrs, files -> [meta, files]}

        final_read_tags = combine_final_tag_files(tags_by_sample)

        tag_bam(merge_bams.out.join(tags_by_sample))

        umi_gene_saturation(final_read_tags)

        pack_images(
            generate_whitelist.out.kneeplot
                .concat(umi_gene_saturation.out.saturation_curve)
                .groupTuple())
    
    emit:

        final_read_tags = final_read_tags
        plots = pack_images.out.collect{it -> it[1]}.collect()
        white_list = generate_whitelist.out.whitelist
        gene_mean_expression = process_matrix.out.meancell
            .filter{it[1] == "gene"}
            .map{it->[it[0], it[2]]}
        raw_gene_expression = process_matrix.out.raw
            .filter{it[1] == "gene"}
            .map{it->[it[0], it[2]]}
        transcript_mean_expression = process_matrix.out.meancell
            .filter{it[1] == "transcript"}
            .map{it->[it[0], it[2]]}
        mitochondrial_expression = process_matrix.out.mitocell
            .filter{it[1] == "gene"}
            .map{it->[it[0], it[2]]}
        umap_matrices = process_matrix.out.umap
            .map{it->[it[0], it[2]]}
            .groupTuple(size:2)
            .map{key, files -> [key, files.flatten()]}
        expression_stats = create_matrix.out.stats
        tagged_transcriptome_bam = merge_tagged_transcriptome_bams.out.tagged_tr_bam

        // Include the tagged transcriptome BAM file in the outputs
        tagged_tr_bam = tag_tr_bam.out.tagged_tr_bam

}