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

process generate_whitelist {
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
    // Note: This is called "uncorrected", but they're actually counts of
    //       high-quality exact matches to longlist. Low-frequency barcodes
    //       are assumed to be false positives. The list is further
    //       filtered by the selected method (basically by abundance).
    // TODO: change this to take precomputed, filtered counts from extract_barcodes
    script:
        // It doesn't make sense to do cell count thresholding of the shortlist for Visium data.
        // A Visium barcode is a tissue coordinate, not a cell.
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

process assign_barcodes {
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
    // Combine all BAMs derived from the initial chunking into per-sample files
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
        samtools merge -@ ${task.cpus - 1} --write-index -o "merged.sorted.bam##idx##merged.sorted.bam.bai" bams/*.bam
        """
}

process cat_tags_by_chrom {
    // Merge per-chunk tags to create per-chromosome tags
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
    # Find the chr column number
    files=(tags/*)
    chr_col=\$(awk -v RS='\t' '/chr/{print NR; exit}' "\${files[0]}")

    # Merge the tags TSVs, keep header from first file and split entries by chromosome
    awk -F'\t' -v chr_col=\$chr_col 'FNR==1{hdr=\$0; next} \
    {if (!seen[\$chr_col]++) \
        print hdr>"chr_tags/"\$chr_col".tsv"; \
        print>"chr_tags/"\$chr_col".tsv"}' tags/*
    """
}

process stringtie {
    label "singlecell"
    cpus params.threads
    // Memory usage for this process is usually less than 3GB, but in some cases, it may go over this.
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
        # Add chromosome label (-l) to generated transcripts
        # so we don't get name collisions during file merge later
        samtools view -h align.bam ${chr}  \
            | tee >(
                stringtie -L ${params.stringtie_opts} -p ${task.cpus} \
                    -G chr.gtf -l "${chr}.stringtie" -o "stringtie.gff" - ) \
            | samtools fastq \
            | bgzip --threads 2 -c > reads.fastq.gz
        # Get transcriptome sequence
        gffread -g ref_genome.fa -w "transcriptome.fa" "stringtie.gff"
        """
}

process align_to_transcriptome {
    label "singlecell"
    cpus params.threads
    memory "32 GB"
    input:
        tuple val(meta),
              val(chr),
              path('transcriptome.fa'),
              path('chr.gtf'),
              path('stringtie.gff'),
              path("reads.fastq.gz")
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
        minimap2 -ax map-ont --eqx -N 100 \
            --cap-kalloc 100m --cap-sw-mem 50m \
            --end-bonus 10 -p 0.9 -t $mm2_threads \
            transcriptome.fa reads.fastq.gz \
        | samtools view -h -@ $view_threads -b -F 2052 - \
        | samtools sort -n -@ $sort_threads --no-PG - > tr_align.bam

        samtools index -@ ${task.cpus} tr_align.bam
        """
}

process assign_features {
    label "singlecell"
    cpus 1
    // This step is performed per chromosome. The tags file per chrom can vary
    // quite widely in size. We don't have a fixed memory size here to get better
    // parallelism on single-host setups.
    memory { 1.0.GB.toBytes() + (tags.size() * 2) }
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
    # gffcompare maps transcript reference IDs to query transcripts.
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

process tag_tr_bam {
    label "singlecell"
    cpus 4
    memory "16 GB"
    publishDir "${params.out_dir}/${meta.alias}", mode: 'copy'
    input:
        tuple val(meta),
              val(chr),
              path('tr_align.bam'),
              path('tr_align.bam.bai'),
              path("feature_assigns.tsv"),
              path("extract_barcodes_with_bc.tsv")
    output:
        tuple val(meta),
              val(chr),
              path("tr_align_tagged.bam"),
              path("tr_align_tagged.bam.bai"),
              emit: tagged_tr_bam
    script:
        """
        # Combine feature assignments and barcodes into a single tags file
        mkdir tags
        paste extract_barcodes_with_bc.tsv feature_assigns.tsv > tags/tags.tsv

        workflow-glue tag_bam \
            tr_align.bam tr_align_tagged.bam tags \
            --threads ${task.cpus}
        samtools index -@ ${task.cpus} tr_align_tagged.bam
        """
}

// Create expression matrices by combining barcode and feature
// tag files. Also outputs the combined tags (per chromosome) to be combined later
process create_matrix {
    label "singlecell"
    cpus 1
    // Benchmarking showed that memory usage was ~15x the size of read_tags input.
    // Set a minimum memory requirement of 1.0GB to allow for overhead.
    memory { 1.0.GB.toBytes() + (read_tags.size() * 20) }
    input:
        tuple val(meta),
              val(chr),
              path("features.tsv"),
              path(read_tags, stageAs: "barcodes.tsv")
    output:
        tuple val(meta),
              val(chr),
              path("summary.tsv"),
              emit: summary
        tuple val(meta),
              val(chr),
              val("gene"),
              path("expression.gene.hdf"),
              emit: gene
        tuple val(meta),
              val(chr),
              val("transcript"),
              path("expression.transcript.hdf"),
              emit: transcript
        tuple val(meta),
              val(chr),
              path("stats.json"),
              emit: stats
    """
    workflow-glue create_matrix \
        ${chr} barcodes.tsv features.tsv \
        --tsv_out summary.tsv \
        --hdf_out expression.hdf \
        --stats stats.json
    """
}

// Combines multiple expression matrices (e.g., from different chromosomes)
// and calculates summary information on the matrix, including UMAPs
process process_matrix {
    label "singlecell"
    cpus 1
    memory "16 GB"
    publishDir "${params.out_dir}/${meta.alias}", mode: 'copy', pattern: "*{mito,umap,raw,processed}*"
    input:
        tuple val(meta),
              val(feature),
              path('inputs/matrix*.hdf')
    output:
        tuple val(meta),
              val(feature),
              path("${feature}_raw_feature_bc_matrix"),
              emit: raw
        tuple val(meta),
              val(feature),
              path("${feature}_processed_feature_bc_matrix"),
              emit: processed
        tuple val(meta),
              val(feature),
              path("${feature}.expression.mean-per-cell.tsv"),
              emit: meancell
        // Mito per cell makes sense only for feature=gene for now.
        tuple val(meta),
              val(feature),
              path("gene.expression.mito-per-cell.tsv"),
              emit: mitocell,
              optional: true
        tuple val(meta),
              val(feature),
              path("${feature}.expression.umap*.tsv"),
              emit: umap
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

// Merge annotated GFFs and transcriptome sequence files
process merge_transcriptome {
    label "singlecell"
    cpus 2
    memory "2 GB"
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
        | bgzip --threads ${task.cpus} -c \
        > "transcriptome.fa.gz"
    find gffs/ -name '*.gff' -exec cat {} + \
        | grep -v '^#' \
        | bgzip --threads ${task.cpus} -c \
        > "transcriptome.gff.gz"
    """
}

process combine_final_tag_files {
    // Create final per-sample read summaries with information from all stages
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
        tuple val(meta),
              path('align.bam'),
              path('align.bam.bai'),
              path('tags/tag_*.tsv')
    output:
        tuple val(meta),
              path("tagged.bam"),
              path('tagged.bam.bai')
    script:
        """
        workflow-glue tag_bam \
            align.bam tagged.bam tags \
            --threads ${task.cpus}
        samtools index -@ ${task.cpus} "tagged.bam"
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
        // Split the GTF by chromosome
        chr_gtf = split_gtf_by_chroms(gtf)
            .flatten()
            .map { fname -> tuple(fname.baseName, fname) }  // [chr, gtf]

        generate_whitelist(high_qual_bc_counts)

        assign_barcodes(
            generate_whitelist.out.whitelist
                .cross(extracted_barcodes)
                .map { it ->
                    meta = it[0][0]
                    whitelist = it[0][1]
                    barcodes = it[1][1]
                    [meta, whitelist, barcodes]
                }
        )

        // Combine the tag chunks to per-chromosome chunks and emit [meta, chr, tags]
        chr_tags = cat_tags_by_chrom(assign_barcodes.out.tags.groupTuple())
            .transpose()
            .map { meta, file -> [meta, file.baseName, file] }

        // Combine the BAM chunks per sample
        merge_bams(bam.groupTuple())

        // Run StringTie per chromosome
        stringtie(
            ref_genome_fasta,
            ref_genome_idx,
            merge_bams.out.merged_bam
                .combine(chr_gtf)
        )

        align_to_transcriptome(stringtie.out.read_tr_map)

        assign_features(
            align_to_transcriptome.out.read_tr_map
                .join(chr_tags, by: [0, 1])
        )

        // Tag the tr_align.bam
        tagged_tr_bam = tag_tr_bam(
            align_to_transcriptome.out.read_tr_map
                .join(assign_features.out.feature_assigns)
                .join(assign_barcodes.out.tags)
        )

        create_matrix(
            assign_features.out.feature_assigns
                // Join on [sample meta, chr]
                .join(chr_tags, by: [0, 1])
        )

        // Aggregate per-chromosome expression matrices to create MEX and UMAP TSVs
        process_matrix(
            create_matrix.out.gene.groupTuple(by: [0, 2])
                .mix(
                    create_matrix.out.transcript.groupTuple(by: [0, 2])
                )
                .map { meta, chroms, feature, hdfs -> [meta, feature, hdfs] }
        )

        merge_transcriptome(
            assign_features.out.annotation.groupTuple()
                .join(stringtie.out.read_tr_map.groupTuple())
                .map {
                    meta, ann_tr_gff, chr, tr_fa, ref_gtf, str_gff, fastq ->
                    [meta, tr_fa, ann_tr_gff]
                }
        )

        tags_by_sample = create_matrix.out.summary
            .groupTuple()
            .map { meta, chrs, files -> [meta, files] }
        final_read_tags = combine_final_tag_files(tags_by_sample)
        tag_bam(merge_bams.out.join(tags_by_sample))

        umi_gene_saturation(final_read_tags)

        pack_images(
            generate_whitelist.out.kneeplot
                .concat(umi_gene_saturation.out.saturation_curve)
                .groupTuple()
        )

    emit:
        final_read_tags = final_read_tags
        plots = pack_images.out.collect { it -> it[1] }.collect()
        white_list = generate_whitelist.out.whitelist
        gene_mean_expression = process_matrix.out.meancell
            .filter { it[1] == "gene" }
            .map { it -> [it[0], it[2]] }
        raw_gene_expression = process_matrix.out.raw
            .filter { it[1] == "gene" }
            .map { it -> [it[0], it[2]] }
        transcript_mean_expression = process_matrix.out.meancell
            .filter { it[1] == "transcript" }
            .map { it -> [it[0], it[2]] }
        mitochondrial_expression = process_matrix.out.mitocell
            .filter { it[1] == "gene" }
            .map { it -> [it[0], it[2]] }
        umap_matrices = process_matrix.out.umap
            .map { it -> [it[0], it[2]] }
            .groupTuple(size: 2)
            .map { key, files -> [key, files.flatten()] }
        expression_stats = create_matrix.out.stats
        tagged_tr_align_bam = tagged_tr_bam.out.tagged_tr_bam.collect()
}
