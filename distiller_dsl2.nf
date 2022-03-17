#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Define help message
def helpMessage() {
    log.info"""
    Distiller is a Hi-C mapping pipeline.
    """.stripIndent()
}

// Show help message
if (params.getOrDefault('help', 'false').toBoolean()) {
    helpMessage()
    exit 0
}

// Include important functions:
include { needsDownloading; checkLeftRightChunk } from './modules/local/functions'

// Check required parameters
MIN_RES = params['bin'].resolutions.collect { it as int }.min() // TODO: move to parameters dictionary
ASSEMBLY_NAME = params['input'].genome.assembly_name // TODO: move to the parameters dictionary

// Include modules and subworkflows
include { DOWNLOAD_TRUNCATE as download_truncate_chunk_fastqs } from './modules/local/download_truncate' addParams( options: [:] )
include { LOCAL_TRUNCATE as local_truncate_chunk_fastqs } from './modules/local/local_truncate' addParams( options: [:] )
include { FASTQC as fastqc } from './modules/local/fastqc' addParams( options: [:] )
include { MAP_PARSE_SORT as map_parse_sort_chunks } from './modules/local/map_parse_sort' addParams( options: [:] )
include { MERGE_DEDUP_SPLITBAM as merge_dedup_splitbam } from './modules/local/merge_dedup_splitbam' addParams( options: [:] )
include { BIN_ZOOM as bin_zoom_library_pairs } from './modules/local/bin_zoom_library_pairs' addParams( options: [:] )
include { MERGE_ZOOM as merge_zoom_library_group_coolers } from './modules/local/merge_zoom_library_group_coolers' addParams( options: [:] )
include { MERGE_STATS as merge_stats_libraries_into_groups } from './modules/local/merge_stats_libraries_into_groups' addParams( options: [:] )
include { MULTIQC } from './modules/local/multiqc/main.nf' addParams( options: [:] )

// Define workflow
workflow DISTILLER {

    // CHROM_SIZES:
    CHROM_SIZES = Channel.from([
              file(params.input.genome.chrom_sizes_path)
                 ])

    // LIBRARY_GROUPS:
    LIBRARY_GROUPS = Channel.from(
              params.input.library_groups.collect{ k, v -> [k, v] }
              )

    // the Channel the location of Raw Data (fastqs):
    RUN_SOURCES = Channel.from(
            params.input.raw_reads_paths.collect{
                k, v -> v.collect{k2, v2 -> (v2.size() == 1)
                                            ? [k,k2]+v2+[null]
                                            : [k,k2]+v2
                                 }
                }.sum()
        ).branch {
            download_truncate:
                        ( needsDownloading(it[2]) || needsDownloading(it[3]) )
            local_truncate:
                        ((params['map'].get('chunksize', 0) > 0) || (params['input'].get('truncate_fastq_reads', 0) > 0))
            local_no_processing:
                        true
        }

    /* Download or truncate input data */
    fastq_download_truncated = download_truncate_chunk_fastqs( RUN_SOURCES.download_truncate ).output
    fastq_local_truncated = local_truncate_chunk_fastqs( RUN_SOURCES.local_truncate ).output

    FASTQ = fastq_download_truncated
        .mix(fastq_local_truncated)
        .transpose()
        .map{[it[0],
              it[1],
              checkLeftRightChunk(it[2], it[3]), // index of the chunk (safety check)
              it[2],
              it[3]]}
        .mix(RUN_SOURCES.local_no_processing.transpose())

    /* Optional quality control */
    if (params.getOrDefault('do_fastqc', 'false').toBoolean()){
        fastqc( FASTQ.map{ it -> [ it[0], it[1], it[2], [1, 2], [it[3], it[4]] ] }.transpose() )
    }

    /* Mapping */
    BWA_INDEX = Channel.from([[
             params.input.genome.bwa_index_wildcard_path
                .split('/')[-1]
                .replaceAll('\\*$', "")
                .replaceAll('\\.$', ""),
             file(params.input.genome.bwa_index_wildcard_path),
            ]])

    INPUT_MAPPING = FASTQ
        .combine(BWA_INDEX)
        .combine(CHROM_SIZES)
        .multiMap{ it ->
            fastq: it[0..4]
            bwa_index: it[5..6]
            chrom_sizes: it[7]
        }
    BAM = map_parse_sort_chunks( INPUT_MAPPING ).output

    /* Merge .pairsams into libraries */
    INPUT_DEDUP = BAM
            .map{lib, run, chunk, pairsam, bam -> [lib, pairsam]}
            .groupTuple()
    PAIRS = merge_dedup_splitbam( INPUT_DEDUP )

//    /* Produce stats with filters: */
//    if (params.get('stats', [:]).get('use_filters', 'false').toBoolean()) {
//        INPUT_DEDUP_LIBRARIES = BAM
//                .map{lib, run, chunk, pairsam, bam -> [lib, pairsam]}
//        PAIRS_LIBRARIES = merge_dedup_splitbam( INPUT_DEDUP_LIBRARIES )
//    }

    /* Bin indexed .pairs into .cool matrices */
    FILTERS = Channel.from( params.bin.filters.collect{ name, expr -> [name, expr] } )
    INPUT_BIN = FILTERS
                .combine( PAIRS.output.map {v -> [v[0], v[1]]} )
                .combine(CHROM_SIZES)
                .multiMap{ it ->
                    filter: it[0..1]
                    pairs: it[2..3]
                    chrom_sizes: it[4]
                }
    COOLERS = bin_zoom_library_pairs( INPUT_BIN ).output

    /* Merge coolers by groups */
    INPUT_COOLERS_MERGE = LIBRARY_GROUPS.combine(COOLERS)
        .filter{ it[1].contains(it[2]) }
        .map {library_group, libraries, lib, filter_name, single_res_clr, multires_clr ->
             [library_group, filter_name, single_res_clr]}
        .groupTuple(by: [0, 1])
        .multiMap{it ->
            filter: it[1] //TODO: pass filter expression
            coolers:  [it[0], it[2]]
        }

    COOLERS_MERGED = merge_zoom_library_group_coolers( INPUT_COOLERS_MERGE ).output

    /* Merge .stats for library groups */
    INPUT_STATS_MERGE = LIBRARY_GROUPS
        .combine( PAIRS.stats )
        .filter{ it[1].contains(it[2]) }
        .map {library_group, libraries, lib, stats -> [library_group, stats]}
        .groupTuple()

    STATS_MERGED = merge_stats_libraries_into_groups( INPUT_STATS_MERGE ).output

    /* MultiQC report */
    REPORTS = MULTIQC( STATS_MERGED.map{it[1]} )

}

workflow {

    DISTILLER ( )

}
workflow.onComplete {
    log.info "Distilled!"
}