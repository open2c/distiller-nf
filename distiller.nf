#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

/*
 * Miscellaneous code for the pipeline
 */

switch(params.compression_format) {
    case 'gz':
        suffix = 'gz'
        decompress_command = 'bgzip -cd -@ 3'
        break
    case 'lz4':
        suffix = 'lz4'
        decompress_command = 'lz4c -cd'
        break
    default:
        suffix = 'gz'
        decompress_command = 'bgzip -cd -@ 3'
        break
}

boolean isSingleFile(object) {    
    object instanceof Path  
}

String getOutDir(output_type) {
    new File(params.output.get('base_dir', ''),
             params.output.dirs.get(output_type, output_type)).getCanonicalPath()
}

String getIntermediateDir(intermediate_type) {
    new File(params.intermediates.get('base_dir', ''),
             params.intermediates.dirs.get(
                intermediate_type, intermediate_type)).getCanonicalPath()
}

Boolean needsDownloading(query) {
    return (
        (query instanceof String) && (
           query.startsWith('sra:') 
           || query.startsWith('http://')
           || query.startsWith('https://')
           || query.startsWith('ftp://') 
       )
   )
}

String checkLeftRightChunk(left_chunk_fname,right_chunk_fname) {
    // checks if the chunk index is the same 
    // both for left and right chunks and returns 
    // that chunk index:
    // left by design:  ${library}.${run}.*.1.fastq.gz
    // right by design: ${library}.${run}.*.2.fastq.gz
    left_chunk_idx  =  left_chunk_fname.toString().tokenize('.')[-4]
    right_chunk_idx = right_chunk_fname.toString().tokenize('.')[-4]
    leftness_of_chunk  =  left_chunk_fname.toString().tokenize('.')[-3]
    rightness_of_chunk = right_chunk_fname.toString().tokenize('.')[-3]
    // assertions (should never happen by design, by just in case):
    // sidedness should be different:
    assert leftness_of_chunk != rightness_of_chunk: ("Sidedness suffix of"
                                                    +"LEFT and RIGHT sides of"
                                                    +"fastq chunks should DIFFER!")
    assert left_chunk_idx == right_chunk_idx: ("Chunk index of"
                                            +"LEFT and RIGHT sides of"
                                            +"fastq chunks should be IDENTICAL!")
    // return Chunk index of fastq chunk:
    return left_chunk_idx
}

// CHROM_SIZES:
// we need 2 copies of this Channel
// for Parsing and Binning:
Channel.from([
          file(params.input.genome.chrom_sizes_path)
             ]).into{CHROM_SIZES_FOR_PARSING;
                     CHROM_SIZES_FOR_BINNING}


// LIBRARY_GROUPS:
// we need 2 copies of this Channel
// for Parsing and Binning:
Channel.from(
          params.input.library_groups.collect{ k, v -> [k, v] }
            ).into{LIBRARY_GROUPS_FOR_COOLER_MERGE;
                   LIBRARY_GROUPS_FOR_STATS_MERGE}


// the Channel the location of Raw Data (fastqs):
LIB_RUN_SOURCES = Channel.from(
    params.input.raw_reads_paths.collect{
        k, v -> v.collect{k2, v2 -> (v2.size() == 1)
                                    ? [k,k2]+v2+[null]
                                    : [k,k2]+v2
                         }
        }.sum()
    )


LIB_RUN_SOURCES_DOWNLOAD_TRUNCATE_CHUNK = Channel.create()
LIB_RUN_SOURCES_LOCAL_TRUNCATE_CHUNK = Channel.create()
LIB_RUN_SOURCES_LOCAL_NO_PROCESSING = Channel.create()
LIB_RUN_SOURCES.choice(
    LIB_RUN_SOURCES_DOWNLOAD_TRUNCATE_CHUNK, 
    LIB_RUN_SOURCES_LOCAL_TRUNCATE_CHUNK, 
    LIB_RUN_SOURCES_LOCAL_NO_PROCESSING) {
    a -> ( ( needsDownloading(a[2]) || needsDownloading(a[3]) ) 
            ? 0
            : ( (    (params['map'].get('chunksize', 0) > 0)
                  || (params['input'].get('truncate_fastq_reads', 0) > 0)
                ) ? 1 : 2
              )
         )
}

LIB_RUN_SOURCES_LOCAL_NO_PROCESSING
    .map{ v -> [v[0], v[1], file(v[2]), file(v[3])]}
    .set{ LIB_RUN_SOURCES_LOCAL_NO_PROCESSING }

LIB_RUN_SOURCES_LOCAL_TRUNCATE_CHUNK
    .map{ v -> [v[0], v[1], file(v[2]), file(v[3])] }
    .set{ LIB_RUN_SOURCES_LOCAL_TRUNCATE_CHUNK }


/*
 * Download and chunk fastqs.
 */


def fastqDumpCmd(file_or_srr, library, run, srr_start=0, srr_end=-1, threads=1) {
    srr_start_flag = (srr_start == 0) ? '' : (' --minSpotId ' + srr_start)
    srr_end_flag = (srr_end == -1) ? '' : (' --maxSpotId ' + srr_end)
    sed_fwd_filter = "\'2,\$s/^@\\(SRR.*length\\)/\\v\t\\1/\'"
    sed_rev_filter = "\"s/^.\\?\\t/@/\""
    cmd = """
        HOME=`readlink -e ./`
        fastq-dump ${file_or_srr} -Z --split-spot ${srr_start_flag} ${srr_end_flag} \
                       | sed ${sed_fwd_filter} \
                       | split -n r/2 -t\$'\\v' --numeric-suffixes=1 --suffix-length 1 \
                         --filter 'sed ${sed_rev_filter} | bgzip -c -@ ${threads} > \$FILE.fastq.gz' - \
                         ${library}.${run}.  """
    return cmd
}


def sraDownloadTruncateCmd(sra_query, library, run, truncate_fastq_reads=0, 
                           chunksize=0, threads=1) {
    srr = ( sra_query =~ /SRR\d+/ )[0]
    srrnum = srr.substring(3)

    srr_start = 0
    if ( sra_query.contains('start=') ) {
        srr_start = ( sra_query =~ /start=(\d+)/ )[0][1] 
    } 
        
    srr_end = -1
    if ( truncate_fastq_reads ) {
        srr_end = srr_start + truncate_fastq_reads
    } else if ( sra_query.contains('end=') ) {
        srr_end = (sra_query =~ /end=(\d+)/)[0][1]
    }

    if ((srr_start > 0) || (srr_end != -1)) {
        cmd = """
            ${fastqDumpCmd(srr, library, run, srr_start, srr_end)}
            if [ -d ./ncbi ]; then rm -Rf ./ncbi; fi
        """

    }
    else {
        cmd = """
            wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR${srrnum.take(3)}/${srr}/${srr}.sra -O ${srr}.sra
            ${fastqDumpCmd(srr+'.sra', library, run, 0, -1)}
            rm ${srr}.sra
        """
    }

    chunk_lines = 4 * chunksize
    if ( (truncate_fastq_reads == 0) && (chunk_lines > 0) ) {
        for (side in 1..2) {
            cmd += """
                zcat ${library}.${run}.${side}.fastq.gz | \
                     split -l ${chunk_lines} --numeric-suffixes=1 \
                     --filter 'bgzip -c -@ ${threads} > \$FILE.${side}.fastq.gz' - \
                     ${library}.${run}.
                rm ${library}.${run}.${side}.fastq.gz 
            """
        }
    } else {
        for (side in 1..2) {
            cmd += """
                mv ${library}.${run}.${side}.fastq.gz ${library}.${run}.0.${side}.fastq.gz
            """
        }
    }

    return cmd
}


String fastqDownloadTruncateCmd(query, library, run, side, 
                                truncate_fastq_reads=0, chunksize=0, threads=1) {
    truncate_lines = 4 * truncate_fastq_reads
    chunk_lines = 4 * chunksize

    if (truncate_lines > 0) {
        cmd = """head -n ${truncate_lines} < <(wget ${query} -O - | gunzip -cd )\
                 | bgzip -c -@ ${threads} \
                 > ${library}.${run}.0.${side}.fastq.gz
              """
    } else if (chunk_lines > 0) {
        cmd = """wget ${query} -O - \
                 | gunzip -cd \
                 | split -l ${chunk_lines} --numeric-suffixes=1 \
                 --filter 'bgzip -c -@ ${threads} > \$FILE.${side}.fastq.gz' - \
                 ${library}.${run}.
             """
    } else {
        cmd = "wget ${query} -O ${library}.${run}.0.${side}.fastq.gz"
    }

    return cmd
}


String fastqLocalTruncateChunkCmd(path, library, run, side, 
                                  truncate_fastq_reads=0, chunksize=0, threads=1) {
    truncate_lines = 4 * truncate_fastq_reads
    chunk_lines = 4 * chunksize

    if (truncate_lines > 0) {
        cmd = """head -n ${truncate_lines} < <( zcat ${path} ) \
                 | bgzip -c -@ ${threads} \
                 >  ${library}.${run}.0.${side}.fastq.gz 
        """
    } else if (chunk_lines > 0) {
        cmd = """
            zcat ${path} | \
            split -l ${chunk_lines} --numeric-suffixes=1 \
            --filter 'bgzip -c -@ ${threads} > \$FILE.${side}.fastq.gz' - \
            ${library}.${run}.
        """
    } else {
        // this line should never be reached, b/c local files should only
        // be processed if truncation or chunking is requested.
        cmd = "mv ${path} ${library}.${run}.0.${side}.fastq.gz"
    }

    return cmd
}


process download_truncate_chunk_fastqs{
    tag "library:${library} run:${run}"
    storeDir getIntermediateDir('processed_fastqs')

    input:
    set val(library), val(run), 
        val(query1), val(query2) from LIB_RUN_SOURCES_DOWNLOAD_TRUNCATE_CHUNK
     
    output:
    set library, run, 
        "${library}.${run}.*.1.fastq.gz", 
        "${library}.${run}.*.2.fastq.gz" into LIB_RUN_CHUNK_DOWNLOADED_PROCESSED
 
    script:

    truncate_fastq_reads = params['input'].get('truncate_fastq_reads',0)
    chunksize = params['map'].get('chunksize', 0) 

    if (query1.startsWith('sra:')) {
        if ( !(( query2 == null) || (! query2.toBoolean())) ) {
            error "Runs defined with SRA should only contain one line"
        }
        
        download_truncate_chunk_cmd1 = sraDownloadTruncateCmd(
            query1, library, run, truncate_fastq_reads, chunksize, task.cpus)

        download_truncate_chunk_cmd2 = ""
    } else {
        download_truncate_chunk_cmd1 = fastqDownloadTruncateCmd(
            query1, library, run, 1, truncate_fastq_reads, chunksize, task.cpus)

        download_truncate_chunk_cmd2 = fastqDownloadTruncateCmd(
            query2, library, run, 2, truncate_fastq_reads, chunksize, task.cpus)
    }


    """
    ${download_truncate_chunk_cmd1}
    ${download_truncate_chunk_cmd2}
    """
}

process local_truncate_chunk_fastqs{
    tag "library:${library} run:${run}"
    storeDir getIntermediateDir('processed_fastqs')

    input:
    set val(library), val(run), 
        file(fastq1), file(fastq2) from LIB_RUN_SOURCES_LOCAL_TRUNCATE_CHUNK
     
    output:
    set library, run, 
        "${library}.${run}.*.1.fastq.gz", 
        "${library}.${run}.*.2.fastq.gz" into LIB_RUN_CHUNK_LOCAL_PROCESSED
 
    script:

    truncate_fastq_reads = params['input'].get('truncate_fastq_reads',0)
    chunksize = params['map'].get('chunksize', 0) 

    truncate_chunk_cmd1 = fastqLocalTruncateChunkCmd(
        fastq1, library, run, 1, truncate_fastq_reads, chunksize, task.cpus)
    truncate_chunk_cmd2 = fastqLocalTruncateChunkCmd(
        fastq2, library, run, 2, truncate_fastq_reads, chunksize, task.cpus)

    """
    ${truncate_chunk_cmd1}
    ${truncate_chunk_cmd2}
    """
}


// use new transpose operator 
// to undo 'groupBy' of 'chunk_fastqs' process:
// https://github.com/nextflow-io/nextflow/issues/440
LIB_RUN_CHUNK_DOWNLOADED_PROCESSED
    .transpose()
    .map{[it[0],
          it[1],
          // index of the chunk (checked for safety):
          checkLeftRightChunk(it[2],it[3]),
          it[2],
          it[3]]}
    .set{ LIB_RUN_CHUNK_FASTQS }

LIB_RUN_CHUNK_LOCAL_PROCESSED
    .transpose()
    .map{[it[0],
          it[1],
          // index of the chunk (checked for safety):
          checkLeftRightChunk(it[2],it[3]),
          it[2],
          it[3]]}
    .mix(LIB_RUN_CHUNK_FASTQS)
    .set{ LIB_RUN_CHUNK_FASTQS }


LIB_RUN_SOURCES_LOCAL_NO_PROCESSING
    .map{[it[0],
          it[1],
          // index of the non-chunked is 0:
          0,
          it[2],
          it[3]]}
    .mix(LIB_RUN_CHUNK_FASTQS)
    .set{LIB_RUN_CHUNK_FASTQS}


/*
 * FastQC the input files.
 */

LIB_RUN_CHUNK_FASTQS_FOR_QC = Channel.create()
LIB_RUN_CHUNK_FASTQS
    .tap(LIB_RUN_CHUNK_FASTQS_FOR_QC)
    .set{LIB_RUN_CHUNK_FASTQS}

LIB_RUN_CHUNK_FASTQS_FOR_QC
    .filter { it -> params.get('do_fastqc', 'false').toBoolean() }
    .map{ v -> [v[0], v[1], v[2], [[1,file(v[3])], [2,file(v[4])]]]} 
    .flatMap{ 
        vs -> vs[3].collect{ 
            it -> [vs[0],
                   vs[1], 
                   vs[2], 
                   it[0],
                   it[1]] } }
    .set {LIB_RUN_CHUNK_SIDE_FASTQS_FOR_QC}

process fastqc{

    tag "library:${library} run:${run} chunk:${chunk} side:${side}"
    publishDir path: getOutDir('fastqc'), mode:"copy"

    input:
    set val(library), val(run), val(chunk), val(side), 
        file(fastq) from LIB_RUN_CHUNK_SIDE_FASTQS_FOR_QC

    output:
    set library, run, chunk, side,  
        "${library}.${run}.${chunk}.${side}_fastqc.html", 
        "${library}.${run}.${chunk}.${side}_fastqc.zip" into LIB_RUN_CHUNK_SIDE_QCS

    """
    mkdir -p ./temp_fastqc/
    ln -s \"\$(readlink -f ${fastq})\" ./temp_fastqc/${library}.${run}.${chunk}.${side}.fastq.gz
    fastqc --threads ${task.cpus} -o ./ -f fastq ./temp_fastqc/${library}.${run}.${chunk}.${side}.fastq.gz
    rm -r ./temp_fastqc/
    """
          
}



BWA_INDEX = Channel.from([[
             params.input.genome.bwa_index_wildcard
                .split('/')[-1]
                .replaceAll('\\*$', "")
                .replaceAll('\\.$', ""),
             file(params.input.genome.bwa_index_wildcard),
            ]])

/*
 * Map fastq files
 */
process map_parse_sort_chunks {
    tag "library:${library} run:${run} chunk:${chunk}"
    storeDir getIntermediateDir('mapped_parsed_sorted_chunks')
 
    input:
    set val(library), val(run), val(chunk), file(fastq1), file(fastq2) from LIB_RUN_CHUNK_FASTQS
    set val(bwa_index_base), file(bwa_index_files) from BWA_INDEX.first()
    file(chrom_sizes) from CHROM_SIZES_FOR_PARSING.first()
     
    output:
    set library, run, chunk,
        "${library}.${run}.${chunk}.pairsam.${suffix}",
        "${library}.${run}.${chunk}.bam" into LIB_RUN_CHUNK_PAIRSAMS

    script:
    // additional mapping options or empty-line
    mapping_options = params['map'].get('mapping_options','')
    dropsam_flag = params['map'].get('drop_sam','false').toBoolean() ? '--drop-sam' : ''
    dropreadid_flag = params['map'].get('drop_readid','false').toBoolean() ? '--drop-readid' : ''
    dropseq_flag = params['map'].get('drop_seq','false').toBoolean() ? '--drop-seq' : ''
    keep_unparsed_bams_command = ( 
        params['map'].get('keep_unparsed_bams','false').toBoolean() ? 
        "| tee >(samtools view -bS > ${library}.${run}.${chunk}.bam)" : "" )

    """
    mkdir ./tmp4sort
    touch ${library}.${run}.${chunk}.bam 
    bwa mem -t ${task.cpus} ${mapping_options} -SP ${bwa_index_base} ${fastq1} ${fastq2} \
        ${keep_unparsed_bams_command} \
        | pairtools parse ${dropsam_flag} ${dropreadid_flag} ${dropseq_flag} \
            --add-columns mapq \
            -c ${chrom_sizes} \
            | pairtools sort --nproc ${task.cpus} \
                             -o ${library}.${run}.${chunk}.pairsam.${suffix} \
                             --tmpdir ./tmp4sort \
            | cat

    rm -rf ./tmp4sort

    """        

}

/*
 * Merge .pairsams into libraries
 */

LIB_RUN_CHUNK_PAIRSAMS
    .map {library, run, chunk, pairsam, bam -> tuple(library, pairsam)}
    .groupTuple()
    .set {LIB_PAIRSAMS_TO_MERGE}

process merge_dedup_splitbam {
    tag "library:${library}"
    storeDir getIntermediateDir('pairs_library')
 
    input:
    set val(library), file(run_pairsam) from LIB_PAIRSAMS_TO_MERGE
     
    output:
    set library, "${library}.nodups.pairs.gz", 
                 "${library}.nodups.pairs.gz.px2", 
                 "${library}.nodups.bam",
                 "${library}.dups.pairs.gz", "${library}.dups.bam", 
                 "${library}.unmapped.pairs.gz", 
                 "${library}.unmapped.bam" into LIB_PAIRS_BAMS
    set library, "${library}.dedup.stats" into LIB_DEDUP_STATS
 
    script:
    dropsam = params['map'].get('drop_sam','false').toBoolean()
    merge_command = ( 
        isSingleFile(run_pairsam) ?
        "${decompress_command} ${run_pairsam}" : 
        "pairtools merge ${run_pairsam} --nproc ${task.cpus} --tmpdir ./tmp4sort"
    )

    if(dropsam) 
        """
        mkdir ./tmp4sort

        ${merge_command} | pairtools dedup \
            --max-mismatch ${params.filter.pcr_dups_max_mismatch_bp} \
            --mark-dups \
            --output ${library}.nodups.pairs.gz \
            --output-unmapped ${library}.unmapped.pairs.gz \
            --output-dups ${library}.dups.pairs.gz \
            --output-stats ${library}.dedup.stats \
            | cat

        touch ${library}.unmapped.bam
        touch ${library}.nodups.bam
        touch ${library}.dups.bam

        rm -rf ./tmp4sort
        pairix ${library}.nodups.pairs.gz
        """
    else 
        """
        mkdir ./tmp4sort

        ${merge_command} | pairtools dedup \
            --max-mismatch ${params.filter.pcr_dups_max_mismatch_bp} \
            --mark-dups \
            --output \
                >( pairtools split \
                    --output-pairs ${library}.nodups.pairs.gz \
                    --output-sam ${library}.nodups.bam \
                 ) \
            --output-unmapped \
                >( pairtools split \
                    --output-pairs ${library}.unmapped.pairs.gz \
                    --output-sam ${library}.unmapped.bam \
                 ) \
            --output-dups \
                >( pairtools split \
                    --output-pairs ${library}.dups.pairs.gz \
                    --output-sam ${library}.dups.bam \
                 ) \
            --output-stats ${library}.dedup.stats \
            | cat

        rm -rf ./tmp4sort
        pairix ${library}.nodups.pairs.gz
        """
}

LIB_PAIRS_BAMS
    .map {v -> tuple(v[0], v[1], v[2])}
    .set {LIB_IDX_PAIRS}

/*
 * Bin indexed .pairs into .cool matrices.
 */ 

process bin_library_pairs{
    tag "library:${library} resolution:${res}"
    publishDir path: getOutDir('coolers_library'), mode:"copy", saveAs: {"${library}.${res}.cool"}


    input:
        set val(library), file(pairs_lib), file(pairs_index_lib) from LIB_IDX_PAIRS
        each res from params['bin'].resolutions
        file(chrom_sizes) from CHROM_SIZES_FOR_BINNING.first()

    output:
        set library, res, "${library}.${res}.cool" into LIB_RES_COOLERS, LIB_RES_COOLERS_TO_ZOOM

    script:

    // get any additional balancing options, if provided
    balance_options = params['bin'].get('balance_options','')
    // balancing command if it's requested
    balance_command = ( params['bin'].get('balance','false').toBoolean() ? 
        "cooler balance --nproc ${task.cpus} ${balance_options} ${library}.${res}.cool" : "" )

    """
    cooler cload pairix \
        --nproc ${task.cpus} \
        --assembly ${params.input.genome.assembly} \
        ${chrom_sizes}:${res} ${pairs_lib} ${library}.${res}.cool

    ${balance_command}
    """
}


/*
 * Zoomify .cool matrices with highest resolution for libraries (when requested).
 */ 

// use library-cooler file with the highest resolution (smallest bin size) to zoomify:
LIB_RES_COOLERS_TO_ZOOM
    .map{ library,res,cool -> [library,[res,cool]] }
    .groupTuple( by: 0, sort: {res,cool -> res} )
    // after grouping by library get res_cool_list with smallest res (highest resoution)
    .map{ library,res_cool_list -> [library,res_cool_list[0]].flatten() }
    .set{LIB_RES_COOLERS_TO_ZOOM}


process zoom_library_coolers{
    tag "library:${library} zoom"
    publishDir path: getOutDir('zoom_coolers_library'), mode:"copy", saveAs: {"${library}.${res}.multires.cool"}


    input:
        set val(library), val(res), file(cool) from LIB_RES_COOLERS_TO_ZOOM

    output:
        set library, res, "${library}.${res}.multires.cool" into LIB_RES_COOLERS_ZOOMED

    // run ot only if requested
    when:
    params['bin'].get('zoomify','false').toBoolean()


    script:

    // additional balancing options as '--balance-args' or empty-line
    balance_options = params['bin'].get('balance_options','')
    balance_options = ( balance_options ? "--balance-args \"${balance_options}\"": "")
    // balancing flag if it's requested
    balance_flag = ( params['bin'].get('balance','false').toBoolean() ? "--balance ${balance_options}" : "--no-balance" )

    """
    cooler zoomify \
        --nproc ${task.cpus} \
        --out ${library}.${res}.multires.cool \
        ${balance_flag} \
        ${cool}
    """
}



/*
 * Merge .cool matrices for library groups.
 */ 


LIBRARY_GROUPS_FOR_COOLER_MERGE
    .combine(LIB_RES_COOLERS)
    .filter{ it[1].contains(it[2]) } 
    .map {library_group, libraries, library, res, file -> tuple(library_group, res, file)}
    .groupTuple(by: [0, 1])
    .set { LIBGROUP_RES_COOLERS_TO_MERGE }

process make_library_group_coolers{
    tag "library_group:${library_group} resolution:${res}"
    publishDir path: getOutDir('coolers_library_group'), mode:"copy", saveAs: {"${library_group}.${res}.cool"}

    input:
        set val(library_group), val(res), file(coolers) from LIBGROUP_RES_COOLERS_TO_MERGE

    output:
        set library_group, res, "${library_group}.${res}.cool" into LIBGROUP_RES_COOLERS, LIBGROUP_RES_COOLERS_TO_ZOOM

    script:

    // get any additional balancing options, if provided
    balance_options = params['bin'].get('balance_options','')
    // balancing command if it's requested
    balance_command = ( params['bin'].get('balance','false').toBoolean() ? 
        "cooler balance --nproc ${task.cpus} ${balance_options} ${library_group}.${res}.cool" : "" )

    if( isSingleFile(coolers))
        // .cool is already balanced in such case:
        """
        ln -s \$(readlink -f ${coolers}) ${library_group}.${res}.cool
        """
    else
        // 'weight' column is gone after merging, so balance if requested:
        """
        cooler merge ${library_group}.${res}.cool ${coolers}

        ${balance_command}
        """
}



/*
 * Zoomify .cool matrices with highest resolution for library groups (when requested).
 */ 

// use library-group-cooler file with the highest resolution (smallest bin size) to zoomify:
LIBGROUP_RES_COOLERS_TO_ZOOM
    .map{ library_group,res,cool -> [library_group,[res,cool]] }
    .groupTuple( by: 0, sort: {res,cool -> res})
    // after grouping by library_group get res_cool_list with smallest bin (highest resoution)
    .map{ library_group,res_cool_list -> [library_group,res_cool_list[0]].flatten() }
    .set{LIBGROUP_RES_COOLERS_TO_ZOOM}


process zoom_library_group_coolers{
    tag "library_group:${library_group} zoom"
    publishDir path: getOutDir('zoom_coolers_library_group'), mode:"copy", saveAs: {"${library_group}.${res}.multires.cool"}


    input:
        set val(library_group), val(res), file(cool) from LIBGROUP_RES_COOLERS_TO_ZOOM

    output:
        set library_group, res, "${library_group}.${res}.multires.cool" into LIBGROUP_RES_COOLERS_ZOOMED

    // run ot only if requested
    when:
    params['bin'].get('zoomify','false').toBoolean()


    script:

    // additional balancing options as '--balance-args' or empty-line
    balance_options = params['bin'].get('balance_options','')
    balance_options = ( balance_options ? "--balance-args \"${balance_options}\"": "")
    // balancing flag if it's requested
    balance_flag = ( params['bin'].get('balance','false').toBoolean() ? "--balance ${balance_options}" : "--no-balance" )

    """
    cooler zoomify \
        --nproc ${task.cpus} \
        --out ${library_group}.${res}.multires.cool \
        ${balance_flag} \
        ${cool}
    """
}




/*
 * Merge .stats for library groups
 */ 



LIBRARY_GROUPS_FOR_STATS_MERGE
    .combine(LIB_DEDUP_STATS)
    .filter{ it[1].contains(it[2]) } 
    .map {library_group, libraries, library, stats -> tuple(library_group, stats)}
    .groupTuple()
    .set { LIBGROUP_STATS_TO_MERGE }


process merge_stats_libraries_into_groups {
    tag "library_group:${library_group}"
    publishDir path: getOutDir('stats_library_group'), pattern: "*.stats", mode:"copy"
 
    input:
    set val(library_group), file(stats) from LIBGROUP_STATS_TO_MERGE
     
    output:
    set library_group, "${library_group}.stats" into LIBGROUP_STATS

    script:
    if( isSingleFile(stats))
        """
        ln -s ${stats} ${library_group}.stats
        """
    else
        """
        pairtools stats --merge ${stats} -o ${library_group}.stats
        """
}


