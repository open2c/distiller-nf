#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

/*
 * Miscellaneous code for the pipeline
 */

MIN_RES = params['bin'].resolutions.collect { it as int }.min()
ASSEMBLY_NAME = params['input'].genome.assembly_name

pairsgz_decompress_command = 'bgzip -cd -@ 3'

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

String getOutputDir(output_type) {
    new File(params.output.dirs.get(output_type, output_type)).getCanonicalPath()
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
    def srr_start_flag = (srr_start == 0) ? '' : (' --minSpotId ' + srr_start)
    def srr_end_flag = (srr_end == -1) ? '' : (' --maxSpotId ' + srr_end)
    def bgzip_threads = Math.max(1,((threads as int)-2).intdiv(2))

    def cmd = """
        HOME=`readlink -e ./`
        fastq-dump ${file_or_srr} -Z --split-spot ${srr_start_flag} ${srr_end_flag} \
                       | pyfilesplit --lines 4 \
                         >(bgzip -c -@{bgzip_threads} > ${library}.${run}.1.fastq.gz) \
                         >(bgzip -c -@{bgzip_threads} > ${library}.${run}.2.fastq.gz) \
                         | cat """

    return cmd
}


def sraDownloadTruncateCmd(sra_query, library, run, truncate_fastq_reads=0,
                           chunksize=0, threads=1) {
    def cmd = ""

    def srr = ( sra_query =~ /SRR\d+/ )[0]
    def srrnum = srr.substring(3)

    def srr_start = 0
    def srr_end = -1

    if ( sra_query.contains('start=') ) {
        srr_start = ( sra_query =~ /start=(\d+)/ )[0][1]
    }

    if ( truncate_fastq_reads ) {
        srr_end = srr_start + truncate_fastq_reads
    } else if ( sra_query.contains('end=') ) {
        srr_end = (sra_query =~ /end=(\d+)/)[0][1]
    }

    wget_url = "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR${srrnum.take(3)}/${srr}/${srr}.sra"
    if ((srr_start > 0) || (srr_end != -1)) {
        cmd = """
            ${fastqDumpCmd(srr, library, run, srr_start, srr_end, threads)}
            if [ -d ./ncbi ]; then rm -Rf ./ncbi; fi
        """
    }
    else {
        cmd = """
            if wget --spider ${wget_url} 2>/dev/null; then
                wget ${wget_url} -qO ${srr}.sra
                ${fastqDumpCmd(srr+'.sra', library, run, 0, -1, threads)}
                rm ${srr}.sra
            else
                echo 'Cannot wget an sra, fall back to fastq-dump'
                ${fastqDumpCmd(srr, library, run, 0, -1, threads)}
                if [ -d ./ncbi ]; then rm -Rf ./ncbi; fi
            fi
        """
    }

    def chunk_lines = 4 * chunksize
    def split_bgzip_threads = Math.max(1, (threads as int)-1)
    if ( (truncate_fastq_reads == 0) && (chunk_lines > 0) ) {
        for (side in 1..2) {
            cmd += """
                zcat ${library}.${run}.${side}.fastq.gz | \
                     split -l ${chunk_lines} --numeric-suffixes=1 \
                     --filter 'bgzip -c -@ ${split_bgzip_threads} > \$FILE.${side}.fastq.gz' - \
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
                                truncate_fastq_reads=0, chunksize=0, threads=2) {
    def cmd = ''

    def truncate_lines = 4 * truncate_fastq_reads
    def chunk_lines = 4 * chunksize
    def bgzip_threads = Math.max(1, (threads as int)-1)

    if (truncate_lines > 0) {
        cmd = """head -n ${truncate_lines} < <(wget ${query} -O - | gunzip -cd )\
                 | bgzip -c -@ ${bgzip_threads} \
                 > ${library}.${run}.0.${side}.fastq.gz
              """
    } else if (chunk_lines > 0) {
        cmd = """wget ${query} -O - \
                 | gunzip -cd \
                 | split -l ${chunk_lines} --numeric-suffixes=1 \
                 --filter 'bgzip -c -@ ${bgzip_threads} > \$FILE.${side}.fastq.gz' - \
                 ${library}.${run}.
             """
    } else {
        cmd = "wget ${query} -O ${library}.${run}.0.${side}.fastq.gz"
    }

    return cmd
}


String fastqLocalTruncateChunkCmd(path, library, run, side,
                                  truncate_fastq_reads=0, chunksize=0, threads=1) {
    def cmd = ""

    def truncate_lines = 4 * truncate_fastq_reads
    def chunk_lines = 4 * chunksize
    def bgzip_threads = Math.max(1, (threads as int)-1)

    if (truncate_lines > 0) {
        cmd = """head -n ${truncate_lines} < <( zcat ${path} ) \
                 | bgzip -c -@ ${bgzip_threads} \
                 >  ${library}.${run}.0.${side}.fastq.gz
        """
    } else if (chunk_lines > 0) {
        cmd = """
            zcat ${path} | \
            split -l ${chunk_lines} --numeric-suffixes=1 \
            --filter 'bgzip -c -@ ${bgzip_threads} > \$FILE.${side}.fastq.gz' - \
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
    storeDir getOutputDir('processed_fastqs')

    input:
    set val(library), val(run),
        val(query1), val(query2) from LIB_RUN_SOURCES_DOWNLOAD_TRUNCATE_CHUNK

    output:
    set library, run,
        "${library}.${run}.*.1.fastq.gz",
        "${library}.${run}.*.2.fastq.gz" into LIB_RUN_CHUNK_DOWNLOADED_PROCESSED

    script:
    def truncate_fastq_reads = params['input'].get('truncate_fastq_reads',0)
    def chunksize = params['map'].get('chunksize', 0)

    def download_truncate_chunk_cmd1 = ""
    def download_truncate_chunk_cmd2 = ""

    if (query1.startsWith('sra:')) {
        if ( !(( query2 == null) || (! query2.toBoolean())) ) {
            error "Runs defined with SRA should only contain one line"
        }

        download_truncate_chunk_cmd1 += sraDownloadTruncateCmd(
            query1, library, run, truncate_fastq_reads, chunksize, task.cpus)

    } else {
        download_truncate_chunk_cmd1 += fastqDownloadTruncateCmd(
            query1, library, run, 1, truncate_fastq_reads, chunksize, task.cpus)

        download_truncate_chunk_cmd2 += fastqDownloadTruncateCmd(
            query2, library, run, 2, truncate_fastq_reads, chunksize, task.cpus)
    }


    """
    ${download_truncate_chunk_cmd1}
    ${download_truncate_chunk_cmd2}
    """
}

process local_truncate_chunk_fastqs{
    tag "library:${library} run:${run}"
    storeDir getOutputDir('processed_fastqs')

    input:
    set val(library), val(run),
        file(fastq1), file(fastq2) from LIB_RUN_SOURCES_LOCAL_TRUNCATE_CHUNK

    output:
    set library, run,
        "${library}.${run}.*.1.fastq.gz",
        "${library}.${run}.*.2.fastq.gz" into LIB_RUN_CHUNK_LOCAL_PROCESSED

    script:

    def truncate_fastq_reads = params['input'].get('truncate_fastq_reads',0)
    def chunksize = params['map'].get('chunksize', 0)

    def truncate_chunk_cmd1 = fastqLocalTruncateChunkCmd(
        fastq1, library, run, 1, truncate_fastq_reads, chunksize, task.cpus)
    def truncate_chunk_cmd2 = fastqLocalTruncateChunkCmd(
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
    storeDir getOutputDir('fastqc')

    input:
    set val(library), val(run), val(chunk), val(side),
        file(fastq) from LIB_RUN_CHUNK_SIDE_FASTQS_FOR_QC

    output:
    set library, run, chunk, side,
        "${library}.${run}.${chunk}.${side}_fastqc.html",
        "${library}.${run}.${chunk}.${side}_fastqc.zip" into LIB_RUN_CHUNK_SIDE_QCS

    """
    TASK_TMP_DIR=\$(mktemp -d -p ${task.distillerTmpDir} distiller.tmp.XXXXXXXXXX)
    ln -s \"\$(readlink -f ${fastq})\" \$TASK_TMP_DIR/${library}.${run}.${chunk}.${side}.fastq.gz
    fastqc --threads ${task.cpus} -o ./ -f fastq \$TASK_TMP_DIR/${library}.${run}.${chunk}.${side}.fastq.gz
    rm -r \$TASK_TMP_DIR
    """

}



BWA_INDEX = Channel.from([[
             params.input.genome.bwa_index_wildcard_path
                .split('/')[-1]
                .replaceAll('\\*$', "")
                .replaceAll('\\.$', ""),
             file(params.input.genome.bwa_index_wildcard_path),
            ]])

/*
 * Map fastq files
 */
process map_parse_sort_chunks {
    tag "library:${library} run:${run} chunk:${chunk}"
    storeDir getOutputDir('mapped_parsed_sorted_chunks')

    input:
    set val(library), val(run), val(chunk), file(fastq1), file(fastq2) from LIB_RUN_CHUNK_FASTQS
    set val(bwa_index_base), file(bwa_index_files) from BWA_INDEX.first()
    file(chrom_sizes) from CHROM_SIZES_FOR_PARSING.first()

    output:
    set library, run, chunk,
        "${library}.${run}.${ASSEMBLY_NAME}.${chunk}.pairsam.${suffix}",
        "${library}.${run}.${ASSEMBLY_NAME}.${chunk}.bam" into LIB_RUN_CHUNK_PAIRSAMS

    script:
    // additional mapping options or empty-line
    def mapping_options = params['map'].get('mapping_options','')
    def trim_options = params['map'].get('trim_options','')

    def dropsam_flag = params['parse'].get('make_pairsam','false').toBoolean() ? '' : '--drop-sam'
    def dropreadid_flag = params['parse'].get('drop_readid','false').toBoolean() ? '--drop-readid' : ''
    def dropseq_flag = params['parse'].get('drop_seq','false').toBoolean() ? '--drop-seq' : ''
    def keep_unparsed_bams_command = (
        params['parse'].get('keep_unparsed_bams','false').toBoolean() ?
        "| tee >(samtools view -bS > ${library}.${run}.${ASSEMBLY_NAME}.${chunk}.bam)" : "" )
    def parsing_options = params['parse'].get('parsing_options','')
    def bwa_threads = (task.cpus as int)
    def sorting_threads = (task.cpus as int)

    def mapping_command = (
        params['map'].get('trim_options','').toBoolean() ?
        "bwa mem \
        -t ${bwa_threads} \
        ${mapping_options} \
        -SP ${bwa_index_base} \
        ${fastq1} ${fastq2} \
        ${keep_unparsed_bams_command}" : " \
        fastp ${trim_options} -i ${fastq1} -I ${fastq2} --stdout | \
        bwa mem \
        -p \
        -t ${sorting_threads} \
        ${mapping_options} \
        -SP ${bwa_index_base} \
        - ${keep_unparsed_bams_command}")


    """
    TASK_TMP_DIR=\$(mktemp -d -p ${task.distillerTmpDir} distiller.tmp.XXXXXXXXXX)
    touch ${library}.${run}.${ASSEMBLY_NAME}.${chunk}.bam

    ${mapping_command} \
    | pairtools parse ${dropsam_flag} ${dropreadid_flag} ${dropseq_flag} \
      ${parsing_options} \
      -c ${chrom_sizes} \
      | pairtools sort --nproc task.cpus \
                     -o ${library}.${run}.${ASSEMBLY_NAME}.${chunk}.pairsam.${suffix} \
                     --tmpdir \$TASK_TMP_DIR \
      | cat

    rm -rf \$TASK_TMP_DIR
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
    storeDir getOutputDir('pairs_library')

    input:
    set val(library), file(run_pairsam) from LIB_PAIRSAMS_TO_MERGE

    output:
    set library, "${library}.${ASSEMBLY_NAME}.nodups.pairs.gz",
                 "${library}.${ASSEMBLY_NAME}.nodups.pairs.gz.px2",
                 "${library}.${ASSEMBLY_NAME}.nodups.bam",
                 "${library}.${ASSEMBLY_NAME}.dups.pairs.gz",
                 "${library}.${ASSEMBLY_NAME}.dups.bam",
                 "${library}.${ASSEMBLY_NAME}.unmapped.pairs.gz",
                 "${library}.${ASSEMBLY_NAME}.unmapped.bam" into LIB_PAIRS_BAMS
    set library, "${library}.${ASSEMBLY_NAME}.dedup.stats" into LIB_DEDUP_STATS

    script:
    def make_pairsam = params['parse'].get('make_pairsam','false').toBoolean()
    def merge_command = (
        isSingleFile(run_pairsam) ?
        "${decompress_command} ${run_pairsam}" :
        "pairtools merge ${run_pairsam} --nproc ${task.cpus} --tmpdir \$TASK_TMP_DIR"
    )

    if(make_pairsam)
        """
        TASK_TMP_DIR=\$(mktemp -d -p ${task.distillerTmpDir} distiller.tmp.XXXXXXXXXX)

        ${merge_command} | pairtools dedup \
            --max-mismatch ${params.dedup.max_mismatch_bp} \
            --mark-dups \
            --output \
                >( pairtools split \
                    --output-pairs ${library}.${ASSEMBLY_NAME}.nodups.pairs.gz \
                    --output-sam ${library}.${ASSEMBLY_NAME}.nodups.bam \
                 ) \
            --output-unmapped \
                >( pairtools split \
                    --output-pairs ${library}.${ASSEMBLY_NAME}.unmapped.pairs.gz \
                    --output-sam ${library}.${ASSEMBLY_NAME}.unmapped.bam \
                 ) \
            --output-dups \
                >( pairtools split \
                    --output-pairs ${library}.${ASSEMBLY_NAME}.dups.pairs.gz \
                    --output-sam ${library}.${ASSEMBLY_NAME}.dups.bam \
                 ) \
            --output-stats ${library}.${ASSEMBLY_NAME}.dedup.stats \
            | cat

        rm -rf \$TASK_TMP_DIR
        pairix ${library}.${ASSEMBLY_NAME}.nodups.pairs.gz
        """
    else
        """
        TASK_TMP_DIR=\$(mktemp -d -p ${task.distillerTmpDir} distiller.tmp.XXXXXXXXXX)

        ${merge_command} | pairtools dedup \
            --max-mismatch ${params.dedup.max_mismatch_bp} \
            --mark-dups \
            --output ${library}.${ASSEMBLY_NAME}.nodups.pairs.gz \
            --output-unmapped ${library}.${ASSEMBLY_NAME}.unmapped.pairs.gz \
            --output-dups ${library}.${ASSEMBLY_NAME}.dups.pairs.gz \
            --output-stats ${library}.${ASSEMBLY_NAME}.dedup.stats \
            | cat

        touch ${library}.${ASSEMBLY_NAME}.unmapped.bam
        touch ${library}.${ASSEMBLY_NAME}.nodups.bam
        touch ${library}.${ASSEMBLY_NAME}.dups.bam

        rm -rf \$TASK_TMP_DIR
        pairix ${library}.${ASSEMBLY_NAME}.nodups.pairs.gz
        """
}

LIB_PAIRS_BAMS
    .map {v -> tuple(v[0], v[1])}
    .set {LIB_PAIRS}
FILTERS = Channel.from(
    params.bin.filters.collect{ name, expr -> [name, expr] } )
FILTERS
    .combine(LIB_PAIRS)
    .set {LIB_FILTER_PAIRS}

/*
 * Bin indexed .pairs into .cool matrices.
 */

process bin_zoom_library_pairs{
    tag "library:${library} filter:${filter_name}"
    storeDir getOutputDir('coolers_library')

    input:
        set val(filter_name), val(filter_expr), val(library), file(pairs_lib) from LIB_FILTER_PAIRS
        file(chrom_sizes) from CHROM_SIZES_FOR_BINNING.first()

    output:
        set library, filter_name, "${library}.${ASSEMBLY_NAME}.${filter_name}.${MIN_RES}.cool",
            "${library}.${ASSEMBLY_NAME}.${filter_name}.${MIN_RES}.mcool" into LIB_FILTER_COOLERS_ZOOMED

    script:

    def res_str = params['bin'].resolutions.join(',')
    // get any additional balancing options, if provided
    def balance_options = params['bin'].get('balance_options','')
    balance_options = ( balance_options ? "--balance-args \"${balance_options}\"": "")
    // balancing flag if it's requested
    def balance_flag = ( params['bin'].get('balance','true').toBoolean() ? "--balance ${balance_options}" : " " )
    def filter_command = (filter_expr == '' ? '' : "| pairtools select '${filter_expr}'")

    """
    ${pairsgz_decompress_command} ${pairs_lib} ${filter_command} | cooler cload pairs \
        -c1 2 -p1 3 -c2 4 -p2 5 \
        --assembly ${ASSEMBLY_NAME} \
        ${chrom_sizes}:${MIN_RES} - ${library}.${ASSEMBLY_NAME}.${filter_name}.${MIN_RES}.cool

    cooler zoomify \
        --nproc ${task.cpus} \
        --out ${library}.${ASSEMBLY_NAME}.${filter_name}.${MIN_RES}.mcool \
        --resolutions ${res_str} \
        ${balance_flag} \
        ${library}.${ASSEMBLY_NAME}.${filter_name}.${MIN_RES}.cool

    """
}

/*
 * Merge .cool matrices for library groups.
 */

LIBRARY_GROUPS_FOR_COOLER_MERGE
    .combine(LIB_FILTER_COOLERS_ZOOMED)
    .filter{ it[1].contains(it[2]) }
    .map {library_group, libraries, library, filter_name, single_res_clr, multires_clr -> tuple(library_group, filter_name, single_res_clr)}
    .groupTuple(by: [0, 1])
    .set { LIBGROUP_FILTER_COOLERS_TO_MERGE }

process merge_zoom_library_group_coolers{
    tag "library_group:${library_group} filter:${filter_name}"
    publishDir path: getOutputDir('coolers_library_group'), mode: "copy"

    input:
        set val(library_group), val(filter_name), file(coolers) from LIBGROUP_FILTER_COOLERS_TO_MERGE

    output:
        set library_group, filter_name,
            "${library_group}.${ASSEMBLY_NAME}.${filter_name}.${MIN_RES}.cool",
            "${library_group}.${ASSEMBLY_NAME}.${filter_name}.${MIN_RES}.mcool" into LIBGROUP_FILTER_RES_COOLERS

    script:

    def res_str = params['bin'].resolutions.join(',')
    def balance_options = params['bin'].get('balance_options','')
    balance_options = ( balance_options ? "--balance-args \"${balance_options}\"": "")
    // balancing flag if it's requested
    def balance_flag = ( params['bin'].get('balance','false').toBoolean() ? "--balance ${balance_options}" : "--no-balance" )

    def merge_command = ""
    if( isSingleFile(coolers))
        merge_command = """
            ln -s \$(readlink -f ${coolers}) ${library_group}.${ASSEMBLY_NAME}.${filter_name}.${MIN_RES}.cool
        """
    else
        merge_command = """
            cooler merge ${library_group}.${ASSEMBLY_NAME}.${filter_name}.${MIN_RES}.cool ${coolers}
        """

    zoom_command = """
    cooler zoomify \
        --nproc ${task.cpus} \
        --out ${library_group}.${ASSEMBLY_NAME}.${filter_name}.${MIN_RES}.mcool \
        --resolutions ${res_str} \
        ${balance_flag} \
        ${library_group}.${ASSEMBLY_NAME}.${filter_name}.${MIN_RES}.cool
    """

    """
    ${merge_command}
    ${zoom_command}
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
    publishDir path: getOutputDir('stats_library_group'), mode: "copy"

    input:
    set val(library_group), file(stats) from LIBGROUP_STATS_TO_MERGE

    output:
    set library_group, "${library_group}.${ASSEMBLY_NAME}.stats" into LIBGROUP_STATS

    script:
    if( isSingleFile(stats))
        """
        ln -s ${stats} ${library_group}.${ASSEMBLY_NAME}.stats
        """
    else
        """
        pairtools stats --merge ${stats} -o ${library_group}.${ASSEMBLY_NAME}.stats
        """
}
