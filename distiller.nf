#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

/*
 * Miscellaneous code for the pipeline
 */

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

LIB_RUN_SOURCES = Channel.from(
    params.input.raw_reads_paths.collect{
        k, v -> v.collect{k2, v2 -> [k,k2]+v2}}.sum())

/*
 * Download fastqs from SRA.
 */

LIB_RUN_FASTQS = Channel.create()
LIB_RUN_SRAS = Channel.create()
LIB_RUN_SOURCES.choice(LIB_RUN_SRAS, LIB_RUN_FASTQS) {
    a -> a.size() == 3 ? 0 : 1
    }



process download_sra {
    tag "$query"
    storeDir getIntermediateDir('downloaded_fastqs')
 
    input:
    set val(library), val(run), val(query) from LIB_RUN_SRAS
     
    output:
    set library, run, 
        "${library}.${run}.1.fastq.gz", 
        "${library}.${run}.2.fastq.gz" into LIB_RUN_FASTQ_SRA
 
    script:
    if( query == null) error "No files provided for library ${library}, run ${run}"

    sra_cli = query.tokenize(':')[-1]
    srr = sra_cli.contains('?') ? sra_cli.tokenize('\\?')[0] : sra_cli
    srrnum = srr.substring(3)
    srrnum_three_digs = srrnum.take(3)
    sra_cli = srr + (sra_cli.contains('?') ? (
        sra_cli.tokenize('\\?')[-1].tokenize('&').collect{
            it.startsWith('start=') 
            ? (' --minSpotId '+it.tokenize('=')[1])
            : (it.startsWith('end=') ? (' --maxSpotId '+it.tokenize('=')[1]) : '')
        }).join(' ') : '')

    if( query.startsWith('sra:') ) {
        if( query.contains('?') ) {
            """
            fastq-dump -F ${sra_cli} --split-files --gzip
            mv ${srr}_1.fastq.gz ${library}.${run}.1.fastq.gz
            mv ${srr}_2.fastq.gz ${library}.${run}.2.fastq.gz
            """
        }
        else {
            """ 
            wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR${srrnum_three_digs}/${srr}/${srr}.sra -O ${srr}.sra
            fastq-dump -F ${srr}.sra --split-files --gzip
            mv ${srr}_1.fastq.gz ${library}.${run}.1.fastq.gz
            mv ${srr}_2.fastq.gz ${library}.${run}.2.fastq.gz
            """
        }
    }
    else {
        error "Runs can be defined with one line only with SRA"
    }

}

LIB_RUN_FASTQS
    .mix(LIB_RUN_FASTQ_SRA)
    .map{ v -> [v[0], v[1], file(v[2]), file(v[3])] }
    .set { LIB_RUN_FASTQS }


/*
 * FastQC the input files
 */

LIB_RUN_FASTQS_FOR_QC = Channel.create()
LIB_RUN_FASTQS
    .tap(LIB_RUN_FASTQS_FOR_QC)
    .set{LIB_RUN_FASTQS}

LIB_RUN_FASTQS_FOR_QC
    .filter { it -> params.get('do_fastqc', 'false').toBoolean() }
    .map{ v -> [v[0], v[1], [[1,file(v[2])], [2,file(v[3])]]]} 
    .flatMap{ 
        vs -> vs[2].collect{ 
            it -> [vs[0],
                   vs[1], 
                   it[0],
                   it[1]] } }
    .set {LIB_RUN_SIDE_FASTQS_FOR_QC}

process fastqc{

    tag "library:${library} run:${run} side:${side}"
    publishDir path: getOutDir('fastqc'), mode:"copy"
    cpus params.fastqc_cpus

    input:
    set val(library), val(run), val(side), file(fastq) from LIB_RUN_SIDE_FASTQS_FOR_QC

    output:
    set library, run, side,  
        "${library}.${run}.${side}_fastqc.html", 
        "${library}.${run}.${side}_fastqc.zip" into LIB_RUN_SIDE_FASTQCS

    """
    mkdir -p ./temp_fastqc/
    ln -s \$(readlink -f ${fastq}) ./temp_fastqc/${library}.${run}.${side}.fastq.gz
    fastqc --threads ${task.cpus} -o ./ -f fastq ./temp_fastqc/${library}.${run}.${side}.fastq.gz
    rm -r ./temp_fastqc/
    """
          
}



/* 
 * Chunk fastqs
 */ 

LIB_RUN_FASTQS_NO_CHUNK = Channel.create()
LIB_RUN_FASTQS_FOR_CHUNK = Channel.create()
LIB_RUN_FASTQS
    .choice(LIB_RUN_FASTQS_NO_CHUNK, LIB_RUN_FASTQS_FOR_CHUNK) {
    it -> params['map'].get('chunksize', 0) == 0 ? 0 : 1
}


process chunk_fastqs {
    tag "library:${library} run:${run}"
    storeDir getIntermediateDir('fastq_chunks')

    cpus params.chunk_cpus

    input:
    set val(library), val(run),file(fastq1), file(fastq2) from LIB_RUN_FASTQS_FOR_CHUNK

    output:
    set library, run, 
        "${library}.${run}.*.1.fastq.gz", 
        "${library}.${run}.*.2.fastq.gz" into LIB_RUN_FASTQ_CHUNKED


    script:
    chunksize_lines = 4 * params['map'].chunksize
   
    """
    zcat ${fastq1} | split -l ${chunksize_lines} -d \
        --filter 'pbgzip -c -n ${task.cpus} > \$FILE.1.fastq.gz' - \
        ${library}.${run}.

    zcat ${fastq2} | split -l ${chunksize_lines} -d \
        --filter 'pbgzip -c -n ${task.cpus} > \$FILE.2.fastq.gz' - \
        ${library}.${run}.
    """
}


LIB_RUN_FASTQ_CHUNKED
    .map{v->[v[0],
             v[1],
             (v[2] instanceof Collection ? v[2] : [v[2]]),
             (v[3] instanceof Collection ? v[3] : [v[3]])
             ] }
    .into {CHUNKS_SIDE_1; CHUNKS_SIDE_2}

CHUNKS_SIDE_1
    .flatMap{ 
        vs -> vs[2].collect{ 
            it -> [vs[0],
                   vs[1], 
                   it.toString().tokenize('.')[-4], 
                   it] }}
    .set{CHUNKS_SIDE_1}

CHUNKS_SIDE_2
    .flatMap{ 
        vs -> vs[3].collect{ 
            it -> [vs[0],
                   vs[1], 
                   it.toString().tokenize('.')[-4],
                   it] }}
    .set{CHUNKS_SIDE_2}

CHUNKS_SIDE_1
    .phase(CHUNKS_SIDE_2) { [it[0], it[1], it[2]] }
    .map{ [it[0][0], it[0][1], it[0][2], it[0][3], it[1][3]] }
    .set{ LIB_RUN_CHUNK_FASTQ }

LIB_RUN_CHUNK_FASTQ
    .mix(
        LIB_RUN_FASTQS_NO_CHUNK
            .map { [it[0], it[1], 0, it[2]] } )
    .set{LIB_RUN_CHUNK_FASTQ}

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
process map_runs {
    tag "library:${library} run:${run} chunk:${chunk}"
    storeDir getIntermediateDir('bam_run')

    cpus params.map_cpus
 
    input:
    set val(library), val(run), val(chunk), file(fastq1), file(fastq2) from LIB_RUN_CHUNK_FASTQ
    set val(bwa_index_base), file(bwa_index_files) from BWA_INDEX.first()
     
    output:
    set library, run, "${library}.${run}.${chunk}.bam" into LIB_RUN_CHUNK_BAMS
 
    """
    bwa mem -t ${task.cpus} -SP ${bwa_index_base} ${fastq1} ${fastq2} \
        | samtools view -bS > ${library}.${run}.${chunk}.bam \
        | cat
    """
}


/*
 * Parse mapped bams
 */

LIB_RUN_CHUNK_BAMS
     .groupTuple(by: [0, 1])
     .set {LIB_RUN_BAMS}

CHROM_SIZES = Channel.from([ file(params.input.genome.chrom_sizes_path) ])

process parse_runs {
    tag "library:${library} run:${run}"
    storeDir getIntermediateDir('pairsam_run')
    publishDir path: getOutDir('stats_run'), pattern: "*.stats", mode:"copy"

    cpus params.parse_cpus
 
    input:
    set val(library), val(run), file(bam) from LIB_RUN_BAMS
    file(chrom_sizes) from CHROM_SIZES.first()
     
    output:
    set library, run, "${library}.${run}.pairsam.gz" into LIB_RUN_PAIRSAMS
    set library, run, "${library}.${run}.stats" into LIB_RUN_STATS
 
    script:
    dropsam_flag = params['map'].get('drop_sam','false').toBoolean() ? '--drop-sam' : ''
    dropreadid_flag = params['map'].get('drop_readid','false').toBoolean() ? '--drop-readid' : ''
    dropseq_flag = params['map'].get('drop_seq','false').toBoolean() ? '--drop-seq' : ''
    stats_command = (params.get('do_stats', 'true').toBoolean() ?
        "pairsamtools stats ${library}.${run}.pairsam.gz -o ${library}.${run}.stats" :
        "touch ${library}.${run}.stats" )
    n_parse_processes = (int)Math.ceil(task.cpus / 2)
    n_parse_processes = n_parse_processes < 1 ? 1 : n_parse_processes

    if( isSingleFile(bam))
        """
        mkdir ./tmp4sort
        pairsamtools parse ${dropsam_flag} ${dropreadid_flag} ${dropseq_flag} \
            -c ${chrom_sizes}  ${bam} \
                | pairsamtools sort --nproc ${task.cpus} \
                                    -o ${library}.${run}.pairsam.gz \
                                    --tmpdir ./tmp4sort \
                | cat

        rm -rf ./tmp4sort

        ${stats_command}
        """
    else 
        """
        mkdir ./tmp4sort
        mkdir ./tmp_pairsam
        parallel -P${n_parse_processes} 'pairsamtools parse \
            ${dropsam_flag} ${dropseq_flag} ${dropreadid_flag} -c ${chrom_sizes} {} \
            | pairsamtools sort --nproc 4 \
                                -o ./tmp_pairsam/{}.pairsam.gz \
                                --tmpdir ./tmp4sort ' ::: ${bam}

        pairsamtools merge ./tmp_pairsam/* --nproc ${task.cpus} -o ${library}.${run}.pairsam.gz
    
        rm -rf ./tmp4sort
        rm -rf ./tmp_pairsam

        ${stats_command}
        """
}



/*
 * Merge .pairsams for runs into libraries
 */

LIB_RUN_PAIRSAMS
    .map {library, run, file -> tuple(library, file)}
    .groupTuple()
    .set {LIB_PAIRSAMS_TO_MERGE}

process merge_runs_into_libraries {
    tag "library:${library}"
    storeDir getIntermediateDir('pairsam_library')

    cpus params.merge_cpus
 
    input:
    set val(library), file(run_pairsam) from LIB_PAIRSAMS_TO_MERGE
     
    output:
    set library, "${library}.pairsam.gz" into LIB_PAIRSAMS

    script:
    if( isSingleFile(run_pairsam))
        """
        ln -s \$(readlink -f ${run_pairsam}) ${library}.pairsam.gz
        """
    else
        """
        pairsamtools merge ${run_pairsam} --nproc ${task.cpus} -o ${library}.pairsam.gz
        """
}

/*
 * Merge .stats for runs into libraries
 */

LIB_RUN_STATS
    .map {library, run, stats -> tuple(library, stats)}
    .groupTuple()
    .set {LIB_STATS_TO_MERGE}

process merge_stats_runs_into_libraries {
    tag "library:${library}"
    publishDir path: getOutDir('stats_library'), pattern: "*.stats", mode:"copy"
 
    input:
    set val(library), file(run_stats) from LIB_STATS_TO_MERGE
     
    output:
    set library, "${library}.stats" into LIB_STATS

    script:
    if( isSingleFile(run_stats))
        """
        ln -s ${run_stats} ${library}.stats
        """
    else
        """
        pairsamtools stats --merge ${run_stats} -o ${library}.stats
        """
}



/*
 * Make pairs bam
 */

process filter_make_pairs {
    tag "library:${library}"
    publishDir path:'.', mode:"copy", saveAs: {
      if( it.endsWith('.nodups.pairs.gz' ))
        return getOutDir("pairs_library") +"/${library}.nodups.pairs.gz"

      if( it.endsWith('.nodups.bam' ))
        return getOutDir("bams_library") +"/${library}.nodups.bam"

      if( it.endsWith('.dups.pairs.gz' ))
        return getOutDir("pairs_library") +"/${library}.dups.pairs.gz"

      if( it.endsWith('.nodups.bam' ))
        return getOutDir("bams_library") +"/${library}.dups.bam"

      if( it.endsWith('.unmapped.pairs.gz' ))
        return getOutDir("pairs_library") +"/${library}.unmapped.pairs.gz"

      if( it.endsWith('.nodups.bam' ))
        return getOutDir("bams_library") +"/${library}.unmapped.bam"

      if( it.endsWith('.dedup.stats' ))
        return getOutDir("stats_library") +"/${library}.dedup.stats.tsv"
    }

    cpus params.filter_make_pairs_cpus
 
    input:
    set val(library), file(pairsam_lib) from LIB_PAIRSAMS
     
    output:
    set library, "${library}.nodups.pairs.gz", "${library}.nodups.bam",
                 "${library}.dups.pairs.gz", "${library}.dups.bam", 
                 "${library}.unmapped.pairs.gz", 
                 "${library}.unmapped.bam" into LIB_PAIRS_BAMS
    set library, "${library}.dedup.stats" into LIB_DEDUP_STATS
 
    """
    pairsamtools select '(pair_type == "CX") or (pair_type == "LL")' \
        ${pairsam_lib} \
        --output-rest >( pairsamtools split \
            --output-pairs ${library}.unmapped.pairs.gz \
            --output-sam ${library}.unmapped.bam \
            ) | \
        pairsamtools dedup \
            --max-mismatch ${params.filter.pcr_dups_max_mismatch_bp} \
            --output \
                >( pairsamtools split \
                    --output-pairs ${library}.nodups.pairs.gz \
                    --output-sam ${library}.nodups.bam \
                 ) \
            --output-dups \
                >( pairsamtools markasdup \
                    | pairsamtools split \
                        --output-pairs ${library}.dups.pairs.gz \
                        --output-sam ${library}.dups.bam \
                 ) \
            --stats-file ${library}.dedup.stats \
            | cat


    """
}


/*
 * Index .pairs.gz files with pairix
 */

LIB_PAIRS_BAMS
    .map {v -> tuple(v[0], v[1])}
    .set {LIB_PAIRS}

process index_pairs{
    tag "library:${library}"
    publishDir path: getOutDir('pairs_library'), mode:"copy", saveAs: {"${library}.nodups.pairs.gz.px2"}

    input:
    set val(library), file(pairs_lib) from LIB_PAIRS
     
    output:
    set val(library), file(pairs_lib), file('*.px2') into LIB_IDX_PAIRS
 
    """
    pairix ${pairs_lib}
    """
}



/*
 * Bin indexed .pairs into .cool matrices.
 */ 

CHROM_SIZES = Channel.from([ file(params.input.genome.chrom_sizes_path) ])

process bin_library_pairs{
    tag "library:${library} resolution:${res}"
    publishDir path: getOutDir('coolers_library'), mode:"copy", saveAs: {"${library}.${res}.cool"}

    cpus params.bin_cpus

    input:
        set val(library), file(pairs_lib), file(pairs_index_lib) from LIB_IDX_PAIRS
        each res from params['bin'].resolutions
        file(chrom_sizes) from CHROM_SIZES.first()

    output:
        set library, res, "${library}.${res}.cool" into LIB_RES_COOLERS

    """
    cooler cload pairix \
        --nproc ${task.cpus} \
        --assembly ${params.input.genome.assembly} \
        ${chrom_sizes}:${res} ${pairs_lib} ${library}.${res}.cool
    """
}


/*
 * Merge .cool matrices for library groups.
 */ 


LIBRARY_GROUPS = Channel.from( params.input.library_groups.collect{ k, v -> [k, v] })

LIBRARY_GROUPS
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
        set library_group, res, "${library_group}.${res}.cool" into LIBGROUP_RES_COOLERS

    """
    cooler merge ${library_group}.${res}.cool ${coolers}
    """
}


/*
 * Merge .stats for library groups
 */ 


LIBRARY_GROUPS = Channel.from( params.input.library_groups.collect{ k, v -> [k, v] })

LIBRARY_GROUPS
    .combine(LIB_STATS.mix(LIB_DEDUP_STATS))
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
        pairsamtools stats --merge ${stats} -o ${library_group}.stats
        """
}


