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
    ln -s \"\$(readlink -f ${fastq})\" ./temp_fastqc/${library}.${run}.${side}.fastq.gz
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
            .map { [it[0], it[1], 0, it[2], it[3]] } )
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
    set library, run, chunk, "${library}.${run}.${chunk}.bam" into LIB_RUN_CHUNK_BAMS
 
    """
    bwa mem -t ${task.cpus} -SP ${bwa_index_base} ${fastq1} ${fastq2} \
        | samtools view -bS > ${library}.${run}.${chunk}.bam \
        | cat
    """
}


/*
 * Parse mapped bams
 */

CHROM_SIZES = Channel.from([ file(params.input.genome.chrom_sizes_path) ])

process parse_runs {
    tag "library:${library} run:${run} chunk:${chunk} parsing"
    storeDir getIntermediateDir('pairsam_chunk')

    cpus params.parse_cpus
 
    input:
    set val(library), val(run), val(chunk), file(bam) from LIB_RUN_CHUNK_BAMS
    file(chrom_sizes) from CHROM_SIZES.first()
     
    output:
    set library, run, chunk, "${library}.${run}.${chunk}.pairsam.gz" into LIB_RUN_PAIRSAM_CHUNKS, LIB_RUN_PAIRSAM_CHUNKS_TO_COUNT
 
    script:
    dropsam_flag = params['map'].get('drop_sam','false').toBoolean() ? '--drop-sam' : ''
    dropreadid_flag = params['map'].get('drop_readid','false').toBoolean() ? '--drop-readid' : ''
    dropseq_flag = params['map'].get('drop_seq','false').toBoolean() ? '--drop-seq' : ''

        """
        mkdir ./tmp4sort
        pairsamtools parse ${dropsam_flag} ${dropreadid_flag} ${dropseq_flag} \
            -c ${chrom_sizes}  ${bam} \
                | pairsamtools sort --nproc ${task.cpus} \
                                    -o ${library}.${run}.${chunk}.pairsam.gz \
                                    --tmpdir ./tmp4sort \
                | cat

        rm -rf ./tmp4sort

        """
}


/*
 * Channels for pre-merging or merging into runs 
 */
LIB_RUN_PAIRSAMS_PREMERGE = Channel.create()
LIB_RUN_PAIRSAMS_MERGE_RUN = Channel.create()


/*
 * Adding total number of chunks to the chunks-channel
 * and reidrecting it either to pre-merging step or straight to merge-run.
 * Everything is explicitly synchronized !
 */
LIB_RUN_PAIRSAM_CHUNKS_TO_COUNT
                    .count()
                    .combine(LIB_RUN_PAIRSAM_CHUNKS)
                    .choice(LIB_RUN_PAIRSAMS_PREMERGE
                        ,LIB_RUN_PAIRSAMS_MERGE_RUN){
                             (it[0] > params.merge.merge_threshold) ? 0 : 1
                           }

/*
 * Get rid of that total number of chunks in pre-merge and merge-run Channels 
 */
LIB_RUN_PAIRSAMS_PREMERGE.map { it[1..-1] }.set{LIB_RUN_PAIRSAMS_PREMERGE}
LIB_RUN_PAIRSAMS_MERGE_RUN.map { it[1..-1] }.set{LIB_RUN_PAIRSAMS_MERGE_RUN}

/*
 * Group pairsam chunks in batches for premerge.
 * k3.sum() - combines chunk id-s of each pairsam in the batch
 * batch content is stochastic, thus resuming from after pre-merge step could be a problem!
 */
LIB_RUN_PAIRSAMS_PREMERGE
     .groupTuple(by: [0, 1], size: params.merge.merge_batch_size, remainder: true)
     .map {k1,k2,k3,v->[k1, k2, k3.sum(), v]}
     .set{LIB_RUN_PAIRSAMS_PREMERGE}

/*
 * Pre-merge .pairsams chunks in batches of merge_batch_size
 * This step is usefull for HPC and big number of chunks
 */
process premerge_pairsam_chunks {
    tag "library:${library} run:${run} batch:${batch} premerging"
    storeDir getIntermediateDir('pairsam_chunk_premerge')

    cpus params.merge_cpus
 
    input:
    set val(library), val(run), val(batch), file(run_chunk_pairsam) from LIB_RUN_PAIRSAMS_PREMERGE

    output:
    set library, run, batch, "${library}.${run}.${batch}.premerged.pairsam.gz" into LIB_RUN_PAIRSAMS_POST_PREMERGE

    script:
    if( isSingleFile(run_chunk_pairsam))
        """
        ln -s \"\$(readlink -f ${run_chunk_pairsam})\" ${library}.${run}.${batch}.premerged.pairsam.gz
        """
    else
        """
        pairsamtools merge ${run_chunk_pairsam} --nproc ${task.cpus} -o ${library}.${run}.${batch}.premerged.pairsam.gz
        """
}


/*
 * Mix (merge) channels with .pairsam chunks
 * one of these channels should be empty: chunks were either pre-merged or remain unchanged
 * now they are to be piped in the channel for merging into runs.
 */
LIB_RUN_PAIRSAMS_POST_PREMERGE.mix(LIB_RUN_PAIRSAMS_MERGE_RUN).set{LIB_RUN_PAIRSAMS_MERGE_RUN}

/*
 * groupby library and run: techincally k3 is no longer needed (consider omitting)
 */
LIB_RUN_PAIRSAMS_MERGE_RUN
     .groupTuple(by: [0, 1])
     .map {k1,k2,k3,v->[k1,k2,k3.sum(),v]}
     .set{LIB_RUN_PAIRSAMS_MERGE_RUN}



/*
 * Merge premerged or not-premerged chunks of .pairsam into runs
 */
process merge_pairsam_into_runs {
    tag "library:${library} run:${run}"
    storeDir getIntermediateDir('pairsam_run')
    publishDir path: getOutDir('stats_run'), pattern: "*.stats", mode:"copy"

    cpus params.merge_cpus
 
    input:
    set val(library), val(run), val(batch), file(run_batch_pairsam) from LIB_RUN_PAIRSAMS_MERGE_RUN
     
    output:
    set library, run, "${library}.${run}.pairsam.gz" into LIB_RUN_PAIRSAMS
    set library, run, "${library}.${run}.stats" into LIB_RUN_STATS
 
    script:
    stats_command = (params.get('do_stats', 'true').toBoolean() ?
        "pairsamtools stats ${library}.${run}.pairsam.gz -o ${library}.${run}.stats" :
        "touch ${library}.${run}.stats" )
    
    if( isSingleFile(run_batch_pairsam))
        """
        ln -s \"\$(readlink -f ${run_batch_pairsam})\" ${library}.${run}.pairsam.gz

        ${stats_command}
        """
    else
        """
        pairsamtools merge ${run_batch_pairsam} --nproc ${task.cpus} -o ${library}.${run}.pairsam.gz

        ${stats_command}
        """

}







/*
 * Grouby .pairsam by library and merge
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
        ln -s \"\$(readlink -f ${run_pairsam})\" ${library}.pairsam.gz
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
    
    script:
    dropsam = params['map'].get('drop_sam','false').toBoolean()
    if(dropsam) 
        """
        pairsamtools select '(pair_type == "CX") or (pair_type == "LL")' \
            ${pairsam_lib} \
            --output-rest >( pairsamtools split \
                --output-pairs ${library}.unmapped.pairs.gz \
                ) | \
            pairsamtools dedup \
                --max-mismatch ${params.filter.pcr_dups_max_mismatch_bp} \
                --output \
                    >( pairsamtools split \
                        --output-pairs ${library}.nodups.pairs.gz \
                     ) \
                --output-dups \
                    >( pairsamtools markasdup \
                        | pairsamtools split \
                            --output-pairs ${library}.dups.pairs.gz \
                     ) \
                --stats-file ${library}.dedup.stats \
                | cat

        touch ${library}.unmapped.bam
        touch ${library}.nodups.bam
        touch ${library}.dups.bam

        """
    else 
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

    script:
    if( isSingleFile(coolers))
        """
        ln -s \$(readlink -f ${coolers}) ${library_group}.${res}.cool
        """
    else
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


