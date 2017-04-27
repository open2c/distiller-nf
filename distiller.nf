#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

boolean isSingleFile(object) {    
    (! [Collection, Object[], nextflow.util.BlankSeparatedList].any { 
        it.isAssignableFrom(object.getClass()) 
    } ) || (object.size() == 1)
}

LIB_RUN_SOURCES = Channel.from(
    params.fastq_paths.collect{
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
    tag { query }
    storeDir "intermediates/downloaded_fastqs"
 
    input:
    set val(library), val(run), val(query) from LIB_RUN_SRAS
     
    output:
    set library, run, 
        "${library}.${run}.1.fastq.gz", 
        "${library}.${run}.2.fastq.gz" into LIB_RUN_FASTQ_SRA
 
    script:
    sra_cli = query.tokenize(':')[-1]
    srr = sra_cli.tokenize('\\?')[0]
    sra_cli = (
        [srr]
        + sra_cli.tokenize('\\?')[-1].tokenize('&').collect{
            it.startsWith('start=') 
            ? (' --minSpotId '+it.tokenize('=')[1])
            : it.startsWith('end=') 
                ? (' --maxSpotId '+it.tokenize('=')[1])
                : ''
        }).join(' ')

    if( query.startsWith('sra:') )
        """
        fastq-dump -F ${sra_cli} --split-files --gzip
        mv ${srr}_1.fastq.gz ${library}.${run}.1.fastq.gz
        mv ${srr}_2.fastq.gz ${library}.${run}.2.fastq.gz
        """
    else
        error "Runs can be defined with one line only with SRA"

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

    tag { "library:${library} run:${run} side:${side}" }
    publishDir path:"fastqc/", mode:"copy"

    input:
    set val(library), val(run), val(side), file(fastq) from LIB_RUN_SIDE_FASTQS_FOR_QC

    output:
    set library, run, side,  
        "${library}.${run}.${side}_fastqc.html", 
        "${library}.${run}.${side}_fastqc.zip" into LIB_RUN_SIDE_FASTQCS

    """
    mkdir -p ./temp_fastqc/
    ln -s \$(readlink -f ${fastq}) ./temp_fastqc/${library}.${run}.${side}.fastq.gz
    fastqc -o ./ -f fastq ./temp_fastqc/${library}.${run}.${side}.fastq.gz
    """
          
}



/* 
 * Chunk fastqs
 */ 

LIB_RUN_FASTQS_NO_CHUNK = Channel.create()
LIB_RUN_FASTQS_FOR_CHUNK = Channel.create()
LIB_RUN_FASTQS
    .choice(LIB_RUN_FASTQS_NO_CHUNK, LIB_RUN_FASTQS_FOR_CHUNK) {
    it -> params.get('chunksize', 0) == 0 ? 0 : 1
}


process chunk_fastqs {
    tag { "library:${library} run:${run}" }
    storeDir "intermediates/fastq_chunks"

    input:
    set val(library), val(run),file(fastq1), file(fastq2) from LIB_RUN_FASTQS_FOR_CHUNK

    output:
    set library, run, 
        "${library}.${run}.*.1.fastq.gz", 
        "${library}.${run}.*.2.fastq.gz" into LIB_RUN_FASTQ_CHUNKED


    script:
    chunksize_lines = 4 * params.chunksize
   
    """
    zcat ${fastq1} | split -l ${chunksize_lines} -d \
        --filter 'gzip > \$FILE.1.fastq.gz' - \
        ${library}.${run}.

    zcat ${fastq2} | split -l ${chunksize_lines} -d \
        --filter 'gzip > \$FILE.2.fastq.gz' - \
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
             params.genome.bwa_index_basepath.split('/')[-1],
             file(params.genome.bwa_index_basepath+'.amb'),
             file(params.genome.bwa_index_basepath+'.ann'),
             file(params.genome.bwa_index_basepath+'.bwt'),
             file(params.genome.bwa_index_basepath+'.pac'),
             file(params.genome.bwa_index_basepath+'.sa')
            ]])

/*
 * Map fastq files
 */
process map_runs {
    tag { "library:${library} run:${run} chunk:${chunk}" }
    storeDir "intermediates/sam/runs"
 
    input:
    set val(library), val(run), val(chunk), file(fastq1), file(fastq2) from LIB_RUN_CHUNK_FASTQ
    set val(bwa_index_base), file(bwa_index_amb), file(bwa_index_ann), 
        file(bwa_index_bwt), file(bwa_index_pac), file(bwa_index_sa) from BWA_INDEX.first()
     
    output:
    set library, run, "${library}.${run}.${chunk}.bam" into LIB_RUN_CHUNK_BAMS
 
    """
    bwa mem -SP ${bwa_index_base} ${fastq1} ${fastq2} \
        | samtools view -bS > ${library}.${run}.${chunk}.bam
    """
}


/*
 * Parse mapped bams
 */

LIB_RUN_CHUNK_BAMS
     .groupTuple(by: [0, 1])
     .set {LIB_RUN_BAMS}


process parse_runs {
    tag { "library:${library} run:${run}" }
    storeDir "intermediates/pairsam/runs"
    publishDir path:"stats/runs/", pattern: "*.stats", mode:"copy"
 
    input:
    set val(library), val(run), file(bam) from LIB_RUN_BAMS
     
    output:
    set library, run, "${library}.${run}.pairsam.gz" into LIB_RUN_PAIRSAMS
    set library, run, "${library}.${run}.stats" into LIB_RUN_STATS
 
    script:
    dropsam_flag = params.get('drop_sam','false').toBoolean() ? '--drop-sam' : ''
    dropreadid_flag = params.get('drop_readid','false').toBoolean() ? '--drop-readid' : ''
    stats_command = (params.get('do_stats', 'true').toBoolean() ?
        "python -m pairsamtools stats ${library}.${run}.pairsam.gz -o ${library}.${run}.stats" :
        "touch ${library}.${run}.stats" )

    if( isSingleFile(bam))
        """
        python -m pairsamtools parse ${dropsam_flag} ${dropreadid_flag} ${bam} \
            | python -m pairsamtools sort -o ${library}.${run}.pairsam.gz
        ${stats_command}
        """
    else 
        """
        cat <( samtools merge - ${bam} | samtools view -H ) \
            <( samtools cat ${bam} | samtools view ) \
            | python -m pairsamtools parse ${dropsam_flag} ${dropreadid_flag} \
            | python -m pairsamtools sort -o ${library}.${run}.pairsam.gz

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
    tag { "library:${library}" }
    storeDir "intermediates/pairsam/libraries"
 
    input:
    set val(library), file(run_pairsam) from LIB_PAIRSAMS_TO_MERGE
     
    output:
    set library, "${library}.pairsam.gz" into LIB_PAIRSAMS

    script:
    if( isSingleFile(run_pairsam))
        """
        ln -s ${run_pairsam} ${library}.pairsam.gz
        """
    else
        """
        python -m pairsamtools merge ${run_pairsam} -o ${library}.pairsam.gz
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
    tag { "library:${library}" }
    publishDir path:"stats/libraries", pattern: "*.stats", mode:"copy"
 
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
        python -m pairsamtools stats --merge ${run_stats} -o ${library}.stats
        """
}



/*
 * Make pairs bam
 */

process make_pairs_bam {
    tag { "library:${library}" }
    publishDir path:'./', saveAs: {
      if( it.endsWith('.nodups.pairs.gz' ))
        return "pairs/libraries/${library}.nodups.pairs.gz"

      if( it.endsWith('.nodups.bam' ))
        return "sam/libraries/${library}.nodups.bam"

      if( it.endsWith('.dups.pairs.gz' ))
        return "pairs/libraries/${library}.dups.pairs.gz"

      if( it.endsWith('.nodups.bam' ))
        return "sam/libraries/${library}.dups.bam"

      if( it.endsWith('.unmapped.pairs.gz' ))
        return "pairs/libraries/${library}.unmapped.pairs.gz"

      if( it.endsWith('.nodups.bam' ))
        return "sam/libraries/${library}.unmapped.bam"

      if( it.endsWith('.dedup.stats' ))
        return "stats/libraries/${library}.dedup.stats.tsv"
    }
 
    input:
    set val(library), file(pairsam_lib) from LIB_PAIRSAMS
     
    output:
    set library, "${library}.nodups.pairs.gz", "${library}.nodups.bam",
                 "${library}.dups.pairs.gz", "${library}.dups.bam", 
                 "${library}.unmapped.pairs.gz", 
                 "${library}.unmapped.bam" into LIB_PAIRS_BAMS
    set library, "${library}.dedup.stats" into LIB_DEDUP_STATS
 
     """
        python -m pairsamtools select '(PAIR_TYPE == "CX") or (PAIR_TYPE == "LL")' \
            ${pairsam_lib} \
            --output-rest >( python -m pairsamtools split \
                --output-pairs ${library}.unmapped.pairs.gz \
                --output-sam ${library}.unmapped.bam \
                ) | \
        python -m pairsamtools dedup \
            --output \
                >( python -m pairsamtools split \
                    --output-pairs ${library}.nodups.pairs.gz \
                    --output-sam ${library}.nodups.bam \
                 ) \
            --output-dups \
                >( python -m pairsamtools markasdup \
                    | python -m pairsamtools split \
                        --output-pairs ${library}.dups.pairs.gz \
                        --output-sam ${library}.dups.bam \
                 ) \
            --stats-file ${library}.dedup.stats

     """
}


/*
 * Index .pairs.gz files with pairix
 */

LIB_PAIRS_BAMS
    .map {v -> tuple(v[0], v[1])}
    .set {LIB_PAIRS}

process index_pairs{
    tag { "library:${library}" }
    publishDir path:"pairs/libraries/", saveAs: {"${library}.nodups.pairs.gz.px2"}

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

CHROM_SIZES = Channel.from([ file(params.genome.chrom_sizes_path) ])

process make_library_coolers{
    tag { "library:${library} resolution:${res}" }
    publishDir path:"coolers/libraries/", saveAs: {"${library}.${res}.cool"}

    input:
        set val(library), file(pairs_lib), file(pairs_index_lib) from LIB_IDX_PAIRS
        each res from params.cooler_resolutions
        file(chrom_sizes) from CHROM_SIZES.first()

    output:
        set library, res, "${library}.${res}.cool" into LIB_RES_COOLERS

    shell:
        """
        cooler cload pairix \
            --assembly ${params.genome.assembly} \
            ${chrom_sizes}:${res} ${pairs_lib} ${library}.${res}.cool
        """
}


/*
 * Merge .cool matrices for library groups.
 */ 


LIBRARY_GROUPS = Channel.from( params.library_groups.collect{ k, v -> [k, v] })

LIBRARY_GROUPS
    .combine(LIB_RES_COOLERS)
    .filter{ it[1].contains(it[2]) } 
    .map {library_group, libraries, library, res, file -> tuple(library_group, res, file)}
    .groupTuple(by: [0, 1])
    .set { LIBGROUP_RES_COOLERS_TO_MERGE }

process make_library_group_coolers{
    tag {"library_group:${library_group} resolution:${res}"}
    publishDir path:"coolers/library_groups/", saveAs: {"${library_group}.${res}.cool"}

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


LIBRARY_GROUPS = Channel.from( params.library_groups.collect{ k, v -> [k, v] })

LIBRARY_GROUPS
    .combine(LIB_STATS.mix(LIB_DEDUP_STATS))
    .filter{ it[1].contains(it[2]) } 
    .map {library_group, libraries, library, stats -> tuple(library_group, stats)}
    .groupTuple()
    .set { LIBGROUP_STATS_TO_MERGE }


process merge_stats_libraries_into_groups {
    tag { "library_group:${library_group}" }
    publishDir path:"stats/library_groups", pattern: "*.stats", mode:"copy"
 
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
        python -m pairsamtools stats --merge ${stats} -o ${library_group}.stats
        """
}


