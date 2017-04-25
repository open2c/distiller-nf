#!/usr/bin/env nextflow

LIBRARY_RUN_FASTQS = Channel.from(
    params.fastq_paths.collect{
        k, v -> v.collect{k2, v2 -> [k,k2]+v2}}.sum())
    .map{ v -> [v[0], v[1], file(v[2]), file(v[3])] }

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
    publishDir path:"sam/runs/", saveAs: {"${library}.${run}.bam"}
 
    input:
    set val(library), val(run), file(fastq1), file(fastq2) from LIBRARY_RUN_FASTQS
    set val(bwa_index_base), file(bwa_index_amb), file(bwa_index_ann), 
        file(bwa_index_bwt), file(bwa_index_pac), file(bwa_index_sa) from BWA_INDEX.first()
     
    output:
    set val(library), val(run), file('mapped.bam') into mapped_bams
 
    """
    bwa mem -SP ${bwa_index_base} ${fastq1} ${fastq2} | samtools view -bS > mapped.bam
    """
}


/*
 * Parse mapped bams
 */
process parse_runs {
    publishDir path:"pairsam/runs/", saveAs: {"${library}.${run}.pairsam.gz"}
 
    input:
    set val(library), val(run), file(mapped_bam) from mapped_bams
     
    output:
    set val(library), val(run), file('pairsam.gz') into pairsam_runs
 
    """
    python -m pairsamtools parse ${mapped_bam} | python -m pairsamtools sort -o pairsam.gz
    """
}



 
/*
 * Merge runs
 */

pairsam_runs
    .map {library, run, file -> tuple(library, file)}
    .groupTuple()
    .set {pairsam_runs_by_library}

process merge_runs_into_libraries {
    publishDir path:"pairsam/libraries/", saveAs: {"${library}.pairsam.gz"}
 
    input:
    set val(library), file(pairsam_per_run) from pairsam_runs_by_library
     
    output:
    set val(library), file('lib_pairsam.gz') into pairsam_libs
 
    """
    python -m pairsamtools merge ${pairsam_per_run} -o lib_pairsam.gz
    """
}


/*
 * Make pairs bam
 */


process make_pairs_bam {
    publishDir path:'./', saveAs: {
      if( it == 'nodups.pairs.gz' )
        return "pairs/libraries/${library}.nodups.pairs.gz"

      if( it == 'nodups.bam' )
        return "sam/libraries/${library}.nodups.bam"

      if( it == 'dups.pairs.gz' )
        return "pairs/libraries/${library}.dups.pairs.gz"

      if( it == 'nodups.bam' )
        return "sam/libraries/${library}.dups.bam"

      if( it == 'unmapped.pairs.gz' )
        return "pairs/libraries/${library}.unmapped.pairs.gz"

      if( it == 'nodups.bam' )
        return "sam/libraries/${library}.unmapped.bam"

      if( it == 'dedup.stats' )
        return "stats/libraries/${library}.dedup.stats.tsv"
    }
 
    input:
    set val(library), file(pairsam_lib) from pairsam_libs
     
    output:
    set val(library), file('nodups.pairs.gz'), file('nodups.bam'),
                      file('dups.pairs.gz'), file('dups.bam'), 
                      file('unmapped.pairs.gz'), file('unmapped.bam'),
                      file('dedup.stats') into pairs_libs
 
     """
        python -m pairsamtools select '(PAIR_TYPE == "CX") or (PAIR_TYPE == "LL")' \
            ${pairsam_lib} \
            --output-rest >( python -m pairsamtools split \
                --output-pairs unmapped.pairs.gz \
                --output-sam unmapped.bam \
                ) | \
        python -m pairsamtools dedup \
            --output \
                >( python -m pairsamtools split \
                    --output-pairs nodups.pairs.gz \
                    --output-sam nodups.bam \
                 ) \
            --output-dups \
                >( python -m pairsamtools markasdup \
                    | python -m pairsamtools split \
                        --output-pairs dups.pairs.gz \
                        --output-sam dups.bam \
                 ) \
            --stats-file dedup.stats

     """
}
