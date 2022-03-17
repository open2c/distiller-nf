// Import generic module functions
include { initOptions; getSoftwareName; getOutputDir } from './functions'
include { isSingleFile } from './functions'

params.options = [:]
options        = initOptions(params.options)
directory = getOutputDir('pairs_library')

ASSEMBLY_NAME = params['input'].genome.assembly_name // TODO: move to the parameters dictionary, and below:


process MERGE_DEDUP_SPLITBAM {
    tag "library:${library} run:${run}"
    label 'process_medium'
    publishDir "${directory}", mode: params.publish_dir_mode

    conda (params.enable_conda ? "bioconda::pairtools" : null)
//        if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
//            container "https://depot.galaxyproject.org/singularity/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2880dd9d8ad0a7b221d4eacda9a818e92983128d-0"
//        } else {
//            container "quay.io/biocontainers/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2880dd9d8ad0a7b221d4eacda9a818e92983128d-0"
//        }

    input:
    tuple val(library), file(run_pairsam)

    output:
    tuple val(library), path("${library}.${ASSEMBLY_NAME}.nodups.pairs.gz"),
                 path("${library}.${ASSEMBLY_NAME}.nodups.pairs.gz.px2"),
                 path("${library}.${ASSEMBLY_NAME}.nodups.bam"),
                 path("${library}.${ASSEMBLY_NAME}.dups.pairs.gz"),
                 path("${library}.${ASSEMBLY_NAME}.dups.bam"),
                 path("${library}.${ASSEMBLY_NAME}.unmapped.pairs.gz"),
                 path("${library}.${ASSEMBLY_NAME}.unmapped.bam"), emit: output

    tuple val(library), path("${library}.${ASSEMBLY_NAME}.dedup.stats"), emit: stats

    path  "*.version.txt"         , emit: version


    script:
    def software = getSoftwareName(task.process)

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


        pairtools --version > ${software}.version.txt

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

        pairtools --version > ${software}.version.txt

        """

}

