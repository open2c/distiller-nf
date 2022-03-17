// Import generic module functions
include { getSoftwareName; initOptions; saveFiles; getOutputDir } from './functions'
include { isSingleFile } from './functions'

params.options = [:]
options        = initOptions(params.options)
directory = getOutputDir('stats_library_group')

ASSEMBLY_NAME = params['input'].genome.assembly_name // TODO: move to the parameters dictionary

process MERGE_STATS {
    tag "library_group:${library_group}"
    label 'process_low'
    publishDir "${directory}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::pairtools" : null)
//        if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
//            container "https://depot.galaxyproject.org/singularity/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2880dd9d8ad0a7b221d4eacda9a818e92983128d-0"
//        } else {
//            container "quay.io/biocontainers/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2880dd9d8ad0a7b221d4eacda9a818e92983128d-0"
//        }

    input:
    tuple val(library_group), file(stats)

    output:
    tuple val(library_group), path("${library_group}.${ASSEMBLY_NAME}.stats"), emit: output
    path  "*.version.txt"         , emit: version


    script:
    def software = getSoftwareName(task.process)

    if( isSingleFile(stats))
        """
        ln -s ${stats} ${library_group}.${ASSEMBLY_NAME}.stats

        echo "NA" > ${software}.version.txt
        """
    else
        """
        pairtools stats --merge ${stats} -o ${library_group}.${ASSEMBLY_NAME}.stats

        pairtools --version > ${software}.version.txt
        """
}


