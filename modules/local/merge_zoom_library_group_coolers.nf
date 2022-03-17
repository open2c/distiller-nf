// Import generic module functions
include { initOptions; getOutputDir; getSoftwareName; isSingleFile } from './functions'

params.options = [:]
options        = initOptions(params.options)
directory      = getOutputDir('coolers_library_group')

MIN_RES = params['bin'].resolutions.collect { it as int }.min() // TODO: move to parameters dictionary
ASSEMBLY_NAME = params['input'].genome.assembly_name // TODO: move to the parameters dictionary
pairsgz_decompress_command = 'bgzip -cd -@ 3'

process MERGE_ZOOM {
    tag "library:${library} filter:${filter_name}"
    label 'process_medium'
    publishDir "${directory}", mode: params.publish_dir_mode

    conda (params.enable_conda ? "bioconda::sra-tools>=2.8.1 bioconda::pbgzip" : null)
//        if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
//            container "https://depot.galaxyproject.org/singularity/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2880dd9d8ad0a7b221d4eacda9a818e92983128d-0"
//        } else {
//            container "quay.io/biocontainers/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2880dd9d8ad0a7b221d4eacda9a818e92983128d-0"
//        }

    input:
    val(filter_name)
    tuple val(library_group), file(coolers)

    output:
    tuple val(library_group), val(filter_name),
        path("${library_group}.${ASSEMBLY_NAME}.${filter_name}.${MIN_RES}.cool"),
        path("${library_group}.${ASSEMBLY_NAME}.${filter_name}.${MIN_RES}.mcool"), emit: output
    path  "*.version.txt"         , emit: version


    script:
    def software = getSoftwareName(task.process)

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

    cooler --version > ${software}.version.txt
    """
}

