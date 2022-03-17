// Import generic module functions
include { getSoftwareName; initOptions; saveFiles; getOutputDir } from './functions'

params.options = [:]
options        = initOptions(params.options)
directory = getOutputDir('coolers_library')

MIN_RES = params['bin'].resolutions.collect { it as int }.min() // TODO: move to parameters dictionary
ASSEMBLY_NAME = params['input'].genome.assembly_name // TODO: move to the parameters dictionary
pairsgz_decompress_command = 'bgzip -cd -@ 3'

process BIN_ZOOM {
    tag "library:${library} filter:${filter_name}"
    label 'process_low'
    publishDir "${directory}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::sra-tools>=2.8.1 bioconda::pbgzip" : null)
//        if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
//            container "https://depot.galaxyproject.org/singularity/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2880dd9d8ad0a7b221d4eacda9a818e92983128d-0"
//        } else {
//            container "quay.io/biocontainers/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2880dd9d8ad0a7b221d4eacda9a818e92983128d-0"
//        }

    input:
    tuple val(filter_name), val(filter_expr)
    tuple val(library), file(pairs_lib)
    file(chrom_sizes)

    output:
    tuple val(library), val(filter_name),
        path("${library}.${ASSEMBLY_NAME}.${filter_name}.${MIN_RES}.cool"),
        path("${library}.${ASSEMBLY_NAME}.${filter_name}.${MIN_RES}.mcool"), emit: output
    path  "*.version.txt"         , emit: version


    script:
    def software = getSoftwareName(task.process)

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

    cooler --version > ${software}.version.txt
    """
}

