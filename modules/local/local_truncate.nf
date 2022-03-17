// Import generic module functions
include { getSoftwareName; initOptions; saveFiles; getOutputDir } from './functions'
include { fastqLocalTruncateChunkCmd } from './functions'

params.options = [:]
options        = initOptions(params.options)
directory = getOutputDir('processed_fastqs')

process LOCAL_TRUNCATE {
    tag "library:${library} run:${run}"
    label 'process_low'
    publishDir "${directory}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::pbgzip" : null)
//        if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
//            container "https://depot.galaxyproject.org/singularity/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2880dd9d8ad0a7b221d4eacda9a818e92983128d-0"
//        } else {
//            container "quay.io/biocontainers/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2880dd9d8ad0a7b221d4eacda9a818e92983128d-0"
//        }

    input:
    tuple val(library), val(run), val(query1), val(query2)

    output:
    tuple val(library), val(run),
        path("${library}.${run}.*.1.fastq.gz"), path("${library}.${run}.*.2.fastq.gz"), emit: output
    path  "*.version.txt"         , emit: version


    script:
    def software = getSoftwareName(task.process)

    def truncate_fastq_reads = params['input'].get('truncate_fastq_reads',0)
    def chunksize = params['map'].get('chunksize', 0)

    def truncate_chunk_cmd1 = fastqLocalTruncateChunkCmd(
        fastq1, library, run, 1, truncate_fastq_reads, chunksize, task.cpus)
    def truncate_chunk_cmd2 = fastqLocalTruncateChunkCmd(
        fastq2, library, run, 2, truncate_fastq_reads, chunksize, task.cpus)

    """
    ${truncate_chunk_cmd1}
    ${truncate_chunk_cmd2}

    bgzip --version | head -n 1 > ${software}.version.txt
    """
}

