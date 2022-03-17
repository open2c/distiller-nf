// Import generic module functions
include { initOptions; getSoftwareName; getOutputDir } from './functions'

params.options = [:]
options        = initOptions(params.options)
directory      = getOutputDir('fastqc')

process FASTQC {
    tag "library:${library} run:${run} chunk:${chunk} side:${side}"
    label 'process_low'
    publishDir "${directory}", mode: params.publish_dir_mode

    conda (params.enable_conda ? "bioconda::sra-tools>=2.8.1 bioconda::pbgzip" : null)
//        if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
//            container "https://depot.galaxyproject.org/singularity/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2880dd9d8ad0a7b221d4eacda9a818e92983128d-0"
//        } else {
//            container "quay.io/biocontainers/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2880dd9d8ad0a7b221d4eacda9a818e92983128d-0"
//        }

    input:
    tuple val(library), val(run), val(chunk), val(side), file(fastq)

    output:
    tuple val(library), val(run), val(chunk), val(side),
        path("${library}.${run}.${chunk}.${side}_fastqc.html"),
        path("${library}.${run}.${chunk}.${side}_fastqc.zip"), emit: output
    path  "*.version.txt"         , emit: version


    script:
    def software = getSoftwareName(task.process)

    """
    TASK_TMP_DIR=\$(mktemp -d -p ${task.distillerTmpDir} distiller.tmp.XXXXXXXXXX)
    ln -s \"\$(readlink -f ${fastq})\" \$TASK_TMP_DIR/${library}.${run}.${chunk}.${side}.fastq.gz
    fastqc --threads ${task.cpus} -o ./ -f fastq \$TASK_TMP_DIR/${library}.${run}.${chunk}.${side}.fastq.gz
    rm -r \$TASK_TMP_DIR

    fastqc --version >> ${software}.version.txt
    """
}

