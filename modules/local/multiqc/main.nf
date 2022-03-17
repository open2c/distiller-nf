// git+https://github.com/open2c/MultiQC.git@pairtools-module

include { initOptions; getOutputDir } from '../functions'

params.options = [:]
options        = initOptions(params.options)
directory      = getOutputDir('mapped_parsed_sorted_chunks')

process MULTIQC {
    label 'process_medium'

//git+https://github.com/open2c/MultiQC.git@pairtools-module
//    conda (params.enable_conda ? 'bioconda::multiqc=1.12' : null)
//    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
//        'https://depot.galaxyproject.org/singularity/multiqc:1.12--pyhdfd78af_0' :
//        'quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0' }"

    publishDir "${directory}", mode: params.publish_dir_mode

    input:
    path multiqc_files

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    multiqc -f -o ./ $args .
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}
