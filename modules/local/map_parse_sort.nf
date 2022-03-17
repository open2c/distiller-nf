// Import generic module functions
include { getSoftwareName; initOptions; saveFiles; getOutputDir } from './functions'
include { sraDownloadTruncateCmd; fastqDownloadTruncateCmd } from './functions'

params.options = [:]
options        = initOptions(params.options)
directory = getOutputDir('mapped_parsed_sorted_chunks')

ASSEMBLY_NAME = params['input'].genome.assembly_name // TODO: move to the parameters dictionary, and below:
switch(params.compression_format) {
    case 'gz':
        suffix = 'gz'
        decompress_command = 'bgzip -cd -@ 3'
        break
    case 'lz4':
        suffix = 'lz4'
        decompress_command = 'lz4c -cd'
        break
    default:
        suffix = 'gz'
        decompress_command = 'bgzip -cd -@ 3'
        break
}

process MAP_PARSE_SORT{
    tag "library:${library} run:${run}"
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
    tuple val(library), val(run), val(chunk), file(fastq1), file(fastq2)
    tuple val(bwa_index_base), file(bwa_index_files)
    file(chrom_sizes)

    output:
    tuple val(library), val(run), val(chunk),
        path("${library}.${run}.${ASSEMBLY_NAME}.${chunk}.pairsam.${suffix}"),
        path("${library}.${run}.${ASSEMBLY_NAME}.${chunk}.bam"), emit: output
    path  "*.version.txt"         , emit: version


    script:
    def software = getSoftwareName(task.process)

    def mapping_options = params['map'].getOrDefault('mapping_options','')
    def trim_options = params['map'].getOrDefault('trim_options','')

    def dropsam_flag = params['parse'].getOrDefault('make_pairsam','false').toBoolean() ? '' : '--drop-sam'
    def dropreadid_flag = params['parse'].getOrDefault('drop_readid','false').toBoolean() ? '--drop-readid' : ''
    def dropseq_flag = params['parse'].getOrDefault('drop_seq','false').toBoolean() ? '--drop-seq' : ''
    def keep_unparsed_bams_command = (
        params['parse'].getOrDefault('keep_unparsed_bams','false').toBoolean() ?
        "| tee >(samtools view -bS > ${library}.${run}.${ASSEMBLY_NAME}.${chunk}.bam)" : "" )
    def parsing_options = params['parse'].getOrDefault('parsing_options','')
    def bwa_threads = (task.cpus as int)
    def sorting_threads = (task.cpus as int)

    def mapping_command = (
        trim_options ?
        "fastp ${trim_options} \
        --json ${library}.${run}.${ASSEMBLY_NAME}.${chunk}.fastp.json \
        --html ${library}.${run}.${ASSEMBLY_NAME}.${chunk}.fastp.html \
        -i ${fastq1} -I ${fastq2} --stdout | \
        bwa mem -p -t ${bwa_threads} ${mapping_options} -SP ${bwa_index_base} \
        - ${keep_unparsed_bams_command}" : \
        \
        "bwa mem -t ${bwa_threads} ${mapping_options} -SP ${bwa_index_base} \
        ${fastq1} ${fastq2} ${keep_unparsed_bams_command}"
        )


    """
    TASK_TMP_DIR=\$(mktemp -d -p ${task.distillerTmpDir} distiller.tmp.XXXXXXXXXX)
    touch ${library}.${run}.${ASSEMBLY_NAME}.${chunk}.bam

    ${mapping_command} \
    | pairtools parse ${dropsam_flag} ${dropreadid_flag} ${dropseq_flag} \
      ${parsing_options} \
      -c ${chrom_sizes} \
      | pairtools sort --nproc ${sorting_threads} \
                     -o ${library}.${run}.${ASSEMBLY_NAME}.${chunk}.pairsam.${suffix} \
                     --tmpdir \$TASK_TMP_DIR \
      | cat

    rm -rf \$TASK_TMP_DIR

    #bwa 2>&1 | head -n 3 | tail -n 1 &> ${software}.version.txt # TODO: fix error message
    pairtools --version >> ${software}.version.txt
    """
}

