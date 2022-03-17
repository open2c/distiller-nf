
params = [:]

/*
 * Miscellaneous code for the pipeline
 */

pairsgz_decompress_command = 'bgzip -cd -@ 3'

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

boolean isSingleFile(object) {
    object instanceof Path
}

String getOutputDir(output_type) {
    new File(params.output.dirs.get(output_type, output_type)).getCanonicalPath()
}

Boolean needsDownloading(query) {
    return (
        (query instanceof String) && (
           query.startsWith('sra:')
           || query.startsWith('http://')
           || query.startsWith('https://')
           || query.startsWith('ftp://')
       )
   )
}

String checkLeftRightChunk(left_chunk_fname,right_chunk_fname) {
    // checks if the chunk index is the same
    // both for left and right chunks and returns
    // that chunk index:
    // left by design:  ${library}.${run}.*.1.fastq.gz
    // right by design: ${library}.${run}.*.2.fastq.gz
    left_chunk_idx  =  left_chunk_fname.toString().tokenize('.')[-4]
    right_chunk_idx = right_chunk_fname.toString().tokenize('.')[-4]
    leftness_of_chunk  =  left_chunk_fname.toString().tokenize('.')[-3]
    rightness_of_chunk = right_chunk_fname.toString().tokenize('.')[-3]
    // assertions (should never happen by design, by just in case):
    // sidedness should be different:
    assert leftness_of_chunk != rightness_of_chunk: ("Sidedness suffix of"
                                                    +"LEFT and RIGHT sides of"
                                                    +"fastq chunks should DIFFER!")
    assert left_chunk_idx == right_chunk_idx: ("Chunk index of"
                                            +"LEFT and RIGHT sides of"
                                            +"fastq chunks should be IDENTICAL!")
    // return Chunk index of fastq chunk:
    return left_chunk_idx
}


/*
 * Download and chunk fastqs.
 */


//def fastqDumpCmd(file_or_srr, library, run, srr_start=0, srr_end=-1, threads=1, use_custom_split=true ) {
def fastqDumpCmd(file_or_srr, library, run, srr_start, srr_end, threads, use_custom_split) {
    def srr_start_flag = (srr_start == 0) ? '' : (' --minSpotId ' + srr_start)
    def srr_end_flag = (srr_end == -1) ? '' : (' --maxSpotId ' + srr_end)
    def bgzip_threads = Math.max(1,((threads as int)-2).intdiv(2))
    def cmd = ''

    if (use_custom_split) {
        cmd = """
            #HOME=`readlink -e ./`
            fastq-dump ${file_or_srr} -Z --split-spot ${srr_start_flag} ${srr_end_flag} \
                        | pyfilesplit --lines 4 \
                            >(bgzip -c -@{bgzip_threads} > ${library}.${run}.1.fastq.gz) \
                            >(bgzip -c -@{bgzip_threads} > ${library}.${run}.2.fastq.gz) \
                            | cat """


    } else {
        cmd = """
            #HOME=`readlink -e ./`
            fastq-dump ${file_or_srr} --gzip --split-spot --split-3 ${srr_start_flag} ${srr_end_flag}
            mv *_1.fastq.gz ${library}.${run}.1.fastq.gz
            mv *_2.fastq.gz ${library}.${run}.2.fastq.gz
        """
    }

    return cmd
}


def sraDownloadTruncateCmd(sra_query, library, run, truncate_fastq_reads, chunksize, threads, use_custom_split) {

//    def args.getOrDefault('truncate_fastq_reads', 0)
//    def args.getOrDefault('chunksize', 0_
//    def args.getOrDefault('threads', 1)
//    def args.getOrDefault('use_custom_split', true)

    def cmd = ""

    def srr = ( sra_query =~ /(?:SRR|ERR)\d+/ )[0]
    def srrnum = srr.substring(3)

    def srr_start = 0
    def srr_end = -1

    if ( sra_query.contains('start=') ) {
        srr_start = ( sra_query =~ /start=(\d+)/ )[0][1]
    }

    if ( truncate_fastq_reads ) {
        srr_end = srr_start + truncate_fastq_reads
    } else if ( sra_query.contains('end=') ) {
        srr_end = (sra_query =~ /end=(\d+)/)[0][1]
    }

    if ((srr_start > 0) || (srr_end != -1)) {
        cmd = """
            ${fastqDumpCmd(srr, library, run, srr_start, srr_end, threads, use_custom_split)}
            if [ -d ./ncbi ]; then rm -Rf ./ncbi; fi
        """
    }
    else {
        cmd = """
	    URL=\$(srapath ${srr})
            if wget --spider \$URL 2>/dev/null; then
		echo 'Dowloading from ', \$URL
                wget \$URL -qO ${srr}.sra
                ${fastqDumpCmd(srr+'.sra', library, run, 0, -1, threads, use_custom_split)}
                rm ${srr}.sra
            else
                echo 'Cannot wget an sra, fall back to fastq-dump'
                ${fastqDumpCmd(srr, library, run, 0, -1, threads, use_custom_split)}
                if [ -d ./ncbi ]; then rm -Rf ./ncbi; fi
            fi
        """
    }

    def chunk_lines = 4 * chunksize
    def split_bgzip_threads = Math.max(1, (threads as int)-1)
    if ( (truncate_fastq_reads == 0) && (chunk_lines > 0) ) {
        for (side in 1..2) {
            cmd += """
                zcat ${library}.${run}.${side}.fastq.gz | \
                     split -l ${chunk_lines} --numeric-suffixes=1 \
                     --filter 'bgzip -c -@ ${split_bgzip_threads} > \$FILE.${side}.fastq.gz' - \
                     ${library}.${run}.
                rm ${library}.${run}.${side}.fastq.gz
            """
        }
    } else {
        for (side in 1..2) {
            cmd += """
                mv ${library}.${run}.${side}.fastq.gz ${library}.${run}.0.${side}.fastq.gz
            """
        }
    }

    return cmd
}


//String fastqDownloadTruncateCmd(query, library, run, side, truncate_fastq_reads=0, chunksize=0, threads=2) {
String fastqDownloadTruncateCmd(query, library, run, side, truncate_fastq_reads, chunksize, threads) {
    def cmd = ''

    def truncate_lines = 4 * truncate_fastq_reads
    def chunk_lines = 4 * chunksize
    def bgzip_threads = Math.max(1, (threads as int)-1)

    if (truncate_lines > 0) {
        cmd = """head -n ${truncate_lines} < <(wget ${query} -O - | gunzip -cd )\
                 | bgzip -c -@ ${bgzip_threads} \
                 > ${library}.${run}.0.${side}.fastq.gz
              """
    } else if (chunk_lines > 0) {
        cmd = """wget ${query} -O - \
                 | gunzip -cd \
                 | split -l ${chunk_lines} --numeric-suffixes=1 \
                 --filter 'bgzip -c -@ ${bgzip_threads} > \$FILE.${side}.fastq.gz' - \
                 ${library}.${run}.
             """
    } else {
        cmd = "wget ${query} -O ${library}.${run}.0.${side}.fastq.gz"
    }

    return cmd
}


//String fastqLocalTruncateChunkCmd(path, library, run, side, truncate_fastq_reads=0, chunksize=0, threads=1) {
String fastqLocalTruncateChunkCmd(path, library, run, side, truncate_fastq_reads, chunksize, threads) {
    def cmd = ""

    def truncate_lines = 4 * truncate_fastq_reads
    def chunk_lines = 4 * chunksize
    def bgzip_threads = Math.max(1, (threads as int)-1)

    if (truncate_lines > 0) {
        cmd = """head -n ${truncate_lines} < <( zcat ${path} ) \
                 | bgzip -c -@ ${bgzip_threads} \
                 >  ${library}.${run}.0.${side}.fastq.gz
        """
    } else if (chunk_lines > 0) {
        cmd = """
            zcat ${path} | \
            split -l ${chunk_lines} --numeric-suffixes=1 \
            --filter 'bgzip -c -@ ${bgzip_threads} > \$FILE.${side}.fastq.gz' - \
            ${library}.${run}.
        """
    } else {
        // this line should never be reached, b/c local files should only
        // be processed if truncation or chunking is requested.
        cmd = "mv ${path} ${library}.${run}.0.${side}.fastq.gz"
    }

    return cmd
}


//
//  Useful utility functions used in nf-core DSL2 module files
//

//
// Extract name of software tool from process name using $task.process
//
def getSoftwareName(task_process) {
    return task_process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()
}

//
// Function to initialise default values and to generate a Groovy Map of available options for nf-core modules
//
def initOptions(Map args) {
    def Map options = [:]
    options.args            = args.args ?: ''
    options.args2           = args.args2 ?: ''
    options.args3           = args.args3 ?: ''
    options.publish_by_meta = args.publish_by_meta ?: []
    options.publish_dir     = args.publish_dir ?: ''
    options.publish_files   = args.publish_files
    options.suffix          = args.suffix ?: ''
    return options
}

//
// Tidy up and join elements of a list to return a path string
//
def getPathFromList(path_list) {
    def paths = path_list.findAll { item -> !item?.trim().isEmpty() }      // Remove empty entries
    paths     = paths.collect { it.trim().replaceAll("^[/]+|[/]+\$", "") } // Trim whitespace and trailing slashes
    return paths.join('/')
}

//
// Function to save/publish module results
//
def saveFiles(Map args) {
    if (!args.filename.endsWith('.version.txt')) {
        def ioptions  = initOptions(args.options)
        def path_list = [ ioptions.publish_dir ?: args.publish_dir ]
        if (ioptions.publish_by_meta) {
            def key_list = ioptions.publish_by_meta instanceof List ? ioptions.publish_by_meta : args.publish_by_meta
            for (key in key_list) {
                if (args.meta && key instanceof String) {
                    def path = key
                    if (args.meta.containsKey(key)) {
                        path = args.meta[key] instanceof Boolean ? "${key}_${args.meta[key]}".toString() : args.meta[key]
                    }
                    path = path instanceof String ? path : ''
                    path_list.add(path)
                }
            }
        }
        if (ioptions.publish_files instanceof Map) {
            for (ext in ioptions.publish_files) {
                if (args.filename.endsWith(ext.key)) {
                    def ext_list = path_list.collect()
                    ext_list.add(ext.value)
                    return "${getPathFromList(ext_list)}/$args.filename"
                }
            }
        } else if (ioptions.publish_files == null) {
            return "${getPathFromList(path_list)}/$args.filename"
        }
    }
}