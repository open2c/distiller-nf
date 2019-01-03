### 0.3.1 (2019-01-03) ###

* bugfixes in cooler, fastqc task and apt.list

### 0.3.0 (2019-01-02) ###

A major update of the pipeline and the config syntax.
* Simplify and speed up the pipeline by merging several groups of processes.
* Greatly reduce the storage requirements down to ~2x of the fastq.gz size.
* Introduce custom pair filtering during binning.
* Make .mcools with custom resolutions.
* Speed up the .sra->.fastq.gz conversion.
* Rename multiple config options.
* Switch from pbgzip to native bgzip multithreading.
* Allow users to change the location of the temp files.
* Post the list of the packages inside the docker on github.

### 0.2.0 (2018-11-18) ###

* report mapq in .pairs

### 0.1.1 (2018-11-09) ###

* remove .sra files when downloading data from SRA.
* use process selectors in configs.

### 0.1.0 (2018-08-04) ###

* Distiller-env: switch to pairtools v0.2.0 and cooler v0.7.10.
