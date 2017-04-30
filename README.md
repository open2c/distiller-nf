# distiller-nf

## A modular Hi-C mapping pipeline for reproducible data analysis.

The `distiller` pipeline aims to provide the following functionality:

- Align the sequences of Hi-C molecules to the reference genome
- Parse .sam alignment and form files with Hi-C pairs
- Filter PCR duplicates
- Aggregate pairs into binned matrices of Hi-C interactions

### Installation

Requirements:

- java 7/8
- [nextflow](https://www.nextflow.io/)
- docker (should be able to run w/o root privileges, 
[tutorial](https://www.digitalocean.com/community/tutorials/how-to-install-and-use-docker-on-ubuntu-16-04))

To setup a new project, execute the following line in the project folder:

```sh
$ nextflow clone mirnylab/distiller
```

This will download the distiller pipeline and the configuration files.

Then:
- configure the location of the input files and other project details
in project.yml 
- configure the hardware/system details in nextflow.config

Launch distiller as:

```sh
$ nextflow distiller.nf -params-file project.yml -C nextflow.config
```
