# Module 1: Quality Control & Assembly

## Requirements

The `assembly/pipeline.py` script requires that the `pandas` library is installed for whichever `python` is invoked by `#!/usr/bin/env python`.

The `assembly/pipeline_qc.sh` script requires the following conda environments:

 - `mash-2.0`
 - `fastqc-0.11.7`
 - `kraken2-2.0.7_beta`
 - `seqtk-1.3`

The `assembly/pipeline_assembly.sh` script requires the following conda environments:

 - `shovill-1.0.1`
 - `quast-4.6.3`
 - `busco-3.0.2`

The `assembly/pipeline_assembly_contaminant.sh` script requires the following conda environments:

 - `bbmap-38.22`
 - `shovill-1.0.1`
 - `quast-4.6.3`
 - `busco-3.0.2`

## Configuration

A configuration file located at `assembly/config.ini` Use this file to configure the install location of these scripts and the locations of several reference databases used. The config file is a simple key-value [ini](https://en.wikipedia.org/wiki/INI_file) format. eg:

```
[scripts]
script-path = /home/dfornika/code/cpo-pipeline/assembly

[databases]
mash-genomedb = /data/ref_databases/mash/refseq.genomes.k21s1000.msh
mash-plasmiddb = /data/ref_databases/mash/refseq.plasmid.k21s1000.msh
kraken2-genomedb = /data/ref_databases/kraken2/2018-09-20_standard
kraken2-plasmiddb = /data/ref_databases/kraken2/2018-09-20_plasmid
busco-db = /data/ref_databases/busco/enterobacteriales_odb9
```

Replace these paths with your system-specific paths.