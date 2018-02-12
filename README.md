## HGSC bcftools Plugins

This repository contains a few little, simple plugins used at the HGSC for manipulating VCF/BCF files.

To use them, simply copy the .c files into the plugins folder of your [bcftools](https://github.com/samtools/bcftools) 1.6 installation and follow the instructions for building bcftools (i.e., run 'make').

A description of how we might chain these plugins with other bcftools commands in a typical post-processing pipeline can be found in FILTERING_AND_FORMATTING.md.

A description of the fields output by the summary stats plugins is can be found in SUMMARY_STATS.md. 

### Copyright

Copyright 2018 Baylor College of Medicine Human Genome Sequencing Center

### License

See LICENSE.txt.

