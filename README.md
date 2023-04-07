# scrHLAtag
Pipeline for processing scrHLA typing data

### Overview

scrHLAtag is command line software meant to process long-read (PacBio) sequencing data generated from 10X Genomics 5 prime libraries.  Using scrHLAtag, you can supply relevant HLA alleles, and obtain a counts file of the number of reads that map to your HLA query and their associated cell-barcodes and umis.  scrHLAtag takes as input a BAM file that has been processed using the single cell IsoSeq3 pipeline through the dedup stage.

### Installation

scrHLAtag is written in Rust and requires two additional tools are available on the command line, minimap2 (ref) and samtools (ref). Once these dependencies are met, scrHLAtag repository is cloned 

### Usage

##### Invocation.

Fist users will generate a simple text files with the HLA alleles of interest. The names of the alleles should be taken from the Anthony Nolan registry (link)


##### Help menu

count reads bam files per cb and umi for scrHLA typing

USAGE:
    scrHLAtag [FLAGS] [OPTIONS] --bam <bam>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information
    -v, --verbose    verbose

OPTIONS:
    -a, --hlaalleles <alleles>    table of hla genes to search for (tsv)
    -b, --bam <bam>               input bam
    -c <cb>                       character to parse cell barcode; default = 'CB='
    -g, --genome <genome>         reference genome
    -o, --out <outfile>           folder for output; default out
    -t, --threads <threads>       threads
    -u <umi>                      character to parse umi; default = 'XM='
