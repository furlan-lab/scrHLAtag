# scrHLAtag
Pipeline for processing scrHLA typing data

<object data="https://github.com/furlan-lab/scrHLAtag/scrHLAtag.pdf" type="application/pdf" width="300px" height="500px">
    <embed src="[http://yoursite.com/the.pdf](https://github.com/furlan-lab/scrHLAtag/scrHLAtag.pdf)">
        <p>This browser does not support PDFs. Please download the PDF to view it: <a href="[http://yoursite.com/the.pdf](https://github.com/furlan-lab/scrHLAtag/scrHLAtag.pdf)">Download PDF</a>.</p>
    </embed>
</object>

### Overview

scrHLAtag is a command line tool written in Rust for processing long-read (PacBio) sequencing data generated from 10X Genomics libraries.  Using scrHLAtag, you can supply relevant HLA alleles, and obtain a counts file of the number of reads that map to your HLA query and their associated cell-barcodes and umis.  scrHLAtag takes as input a BAM file that has been processed using the single cell IsoSeq3 pipeline through the dedup stage.

### Installation

scrHLAtag is written in Rust.  To install the rust compiler go to https://www.rust-lang.org/tools/install.  scrHLAtag requires two additional tools be available on the command line, minimap2 (ref) and samtools (ref). 

To install scrHLAtag:
1. clone the repository by typing `https://github.com/furlan-lab/scrHLAtag.git` from the location you want to build from
2. build by entring `cargo build release`
3. the build process will create a self contained binary executable file in targets/release directory called `scrHLAtag`
4. move this binary elsewhere if desired (ideally somewhere referenced by your PATH environment variable - e.g. `~/.local/bin`)

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
