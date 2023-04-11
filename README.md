<img width="150" alt="image" src="scrHLAtag.png">


#                       scrHLAtag
Pipeline for processing scrHLA typing data


### Overview

scrHLAtag is a command line tool written in Rust for processing long-read (PacBio) sequencing data generated from 10X Genomics libraries.  Using scrHLAtag, you can supply relevant HLA alleles, and obtain a counts file of the number of reads that map to your HLA query and their associated cell-barcodes and umis.  scrHLAtag takes as input a BAM file that has been processed using the single cell IsoSeq3 pipeline through the dedup stage.

### Installation

scrHLAtag is written in Rust.  To install the rust compiler go to https://www.rust-lang.org/tools/install.  scrHLAtag requires two additional tools be available on the command line, minimap2 (https://github.com/lh3/minimap2) and samtools (https://samtools.github.io). 

To install scrHLAtag:
1. clone the repository by typing `https://github.com/furlan-lab/scrHLAtag.git` from the location you want to build from
2. enter the cloned repo by typing `cd scrHLAtag`
3. build by typing `cargo build --release`
4. the build process will create a self contained binary executable file in `targets/release` directory called `scrHLAtag`
5. move this binary elsewhere if desired (ideally somewhere referenced by your PATH environment variable - e.g. `~/.local/bin`)

### Updates

**version 0.1.0** - 4/9/23 - first stable version.  Implements both transcriptomic and genomic options

### Usage

##### Alleles File.

First users will generate a simple 'alleles_file' that lists the HLA alleles to be counted (txt.file). The names of the alleles should be taken from the Anthony Nolan registry (https://www.ebi.ac.uk/ipd/imgt/hla/).  The alleles_file should look something like this:

```sh
A*30:02:01
A*24:02:01
B*57:02:01
B*48:01:01
C*18:02:01
C*08:06
DRB1*03:02:01
DRB1*04:01:01
DRB3*01:62:01
DQB1*04:02:01
```
**Important Notes:**
1. Because scrHLAtag is designed for use with single cell RNA data, ***only up to 3-field HLA nomenclature should be used***.  HLA alleles which only vary in the 4th field from each other (i.e. A\*30:02:01:01, A\*30:02:01:02, A\*30:02:01:03, etc) have been removed from the reference and only the \*:\*:\*:01 allele was retained in the reference.  Code showing how this was accomplished, can be found in `scripts/generation_refs.sh`.  See https://hla.alleles.org/nomenclature/naming.html for more details about HLA nomenclature.
2. The HLA reference file (from Nolan registry) is included in this program and does not need to be supplied during invocation.  *Current version: 3.51.0 - 2023-01*
3. Because the `*` character is not fasta friendly, the `|` character is used in our fasta reference files as a separator.  The alleles_file should still contain \*, however.


##### Invocation

**To run scrHLAtag simply type:**

`scHLAtag -b BAMFILE -a ALLELES_FILE -o OUTFOLDER`



##### Help menu

```sh
scrHLAtag 0.1.0
copyright Scott Furlan
scrHLAtag is a command line tool for aligning and counting long-read sequence specific for HLA alleles in single cell
libraries

USAGE:
    scrHLAtag [FLAGS] [OPTIONS] --alleles <alleles_file> --bam <bam>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information
    -v, --verbose    verbose

OPTIONS:
    -l, --level <align_level>       align to 'genome', 'transcriptome', or 'both'; default is 'both'
    -a, --alleles <alleles_file>    table of hla genes to search for (tsv file)
    -b, --bam <bam>                 input bam
    -c, --cb_tag <cb>               character to parse cell barcode; default = 'CB'
    -o, --out <output_folder>       folder for output; default 'out'
    -t, --threads <threads>         threads
    -u, --umi_tag <umi>             character to parse umi; default = 'XM'
```
 
### Output
 
scrHLAtag will output some combination of the following files (depending on alignment level):
```sh
Aligned_mm2_sorted_gene.bam          = minimap2 output bam file aligned to genome, sorted by read name 
Aligned_mm2_sorted_mRNA.bam          = minimap2 output bam file aligned to transcriptome, sorted by read name 
Aligned_mm2_sorted_gene.bam.bai      = index for above bam
Aligned_mm2_sorted_mRNA.bam.bai      = index for above bam
counts_gene.txt.gz                   = counts file; columns: CB, UMI, allele, read_count
counts._mRNAtxt.gz                   = counts file; columns: CB, UMI, allele, read_count
align_gene.fa                        = fasta reference file used for minimap2 alignment
align_mRNA.fa                        = fasta reference file used for minimap2 alignment
molecules_info_gene.txt.gz           = a file listing alignment metrics for each molecule; columns: CB, UMI, allele, start_pos, mapq, cigar, NM, AS, s1, de)
molecules_info_mRNA.txt.gz           = a file listing alignment metrics for each molecule; columns: CB, UMI, allele, start_pos, mapq, cigar, NM, AS, s1, de)
```
See minimap2 manual (https://lh3.github.io/minimap2/minimap2.html) for a discussion of the molecule_info metrics










