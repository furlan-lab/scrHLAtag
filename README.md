<p align="center"><img src="img/Artboard1.png" alt="" width="400"></a></p>


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
**version 0.1.7** - 5/13/25 - IMGT/HLA is now version 3.60.0 - 2025-04; updated reference fastas to include 4-field alleles that differ at the mRNA level despite nomenclature

**version 0.1.6** - 4/21/25 - added qual string to fastq; fixed a bug that searched for nb tags in the aligned bam and printed long logfiles when it found none.  This had previously been fixed by writing nb in the readname.

**version 0.1.5** - 6/15/23 - added more checking for all components of a count before writing (i.e. will skip a read if cb, umi, or name is not present)

**version 0.1.4** - 6/13/23 - added a feature "hla_sep" which allows the user to define the "string separator" character - i.e. B\*48:01:01 where \* is the "string separator"

**version 0.1.3** - 6/12/23 - if no alleles file is supplied - now runs against all alleles in the database

**version 0.1.2** - 4/26/23 - returns nb tag (read more about nb tag here: https://isoseq.how/umi/isoseq-correct.html) in molecule_info file if present in input bam file

**version 0.1.1** - 4/11/23 - returns name of molecule in molecule_info and optionally sequence

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
1. Because scrHLAtag is designed for use with single cell RNA data, ***only up to 3-field HLA nomenclature should be used***.  HLA alleles which only vary in the 4th field from each other (i.e. A\*30:02:01:01, A\*30:02:01:02, A\*30:02:01:03, etc) have been removed from the reference and only the \*:\*:\*:01 allele was retained in the reference.  Code showing how this was accomplished, can be found in `scripts/generation_refs.sh`.  See https://hla.alleles.org/pages/nomenclature/naming_alleles/ for more details about HLA nomenclature.
2. The HLA reference file (from Nolan registry) is included in this program and does not need to be supplied during invocation.  *Current version: 3.60.0 - 2025-04*
3. Because the `*` character is not fasta friendly, the `|` character is used in our fasta reference files as a separator.  The alleles_file should still contain \*, however.


##### Invocation

**To run scrHLAtag simply type:**

`scrHLAtag -b BAMFILE -a ALLELES_FILE -o OUTFOLDER`



##### Help menu

```sh
scrHLAtag 0.1.6
copyright Scott Furlan
scrHLAtag is a command line tool for aligning and counting long-read sequence specific for HLA alleles in single cell
libraries

USAGE:
    scrHLAtag [FLAGS] [OPTIONS] --bam <bam>

FLAGS:
        --help       Prints help information
    -s, --ret_seq    include this flag to return sequence in molecule_info text file
    -V, --version    Prints version information
    -v, --verbose    verbose

OPTIONS:
    -l, --level <align_level>       align to 'genome', 'transcriptome', or 'both'; default is 'both'
    -a, --alleles <alleles_file>    table of hla genes to search for (tsv file); if not supplied will run against all
                                    alleles in database
    -b, --bam <bam>                 input bam
    -c, --cb_tag <cb>               character to parse cell barcode; default = 'CB'
    -h, --hla_sep <hsep>            character to separate HLA allele names; default = '*'
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
molecules_info_gene.txt.gz           = a file listing alignment metrics for each molecule; columns: CB, nb_tag (if present), UMI, allele, start_pos, mapq, cigar, NM, AS, s1, de); optionally including sequence in the last column
molecules_info_mRNA.txt.gz           = a file listing alignment metrics for each molecule; columns: CB, nb_tag (if present), UMI, allele, start_pos, mapq, cigar, NM, AS, s1, de); optionally including sequence in the last column
```
See minimap2 manual (https://lh3.github.io/minimap2/minimap2.html) for a discussion of the molecule_info metrics

See Isoseq FAQ for a discussion of the nb tag (https://isoseq.how/umi/isoseq-correct.html)










