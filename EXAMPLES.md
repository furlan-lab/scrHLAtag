<img width="75" alt="image" src="scrHLAtag.png">

# Example usage with test data

### Build and install scrHLAtag
```sh
loc=~/develop/scrHLAtag # or location where you have cloned the repository
cd $loc
cargo build --release && cp target/release/scrHLAtag ~/.local/bin
```

### Run on small dataset
```sh
loc=~/develop/scrHLAtag
scrHLAtag -v -b $loc/data/test.bam -a $loc/data/testhla.tsv -o out -l transcriptome -s
```

You should get
```txt
Writing log file: 'out/scrHLAtag.log'
Running with 1 thread(s)!
Checking command line programs!

Making fastq from BAM: '/Users/sfurlan/develop/scrHLAtag/data/test.bam'
Total reads processed: 196
Reads with errors: 0
Matching HLA alleles with known transcriptome reference!
Making partial reference: 'out/align_transcriptome.fa'
Read HLA reference file: '/Users/sfurlan/develop/scrHLAtag/data/HLA_DB_3field_mRNA.fa.gz'
        Found: A*24:02:01
        Found: A*30:02:01
        Found: A*33:03:01
        Found: B*42:01:01
        Found: B*48:01:01
        Found: B*57:02:01
        Found: C*08:06
        Found: C*17:01:02
        Found: C*18:02:01
Could not find: DMG*04:10:01
        Found: DPA1*01:03:01
        Found: DPA1*02:02:02
        Found: DPA1*03:01:01
        Found: DPB1*01:01:01
        Found: DPB1*04:02:01
        Found: DPB1*105:01:01
        Found: DQA1*03:03:01
        Found: DQA1*04:01:01
        Found: DQB1*03:01:01
        Found: DQB1*04:02:01
        Found: DRB1*03:02:01
        Found: DRB1*04:01:01
        Found: DRB3*01:62:01
Aligning reads using minimap2:
[M::mm_idx_gen::0.001*3.00] collected minimizers
[M::mm_idx_gen::0.001*2.48] sorted minimizers
[M::main::0.001*2.42] loaded/built the index for 22 target sequence(s)
[M::mm_mapopt_update::0.001*2.34] mid_occ = 50
[M::mm_idx_stat] kmer size: 19; skip: 19; is_hpc: 0; #seq: 22
[M::mm_idx_stat::0.001*2.26] distinct minimizers: 1040 (52.79% are singletons); average occurrences: 1.806; average spacing: 10.382; total length: 19497
[M::worker_pipeline::0.061*1.02] mapped 196 sequences
[M::main] Version: 2.24-r1150-dirty
[M::main] CMD: minimap2 --cs=long --secondary=no -x map-hifi -Q --MD -a -t 1 -o out/Aligned_mm2_transcriptome.sam out/align_transcriptome.fa out/fastq.fq.gz
[M::main] Real time: 0.061 sec; CPU: 0.062 sec; Peak RSS: 0.004 GB

Minimap2 complete; Running samtools sort

Samtools sort complete; Running samtools view

Samtools view complete; Running samtools index

Counting reads:
Total reads processed: 196
Reads with errors: 0
Unmapped reads: 27
Mapped reads: 169

Cleaning up!

Done!!!
```

### Run on larger dataset
```sh
loc=~/develop/scrHLAtag
scrHLAtag -v -b $loc/data/full_dedup.bam -a $loc/data/testhla.tsv -o out
```

You should get
```txt
Writing log file: 'out/scrHLAtag.log'
Running with 1 thread(s)!
Checking command line programs!

Making fastq from BAM: '/Users/sfurlan/develop/scrHLAtag/data/full_dedup.bam'
Total reads processed: 261690
Reads with errors: 0
Matching HLA alleles with known genome reference!
Making partial reference: 'out/align_genome.fa'
Read HLA reference file: '/Users/sfurlan/develop/scrHLAtag/data/HLA_DB_3field_gene.fa.gz'
        Found: A*24:02:01
        Found: A*30:02:01
        Found: A*33:03:01
        Found: B*42:01:01
        Found: B*48:01:01
        Found: B*57:02:01
        Found: C*08:06
Could not find: C*17:01:02
        Found: C*18:02:01
Could not find: DMG*04:10:01
        Found: DPA1*01:03:01
        Found: DPA1*02:02:02
        Found: DPA1*03:01:01
        Found: DPB1*01:01:01
        Found: DPB1*04:02:01
        Found: DPB1*105:01:01
        Found: DQA1*03:03:01
        Found: DQA1*04:01:01
        Found: DQB1*03:01:01
        Found: DQB1*04:02:01
        Found: DRB1*03:02:01
        Found: DRB1*04:01:01
        Found: DRB3*01:62:01
Aligning reads using minimap2:
[M::mm_idx_gen::0.003*1.45] collected minimizers
[M::mm_idx_gen::0.004*1.27] sorted minimizers
[M::main::0.004*1.27] loaded/built the index for 21 target sequence(s)
[M::mm_mapopt_update::0.004*1.25] mid_occ = 50
[M::mm_idx_stat] kmer size: 19; skip: 19; is_hpc: 0; #seq: 21
[M::mm_idx_stat::0.005*1.24] distinct minimizers: 10080 (63.05% are singletons); average occurrences: 1.565; average spacing: 10.029; total length: 158202
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/184702|BARCODE=ACAGCTATCAAAGACA_TCTTGTTAGC'. Continue anyway.
[M::worker_pipeline::12.204*1.00] mapped 14611 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/92898|BARCODE=ACGATACGTTGAGTTC_AGATTATGCA'. Continue anyway.
[M::worker_pipeline::16.675*1.00] mapped 5347 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/142954|BARCODE=ACGCCGAAGCTTATCG_CGTCAAGTGG'. Continue anyway.
[M::worker_pipeline::18.150*1.00] mapped 1716 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/212528|BARCODE=AGCAGCCCACCTTGTC_GTATGGGTCC'. Continue anyway.
[M::worker_pipeline::27.774*1.00] mapped 11421 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/153673|BARCODE=ATCCACCAGCGTTGCC_ACCTCAACAC'. Continue anyway.
[M::worker_pipeline::39.696*1.00] mapped 13909 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/73150|BARCODE=ATTGGACTCATAAAGG_TTTATGACGG'. Continue anyway.
[M::worker_pipeline::45.032*1.00] mapped 6418 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/236409|BARCODE=CAGCAGCAGTCGTTTG_ATTTGCGAGC'. Continue anyway.
[M::worker_pipeline::57.069*1.00] mapped 14366 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/82969|BARCODE=CAGGTGCAGAATGTGT_GTGACGAGAA'. Continue anyway.
[M::worker_pipeline::59.588*1.00] mapped 3043 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/130102|BARCODE=CCACTACGTGCAGACA_TATTCCCCAA'. Continue anyway.
[M::worker_pipeline::68.129*1.00] mapped 10400 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/136505|BARCODE=CCAGCGAAGCTAGGCA_TCGGAACTAG'. Continue anyway.
[M::worker_pipeline::68.278*1.00] mapped 189 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/6728|BARCODE=CCGTGGAAGACGACGT_CTACCAGTGT'. Continue anyway.
[M::worker_pipeline::71.987*1.00] mapped 4330 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/138996|BARCODE=CGCGGTATCATTTGGG_TAATTGATCC'. Continue anyway.
[M::worker_pipeline::83.318*1.00] mapped 13663 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/243789|BARCODE=CGTTCTGCAGGACCCT_GGCAAACGCC'. Continue anyway.
[M::worker_pipeline::92.744*1.00] mapped 11052 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/171803|BARCODE=CTACGTCAGAGCTTCT_CAATTATGCC'. Continue anyway.
[M::worker_pipeline::97.572*1.00] mapped 5650 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/244761|BARCODE=CTAGTGACAGCTGTAT_TTTAACTGCC'. Continue anyway.
[M::worker_pipeline::99.697*1.00] mapped 2519 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/225366|BARCODE=GAAGCAGTCGAGAGCA_AGCAACTGGG'. Continue anyway.
[M::worker_pipeline::114.811*1.00] mapped 17787 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/96883|BARCODE=GACCAATAGCAGATCG_TATGCTGACG'. Continue anyway.
[M::worker_pipeline::117.439*1.00] mapped 3164 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/168365|BARCODE=GACGTTAAGCCCAATT_CGTATACACT'. Continue anyway.
[M::worker_pipeline::119.466*1.00] mapped 2371 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/213724|BARCODE=GACTACAGTTGATTGC_CGGTTGTCGT'. Continue anyway.
[M::worker_pipeline::120.843*1.00] mapped 1563 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/97177|BARCODE=GATGAGGGTCATGCCG_GGTCTCCATC'. Continue anyway.
[M::worker_pipeline::125.467*1.00] mapped 5395 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/148640|BARCODE=GCGCAGTCAACACCCG_ACTGTCGCAG'. Continue anyway.
[M::worker_pipeline::134.151*1.00] mapped 10047 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/184091|BARCODE=GGACGTCCACGCGAAA_GCTAAAGGCT'. Continue anyway.
[M::worker_pipeline::142.611*1.00] mapped 10148 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/199008|BARCODE=GTTCATTAGATATACG_TCACTTCCAG'. Continue anyway.
[M::worker_pipeline::166.347*1.00] mapped 28008 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/243571|BARCODE=GTTCTCGCAAACTGCT_GGGTGGCGTGT'. Continue anyway.
[M::worker_pipeline::167.618*1.00] mapped 1430 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/41639|BARCODE=TACTTGTTCCGTTGTC_CGGTCCCGAC'. Continue anyway.
[M::worker_pipeline::176.093*1.00] mapped 10088 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/251905|BARCODE=TCACAAGGTTTGTTGG_GCAGGAACAT'. Continue anyway.
[M::worker_pipeline::182.634*1.00] mapped 7710 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/94023|BARCODE=TCATTACAGCGTCTAT_CCGCCATTAC'. Continue anyway.
[M::worker_pipeline::185.206*1.00] mapped 2984 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/196395|BARCODE=TCGCGTTCAATGTTGC_GGAGTGCCAC'. Continue anyway.
[M::worker_pipeline::188.061*1.00] mapped 3439 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/235311|BARCODE=TGGACGCTCCGCAAGC_TTCAGGTACG'. Continue anyway.
[M::worker_pipeline::203.060*1.00] mapped 17803 sequences
[M::worker_pipeline::220.754*1.00] mapped 21081 sequences
[M::main] Version: 2.24-r1150-dirty
[M::main] CMD: minimap2 --cs=long --secondary=no -x map-hifi -Q --MD -a -t 1 -o out/Aligned_mm2_genome.sam out/align_genome.fa out/fastq.fq.gz
[M::main] Real time: 220.755 sec; CPU: 220.383 sec; Peak RSS: 0.129 GB

Minimap2 complete; Running samtools sort

Samtools sort complete; Running samtools view

Samtools view complete; Running samtools index

Counting reads:
Total reads processed: 289488
Reads with errors: 0
Unmapped reads: 39902
Mapped reads: 249586

Cleaning up!

Matching HLA alleles with known transcriptome reference!
Making partial reference: 'out/align_transcriptome.fa'
Read HLA reference file: '/Users/sfurlan/develop/scrHLAtag/data/HLA_DB_3field_mRNA.fa.gz'
        Found: A*24:02:01
        Found: A*30:02:01
        Found: A*33:03:01
        Found: B*42:01:01
        Found: B*48:01:01
        Found: B*57:02:01
        Found: C*08:06
        Found: C*17:01:02
        Found: C*18:02:01
Could not find: DMG*04:10:01
        Found: DPA1*01:03:01
        Found: DPA1*02:02:02
        Found: DPA1*03:01:01
        Found: DPB1*01:01:01
        Found: DPB1*04:02:01
        Found: DPB1*105:01:01
        Found: DQA1*03:03:01
        Found: DQA1*04:01:01
        Found: DQB1*03:01:01
        Found: DQB1*04:02:01
        Found: DRB1*03:02:01
        Found: DRB1*04:01:01
        Found: DRB3*01:62:01
Aligning reads using minimap2:
[M::mm_idx_gen::0.001*3.11] collected minimizers
[M::mm_idx_gen::0.001*2.52] sorted minimizers
[M::main::0.001*2.51] loaded/built the index for 22 target sequence(s)
[M::mm_mapopt_update::0.001*2.42] mid_occ = 50
[M::mm_idx_stat] kmer size: 19; skip: 19; is_hpc: 0; #seq: 22
[M::mm_idx_stat::0.001*2.34] distinct minimizers: 1040 (52.79% are singletons); average occurrences: 1.806; average spacing: 10.382; total length: 19497
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/184702|BARCODE=ACAGCTATCAAAGACA_TCTTGTTAGC'. Continue anyway.
[M::worker_pipeline::4.253*1.00] mapped 14611 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/92898|BARCODE=ACGATACGTTGAGTTC_AGATTATGCA'. Continue anyway.
[M::worker_pipeline::5.824*1.00] mapped 5347 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/142954|BARCODE=ACGCCGAAGCTTATCG_CGTCAAGTGG'. Continue anyway.
[M::worker_pipeline::6.381*1.00] mapped 1716 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/212528|BARCODE=AGCAGCCCACCTTGTC_GTATGGGTCC'. Continue anyway.
[M::worker_pipeline::9.692*1.00] mapped 11421 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/153673|BARCODE=ATCCACCAGCGTTGCC_ACCTCAACAC'. Continue anyway.
[M::worker_pipeline::13.817*1.00] mapped 13909 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/73150|BARCODE=ATTGGACTCATAAAGG_TTTATGACGG'. Continue anyway.
[M::worker_pipeline::15.699*1.00] mapped 6418 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/236409|BARCODE=CAGCAGCAGTCGTTTG_ATTTGCGAGC'. Continue anyway.
[M::worker_pipeline::19.769*1.00] mapped 14366 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/82969|BARCODE=CAGGTGCAGAATGTGT_GTGACGAGAA'. Continue anyway.
[M::worker_pipeline::20.621*1.00] mapped 3043 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/130102|BARCODE=CCACTACGTGCAGACA_TATTCCCCAA'. Continue anyway.
[M::worker_pipeline::23.552*1.00] mapped 10400 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/136505|BARCODE=CCAGCGAAGCTAGGCA_TCGGAACTAG'. Continue anyway.
[M::worker_pipeline::23.606*1.00] mapped 189 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/6728|BARCODE=CCGTGGAAGACGACGT_CTACCAGTGT'. Continue anyway.
[M::worker_pipeline::24.833*1.00] mapped 4330 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/138996|BARCODE=CGCGGTATCATTTGGG_TAATTGATCC'. Continue anyway.
[M::worker_pipeline::28.747*1.00] mapped 13663 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/243789|BARCODE=CGTTCTGCAGGACCCT_GGCAAACGCC'. Continue anyway.
[M::worker_pipeline::31.991*1.00] mapped 11052 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/171803|BARCODE=CTACGTCAGAGCTTCT_CAATTATGCC'. Continue anyway.
[M::worker_pipeline::33.558*1.00] mapped 5650 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/244761|BARCODE=CTAGTGACAGCTGTAT_TTTAACTGCC'. Continue anyway.
[M::worker_pipeline::34.241*1.00] mapped 2519 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/225366|BARCODE=GAAGCAGTCGAGAGCA_AGCAACTGGG'. Continue anyway.
[M::worker_pipeline::39.462*1.00] mapped 17787 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/96883|BARCODE=GACCAATAGCAGATCG_TATGCTGACG'. Continue anyway.
[M::worker_pipeline::40.307*1.00] mapped 3164 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/168365|BARCODE=GACGTTAAGCCCAATT_CGTATACACT'. Continue anyway.
[M::worker_pipeline::41.000*1.00] mapped 2371 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/213724|BARCODE=GACTACAGTTGATTGC_CGGTTGTCGT'. Continue anyway.
[M::worker_pipeline::41.508*1.00] mapped 1563 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/97177|BARCODE=GATGAGGGTCATGCCG_GGTCTCCATC'. Continue anyway.
[M::worker_pipeline::43.126*1.00] mapped 5395 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/148640|BARCODE=GCGCAGTCAACACCCG_ACTGTCGCAG'. Continue anyway.
[M::worker_pipeline::46.054*1.00] mapped 10047 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/184091|BARCODE=GGACGTCCACGCGAAA_GCTAAAGGCT'. Continue anyway.
[M::worker_pipeline::48.922*1.00] mapped 10148 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/199008|BARCODE=GTTCATTAGATATACG_TCACTTCCAG'. Continue anyway.
[M::worker_pipeline::56.909*1.00] mapped 28008 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/243571|BARCODE=GTTCTCGCAAACTGCT_GGGTGGCGTGT'. Continue anyway.
[M::worker_pipeline::57.369*1.00] mapped 1430 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/41639|BARCODE=TACTTGTTCCGTTGTC_CGGTCCCGAC'. Continue anyway.
[M::worker_pipeline::60.285*1.00] mapped 10088 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/251905|BARCODE=TCACAAGGTTTGTTGG_GCAGGAACAT'. Continue anyway.
[M::worker_pipeline::62.570*1.00] mapped 7710 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/94023|BARCODE=TCATTACAGCGTCTAT_CCGCCATTAC'. Continue anyway.
[M::worker_pipeline::63.474*1.00] mapped 2984 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/196395|BARCODE=TCGCGTTCAATGTTGC_GGAGTGCCAC'. Continue anyway.
[M::worker_pipeline::64.418*1.00] mapped 3439 sequences
[WARNING] failed to parse the FASTA/FASTQ record next to 'molecule/235311|BARCODE=TGGACGCTCCGCAAGC_TTCAGGTACG'. Continue anyway.
[M::worker_pipeline::69.541*1.00] mapped 17803 sequences
[M::worker_pipeline::75.400*1.00] mapped 21081 sequences
[M::main] Version: 2.24-r1150-dirty
[M::main] CMD: minimap2 --cs=long --secondary=no -x map-hifi -Q --MD -a -t 1 -o out/Aligned_mm2_transcriptome.sam out/align_transcriptome.fa out/fastq.fq.gz
[M::main] Real time: 75.401 sec; CPU: 75.328 sec; Peak RSS: 0.099 GB

Minimap2 complete; Running samtools sort

Samtools sort complete; Running samtools view

Samtools view complete; Running samtools index

Counting reads:
Total reads processed: 261687
Reads with errors: 0
Unmapped reads: 40357
Mapped reads: 221330

Cleaning up!

Done!!!
```



## Thank you for trying out scrHLAtag!

