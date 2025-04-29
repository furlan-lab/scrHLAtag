

# Run scrHLAtag on bams with qual data for eventual variant calling
```sh
sFH2
cd ~/develop/scrHLAtag/data
cp /fh/fast/furlan_s/user/skanaan/scrHLAtyping/Manuscript/2nd_Draft/geo_upload/scrHLAtyping/AML8_34_hifi_reads.bam .
cp /fh/fast/furlan_s/user/skanaan/scrHLAtyping/Manuscript/2nd_Draft/geo_upload/scrHLAtyping/AML8_BM_hifi_reads.bam .
MM
micromamba activate isoseq4
export ref_dir=/fh/fast/furlan_s/grp/refs/long_read_refs/pacbio
export primers_5p=$ref_dir/5p_10x_primers.fasta
export cbc_include=/fh/fast/furlan_s/grp/refs/long_read_refs/pacbio/737K-august-2016.txt
export hg38=/fh/fast/furlan_s/grp/refs/GRCh38/GRCh38.p13.genome.fa
export annotation=/fh/fast/furlan_s/grp/refs/long_read_refs/annotation.sorted.gtf
export pA=$ref_dir/polyA.list
export cage=$ref_dir/refTSS_v3.1_human_coordinate.hg38.bed
export threads=30
export BAM=AML8_34_hifi_reads.bam
export BAM=AML8_BM_hifi_reads.bam
ls -asl $BAM
export bam=$BAM
export out_subdir=AML8
export prefix=$out_subdir/AML8_BM
echo $prefix
mkdir -p $out_subdir
sbatch -n 1 -c 30 -p campus-new -M gizmo --mem-per-cpu=21000MB --wrap='lima \
    --log-level INFO \
    --log-file $out_subdir/lima.log \
    --isoseq \
    -j $threads \
    $bam \
    $primers_5p \
    $prefix.bam &&
isoseq tag \
    --log-level INFO \
    --log-file $out_subdir/tag.log \
    -j $threads \
    --design 16B-10U-T \
    $prefix.5p--3p.bam \
    $prefix.tagged.bam &&
isoseq refine \
    --log-level INFO \
    --log-file $out_subdir/refine.log \
    -j $threads \
    --require-polya \
    $prefix.tagged.bam \
    $primers_5p \
    $prefix.refined.bam &&
isoseq correct \
    --log-level INFO \
    --log-file $out_subdir/correct.log \
    -j $threads \
    --barcodes $cbc_include \
    $prefix.refined.bam \
    $prefix.corrected.bam'

cat > aml8_alleles.tsv << EOL
A*11:01:01
A*24:02:01
B*07:06:01
B*15:02:01
C*07:02:01
C*08:01:01
DMA*01:01:01
DMA*01:02:01
DMB*01:01:01
DMB*01:02:01
DOA*01:01:02
DOA*01:01:04
DOA*01:10
DOB*01:01:01
DPA1*01:03:01
DPA1*02:01:01
DPA1*02:02:02
DPB1*01:01:01
DPB1*13:01:01
DPB1*21:01
DPB2*02:01
DPB2*03:01:01
DQA1*01:03:01
DQA1*01:04:01
DQA1*06:01:01
DQB1*03:01:01
DQB1*06:01:01
DQB1*06:02:01
DRA*01:01:01
DRA*01:02:02
DRB1*11:05
DRB1*12:02:01
DRB1*15:01:01
DRB3*02:02:01
DRB3*03:01:03
DRB5*01:01:01
E*01:01:01
E*01:03:01
E*01:03:02
F*01:01:01
F*01:01:02
H*02:07:01
MICA*008:01:03
MICA*008:04:02
MICA*019:01:01
MICA*027:01:01
MICB*005:02:01
MICB*005:03:01
MICB*008:01:01
MICB*014:01:01
TAP1*01:01:01
TAP1*03:01
TAP2*02:01:02
EOL


ml minimap2/2.26-GCCcore-12.3.0
ml SAMtools
cargo build --release && cp ../target/release/scrHLAtag ..

sbatch -n 1 -c 30 -p campus-new -M gizmo --mem-per-cpu=21000MB --wrap='../scrHLAtag -v -b AML8/AML8_34.corrected.bam -a aml8_alleles.tsv -o AML8 -l transcriptome -s -t 30'

../scrHLAtag -v -b AML8_34.corrected.bam -a aml8_alleles.tsv -o AML8_34 -l transcriptome -s -t 16

../scrHLAtag -v -b AML8_BM.corrected.bam -a aml8_alleles.tsv -o AML8_BM -l transcriptome -s -t 16
```

```R
#pull up reads in IGV
#grab some and copy here
R
library(DECIPHER)
library(magrittr)
library(seqinr)
library(msa)
seqs <- read.fasta("align_mRNA.fa")
A11_ref <- DNAStringSet(paste0(seqs[[1]], collapse=""))
A24_ref <- DNAStringSet(paste0(seqs[[2]], collapse=""))
mA11 <- "ATGGCCGTCATGGCGCCCCGAACCCTCCTCCTGCTACTCTCGGGGGCCCTGGCCCTGACCCAGACCTGGGCGGGCTCCCACTCCATGAGGTATTTCTACACCTCCGTGTCCCGGCCCGGCCGCGGGGAGCCCCGCTTCATCGCCGTGGGCTACGTGGACGACACGCAGTTCGTGCGGTTCGACAGCGACGCCGCGAGCCAGAGGATGGAGCCGCGGGCGCCGTGGATAGAGCAGGAGGGGCCGGAGTATTGGGACCAGGAGACACGGAATGTGAAGGCCCAGTCACAGACTGACCGAGTGGACCTGGGGACCCTGCGCGGCTACTACAACCAGAGCGAGGACGGTTCTCACACCATCCAGATAATGTATGGCTGCGACGTGGGGCCGGACGGGCGCTTCCTCCGCGGGTACCGGCAGGACGCCTACGACGGCAAGGATTACATCGCCCTGAACGAGGACCTGCGCTCTTGGACCGCGGCGGACATGGCAGCTCAGATCACCAAGCGCAAGTGGGAGGCGGCCCATGCGGCGGAGCAGCAGAGAGCCTACCTGGAGGGCCGGTGCGTGGAGTGGCTCCGCAGATACCTGGAGAACGGGAAGGAGACGCTGCAGCGCACGGACCCCCCCCAAGACACATATGACCCACCACCCCATCTCTGACCATGAGGCCACCCTGAGGTGCTGGGCCCTGGGCTTCTACCCTGCGGAGATCACACTGACCTGGCAGCGGGATGGGGAGGATCAGACCCAGGACACGGAGCTCGTGGAGACCAGGCCTGCAGGGGATGGAACCTTCCAGAAGTGGGCGGCTGTGGTGGTGCCTTCTGGAGAGGAGCAGAGATACACCTGCCATGTGCAGCATGAGGGTCTGCCCAAGCCCCTCACCCTGAGATGGGAGCTGTCTTCCCAGCCCACCATCCCCATCGTGGGCATCATTGCTGGCCTGGTTCTCCTTGGAGCTGTGATCACTGGAGCTGTGGTCGCTGCCGTGATGTGGAGGAGGAAGAGCTCAGGTGGAGAAGGGGTGAAAGATAGAAAAGGAGGGAGTTACACTCAGGCTGCAAGCAGTGACAGTGCCCAGGGCTCTGATGTGTCTCTCACAGCTTGTAAAGTGTGAGACAGCTGCCTTGTGTGGGACTGAGAGGCAGGAGTTGTTCCTGCCCTTCCCTTTGTGACTTGAAGAACCCTGACTTTGTTTCTGCAAAGGCACCTGCATGTGTCTGTGTTCGTGTAGGCATAATGTGAGGAGGTGGGGAGACCACCCCACCCCGATGTCCACCATGACCCTCTTCCCACGCTGACCTGTGCTCCCTCTCCAATCATCTTTCCTGTTCCAGAGAGGTGGGGCTGAGGTGTCTCCATCTCTGTCTCAACTTCATGGTGCACTGAGCTGTAACTTCTTCCTTCCCTATTAAAATTAGAACCTG"
mmA11<-"CTCTCGGAGGCCCTGGCCCTTACCCAGACCTGGGCGGGCTCCCACTCCTTGAAGTATTTCCACACTTCCGTGTCCCGGCCCGGCCGCGGGGAGCCCCGCTTCATCTCTGTGGGCTACGTGGACGACACCCAGTTCGTGCGCTTCGACAACGACGCCGCGAGTCCGAGGATGGTGCCGCGGGCGCCGTGGATGGAGCAGGAGGGGTCAGAGTATTGGGACCGGGAGACACGGAGCGCCAGGGACACCGCACAGATTTTCCGAGTGAACCTGCGGACGCTGCGCGGCTACTACAATCAGAGCGAGGCCGGGTCTCACACCCTGCAGTGGATGCATGGCTGCGAGCTGGGGCCCGACGGGCGCTTCCTCCGCGGGTATGAACAGTTCGCCTACGACGGCAAGGATTATCTCACCCTGAATGAGGACCTGCGCTCCTGGACCGCGGTGGACACGGCGGCTCAGATCTCCGAGCAAAAGTCAAATGATGCCTCTGAGGCGGAGCACCAGAGAGCCTACCTGGAAGACACATGCGTGGAGTGGCTCCACAAATACCTGGAGAAGGGGAAGGAGACGCTGCTTCACCTGGAGCCCCCAAAGACACACGTGACTCACCACCCCATCTCTGACCATGAGGCCACCCTGAGGTGCTGGGCCCTGGGCTTCTACCCTGCGGAGATCACACTGACCTGGCAGCAGGATGGGGAGGGCCATACCCAGGACACGGAGCTCGTGGAGACCAGGCCTGCAGGGGATGGAACCTTCCAGAAGTGGGCAGCTGTGGTGGTGCCTTCTGGAGAGGAGCAGAGATACACGTGCCATGTGCAGCATGAGGGGCTACCCGAGCCCGTCACCCTGAGATGGAAGCCGGCTTCCCAGCCCACCATCCCCATCGTGGGCATCATTGCTGGCCTGGTTCTCCTTGGATCTGTGGTCTCTGGAGCTGTGGTTGCTGCTGTGATATGGAGGAAGAAGAGCTCAGGTGGAAAAGGAGGGAGCTACTCTAAGGCTGAGTGGAGCGACAGTGCCCAGGGGTCTGAGTCTCACAGCTTGTAAAGCCTGAGACAGCTGCCTTGTGTGCGACTGAGATGCACAGCTGCCTTGTGTGCGACTGAGATGCAGGATTTCCTCACGCCTCCCCTATGTGTCTTAGGGGACTCTGGCTTCTCTTTTTGCAAGGGCCTCTGAATCTGTCTGTGTCCCTGTTAGCACAATGTGAGGAGGTAGAGAAACAGTCCACCTCTGTGTCTACCATGACCCCCTTCCTCACACTGACCTGTGTTCCTTCCCTGTTCTCTTTTCTATTAAAAATAAGAACCTGGGCAGAGTGCGGCAGCTCATGCCTGTAATCCCAGCACTTAGGGAGGCCGAGGAGGGCAGATCACGAGGTCAGGAGATCGAAACCATCCTGGCTAACACGGTGAAACCCCGTCTCTACTAAAAAATACAAAAAATTAGCTGGGCGCAGAGGCACGGGCCTGTAGTCCCAGCTACTCAGGAGGCGGAGGCAGGAGAATGGCGTCAACCCGGGAGGCGGAGGTTGCAGTGAGCCAGGATTGTGCGACTGCACTCCAGCCTGGGTGACAGGGTGAAACGCCATCTC"
reads<-c(DNAStringSet(mmA11), DNAStringSet(mA11), A11_ref, A24_ref)
names(reads)<-c("MM", "M", "A11ref", "A24ref")
A <-msa(reads)
BrowseSeqs(A@unmasked)

#mmA11 blasts as HLA-E; other clear outliers in IGV match to HLA-F
```

## longcallr doesn't work well.  lets try deepvariant

```sh

BIN_VERSION=1.8.0
cd /fh/fast/furlan_s/grp/sifs
export SINGULARITY_TMPDIR=/fh/fast/furlan_s/grp/sifs/tmp_dir
singularity pull docker://google/deepvariant:"${BIN_VERSION}-gpu"

sFH2
ml SAMtools/1.11-GCC-10.2.0
ml Singularity/3.5.3
export BIN_VERSION=1.8.0
export INPUT_DIR=/home/sfurlan/develop/scrHLAtag/data/AML8
export OUTPUT_DIR=/home/sfurlan/develop/scrHLAtag/data/AML8/deepvariant

cd ${INPUT_DIR}
export THREADS="1"
export REF_FASTA=align_mRNA.fa
samtools faidx $REF_FASTA
export SINGULARITY_BINDPATH="/fh/scratch,/fh/fast,/shared,/fh/working"
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${OUTPUT_DIR}/logs"
mkdir -p "${OUTPUT_DIR}/intermediate_results_dir"
export TMPDIR="${OUTPUT_DIR}/tmp_dir"
mkdir -p "${TMPDIR}"

export IBAM="Aligned_mm2_sorted_mRNA.bam"

sbatch -n 1 -c ${THREADS} --gpus=1 -p campus-new --mem-per-cpu=8GB --wrap='
    singularity exec --bind "${INPUT_DIR}":"/input","${OUTPUT_DIR}":"/output",/usr/lib/locale/:/usr/lib/locale/,"${TMPDIR},/tmp" \
    /fh/fast/furlan_s/grp/sifs/deepvariant_1.8.0-gpu.sif \
    run_deepvariant \
    --model_type="MASSEQ" \
    --ref=/input/${REF_FASTA} \
    --reads="/input/${IBAM}" \
    --output_vcf="/output/deepvariant_calls.vcf.gz" \
    --output_gvcf="/output/deepvariant_calls.g.vcf.gz" \
    --num_shards ${THREADS} \
    --make_examples_extra_args="vsc_min_count_indels=1,vsc_min_fraction_indels=0.01,realign_reads=false,min_mapping_quality=1,min_base_quality=1" \
    --postprocess_variants_extra_args="debug_output_all_candidates=ALT" \
    --intermediate_results_dir="/output/intermediate_results_dir"'


```


```sh
cd AML8_34
../../addback -b Aligned_mm2_sorted_mRNA.bam
samtools index out.bam

cd ../AML8_BM
../../addback -b Aligned_mm2_sorted_mRNA.bam
samtools index out.bam

cat > variants.csv << EOL
seq,start,ref_nt,query_nt,name
A|11:01:01,1011,A,AGGTGGAGAAGGGGTGAA,HLA_indel1
A|11:01:01,1013,A,GTGGAGAAGGGGTGAA,HLA_indel2
EOL
sed -E 's/("([^"]*)")?,/\2\t/g' variants.csv > variants.tsv
cat variants.tsv

~/develop/mutCaller/mutcaller ALIGNED --help
~/develop/mutCaller/mutcaller ALIGNED -b out.bam -s variants.tsv
```

