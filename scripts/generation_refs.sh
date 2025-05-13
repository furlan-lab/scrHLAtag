---
output: html_document
editor_options: 
  chunk_output_type: console
---
## download IMGT-HLA database
```bash
cd ~/develop
git clone https://github.com/ANHIG/IMGTHLA.git
cd IMGTHLA
cat release_version.txt #version 3.60.0 date: 2025-04-09
cat Allelelist.txt
unzip hla_gen.fasta.zip #for some reason this file is zipped in this version
y
```

## create a list of 3-fied HLA tyoes appropriate for RNA sequencing (mRNA level)

```r
R
library(colorout)
library("Biostrings")
library(tidyr)
library(Rgb)
library(seqinr)
ROOT<-"~/develop/IMGTHLA"
fa<-seqinr::read.fasta(file.path(ROOT, "hla_nuc.fasta"), forceDNAtolower = FALSE)
ns<-lapply(fa, function(l) {attr(l, "Annot")})
length(ns) #42579
n<-strsplit(as.character(unlist(ns)), " ") %>% sapply("[[", 2)
tk<-which(unlist(lapply(strsplit(n, ":"), function(l) length(l)<4)))
tf<-which(unlist(lapply(strsplit(n, ":"), function(l) length(l)==4)))
o<-which(unlist(lapply(strsplit(n, ":"), function(l) length(l)>4)))

n[tf] %>% length #8255 that need to be collapsed

names<-unique(unlist(lapply(strsplit(n[tf], ":"), function(l) paste(l[1], l[2], l[3], sep=":"))))
cum<-vector("list", length(names))
names(cum)<-names
for(entry in n[tf]){
    #check to see if there is an entry already
    cat<-unlist(lapply(strsplit(entry, ":"), function(l) paste(l[1], l[2], l[3], sep=":")))
    cum[[cat]]<-c(cum[[cat]], entry)
}

length(cum) #768 categories

tg<-as.character(sapply(cum, "[[", 1))
ntk<-c(tk, match(tg, n))
ntk<-ntk[order(ntk)]
length(ntk) #35092

fa2<-fa[ntk]
nst<-lapply(fa2, function(l) {attr(l, "Annot")})
nt<-strsplit(as.character(unlist(nst)), " ") %>% sapply("[[", 2)
length(nt)==length(unique(nt))
newnames<-gsub(">", "", as.character(unlist(nst)))
newnamesl<-strsplit(newnames, " ")
#change 4 field to 3
f3f<- sapply(newnamesl, "[[", 2)
f3f_w<-strsplit(f3f, ":")
fixed_3f<-lapply(f3f_w, function(entry) {
    if(length(entry)<4){
        paste0(entry, collapse=":")
        }else{
            paste0(entry[1:3], collapse=":")
        }
    }) %>% unlist()
newnames<-paste(fixed_3f, sapply(newnamesl, "[[", 1), sapply(newnamesl, "[[", 3), sapply(newnamesl, "[[", 4), sep=" ")
#ok
newnames<-gsub("\\*", "|", newnames)
# newnames<-gsub("\\:", "_", newnames)
newnames
length(fa2)

## the fairly common DRB4*01:03N (a.k.a DRB4*01:03:01:02N) in association with DRB1*07s is missing from fa2
##let's investigate HLA nomunclature of DRB4*01:03N vs DRB4*01:03:01:02N
ref_align_og <- readDNAStringSet(file.path(ROOT, "hla_nuc.fasta")) # mRNA alignment
allele_3field <- sapply(
  strsplit(sapply(strsplit(ref_align_og@ranges@NAMES, " "), '[', 2), ":"), 
  function(x) paste(na.omit(x[1:3]), collapse = ":")
)
allele_3field <- unique(allele_3field)
allele_3field <- paste0(" ", allele_3field)
allele_3field <- gsub("\\*", "\\\\*", allele_3field)
length(allele_3field) == length(fa2) # should be true

## see if we are correctly naming the 4-field alleles by their 3-field name in 'nuc' space (because the 4th field is supposedly invisible in RNA)
correctly_named <- pbmcapply::pbmclapply(allele_3field, function(x) {
  x <- paste0(x, "(?=[: ])")
  x <- grep(x, ref_align_og@ranges@NAMES, value = TRUE, perl  = TRUE) # value=T to return the actual strings instead of index, perl=T to support look-ahead and look-behind constructs
  test <- ref_align_og[ref_align_og@ranges@NAMES %in% x]
  if (length(test) > 1) {
    algn <- DECIPHER::AlignSeqs(test, verbose = F)
    cm <- Biostrings::consensusMatrix(algn)
    as.they.should <- all(colSums(cm > 0) == 1)
  } else {
    as.they.should <- TRUE
  }
  as.they.should
}, mc.cores = parallel::detectCores()) %>% unlist()

allele_3field <- data.frame(allele_3f = allele_3field, correctly_named = correctly_named)
table(allele_3field$correctly_named) # 35067/35092 are correctly reduced to their 3-field counterpart, but not for 25 of the 35092
allele_3false <- allele_3field[which(allele_3field$correctly_named == F),] # this includes DRB4*01:03:01!!

alleles.to.reintroduce <- pbmcapply::pbmclapply(allele_3false$allele_3f, function(x) {
  test <- ref_align_og[ref_align_og@ranges@NAMES %like% x]
  #grab the names of the outlier seqs:
  algn <- DECIPHER::AlignSeqs(test, verbose = F)
  cm <- Biostrings::consensusMatrix(algn, baseOnly = FALSE)
  cons_vec <- apply(cm, 2, function(col) {
    # in case of tie, which.max() just picks the first line for a quick consensus
    rownames(cm)[which.max(col)]
  })
  consensus <- paste0(cons_vec, collapse = "")
  #turn aligned sequences into a character matrix
  mat <- do.call(rbind, strsplit(as.character(algn), split = ""))
  rownames(mat) <- names(algn)
  #flag any row with a mismatch
  bad_nomenclature <- rownames(mat)[apply(mat != matrix(cons_vec,
                                                        nrow = nrow(mat),
                                                        ncol = length(cons_vec),
                                                        byrow = TRUE),
                                          1, any)]
  bad_seqs <- apply(mat[bad_nomenclature, , drop = FALSE], 1, paste0, collapse = "")
  bad_nomenclature <- bad_nomenclature[!duplicated(bad_seqs)] #duplicated() flags any sequence matching one we’ve already seen and keeps the first of each.
}, mc.cores = parallel::detectCores()) %>% unlist()

# reintroduce those alleles, as de facto they don't have 3-field equivalents
alleles.to.reintroduce <- stringr::str_split_fixed(alleles.to.reintroduce, " ", 4) %>% as.data.frame()
fa3 <- fa[which(names(fa) %in% alleles.to.reintroduce$V1)]
newnames2 <- unlist(lapply(fa3, function(l) {attr(l, "Annot")}))
newnames2 <- gsub("\\*", "|", newnames2)
newnames2 <- gsub(">", "", as.character(unlist(newnames2)))
newnames2 <- stringr::str_split_fixed(newnames2, " ", 4) %>% as.data.frame()
newnames2 <- newnames2[, c(2,1,3,4)]
newnames2 <- apply(newnames2, 1, paste, collapse = " ")

fa3 <- c(fa2, fa3)
newnames2 <- c(newnames, newnames2)

seqinr::write.fasta(fa3, file.out="~/develop/scrHLAtag/data/HLA_DB_3field_mRNA.fa", names=newnames2)
#35092+49 seqs
all_alleles_mRNA <- gsub("\\|", "*", newnames2)
all_alleles_mRNA <- stringr::str_split_fixed(all_alleles_mRNA, " ", 2)[ ,1]
```

## create a list of 3-fied HLA tyoes appropriate for RNA sequencing (GENE level)

```r
R
library(colorout)
library("Biostrings")
library(tidyr)
library(Rgb)
library(seqinr)
ROOT<-"~/develop/IMGTHLA"
fa<-seqinr::read.fasta(file.path(ROOT, "hla_gen.fasta"), forceDNAtolower = FALSE)
ns<-unlist(lapply(fa, function(l) {attr(l, "Annot")}))
length(ns) #24355
n<-strsplit(as.character(unlist(ns)), " ") %>% sapply("[[", 2)
tk<-which(unlist(lapply(strsplit(n, ":"), function(l) length(l)<4)))
tf<-which(unlist(lapply(strsplit(n, ":"), function(l) length(l)==4)))
o<-which(unlist(lapply(strsplit(n, ":"), function(l) length(l)>4)))

n[tf] %>% length #8243 that need to be collapsed

names<-unique(unlist(lapply(strsplit(n[tf], ":"), function(l) paste(l[1], l[2], l[3], sep=":"))))
cum<-vector("list", length(names))
names(cum)<-names
for(entry in n[tf]){
    #check to see if there is an entry already
    cat<-unlist(lapply(strsplit(entry, ":"), function(l) paste(l[1], l[2], l[3], sep=":")))
    cum[[cat]]<-c(cum[[cat]], entry)
}

length(cum) #767 categories

tg<-as.character(sapply(cum, "[[", 1))
ntk<-c(tk, match(tg, n))
ntk<-ntk[order(ntk)]

fa2<-fa[ntk]
length(fa2)
nst<-lapply(fa2, function(l) {attr(l, "Annot")})
nt<-strsplit(as.character(unlist(nst)), " ") %>% sapply("[[", 2)
length(nt)==length(unique(nt))
newnames<-gsub(">", "", as.character(unlist(nst)))
newnamesl<-strsplit(newnames, " ")
#change 4 field to 3
f3f<- sapply(newnamesl, "[[", 2)
f3f_w<-strsplit(f3f, ":")
fixed_3f<-lapply(f3f_w, function(entry) {
    if(length(entry)<4){
        paste0(entry, collapse=":")
        }else{
            paste0(entry[1:3], collapse=":")
        }
    }) %>% unlist()
newnames<-paste(fixed_3f, sapply(newnamesl, "[[", 1), sapply(newnamesl, "[[", 3), sapply(newnamesl, "[[", 4), sep=" ")
#ok
newnames<-gsub("\\*", "|", newnames)
# newnames<-gsub("\\:", "_", newnames)
newnames
length(fa2)

## reintroduce the same list of alleles that erroneously don't have 3-field equivalents
fa3 <- fa[which(names(fa) %in% alleles.to.reintroduce$V1)]
newnames2 <- unlist(lapply(fa3, function(l) {attr(l, "Annot")}))
newnames2 <- gsub("\\*", "|", newnames2)
newnames2 <- gsub(">", "", as.character(unlist(newnames2)))
newnames2 <- stringr::str_split_fixed(newnames2, " ", 4) %>% as.data.frame()
newnames2 <- newnames2[, c(2,1,3,4)]
newnames2 <- apply(newnames2, 1, paste, collapse = " ")


fa3 <- c(fa2, fa3)
newnames2 <- c(newnames, newnames2)

seqinr::write.fasta(fa3, file.out="~/develop/scrHLAtag/data/HLA_DB_3field_gene.fa", names=newnames2)
#16879+37 seqs
all_alleles_gene <- gsub("\\|", "*", newnames2)
all_alleles_gene <- stringr::str_split_fixed(all_alleles_gene, " ", 2)[ ,1]

#if you still have `all_alleles_mRNA` from previous chunck:
all_alleles <- c(all_alleles_gene, all_alleles_mRNA) %>% unique() %>% sort()
write(all_alleles, file.path("~/develop/scrHLAtag/data/all_alleles.tsv"))
```



## check seq lengths
```r
R
library(colorout)
library("Biostrings")
library(tidyr)
library(Rgb)
library(seqinr)
ROOT<-"~/develop/scrHLAtag/data/"
fa_mRNA<-seqinr::read.fasta(file.path(ROOT, "HLA_DB_3field_mRNA.fa"), forceDNAtolower = FALSE)
fa_gene<-seqinr::read.fasta(file.path(ROOT, "HLA_DB_3field_gene.fa"), forceDNAtolower = FALSE)

length(fa_mRNA[[1]])
attr(fa_mRNA[[1]], "Annot")
length(fa_gene[[1]])
attr(fa_gene[[1]], "Annot")

ns_mRNA<-lapply(fa_mRNA, function(l) {attr(l, "Annot")}) %>% unlist %>% as.character() %>% strsplit(" ") %>% sapply("[[", 2)
head(ns_mRNA)
tail(ns_mRNA)

ns_gene<-lapply(fa_gene, function(l) {attr(l, "Annot")}) %>% unlist %>% as.character() %>% strsplit(" ") %>% sapply("[[", 2)
head(ns_gene)
tail(ns_gene)

#all(ns_gene %in% ns_mRNA)

```


```bash
gzip ~/develop/scrHLAtag/data/HLA_DB_3field_gene.fa
```


```bash
gzip ~/develop/scrHLAtag/data/HLA_DB_3field_mRNA.fa
```
