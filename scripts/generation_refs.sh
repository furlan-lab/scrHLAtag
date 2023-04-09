## download IMGT-HLA database
```bash
cd ~/develop
git clone https://github.com/ANHIG/IMGTHLA.git #version 3.51.0 4/8/23
cd IMGTHLA
cat Allelelist.txt
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
fa<-read.fasta(file.path(ROOT, "hla_gen.fasta"), forceDNAtolower = FALSE)
ns<-unlist(lapply(fa, function(l) {attr(l, "Annot")}))
length(ns) #19184
n<-strsplit(as.character(unlist(ns)), " ") %>% sapply("[[", 2)
tk<-which(unlist(lapply(strsplit(n, ":"), function(l) length(l)<4)))
tf<-which(unlist(lapply(strsplit(n, ":"), function(l) length(l)==4)))
o<-which(unlist(lapply(strsplit(n, ":"), function(l) length(l)>4)))

n[tf] #6822 that need to be collapsed

names<-unique(unlist(lapply(strsplit(n[tf], ":"), function(l) paste(l[1], l[2], l[3], sep=":"))))
cum<-vector("list", length(names))
names(cum)<-names
for(entry in n[tf]){
    #check to see if there is an entry already
    cat<-unlist(lapply(strsplit(entry, ":"), function(l) paste(l[1], l[2], l[3], sep=":")))
    cum[[cat]]<-c(cum[[cat]], entry)
}

length(cum) #677 categories

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

write.fasta(fa2, file.out="~/develop/scrHLAtag/data/HLA_DB_3field_gene.fa", names=newnames)
#13039 seqs
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
fa<-read.fasta(file.path(ROOT, "hla_nuc.fasta"), forceDNAtolower = FALSE)
ns<-lapply(fa, function(l) {attr(l, "Annot")})
length(ns) #36620
n<-strsplit(as.character(unlist(ns)), " ") %>% sapply("[[", 2)
tk<-which(unlist(lapply(strsplit(n, ":"), function(l) length(l)<4)))
tf<-which(unlist(lapply(strsplit(n, ":"), function(l) length(l)==4)))
o<-which(unlist(lapply(strsplit(n, ":"), function(l) length(l)>4)))

n[tf] #6832 that need to be collapsed

names<-unique(unlist(lapply(strsplit(n[tf], ":"), function(l) paste(l[1], l[2], l[3], sep=":"))))
cum<-vector("list", length(names))
names(cum)<-names
for(entry in n[tf]){
    #check to see if there is an entry already
    cat<-unlist(lapply(strsplit(entry, ":"), function(l) paste(l[1], l[2], l[3], sep=":")))
    cum[[cat]]<-c(cum[[cat]], entry)
}

length(cum) #677 categories

tg<-as.character(sapply(cum, "[[", 1))
ntk<-c(tk, match(tg, n))
ntk<-ntk[order(ntk)]
length(ntk) #30465

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
write.fasta(fa2, file.out="~/develop/scrHLAtag/data/HLA_DB_3field_mRNA.fa", names=newnames)
#13039 seqs
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
fa_mRNA<-read.fasta(file.path(ROOT, "HLA_DB_3field_mRNA.fa"), forceDNAtolower = FALSE)
fa_gene<-read.fasta(file.path(ROOT, "HLA_DB_3field_gene.fa"), forceDNAtolower = FALSE)

length(fa_mRNA[[1]])
attr(fa_mRNA[[1]], "Annot")
length(fa_gene[[1]])
attr(fa_gene[[1]], "Annot")

ns_mRNA<-lapply(fa_mRNA, function(l) {attr(l, "Annot")}) %>% unlist %>% as.character() %>% strsplit(" ") %>% sapply("[[", 2)

ns_gene<-lapply(fa_gene, function(l) {attr(l, "Annot")}) %>% unlist %>% as.character() %>% strsplit(" ") %>% sapply("[[", 2)

#all(ns_gene %in% ns_mRNA)

```


```bash
gzip ~/develop/scrHLAtag/data/HLA_DB_3field_gene.fa
```


```bash
gzip ~/develop/scrHLAtag/data/HLA_DB_3field_mRNA.fa
```
