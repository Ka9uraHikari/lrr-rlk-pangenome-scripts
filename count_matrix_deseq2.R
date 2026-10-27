suppressPackageStartupMessages(library(DESeq2))
a <- commandArgs(trailingOnly=TRUE)
x <- as.matrix(read.delim(a[1],row.names=1,check.names=FALSE));md<-read.delim(a[2],check.names=FALSE)
dir.create(a[3],recursive=TRUE,showWarnings=FALSE)
for(trt in unique(md$treatment[md$experiment=='disease'])){
 m<-md[md$experiment=='disease' & md$treatment==trt,]
 m$timepoint<-factor(m$timepoint,levels=c('0d','0h',setdiff(unique(m$timepoint),c('0d','0h'))));m$timepoint<-droplevels(m$timepoint)
 rownames(m)<-m$run_accession_CRR
 y<-round(x[,m$run_accession_CRR]);y<-y[rowSums(y>=10)>=3,,drop=FALSE]
 d<-DESeqDataSetFromMatrix(y,m,~timepoint);d<-DESeq(d,quiet=TRUE)
 ref<-levels(m$timepoint)[1]
 for(tm in levels(m$timepoint)[-1]){
  r<-as.data.frame(results(d,contrast=c('timepoint',tm,ref)));r<-cbind(Gene=rownames(r),r)
  write.table(r,file.path(a[3],paste0(trt,'_',tm,'.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
 }
 cat(trt,'completed\n')
}
m<-md[md$experiment=='abiotic',];rownames(m)<-m$run_accession_CRR
print(unique(m[,c('treatment','timepoint','condition_code')]))
controls<-unique(m$condition_code[!m$treatment %in% c('Cold','Heat','Salt','cold','heat','salt')]);stopifnot(length(controls)==1)
m$group<-relevel(factor(m$condition_code),controls)
y<-round(x[,m$run_accession_CRR]);y<-y[rowSums(y)>0,,drop=FALSE]
d<-DESeqDataSetFromMatrix(y,m,~group);d<-DESeq(d,quiet=TRUE)
for(g in setdiff(levels(m$group),controls)){
 r<-as.data.frame(results(d,contrast=c('group',g,controls)));r<-cbind(Gene=rownames(r),r)
 write.table(r,file.path(a[3],paste0(g,'.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
}
writeLines(c(R.version.string,paste('DESeq2',packageVersion('DESeq2'))),file.path(a[3],'session.txt'))
