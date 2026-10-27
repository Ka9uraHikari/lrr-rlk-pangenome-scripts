suppressPackageStartupMessages(library(DESeq2))
a<-commandArgs(trailingOnly=TRUE)
if(length(a)!=4) stop('Usage: Rscript condition_adjusted_correlations.R counts.tsv metadata.tsv markers.tsv results.tsv')
x<-as.matrix(read.delim(a[1],row.names=1,check.names=FALSE));md<-read.delim(a[2],check.names=FALSE);markers<-read.delim(a[3],check.names=FALSE)
candidates<-c(CR133='CsaV4_7G002666',CR132='CsaV4_7G002654',CR018='CsaV4_7G000704',SC013='CsaV4_2G003347',CR094='CsaV4_5G000526',CR042='CsaV4_2G002663',SC022='CsaV4_3G003305',CR085='CsaV4_4G002864',CR095='CsaV4_5G000549')
rows<-list()
for(stratum in c('disease','abiotic')){
 m<-md[md$experiment==stratum,];rownames(m)<-m$run_accession_CRR;m$group<-factor(m$condition_code)
 d<-DESeqDataSetFromMatrix(round(x[,m$run_accession_CRR]),m,~group);d<-estimateSizeFactors(d)
 targets<-unique(c(unname(candidates),markers$Gene));y<-t(log2(counts(d,normalized=TRUE)[targets,,drop=FALSE]+1))
 X<-model.matrix(~0+group,m);residual<-lm.fit(X,y)$residuals;n<-nrow(m);rk<-qr(X)$rank;df<-n-rk-1
 for(cand in names(candidates))for(i in seq_len(nrow(markers))){
  r<-cor(residual[,candidates[cand]],residual[,markers$Gene[i]])
  p<-2*pt(-abs(r)*sqrt(df/(1-r*r)),df)
  r0<-cor(y[,candidates[cand]],y[,markers$Gene[i]]);p0<-2*pt(-abs(r0)*sqrt((n-2)/(1-r0*r0)),n-2)
  rows[[length(rows)+1]]<-data.frame(Candidate=cand,Candidate_gene=candidates[cand],Marker=markers$Marker[i],Marker_gene=markers$Gene[i],Stratum=stratum,N=n,Design_rank=rk,Residual_test_df=df,Unadjusted_r=r0,Unadjusted_P=p0,Condition_adjusted_r=r,Condition_adjusted_P=p)
 }
}
z<-do.call(rbind,rows);z$Unadjusted_BH_q<-ave(z$Unadjusted_P,z$Stratum,FUN=function(v)p.adjust(v,'BH'));z$Adjusted_BH_q<-ave(z$Condition_adjusted_P,z$Stratum,FUN=function(v)p.adjust(v,'BH'))
write.table(z,a[4],sep='\t',row.names=FALSE,quote=FALSE)
