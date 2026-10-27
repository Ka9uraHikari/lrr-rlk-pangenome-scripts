#!/usr/bin/env Rscript
suppressMessages({
  library(DESeq2)
})

BASE <- '/path/to/project/18_transtome_raw'
QUANT <- file.path(BASE, 'quant/salmon')

crrs <- unique(read.table(file.path(BASE,'metadata/sample_table.tsv'), header=T, sep='\t', stringsAsFactors=F)$CRR)
counts <- NULL
samples_done <- c()
for (crr in sort(crrs)) {
  qf <- file.path(QUANT, crr, 'quant.sf')
  if (!file.exists(qf)) next
  d <- read.delim(qf)
  d$gene <- sub('\\..*$', '', d$Name)
  g <- aggregate(NumReads ~ gene, data=d, sum)
  if (is.null(counts)) {
    counts <- g
    colnames(counts)[2] <- crr
  } else {
    counts[[crr]] <- g$NumReads[match(counts$gene, g$gene)]
  }
  samples_done <- c(samples_done, crr)
}
counts[is.na(counts)] <- 0
rownames(counts) <- counts$gene
counts$gene <- NULL
cat("count矩阵:", nrow(counts), "基因 x", ncol(counts), "样本\n")

meta <- read.table(file.path(BASE,'metadata/sample_table.tsv'), header=T, sep='\t', stringsAsFactors=F)
meta <- meta[!duplicated(meta$CRR),]
rownames(meta) <- meta$CRR
meta <- meta[samples_done,]
meta$condition <- paste0(meta$pathogen, '_', meta$timepoint)
meta$pathogen <- factor(meta$pathogen)
meta$timepoint <- factor(meta$timepoint, levels=c('0d','1d','2d','3d','4d','5d','0h','6h','12h','24h','48h','72h'))

outdir <- file.path(BASE, 'results', 'DEG')
dir.create(outdir, recursive=T, showWarnings=F)
for (patho in c('PM','ALS','CC','GM')) {
  idx <- which(meta$pathogen == patho)
  sub_counts <- counts[, idx, drop=F]
  sub_meta <- meta[idx,]
  t0 <- ifelse(patho %in% c('PM','ALS'), '0d', '0h')
  sub_meta$timepoint <- relevel(factor(sub_meta$timepoint), ref=t0)
  keep <- rowSums(sub_counts >= 10) >= 3
  sub_counts <- sub_counts[keep,]
  dds <- DESeqDataSetFromMatrix(countData=round(sub_counts), colData=sub_meta, design=~timepoint)
  dds <- DESeq(dds, quiet=T)
  res_all <- NULL
  for (tp in setdiff(levels(sub_meta$timepoint), t0)) {
    res <- results(dds, contrast=c('timepoint', tp, t0))
    res <- as.data.frame(res)
    res$gene <- rownames(res)
    res$contrast <- paste0(patho, '_', tp, 'vs', t0)
    res_all <- rbind(res_all, res)
  }
  out <- file.path(outdir, paste0(patho, '_DEG_results.csv'))
  write.csv(res_all, out, row.names=F)
  sig <- res_all[!is.na(res_all$padj) & res_all$padj < 0.05 & abs(res_all$log2FoldChange) >= 1,]
  cat(patho, ": 显著DEG数 =", nrow(sig), "\n")
}
cat("DESeq2完成，结果在", outdir, "\n")
