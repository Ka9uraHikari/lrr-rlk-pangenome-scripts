#!/usr/bin/env Rscript
suppressMessages({})
tpm <- read.csv('quant/gene_tpm_matrix.csv', row.names=1, check.names=F)
cores <- c('CsaV4_7G002666','CsaV4_5G000526','CsaV4_2G002661','CsaV4_4G002866','CsaV4_7G002654','CsaV4_4G002864')
groups <- list(PM=c('0d','1d','3d','5d'), ALS=c('0d','1d','2d','4d'), CC=c('0h','6h','24h','72h'), GM=c('0h','6h','12h','48h'))
repcols <- function(p) colnames(tpm)[startsWith(colnames(tpm), paste0(p,'_'))]
cat("基因 | 处理 | 时间点表达(3重复均值TPM) | 动态模式\n")
for (g in cores) {
  if (!g %in% rownames(tpm)) next
  for (p in names(groups)) {
    tps <- groups[[p]]
    vals <- sapply(tps, function(tp) mean(as.numeric(tpm[g, repcols(paste0(p,'_',tp))])))
    peak <- which.max(vals)
    if (vals[1] > max(vals[-1])) pat <- '持续下降'
    else if (peak == length(vals)) pat <- '持续诱导'
    else if (peak == 2) pat <- '早期响应'
    else pat <- '中期响应'
    cat(sprintf("%s | %s | %s | %s\n", g, p, paste(sprintf("%.1f",vals),collapse=','), pat))
  }
}
