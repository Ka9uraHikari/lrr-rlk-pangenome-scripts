#!/usr/bin/env Rscript
suppressMessages({library(WGCNA); allowWGCNAThreads()})
tpm <- read.csv('/path/to/project/18_transtome_raw/quant/gene_tpm_matrix.csv', row.names=1, check.names=F)
keep <- rowMeans(tpm) > 1
tpm <- tpm[keep,]
datExpr <- t(tpm)
datExpr <- datExpr[order(rownames(datExpr)),]
grp <- sub('_R[0-9]$','',rownames(datExpr))
vars <- apply(datExpr, 2, var)
datExpr <- datExpr[, order(vars, decreasing=T)[1:min(5000,ncol(datExpr))]]
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector=powers, verbose=0)
softPower <- sft$powerEstimate
if (is.na(softPower)) softPower <- 6
cat("软阈值:", softPower, "\n")
net <- blockwiseModules(datExpr, power=softPower, networkType='signed',
  TOMType='signed', minModuleSize=30, mergeCutHeight=0.25,
  numericLabels=T, saveTOMs=F, verbose=0)
mods <- net$colors
names(mods) <- colnames(datExpr)
tab <- table(mods)
cat("模块数:", length(tab), "\n")
print(head(tab, 10))
write.table(mods, '/path/to/project/18_transtome_raw/results/WGCNA_modules.csv', sep=',', col.names=F)
grp_f <- factor(grp)
mod_mean <- sapply(unique(mods), function(m) rowMeans(datExpr[, mods==m, drop=F]))
pvals <- apply(mod_mean, 2, function(x) anova(lm(x ~ grp_f))$'Pr(>F)'[1])
res <- data.frame(module=colnames(mod_mean), anova_p=pvals)
res$module <- paste0('M', res$module)
write.csv(res, '/path/to/project/18_transtome_raw/results/WGCNA_module_trait.csv', row.names=F)
cat("WGCNA完成，结果在 results/WGCNA_*.csv\n")
