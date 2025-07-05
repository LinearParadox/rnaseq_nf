#!/usr/bin/env Rscript

library(edgeR)
library(clusterProfiler)
library(dplyr)
library(biomaRt)
library(limma)
library(msigdbr)
library(tximport)
writeLines(capture.output(sessionInfo()), "sessionInfo_deg.txt")

#4 organism
#5 design_matrix
#6 contrast matrix
#7 min_gs size
#8 max_gs size



args = commandArgs(trailingOnly = T)
organism=args[[1]]
design = as.matrix(read.csv(args[[2]], header=T, row.names=1))
coef <- as.matrix(read.csv(args[[3]], header=T, row.names=1))
min_size = as.numeric(args[[4]])
max_size = as.numeric(args[[5]])
if (organism == "human") {
  library(org.Hs.eg.db)
  orgdb = org.Hs.eg.db
  mart <- useEnsembl(dataset = "hsapiens_gene_ensembl", biomart='ensembl')
  symbol_key <- "hgnc_symbol"
} else {
  library(org.Mm.eg.db)
  orgdb = org.Mm.eg.db
  mart <- useEnsembl(dataset = "mmusculus_gene_ensembl", biomart='ensembl')
  symbol_key <- "mgi_symbol"
}
files <- paste0(rownames(design), "/", 'quant.genes.sf')
names(files) <- rownames(design)
txi <- tximport(files, type="salmon", txIn=F, geneIdCol="Name", lengthCol="EffectiveLength", abundanceCol="TPM",countsCol="NumReads")
cts <- txi$counts
normMat <- txi$length
normMat <- normMat/exp(rowMeans(log(normMat)))
normCts <- cts/normMat
eff.lib <- calcNormFactors(normCts) * colSums(normCts)
normMat <- sweep(normMat, 2, eff.lib, "*")
normMat <- log(normMat)
y <- DGEList(cts, genes=rownames(cts))
y <- scaleOffset(y, normMat)

y$samples$lib.size <- colSums(y$counts)
keep <- filterByExpr(y, design = design)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- normLibSizes(y)
annot <-biomaRt::select(mart, keys=rownames(y$counts), keytype="ensembl_gene_id_version", columns=c("ensembl_gene_id_version", symbol_key))
annot <- annot[!duplicated(annot$ensembl_gene_id_version),]
rownames(annot) <- annot[,1]
y$genes$SYMBOL <- annot[[symbol_key]]
colnames(y$genes) <- c("gene_id", "SYMBOL")
dir.create("figs", showWarnings = FALSE)
png("figs/MDS_plot.png")
plotMDS(y)
dev.off()
fit <- glmQLFit(y, design=design, robust=T)
png("figs/edgeR_counts_fit.png")
plotQLDisp(fit)
dev.off()
dir.create("csv", showWarnings = FALSE)
dir.create("csv/gsea", showWarnings = FALSE)
dir.create("csv/de", showWarnings = FALSE)
test<-lapply(colnames(coef), FUN=function(x){
  qlf <- glmQLFTest(fit, contrast=coef[,x])
  df<-data.frame(topTags(qlf, n=Inf, adjust.method="BH", p.value=1)) %>%
    dplyr::select(c("SYMBOL", "gene_id", "logFC", "logCPM", "F", "PValue", "FDR"))
  df$gene_id <- gsub("\\..*","", df$gene_id)
  write.csv(df, paste0("csv/de/", x, ".csv"), row.names = F)
  ranked_list <- df$logFC * -log10(df$PValue)
  names(ranked_list) <- df$gene_id
  ranked_list <- ranked_list[!is.na(ranked_list)]
  ranked_list <- ranked_list[!is.infinite(ranked_list)]
  ranked_list <- ranked_list[order(ranked_list, decreasing = T)]
  hallmark <- msigdbr(species=organism, collection="H") %>%
    dplyr::select(c(gs_name, ensembl_gene))
  dir.create(paste0("csv/gsea/", x), showWarnings = FALSE)
  gsea <- clusterProfiler::GSEA(ranked_list, TERM2GENE = hallmark)
  write.csv(gsea[], paste0("csv/gsea/", x, "/hallmarks.csv"))
  gsea<-gseGO(geneList     = ranked_list,

        OrgDb        = orgdb,
        keyType      = "ENSEMBL",
        ont          = "CC",
        minGSSize = min_size,
        maxGSSize    = max_size,
        pvalueCutoff = 0.05,
        verbose      = FALSE)
  write.csv(gsea[], paste0("csv/gsea/", x, "/gocc.csv"))

  gsea<-gseGO(geneList     = ranked_list,
        OrgDb        = orgdb,
        ont          = "MF",
        keyType      = "ENSEMBL",
        minGSSize = min_size,
        maxGSSize    = max_size,
        pvalueCutoff = 0.05,
        verbose      = FALSE)
  write.csv(gsea[], paste0("csv/gsea/", x, "/gomf.csv"))
  gsea<-gseGO(geneList     = ranked_list,

        OrgDb        = orgdb,
        ont          = "BP",
        keyType      = "ENSEMBL",
        minGSSize = min_size,
        maxGSSize    = max_size,
        pvalueCutoff = 0.05,
        verbose      = FALSE)
  write.csv(gsea[], paste0("csv/gsea/", x, "/gobp.csv"))
})
