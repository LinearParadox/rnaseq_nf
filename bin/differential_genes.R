#!/usr/bin/env Rscript

library(edgeR)
library(clusterProfiler)
library(dplyr)
library(limma)
library(msigdbr)
library(tximeta)

writeLines(capture.output(sessionInfo()), "sessionInfo_deg.txt")

#2 fasta
#3 gtf
#4 organism
#5 design_matrix
#6 contrast matrix
#7 min_gs size
#8 max_gs size



args = commandArgs(trailingOnly = T)
fasta = args[[1]]
gtf = args[[2]]
organism=args[[3]]
design = as.matrix(read.csv(args[[4]], header=T, row.names=1))
coef <- as.matrix(read.csv(args[[5]], header=T, row.names=1))
min_size = as.numeric(args[[6]])
max_size = as.numeric(args[[7]])
if (organism == "human") {
  library(org.Hs.eg.db)
  orgdb = org.Hs.eg.db
  makeLinkedTxome(indexDir="salmon_index", source="Gencode", organism="Homo sapiens",
                  release="--", genome="--", fasta=fasta, gtf=gtf, write=FALSE)
} else {
  library(org.Mm.eg.db)
  makeLinkedTxome(indexDir="salmon_index", source="Gencode", organism="Mus musculus",
                  release="--", genome="", fasta=fasta, gtf=gtf, write=FALSE)
  orgdb = org.Mm.eg.db
}
se <- tximeta(rownames(design), useHub = T)
se <- summarizeToGene(se, assignRanges="abundant")
gse <- addIds(se, "SYMBOL", gene=TRUE)
y <- makeDGEList(gse)
y$samples$lib.size <- colSums(y$counts)
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- normLibSizes(y)
dir.create("figs", showWarnings = FALSE)
png("figs/MDS_plot.png")
plotMDS(y)
dev.off()
print(design)
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
  write.csv(df, paste0("csv/de/", x, ".csv"))
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
  gseGO(geneList     = ranked_list,
        OrgDb        = orgdb,
        keyType      = "ENSEMBL",
        ont          = "CC",
        minGSSize = min_size,
        maxGSSize    = max_size,
        pvalueCutoff = 0.05,
        verbose      = FALSE)
  write.csv(gsea[], paste0("csv/gsea/", x, "/gocc.csv"))
  gseGO(geneList     = ranked_list,
        OrgDb        = orgdb,
        ont          = "MF",
        keyType      = "ENSEMBL",
        minGSSize = min_size,
        maxGSSize    = max_size,
        pvalueCutoff = 0.05,
        verbose      = FALSE)
  write.csv(gsea[], paste0("csv/gsea/", x, "/gomf.csv"))
  gseGO(geneList     = ranked_list,
        OrgDb        = orgdb,
        ont          = "BP",
        keyType      = "ENSEMBL",
        minGSSize = min_size,
        maxGSSize    = max_size,
        pvalueCutoff = 0.05,
        verbose      = FALSE)
  write.csv(gsea[], paste0("csv/gsea/", x, "/gomf.csv"))
})
