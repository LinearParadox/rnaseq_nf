#!/usr/bin/env Rscript

library(edgeR)
library(dplyr)
library(biomaRt)

writeLines(capture.output(sessionInfo()), "sessionInfo_det.txt")
args = commandArgs(trailingOnly = T)
org <- args[[1]]
design <- as.matrix(read.csv(args[[2]], header=T, row.names=1))
coef <- as.matrix(read.csv(args[[3]], header=T, row.names=1))

if (org == "human") {
  library(org.Hs.eg.db)
  orgdb = org.Hs.eg.db
  mart <- useEnsembl(dataset = "hsapiens_gene_ensembl", biomart='ensembl')
} else {
  library(org.Mm.eg.db)
  orgdb = org.Mm.eg.db
  mart <- useEnsembl(dataset = "mmusculus_gene_ensembl", biomart='ensembl')
}
catch <- catchSalmon(rownames(design))
divided.counts <- catch$counts/catch$annotation$Overdispersion
colnames(divided.counts)<-rownames(design)[match(colnames(divided.counts), rownames(design))]
y <- DGEList(counts = divided.counts,
             genes = catch$annotation)
y$samples$lib.size <- colSums(y$counts)
keep<-filterByExpr(y, design = design)
y <-y[keep,,keep.lib.sizes=FALSE]
y <- normLibSizes(y)
annot <-biomaRt::select(mart, keys=rownames(y$genes), keytype="ensembl_transcript_id_version", columns=c("ensembl_transcript_id_version", "hgnc_symbol","transcript_biotype", "transcript_is_canonical"))
annot<-annot[!duplicated(annot$ensembl_transcript_id_version),]
annot[!is.na(annot$transcript_is_canonical), "transcript_is_canonical"] <- "TRUE"
annot[is.na(annot$transcript_is_canonical), "transcript_is_canonical"] <- "FALSE"
rownames(annot) <- annot[,1]
y$genes$SYMBOL <- annot$hgnc_symbol
y$genes$TYPE <- annot$transcript_biotype
y$genes$is_canonical <- annot$transcript_is_canonical
dir.create("figs", showWarnings = FALSE)
png("figs/MDS_plot.png")
plotMDS(y)
dev.off()
fit <- glmQLFit(y,design, robust=T)
png("figs/edgeR_counts_fit.png")
plotQLDisp(fit)
dev.off()
dir.create("csv", showWarnings = FALSE)
dir.create("csv/de", showWarnings = FALSE)
test<-lapply(colnames(coef), FUN=function(x){
  qlf <- glmQLFTest(fit, contrast=coef[,x])
  df<-data.frame(topTags(qlf, n=Inf, adjust.method="BH", p.value=1)) %>%
    dplyr::select(c("SYMBOL", "TYPE", "is_canonical","logFC", "logCPM", "F", "PValue", "FDR"))
  write.csv(df, paste0("csv/de/", x, ".csv"))
})
