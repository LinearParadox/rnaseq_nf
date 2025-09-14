#!/usr/bin/env Rscript

library(edgeR)
library(dplyr)
library(biomaRt)
library(tibble)
library(argparser)

get_samples <- function(contrast, design.mat){
  samps <- names(contrast)[contrast != 0]
  samp_names <- design[, samps]
  keep <- rowSums(samp_names != 0) > 0
  return(rownames(samp_names[keep, , drop=F]))
}


writeLines(capture.output(sessionInfo()), "sessionInfo_det.txt")

p <- arg_parser("Differential expression analysis with edgeR and GSEA")
p <- add_argument(p, "--organism", help="organism: human or mouse")
p <- add_argument(p, "--design", help="design matrix csv file")
p <- add_argument(p, "--coef", help="contrast matrix csv file")

args <- parse_args(p)



org <- args$organism
design <- as.matrix(read.csv(args$design, header=T, row.names=1))
coef <- as.matrix(read.csv(args$coef, header=T, row.names=1))

if (org == "human") {
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
files <- rownames(design)
catch <- catchSalmon(files)
divided.counts <- catch$counts/catch$annotation$Overdispersion
y <- DGEList(counts = divided.counts,
             genes = catch$annotation)
y$samples$lib.size <- colSums(y$counts)
keep<-filterByExpr(y, design = design)
y <-y[keep,,keep.lib.sizes=FALSE]
y <- normLibSizes(y)

annot <-biomaRt::select(mart, keys=rownames(y$genes), keytype="ensembl_transcript_id_version", columns=c("ensembl_transcript_id_version", symbol_key,"transcript_biotype", "transcript_is_canonical"))
annot<-annot[!duplicated(annot$ensembl_transcript_id_version),]
annot[!is.na(annot$transcript_is_canonical), "transcript_is_canonical"] <- "TRUE"
annot[is.na(annot$transcript_is_canonical), "transcript_is_canonical"] <- "FALSE"
rownames(annot) <- annot[,1]
annot <- annot[rownames(y),]
y$genes$SYMBOL <- annot[[symbol_key]]
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
counts <- cpm(y, offset=y$offset)
log_counts <- cpm(y, log=T, offset=y$offset)
write.csv(counts, paste0("csv/cpm.csv"))
write.csv(counts, paste0("csv/logcpm.csv"))
test<-lapply(colnames(coef), FUN=function(x){
  qlf <- glmQLFTest(fit, contrast=coef[,x])
  df<-data.frame(topTags(qlf, n=Inf, adjust.method="BH", p.value=1)) %>%
    dplyr::select(c("SYMBOL", "TYPE", "is_canonical","logFC", "logCPM", "F", "PValue", "FDR"))
  sample.names <- get_samples(coef[,x], design)
  counts.keep <- counts[rownames(df), sample.names]
  df <- cbind(df, counts.keep[rownames(df),])
  write.csv(df, paste0("csv/de/", x, ".csv"))
})
