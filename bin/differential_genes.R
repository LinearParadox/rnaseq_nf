#!/usr/bin/env Rscript

library(edgeR)
library(clusterProfiler)
library(dplyr)
library(biomaRt)
library(tibble)
library(limma)
library(msigdbr)
library(tximport)
library(argparser, quietly=TRUE)
writeLines(capture.output(sessionInfo()), "sessionInfo_deg.txt")

get_samples <- function(contrast, design.mat){
  samps <- names(contrast)[contrast != 0]
  samp_names <- design[, samps]
  keep <- rowSums(samp_names != 0) > 0
  return(rownames(samp_names[keep, , drop=F]))
}

p <- arg_parser("Differential expression analysis with edgeR and GSEA")
p <- add_argument(p, "--organism", help="organism: human or mouse")
p <- add_argument(p, "--design", help="design matrix csv file")
p <- add_argument(p, "--coef", help="contrast matrix csv file")
p <- add_argument(p, "--min_size", help="minimum gene set size")
p <- add_argument(p, "--max_size", help="maximum gene set size")
#p <- add_argument(p, "--formula", help="design formula")
#p <- add_argument(p, "--ref", help="reference level for formula")
#p <- add_argument(p, "--meta", help="sample level metadata. Required if using a formula")

#TODO: ADD PARSING VIA META



args <- parse_args(p)
organism <- args$organism
design <- as.matrix(read.csv(args$design, header=T, row.names=1))
coef <- as.matrix(read.csv(args$coef, header=T, row.names=1))
min_size = as.numeric(args$min_size)
max_size = as.numeric(args$max_size)
if (organism == "human") {
  library(org.Hs.eg.db)
  orgdb = org.Hs.eg.db
  mart <- useEnsembl(dataset = "hsapiens_gene_ensembl", biomart='ensembl')
  symbol_key <- "hgnc_symbol"
  species <- "Homo sapiens"
  prefix <- "^ENSG"
} else {
  library(org.Mm.eg.db)
  orgdb = org.Mm.eg.db
  mart <- useEnsembl(dataset = "mmusculus_gene_ensembl", biomart='ensembl')
  symbol_key <- "mgi_symbol"
  species <- "Mus musculus"
  prefix <- "^ENSMUSG"

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
y <- DGEList(cts, genes=rownames(cts), samples=design)
y <- scaleOffset(y, normMat)
y$samples$lib.size <- colSums(y$counts)
keep <- filterByExpr(y, design = design)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- normLibSizes(y)
annot <-biomaRt::select(mart, keys=rownames(y$counts), keytype="ensembl_gene_id_version", columns=c("ensembl_gene_id_version", symbol_key))
annot <- annot[!duplicated(annot$ensembl_gene_id_version),]
transcripts <- rownames(y)[!(rownames(y) %in% annot$ensembl_gene_id_version)]
df <- data.frame(transcripts, rep("", length(transcripts)))
colnames(df) <- c("ensembl_gene_id_version", symbol_key)
annot <- rbind(annot, df)
colnames(annot) <- c("gene_id", "SYMBOL")
annot <- annot %>% mutate(SYMBOL=case_when(
  SYMBOL == "" ~ gsub("\\..*","", gene_id),
  .default = SYMBOL
))
rownames(annot) <- annot[,1]
y$genes$SYMBOL <- annot[rownames(y), "SYMBOL"]
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
counts <- cpm(y, offset=y$offset)
log_counts <- cpm(y, log=T, offset=y$offset)
write.csv(counts, paste0("csv/cpm.csv"))
write.csv(counts, paste0("csv/logcpm.csv"))
test<-lapply(colnames(coef), FUN=function(x){
  qlf <- glmQLFTest(fit, contrast=coef[,x])
  qlf <- qlf[grepl(prefix, rownames(qlf)),]
  df<-data.frame(topTags(qlf, n=Inf, adjust.method="BH", p.value=1)) %>%
    dplyr::select(c("SYMBOL", "gene_id", "logFC", "logCPM", "F", "PValue", "FDR"))
  sample.names <- get_samples(coef[,x], design)
  counts.keep <- counts[rownames(df), sample.names]
  df <- cbind(df, counts.keep[rownames(df),])
  df$gene_id <- gsub("\\..*","", df$gene_id)
  write.csv(df, paste0("csv/de/", x, ".csv"), row.names = F)
  ranked_list <- df$logFC * -log10(df$PValue)
  names(ranked_list) <- df$gene_id
  ranked_list <- ranked_list[!is.na(ranked_list)]
  ranked_list <- ranked_list[!is.infinite(ranked_list)]
  ranked_list <- ranked_list[order(ranked_list, decreasing = T)]
  hallmark <- msigdbr(species=species, collection="H") %>%
    dplyr::select(all_of(c("gs_name","ensembl_gene")))
  dir.create(paste0("csv/gsea/", x), showWarnings = FALSE)
  gsea <- clusterProfiler::GSEA(ranked_list, TERM2GENE = hallmark)
  write.csv(gsea[], paste0("csv/gsea/", x, "/hallmarks.csv"))
  biocarta <- msigdbr(species=species, subcollection="CP:BIOCARTA") %>%
    dplyr::select(all_of(c("gs_name","ensembl_gene")))
  gsea <- clusterProfiler::GSEA(ranked_list, TERM2GENE = biocarta)
  write.csv(gsea[], paste0("csv/gsea/", x, "/biocarta.csv"))
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
