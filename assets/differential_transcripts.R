library(edgeR)
library(clusterProfiler)
library(dplyr)
library(limma)
library(msigdbr)
library(tximeta)

#1 samples
#2 conditions
#3 fasta
#4 gtf
#5 organism
#6 contrasts
#7 design formula
#8 reference levels (optional)


args = commandArgs(trailingOnly = T)
samples = read.csv(args[[1]], header=F)
cond = read.csv(args[[2]], header=T)
#samples = read.csv("gene_out_paths.csv", header=F)
contrasts <- read.csv(args[[6]], header=F)
contrasts[,2] <- gsub(":", "..", contrasts[,2])
colnames(samples) <- c("names", "files")
#cond = read.csv("cond_out.csv", header=T)
colnames(cond)[1] <- "names"
merged_df <- dplyr::left_join(cond, samples, by="names")
#param3 fasta
#param4 gtf
#params5 == organism
if (args[[5]] == "human") {
  library(org.Hs.eg.db)
  orgdb = org.Hs.eg.db
  makeLinkedTxome(indexDir="quant_files/salmon_index", source="Gencode", organism="Homo sapiens",
                  release="--", genome="", fasta=args[[3]], gtf=args[[4]], write=FALSE)
  organism="human"
} else {
  library(org.Mm.eg.db)
  makeLinkedTxome(indexDir="quant_files/salmon_index", source="Gencode", organism="Mus musculus",
                  release="--", genome="", fasta=args[[3]], gtf=args[[4]], write=FALSE)
  organism="mouse"
  orgdb = org.Mm.eg.db
}
se <- tximeta(merged_df)
se <- summarizeToGene(se, assignRanges="abundant")
gse <- addIds(se, "SYMBOL", gene=TRUE)
y <- makeDGEList(gse)
y$samples$lib.size <- colSums(y$counts)
y <- normLibSizes(y)
keep<-filterByExpr(y)
y <-y[keep,,keep.lib.sizes=FALSE]
dir.create("figs", showWarnings = FALSE)
png("figs/MDS_plot.png")
plotMDS(y)
dev.off()
for(x in colnames(cond[,-c(1)])){
  y$samples[,x] <- factor(y$samples[,x])
}
#arg[-1] references
if(length(args) == 8){
  relevel_df <- read.csv(args[[8]], row.names=1, header=F)
  print(relevel_df)
  #relevel_df <- read.csv("relevel.csv", row.names=1, header=F)
  for(x in rownames(relevel_df)){
    y$samples[,x] <- relevel(y$samples[,x], ref=relevel_df[x,1])
  }
}
design <- model.matrix(as.formula(args[[7]]), data=y$samples)
print(contrasts)
print(contrasts[,2])
my.contrasts <- makeContrasts(contrasts=contrasts[,2], levels=design)
colnames(my.contrasts) <- contrasts[,1]
fit <- glmQLFit(y,design, robust=T)
png("figs/edgeR_counts_fit.png")
plotQLDisp(fit)
dev.off()
dir.create("csv", showWarnings = FALSE)
dir.create("csv/gsea", showWarnings = FALSE)
dir.create("csv/de", showWarnings = FALSE)
test<-lapply(colnames(my.contrasts), FUN=function(x){
  qlf <- glmQLFTest(fit, contrast=my.contrasts[,x])
  df<-data.frame(topTags(qlf, n=Inf, adjust.method="BH", p.value=1)) %>%
    dplyr::select(c("SYMBOL", "gene_id", "logFC", "logCPM", "F", "PValue", "FDR"))
  df$gene_id <- gsub("\\..*","", df$gene_id)
  write.csv(df, paste0("csv/de/", x, ".csv"))
  ranked_list <- df$logFC * -log10(df$PValue)
  names(ranked_list) <- df$gene_id
  ranked_list <- ranked_list[order(ranked_list, decreasing = T)]
  hallmark <- msigdbr(species=organism, category="H") %>%
    dplyr::select(c(gs_name, ensembl_gene))
  gsea <- clusterProfiler::GSEA(ranked_list, TERM2GENE = hallmark)
  write.csv(gsea[], paste0("csv/gsea/", x, "_hallmarks.csv"))
  gseGO(geneList     = ranked_list,
        OrgDb        = orgdb,
        keyType      = "ENSEMBL",
        ont          = "CC",
        maxGSSize    = 500,
        pvalueCutoff = 0.05,
        verbose      = FALSE)
  write.csv(gsea[], paste0("csv/gsea/", x, "_gocc.csv"))
  gseGO(geneList     = ranked_list,
        OrgDb        = orgdb,
        ont          = "MF",
        keyType      = "ENSEMBL",
        maxGSSize    = 500,
        pvalueCutoff = 0.05,
        verbose      = FALSE)
  write.csv(gsea[], paste0("csv/gsea/", x, "_gomf.csv"))
  gseGO(geneList     = ranked_list,
        OrgDb        = orgdb,
        ont          = "BP",
        keyType      = "ENSEMBL",
        maxGSSize    = 500,
        pvalueCutoff = 0.05,
        verbose      = FALSE)
  write.csv(gsea[], paste0("csv/gsea/", x, "_gomf.csv"))
})
