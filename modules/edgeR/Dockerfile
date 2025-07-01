FROM rocker/r2u
RUN R -e 'install.packages(c("dplyr", "BiocManager", "msigdbr")); BiocManager::install(c("clusterProfiler", "org.Mm.eg.db", "org.Hs.eg.db", "edgeR", "tximeta"))'
WORKDIR /working