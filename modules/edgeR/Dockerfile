FROM rocker/r2u
RUN R -e 'install.packages(c("dplyr", "BiocManager", "msigdbr", "readr")); BiocManager::install(c("clusterProfiler", "org.Mm.eg.db", "org.Hs.eg.db", "edgeR", "tximeta"))'
WORKDIR /working