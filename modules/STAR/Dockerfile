FROM ubuntu:24.04

RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    zlib1g-dev \
    libncurses5-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    xxd \ 
    grep \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /usr/local/src
RUN wget https://github.com/samtools/samtools/releases/download/1.22/samtools-1.22.tar.bz2 && \
    tar -xjf samtools-1.22.tar.bz2 && \
    cd samtools-1.22 && \
    make && \
    ./configure && \
    make && \
    make install && \
    mv samtools /usr/local/bin && \
    cd .. && \
    rm -rf samtools*
SHELL ["/bin/bash", "-lc"]
RUN wget -O STAR-2.7.11b.tar.gz https://github.com/alexdobin/STAR/archive/2.7.11b.tar.gz && \
    tar -xzf STAR-2.7.11b.tar.gz && \
    cd STAR-2.7.11b/source && \
    make STAR && \
    mv STAR /usr/local/bin && \
    cd /usr/local/src && \
    rm -rf STAR-2.7.11b.tar.gz STAR-2.7.11b
SHELL ["/bin/bash", "-c"]
WORKDIR /working
