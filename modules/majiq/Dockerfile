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
    python3-venv \
    python3-pip \
    git && \
    rm -rf /var/lib/apt/lists/*
WORKDIR /usr/local/src
RUN wget https://github.com/samtools/htslib/releases/download/1.22/htslib-1.22.tar.bz2 && \
    tar -xjf htslib-1.22.tar.bz2 && \
    cd htslib-1.22 && \
    make && \
    ./configure && \
    make && \
    make install && \
    mv htslib /usr/local/bin && \
    cd .. && \
    rm -rf htslib*
RUN python3 -m venv /opt/venv

# Add venv to path
ENV PATH="/opt/venv/bin:$PATH"
RUN apt-get update && apt-get install -y libzstd-dev
RUN pip install --upgrade pip setuptools
RUN pip install --no-cache-dir git+https://bitbucket.org/biociphers/moccasin@new_moccasin git+https://bitbucket.org/biociphers/majiq_academic.git
SHELL ["/bin/bash", "-c"]
WORKDIR /working