ARG VARIANT="ubuntu-22.04"
FROM mcr.microsoft.com/vscode/devcontainers/base:0-${VARIANT}

# Install core dev tools, Java for Nextflow, Docker CLI
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       build-essential \
       curl \
       git \
       zip \
       unzip \
       ca-certificates \
       python3 python3-pip

RUN apt-get update \
	&& apt-get install -y ca-certificates curl \
	&& install -m 0755 -d /etc/apt/keyrings \
	&& curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc \
	&& chmod a+r /etc/apt/keyrings/docker.asc \
	&& echo \
		"deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] https://download.docker.com/linux/ubuntu \
		$(. /etc/os-release && echo "${UBUNTU_CODENAME:-$VERSION_CODENAME}") stable" | \
		tee /etc/apt/sources.list.d/docker.list > /dev/null \
	&& apt-get update \
	&& apt-get install -y docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin \   
	&& rm -rf /var/lib/apt/lists/*



RUN mkdir -p /home/vscode/.nextflow \
    && chown -R vscode:vscode /home/vscode

USER vscode
RUN sudo usermod -aG docker vscode
ENV HOME=/home/vscode 
ENV NXF_HOME=/home/vscode/.nextflow
WORKDIR /home/vscode

#Install SDKMAN and Java
SHELL ["/bin/bash", "-lc"]

# Installing SDK+Java+Nextflow
RUN curl -s "https://get.sdkman.io?ci=true" | bash \
 && source "$HOME/.sdkman/bin/sdkman-init.sh" \
 && sdk install java 17.0.10-tem \
 && sdk use java 17.0.10-tem \
 && curl -fsSL https://get.nextflow.io | bash

RUN sudo mv nextflow /usr/local/bin && sudo chmod +x /usr/local/bin/nextflow

SHELL ["/bin/sh", "-c"]


# Install nf-core tools for pipeline linting
RUN pip3 install --no-cache-dir nf-core

WORKDIR /workspace
