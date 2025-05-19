FROM ubuntu:latest AS base
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends  libcurl4-openssl-dev\
    build-essential wget curl ca-certificates bzip2 gnupg git libssl-dev libffi-dev libz-dev libbz2-dev && rm -rf /var/lib/apt/lists/* && \
    mkdir -p /ref/ /script/ /raw_data/ /outdir/ && \
    wget -q -O /miniforge.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh && \
    bash /miniforge.sh -b -p /opt/conda && rm /miniforge.sh

FROM base AS pangolin
RUN git clone https://github.com/cov-lineages/pangolin.git && cd pangolin/ &&  \
    /opt/conda/bin/mamba env create -f environment.yml --name pangolin &&  \
    /opt/conda/envs/pangolin/bin/pip install . &&  rm -rf /opt/pangolin/

FROM pangolin AS nextclade
RUN cd /bin/ && wget -O nextclade https://github.com/nextstrain/nextclade/releases/latest/download/nextclade-x86_64-unknown-linux-musl && chmod u+x ./nextclade && \
    nextclade dataset get --name='sars-cov-2' --output-dir='/root/nextclade/'

FROM nextclade AS Freyja
RUN /opt/conda/bin/mamba create --name Freyja --channel conda-forge --channel bioconda --channel defaults freyja==1.5.3 && \
    cd /opt/conda/envs/Freyja && git clone https://github.com/andersen-lab/Freyja-data.git


