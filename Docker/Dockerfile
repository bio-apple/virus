FROM ubuntu:latest AS base
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends  libcurl4-openssl-dev  \
    build-essential wget curl ca-certificates bzip2 gnupg git libssl-dev libffi-dev libz-dev libbz2-dev && rm -rf /var/lib/apt/lists/* && \
    mkdir -p /ref/ /script/ /raw_data/ /outdir/ && \
    wget -q -O /miniforge.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh && \
    bash /miniforge.sh -b -p /opt/conda && rm /miniforge.sh

FROM base AS pangolin
RUN git clone https://github.com/cov-lineages/pangolin.git && cd pangolin/ &&  \
    /opt/conda/bin/mamba env create -f environment.yml --name pangolin &&  \
    /opt/conda/envs/pangolin/bin/pip install . &&  rm -rf /opt/pangolin/

FROM pangolin AS nextclade
RUN cd /bin/ && wget -O nextclade https://github.com/nextstrain/nextclade/releases/latest/download/nextclade-x86_64-unknown-linux-musl && chmod u+x ./nextclade

FROM nextclade AS freyja
RUN /opt/conda/bin/mamba create --name Freyja --channel conda-forge --channel bioconda --channel defaults freyja==1.5.3 && \
    cd /opt/conda/envs/Freyja && git clone https://github.com/andersen-lab/Freyja-data.git

FROM freyja AS virus
RUN /opt/conda/bin/mamba install --channel conda-forge --channel bioconda --channel defaults ivar=1.3 trimmomatic  \
    bowtie2 bbmap fastp seqtk samtools bedtools bcftools bwa prinseq cutadapt drep jvarkit entrez-direct ncbi-datasets-cli &&  \
    /opt/conda/bin/pip3 install seaborn matplotlib numpy pysam pandas && /opt/conda/bin/conda clean -a -y
