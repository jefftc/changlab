# Heavily modified from rocker-versioned Dockerfile.


# https://icbi-lab.github.io/immunedeconv/articles/immunedeconv.html

# /usr/local/CIBERSORT/
#   CIBERSORT.R
#   LM22.txt
# /usr/local/changlab/Rlib/


FROM debian:stretch

LABEL maintainer="Jeffrey Chang <jeffrey.t.chang@uth.tmc.edu>"

ENV BUILD_DATE="$(TZ='America/Los_Angeles' date -I)" \
    TERM=xterm
# Silences warnings with installing locales.
# debconf: unable to initialize frontend: Dialog
ENV DEBIAN_FRONTEND=noninteractive


# Basic stuff.
RUN apt-get update
RUN apt-get install -y --no-install-recommends \
    apt-utils 
RUN apt-get install -y --no-install-recommends \
    locales
ENV LC_ALL=en_US.UTF-8 \
    LANG=en_US.UTF-8
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen && \
    /usr/sbin/locale-gen en_US.UTF-8 && \
    /usr/sbin/update-locale LANG=en_US.UTF-8

# locales needs to be installed before this can be installed.
RUN apt-get install -y --no-install-recommends \
    ca-certificates

# Miscellaneous useful things.
# Don't know why, but apt-get update here prevents 404.
RUN apt-get update && apt-get install -y --no-install-recommends \
    bash-completion \
    file \
    less \
    openssl \
    libssl-dev \
    curl \
    libcurl4-openssl-dev


# Install miniconda.

RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    bzip2 \
    libglib2.0-0 \
    libxext6 \
    libsm6 \
    libxrender1

RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    mercurial \
    subversion

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda2-4.7.12-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc && \
    find /opt/conda/ -follow -type f -name '*.a' -delete && \
    find /opt/conda/ -follow -type f -name '*.js.map' -delete && \
    /opt/conda/bin/conda clean -afy

ENV PATH /opt/conda/bin:$PATH


# Install immunedconv.

RUN conda install -c bioconda -c conda-forge r-immunedeconv


# Install files for CIBERSORT.

RUN mkdir -p /usr/local/CIBERSORT/
COPY CIBERSORT.R LM22.txt /usr/local/CIBERSORT/


# Install Rlib for convenience functions.

RUN mkdir -p /usr/local/changlab/
COPY Rlib /usr/local/changlab/Rlib




CMD ["R"]
