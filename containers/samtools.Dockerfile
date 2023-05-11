# Heavily modified from rocker-versioned Dockerfile.

FROM debian:stretch

LABEL maintainer="Jeffrey Chang <jeffrey.t.chang@uth.tmc.edu>"

ENV BUILD_DATE="$(TZ='America/Los_Angeles' date -I)" \
    LC_ALL=en_US.UTF-8 \
    LANG=en_US.UTF-8 \
    TERM=xterm


# Basic stuff.
# Don't know why, but apt-get update here prevents 404.
RUN apt-get update && apt-get install -y --no-install-recommends \
    apt-utils \
    bash-completion \
    ca-certificates \
    locales
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen && \
    /usr/sbin/locale-gen en_US.UTF-8 && \
    /usr/sbin/update-locale LANG=en_US.UTF-8


# Miscellaneous useful things.
# Don't know why, but apt-get update here prevents 404.
RUN apt-get update && apt-get install -y --no-install-recommends \
    file \
    less \
    vim

RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    curl


# Install samtools, bgzip, tabix.

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gcc \
    make \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev

RUN apt-get install -y --no-install-recommends \
    bzip2 \
    libncurses5-dev

ENV HTSLIB_VERSION=1.14
RUN cd /usr/local/src && \
    curl -L -O https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 && \
    bzcat htslib-${HTSLIB_VERSION}.tar.bz2 | tar xf - && \
    rm -f htslib-${HTSLIB_VERSION}.tar.bz2

RUN cd /usr/local/src/htslib-${HTSLIB_VERSION} && \
    ./configure && \
    make && \
    make install

ENV SAMTOOLS_VERSION=1.14
RUN cd /usr/local/src && \
    curl -L -O https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    bzcat samtools-${SAMTOOLS_VERSION}.tar.bz2 | tar xf - && \
    rm -f samtools-${SAMTOOLS_VERSION}.tar.bz2

RUN cd /usr/local/src/samtools-${SAMTOOLS_VERSION} && \
    ./configure && \
    make && \
    make install



  
# Clean up apt-get.
RUN cd / && \
    apt-get remove --purge -y $BUILDDEPS &&\
    apt-get autoremove -y && \
    apt-get autoclean -y && \
    rm -rf /var/lib/apt/lists/*

CMD ["samtools"]
