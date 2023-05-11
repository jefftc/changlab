# Heavily modified from rocker-versioned Dockerfile.

# /usr/local/bin/Platypus.py
# /usr/local/src/Platypus_0.8.1/


FROM debian:stretch

LABEL maintainer="Jeffrey Chang <jeffrey.t.chang@uth.tmc.edu>"

ENV BUILD_DATE="$(TZ='America/Los_Angeles' date -I)" \
    LC_ALL=en_US.UTF-8 \
    LANG=en_US.UTF-8 \
    TERM=xterm


# Basic stuff.
RUN apt-get update
# Don't know why, but apt-get update here prevents 404.
RUN apt-get update && apt-get install -y --no-install-recommends \
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
    vim \
    openssl \
    libssl-dev \
    libcurl4-openssl-dev \
    curl


# Developer stuff.
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gcc \
    make \
    git

# Python
RUN apt-get update && apt-get install -y --no-install-recommends \
    python2.7 \
    python2.7-dev \
    python-pip



# HTSLIB

RUN apt-get update && apt-get install -y --no-install-recommends \
    libbz2-dev \
    zlib1g-dev

RUN apt-get update && apt-get install -y --no-install-recommends \
    liblzma-dev

ENV HTSLIB_VERSION=1.10.2
RUN cd /usr/local/src && \
    curl -L -O https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 && \
    bzcat htslib-${HTSLIB_VERSION}.tar.bz2 | tar xf - && \
    rm -f htslib-${HTSLIB_VERSION}.tar.bz2

# Install in /usr/lib/ so Platypus can find it.
RUN cd /usr/local/src/htslib-${HTSLIB_VERSION} && \
    ./configure --prefix=/usr && \
    make && \
    make install



# Platypus
# Python 2.6+
# htslib 1.2.1+
# https://www.well.ox.ac.uk/research/research-groups/lunter-group/lunter-group/platypus-a-haplotype-based-variant-caller-for-next-generation-sequence-data

ENV PLATYPUS_VERSION=0.8.1
RUN cd /usr/local/src && \
    curl -L https://www.rdm.ox.ac.uk/resolveuid/599a7efc8ec04059a101c59714353209 > Platypus_${PLATYPUS_VERSION}.tgz && \
    tar xvfz Platypus_${PLATYPUS_VERSION}.tgz && \
    rm -f Platypus_${PLATYPUS_VERSION}.tgz

RUN cd /usr/local/src/Platypus_${PLATYPUS_VERSION} && \
    ./buildPlatypus.sh

RUN ln -s /usr/local/src/Platypus_${PLATYPUS_VERSION}/Platypus.py \
  /usr/local/bin


CMD ["Platypus.py"]
