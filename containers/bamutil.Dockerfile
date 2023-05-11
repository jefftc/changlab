# Heavily modified from rocker-versioned Dockerfile.

# /usr/local/bin/
#   bam

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


# Developer stuff.
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gcc \
    make \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    libncurses5-dev



# libStatGen

RUN cd /usr/local/src && \
    curl -L -O https://github.com/statgen/libStatGen/archive/v1.0.14.tar.gz &&\
    tar xvfz v1.0.14.tar.gz
RUN cd /usr/local/src/libStatGen-1.0.14 && \
    make all
RUN ln -s /usr/local/src/libStatGen-1.0.14 /usr/local/src/libStatGen


RUN cd /usr/local/src && \
    curl -L -O https://github.com/statgen/bamUtil/archive/v1.0.14.tar.gz && \
    tar xvfz v1.0.14.tar.gz
RUN cd /usr/local/src/bamUtil-1.0.14 && \
    make && \
    make install






# Clean up apt-get.
RUN cd / && \
    apt-get remove --purge -y $BUILDDEPS &&\
    apt-get autoremove -y && \
    apt-get autoclean -y && \
    rm -rf /var/lib/apt/lists/*

CMD ["bam"]
