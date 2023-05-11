# Heavily modified from rocker-versioned Dockerfile.

# /usr/local/strelka_workflow/


FROM ubuntu:16.04

# Somehow, Strelka doesn't build on debian.  Missing some C++ STL
# stuff.
#FROM debian:stretch

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



# Install Strelka.

RUN apt-get update && apt-get install -y --no-install-recommends \
  zlib1g-dev \
  rsync \
  python2.7
RUN ln -s /usr/bin/python2.7 /usr/bin/python


RUN cd /usr/local/src && \
  curl -u strelka: -O ftp://ftp.illumina.com/v1-branch/v1.0.15/strelka_workflow-1.0.15.tar.gz

RUN cd /usr/local/src && \
  tar xvfz strelka_workflow-1.0.15.tar.gz && \
  cd strelka_workflow-1.0.15 && \
  ./configure --prefix=/usr/local/strelka_workflow && \
  make
