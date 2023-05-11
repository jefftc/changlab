# Heavily modified from rocker-versioned Dockerfile.

# /usr/local/bin/
#   bowtie2
#   STAR
#   rsem-*

FROM debian:buster

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


# Developer tools needed to build RSEM.
RUN apt-get update && apt-get install -y --no-install-recommends \
    make \
    g++ \
    zlib1g-dev






# Install programs to uncompress FASTQ files.

RUN apt-get update && apt-get install -y --no-install-recommends \
    bash-completion

RUN apt-get update && apt-get install -y --no-install-recommends \
    bzip2 \
    xz-utils

RUN ln -s /bin/zcat /usr/bin/gzcat



# Install bowtie.

RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    unzip

ENV A=bowtie-1.3.0-linux-x86_64
RUN cd /usr/local/src && \
  wget --default-page=${A}.zip \
    http://sourceforge.net/projects/bowtie-bio/files/bowtie/1.3.0/${A}.zip/ &&\
  unzip ${A}.zip && \
  rm ${A}.zip

RUN cp -p /usr/local/src/${A}/bowtie* /usr/local/bin



# Install bowtie2.

ENV A=bowtie2-2.4.2-linux-x86_64
RUN cd /usr/local/src && \
    wget --default-page=${A}.zip \
    http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.2/${A}.zip/&&\
    unzip ${A}.zip && \
    rm ${A}.zip

RUN cp -p /usr/local/src/${A}/bowtie* /usr/local/bin/


# Install STAR.

RUN apt-get update && apt-get install -y --no-install-recommends \
    git

RUN cd /usr/local/src && \
    git clone https://github.com/alexdobin/STAR.git && \
    cd STAR/source && \
    make STAR && \
    cp -Rp /usr/local/src/STAR/bin/Linux_x86_64_static/STAR* /usr/local/bin/&&\
    cd && \
    rm -rf /usr/local/src/STAR



# Install RSEM.

# Need Python for --gff3 option of rsem-prepare-reference.
# System PERL is broken.  Get error about missing POD/Usage.pm.
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    perl

RUN ln -s /usr/bin/python3 /usr/bin/python

RUN cd /usr/local/src && \
    curl -L -O https://github.com/deweylab/RSEM/archive/v1.3.3.tar.gz && \
    tar xvfz v1.3.3.tar.gz && \
    rm -f v1.3.3.tar.gz

RUN cd /usr/local/src/RSEM-1.3.3 && \
    make install



# Use bash so that unnamed pipes are evaluated correctly.  (sh does
# not support unnamed pipes).
SHELL ["/bin/bash", "-c"]
CMD ["/bin/bash"]



