#FROM ubuntu:20.04  # doesn't work.  compile errors in joinx
#FROM ubuntu:16.04  # doesn't work.  compile errors in joinx
FROM ubuntu:14.04

LABEL maintainer="Jeffrey Chang <jeffrey.t.chang@uth.tmc.edu>"

ENV BUILD_DATE="$(TZ='America/Los_Angeles' date -I)" \
    LC_ALL=C \
    TERM=xterm

# Basic stuff.
RUN apt-get update
RUN apt-get update && apt-get install -y --no-install-recommends \
    apt-utils
RUN DEBIAN_FRONTEND=noninteractive TZ="America/New_York" \
    apt-get -y install tzdata
RUN apt-get install -y --no-install-recommends \
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
    openssl \
    libssl-dev \
    curl \
    libcurl4-openssl-dev



# Install prerequisites.

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    git \
    cmake \
    curl \
    cpanminus \
    libbz2-dev \
    libgtest-dev \
    libbam-dev \
    zlib1g-dev 



# Install samtools.

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


RUN cd /usr/local/src && \
    curl -k -L https://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2/download > samtools-0.1.19.tar.bz2

RUN cd /usr/local/src && \
    tar jxf samtools-0.1.19.tar.bz2 && \
    cd samtools-0.1.19 && \
    make && \
    mv samtools /usr/local/bin/



# Install calcRoiCovg.

RUN cd /usr/local/src && \
    git clone https://github.com/Beifang/calcRoiCovg.git
ENV SAMTOOLS_ROOT=/usr/local/src/samtools-0.1.19
RUN cd /usr/local/src/calcRoiCovg && \
    make && \
    mv calcRoiCovg /usr/local/bin/


# Install bedtools.

RUN apt-get update && apt-get install -y --no-install-recommends \
    python2.7 && \
    ln -s /usr/bin/python2.7 /usr/bin/python

RUN cd /usr/local/src && \
    curl -L -O https://github.com/arq5x/bedtools2/archive/v2.27.1.tar.gz

RUN cd /usr/local/src && \
    tar -zxvf v2.27.1.tar.gz && \
    cd bedtools2-2.27.1 && \
    make && \
    mv /usr/local/src/bedtools2-2.27.1/bin/* /usr/local/bin/


# Install joinx.

RUN cd /usr/local/src && \
    git clone --recursive https://github.com/genome/joinx.git

# Fix bug.
RUN cd /usr/local/src/joinx/src/lib/io && \
    perl -i -p -e 's/return std::getline\(_in, line\);/std::getline\(_in, line\); return true;/' StreamLineSource.cpp

RUN cd /usr/local/src/joinx && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make deps && \
    make && \
    make install


# Install PERL modules

RUN cpanm Test::Most && \
    cpanm Statistics::Descriptive && \
    cpanm Statistics::Distributions && \
    cpanm Bit::Vector



# Install Genome MuSiC2.
# http://gmt.genome.wustl.edu/packages/genome-music/install.html
# https://github.com/ding-lab/MuSiC2

RUN cd /usr/local/src && \
    git clone https://github.com/ding-lab/MuSiC2 && \
    cd /usr/local/src/MuSiC2 && \
    cpanm MuSiC2-0.2.tar.gz





# Clean up apt-get.

RUN cd / && \
    apt-get remove --purge -y $BUILDDEPS &&\
    apt-get autoremove -y && \
    apt-get autoclean -y && \
    rm -rf /var/lib/apt/lists/*

CMD ["music2"]

