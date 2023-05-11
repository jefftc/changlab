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


RUN apt-get update && apt-get install -y --no-install-recommends \
#    build-essential \
#    gcc \
#    make \
    libbz2-dev \
    zlib1g-dev

#    liblzma-dev \
#RUN apt-get install -y --no-install-recommends \
#    bzip2 \
#    libncurses5-dev


RUN cd /usr/local/src && \
  curl -L -O https://github.com/freebayes/freebayes/releases/download/v1.3.4/freebayes-1.3.4-linux-static-AMD64.gz && \
  gunzip freebayes-1.3.4-linux-static-AMD64.gz && \
  chmod 755 freebayes-1.3.4-linux-static-AMD64 && \
  ln -s /usr/local/src/freebayes-1.3.4-linux-static-AMD64 \
    /usr/local/bin/freebayes


## Samtools
#
#ENV HTSLIB_VERSION=1.10.2
#RUN cd /usr/local/src && \
#    curl -L -O https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 && \
#    bzcat htslib-${HTSLIB_VERSION}.tar.bz2 | tar xf - && \
#    rm -f htslib-${HTSLIB_VERSION}.tar.bz2
#
#RUN cd /usr/local/src/htslib-${HTSLIB_VERSION} && \
#    ./configure && \
#    make && \
#    make install
#
#ENV SAMTOOLS_VERSION=1.10
#RUN cd /usr/local/src && \
#    curl -L -O https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
#    bzcat samtools-${SAMTOOLS_VERSION}.tar.bz2 | tar xf - && \
#    rm -f samtools-${SAMTOOLS_VERSION}.tar.bz2
#
#RUN cd /usr/local/src/samtools-${SAMTOOLS_VERSION} && \
#    ./configure && \
#    make && \
#    make install
#
#
## vcftools
#
#RUN apt-get install -y --no-install-recommends \
#    git 
#
#RUN apt-get install -y --no-install-recommends \
#    libvcflib-tools
#
#RUN apt-get install -y --no-install-recommends \
#    pkg-config \
#    autoconf \
#    automake \
#    libtool
#
#RUN cd /usr/local/src && \
#    git clone https://github.com/vcftools/vcftools.git && \
#    cd vcftools && \
#    ./autogen.sh && \
#    ./configure && \
#    make && \
#    make install
#
#
#
### tabixpp
##
##RUN cd /usr/local/src && \
##  git clone --recursive https://github.com/ekg/tabixpp.git
### Bug: -lcurl needs to be included in LIBS or will get errors.
##RUN cd /usr/local/src/tabixpp && \
##  LIBS="-lhts -lpthread -lm -lbz2 -llzma -lz -lcurl" make
#
#
## vcflib
#
#RUN apt-get install -y --no-install-recommends \
#    libtabixpp-dev
#
#RUN apt-get install -y --no-install-recommends \
#    cmake
#
#RUN cd /usr/local/src && \
#  git clone --recursive https://github.com/vcflib/vcflib.git
#
#RUN cd /usr/local/src/vcflib && \
#  mkdir -p build && cd build && \
#  cmake ..
#RUN cd /usr/local/src/vcflib/build && \
#  cmake --build .
#RUN cd /usr/local/src/vcflib/build && \
#  cmake --install .
#  
#
## freebayes
#
#RUN apt-get install -y --no-install-recommends \
#    bc \
#    parallel \
#    meson \
#    ninja-build
#
#
#https://github.com/freebayes/freebayes/releases/download/v1.3.4/freebayes-1.3.4-linux-static-AMD64.gz
#gunzip freebayes-1.3.4-linux-static-AMD64.gz
#
#RUN cd /usr/local/src && \
#    git clone --recursive git://github.com/ekg/freebayes.git
#
#RUN apt-get install -y --no-install-recommends \
#    libvcflib-dev
#
##RUN cd /usr/local/src/freebayes && \
##    meson build/ --buildtype debug
##    make && \
##    make install
#
#
##Dependency libvcflib found: NO
##Dependency libseqlib found: NO


  




CMD ["freebayes"]
