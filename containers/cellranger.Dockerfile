# Heavily modified from rocker-versioned Dockerfile.

# /usr/local/cellranger/
#   cellranger
# /usr/local/bcl2fastq/


# No.  This makes the container too big.
## /usr/local/genomes/
##   refdata-cellranger-GRCh38-3.0.0/
##   refdata-cellranger-GRCh38-3.0.0_premrna/
##   refdata-cellranger-GRCh38-and-mm10-3.1.0/
##   refdata-cellranger-GRCh38-and-mm10-3.1.0_premrna/
##   refdata-cellranger-hg19-3.0.0/
##   refdata-cellranger-hg19-3.0.0_premrna/
##   refdata-gex-GRCh38-2020-A/
##   refdata-gex-GRCh38-and-mm10-2020-A/
##   refdata-gex-mm10-2020-A/



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



## Copy the genome files.
#
#RUN mkdir /usr/local/genomes
#COPY refdata-cellranger-hg19-3.0.0 /usr/local/genomes
#COPY refdata-cellranger-hg19-3.0.0_premrna /usr/local/genomes
#COPY refdata-cellranger-GRCh38-3.0.0 /usr/local/genomes
#COPY refdata-cellranger-GRCh38-3.0.0_premrna /usr/local/genomes
#COPY refdata-cellranger-GRCh38-and-mm10-3.1.0 /usr/local/genomes
#COPY refdata-cellranger-GRCh38-and-mm10-3.1.0_premrna /usr/local/genomes
#COPY refdata-gex-GRCh38-2020-A /usr/local/genomes
#COPY refdata-gex-GRCh38-and-mm10-2020-A /usr/local/genomes
#COPY refdata-gex-mm10-2020-A /usr/local/genomes



# Install Illumina bcl2fastq.

RUN apt-get update && apt-get install -y --no-install-recommends \
  unzip
  
RUN apt-get update && apt-get install -y --no-install-recommends \
  gcc \
  cmake \
  libboost-dev \
  zlib1g-dev \
  libpthread-stubs0-dev

COPY bcl2fastq2-v2-20-0-tar.zip /usr/local/src
RUN cd /usr/local/src && \
  unzip bcl2fastq2-v2-20-0-tar.zip

RUN apt-get update && apt-get install -y --no-install-recommends \
  g++ \
  make \
  bzip2

ENV TMP=/tmp
ENV SOURCE=${TMP}/bcl2fastq \
  BUILD=${TMP}/bcl2fastq2-v2-20-0-build \
  INSTALL_DIR=/usr/local/bcl2fastq2-v2-20-0

# Configuration gets error:
# CMake Error at cmake/macros.cmake:80 (message):
#   Required header sys/stat.h not found.
# Fix it by linking file to expected place.
RUN ln -s /usr/include/x86_64-linux-gnu/sys /usr/include/

RUN cd ${TMP} && \
  tar xfz /usr/local/src/bcl2fastq2-v2.20.0.422-Source.tar.gz && \
  mkdir ${BUILD} && \
  cd ${BUILD} && \
  ${SOURCE}/src/configure --prefix=${INSTALL_DIR}

RUN cd ${BUILD} && \
  make && \
  make install



# Install Cell Ranger.

COPY cellranger-6.0.2.tar.gz /usr/local/src

RUN cd /usr/local/src && \
  tar xfz cellranger-6.0.2.tar.gz

RUN ln -s /usr/local/src/cellranger-6.0.2 /usr/local/cellranger
RUN /usr/local/cellranger/cellranger sitecheck > /usr/local/sitecheck.txt

# TODO: Move this up.
RUN ln -s /usr/local/bcl2fastq2-v2-20-0 /usr/local/bcl2fastq

ENV PATH $PATH:/usr/local/cellranger:/usr/local/bcl2fastq/bin



# Clean up apt-get.
RUN cd / && \
    apt-get remove --purge -y $BUILDDEPS &&\
    apt-get autoremove -y && \
    apt-get autoclean -y && \
    rm -rf /var/lib/apt/lists/*

CMD ["cellranger"]
