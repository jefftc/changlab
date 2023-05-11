# Heavily modified from rocker-versioned Dockerfile.

# /usr/local/NGSCheckMate/
#   ncm.py
#   ncm.conf
#   SNP/
#     SNP_GRCh37_hg19_wChr.bed
#     SNP_GRCh37_hg19_woChr.bed
#     SNP_GRCh38_hg38_wChr.bed

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



# samtools

RUN apt-get install -y --no-install-recommends \
    bzip2 \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev

ENV HTSLIB_VERSION=1.10.2
RUN cd /usr/local/src && \
    curl -L -O https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 && \
    bzcat htslib-${HTSLIB_VERSION}.tar.bz2 | tar xf - && \
    rm -f htslib-${HTSLIB_VERSION}.tar.bz2

RUN cd /usr/local/src/htslib-${HTSLIB_VERSION} && \
    ./configure && \
    make && \
    make install


RUN apt-get install -y --no-install-recommends \
    libncurses5-dev

ENV SAMTOOLS_VERSION=1.10
RUN cd /usr/local/src && \
    curl -L -O https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    bzcat samtools-${SAMTOOLS_VERSION}.tar.bz2 | tar xf - && \
    rm -f samtools-${SAMTOOLS_VERSION}.tar.bz2

RUN cd /usr/local/src/samtools-${SAMTOOLS_VERSION} && \
    ./configure && \
    make && \
    make install

ENV BCFTOOLS_VERSION=1.10.2
RUN cd /usr/local/src && \
    curl -L -O https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    bzcat bcftools-${HTSLIB_VERSION}.tar.bz2 | tar xf - && \
    rm -f bcftools-${HTSLIB_VERSION}.tar.bz2
RUN cd /usr/local/src/bcftools-${HTSLIB_VERSION} && \
    ./configure && \
    make && \
    make install



# Install R 3.3.3.
RUN apt-get update && apt-get install -y --no-install-recommends \
    r-base
RUN apt-get install -y --no-install-recommends \
    libreadline-dev \
    libblas-dev \
    liblapack-dev \
    libicu-dev


# Add a default CRAN mirror
RUN mkdir -p /usr/lib/R/etc && \
    echo "options(repos=c(CRAN='https://cran.rstudio.com/'), download.file.method='libcurl')" \
    >> /usr/lib/R/etc/Rprofile.site

# Add a library directory (for user-installed packages)
RUN mkdir -p /usr/lib/R/site-library && \
    chown root:staff /usr/lib/R/site-library && \
    chmod g+wx /usr/lib/R/site-library

# Fix library path
# spython messes up when parsing the string "{U}SER".  Will interpret as
# the {U}SER variable.  So need to munge the string.
ENV A=USE
RUN echo "R_LIBS_${A}R='/usr/lib/R/site-library'" >> /usr/lib/R/etc/Renviron && \
  echo "R_LIBS=\${R_LIBS-'/usr/lib/R/site-library:/usr/lib/R/library'}" >> /usr/lib/R/etc/Renviron



# NGSCheckMate

RUN cd /usr/local/ && \
    git clone https://github.com/parklab/NGSCheckMate.git
ENV NCM_HOME=/usr/local/NGSCheckMate

#RUN ln -s /usr/local/NGSCheckMate/ncm.py /usr/local/bin/ && \
#  ln -s /usr/local/NGSCheckMate/ncm_fastq.py /usr/local/bin/ && \
#  ln -s /usr/local/NGSCheckMate/ngscheckmate_fastq /usr/local/bin/


# To use BAM module, need to set in /usr/local/NGSCheckMate/ncm.conf:
# REF=<path for the reference FASTA file >  




# Clean up apt-get.
RUN cd / && \
    apt-get remove --purge -y $BUILDDEPS &&\
    apt-get autoremove -y && \
    apt-get autoclean -y && \
    rm -rf /var/lib/apt/lists/*

#CMD ["ncm.py"]
