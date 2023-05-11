# Heavily modified from rocker-versioned Dockerfile.

# /usr/local/MISO/
#   hg19_v1/
#     SE.hg19.gff3
#     [...]
#   mm10_v1/
#     [...]
#   hg19_v2/
#     [...]
#   mm10_v2/
#     [...]
# /usr/local/bin/
#   gtf2gff3.pl
#   summarize_miso
#   compare_miso
#   sashimi_plot
#   pe_utils
#   exon_utils
#   index_gff


# TODO:
# - chmod +x gtf2gff3.pl
# - gtf2gff3.pl gets warning (but still runs):
#   Config::Std not installed. Will use built-in defaults.

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


# Download the MISO event files.
# http://miso.readthedocs.io/en/fastmiso/annotation.html

RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    curl \
    unzip


ENV F=miso_annotations_mm10_v1.zip
RUN mkdir -p /usr/local/MISO && \
  cd /usr/local/MISO && \
  curl -O http://hollywood.mit.edu/burgelab/miso/annotations/${F} && \
  unzip ${F} && \
  mv mm10 mm10_v1 && \
  rm ${F}

ENV F=miso_annotations_hg19_v1.zip
RUN mkdir -p /usr/local/MISO && \
  cd /usr/local/MISO && \
  curl -O http://hollywood.mit.edu/burgelab/miso/annotations/${F} && \
  unzip ${F} && \
  mv hg19 hg19_v1 && \
  rm ${F}


ENV F=miso_annotations_mm10_v2.zip
RUN mkdir -p /usr/local/MISO && \
  cd /usr/local/MISO && \
  curl -O http://hollywood.mit.edu/burgelab/miso/annotations/ver2/${F} && \
  unzip ${F} && \
  mv mm10 mm10_v2 && \
  rm ${F}

ENV F=miso_annotations_hg19_v2.zip
RUN mkdir -p /usr/local/MISO && \
  cd /usr/local/MISO && \
  curl -O http://hollywood.mit.edu/burgelab/miso/annotations/ver2/${F} && \
  unzip ${F} && \
  mv hg19 hg19_v2 && \
  rm ${F}




# Python
RUN apt-get update && apt-get install -y --no-install-recommends \
    python2.7-dev \
    python-pip \
    python-tk     # needed for matplotlib

RUN pip install setuptools wheel
RUN pip install numpy

# PySam is a bit complicated and has a lot of dependencies.
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gcc \
    make \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev

# Cython needs to be installed before pysam
RUN pip install cython
RUN pip install pysam

RUN pip install misopy

# Must be installed after gcc.
RUN pip install matplotlib

# misopy requires tagBam from bedtools.
RUN apt-get install -y --no-install-recommends \
    git

RUN cd /usr/local/src && \
    git clone https://github.com/arq5x/bedtools2.git

RUN cd /usr/local/src/bedtools2 && \
    make && \
    make install


# tagBam requires samtools.
RUN apt-get install -y --no-install-recommends \
    bzip2 \
    libncurses5-dev

ENV HTSLIB_VERSION=1.10.2
RUN cd /usr/local/src && \
    curl -L -O https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 && \
    bzcat htslib-${HTSLIB_VERSION}.tar.bz2 | tar xf - && \
    rm -f htslib-${HTSLIB_VERSION}.tar.bz2

RUN cd /usr/local/src/htslib-${HTSLIB_VERSION} && \
    ./configure && \
    make && \
    make install

ENV SAMTOOLS_VERSION=1.10
RUN cd /usr/local/src && \
    curl -L -O https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    bzcat samtools-${SAMTOOLS_VERSION}.tar.bz2 | tar xf - && \
    rm -f samtools-${SAMTOOLS_VERSION}.tar.bz2

RUN cd /usr/local/src/samtools-${SAMTOOLS_VERSION} && \
    ./configure && \
    make && \
    make install


# Install gtf2gff3.pl.
RUN cd /usr/local/bin && \
  curl -O http://hollywood.mit.edu/burgelab/miso/scripts/gtf2gff3.pl


CMD ["/usr/local/bin/miso"]
