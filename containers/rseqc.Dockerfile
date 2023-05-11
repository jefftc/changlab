# Heavily modified from rocker-versioned Dockerfile.
FROM debian:stretch


# Paths:
# /usr/local/RSeQC/     Contains downloaded data files
# /usr/local/bin/
#   <scripts>


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


# Download data files for RSeQC.
# Do this early, so we can change installation below without
# re-downloading these files.

RUN mkdir /usr/local/RSeQC

## Download testing data sets.
#
##Pair-end strand specific (Illumina). BAM file md5sum=fbd1fb1c153e3d074524ec70e6e21fb9
#RUN curl -O http://dldcc-web.brc.bcm.edu/lilab/liguow/RSeQC/dat/Pairend_StrandSpecific_51mer_Human_hg19.bam
#RUN curl -O http://dldcc-web.brc.bcm.edu/lilab/liguow/RSeQC/dat/Pairend_StrandSpecific_51mer_Human_hg19.bam.bai
#
##Pair-end  non-strand specific (Illumina). BAM file md5sum=ba014f6b397b8a29c456b744237a12de
#RUN curl -O http://dldcc-web.brc.bcm.edu/lilab/liguow/RSeQC/dat/Pairend_nonStrandSpecific_36mer_Human_hg19.bam
#RUN curl -O http://dldcc-web.brc.bcm.edu/lilab/liguow/RSeQC/dat/Pairend_nonStrandSpecific_36mer_Human_hg19.bam.bai
#
##Single-end strand specific (SOLiD). BAM file md5sum=b39951a6ba4639ca51983c2f0bf5dfce
#RUN curl -O http://dldcc-web.brc.bcm.edu/lilab/liguow/RSeQC/dat/SingleEnd_StrandSpecific_50mer_Human_hg19.bam
#RUN curl -O http://dldcc-web.brc.bcm.edu/lilab/liguow/RSeQC/dat/SingleEnd_StrandSpecific_50mer_Human_hg19.bam.bai

# Download gene models (update to 06/27/2013)

ENV BCM=http://dldcc-web.brc.bcm.edu/lilab/liguow/RSeQC/dat

#human (hg19/GRCh37)
RUN cd /usr/local/RSeQC && \
  curl -O ${BCM}/hg19_RefSeq.bed.gz && \
  curl -O ${BCM}/hg19_Ensembl.bed.gz && \
  curl -O ${BCM}/hg19_GENCODE_v14.bed.gz && \
  curl -O ${BCM}/hg19_GENCODE_v12.bed.gz && \
  curl -O ${BCM}/hg19_UCSC_knownGene.bed.gz && \
  curl -O ${BCM}/hg19_Vega.bed.gz && \
  curl -O ${BCM}/hg19_AceView.bed.gz

#Mouse (mm9)
RUN cd /usr/local/RSeQC && \
  curl -O ${BCM}/mm9_NCBI37_Ensembl.bed.gz && \
  curl -O ${BCM}/mm9_NCBI37_MGC.bed.gz && \
  curl -O ${BCM}/mm9_NCBI37_Refseq.bed.gz

#Mouse (mm10)
RUN cd /usr/local/RSeQC && \
  curl -O ${BCM}/GRCm38_mm10_Ensembl.bed.gz && \
  curl -O ${BCM}/GRCm38_mm10_MGC.bed.gz && \
  curl -O ${BCM}/GRCm38_mm10_RefSeq.bed.gz

#Fly (D. melanogaster) (BDGP R5/dm3)
RUN cd /usr/local/RSeQC && \
  curl -O ${BCM}/fly_dm3_EnsemblGene.bed.gz && \
  curl -O ${BCM}/fly_dm3_RefSeq.bed.gz && \
  curl -O ${BCM}/fly_dm3_flyBaseGene.bed.gz

# Download ribosome RNA (update to 08/17/2012)
ENV GOOGLE=https://sites.google.com/site/liguowangspublicsite/home
RUN cd /usr/local/RSeQC && \
  curl -O ${GOOGLE}/hg19_rRNA.bed && \
  curl -O ${GOOGLE}/mm10_rRNA.bed && \
  curl -O ${GOOGLE}/mm9_rRNA.bed


# Fetch chromosome sizes.

RUN cd /usr/local/bin && \
  curl -O http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes && \
  chmod 755 fetchChromSizes

RUN cd /usr/local/RSeQC && \
  fetchChromSizes hg19 > hg19.chrom.sizes && \
  fetchChromSizes mm9 > mm9.chrom.sizes && \
  fetchChromSizes mm10 > mm10.chrom.sizes && \
  fetchChromSizes danRer7 > danRer7.chrom.sizes 





# Developer stuff.
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gcc \
    make


# Python
RUN apt-get update && apt-get install -y --no-install-recommends \
    python2.7-dev \
    python-pip


RUN pip install setuptools wheel
RUN pip install numpy bx-python

# PySam is a bit complicated.
RUN apt-get update && apt-get install -y --no-install-recommends \
    libbz2-dev    \
    liblzma-dev \
    zlib1g-dev
# Cython needs to be installed before pysam.
RUN pip install Cython && \
  pip install pysam

# Install RSeQC
RUN pip install RSeQC


#CMD ["python"]
