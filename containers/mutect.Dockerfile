# Heavily modified from rocker-versioned Dockerfile.


# /usr/lib/jvm/jre1.6.0_45/
#   bin/java
# /opt/muTect/
#   muTect-1.1.4.jar
#   data.hg19/
#     dbsnp_132_b37.leftAligned.vcf
#     dbsnp_132_b37.leftAligned.vcf.idx
#     hg19_cosmic_v54_120711.vcf
#     hg19_cosmic_v54_120711.vcf.idx
##   data.numeric.hg19/            Sorted numerically
##     dbsnp_132_b37.leftAligned.vcf
##     dbsnp_132_b37.leftAligned.vcf.idx
##     hg19_cosmic_v54_120711.vcf
##     hg19_cosmic_v54_120711.vcf.idx
##   data.mm9/
##     dbsnp_128_mm9.vcf           MISSING.  Can't find online.
#   muTect-1.1.4-bin.zip



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
    vim \
    openssl \
    libssl-dev \
    curl \
    libcurl4-openssl-dev


# Compression
RUN apt-get update && apt-get install -y --no-install-recommends \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    unzip \
    zip


# Install muTect.

ENV MUTECT_DIR=/opt/muTect
COPY muTect-1.1.4-bin.zip ${MUTECT_DIR}/
RUN cd ${MUTECT_DIR} && \
    unzip muTect-1.1.4-bin.zip


# Install muTect files.

COPY b37_cosmic_v54_120711.vcf.* ${MUTECT_DIR}/data.hg19/
COPY dbsnp_132_b37.leftAligned.vcf.* ${MUTECT_DIR}/data.hg19/
RUN cd ${MUTECT_DIR}/data.hg19 && \
    gunzip *.gz



# Requires Java 6.

## # Only OpenJDK-8 available.  Install it.
## RUN apt-get update && apt-get install -y --no-install-recommends \
##     openjdk-8-jdk \
##     ant
## ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/
## RUN export JAVA_HOME=/usr/lib/jvm/java-8-openjdk-amd64/

COPY jre-6u45-linux-x64.bin /usr/lib/jvm/
RUN cd /usr/lib/jvm && \
    sh jre-6u45-linux-x64.bin && \
    ln -s /usr/lib/jvm/jre1.6.0_45/bin/java /usr/local/bin/java && \
    export JAVA_HOME=/usr/lib/jvm/jre1.6.0_45/
ENV JAVA_HOME /usr/lib/jvm/jre1.6.0_45/


## # Sort the VCF files.
## 
## RUN mkdir -p ${MUTECT_DIR}/data.numeric.hg19
## 
## ENV F1=${MUTECT_DIR}/data.hg19/dbsnp_132_b37.leftAligned.vcf \
##   F2=${MUTECT_DIR}/data.numeric.hg19/dbsnp_132_b37.leftAligned.vcf
## RUN grep "^#" ${F1} > ${F2} && \
##   grep -v "^#" ${F1} | sort -V -k1,1 -k2,2n >> ${F2}
## 
## ENV F1=${MUTECT_DIR}/data.hg19/hg19_cosmic_v54_120711.vcf \
##   F2=${MUTECT_DIR}/data.numeric.hg19/hg19_cosmic_v54_120711.vcf
## RUN grep "^#" ${F1} > ${F2} && \
##   grep -v "^#" ${F1} | sort -V -k1,1 -k2,2n >> ${F2}


# Helpful for copying files out.
RUN apt-get update && apt-get install -y --no-install-recommends \
    rsync

