# Heavily modified from rocker-versioned Dockerfile.

# /usr/local/bin/
#   bedGraphToBigWig
#   bedSort
#   bedToBigBed
#   fetchChromSizes
#   gtfToGenePred

FROM debian:buster

LABEL maintainer="Jeffrey Chang <jeffrey.t.chang@uth.tmc.edu>"

ENV BUILD_DATE="$(TZ='America/Los_Angeles' date -I)" \
    LC_ALL=en_US.UTF-8 \
    LANG=en_US.UTF-8 \
    TERM=xterm

# Basic stuff.
RUN apt-get update
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
    bash-completion \
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



# kentUtils

#RUN apt-get update && apt-get install -y --no-install-recommends \
#    rsync

#RUN rsync -azvP \
#  rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/ \
#  /usr/local/bin/

RUN apt-get install -y --no-install-recommends \
    libpng-dev \
    uuid-dev
    #default-libmysqlclient-dev

RUN apt-get install -y --no-install-recommends \
    mariadb-server

ENV MACHTYPE=x86_64
RUN mkdir -p /root/bin/${MACHTYPE}

RUN cd /usr/local/src && \
    curl -k -L -O https://github.com/ucscGenomeBrowser/kent/archive/v377_base.tar.gz && \
    tar xvfz v377_base.tar.gz && \
    rm -f v377_base.tar.gz

RUN cd /usr/local/src/kent-377_base/src/lib && \
    make
RUN cd /usr/local/src/kent-377_base/src/jkOwnLib && \
    make
RUN cd /usr/local/src/kent-377_base/src/htslib && \
    make

# Need jkhgap.a for bedSort.
# This requires mysql_config for some reason.
RUN apt-get install -y --no-install-recommends \
    libmariadbclient-dev
RUN apt-get install -y --no-install-recommends \
    default-libmysqlclient-dev
RUN cd /usr/local/src/kent-377_base/src && \
    make libs

RUN cd /usr/local/src/kent-377_base/src/utils/bedToBigBed && \
    make && \
    cp -p /root/bin/${MACHTYPE}/bedToBigBed /usr/local/bin/
RUN cd /usr/local/src/kent-377_base/src/utils/bedGraphToBigWig && \
    make && \
    cp -p /root/bin/${MACHTYPE}/bedGraphToBigWig /usr/local/bin/
RUN cd /usr/local/src/kent-377_base/src/hg/bedSort && \
    make && \
    cp -p /root/bin/${MACHTYPE}/bedSort /usr/local/bin/
RUN cd /usr/local/src/kent-377_base/src/hg/utils/gtfToGenePred && \
    L="${L} -lgnutls" make && \
    cp -p /root/bin/${MACHTYPE}/gtfToGenePred /usr/local/bin/
RUN cd /usr/local/src/kent-377_base/src/utils/userApps/ && \
    make && \
    cp -p /root/bin/${MACHTYPE}/fetchChromSizes /usr/local/bin/




# Save space.
RUN rm -rf /usr/local/src/kent-377_base

