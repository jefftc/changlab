# Modified from:
# https://github.com/GenomicParisCentre/dockerfiles/
#   blob/master/tophat2/2.0.14/Dockerfile
FROM ubuntu:14.04

LABEL maintainer="Jeffrey Chang <jeffrey.t.chang@uth.tmc.edu>"

ENV BUILD_DATE="$(TZ='America/Los_Angeles' date -I)" \
    LC_ALL=en_US.UTF-8 \
    LANG=en_US.UTF-8 \
    TERM=xterm

# Basic stuff.
RUN apt-get update
RUN apt-get install -y --no-install-recommends \
    bash-completion \
    ca-certificates \
    locales
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen && \
    /usr/sbin/locale-gen en_US.UTF-8 && \
    /usr/sbin/update-locale LANG=en_US.UTF-8


# Install tophat2.

RUN apt-get install -y --no-install-recommends \
    make \
    g++ \
    wget \
    unzip \
    python


RUN cd /usr/local/src && \
    wget http://ccb.jhu.edu/software/tophat/downloads/tophat-2.0.14.Linux_x86_64.tar.gz && \
    tar zxvf tophat-2.0.14.Linux_x86_64.tar.gz && \
    rm tophat-2.0.14.Linux_x86_64.tar.gz

RUN ln -s /usr/local/src/tophat-2.0.14.Linux_x86_64/* /usr/local/bin/ && \
    rm -f /usr/local/bin/AUTHORS && \
    rm -f /usr/local/bin/COPYING && \
    rm -f /usr/local/bin/README


# Install bowtie2

RUN cd /usr/local/src && \
    wget --default-page=bowtie2-2.2.3-linux-x86_64.zip http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.3/bowtie2-2.2.3-linux-x86_64.zip/ && \
    unzip bowtie2-2.2.3-linux-x86_64.zip && \
    rm bowtie2-2.2.3-linux-x86_64.zip

RUN ln -s /usr/local/src/bowtie2-2.2.3/bowtie* /usr/local/bin/ && \
    ln -s /usr/local/src/bowtie2-2.2.3/bowtie2 /usr/local/bin/bowtie

#ENV PATH $PATH:/bin/tophat-2.0.14.Linux_x86_64
#ENV PATH $PATH:/bin/bowtie2-2.2.3



# Install gffread

RUN cd /usr/local/src && \
    wget http://ccb.jhu.edu/software/stringtie/dl/gffread-0.11.6.tar.gz && \
    tar xvfz gffread-0.11.6.tar.gz && \
    rm -f gffread-0.11.6.tar.gz && \
    cd gffread-0.11.6 && \
    make

RUN ln -s /usr/local/src/gffread-0.11.6/gffread /usr/local/bin

# NEED TO INCLUDE /usr/local/bin IN PATH.
