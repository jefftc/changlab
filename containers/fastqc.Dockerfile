# Heavily modified from rocker-versioned Dockerfile.

# /usr/local/bin/
#   fastqc


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
    file \
    less \
    vim

# \
#    libcurl4-openssl-dev \
#    curl

RUN apt-get install -y --no-install-recommends \
    unzip \
    gnupg \
    wget




# Install OpenJDK-8 for Buster.
# Need to add the adoptopenjdk repository.

# For add-apt-repository.
RUN apt-get install -y --no-install-recommends \
    software-properties-common
RUN wget -qO - https://adoptopenjdk.jfrog.io/adoptopenjdk/api/gpg/key/public \
   | apt-key add -
RUN add-apt-repository --yes https://adoptopenjdk.jfrog.io/adoptopenjdk/deb/
RUN apt-get update && apt-get install -y --no-install-recommends \
    adoptopenjdk-8-hotspot


# FastQC
# https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc

ENV JAVA_HOME /usr/lib/jvm/adoptopenjdk-8-hotspot-amd64/
RUN export JAVA_HOME=/usr/lib/jvm/java-8-openjdk-vamd64/

RUN cd /usr/local/src && \
    wget --max-redirect 4 https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip && \
    unzip fastqc_v0.11.9.zip && \
    chmod 755 FastQC/fastqc && \
    ln -s /usr/local/src/FastQC/fastqc /usr/local/bin/ && \
    rm -rf fastqc_v0.11.9.zip


# Perl libraries.

RUN apt-get install -y --no-install-recommends \
    perl
RUN cpan FindBin


RUN apt-get install -y --no-install-recommends \
    libfontconfig1



CMD ["/usr/local/bin/fastqc"]
