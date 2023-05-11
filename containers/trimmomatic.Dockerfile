# Heavily modified from rocker-versioned Dockerfile.


# /usr/local/bin/trimmomatic-0.39.jar
# /usr/local/src/Trimmomatic-0.39/adapters/


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
    apt-utils \
    bash-completion \
    ca-certificates \
    locales
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen && \
    /usr/sbin/locale-gen en_US.UTF-8 && \
    /usr/sbin/update-locale LANG=en_US.UTF-8


# Miscellaneous useful things.
RUN apt-get update && apt-get install -y --no-install-recommends \
    file \
    less \
    vim \
    openssl \
    libssl-dev \
    libcurl4-openssl-dev \
    curl


# Install OpenJDK-8
RUN apt-get update && apt-get install -y --no-install-recommends \
    openjdk-8-jdk

ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/
RUN export JAVA_HOME=/usr/lib/jvm/java-8-openjdk-amd64/



# Trimmomatic.

RUN apt-get update && apt-get install -y --no-install-recommends \
    unzip

RUN cd /usr/local/src && \
    curl -O http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip && \
    unzip Trimmomatic-0.39.zip && \
    cp -p /usr/local/src/Trimmomatic-0.39/trimmomatic-0.39.jar /usr/local/bin
RUN ln -s /usr/local/bin/trimmomatic-0.39.jar /usr/local/bin/trimmomatic.jar



# Uncompression programs

RUN apt-get update && apt-get install -y --no-install-recommends \
    bzip2 \
    xz-utils
