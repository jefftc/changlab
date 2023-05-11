# Heavily modified from rocker-versioned Dockerfile.

# /usr/local/snpEff/
#   snpEff.jar
#   databases.txt
#   data/
#     hg19/
#     mm10/
#     GRCh37.75/


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




# Install OpenJDK-8
RUN apt-get install -y --no-install-recommends \
    openjdk-8-jdk \
    ant

ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/
RUN export JAVA_HOME=/usr/lib/jvm/java-8-openjdk-amd64/


# SnpEff

RUN apt-get install -y --no-install-recommends \
  unzip

RUN cd /usr/local && \
  curl -L -O https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip && \
  unzip snpEff_latest_core.zip && \
  rm -f snpEff_latest_core.zip

RUN cd /usr/local/snpEff && \
  java -jar snpEff.jar databases > databases.txt

RUN cd /usr/local/snpEff && \
  java -jar snpEff.jar download hg19
RUN cd /usr/local/snpEff && \
  java -jar snpEff.jar download mm10
# Download genome used by Radia.
RUN cd /usr/local/snpEff && \
  java -jar snpEff.jar download GRCh37.75

  
# Clean up apt-get.
RUN cd / && \
    apt-get remove --purge -y $BUILDDEPS &&\
    apt-get autoremove -y && \
    apt-get autoclean -y && \
    rm -rf /var/lib/apt/lists/*

#CMD ["R"]



