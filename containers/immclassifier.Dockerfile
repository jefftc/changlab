# Heavily modified from rocker-versioned Dockerfile.

# /bin/
#   zcat
#   bzcat
# /usr/bin/
#   xzcat
# /usr/local/bin/
#   bwa


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



# Download BWA.
ENV A=bwa-0.7.17.tar.bz2
RUN cd /usr/local/src && \
  curl -L https://sourceforge.net/projects/bio-bwa/files/${A}/download > ${A} \
    && \
  bzcat ${A} | tar xf - && \
  rm ${A}


# BWA dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    zlib1g-dev


# Compile and install
RUN cd /usr/local/src/bwa-0.7.17 && \
  make && \
  cp -p bwa /usr/local/bin/
  

# Clean up apt-get.
RUN cd / && \
    apt-get remove --purge -y $BUILDDEPS &&\
    apt-get autoremove -y && \
    apt-get autoclean -y && \
    rm -rf /var/lib/apt/lists/*

CMD ["bwa"]
