# Heavily modified from rocker-versioned Dockerfile.

# /usr/local/bin/bam-somaticsniper


FROM debian:bullseye

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


# Install build tools.

RUN apt-get update && apt-get install -y --no-install-recommends \
  build-essential \
  gcc \
  make \
  cmake

RUN apt-get update && apt-get install -y --no-install-recommends \
  zlib1g-dev \
  libncurses5-dev

RUN apt-get update && apt-get install -y --no-install-recommends \
  git

# Unit tests are written in python.
RUN apt-get update && apt-get install -y --no-install-recommends \
  python2.7
RUN ln -s /usr/bin/python2.7 /usr/bin/python


# Install SomaticSniper.

RUN cd /usr/local/src/ && \
  git clone git://github.com/genome/somatic-sniper.git
RUN cd /usr/local/src/somatic-sniper && \
  mkdir -p build && \
  cd build && \
  cmake ../ && \
  make deps && \
  make -j

RUN cd /usr/local/src/somatic-sniper/build && \
  make test

RUN cd /usr/local/src/somatic-sniper/build && \
  make install


# Clean up apt-get.
RUN cd / && \
    apt-get remove --purge -y $BUILDDEPS &&\
    apt-get autoremove -y && \
    apt-get autoclean -y && \
    rm -rf /var/lib/apt/lists/*

CMD ["bam-somaticsniper"]

