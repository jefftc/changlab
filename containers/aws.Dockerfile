# Heavily modified from rocker-versioned Dockerfile.

# /usr/local/bin/
#   aws

FROM debian:bullseye

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
    openssl \
    libssl-dev \
    libcurl4-openssl-dev \
    curl

# AWS CLI.

RUN apt-get update && apt-get install -y --no-install-recommends \
  unzip

RUN cd /usr/local/src && \
  curl -L -O "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" && \
  unzip awscli-exe-linux-x86_64.zip && \
  ./aws/install && \
  rm -f awscli-exe-linux-x86_64.zip && \
  rm -rf aws


# Needed for aws help
RUN apt-get update && apt-get install -y --no-install-recommends \
  groff

RUN apt-get update && apt-get install -y --no-install-recommends \
    less \
    file \
    vim

# Clean up apt-get.
RUN cd / && \
    apt-get remove --purge -y $BUILDDEPS &&\
    apt-get autoremove -y && \
    apt-get autoclean -y && \
    rm -rf /var/lib/apt/lists/*


CMD ["aws"]
