# Heavily modified from rocker-versioned Dockerfile.

# /usr/local/illumina/manifest

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
    bash-completion \
    ca-certificates \
    locales
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen && \
    /usr/sbin/locale-gen en_US.UTF-8 && \
    /usr/sbin/update-locale LANG=en_US.UTF-8


RUN apt-get update && apt-get install -y --no-install-recommends \
    python \
    rsync

## Miscellaneous useful things.
## Don't know why, but apt-get update here prevents 404.
#RUN apt-get update && apt-get install -y --no-install-recommends \
#    file \
#    less \
#    vim \
#    openssl \
#    libssl-dev \
#    libcurl4-openssl-dev \
#    curl


# Install illumina manifests.

COPY illumina /usr/local/illumina/manifest


CMD ["bash"]
