# Heavily modified from rocker-versioned Dockerfile.

# /usr/local/bin/
#   bowtie
#   bowtie-align-l
#   bowtie-align-l-debug
#   bowtie-align-s
#   bowtie-align-s-debug
#   bowtie-build
#   bowtie-build-l
#   bowtie-build-l-debug
#   bowtie-build-s
#   bowtie-build-s-debug
#   bowtie-inspect
#   bowtie-inspect-l
#   bowtie-inspect-l-debug
#   bowtie-inspect-s
#   bowtie-inspect-s-debug

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
    vim

RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    curl

RUN apt-get install -y --no-install-recommends \
    unzip



# Install bowtie.

ENV A=bowtie-1.3.1-linux-x86_64
RUN cd /usr/local/src && \
    curl -L https://sourceforge.net/projects/bowtie-bio/files/${A}.zip/download > ${A}.zip
RUN cd /usr/local/src && \
    unzip ${A}.zip && \
    rm ${A}.zip
RUN cp -p /usr/local/src/${A}/bowtie* /usr/local/bin


# Requires python3
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3


# Clean up apt-get.
RUN cd / && \
    apt-get remove --purge -y $BUILDDEPS &&\
    apt-get autoremove -y && \
    apt-get autoclean -y && \
    rm -rf /var/lib/apt/lists/*

CMD ["bowtie"]
