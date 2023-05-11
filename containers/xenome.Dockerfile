# Heavily modified from rocker-versioned Dockerfile.
FROM debian:stretch

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

# Install application specific stuff.
RUN apt-get install -y --no-install-recommends \
    make \
    g++ \
    cmake \
    pandoc \
    libboost-all-dev \
    zlib1g-dev \
    libbz2-dev \
    libsqlite3-dev \
    git

RUN cd /usr/local/src && \
    git clone https://github.com/data61/gossamer.git

RUN cd /usr/local/src/gossamer && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    make test && \
    make install

#  ldd /usr/local/bin/xenome

## Clean up apt-get.
#RUN cd / && \
#    apt-get remove --purge -y $BUILDDEPS &&\
#    apt-get autoremove -y && \
#    apt-get autoclean -y && \
#    rm -rf /var/lib/apt/lists/*

ENTRYPOINT ["/usr/local/bin/xenome"]
