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


# Miscellaneous useful things.
# Don't know why, but apt-get update here prevents 404.
RUN apt-get update && apt-get install -y --no-install-recommends \
    file \
    less \
    openssl \
    libssl-dev \
    curl \
    libcurl4-openssl-dev



RUN apt-get install -y --no-install-recommends \
  python-dev python-pip python-setuptools
RUN pip install wheel

# Development tools.
RUN apt-get update && apt-get install -y --no-install-recommends \
    g++ \
    make \
    cmake \
    autoconf 

# Version >= 0.32.0 requires python3
RUN pip install llvmlite==0.31.0

RUN pip install scrublet

RUN apt-get install -y --no-install-recommends \
  python-tk



# Clean up apt-get.
RUN cd / && \
    apt-get remove --purge -y $BUILDDEPS &&\
    apt-get autoremove -y && \
    apt-get autoclean -y && \
    rm -rf /var/lib/apt/lists/*


CMD ["python"]

