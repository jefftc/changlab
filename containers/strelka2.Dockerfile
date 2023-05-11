# Heavily modified from rocker-versioned Dockerfile.

# /usr/local/strelka/bin/
#   configureStrelkaGermlineWorkflow.py
#   configureStrelkaSomaticWorkflow.py
# /usr/local/manta/bin/
#   configManta.py

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
    libcurl4-openssl-dev \
    curl


# Python
RUN apt-get update && apt-get install -y --no-install-recommends \
    python2.7-dev
RUN ln -s /usr/bin/python2.7 /usr/bin/python2 && \
    ln -s /usr/bin/python2.7 /usr/bin/python

# Compression
RUN apt-get update && apt-get install -y --no-install-recommends \
    zlib1g-dev \
    bzip2

# Other useful things.
RUN apt-get update && apt-get install -y --no-install-recommends \
    rsync \
    less

# Install Strelka2.

ENV V=2.9.10
RUN cd /usr/local/src && \
    curl -L -O https://github.com/Illumina/strelka/releases/download/v${V}/strelka-${V}.centos6_x86_64.tar.bz2 && \
    bzcat strelka-${V}.centos6_x86_64.tar.bz2 | tar xf -
#RUN rsync -aP /usr/local/src/strelka-${V}.centos6_x86_64/ /usr/local/
RUN ln -s /usr/local/src/strelka-${V}.centos6_x86_64 /usr/local/strelka
ENV PATH $PATH:/usr/local/strelka/bin/


# Install Manta

ENV V=1.6.0
RUN cd /usr/local/src && \
    curl -L -O https://github.com/Illumina/manta/releases/download/v${V}/manta-${V}.centos6_x86_64.tar.bz2 && \
    bzcat manta-${V}.centos6_x86_64.tar.bz2 | tar xf -
#RUN rsync -aP /usr/local/src/manta-${V}.centos6_x86_64/ /usr/local/
RUN ln -s /usr/local/src/manta-${V}.centos6_x86_64 /usr/local/manta

ENV PATH $PATH:/usr/local/manta/bin/





# Clean up apt-get.
RUN cd / && \
    apt-get remove --purge -y $BUILDDEPS &&\
    apt-get autoremove -y && \
    apt-get autoclean -y && \
    rm -rf /var/lib/apt/lists/*

#CMD ["configureStrelkaSomaticWorkflow.py"]
