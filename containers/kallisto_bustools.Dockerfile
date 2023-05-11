# /whitelist/
#   10xv2_whitelist.txt
#   10xv3_whitelist.txt


FROM debian:buster

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


# Install Kallisto

RUN apt-get install -y --no-install-recommends \
    openssl \
    curl \
    libcurl4-openssl-dev

RUN cd /usr/local/src && \
    curl -O https://codeload.github.com/pachterlab/kallisto/tar.gz/v0.46.1 && \
    tar xvfz v0.46.1 && \
    rm -f v0.46.1

RUN apt-get install -y --no-install-recommends \
    file \
    less \
    g++ \
    make \
    cmake \
    autoconf \
    libhdf5-dev \
    zlib1g-dev \
    libhts-dev

RUN cd /usr/local/src/kallisto-0.46.1 && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    make install



# Install BUStools

RUN apt-get install -y --no-install-recommends \
    git

RUN cd /usr/local/src && \
  git clone https://github.com/BUStools/bustools.git && \
  cd bustools && \
  mkdir build && \
  cd build && \
  cmake .. && \
  make && \
  make install



# Install scanpy.

RUN apt-get install -y --no-install-recommends \
    pkg-config \
    libhdf5-dev \
    libfreetype6-dev \
    libpng-dev \
    libxml2-dev \
    python3 \
    python3-dev \
    python3-setuptools \
    python3-pip

RUN ln -s /usr/bin/python3 /usr/bin/python

# Somehow the apt-get libraries (numpy et al.) aren't compatible with
# the scanpy from pip.  The pip scanpy seems to require later versions
# of these libraries.  Install them all from pip.
RUN pip3 install \
    numpy \
    scipy \
    matplotlib \
    pandas \
    python-igraph \
    louvain \
    scanpy



# Install t2g.py.

RUN cd /usr/local/bin && \
    curl -L -O https://github.com/BUStools/getting_started/releases/download/getting_started/t2g.py && \
    chmod +x /usr/local/bin/t2g.py

# t2g.py --use-version is broken, as documented at:
# https://groups.google.com/forum/?utm_medium=email&utm_source=footer#!msg/kallisto-and-applications/v3oejhRUcVo/39QVurrbAAAJ
# But the fix is incorrect as well.
# Need to change:
#               if 'transcript_version' not in d or 'gene_version' not in d:
#                    continue
#               # tid += '.' + d['transcript_version']
#               # gid += '.' + d['gene_version']
# to:
#               if False:
#                    continue
#               tid = d['transcript_id']
#               gid = d['gene_id']

#cp /usr/local/bin/t2g.py /usr/local/bin/t2g.py.old
#cp /usr/local/bin/t2g.py.old /usr/local/bin/t2g.py
#diff /usr/local/bin/t2g.py.old /usr/local/bin/t2g.py

RUN perl -i -p -e \
        "s/if 'transcript_version' not in d or 'gene_version' .*$/if False:/" \
        /usr/local/bin/t2g.py && \
    perl -i -p -e "s/# tid \+= '\.' .*$/tid = d['transcript_id']/" \
        /usr/local/bin/t2g.py && \
    perl -i -p -e "s/# gid \+= '\.' .*$/gid = d['gene_id']/" \
        /usr/local/bin/t2g.py


# Install whitelist, etc.

RUN mkdir /whitelist
RUN cd /whitelist && \
    curl -L -O https://github.com/BUStools/getting_started/releases/download/getting_started/10xv2_whitelist.txt && \
    curl -L -O https://github.com/BUStools/getting_started/releases/download/species_mixing/10xv3_whitelist.txt




  




## Clean up apt-get.
#RUN cd / && \
#    apt-get remove --purge -y $BUILDDEPS &&\
#    apt-get autoremove -y && \
#    apt-get autoclean -y && \
#    rm -rf /var/lib/apt/lists/*
