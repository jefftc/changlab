# Heavily modified from rocker-versioned Dockerfile.

# /usr/local/bin/povray
# /usr/local/etc/povray/3.6/
#   povray.conf
# /usr/local/genomicode/bin/
#   pcaplot.py
#   pybinreg.py
#   scoresig.py
#   [...]


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
    vim

RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    curl


# Developer stuff.

#RUN apt-get update && apt-get install -y --no-install-recommends \
#    build-essential \
#    gcc \
#    make
#    bzip2 \
#    libbz2-dev \
#    libncurses5-dev
#    liblzma-dev \

# Prerequisites

#RUN apt-get update && apt-get install -y --no-install-recommends \
#  libboost-dev \
#  zlib1g-dev \
#  libpng-dev \
#  libjpeg62-turbo-dev \
#  libtiff5-dev \
#  libopenexr-dev


# Install povray.

RUN cd /usr/local/src && \
  curl -L -O http://www.povray.org/redirect/www.povray.org/ftp/pub/povray/Old-Versions/Official-3.62/Linux/povlinux-3.6.tgz && \
  tar xvfz povlinux-3.6.tgz && \
  cd povray-3.6 && \
  ./install -no-arch-check


# Configure security.  Set to none, otherwise may get errors like:
# Scene File Parser Initialization Error: Reading from
#  '/temp00002/tmpjAb_oq.pov' is not permitted.  Check the configuration
#  in '/usr/local/etc/povray/3.6/povray.conf'.
#
# [File I/O Security]
# ;none       ; all read and write operations on files are allowed.
# ;read-only  ; uses the "read+write" directories for writing (see below).
# restricted  ; uses _only_ "read" and "read+write" directories for file I/O.

RUN perl -i -p -e 's/^restricted/none/' /usr/local/etc/povray/3.6/povray.conf




# genomicode requirements

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gcc \
    make

RUN apt-get update && apt-get install -y --no-install-recommends \
    python2.7 \
    python2.7-dev \
    python-pip
RUN pip install numpy
RUN pip install h5py

# genomicode scripts


COPY changlab /usr/local/src/changlab

ENV USER=root
RUN cd /usr/local/src/changlab && \
  python setup.py install
# rmdir /changlab

RUN mkdir -p /usr/local/genomicode/bin && \
  cp -Rp /usr/local/lib/python2.7/dist-packages/genomicode/bin/* /usr/local/genomicode/bin
ENV PATH /usr/local/genomicode/bin:$PATH



# Clean up apt-get.
RUN cd / && \
    apt-get remove --purge -y $BUILDDEPS &&\
    apt-get autoremove -y && \
    apt-get autoclean -y && \
    rm -rf /var/lib/apt/lists/*


CMD ["povray"]

