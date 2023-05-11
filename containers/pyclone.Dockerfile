# Heavily modified from rocker-versioned Dockerfile.

#FROM debian:stretch
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
#RUN apt-get update && apt-get install -y --no-install-recommends \
#    file \
#    less \
#    vim \
#    openssl \
#    libssl-dev \
#    libcurl4-openssl-dev \
#    curl


# Install miniconda for PyClone.

RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    bzip2 \
    libglib2.0-0 \
    libxext6 \
    libsm6 \
    libxrender1

RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    mercurial \
    subversion

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda2-4.7.12-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc && \
    find /opt/conda/ -follow -type f -name '*.a' -delete && \
    find /opt/conda/ -follow -type f -name '*.js.map' -delete && \
    /opt/conda/bin/conda clean -afy

ENV PATH /opt/conda/bin:$PATH


# Install PyClone.
# https://bitbucket.org/aroth85/pyclone/wiki/Home

RUN conda install pyclone -c aroth85

# Server is headless, and matplotlib code will generate errors without
# graphics driver.
# http://matplotlib.org/faq/
#   howto_faq.html#generate-images-without-having-a-window-appear
#
# ImportError: Failed to import any qt binding
# >>> import PyQt5
# >>> import matplotlib.pyplot
#
# Add before "import matplotlib" in:
# /opt/conda/lib/python2.7/site-packages/pyclone/post_process/plot/
#   clusters.py
#   loci.py
#   _scatter.py
#   utils.py
#
# import matplotlib          <<< add
# matplotlib.use('Agg')      <<< add

# In _scatter, needs to go after from __future__.
RUN cd /opt/conda/lib/python2.7/site-packages/pyclone/post_process/plot/ && \
  perl -pi -e 'print "import matplotlib\nmatplotlib.use(\"Agg\")\n" if $. == 1' clusters.py && \
  perl -pi -e 'print "import matplotlib\nmatplotlib.use(\"Agg\")\n" if $. == 1' loci.py && \
  perl -pi -e 'print "import matplotlib\nmatplotlib.use(\"Agg\")\n" if $. == 7' _scatter.py && \
  perl -pi -e 'print "import matplotlib\nmatplotlib.use(\"Agg\")\n" if $. == 1' utils.py


CMD ["PyClone"]
