# Heavily modified from rocker-versioned Dockerfile.

# /usr/local/genomicode/bin/
#   slice_svm.py
#   [...]
# /usr/local/genomidata/
#   psid2platform/
#   convert_platform/
# /usr/local/bin/
#   cluster            Cluster 3.0
# /usr/local/lib/python2.7/dist-packages/genomicode/
#   bin/
#     [... scripts]
#   Rlib/
#   Verdana.ttf
#   Verdana Bold.ttf
#   MS PGothic.ttf
#   [... python libraries]
# /usr/local/src/changlab/Betsy/scripts/
#   summarize_alignments_at_close_variants.py
#
# Deprecate this use:
# /usr/local/src/changlab/scripts/
#   slice_svm.py

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


# Miscellaneous useful things.
# Don't know why, but apt-get update here prevents 404.
RUN apt-get update && apt-get install -y --no-install-recommends \
    file \
    less \
    vim \
    openssl \
    libssl-dev \
    libcurl4-openssl-dev \
    curl


# Developer stuff.
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gcc \
    make \
    git

# Python
RUN apt-get update && apt-get install -y --no-install-recommends \
    python2.7 \
    python2.7-dev \
    python-pip


# Genomicode requirements.
RUN pip install setuptools wheel
RUN pip install \
  xlwt openpyxl xlrd \
  psutil

# Install libraqm for Pillow.
RUN apt-get update && apt-get install -y --no-install-recommends \
    libfribidi-dev \
    libfreetype6-dev \
    libpangocairo-1.0-0 \
    libpango1.0-dev \
    gtk-doc-tools

RUN cd /usr/local/src && \
    curl -L -O https://github.com/harfbuzz/harfbuzz/releases/download/2.8.0/harfbuzz-2.8.0.tar.xz && \
    xzcat harfbuzz-2.8.0.tar.xz | tar xvf - && \
    cd harfbuzz-2.8.0 && \
    ./configure && \
    make && \
    make install

RUN cd /usr/local/src && \
    curl -L -O https://github.com/HOST-Oman/libraqm/releases/download/v0.7.1/raqm-0.7.1.tar.gz && \
    tar xvfz raqm-0.7.1.tar.gz && \
    cd raqm-0.7.1 && \
    ./configure && \
    make && \
    make install

# Pillow 6.2.2 is broken and gives weird errors (e.g. IOError: invalid
# face handle).  Use an earlier version.
# 6.2.2 broken
# 5.4.1 broken
# 5.0.0 broken
# 4.3.0 works
# 3.4.2 works
RUN pip install Pillow==4.3.0

RUN pip install \
  numpy pandas matplotlib scipy
RUN pip install h5py


RUN apt-get update && apt-get install -y --no-install-recommends \
    graphviz-dev
RUN pip install pygraphviz

# PySam is a bit complicated.
RUN apt-get update && apt-get install -y --no-install-recommends \
    libbz2-dev    \
    liblzma-dev
# Cython needs to be installed before pysam.
RUN pip install Cython && \
  pip install pysam

# arial10
RUN cd /usr/local/lib/python2.7/dist-packages && \
  curl -O https://raw.githubusercontent.com/juanpex/django-model-report/master/model_report/arial10.py


# More requirements

# Install miscellaneous tools needed in programs.
# Need ghostscript so imagemagick can read PDF files.
# GDBM needed for geneidlib to read psid2platform.
RUN apt-get install -y --no-install-recommends \
    libgs-dev \
    ghostscript \
    python-gdbm \
    pigz \
    wget


# Install Cluster 3.0.
RUN cd /usr/local/src && \
  curl -O http://bonsai.hgc.jp/~mdehoon/software/cluster/cluster-1.59.tar.gz && \
  tar xvfz cluster-1.59.tar.gz && \
  cd cluster-1.59 && \
  ./configure --without-x && \
  make && \
  make install


# Install imagemagick.
# apt-get version doesn't use ghostscript and can't read PDF.
# Install from source.

ENV A=ImageMagick-6.9.9-37
#ENV A=ImageMagick-7.0.7-25
RUN cd /usr/local/src && \
    wget --default-page=${A}.tar.gz \
     https://sourceforge.net/projects/imagemagick/files/im6-src/${A}.tar.gz/&&\
    tar zxvf ${A}.tar.gz && \
    rm ${A}.tar.gz
RUN cd /usr/local/src/${A} && \
    ./configure --prefix=/usr --with-gslib=yes && \
    make && \
    make install && \
    rm -rf /usr/local/src/${A}


# Not sure this is needed?
#POV-Ray
# sudo apt-get install libboost-dev zlib1g-dev libpng12-dev libjpeg8-dev libtiff5-dev libopenexr-dev
#wget http://www.povray.org/redirect/www.povray.org/ftp/pub/povray/Official/Unix/povray-3.6.tar.gz
#tar -zxvf povray-3.6.tar.gz
#cd povray-3.6.1
#./configure COMPILED_BY="Jeffrey Chang <jeffrey.t.chang@uth.tmc.edu>"
#make
#sudo make install





# samtools library requires samtools.  Just install it there.

RUN apt-get install -y --no-install-recommends \
    bzip2

ENV HTSLIB_VERSION=1.10.2
RUN cd /usr/local/src && \
    curl -L -O https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 && \
    bzcat htslib-${HTSLIB_VERSION}.tar.bz2 | tar xf - && \
    rm -f htslib-${HTSLIB_VERSION}.tar.bz2

RUN cd /usr/local/src/htslib-${HTSLIB_VERSION} && \
    ./configure && \
    make && \
    make install


RUN apt-get install -y --no-install-recommends \
    libncurses5-dev

ENV SAMTOOLS_VERSION=1.10
RUN cd /usr/local/src && \
    curl -L -O https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    bzcat samtools-${SAMTOOLS_VERSION}.tar.bz2 | tar xf - && \
    rm -f samtools-${SAMTOOLS_VERSION}.tar.bz2

RUN cd /usr/local/src/samtools-${SAMTOOLS_VERSION} && \
    ./configure && \
    make && \
    make install



# Install rpy2.

# Install R 3.3.3.
RUN apt-get install -y --no-install-recommends \
    r-base
RUN apt-get install -y --no-install-recommends \
    libreadline-dev \
    libblas-dev \
    liblapack-dev \
    libicu-dev

RUN pip install rpy2==2.4.4




# Install biomaRt library for geneidlib.  We will be able to run
# annotate_matrix.py in a container so the user doesn't need to
# install it.

# Add a default CRAN mirror
RUN mkdir -p /usr/lib/R/etc && \
    echo "options(repos=c(CRAN='https://cran.rstudio.com/'), download.file.method='libcurl')" \
    >> /usr/lib/R/etc/Rprofile.site

# Add a library directory (for user-installed packages)
RUN mkdir -p /usr/lib/R/site-library && \
    chown root:staff /usr/lib/R/site-library && \
    chmod g+wx /usr/lib/R/site-library

# Fix library path
# spython messes up when parsing the string "{U}SER".  Will interpret as
# the {U}SER variable.  So need to munge the string.
ENV A=USE
RUN echo "R_LIBS_${A}R='/usr/lib/R/site-library'" >> /usr/lib/R/etc/Renviron && \
  echo "R_LIBS=\${R_LIBS-'/usr/lib/R/site-library:/usr/lib/R/library'}" >> /usr/lib/R/etc/Renviron

RUN apt-get update && apt-get install -y --no-install-recommends \
    libxml2-dev \
    libcurl4-openssl-dev

# Installs biomaRt_2.30.0 and RCurl 1.95-0.1.2.
# BiocManager will try to install version of RCurl (1.98-1.2), which
# requires R >= 3.4.0.
# Need to install an older version of XML library first.  3.99+
# requires R 4.0.
RUN Rscript -e "install.packages('devtools')"
RUN Rscript -e \
  'require(devtools); install_version("XML", version="3.98-1.20")'
RUN Rscript -e \
  'source("https://bioconductor.org/biocLite.R"); biocLite("biomaRt")'

#RUN Rscript -e 'if(!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager"); BiocManager::install()'
#
#RUN Rscript -e "BiocManager::install(c('biomaRt'))"



## Install R library edgeR.
## Needed for slice_matrix.py --convert_counts_to_cpm
## No.  Just have it run it from the R library without rpy2.
#
#RUN Rscript -e \
#  'require(devtools); install_version("locfit", version="1.5-9.2")'
#RUN Rscript -e \
#  'source("https://bioconductor.org/biocLite.R"); biocLite("edgeR")'




RUN pip install XlsxWriter
RUN pip install python-pptx




# Install genomicode.

COPY changlab /usr/local/src/changlab

ENV USER=root
RUN cd /usr/local/src/changlab && \
  python setup.py install
# rmdir /changlab

RUN mkdir -p /usr/local/genomicode/bin && \
  cp -Rp /usr/local/lib/python2.7/dist-packages/genomicode/bin/* /usr/local/genomicode/bin
ENV PATH /usr/local/genomicode/bin:$PATH

RUN mkdir -p /usr/local/genomidata/
COPY psid2platform /usr/local/genomidata/psid2platform
COPY convert_platform /usr/local/genomidata/convert_platform



# Configure .genomicoderc?

# Have error:
# /bin/bash: /usr/bin/lesspipe.sh: No such file or directory
# This fixes things in Docker, but breaks "less" in Singularity.
#RUN ln -s /bin/lesspipe /usr/bin/lesspipe.sh
# Just ignore this.  LESSOPEN environment variable set incorrectly for
# container.


CMD ["python"]
