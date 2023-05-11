# Heavily modified from rocker-versioned Dockerfile.

# TODO: Make sure R doesn't try to load any of the user's libraries.
# Need to make sure .libPaths is:
# .libPaths(c("/usr/local/lib/R/site-library", "/usr/local/lib/R/library"))

# /usr/local/changlab/Rlib/
# /usr/local/ExperimentHub/     Cache for SingleR files

FROM debian:bullseye

LABEL maintainer="Jeffrey Chang <jeffrey.t.chang@uth.tmc.edu>"

ENV R_VERSION=4.0.3 \
    BUILD_DATE="$(TZ='America/Los_Angeles' date -I)" \
    LC_ALL=en_US.UTF-8 \
    LANG=en_US.UTF-8 \
    TERM=xterm

RUN apt-get update
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
    vim \
    openssl \
    libssl-dev \
    libcurl4-openssl-dev \
    curl


# Developer stuff.
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gcc \
    gfortran \
    make \
    git


# Install R.

RUN apt-get update && apt-get install -y --no-install-recommends \
    libreadline-dev \
    libblas-dev \
    liblapack-dev \
    libicu-dev

RUN apt-get update && apt-get install -y --no-install-recommends \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libpcre2-dev


# Download, configure, and make R.
RUN cd /usr/local/src && \
    curl -O https://cran.r-project.org/src/base/R-4/R-${R_VERSION}.tar.gz  && \
    tar xfz R-${R_VERSION}.tar.gz && \
    rm -f R-${R_VERSION}.tar.gz && \
    cd /usr/local/src/R-${R_VERSION} && \
    R_PAPERSIZE=letter \
    R_BATCHSAVE="--no-save --no-restore" \
    R_BROWSER=xdg-open \
    PAGER=/usr/bin/pager \
    PERL=/usr/bin/perl \
    R_UNZIPCMD=/usr/bin/unzip \
    R_ZIPCMD=/usr/bin/zip \
    R_PRINTCMD=/usr/bin/lpr \
    LIBnn=lib \
    AWK=/usr/bin/awk \
    CFLAGS="-g -O2 -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g" \
    CXXFLAGS="-g -O2 -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g" \
    ./configure --enable-R-shlib \
                --enable-memory-profiling \
                --with-readline \
                --with-blas \
                --with-tcltk \
                --disable-nls \
                --with-x=no \
                --with-recommended-packages && \
    make && \
    make install && \
    rm -rf /usr/local/src/R-${R_VERSION}



# Add a default CRAN mirror
RUN mkdir -p /usr/local/lib/R/etc && \
    echo "options(repos=c(CRAN='https://cran.rstudio.com/'), download.file.method='libcurl')" \
    >> /usr/local/lib/R/etc/Rprofile.site

# Add a library directory (for user-installed packages)
RUN mkdir -p /usr/local/lib/R/site-library && \
    chown root:staff /usr/local/lib/R/site-library && \
    chmod g+wx /usr/local/lib/R/site-library

# Fix library path
# spython messes up when parsing the string "{U}SER".  Will interpret as
# the {U}SER variable.  So need to munge the string.
ENV A=USE
RUN echo "R_LIBS_${A}R='/usr/local/lib/R/site-library'" >> /usr/local/lib/R/etc/Renviron && \
  echo "R_LIBS=\${R_LIBS-'/usr/local/lib/R/site-library:/usr/local/lib/R/library:/usr/lib/R/library'}" >> /usr/local/lib/R/etc/Renviron

# Use littler installation scripts
# Installation won't work if I'm not in a writeable directory.
RUN Rscript -e "install.packages(c('littler', 'docopt'))" && \
  ln -s /usr/local/lib/R/site-library/littler/examples/install2.r \
    /usr/local/bin/install2.r && \
  ln -s /usr/local/lib/R/site-library/littler/examples/installGithub.r \
    /usr/local/bin/installGithub.r && \
  ln -s /usr/local/lib/R/site-library/littler/bin/r /usr/local/bin/r


# Install python 2.
# Needed for Seurat.
RUN apt-get install -y --no-install-recommends \
  python-dev python-setuptools
RUN cd /usr/local/src && \
  curl -O https://bootstrap.pypa.io/pip/2.7/get-pip.py && \
  python2 get-pip.py
RUN pip install wheel
RUN apt-get install -y --no-install-recommends \
  libhdf5-dev 
RUN pip install h5py


# Seurat
# umap-learn requires llvmlite
# Version >= 0.32.0 requires python3
RUN pip install llvmlite==0.31.0
# umap-learn >= 0.4.0 requires python3.
RUN pip install umap-learn==0.3.9
RUN pip install virtualenv   # needed for reticulate
RUN install2.r -e hdf5r
RUN apt-get update && apt-get install -y --no-install-recommends \
  libpng-dev 
RUN install2.r -e Seurat


# Install Bioconductor
RUN Rscript -e 'if(!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager"); BiocManager::install()'
RUN apt-get install -y --no-install-recommends \
  libxml2-dev
RUN Rscript -e "BiocManager::install('BiocParallel')"
RUN Rscript -e "BiocManager::install('XVector')"
RUN Rscript -e "BiocManager::install('GenomicFeatures')"
RUN Rscript -e "BiocManager::install('ensembldb')"

# SingleR
# Should install SingleR after Seurat.
# Need to install tibble and dbplyr manually before scRNAseq or there
# will be errors.
RUN install2.r -e tibble
RUN install2.r -e dbplyr
# Install scRNAseq package for training data sets.
RUN Rscript -e "BiocManager::install('scRNAseq')"
# Github version is obsolete.  Now on Bioconductor.
#RUN Rscript -e "devtools::install_github('dviraran/SingleR')"
# Devel version requires R 4.
#RUN Rscript -e "BiocManager::install(version='devel'); BiocManager::install('SingleR')"
# BiocNeighbors (required by SingleR) is hard to install.
# - Use Bioconductor release 3.10 for R 3.6.
#   https://www.bioconductor.org/about/release-announcements/
# - Release 3.10 requires RcppAnnoy version 0.0.16.
#   Incompatible with 0.0.17.
RUN install2.r -e devtools
#RUN Rscript -e "devtools::install_version('RcppAnnoy', version='0.0.16')"
RUN Rscript -e "devtools::install_version('RcppAnnoy')"
RUN Rscript -e "BiocManager::install('BiocNeighbors')"
RUN Rscript -e "BiocManager::install('SingleR')"

# Cache the reference data for SingleR.
RUN Rscript -e "BiocManager::install('celldex')"
#RUN Rscript -e "library(celldex) ; \
#  HumanPrimaryCellAtlasData();"
#  library(SingleR) ; \

RUN mkdir -p /usr/local/ExperimentHub
RUN Rscript -e "library(ExperimentHub); \
  setExperimentHubOption(\"CACHE\", \"/usr/local/ExperimentHub\"); \
  library(celldex); \
  HumanPrimaryCellAtlasData(); \
  BlueprintEncodeData(); \
  DatabaseImmuneCellExpressionData(); \
  NovershternHematopoieticData(); \
  MonacoImmuneData(); \
  ImmGenData(); \
  MouseRNAseqData()"

# Not quite the same thing as in ExperimentHub.
#ENV A=dviraran/SingleR/blob/master/data
#RUN mkdir -p /usr/local/SingleR/data && \
#  cd /usr/local/SingleR/data && \
#  curl -L -o blueprint_encode.rda \
#    https://github.com/${A}/blueprint_encode.rda?raw=true
## load("/usr/local/SingleR/data/blueprint_encode.rda")


# Install genomicode.

COPY changlab /usr/local/src/changlab

# /usr/local/changlab/Rlib/
RUN mkdir -p /usr/local/changlab/
COPY changlab/Rlib /usr/local/changlab/Rlib


# Clean up apt-get.
RUN cd / && \
    apt-get remove --purge -y $BUILDDEPS &&\
    apt-get autoremove -y && \
    apt-get autoclean -y && \
    rm -rf /var/lib/apt/lists/*

CMD ["R"]
