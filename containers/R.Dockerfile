# Heavily modified from rocker-versioned Dockerfile.

# 221019 Updated to Debian 11 (bullseye).  buster not compiling anymore.
#        Updated to R 4.2.1.
# 201218 Updated to Debian 10 (buster).  devtools (gert) not installing in
#        Debian 9 (stretch).
# 201218 Updated from R 3.6.3 to R 4.0.3.  Stuff from Bioconductor
#        was becoming incompatible.


# TODO: Make sure R doesn't try to load any of the user's libraries.
# Need to make sure .libPaths is:
# .libPaths(c("/usr/local/lib/R/site-library", "/usr/local/lib/R/library"))


# /usr/local/changlab/Rlib/
# /usr/local/genomicode/bin/
#   slice_svm.py
#   [...]
# /usr/local/genomidata/
#   affymetrix/
#   illumina/
#   psid2platform/
#   convert_platform/



#FROM debian:stretch
#FROM debian:buster
FROM debian:bullseye

LABEL maintainer="Jeffrey Chang <jeffrey.t.chang@uth.tmc.edu>"

ENV R_VERSION=4.2.1 \
    BUILD_DATE="$(TZ='America/Los_Angeles' date -I)" \
    LC_ALL=en_US.UTF-8 \
    LANG=en_US.UTF-8 \
    TERM=xterm

# Fail now if GenePattern can't be found.  Will install later.
COPY GenePattern_1.0.2.tar.gz /usr/local/src


RUN apt-get update

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    locales \
    libreadline8 \
    bash-completion \
    file

RUN apt-get update && apt-get install -y --no-install-recommends \
    make \
    gfortran \
    g++

RUN apt-get update && apt-get install -y --no-install-recommends \
    libpcre2-dev \
    libpcre3

RUN apt-get update && apt-get install -y --no-install-recommends \
    libbz2-1.0 \
    libcurl4 \
    #libblas-dev \
    #libicu57 \
    liblzma5 \
    unzip \
    zip \
    zlib1g

RUN apt-get update && apt-get install -y --no-install-recommends \
    fonts-texgyre \
    gsfonts \
    libjpeg62-turbo \
    libpangocairo-1.0-0 \
    libpng16-16 \
    libtiff5

# OpenBLAS >= 0.3.4 generates problems with R affy library:
# ERROR; return code from pthread_create() is 22
# Use BLAS instead.
# Check BLAS version:
# R> SessionInfo()
RUN apt-get install -y --no-install-recommends \
    libblas-dev \
    liblapack-dev
#    libopenblas-dev \


RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen && \
    /usr/sbin/locale-gen en_US.UTF-8 && \
    /usr/sbin/update-locale LANG=en_US.UTF-8

RUN apt-get install -y --no-install-recommends \
    curl \
    default-jdk \
    libbz2-dev \
    libcairo2-dev \
    libcurl4-openssl-dev \
    libpango1.0-dev \
    libjpeg-dev \
    libicu-dev \
    libpcre3-dev \
    libpng-dev \
    libreadline-dev \
    libtiff5-dev \
    liblzma-dev \
    libx11-dev \
    libxt-dev \
    perl \
    tcl8.6-dev \
    tk8.6-dev \
    texinfo \
    texlive-extra-utils \
    texlive-fonts-recommended \
    texlive-fonts-extra \
    texlive-latex-recommended \
    x11proto-core-dev \
    xauth \
    xfonts-base \
    xvfb \
    zlib1g-dev


# Install Ghostscript for plotting PDF.
RUN apt-get install -y --no-install-recommends \
    ghostscript libgs-dev


# Install other packages that are convenient to work with.
RUN apt-get install -y --no-install-recommends \
    libssl-dev wget \
    libxml2-dev libatlas3-base \
    vim less



# Download R.
RUN cd /usr/local/src && \
    curl -O https://cran.r-project.org/src/base/R-4/R-${R_VERSION}.tar.gz  && \
    tar xfz R-${R_VERSION}.tar.gz && \
    rm -f R-${R_VERSION}.tar.gz




# Configure and make
RUN cd /usr/local/src/R-${R_VERSION} && \
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
                --with-recommended-packages && \
    make && \
    make install

# Doesn't work.  Cannot delete because directory not empty.
# rm -rf /usr/local/src/R-${R_VERSION}


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
# Needed for scripts (find_diffexp_genes.py).
# Needed for Seurat.
RUN apt-get update && apt-get install -y --no-install-recommends \
  python-dev python-gdbm python-setuptools
# python-pip no longer supported.  Install it manually.
RUN cd /usr/local/src && \
  curl https://bootstrap.pypa.io/pip/2.7/get-pip.py -o get-pip.py && \
  python2.7 get-pip.py

RUN pip install wheel

RUN apt-get install -y --no-install-recommends \
  libhdf5-dev
RUN pip install h5py

# Actually, not needed anymore.
## Install rpy2.
## Needed for find_diffexp_genes.py.
#RUN pip install rpy2==2.4.4 && \
#  echo "/usr/local/lib/R/lib" >> /etc/ld.so.conf.d/libR.conf && \
#  ldconfig

## Install rpy2 manually since we need to change some setup parameters.
#RUN cd /usr/local/src && \
#  pip download rpy2==2.4.4 && \
#  tar xvfz rpy2-2.4.4.tar.gz && \
#  cd rpy2-2.4.4 && \
#  # Add library path after line:  \
#  # extra_link_args = []          \
#  sed -i '/extra_link_args = \[\]/ a \    extra_link_args.append("-Wl,-rpath,/usr/local/lib/R/lib")' setup.py && \
#  sed -i '/extra_link_args = \[\]/ a \    extra_link_args.append("-L/usr/local/lib/R/lib")' setup.py && \
#  python setup.py build && \
#  python setup.py install




# Install R libraries.
RUN apt-get install -y --no-install-recommends \
    libgit2-dev
RUN install2.r -e curl rstudioapi openssl
RUN install2.r -e gert
RUN install2.r -e usethis
RUN install2.r -e R.utils devtools

# GenePattern.
RUN install2.r -e rJava
RUN cd /usr/local/src && \
  R CMD INSTALL GenePattern_1.0.2.tar.gz && \
  rm /usr/local/src/GenePattern_1.0.2.tar.gz

# Need to install an older version of XML library first.  3.99+
# requires R 4.0.
RUN Rscript -e \
  'require(devtools); install_version("XML", version="3.98-1.20")'
# Install Cairo to prevent message:
# 1: In plot.xy(xy.coords(x, y), type = type, ...) :
# semi-transparency is not supported on this device: reported only
# once per page
RUN install2.r -e Cairo

# ERROR: dependency ‘readr’ is not available for package ‘haven’
# Gigantion runs out of /tmp disk space if I try to install everything
# at once.  Break it up into smaller installs.  Install some basic
# packages here so they won't need to be reinstalled later.
RUN install2.r -e readr
RUN install2.r -e shiny
RUN install2.r -e ggplot2
RUN install2.r -e plotly
RUN install2.r -e tidyr
RUN install2.r -e plyr
RUN install2.r -e base64
RUN install2.r -e reshape2
RUN install2.r -e snow
RUN install2.r -e iterators
RUN install2.r -e formatR

RUN install2.r -e WriteXLS MASS
RUN install2.r -e beeswarm venneuler
RUN install2.r -e SparseM mvtnorm pwr 
RUN install2.r -e Matrix MatrixModels RcppEigen
# CMake needed for nloptr
RUN apt-get update && apt-get install -y --no-install-recommends \
    cmake
RUN install2.r -e nloptr
RUN install2.r -e car
RUN install2.r -e mclust cluster sigclust dbscan
RUN install2.r -e e1071 LiblineaR randomForest naivebayes
RUN install2.r -e GSA

# Install Bioconductor
RUN Rscript -e 'if(!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager"); BiocManager::install()'
RUN Rscript -e "BiocManager::install('BiocParallel')"
RUN Rscript -e "BiocManager::install('XVector')"
RUN Rscript -e "BiocManager::install('GenomicFeatures')"
RUN Rscript -e "BiocManager::install('Rsamtools')"
RUN Rscript -e "BiocManager::install('AnnotationForge')"

RUN Rscript -e "BiocManager::install(c('biomaRt'))"
RUN Rscript -e "BiocManager::install(c('affy', 'marray', 'annotate'))"
# Latest version doesn't work.  Breaks oligoClasses.  Install an older
# version.  Needs to be installed right before oligoClasses.
# Otherwise, something else will update it.
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/Archive/ff/ff_2.2-14.tar.gz', repos=NULL)"
RUN Rscript -e "BiocManager::install('oligoClasses')"
RUN Rscript -e "BiocManager::install('oligo')"
RUN Rscript -e "BiocManager::install(c('lumi', 'methylumi'))"

# Not available for R 3.6.0.
#RUN Rscript -e "BiocManager::install(c('IlluminaHumanMethylation450k.db'))"
RUN Rscript -e "BiocManager::install(c('IlluminaHumanMethylation450kprobe'))"
RUN Rscript -e "BiocManager::install(c('IlluminaHumanMethylation27k.db'))"

RUN Rscript -e "BiocManager::install(c('makecdfenv'))"
RUN Rscript -e "BiocManager::install(c('hgu95av2', 'hgu95av2cdf'))"
RUN Rscript -e "BiocManager::install(c('hgu133a2cdf', 'hgu133plus2cdf', 'hgu133acdf'))"
RUN Rscript -e "BiocManager::install(c('hthgu133acdf', 'hthgu133bcdf', 'u133aaofav2cdf'))"
# Not available for R 3.6.0.
#RUN Rscript -e "BiocManager::install(c('hta20cdf'))"
# 356 Mb.  Not necessary.
#RUN Rscript -e "BiocManager::install(c('pd.hta.2.0'))"

RUN Rscript -e "BiocManager::install(c('limma', 'edgeR', 'DESeq2'))"
RUN Rscript -e "BiocManager::install(c('impute'))"
RUN Rscript -e "BiocManager::install(c('GSVA'))"
RUN Rscript -e "BiocManager::install(c('copynumber'))"
RUN Rscript -e "BiocManager::install(c('LPE'))"
RUN Rscript -e "BiocManager::install(c('bioDist'))"
# Need to install samr after impute installed.
RUN install2.r -e samr

# minfi library for Illumina methylation microarrays.
RUN Rscript -e "BiocManager::install('minfi')"

RUN Rscript -e "BiocManager::install(c('org.Hs.eg.db', 'org.Mmu.eg.db'))"
RUN Rscript -e "BiocManager::install('clusterProfiler')"




# Seurat
# umap-learn requires llvmlite
# Version >= 0.32.0 requires python3
RUN pip install llvmlite==0.31.0
# umap-learn >= 0.4.0 requires python3.
RUN pip install umap-learn==0.3.9
RUN pip install virtualenv   # needed for reticulate
RUN apt-get update && apt-get install -y --no-install-recommends \
    libgeos-dev
RUN install2.r -e rgeos   # required for SeuratObject
RUN install2.r -e hdf5r Seurat

# umap package
RUN install2.r -e umap


# maftools
RUN Rscript -e "BiocManager::install('maftools')"


# RcppAnnoy
# Need this to load Jinfeng's RDS files.
RUN Rscript -e "install.packages('RcppAnnoy')"


# VennDiagram
RUN Rscript -e "devtools::install_version('VennDiagram', version='1.6.9')"



# No.  Out of date, and can't download manifest files anymore.
## lumidat
## Get warning:
## Skipping 3 packages not available: lumi, limma, Biobase
## But they should be installed already?
#RUN Rscript -e "devtools::install_github('drmjc/lumidat')"


# Satija Signac library.
RUN Rscript -e "install.packages('Signac')"


# Install the qs library.
RUN Rscript -e "install.packages('qs')"


# Install glmGamPoi for Seurat sctransform.
RUN Rscript -e "BiocManager::install('glmGamPoi')"



# Install genomicode.

COPY changlab /usr/local/src/changlab

ENV USER=root
RUN cd /usr/local/src/changlab && \
  python setup.py install
# rmdir /changlab

# /usr/local/changlab/Rlib/
RUN mkdir -p /usr/local/changlab/
COPY changlab/Rlib /usr/local/changlab/Rlib

# /usr/local/genomicode/bin/
RUN mkdir -p /usr/local/genomicode/bin && \
  cp -Rp /usr/local/lib/python2.7/dist-packages/genomicode/bin/* /usr/local/genomicode/bin
ENV PATH /usr/local/genomicode/bin:$PATH

# /usr/local/genomidata/
RUN mkdir -p /usr/local/genomidata/
COPY psid2platform /usr/local/genomidata/psid2platform
COPY convert_platform /usr/local/genomidata/convert_platform
COPY affymetrix /usr/local/genomidata/affymetrix
COPY illumina /usr/local/genomidata/illumina


# Clean up apt-get.
RUN cd / && \
    apt-get remove --purge -y $BUILDDEPS &&\
    apt-get autoremove -y && \
    apt-get autoclean -y && \
    rm -rf /var/lib/apt/lists/*

CMD ["R"]
