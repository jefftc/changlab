# Heavily modified from rocker-versioned Dockerfile.

#FROM debian:stretch
FROM debian:buster

LABEL maintainer="Jeffrey Chang <jeffrey.t.chang@uth.tmc.edu>"

ENV R_VERSION=4.0.3 \
    BUILD_DATE="$(TZ='America/Los_Angeles' date -I)" \
    LC_ALL=en_US.UTF-8 \
    LANG=en_US.UTF-8 \
    TERM=xterm


RUN apt-get update

RUN apt-get install -y --no-install-recommends \
    bash-completion \
    ca-certificates \
    file \
    fonts-texgyre \
    g++ \
    gfortran \
    gsfonts \
    libblas-dev \
    libbz2-1.0 \
    libcurl4 \
    libjpeg62-turbo \
    libopenblas-dev \
    libpangocairo-1.0-0 \
    libpcre2-dev \
    libpcre3 \
    libpng16-16 \
    libreadline7 \
    libtiff5 \
    liblzma5 \
    locales \
    make \
    unzip \
    zip \
    zlib1g

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




## # Install python 2.
## # Needed for scripts (find_diffexp_genes.py).
## # Needed for Seurat.
## RUN apt-get install -y --no-install-recommends \
##   python-dev python-gdbm python-pip python-setuptools
## RUN pip install wheel
## 
## RUN apt-get install -y --no-install-recommends \
##   libhdf5-dev 
## RUN pip install h5py
## 
## # Actually, not needed anymore.
## ## Install rpy2.
## ## Needed for find_diffexp_genes.py.
## #RUN pip install rpy2==2.4.4 && \
## #  echo "/usr/local/lib/R/lib" >> /etc/ld.so.conf.d/libR.conf && \
## #  ldconfig
## 
## ## Install rpy2 manually since we need to change some setup parameters.
## #RUN cd /usr/local/src && \
## #  pip download rpy2==2.4.4 && \
## #  tar xvfz rpy2-2.4.4.tar.gz && \
## #  cd rpy2-2.4.4 && \
## #  # Add library path after line:  \
## #  # extra_link_args = []          \
## #  sed -i '/extra_link_args = \[\]/ a \    extra_link_args.append("-Wl,-rpath,/usr/local/lib/R/lib")' setup.py && \
## #  sed -i '/extra_link_args = \[\]/ a \    extra_link_args.append("-L/usr/local/lib/R/lib")' setup.py && \
## #  python setup.py build && \
## #  python setup.py install
## 
## 
## 
## 


# Install R libraries.
RUN apt-get install -y --no-install-recommends \
    libgit2-dev
RUN install2.r -e curl rstudioapi openssl
RUN install2.r -e gert
RUN install2.r -e usethis
RUN install2.r -e R.utils devtools


## 
## # GenePattern.
## RUN install2.r -e rJava
## RUN cd /usr/local/src && \
##   R CMD INSTALL GenePattern_1.0.2.tar.gz && \
##   rm /usr/local/src/GenePattern_1.0.2.tar.gz
## 
## # Need to install an older version of XML library first.  3.99+
## # requires R 4.0.
## RUN Rscript -e \
##   'require(devtools); install_version("XML", version="3.98-1.20")'
## # Install Cairo to prevent message:
## # 1: In plot.xy(xy.coords(x, y), type = type, ...) :
## # semi-transparency is not supported on this device: reported only
## # once per page
## RUN install2.r -e Cairo
## 
## # ERROR: dependency ‘readr’ is not available for package ‘haven’
## # Gigantion runs out of /tmp disk space if I try to install everything
## # at once.  Break it up into smaller installs.  Install some basic
## # packages here so they won't need to be reinstalled later.
## RUN install2.r -e readr
## RUN install2.r -e shiny
## RUN install2.r -e ggplot2
## RUN install2.r -e plotly
## RUN install2.r -e tidyr
## RUN install2.r -e plyr
## RUN install2.r -e base64
## RUN install2.r -e reshape2
## RUN install2.r -e snow
## RUN install2.r -e iterators
## RUN install2.r -e formatR
## 
## RUN install2.r -e WriteXLS MASS
## RUN install2.r -e beeswarm venneuler
## RUN install2.r -e SparseM mvtnorm pwr car
## RUN install2.r -e mclust cluster sigclust dbscan
## RUN install2.r -e e1071 LiblineaR randomForest naivebayes
## RUN install2.r -e GSA
## 


# Install Bioconductor
RUN Rscript -e 'if(!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager"); BiocManager::install()'

RUN Rscript -e "BiocManager::install('BiocGenerics')"
RUN Rscript -e "BiocManager::install('DelayedArray')"
RUN Rscript -e "BiocManager::install('DelayedMatrixStats')"
RUN Rscript -e "BiocManager::install('limma')"
RUN Rscript -e "BiocManager::install('S4Vectors')"
RUN Rscript -e "BiocManager::install('SingleCellExperiment')"
RUN Rscript -e "BiocManager::install('SummarizedExperiment')"
RUN Rscript -e "BiocManager::install('batchelor')"
RUN Rscript -e "BiocManager::install('Matrix.utils')"


RUN Rscript -e "install.packages('devtools')"
RUN Rscript -e "devtools::install_github('cole-trapnell-lab/leidenbase')"
RUN apt-get install -y --no-install-recommends \
    libudunits2-dev
RUN Rscript -e "install.packages('units')"

RUN apt-get install -y --no-install-recommends \
    libgdal-dev 
RUN install2.r -e sf

RUN Rscript -e "devtools::install_github('cole-trapnell-lab/monocle3')"



## RUN Rscript -e "BiocManager::install(c('biomaRt'))"
## RUN Rscript -e "BiocManager::install(c('affy', 'marray', 'annotate'))"
## # Latest version doesn't work.  Breaks oligoClasses.  Install an older
## # version.  Needs to be installed right before oligoClasses.
## # Otherwise, something else will update it.
## RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/Archive/ff/ff_2.2-14.tar.gz', repos=NULL)"
## RUN Rscript -e "BiocManager::install('oligoClasses')"
## RUN Rscript -e "BiocManager::install('oligo')"
## RUN Rscript -e "BiocManager::install(c('lumi', 'methylumi'))"
## 
## # Not available for R 3.6.0.
## #RUN Rscript -e "BiocManager::install(c('IlluminaHumanMethylation450k.db'))"
## RUN Rscript -e "BiocManager::install(c('IlluminaHumanMethylation450kprobe'))"
## RUN Rscript -e "BiocManager::install(c('IlluminaHumanMethylation27k.db'))"
## 
## RUN Rscript -e "BiocManager::install(c('makecdfenv'))"
## RUN Rscript -e "BiocManager::install(c('hgu95av2', 'hgu95av2cdf'))"
## RUN Rscript -e "BiocManager::install(c('hgu133a2cdf', 'hgu133plus2cdf', 'hgu133acdf'))"
## RUN Rscript -e "BiocManager::install(c('hthgu133acdf', 'hthgu133bcdf', 'u133aaofav2cdf'))"
## # Not available for R 3.6.0.
## #RUN Rscript -e "BiocManager::install(c('hta20cdf'))"
## # 356 Mb.  Not necessary.
## #RUN Rscript -e "BiocManager::install(c('pd.hta.2.0'))"
## 
## RUN Rscript -e "BiocManager::install(c('limma', 'edgeR', 'DESeq2'))"
## RUN Rscript -e "BiocManager::install(c('impute'))"
## RUN Rscript -e "BiocManager::install(c('GSVA'))"
## RUN Rscript -e "BiocManager::install(c('copynumber'))"
## RUN Rscript -e "BiocManager::install(c('LPE'))"
## RUN Rscript -e "BiocManager::install(c('bioDist'))"
## # Need to install samr after impute installed.
## RUN install2.r -e samr
## 
## # minfi library for Illumina methylation microarrays.
## RUN Rscript -e "BiocManager::install('minfi')"
## 
## RUN Rscript -e "BiocManager::install(c('org.Hs.eg.db', 'org.Mmu.eg.db'))"
## RUN Rscript -e "BiocManager::install('clusterProfiler')"
## 
## 
## 
## 
## # Seurat
## # umap-learn requires llvmlite
## # Version >= 0.32.0 requires python3
## RUN pip install llvmlite==0.31.0
## # umap-learn >= 0.4.0 requires python3.
## RUN pip install umap-learn==0.3.9
## RUN pip install virtualenv   # needed for reticulate
## RUN install2.r -e hdf5r Seurat
##  
## # umap package
## RUN install2.r -e umap
##  
##  
## # SingleR
## # Should install SingleR after Seurat.
## # Need to install tibble and dbplyr manually before scRNAseq or there
## # will be errors.
## RUN install2.r -e tibble
## RUN install2.r -e dbplyr
## # Install scRNAseq package for training data sets.
## RUN Rscript -e "BiocManager::install('scRNAseq')"
## # Github version is obsolete.  Now on Bioconductor.
## #RUN Rscript -e "devtools::install_github('dviraran/SingleR')"
## # Devel version requires R 4.
## #RUN Rscript -e "BiocManager::install(version='devel'); BiocManager::install('SingleR')"
## # BiocNeighbors (required by SingleR) is hard to install.
## # - Use Bioconductor release 3.10 for R 3.6.
## #   https://www.bioconductor.org/about/release-announcements/
## # - Release 3.10 requires RcppAnnoy version 0.0.16.
## #   Incompatible with 0.0.17.
## #RUN Rscript -e "devtools::install_version('RcppAnnoy', version='0.0.16')"
## RUN Rscript -e "devtools::install_version('RcppAnnoy')"
## RUN Rscript -e "BiocManager::install('BiocNeighbors')"
## RUN Rscript -e "BiocManager::install('SingleR')"
## 
## # Cache the reference data for SingleR.
## RUN Rscript -e "BiocManager::install('celldex')"
## #RUN Rscript -e "library(celldex) ; \
## #  HumanPrimaryCellAtlasData();"
## #  library(SingleR) ; \
## 
## RUN mkdir -p /usr/local/ExperimentHub
## RUN Rscript -e "library(ExperimentHub); \
##   setExperimentHubOption(\"CACHE\", \"/usr/local/ExperimentHub\"); \
##   library(celldex); \
##   HumanPrimaryCellAtlasData(); \
##   BlueprintEncodeData(); \
##   DatabaseImmuneCellExpressionData(); \
##   NovershternHematopoieticData(); \
##   MonacoImmuneData(); \
##   ImmGenData(); \
##   MouseRNAseqData()"
## 
## # Not quite the same thing as in ExperimentHub.
## #ENV A=dviraran/SingleR/blob/master/data
## #RUN mkdir -p /usr/local/SingleR/data && \
## #  cd /usr/local/SingleR/data && \
## #  curl -L -o blueprint_encode.rda \
## #    https://github.com/${A}/blueprint_encode.rda?raw=true
## ## load("/usr/local/SingleR/data/blueprint_encode.rda")
## 
## 
## # Harmony
## RUN Rscript -e "devtools::install_github('immunogenomics/harmony')"
## 
## # maftools
## RUN Rscript -e "BiocManager::install('maftools')"
## 
## 
## ## RcppAnnoy
## ## Need this to load Jinfeng's RDS files.
## ## Actually, already installed above.
## #RUN Rscript -e "install.packages('RcppAnnoy')"
## 
## 
## # VennDiagram
## RUN Rscript -e "devtools::install_version('VennDiagram', version='1.6.9')"
## 
## 
## # Install genomicode.
## 
## COPY changlab /usr/local/src/changlab
## 
## ENV USER=root
## RUN cd /usr/local/src/changlab && \
##   python setup.py install
## # rmdir /changlab
## 
## # /usr/local/changlab/Rlib/
## RUN mkdir -p /usr/local/changlab/
## COPY changlab/Rlib /usr/local/changlab/Rlib
## 
## # /usr/local/genomicode/bin/
## RUN mkdir -p /usr/local/genomicode/bin && \
##   cp -Rp /usr/lib64/python2.7/dist-packages/genomicode/bin/* /usr/local/genomicode/bin
## ENV PATH /usr/local/genomicode/bin:$PATH
## 
## # /usr/local/genomidata/
## RUN mkdir -p /usr/local/genomidata/
## COPY psid2platform /usr/local/genomidata/psid2platform
## COPY convert_platform /usr/local/genomidata/convert_platform
## 
##   
## ## # Clean up apt-get.
## ## RUN cd / && \
## ##     apt-get remove --purge -y $BUILDDEPS &&\
## ##     apt-get autoremove -y && \
## ##     apt-get autoclean -y && \
## ##     rm -rf /var/lib/apt/lists/*
## 

CMD ["R"]
