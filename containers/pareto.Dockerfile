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


# Prerequisites for R.

# Development tools.
RUN apt-get update && apt-get install -y --no-install-recommends \
    g++ \
    gfortran \
    make \
    cmake \
    autoconf 

# Languages
RUN apt-get install -y --no-install-recommends \
    perl \
    tcl8.6-dev \
    tk8.6-dev 

# Compression
RUN apt-get install -y --no-install-recommends \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    unzip \
    zip

# Math
RUN apt-get install -y --no-install-recommends \
    libblas-dev \
    libatlas3-base 

# Plotting stuff
RUN apt-get install -y --no-install-recommends \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev 

# Ghostscript, for making PDFs.
RUN apt-get install -y --no-install-recommends \
    texinfo \
    texlive-extra-utils \
    texlive-fonts-recommended \
    texlive-fonts-extra \
    texlive-latex-recommended \
    ghostscript libgs-dev \
    gsfonts \
    fonts-texgyre

# Miscellaneous libraries
RUN apt-get install -y --no-install-recommends \
    libreadline-dev \
    libpcre3-dev


# Install R

ENV R_VERSION=3.6.0 

# Download R.
RUN cd /usr/local/src && \
    curl -O https://cran.r-project.org/src/base/R-3/R-${R_VERSION}.tar.gz  && \
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

# Doesn't work.  Get "Directory not empty" errors.
#RUN rm -rf /usr/local/src/R-${R_VERSION}


# Add a default CRAN mirror
RUN mkdir -p /usr/local/lib/R/etc && \
    echo "options(repos=c(CRAN='https://cran.rstudio.com/'), download.file.method='libcurl')" \
    >> /usr/local/lib/R/etc/Rprofile.site

## Add a library directory (for user-installed packages)
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
    ln -s /usr/local/lib/R/site-library/littler/examples/install2.r /usr/local/bin/install2.r && \
    ln -s /usr/local/lib/R/site-library/littler/examples/installGithub.r /usr/local/bin/installGithub.r && \
    ln -s /usr/local/lib/R/site-library/littler/bin/r /usr/local/bin/r



# Required or useful for vitkl/ParetoTI.
RUN install2.r -e remotes umap Rtsne


# Install bioconductor libraries
RUN Rscript -e 'if(!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager"); BiocManager::install()'

# Pre-requisite for ParetoTI.
RUN apt-get install -y --no-install-recommends \
    libxml2-dev
RUN Rscript -e 'BiocManager::install("vitkl/ParetoTI", \
  dependencies=c("Depends", "Imports", "LinkingTo"))'




# Needed for anaconda.
RUN apt-get install -y --no-install-recommends \
    wget \
    bzip2

# Install anaconda.
RUN wget --quiet --no-check-certificate https://repo.anaconda.com/archive/Anaconda3-2019.03-Linux-x86_64.sh -O ~/anaconda.sh && \
    /bin/bash ~/anaconda.sh -b -p /opt/conda && \
    rm ~/anaconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc && \
    find /opt/conda/ -follow -type f -name '*.a' -delete && \
    find /opt/conda/ -follow -type f -name '*.js.map' -delete && \
    /opt/conda/bin/conda clean -afy

ENV PATH /opt/conda/bin:$PATH

RUN conda create --yes -n reticulate_PCHA python=3.7.3 pip
RUN conda run -n reticulate_PCHA pip install --upgrade \
  py_pcha scipy keras geosketch pydot umap-learn \
  numpy datetime pandas h5py sklearn \
  tensorflow tensorflow-probability



# Install scrublet for doublet removal.  Recommended in ParetoTI
# documentation.
# Install in base conda environment.

#RUN conda run -n reticulate_PCHA pip install --upgrade \
#  scrublet
RUN pip install --upgrade scrublet



# More generally useful libraries for ParetoTI analysis.
RUN install2.r -e cowplot
RUN Rscript -e 'BiocManager::install("SingleCellExperiment")'


# Clean up apt-get.
RUN cd / && \
    apt-get remove --purge -y $BUILDDEPS &&\
    apt-get autoremove -y && \
    apt-get autoclean -y && \
    rm -rf /var/lib/apt/lists/*

CMD ["R"]

