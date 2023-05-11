# Heavily modified from rocker-versioned Dockerfile.

# TODO: Make sure R doesn't try to load any of the user's libraries.
# Need to make sure .libPaths is:
# .libPaths(c("/usr/local/lib/R/site-library", "/usr/local/lib/R/library"))


# /usr/local/changlab/Rlib/

FROM debian:bullseye

LABEL maintainer="Jeffrey Chang <jeffrey.t.chang@uth.tmc.edu>"

ENV R_VERSION=4.1.2 \
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




# Install R libraries.

RUN apt-get install -y --no-install-recommends \
    libxml2-dev
RUN install2.r -e R.utils 
RUN install2.r -e devtools


# Install Bioconductor
RUN Rscript -e 'if(!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager"); BiocManager::install()'


# Install sva package.

RUN Rscript -e 'BiocManager::install("locfit")'
RUN Rscript -e 'BiocManager::install("edgeR")'
RUN apt-get update && apt-get install -y --no-install-recommends \
    libpng-dev
RUN Rscript -e 'BiocManager::install("png")'
RUN Rscript -e 'BiocManager::install("annotate")'
RUN Rscript -e 'BiocManager::install("genefilter")'
RUN Rscript -e 'BiocManager::install("sva")'

# Install for testing purposes.
RUN Rscript -e 'BiocManager::install("bladderbatch")'


# Install genomicode.

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
