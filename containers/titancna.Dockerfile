# Heavily modified from rocker-versioned Dockerfile.

# /usr/local/bin/
#   gcCounter        From HMMcopy
#   mapCounter       From HMMcopy
#   readCounter      From HMMcopy
# /usr/local/data/
#   gc_hg19.wig      From HMMcopy
#   gc.wig           From HMMcopy
#   map_hg19.wig     From HMMcopy
#   map.wig          From HMMcopy
#   
# /usr/local/src/hmmcopy_utils/

FROM debian:buster

LABEL maintainer="Jeffrey Chang <jeffrey.t.chang@uth.tmc.edu>"

ENV R_VERSION=4.1.2 \
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




# Install TitanCNA.

RUN Rscript -e 'if(!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")'
RUN Rscript -e "BiocManager::install('TitanCNA')"

#RUN Rscript -e "install.packages('devtools')"
#RUN Rscript -e "remotes::install_github('gavinha/TitanCNA')"



# Install HMMcopy.

RUN apt-get --allow-releaseinfo-change update && \
  apt-get install -y --no-install-recommends \
  make \
  cmake

RUN apt-get update && apt-get install -y --no-install-recommends \
  git

RUN Rscript -e "install.packages('devtools')"
RUN Rscript -e "BiocManager::install('HMMcopy')"

# Install binary programs, e.g. readCounter, gcCounter.
RUN cd /usr/local/src && \
    git clone https://github.com/shahcompbio/hmmcopy_utils.git
RUN cd /usr/local/src/hmmcopy_utils && \
    cmake . && \
    make

RUN cp -p /usr/local/src/hmmcopy_utils/bin/* /usr/local/bin

## Install wig files.
#RUN mkdir -p /usr/local/data && \
#    cp -p /usr/local/src/hmmcopy_utils/data/*.wig /usr/local/data/





# Install ichorCNA.

RUN Rscript -e "remotes::install_github('broadinstitute/ichorCNA')"


# Install bowtie.
# Used by HMMCopy.

ENV A=bowtie-1.3.1-linux-x86_64
RUN cd /usr/local/src && \
    curl -L https://sourceforge.net/projects/bowtie-bio/files/${A}.zip/download > ${A}.zip
RUN cd /usr/local/src && \
    unzip ${A}.zip && \
    rm ${A}.zip
RUN cp -p /usr/local/src/${A}/bowtie* /usr/local/bin



# For parallel computing.

RUN Rscript -e "install.packages('doMC')"
RUN Rscript -e "install.packages('doMPI')"



# Install changlab R libraries

# /usr/local/changlab/Rlib/
RUN mkdir -p /usr/local/changlab/
COPY changlab/Rlib /usr/local/changlab/Rlib






## # Clean up src code.
## RUN rm -rf /usr/local/src/*


# Clean up apt-get.
RUN cd / && \
    apt-get remove --purge -y $BUILDDEPS &&\
    apt-get autoremove -y && \
    apt-get autoclean -y && \
    rm -rf /var/lib/apt/lists/*

CMD ["R"]
