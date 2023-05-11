# Heavily modified from rocker-versioned Dockerfile.

# /usr/local/beast/
# /usr/local/changlab/Rlib/

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


# Miscellaneous basic stuff
RUN Rscript -e "install.packages('data.table')"
RUN Rscript -e "install.packages('dplyr')"
RUN Rscript -e "install.packages('ggplot2')"
RUN Rscript -e "install.packages('TreeTools')"


# Install ggtree.

RUN Rscript -e "install.packages('BiocManager')"
RUN Rscript -e "BiocManager::install('ggtree')"
RUN Rscript -e "install.packages('ggsci')"
RUN Rscript -e "install.packages('ggnewscale')"


# Other packages

RUN Rscript -e "install.packages('matrixStats')"
RUN Rscript -e "install.packages('geomnet')"
RUN Rscript -e "install.packages('phytools')"



# Install BEAGLE.

RUN apt-get install -y --no-install-recommends \
    build-essential \
    automake \
    autoconf \
    libtool

RUN apt-get update && apt-get install -y --no-install-recommends \
    git

# Somehow, this causes problems with OpenCL.
# OpenCL error: Unknown error from file <GPUInterfaceOpenCL.cpp>, line 115.
#RUN apt-get update && apt-get install -y --no-install-recommends \
#    ocl-icd-opencl-dev

RUN cd /usr/local/src && \
  git clone --depth=1 https://github.com/beagle-dev/beagle-lib.git && \
  cd beagle-lib && \
  ./autogen.sh && \
  ./configure --prefix=/usr && \
  make install


# Babette imports
RUN Rscript -e "install.packages('phangorn')"
RUN Rscript -e "install.packages('remotes')"
RUN Rscript -e "install.packages('stringr')"
RUN Rscript -e "install.packages('testit')"

# Babette suggests
RUN Rscript -e "install.packages('ape')"
RUN Rscript -e "install.packages('hunspell')"
RUN Rscript -e "install.packages('lintr')"
RUN Rscript -e "install.packages('nLTT')"
RUN Rscript -e "install.packages('rmarkdown')"
RUN Rscript -e "install.packages('spelling')"
RUN Rscript -e "install.packages('plotrix')"

# BEAST2 / Babette
# RUN Rscript -e "install.packages('beautier')"
# RUN Rscript -e "install.packages('beastier')"
# RUN Rscript -e "install.packages('mauricer')"
# RUN Rscript -e "install.packages('tracerer')"
# RUN Rscript -e "install.packages('babette')"

RUN Rscript -e "install.packages('devtools')"

# force refresh github
ENV GITHUB_DATE=210815

# Install beastier.

# When running large trees, will get an error:
#   java.lang.reflect.InvocationTargetException
#   [...]
#   Caused by: java.lang.StackOverflowError
# Triggered in:
#   beastier::is_beast2_input_file
#   create_beast2_validate_cmd
#   create_beast2_validate_cmd_jar
#   /usr/lib/jvm/java-11-openjdk-amd64/bin/java \
#     -cp /usr/local/beast/lib/launcher.jar \
#     beast.app.beastapp.BeastLauncher -validate <file>
# To fix, increase the default stack size for java.

RUN cd /usr/local/src && \
  curl -L -O \
    https://github.com/ropensci/beastier/archive/refs/heads/master.zip && \
  unzip master.zip
# Add -Xss2m before -cp.
RUN perl -i -p -e 's/"-cp"/"-Xss2m", "-cp"/' \
  /usr/local/src/beastier-master/R/create_beast2_validate_cmd.R
RUN perl -i -p -e 's/"-cp"/"-Xss2m", "-cp"/' \
  /usr/local/src/beastier-master/R/create_beast2_continue_cmd_from_options.R


RUN Rscript -e "devtools::install('/usr/local/src/beastier-master')"
# RUN Rscript -e "remotes::install_github('ropensci/beastier')"
RUN Rscript -e "remotes::install_github('ropensci/beautier')"
RUN Rscript -e "remotes::install_github('ropensci/mauricer')"
RUN Rscript -e "remotes::install_github('ropensci/tracerer')"
RUN Rscript -e "remotes::install_github('ropensci/babette')"




# Install BEAST2 in "/usr/local".

#RUN cd /usr/local/src && \
#    curl -O -L https://github.com/CompEvol/beast2/releases/download/v2.6.3/BEAST_with_JRE.v2.6.3.Linux.tgz && \
#  tar xfz BEAST_with_JRE.v2.6.3.Linux.tgz && \
#  rm -rf BEAST_with_JRE.v2.6.3.Linux.tgz

# beastier::get_default_base2_folder calls
#   rappdirs::user_data_dir
#   which checks USER_DATA_DIR
env R_USER_DATA_DIR=/usr/local
RUN Rscript -e 'remotes::install_github("richelbilderbeek/beastierinstall")'
RUN Rscript -e "library(beastier); beastierinstall::install_beast2()"
RUN rm -f /usr/local/BEAST.*.tgz



# TODO: Move these up somewhere.
RUN Rscript -e "install.packages('tidytree')"



# Install changlab R libraries

# /usr/local/changlab/Rlib/
RUN mkdir -p /usr/local/changlab/
COPY changlab/Rlib /usr/local/changlab/Rlib



# Clean up src code.
RUN rm -rf /usr/local/src/*


# Clean up apt-get.
RUN cd / && \
    apt-get remove --purge -y $BUILDDEPS &&\
    apt-get autoremove -y && \
    apt-get autoclean -y && \
    rm -rf /var/lib/apt/lists/*

CMD ["R"]
