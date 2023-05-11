# Heavily modified from rocker-versioned Dockerfile.

# /usr/local/src/gatk-4.1.8.0/
# /usr/local/bin/gatk
# /usr/local/picard.jar

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
    vim \
    openssl \
    libssl-dev \
    curl \
    libcurl4-openssl-dev


# Requires Python 2.6+.

RUN apt-get update && apt-get install -y --no-install-recommends \
    python2.7-dev \
    python-pip


# Requires Java 8.

# Install OpenJDK-8
RUN apt-get install -y --no-install-recommends \
    openjdk-8-jdk \
    ant

ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/
RUN export JAVA_HOME=/usr/lib/jvm/java-8-openjdk-amd64/





# GATK requires R and some R libraries.

# Install Ghostscript for plotting PDF.
RUN apt-get update && apt-get install -y --no-install-recommends \
    gsfonts \
    libjpeg62-turbo \
    libcairo2-dev \
    libpango1.0-dev \
    libpangocairo-1.0-0 \
    libpng-dev \
    libtiff5-dev \
    ghostscript \
    libgs-dev \
    libx11-dev \
    x11proto-core-dev \
    libxt-dev

RUN apt-get update && apt-get install -y --no-install-recommends \
    libreadline-dev \
    libblas-dev \
    liblapack-dev \
    libicu-dev \
    liblzma5

# Installs R 3.3.3.
RUN apt-get update && apt-get install -y --no-install-recommends \
    r-base


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
    build-essential \
    gcc \
    make

RUN apt-get update && apt-get install -y --no-install-recommends \
    zlib1g-dev
RUN Rscript -e "install.packages('git2r')"
RUN apt-get update && apt-get install -y --no-install-recommends \
    libxml2-dev
RUN Rscript -e "install.packages('roxygen2')"
RUN Rscript -e "install.packages('rversions')"
RUN Rscript -e "install.packages('usethis')"
RUN Rscript -e "install.packages('devtools')"
RUN Rscript -e "require(devtools); install_version('caTools', version='1.17.1.4')"
RUN Rscript -e "install.packages('gplots')"
RUN Rscript -e "install.packages('ggplot2')"
RUN Rscript -e "install.packages('gsalib')"
RUN Rscript -e "install.packages('reshape')"


# Other miscellaneous R libraries.

RUN Rscript -e "install.packages('optparse')"
# Version 1.13.0 gives error about not able to load datatable.so.
RUN Rscript -e "require(devtools); install_version('data.table', version='1.12.8')"



# Install gatk 4.
RUN apt-get update && apt-get install -y --no-install-recommends \
    unzip

RUN mkdir -p /usr/local/src && \
    cd /usr/local/src && \
    curl -L -O https://github.com/broadinstitute/gatk/releases/download/4.1.8.0/gatk-4.1.8.0.zip && \
    unzip gatk-4.1.8.0.zip && \
    rm -f gatk-4.1.8.0.zip

RUN ln -s /usr/local/src/gatk-4.1.8.0/gatk /usr/local/bin


# Install Picard, too.

RUN cd /usr/local && \
    curl -L -O https://github.com/broadinstitute/picard/releases/download/2.25.5/picard.jar





CMD ["gatk"]
