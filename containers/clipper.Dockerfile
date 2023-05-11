# Heavily modified from rocker-versioned Dockerfile.


# /usr/local/src/merge_peaks/bin/
# /usr/local/bin/trimmomatic-0.39.jar
# /usr/local/src/Trimmomatic-0.39/adapters/


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
    make \
    git



# Install R 3.3.3.
# Need to install R before PERL.

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



# Python 2
RUN apt-get update && apt-get install -y --no-install-recommends \
    python2.7 \
    python2.7-dev \
    python-pip

RUN pip install setuptools wheel


# Install samtools and pysam.

# PySam is a bit complicated.
RUN apt-get update && apt-get install -y --no-install-recommends \
    libbz2-dev    \
    liblzma-dev   \
    zlib1g-dev    \
    bzip2

# Cython needs to be installed before pysam.
RUN pip install Cython && \
    pip install pysam

ENV HTSLIB_VERSION=1.10.2
RUN cd /usr/local/src && \
    curl -L -O https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 && \
    bzcat htslib-${HTSLIB_VERSION}.tar.bz2 | tar xf - && \
    rm -f htslib-${HTSLIB_VERSION}.tar.bz2 && \
    cd htslib-${HTSLIB_VERSION} && \
    ./configure && \
    make && \
    make install && \
    rm -rf htslib-${HTSLIB_VERSION}


RUN apt-get install -y --no-install-recommends \
    libncurses5-dev

ENV SAMTOOLS_VERSION=1.10
RUN cd /usr/local/src && \
    curl -L -O https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    bzcat samtools-${SAMTOOLS_VERSION}.tar.bz2 | tar xf - && \
    rm -f samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION} && \
    ./configure && \
    make && \
    make install && \
    cd && \
    rm -rf /usr/local/src/samtools-${SAMTOOLS_VERSION}


# Python 2 libraries and programs.

RUN pip install \
    numpy \
    pandas \
    scipy \
    matplotlib \
    sklearn \
    bx-python

RUN pip install gffutils
RUN pip install pybedtools
RUN pip install statsmodels
RUN pip install HTSeq


# Python 3 libraries and programs.

# apt-get get Python 3.5.3.  However, Clipper requires scikit-learn
# 0.23.2, which requires Python 3.6.  Install a recent version of
# Python from source.

RUN apt-get update && apt-get install -y --no-install-recommends \
    libreadline-gplv2-dev\
    libncursesw5-dev \
    libssl-dev \
    libsqlite3-dev \
    tk-dev \
    libgdbm-dev \
    libc6-dev \
    libbz2-dev \
    libffi-dev \
    zlib1g-dev

RUN cd /usr/local/src && \
    curl -k -L -O https://www.python.org/ftp/python/3.8.6/Python-3.8.6.tgz && \
    tar xzf Python-3.8.6.tgz && \
    rm -f Python-3.8.6.tgz && \
    cd Python-3.8.6 && \
    ./configure --enable-optimizations && \
    make && \
    make altinstall
# TODO: Merge this with previous part.  No make altinstall
RUN cd /usr/local/src/Python-3.8.6 && \
    make install

RUN cd /usr/local/src && \
    curl -k -L https://bootstrap.pypa.io/get-pip.py -o get-pip.py && \
    python3 get-pip.py

RUN pip3 install setuptools wheel
RUN pip3 install cutadapt
RUN pip3 install umi_tools



# PERL

# Need PERL 5.10.1.  According to protocol:
#   changes to sorting in 5.22 may work but cause slightly different
#   peak output
# make test fails.  Seems like Errno in 5.10 is broken.

#ENV V=5.10.1
#ENV V=5.14.1
ENV V=5.32.0
RUN cd /usr/local/src && \
    curl -L -O https://www.cpan.org/src/5.0/perl-${V}.tar.gz && \
    zcat perl-${V}.tar.gz | tar xf - && \
    rm -f perl-${V}.tar.gz && \
    cd perl-${V} && \
    ./Configure -des -Dcc=gcc -Dglibpth='/usr/lib/x86_64-linux-gnu' -Dplibpth='/usr/lib/x86_64-linux-gnu' -Dlibpth='/usr/lib/x86_64-linux-gnu'
RUN cd /usr/local/src/perl-${V} && \
    make && \
    make install && \
    cd && \
    rm -rf /usr/local/src/perl-${V}


# FastQC
# https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc

# Install OpenJDK-8
RUN apt-get install -y --no-install-recommends \
    openjdk-8-jdk \
    ant

ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/
RUN export JAVA_HOME=/usr/lib/jvm/java-8-openjdk-amd64/

RUN cd /usr/local/src && \
    curl -L -O https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip && \
    unzip fastqc_v0.11.9.zip && \
    chmod 755 FastQC/fastqc && \
    ln -s /usr/local/src/FastQC/fastqc /usr/local/bin/ && \
    rm -rf fastqc_v0.11.9.zip



# STAR

RUN apt-get update && apt-get install -y --no-install-recommends \
    git

RUN cd /usr/local/src && \
    git clone https://github.com/alexdobin/STAR.git && \
    cd STAR/source && \
    make STAR && \
    cp -Rp /usr/local/src/STAR/bin/Linux_x86_64_static/STAR* /usr/local/bin/&&\
    cd && \
    rm -rf /usr/local/src/STAR



# fastq-tools

RUN apt-get install -y --no-install-recommends \
    libpcre3-dev

RUN cd /usr/local/src && \
    curl -k -L -O http://homes.cs.washington.edu/~dcjones/fastq-tools/fastq-tools-0.8.tar.gz && \
    tar xvfz fastq-tools-0.8.tar.gz && \
    cd fastq-tools-0.8 && \
    ./configure && \
    make && \
    make install && \
    rm -rf fastq-tools-0.8.tar.gz && \
    rm -rf fastq-tools-0.8


# Bedtools

RUN cd /usr/local/src && \
    curl -k -L -O https://github.com/arq5x/bedtools2/releases/download/v2.29.2/bedtools.static.binary && \
    chmod 755 bedtools.static.binary && \
    mv bedtools.static.binary /usr/local/bin/bedtools

RUN cd /usr/local/src && \
    curl -k -L -O https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz && \
    tar -zxvf bedtools-2.29.1.tar.gz && \
    cd bedtools2 && \
    make && \
    cp -p bin/* /usr/local/bin/


# bedToBigBed

RUN apt-get install -y --no-install-recommends \
    libpng-dev \
    uuid-dev \
    default-libmysqlclient-dev

ENV MACHTYPE=x86_64
RUN mkdir -p /root/bin/${MACHTYPE}

RUN cd /usr/local/src && \
    curl -k -L -O https://github.com/ucscGenomeBrowser/kent/archive/v377_base.tar.gz && \
    tar xvfz v377_base.tar.gz && \
    rm -f v377_base.tar.gz && \
    cd /usr/local/src/kent-377_base/src/lib && \
    make && \
    cd /usr/local/src/kent-377_base/src/jkOwnLib && \
    make && \
    cd /usr/local/src/kent-377_base/src/htslib && \
    make
RUN cd /usr/local/src/kent-377_base/src/utils/bedToBigBed && \
    make && \
    cp -p /root/bin/${MACHTYPE}/bedToBigBed /usr/local/bin/
RUN cd /usr/local/src/kent-377_base/src/utils/bedGraphToBigWig && \
    make && \
    cp -p /root/bin/${MACHTYPE}/bedGraphToBigWig /usr/local/bin/
# Need jkhgap.a for bedSort.
RUN cd /usr/local/src/kent-377_base/src && \
    make libs
RUN cd /usr/local/src/kent-377_base/src/hg/bedSort && \
    make && \
    cp -p /root/bin/${MACHTYPE}/bedSort /usr/local/bin/




# Trimmomatic.
RUN cd /usr/local/src && \
    curl -O http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip && \
    unzip Trimmomatic-0.39.zip && \
    cp -p /usr/local/src/Trimmomatic-0.39/trimmomatic-0.39.jar /usr/local/bin
RUN ln -s /usr/local/bin/trimmomatic-0.39.jar /usr/local/bin/trimmomatic.jar



# Yeo Lab tools

# CWL

RUN apt-get install -y --no-install-recommends \
    cwltool

RUN cd /usr/local/src && \
    git clone https://github.com/YeoLab/makebigwigfiles && \
    cd makebigwigfiles && \
    python setup.py build && \
    python setup.py install

# Clipper seems to require Python3 now.
RUN pip3 install Cython
# Clipper requires scikit-learn 0.23.2, which requires Python 3.6, but
# apt-get only has Python 3.5.  Need to install Python 3.6.
RUN cd /usr/local/src && \
    git clone https://github.com/YeoLab/clipper

# Copy Han's rRNA data into the container.
COPY dat01/rRNA.AS.STRUCTURE.COMPILED.gff /usr/local/src/clipper/clipper/data
COPY dat01/rRNA_exons.bed /usr/local/src/clipper/clipper/data/regions
COPY dat01/rRNA_genes.bed /usr/local/src/clipper/clipper/data/regions

RUN cd /usr/local/src/clipper && \
    python3 setup.py build && \
    python3 setup.py install

RUN cd /usr/local/src && \
    git clone https://github.com/YeoLab/eclipdemux && \
    cd eclipdemux && \
    python setup.py build && \
    python setup.py install


# Install merge_peaks and dependencies.

# Set variable to configure CPAN with default options.
ENV PERL_MM_USE_DEFAULT=1

RUN cpan Test::Regexp
RUN cpan Regexp::Common
RUN cpan YAML
RUN cpan IO::Pty
RUN cpan inc::latest

# IPC::Run (required for Statistics::R) has errors.  Not sure how to
# fix.  Need to force install.
#RUN cpan IPC::Run
# 20200505.0 fails with:
# No such file or directory: getprotobyname('tcp') at t/win32_compile.t line 90.
# Compilation failed in require at t/win32_compile.t line 90.
# BEGIN failed--compilation aborted at t/win32_compile.t line 90.
# Makefile:878: recipe for target 'test_dynamic' failed
#RUN cpan TODDR/IPC-Run-20180523.0.tar.gz    # fails
#RUN cpan TODDR/IPC-Run-0.99.tar.gz          # fails (oldest version)
RUN cpan -f IPC::Run
# perldoc -l IPC::Run    # Check if installed

RUN perl -MCPAN -e 'install Statistics::Basic'
RUN perl -MCPAN -e 'install Statistics::Distributions'
RUN cpan FANGLY/Statistics-R-0.34.tar.gz



RUN pip3 install \
    setuptools \
    numpy \
    matplotlib
RUN pip3 install cwl

RUN cd /opt/ && \
    curl -L -O https://github.com/nboley/idr/archive/2.0.2.zip && \
    unzip /opt/2.0.2.zip
RUN cd /opt/idr-2.0.2/ && \
    python3 setup.py install

RUN cd /usr/local/src && \
    git clone https://github.com/YeoLab/merge_peaks
RUN for i in /usr/local/src/merge_peaks/bin/perl/*.pl; do \
      ln -s $i /usr/local/bin/; \
    done


# Save space.
RUN rm -rf /usr/local/src/Python-3.8.6
RUN rm -rf /usr/local/src/kent-377_base
RUN rm /usr/local/src/bedtools-2.29.1.tar.gz && \
    rm -rf /usr/local/src/bedtools2
RUN rm -rf /usr/local/src/clipper

