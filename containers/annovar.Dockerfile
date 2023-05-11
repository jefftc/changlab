# Heavily modified from rocker-versioned Dockerfile.

# Need to download Annovar from:
# https://annovar.openbioinformatics.org/en/latest/
# /usr/local/annovar/
#   annotate_variation.pl
#   humandb/


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


# PERL

RUN apt-get update && apt-get install -y --no-install-recommends \
  build-essential \
  make \
  gcc

# Somehow, the system PERL doesn't work.  Get error:
# 'Can't find any loadable formatter class in Pod::Perldoc::Toterm
# Pod::Perldoc::Toterm Pod::Perldoc::ToTerm Pod::Perldoc::ToTERM
# Pod::Simple::term Pod::Simple::term Pod::Simple::Term
# Pod::Simple::TERM Pod::term Pod::term Pod::Term Pod::TERM
# Pod::Perldoc::Totext Pod::Perldoc::Totext Pod::Perldoc::ToText
# Pod::Perldoc::ToTEXT Pod::Simple::text Pod::Simple::text
# Pod::Simple::Text Pod::Simple::TEXT Pod::text Pod::text Pod::Text
# Pod::TEXT Pod::Perldoc::ToPod?!

# Install my own PERL.

#RUN apt-get update && apt-get install -y --no-install-recommends \
#  perl \
#  perl-doc \
#  libpath-tiny-perl

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


RUN cpan Pod::Simple
RUN cpan Pod::Perldoc
RUN cpan Pod::Usage


# Annovar

RUN mkdir -p /usr/local/src
COPY annovar.latest.tar.gz /usr/local/src
RUN cd /usr/local/src && \
  tar xfz annovar.latest.tar.gz && \
  mv annovar /usr/local

# Download Annovar databases.

RUN apt-get update && apt-get install -y --no-install-recommends \
  wget \
  unzip

RUN cd /usr/local/annovar && \
  ./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene \
  humandb/
RUN cd /usr/local/annovar && \
  ./annotate_variation.pl -buildver hg19 -downdb cytoBand humandb/
RUN cd /usr/local/annovar && \
  ./annotate_variation.pl -buildver hg19 -downdb genomicSuperDups humandb/
RUN cd /usr/local/annovar && \
  ./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar \
  esp6500siv2_all humandb/
RUN cd /usr/local/annovar && \
  ./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar \
  1000g2015aug humandb/
RUN cd /usr/local/annovar && \
  ./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 \
  humandb/
RUN cd /usr/local/annovar && \
  ./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar avsnp147 \
  humandb/
RUN cd /usr/local/annovar && \
  ./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp30a \
  humandb/
RUN cd /usr/local/annovar && \
  ./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar snp138 \
  humandb/
RUN cd /usr/local/annovar && \
  ./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ljb26_all \
  humandb/

ENV PATH /usr/local/annovar:$PATH

  
# Clean up apt-get.
RUN cd / && \
    apt-get remove --purge -y $BUILDDEPS &&\
    apt-get autoremove -y && \
    apt-get autoclean -y && \
    rm -rf /var/lib/apt/lists/*

CMD ["annotate_variation.pl"]



