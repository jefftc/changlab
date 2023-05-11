FROM ubuntu:16.04

LABEL maintainer="Jeffrey Chang <jeffrey.t.chang@uth.tmc.edu>"

ENV BUILD_DATE="$(TZ='America/Los_Angeles' date -I)" \
    LC_ALL=en_US.UTF-8 \
    LANG=en_US.UTF-8 \
    TERM=xterm

# Install some basic packages.
RUN apt-get update --fix-missing && \
    apt-get install -y --no-install-recommends \
    build-essential

# Miscellaneous useful things.
RUN apt-get install -y --no-install-recommends \
    bash-completion \
    ca-certificates \
    locales \
    less \
    file \
    git

# Needed to fix singularity image.
RUN apt-get install -y --no-install-recommends \
    vim-tiny

RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen && \
    /usr/sbin/locale-gen en_US.UTF-8 && \
    /usr/sbin/update-locale LANG=en_US.UTF-8


# Pre-requisites for R.
RUN apt-get install -y --no-install-recommends \
    # R
    zlib1g-dev \
    libpng-dev \
    # For plotting PDF.
    ghostscript \
    libgs-dev

# Install perl.
RUN apt-get install -y --no-install-recommends \
    perl

# Install R.
RUN apt-get install -y --no-install-recommends \
    r-base

# Install OpenJDK-8
RUN apt-get install -y --no-install-recommends \
    openjdk-8-jdk \
    ant

## Fix certificate issues
#RUN apt-get update && \
#    apt-get install ca-certificates-java && \
#    apt-get clean && \
#    update-ca-certificates -f;

# Setup JAVA_HOME -- useful for docker commandline
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/
RUN export JAVA_HOME=/usr/lib/jvm/java-8-openjdk-amd64/



# Install VarDictJava.
RUN cd /usr/src && \
  git clone --recursive \
    https://github.com/AstraZeneca-NGS/VarDictJava.git && \
  cd /usr/src/VarDictJava && \
  ./gradlew clean installDist && \
  mv /usr/src/VarDictJava/build/install/VarDict /usr/local

# Doesn't work.
#  ./gradlew clean javadoc

ENV PATH /usr/local/VarDict/bin:$PATH
RUN export PATH

# /usr/local/VarDict/bin/
#   VarDict
#   VarDict.bat
#   testsomatic.R
#   teststrandbias.R
#   var2vcf_paired.pl
#   var2vcf_valid.pl


## Clean up apt-get.
#RUN cd / && \
#    apt-get remove --purge -y $BUILDDEPS &&\
#    apt-get autoremove -y && \
#    apt-get autoclean -y && \
#    rm -rf /var/lib/apt/lists/*


#ENTRYPOINT ["VarDict"]
