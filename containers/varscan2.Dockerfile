# Heavily modified from rocker-versioned Dockerfile.

# /usr/local/bin/
#   bam-readcount
#   fpfilter.pl
#   mergeSegments.pl
#   VarScan.v2.3.9.jar
#   VarScan.jar
# /usr/local/genomicode/bin/
#   pyvarscan.py
#   [...]


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




# Requires Python.

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



# bam-readcount

RUN apt-get install -y --no-install-recommends \
    cmake \
    make \
    g++ \
    git

RUN apt-get install -y --no-install-recommends \
    patch


RUN cd /usr/local/src/ && \
  git clone https://github.com/genome/bam-readcount.git
RUN cd /usr/local/src/bam-readcount && \
  cmake . && \
  make && \
  cp bin/bam-readcount /usr/local/bin



# Varscan2

ENV A=VarScan.v2.3.9
RUN cd /usr/local/bin && \
  curl -L https://sourceforge.net/projects/varscan/files/${A}.jar/download > \
    ${A}.jar && \
  ln -s /usr/local/bin/${A}.jar /usr/local/bin/VarScan.jar

# fpfilter.pl

# Use the old version of fpfilter.pl from:
# https://sourceforge.net/projects/varscan/files/scripts/
# (The new version does not work with VarScan files.)

ENV A=sourceforge.net/projects/varscan/files
RUN cd /usr/local/bin && \
  curl -L https://${A}/scripts/fpfilter.pl/download > fpfilter.pl && \
  chmod 755 fpfilter.pl && \
  curl -L https://${A}/scripts/fpfilter-help.txt/download > \
    fpfilter-help.txt && \
  curl -L https://${A}/scripts/mergeSegments.pl/download > mergeSegments.pl &&\
  chmod 755 mergeSegments.pl && \
  curl -L https://${A}/scripts/mergeSegments-help.txt/download > \
    mergeSegments-help.txt

# Need to install PERL libraries??
#RUN cd /usr/local/src/ && \
#  git clone https://github.com/genome/fpfilter-tool && \
#  cp /usr/local/src/fpfilter-tool/fpfilter.pl /usr/local/bin




# genomicode / pyvarscan

RUN pip install setuptools wheel
RUN pip install numpy
RUN pip install psutil

COPY changlab /usr/local/src/changlab

ENV USER=root
RUN cd /usr/local/src/changlab && \
  python setup.py install
# rmdir /changlab

# /usr/local/genomicode/bin/
RUN mkdir -p /usr/local/genomicode/bin && \
  cp -Rp /usr/local/lib/python2.7/dist-packages/genomicode/bin/* /usr/local/genomicode/bin
ENV PATH /usr/local/genomicode/bin:$PATH

  
# Clean up apt-get.
RUN cd / && \
    apt-get remove --purge -y $BUILDDEPS &&\
    apt-get autoremove -y && \
    apt-get autoclean -y && \
    rm -rf /var/lib/apt/lists/*

#CMD ["R"]



