# Stuff for running subread/featureCounts.

# /usr/local/bin/
#   exactSNP
#   subindel
#   sublong
#   subread-buildindex
#   featureCounts
#   subjunc
#   subread-align
#   utilities/

FROM debian:stretch

LABEL maintainer="Jeffrey Chang <jeffrey.t.chang@uth.tmc.edu>"

ENV BUILD_DATE="$(TZ='America/Los_Angeles' date -I)" \
    LC_ALL=en_US.UTF-8 \
    LANG=en_US.UTF-8 \
    TERM=xterm


# Basic stuff.
# Don't know why, but apt-get update here prevents 404.
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
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
    openssl \
    libssl-dev \
    curl \
    libcurl4-openssl-dev



# Download subread.

ENV A=subread \
    B=subread-2.0.3-Linux-x86_64.tar.gz
RUN cd /usr/local/src && \
  curl -k -L https://sourceforge.net/projects/${A}/files/${B}/download > ${B}

RUN cd /usr/local/src && \
  tar xfz ${B} && \
  cd subread-2.0.3-Linux-x86_64 && \
  cp -rp bin/* /usr/local/bin/ && \
  rm -f ${B}


## Clean up apt-get.
#RUN cd / && \
#    apt-get remove --purge -y $BUILDDEPS &&\
#    apt-get autoremove -y && \
#    apt-get autoclean -y && \
#    rm -rf /var/lib/apt/lists/*

#ENTRYPOINT ["velocyto"]
