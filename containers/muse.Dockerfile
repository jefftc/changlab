# Heavily modified from rocker-versioned Dockerfile.

# /usr/local/bin/MuSE
# /opt/MuSE/
#   data.hg19/
#     dbsnp_132_b37.leftAligned.vcf.gz
#     dbsnp_132_b37.leftAligned.vcf.gz.tbi


FROM debian:bullseye

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


# Install MuSE.

ENV MUSE=MuSEv1.0rc_submission_c039ffa
RUN cd /usr/local/src && \
  curl -L -O http://bioinformatics.mdanderson.org/Software/MuSE/${MUSE}
RUN chmod 755 /usr/local/src/${MUSE} && \
  ln -s /usr/local/src/${MUSE} /usr/local/bin/MuSE


# Download and process a dbSNP VCF file.
# Should be bgzip compressed and tabix indexed.

ENV MUSE_DIR=/opt/MuSE
RUN mkdir -p ${MUSE_DIR}/data.hg19
COPY dbsnp_132_b37.leftAligned.vcf.gz ${MUSE_DIR}/data.hg19/
COPY dbsnp_132_b37.leftAligned.vcf.gz.tbi ${MUSE_DIR}/data.hg19/


# Clean up apt-get.
RUN cd / && \
    apt-get remove --purge -y $BUILDDEPS &&\
    apt-get autoremove -y && \
    apt-get autoclean -y && \
    rm -rf /var/lib/apt/lists/*

CMD ["MuSE"]

