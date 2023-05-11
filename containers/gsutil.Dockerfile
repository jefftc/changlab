# Heavily modified from rocker-versioned Dockerfile.

# /usr/local/google-cloud-sdk/bin/
#   gsutil

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


# Requires Python with crcmod.

RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    python3-dev \
    python3-pip \
    python3-setuptools

RUN pip3 install wheel
# Error.  Not installed.
#RUN pip3 uninstall crcmod
RUN pip3 install --no-cache-dir -U crcmod


# Install gsutil

#curl -o install.sh https://sdk.cloud.google.com
#RUN curl https://sdk.cloud.google.com | bash
RUN curl -O https://dl.google.com/dl/cloudsdk/channels/rapid/install_google_cloud_sdk.bash && \
  sh install_google_cloud_sdk.bash --install-dir=/usr/local/

RUN ln -s /usr/local/google-cloud-sdk/bin/gsutil /usr/local/bin

CMD ["/usr/local/google-cloud-sdk/bin/gsutil"]
