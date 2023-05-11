FROM ubuntu:16.04

LABEL maintainer="Jeffrey Chang <jeffrey.t.chang@uth.tmc.edu>"

# Install some required packages.
RUN apt-get update --fix-missing && \
    apt-get install -y --no-install-recommends \
    build-essential \
    # General use
    less \
    file \
    # For Anaconda
    wget \
    bzip2 \
    # Bioconductor
    zlib1g-dev \
    libpng-dev \
    # For plotting PDF.
    ghostscript \
    libgs-dev


# Install anaconda.
RUN wget --quiet --no-check-certificate https://repo.anaconda.com/archive/Anaconda3-2019.03-Linux-x86_64.sh -O ~/anaconda.sh && \
    /bin/bash ~/anaconda.sh -b -p /opt/conda && \
    rm ~/anaconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc && \
    find /opt/conda/ -follow -type f -name '*.a' -delete && \
    find /opt/conda/ -follow -type f -name '*.js.map' -delete && \
    /opt/conda/bin/conda clean -afy

ENV PATH /opt/conda/bin:$PATH

# Update Anaconda packages
RUN conda update conda && \
    conda update anaconda && \
    conda update --all

# conda search r-base
RUN conda install \
    # velocyto.py
    numpy \
    scipy \
    cython \
    numba \
    matplotlib \
    scikit-learn \
    h5py \
    click \
    # velocyto.R
    r-base \
    libboost \
    # Pagoda2
    r-kernsmooth

# Add a default CRAN mirror
RUN mkdir -p /opt/conda/lib/R/etc && \
    echo "options(repos=c(CRAN='https://cran.rstudio.com/'), download.file.method='libcurl')" \
    >> /opt/conda/lib/R/etc/Rprofile.site



# Install velocyto.py
RUN pip install velocyto



# Install velocyto.R

# Install Bioconductor
RUN Rscript -e 'if(!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager"); BiocManager::install()' && \
    Rscript -e 'BiocManager::install(c("AnnotationDbi", "BioGenerics", "GO.db", "pcaMethods"))'

# Install devtools
RUN ln -s /bin/tar /bin/gtar && \
    Rscript -e "install.packages('git2r', configure.vars='LIBS=-L/usr/lib/x86_64-linux-gnu')" && \
    Rscript -e "install.packages('devtools')"

# Install velocyto.
RUN Rscript -e "devtools::install_github('velocyto-team/velocyto.R')"
#RUN Rscript -e "library(devtools); install_github('velocyto-team/velocyto.R')"

# Install Pagoda2 (for analyzing velocyto data).
RUN Rscript -e "devtools::install_github('hms-dbmi/pagoda2')"


# So R doesn't generate warning messages from Singularity image.
# Doesn't seem to be a problem from Docker.
ENV LANG=C.UTF-8 \
  LC_COLLATE=C.UTF-8 \
  LC_CTYPE=C.UTF-8 \
  LC_MESSAGES=C.UTF-8 \
  LC_MONETARY=C.UTF-8 \
  LC_NUMERIC=C.UTF-8 \
  LC_TIME=C.UTF-8 \
  LC_ALL=C.UTF-8 


## Clean up apt-get.
#RUN cd / && \
#    apt-get remove --purge -y $BUILDDEPS &&\
#    apt-get autoremove -y && \
#    apt-get autoclean -y && \
#    rm -rf /var/lib/apt/lists/*


#RUN useradd -m user
#USER user
#WORKDIR "/home/user"

#ENTRYPOINT ["velocyto"]
