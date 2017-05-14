FROM ubuntu:16.04

MAINTAINER mirnylib

RUN apt-get update --fix-missing && \
  apt-get install -q -y wget curl bzip2 libbz2-dev git build-essential zlib1g-dev locales vim


# Set the locale
RUN locale-gen en_US.UTF-8  
ENV LANG en_US.UTF-8  
ENV LANGUAGE en_US:en  
ENV LC_ALL en_US.UTF-8     

# Install conda
RUN curl -LO http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda3 -b && \
    rm Miniconda3-latest-Linux-x86_64.sh 
ENV PATH=/miniconda3/bin:${PATH}

# Install conda dependencies
ADD environment.yml /
RUN pwd
RUN conda config --set always_yes yes --set changeps1 no && \
    conda config --add channels conda-forge && \
    conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels golobor && \ 
    conda config --get && \
    conda update -q conda && \
    conda info -a && \
    conda env update -q -n root --file environment.yml && \
    conda clean --tarballs --index-cache --lock

