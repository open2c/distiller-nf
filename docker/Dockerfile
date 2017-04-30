FROM ubuntu:16.04

MAINTAINER mirnylib

RUN apt-get update --fix-missing && \
  apt-get install -q -y wget curl bzip2 git build-essential zlib1g-dev locales vim

# Set the locale
RUN locale-gen en_US.UTF-8  
ENV LANG en_US.UTF-8  
ENV LANGUAGE en_US:en  
ENV LC_ALL en_US.UTF-8     

# Install conda
RUN curl -LO http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda3 -b
RUN rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH=/miniconda3/bin:${PATH}

# Install conda dependencies
RUN conda config --set always_yes yes --set changeps1 no
RUN conda config --add channels conda-forge               
RUN conda config --add channels defaults                  
RUN conda config --add channels bioconda                  
RUN conda config --add channels golobor
RUN conda config --get                                    

RUN conda update -q conda
# Useful for debugging any issues with conda
RUN conda info -a

ADD environment.yml /
RUN conda env update -q -n root --file environment.yml

