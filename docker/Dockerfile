FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive

LABEL author="Edoardo Giacopuzzi"
LABEL contact="edoardo.giacopuzzi@fht.org"
LABEL image_version="0.5"
LABEL regenie_version="3.2.7"
LABEL bgenix_version="1.1.7"
LABEL R_version="4.1.0"

# Install software
WORKDIR "/opt"
RUN apt-get update && \
    apt-get install -y \
    libcairo2-dev \
    libxt-dev \
    pandoc \
    git \
    unzip \
    gfortran \
    build-essential \
    automake \
    python3 \
    python3-pip \
    zlib1g-dev \
    libgomp1 \
    procps \
    libx11-6 \
    openjdk-11-jdk \
    libjpeg-dev \
    libpcre2-dev pcre2-utils \
    libxml2-dev \
    libopenblas-dev \
    libblas-dev liblapack-dev \
    autoconf ca-certificates wget libbz2-dev \
    libc6-dev libcurl4-openssl-dev libfreetype6 libgsl-dev liblzma-dev \
    libncurses5-dev libperl-dev libreadline-dev libssl-dev libz-dev \
&& ln -s /usr/bin/python3 /usr/bin/python \
&& wget http://security.ubuntu.com/ubuntu/pool/main/i/icu/libicu60_60.2-3ubuntu3.2_amd64.deb \
&& apt-get install ./libicu60_60.2-3ubuntu3.2_amd64.deb \
&& apt-get -y autoremove \
&& apt-get -y clean all \
&& rm -rf /var/cache 

# Install regenie
WORKDIR /opt
ENV REGENIE_VERSION "3.2.7"
RUN mkdir regenie && cd regenie && \
    wget https://github.com/rgcgithub/regenie/releases/download/v${REGENIE_VERSION}/regenie_v${REGENIE_VERSION}.gz_x86_64_Linux.zip && \
    unzip -q regenie_v${REGENIE_VERSION}.gz_x86_64_Linux.zip && \
    rm regenie_v${REGENIE_VERSION}.gz_x86_64_Linux.zip && \
    mv regenie_v${REGENIE_VERSION}.gz_x86_64_Linux regenie && \
    chmod +x regenie
ENV PATH="/opt/regenie/:${PATH}"

#Install plink
WORKDIR /opt
RUN mkdir plink && cd plink \
    && wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20220402.zip \
    && unzip -q plink_linux_x86_64_20220402.zip \
    && rm plink_linux_x86_64_20220402.zip \
    && chmod a+x plink \
    && mv plink /usr/local/bin/plink

# Install plink2
WORKDIR /opt
RUN wget https://s3.amazonaws.com/plink2-assets/alpha3/plink2_linux_avx2_20220814.zip \
    && unzip -q plink2_linux_avx2_20220814.zip \
    && rm plink2_linux_avx2_20220814.zip \
    && chmod a+x plink2 \
    && mv plink2 /usr/local/bin/plink2

# Install htslib
WORKDIR "/opt"
RUN wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2 \
    && tar -xjf htslib-1.16.tar.bz2 \
    && cd htslib-1.16 \
    && autoreconf -i \
    && ./configure \
    && make \
    && make install \
    && cd .. \
    && rm -rf htslib-1.16*

# Install bedtools
WORKDIR "/opt"
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary \
    && mv bedtools.static.binary bedtools \
    && chmod a+x bedtools \
    && ln -s /opt/bedtools /usr/local/bin/bedtools

# Install bgen tools
WORKDIR "/opt"
RUN wget http://code.enkre.net/bgen/tarball/release/bgen.tgz && \
    tar -zxvf bgen.tgz && \
    cd bgen.tgz && \
    ./waf configure && \
    ./waf && \
    ./waf install && \
    cd .. && \
    rm -rf bgen.tgz

# Prepare R environment
WORKDIR "/opt"
ENV R_LIBS /opt/R/libs
ENV R_LIBS_USER /opt/R/user_libs
ENV BIOC_VERSION 3.14

RUN mkdir -p /home/ruser/R/libs \
    && mkdir -p /home/ruser/R/user_libs \
    && mkdir -p /root/R_groundhog \
    && mkdir -p /opt/R/user_libs

# Install R from source
RUN wget https://cloud.r-project.org/src/base/R-4/R-4.1.0.tar.gz \
    && tar -zxvf R-4.1.0.tar.gz \
    && rm R-4.1.0.tar.gz \
    && cd R-4.1.0 \
    && ./configure --enable-R-shlib \
    && make \
    && make install \
    && cd .. \
    && rm -rf R-4.1.0

# Install R packages
RUN R -e "install.packages(c('data.table', 'optparse'), repos = c(CRAN = 'https://cloud.r-project.org'))"

# Install packages in the root env
RUN pip install -y pandas==1.5.3 statsmodels==0.13.5

# Install miniconda and create a python2 environment
WORKDIR /opt
ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda
ENV PATH=$CONDA_DIR/bin:$PATH
RUN conda create -n python2 "python<3" 

# Use python2 to install qctool
WORKDIR /opt
RUN wget https://code.enkre.net/qctool/zip/release/qctool.tgz && \
    unzip qctool.tgz && \
    rm qctool.tgz && \
    cd qctool && \
    conda run -n python2 ./waf configure && \
    conda run -n python2 ./waf && \
    mv build/release/apps/qctool_v2.2.0 /usr/bin/qctool && \
    cd .. && \
    rm -rf qctool