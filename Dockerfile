FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive

LABEL author="Edoardo Giacopuzzi"
LABEL contact="edoardo.giacopuzzi@fht.org"
LABEL image_version="0.2"
LABEL regenie_version="3.1.3"
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
    zlib1g-dev \
    libgomp1 \
    procps \
    libx11-6 \
    openjdk-11-jdk \
    python3 \
    libjpeg-dev \
    libpcre2-dev pcre2-utils \
    libxml2-dev \
    autoconf ca-certificates wget libbz2-dev \
    libc6-dev libcurl4-openssl-dev libfreetype6 libgsl-dev liblzma-dev \
    libncurses5-dev libperl-dev libreadline-dev libssl-dev libz-dev \
&& ln -s /usr/bin/python3 /usr/bin/python \
&& wget http://security.ubuntu.com/ubuntu/pool/main/i/icu/libicu60_60.2-3ubuntu3.2_amd64.deb \
&& apt-get install ./libicu60_60.2-3ubuntu3.2_amd64.deb \
&& apt-get -y autoremove \
&& apt-get -y clean all \
&& rm -rf /var/cache 

# Install regenie (not as conda package available)
WORKDIR "/opt"
RUN mkdir regenie && cd regenie && \
    wget https://github.com/rgcgithub/regenie/releases/download/v3.1.3/regenie_v3.1.3.gz_x86_64_Linux.zip && \
    unzip -q regenie_v3.*.gz_x86_64_Linux.zip && \
    rm regenie_v3.*.gz_x86_64_Linux.zip && \
    mv regenie_v3.*.gz_x86_64_Linux regenie && \
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

# Install jbang (not as conda package available)
WORKDIR "/opt"
RUN wget https://github.com/jbangdev/jbang/releases/download/v0.91.0/jbang-0.91.0.zip && \
    unzip -q jbang-*.zip && \
    mv jbang-0.91.0 jbang  && \
    rm jbang*.zip
ENV PATH="/opt/jbang/bin:${PATH}"

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
RUN R -e "install.packages(c('BiocManager', 'groundhog'), repos = c(CRAN = 'https://cloud.r-project.org'))"
COPY install_R_pkgs.R .
RUN Rscript install_R_pkgs.R
