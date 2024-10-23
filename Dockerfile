FROM rocker/tidyverse:4.4.0

LABEL maintainer = "Ammar Naqvi (naqvia@chop.edu)"

#########################################
RUN apt-get update && apt-get install -y --no-install-recommends apt-utils dialog

# Add curl, bzip2 and some dev libs
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    bedtools \
    cpanminus \
    curl \
    bzip2 \
    libbz2-dev \
    liblzma-dev \
    libreadline-dev \
    libgdal-dev \
    libudunits2-dev \
    libmagick++-dev \
    zlib1g

# Define environment variable for cleaner Dockerfile
ENV UCSC_USER_UTILS="http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64"
ENV PATH="/usr/local/ucsc-tools:${PATH}"

# Download bigBedToBed
RUN mkdir -p /usr/local/ucsc-tools && \
    wget "${UCSC_USER_UTILS}/bigBedToBed" -O /usr/local/ucsc-tools/bigBedToBed && \
    chmod +x /usr/local/ucsc-tools/bigBedToBed

# install perl packages
RUN cpanm install DBI
RUN cpanm install Statistics::Lite

# Set the Bioconductor repository as the primary repository
RUN R -e "options(repos = BiocManager::repositories())"

# Install BiocManager and the desired version of Bioconductor
RUN R -e "install.packages('BiocManager', dependencies=TRUE)"
RUN R -e "BiocManager::install(version = '3.19')"

# Install packages
RUN R -e 'BiocManager::install(c( \
  "biomaRt", \
  "BSgenome.Hsapiens.UCSC.hg38", \
  "data.table", \
  "DBI", \
  "GenomicRanges", \
  "ggthemes", \
  "optparse", \
  "RColorBrewer", \
  "RSQLite", \
  "R.utils", \
  "tidytext" \
))'
  
	
RUN R -e "remotes::install_github('clauswilke/colorblindr', ref = '1ac3d4d62dad047b68bb66c06cee927a4517d678', dependencies = TRUE)"
RUN R -e "remotes::install_github('thomasp85/patchwork', ref = '1cb732b129ed6a65774796dc1f618558c7498b66')"
	
# Install python and python packages
# Install pip3 and low-level python installation reqs
RUN apt-get -y --no-install-recommends install \
    python3-pip  python3-dev
RUN ln -s /usr/bin/python3 /usr/bin/python 
RUN python3 -m pip install --upgrade pip

RUN pip3 install \
    "rpg==2.0.3"

WORKDIR /rocker-build/

ADD Dockerfile . 

