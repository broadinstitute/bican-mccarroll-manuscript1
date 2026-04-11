#!/usr/bin/env bash
# MIT License
#
# Copyright 2026 Broad Institute
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

set -x
set -euo pipefail

# Set noninteractive mode for apt-get to avoid prompts during package installation.

export DEBIAN_FRONTEND=noninteractive

# Load the DISTRIB_CODENAME and DISTRIB_RELEASE variables.

source /etc/lsb-release

# Install dependencies needed to update apt sources.

apt-get -qq update
apt-get -qq install wget

# Setup apt sources for R and R version info for r-base-core.

wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc >> /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
echo "deb https://cloud.r-project.org/bin/linux/ubuntu ${DISTRIB_CODENAME}-cran40/" >> /etc/apt/sources.list

r_version=4.5.3
r_base_version=${r_version}-1.${DISTRIB_RELEASE/./}.0

# Install other dependencies.

apt-get -qq update
apt-get -qq install --no-install-recommends \
    automake \
    bcftools \
    build-essential \
    cmake \
    curl \
    file \
    gfortran \
    git \
    libblas-dev \
    libbz2-dev \
    libcurl4-openssl-dev \
    libgeos-dev \
    libglpk40 \
    libgsl-dev \
    libhdf5-dev \
    libjpeg-dev \
    liblapack-dev \
    liblzma-dev \
    libopenblas-dev \
    libpng-dev \
    libpq-dev \
    libssl-dev \
    libtool \
    libwebp-dev \
    libxml2-dev \
    libz-dev \
    pkg-config \
    tk \
    zlib1g-dev \
    r-base-core="${r_base_version}"

# Cleanup.
rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
