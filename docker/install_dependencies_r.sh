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

# Based on / hat tip:
#   - https://github.com/broadinstitute/Drop-seq/blob/master/src/docker/R/common/install.sh
#   - https://cloud.r-project.org/bin/linux/ubuntu/
#   - https://github.com/rocker-org/rocker-versioned2/blob/77e4b1c/scripts/install_R_source.sh
#   - https://github.com/rocker-org/rocker-versioned2/blob/R4.5.2/scripts/setup_R.sh#L33-L40
#   - https://packagemanager.posit.co/client/#/repos/cran/setup?distribution=ubuntu-24.04&r_environment=other
#   - https://forum.posit.co/t/unable-to-install-binary-packages-from-packagemanager-rstudio-com-on-linux/82161
#   - https://docs.posit.co/rspm/admin/serving-binaries.html#binary-user-agents
set -x
set -euo pipefail

# Get the path to the top-level directory of the repository.

src_dir=$(realpath "$(dirname "${BASH_SOURCE[0]}")"/..)

# Load the DISTRIB_CODENAME and DISTRIB_RELEASE variables.

source /etc/lsb-release

# Install R dependencies.

## Get the R home directory location.
RHOME=$(R RHOME)

## Configure R to use the Posit Public Package Manager (P3M) for binary packages.
cat <<EOF >>"${RHOME}/etc/Rprofile.site"
options(repos = c(CRAN = sprintf("https://packagemanager.posit.co/cran/latest/bin/linux/$DISTRIB_CODENAME-%s/%s", R.version["arch"], substr(getRversion(), 1, 3))))
options(HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(), paste(getRversion(), R.version["platform"], R.version["arch"], R.version["os"])))
EOF

## Install the packages from the local checkout.
Rscript -e '
source("'"$src_dir"'/R/tools/install_packages.R");
install_all_packages(path="'"$src_dir"'", dependencies_only=TRUE);
'

# Cleanup.
rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
