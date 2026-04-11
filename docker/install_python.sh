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

# Get the path to the top-level directory of the repository.

src_dir=$(realpath "$(dirname "${BASH_SOURCE[0]}")"/..)

# Install python packages.

uv_dir=/tmp/uv
envs_dir="$src_dir"/python
base_dir=/usr/local/bican-mccarroll-manuscript1
pythons_dir="$base_dir"/pythons
venv_dir=$base_dir/venv

export UV_PYTHON_INSTALL_DIR=$pythons_dir
export UV_TOOL_DIR=$venv_dir
export UV_TOOL_BIN_DIR=$base_dir
export PATH=$PATH:$uv_dir:$base_dir

## Install uv.
curl -LsSf https://astral.sh/uv/install.sh | env UV_UNMANAGED_INSTALL=$uv_dir sh

## Install the python environments from the local checkout.
for pyproject_toml in "$envs_dir"/*/pyproject.toml; do
  env_name=$(basename "$(dirname "$pyproject_toml")")
  env_dir="$envs_dir"/"$env_name"
  echo "Building $env_name"
  uv tool install "$env_dir"
done

# Cleanup.
rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
