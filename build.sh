#!/usr/bin/bash

set -euo pipefail

project_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
build_dir="${project_dir}/build"

mkdir -p "${build_dir}"

cmake -S "${project_dir}" -B "${build_dir}" \
    -DCMAKE_BUILD_TYPE=Debug \
    -DCMAKE_CXX_FLAGS="-Wall -Werror" \
    -DCMAKE_CXX_COMPILER=clang++
