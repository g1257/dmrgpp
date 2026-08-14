#!/usr/bin/bash

set -euo pipefail

source_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
build_dir="${source_dir}/build"

mkdir -p "${build_dir}"

cmake_options=(
    -DCMAKE_BUILD_TYPE=Release
    -DCMAKE_COMPILE_WARNING_AS_ERROR=ON
)

if [[ $(uname -s) == Darwin ]]; then
    cmake_options+=(
        -DKokkosKernels_ENABLE_TPL_BLAS=ON
        -DKokkosKernels_ENABLE_TPL_LAPACK=ON
    )
else
    cmake_options+=(
        -DBUILD_SHARED_LIBS=OFF
    )
fi

cmake -S "${source_dir}" -B "${build_dir}" "${cmake_options[@]}"
