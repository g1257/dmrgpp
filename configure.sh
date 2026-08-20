#!/usr/bin/bash

set -euo pipefail

computeProcsForMake()
{
    local nprocs

    if [[ $(uname -s) == Darwin ]]; then
        nprocs=$(sysctl -n hw.ncpu)
    else
        nprocs=$(nproc)
    fi

    if (( nprocs > 3 )); then
        echo $((nprocs - 2))
    else
        echo 1
    fi
}

source_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
build_dir="${source_dir}/build"

mkdir -p "${build_dir}"

cmake_options=(
    -DCMAKE_BUILD_TYPE=Release
    -DCMAKE_COMPILE_WARNING_AS_ERROR=ON
)

if [[ $(uname -s) != Darwin ]]; then
    cmake_options+=(
        -DBUILD_SHARED_LIBS=OFF
    )
fi

cmake -S "${source_dir}" -B "${build_dir}" "${cmake_options[@]}"
cmake --build "${build_dir}" -j "$(computeProcsForMake)"
