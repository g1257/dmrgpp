#!/bin/bash

set -euo pipefail
set -x

step="${1:-}"
if [[ -z "${step}" ]]; then
  echo "Usage: $0 {configure|build|test|install}" >&2
  exit 2
fi

# This workflow runs on a self-hosted GPU VM.  Dependencies are expected to be
# provided by a prebuilt, host-local Spack environment rather than rebuilt for
# each CI run.  Override these paths with GitHub Actions variables or runner
# environment variables if the VM layout changes.
spack_setup="${DMRGPP_CUDA_SPACK_SETUP:-/scratch/spack/share/spack/setup-env.sh}"
spack_env="${DMRGPP_CUDA_SPACK_ENV:-/scratch/spack_envs/dmrgpp}"
cuda_root="${DMRGPP_CUDA_ROOT:-/usr/local/cuda-12.6}"

if [[ -r "${spack_setup}" ]]; then
  # Spack's shell setup may reference unset variables internally.
  set +u
  # shellcheck source=/dev/null
  source "${spack_setup}"
  set -u
elif ! command -v spack >/dev/null 2>&1; then
  echo "Spack setup script '${spack_setup}' is not readable and spack is not on PATH" >&2
  exit 1
fi

if [[ ! -d "${spack_env}" ]]; then
  echo "Spack environment '${spack_env}' does not exist" >&2
  exit 1
fi

if declare -F spack >/dev/null 2>&1; then
  set +u
  spack env activate -d "${spack_env}"
  set -u
else
  eval "$(spack env activate --sh "${spack_env}")"
fi

if [[ -x "${cuda_root}/bin/nvcc" ]]; then
  export CUDA_HOME="${cuda_root}"
  export CUDACXX="${cuda_root}/bin/nvcc"
  export PATH="${cuda_root}/bin:${PATH}"
fi

source_dir="${DMRGPP_SOURCE_DIR:-${GITHUB_WORKSPACE:-$(pwd)}/source}"
build_dir="${DMRGPP_BUILD_DIR:-${GITHUB_WORKSPACE:-$(pwd)}/build}"
install_dir="${DMRGPP_INSTALL_DIR:-${GITHUB_WORKSPACE:-$(pwd)}/install}"

case "${step}" in
  configure)
    nvidia-smi
    nvcc --version
    /usr/bin/g++-13 --version
    spack find

    if [[ -d "${build_dir}" ]]; then
      rm -rf "${build_dir}"
    fi

    cuda_root_args=()
    if [[ -n "${CUDA_HOME:-}" ]]; then
      cuda_root_args+=("-DCUDAToolkit_ROOT=${CUDA_HOME}")
    fi

    cmake -S "${source_dir}" \
          -B "${build_dir}" \
          -GNinja \
          -DCMAKE_BUILD_TYPE=RelWithDebInfo \
          -DCMAKE_CXX_COMPILER=/usr/bin/g++-13 \
          -DCMAKE_COMPILE_WARNING_AS_ERROR=ON \
          -DFETCHCONTENT_TRY_FIND_PACKAGE_MODE=NEVER \
          -DKokkos_ENABLE_CUDA=ON \
          -DKokkos_ARCH_AMPERE80=ON \
          "${cuda_root_args[@]}"
    ;;

  build)
    cmake --build "${build_dir}" --parallel "$(nproc)"
    ;;

  test)
    CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}" \
      ctest --test-dir "${build_dir}" -LE Nightly --output-on-failure
    ;;

  install)
    cmake --install "${build_dir}" --prefix "${install_dir}"
    ;;

  *)
    echo "Invalid step '${step}'" >&2
    exit 2
    ;;
esac
