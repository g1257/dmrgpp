#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: run_case.sh CASE [RANKS]

CASE is one of: nb7_u2p5, nb11_u2p0, nb11_u2p5, nb11_u3p5, nb15_u2p5

Environment:
  CINCUENTA_EXE  MPI-enabled cincuenta executable (default: <repo>/build/cincuenta/src/cincuenta)
  RUN_ROOT       parent output directory (default: <repo>/runs/garcia2004_star)
  RUN_NAME       output directory basename (default: CASE_UTC-TIMESTAMP)
  MPIEXEC        MPI launcher outside Slurm (default: mpiexec)
  MPIEXEC_NFLAG  rank-count flag outside Slurm (default: -n)
  EXTRA_ARGS     extra whitespace-separated cincuenta arguments

Inside a Slurm allocation the script uses srun. Outside Slurm it uses mpiexec.
Every invocation creates a new directory and refuses to overwrite an old one.
EOF
}

if [[ $# -lt 1 || $# -gt 2 || $1 == -h || $1 == --help ]]; then
    usage
    [[ $# -ge 1 && ( $1 == -h || $1 == --help ) ]] && exit 0
    exit 2
fi

case_name=${1#input_}
case_name=${case_name%.ain}
case "$case_name" in
    nb7_u2p5)  default_ranks=16 ;;
    nb11_u2p0|nb11_u2p5|nb11_u3p5) default_ranks=32 ;;
    nb15_u2p5) default_ranks=64 ;;
    *) echo "Unknown case: $case_name" >&2; usage >&2; exit 2 ;;
esac
ranks=${2:-$default_ranks}
[[ $ranks =~ ^[1-9][0-9]*$ ]] || { echo "RANKS must be a positive integer" >&2; exit 2; }

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(git -C "$script_dir" rev-parse --show-toplevel)
input_source="$script_dir/inputs/input_${case_name}.ain"
exe=${CINCUENTA_EXE:-$repo_root/build/cincuenta/src/cincuenta}
run_root=${RUN_ROOT:-$repo_root/runs/garcia2004_star}
run_name=${RUN_NAME:-${case_name}_$(date -u +%Y%m%dT%H%M%SZ)}
run_dir="$run_root/$run_name"

[[ -f $input_source ]] || { echo "Missing input: $input_source" >&2; exit 1; }
[[ -x $exe ]] || {
    echo "Cincuenta executable is not executable: $exe" >&2
    echo "Set CINCUENTA_EXE to the MPI-enabled cluster build." >&2
    exit 1
}
mkdir -p "$run_root"
if ! mkdir "$run_dir"; then
    echo "Run directory already exists; refusing to overwrite: $run_dir" >&2
    exit 1
fi
cp "$input_source" "$run_dir/input.ain"

if [[ -n ${SLURM_JOB_ID:-} ]]; then
    launcher=(srun --ntasks="$ranks")
else
    launcher=("${MPIEXEC:-mpiexec}" "${MPIEXEC_NFLAG:--n}" "$ranks")
fi
extra_args=()
if [[ -n ${EXTRA_ARGS:-} ]]; then
    # EXTRA_ARGS is intentionally simple; arguments containing spaces are unsupported.
    read -r -a extra_args <<<"$EXTRA_ARGS"
fi
command=("${launcher[@]}" "$exe" -f input.ain -p 12 -l cincuenta.log "${extra_args[@]}")

export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
export MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}

{
    echo "case=$case_name"
    echo "ranks=$ranks"
    echo "utc_start=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "hostname=$(hostname)"
    echo "repo_root=$repo_root"
    echo "git_commit=$(git -C "$repo_root" rev-parse HEAD)"
    echo "git_describe=$(git -C "$repo_root" describe --always --dirty)"
    echo "executable=$exe"
    sha256sum "$exe" "$input_source"
    printf 'command='
    printf '%q ' "${command[@]}"
    printf '\n'
    env | grep -E '^(SLURM|SPACK|OMP|MKL|OPENBLAS|CUDA|ROCR|HIP)_' | sort || true
} >"$run_dir/manifest.txt"

printf 'Run directory: %s\n' "$run_dir"
printf 'Launching: '
printf '%q ' "${command[@]}"
printf '\n'

start_epoch=$(date +%s)
set +e
(
    cd "$run_dir"
    "${command[@]}"
) 2>&1 | tee "$run_dir/launcher.log"
status=${PIPESTATUS[0]}
set -e
end_epoch=$(date +%s)
{
    echo "utc_end=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "elapsed_seconds=$((end_epoch - start_epoch))"
    echo "exit_status=$status"
} >>"$run_dir/manifest.txt"

if (( status != 0 )); then
    echo "Cincuenta failed with status $status; outputs retained in $run_dir" >&2
    exit "$status"
fi

echo "Completed successfully: $run_dir"
