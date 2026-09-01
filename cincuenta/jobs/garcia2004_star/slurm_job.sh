#!/usr/bin/env bash
#SBATCH --job-name=cincuenta-star
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --output=cincuenta-star-%j.out
#SBATCH --error=cincuenta-star-%j.err

set -euo pipefail

# Override these with sbatch --export or edit a private cluster-specific copy.
: "${CASE:=nb7_u2p5}"
: "${RANKS:=${SLURM_NTASKS:-16}}"

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
exec "$script_dir/run_case.sh" "$CASE" "$RANKS"
