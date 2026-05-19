#!/bin/bash
#SBATCH --job-name=scProtXpress
#SBATCH --output=logs/scprotxpress_%j.out
#SBATCH --error=logs/scprotxpress_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=120G
#SBATCH --partition=large

set -euo pipefail

# ------------------------------------------------------------
# Edit before submitting:
# ------------------------------------------------------------
export SCPROTXPRESS_DATA_ROOT="${SCPROTXPRESS_DATA_ROOT:-/scratch/USER/data}"
export SCPROTXPRESS_OUT_ROOT="${SCPROTXPRESS_OUT_ROOT:-/scratch/USER/results}"
export SCPROTXPRESS_RUN_NAME="${SCPROTXPRESS_RUN_NAME:-run_$(date +%Y%m%d_%H%M%S)}"

# Pin BLAS threads -- foreach workers also pin via RhpcBLASctl inside R.
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

MODALITY="${1:-CyTOF}"

case "$MODALITY" in
  CyTOF)
    DRIVER="$(Rscript -e 'cat(system.file("drivers", "cytof_driver.R", package = "scProtXpress"))')"
    ;;
  CITE-seq|CITEseq)
    DRIVER="$(Rscript -e 'cat(system.file("drivers", "citeseq_driver.R", package = "scProtXpress"))')"
    ;;
  *)
    echo "Unknown modality: $MODALITY (expected CyTOF or CITE-seq)" >&2
    exit 1
    ;;
esac

mkdir -p logs
echo "[$(date -Iseconds)] starting scProtXpress driver: $DRIVER"
echo "[$(date -Iseconds)] data_root=$SCPROTXPRESS_DATA_ROOT"
echo "[$(date -Iseconds)] out_root=$SCPROTXPRESS_OUT_ROOT"
echo "[$(date -Iseconds)] run_name=$SCPROTXPRESS_RUN_NAME"
echo "[$(date -Iseconds)] cpus_per_task=$SLURM_CPUS_PER_TASK"

Rscript "$DRIVER"
echo "[$(date -Iseconds)] driver finished"
