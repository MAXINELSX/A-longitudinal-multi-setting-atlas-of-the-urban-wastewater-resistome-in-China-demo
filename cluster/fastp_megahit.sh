set -euo pipefail
set -x

echo "[BOOT] $(date) host=$(hostname) job=$SLURM_JOB_ID task=${SLURM_ARRAY_TASK_ID:-NA} pwd=$(pwd)"

source /appsnew/usr/python/Anaconda3-2025.06-1/etc/profile.d/conda.sh
conda activate qc_asm

TAG="${TAG:-J007_01}"
DATA_DIR="${DATA_DIR:-/home/tangs/tangs_lisx/lustre1/usb_copy/usb3/meishiqiang1113076/X101SC25052585-Z02-J008/00.CleanData}"

USB_ROOT="${USB_ROOT:-/home/tangs/tangs_lisx/lustre1/usb_copy}"
WORK_ROOT="${WORK_ROOT:-${USB_ROOT}/megahit}"
FINAL_ROOT="${FINAL_ROOT:-${USB_ROOT}/ass}"
mkdir -p "$WORK_ROOT" "$FINAL_ROOT/logs"

SAMPLE_LIST="${WORK_ROOT}/sample_list.${TAG}.txt"

if [[ "${SLURM_ARRAY_TASK_ID:-0}" -eq 0 ]]; then
  if [[ ! -s "$SAMPLE_LIST" ]]; then
    find "$DATA_DIR" -mindepth 1 -maxdepth 1 -type d -printf "%f\n" | sort | awk 'NF{print}' > "$SAMPLE_LIST"
  fi
fi

for i in {1..60}; do
  [[ -s "$SAMPLE_LIST" ]] && break
  sleep 1
done
[[ -s "$SAMPLE_LIST" ]] || { echo "ERROR: sample_list not ready: $SAMPLE_LIST"; exit 10; }

echo "[INFO] TAG=$TAG"
echo "[INFO] DATA_DIR=$DATA_DIR"
echo "[INFO] SAMPLE_LIST=$SAMPLE_LIST"
echo "[INFO] total samples: $(grep -c . "$SAMPLE_LIST")"
echo "[INFO] SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-NA}"
echo "[INFO] PATH=$PATH"
which fastp
which megahit

export TAG DATA_DIR USB_ROOT WORK_ROOT FINAL_ROOT SAMPLE_LIST
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
bash "$SCRIPT_DIR/pipeline_one_sample.sh"
