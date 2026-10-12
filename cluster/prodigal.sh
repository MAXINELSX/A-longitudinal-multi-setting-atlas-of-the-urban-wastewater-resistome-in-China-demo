set -euo pipefail
set -x

echo "[BOOT] $(date) host=$(hostname) job=$SLURM_JOB_ID task=$SLURM_ARRAY_TASK_ID"

source /appsnew/usr/python/Anaconda3-2025.06-1/etc/profile.d/conda.sh
conda activate qc_asm
command -v prodigal >/dev/null || { echo "ERROR: prodigal not found"; exit 1; }

LISTFILE="${LISTFILE:?need LISTFILE via --export LISTFILE=/path/to/list.txt}"

OUTROOT="/home/tangs/tangs_lisx/lustre1/usb_copy/prodigal/meta_gff"
LOGDIR="/home/tangs/tangs_lisx/lustre1/usb_copy/prodigal/logs"
mkdir -p "$OUTROOT" "$LOGDIR"

CLOSED_OPT=()
if prodigal -h 2>&1 | grep -q -- '--closed-ends'; then
  CLOSED_OPT=(--closed-ends)
fi

LINE=$((SLURM_ARRAY_TASK_ID + 1))
FILE="$(sed -n "${LINE}p" "$LISTFILE" | sed 's/\r$//')"

if [[ -z "${FILE:-}" ]]; then
  echo "[SKIP] empty line at $LINE"
  exit 0
fi
if [[ ! -f "$FILE" ]]; then
  echo "[WARN] missing file: $FILE"
  exit 0
fi

BN="$(basename "$FILE")"
BN="${BN%.gz}"
BN="${BN%.fa}"
BN="${BN%.fna}"
BN="${BN%.fasta}"

OUTDIR="$OUTROOT/$BN"
mkdir -p "$OUTDIR"

if [[ -s "$OUTDIR/proteins.faa" && -s "$OUTDIR/genes.fna" && -s "$OUTDIR/genes.gff" ]]; then
  echo "[SKIP] done: $BN"
  exit 0
fi

echo "[RUN ] $BN <= $FILE"
prodigal -i "$FILE" -p meta "${CLOSED_OPT[@]}" -g 11 \
  -a "$OUTDIR/proteins.faa" \
  -d "$OUTDIR/genes.fna" \
  -f gff \
  -o "$OUTDIR/genes.gff"

echo "[DONE] $BN"
