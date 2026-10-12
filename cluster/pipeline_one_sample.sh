set -euo pipefail
shopt -s nullglob

TAG="${TAG:-J007_01}"

DATA_DIR="${DATA_DIR:-/home/tangs/tangs_lisx/lustre1/usb_copy/usb4/yaoqiao1103036/X101SC25052585-Z02-J007_multipath/X101SC25052585-Z02-J007_01/00.CleanData}"

USB_ROOT="${USB_ROOT:-/home/tangs/tangs_lisx/lustre1/usb_copy}"
WORK_ROOT="${WORK_ROOT:-${USB_ROOT}/megahit}"
FINAL_ROOT="${FINAL_ROOT:-${USB_ROOT}/ass}"

SAMPLE_LIST="${SAMPLE_LIST:-${WORK_ROOT}/sample_list.${TAG}.txt}"

CPUS="${SLURM_CPUS_PER_TASK:-64}"
QUAL=20
LEN=50
MIN_CONTIG_LEN=500

TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
LINE=$((TASK_ID + 1))
SAMPLE_NAME="$(sed -n "${LINE}p" "$SAMPLE_LIST" | tr -d '\r\n')"
[[ -n "$SAMPLE_NAME" ]] || { echo "ERROR: empty sample for task $TASK_ID (line $LINE)"; exit 1; }

SAMPLE_DIR="${DATA_DIR}/${SAMPLE_NAME}"
[[ -d "$SAMPLE_DIR" ]] || { echo "ERROR: sample dir not found: $SAMPLE_DIR"; exit 1; }

command -v fastp >/dev/null || { echo "ERROR: fastp not in PATH (check conda env)"; exit 1; }
command -v megahit >/dev/null || { echo "ERROR: megahit not in PATH (check conda env)"; exit 1; }

r1_arr=("${SAMPLE_DIR}"/*_1.clean.fq.gz)
r2_arr=("${SAMPLE_DIR}"/*_2.clean.fq.gz)
if [[ ${#r1_arr[@]} -ne 1 || ${#r2_arr[@]} -ne 1 ]]; then
  echo "ERROR: cannot uniquely find *_1.clean.fq.gz/*_2.clean.fq.gz in $SAMPLE_DIR"
  echo "R1=${#r1_arr[@]} R2=${#r2_arr[@]}"
  exit 1
fi
IN_R1="${r1_arr[0]}"
IN_R2="${r2_arr[0]}"

prefix="$(basename "$IN_R1")"
prefix="${prefix%_1.clean.fq.gz}"

WORK_SAMPLE="${WORK_ROOT}/${TAG}/${SAMPLE_NAME}"
CLEAN2_DIR="${WORK_SAMPLE}/clean2"
ASM_DIR="${WORK_SAMPLE}/megahit_out"
TMP_DIR="${WORK_SAMPLE}/tmp"
mkdir -p "$CLEAN2_DIR"

MH_LOG="${WORK_SAMPLE}/${SAMPLE_NAME}.megahit.log"

FINAL_SAMPLE="${FINAL_ROOT}/${TAG}/${SAMPLE_NAME}"
FINAL_FA="${FINAL_SAMPLE}/final.contigs.fa"
mkdir -p "$FINAL_SAMPLE"

OUT_R1="${CLEAN2_DIR}/${prefix}_1.clean2.fq.gz"
OUT_R2="${CLEAN2_DIR}/${prefix}_2.clean2.fq.gz"
FP_HTML="${CLEAN2_DIR}/${prefix}.fastp.html"
FP_JSON="${CLEAN2_DIR}/${prefix}.fastp.json"
FP_LOG="${CLEAN2_DIR}/${prefix}.fastp.log"

echo "[INFO] TAG=$TAG SAMPLE=$SAMPLE_NAME CPUS=$CPUS"
echo "[INFO] IN_R1=$IN_R1"
echo "[INFO] IN_R2=$IN_R2"
echo "[INFO] ASM_DIR=$ASM_DIR"
echo "[INFO] FINAL_FA=$FINAL_FA"
echo

if [[ -s "$FINAL_FA" ]]; then
  echo "[SKIP] final exists: $FINAL_FA"
  exit 0
fi

if [[ -s "$ASM_DIR/final.contigs.fa" ]]; then
  echo "[INFO] found existing final in WORK -> copy to ass"
  cp -f "$ASM_DIR/final.contigs.fa" "$FINAL_FA"
  ( cd "$FINAL_SAMPLE" && md5sum final.contigs.fa > final.contigs.fa.md5 )
  exit 0
fi

if [[ -s "$OUT_R1" && -s "$OUT_R2" ]]; then
  echo "[INFO] clean2 exists, skip fastp"
else
  echo "[RUN ] fastp (2nd round)"
  rm -f "$OUT_R1" "$OUT_R2" "$FP_HTML" "$FP_JSON" "$FP_LOG"
  fastp \
    -i "$IN_R1" -I "$IN_R2" \
    -o "$OUT_R1" -O "$OUT_R2" \
    -w "$CPUS" \
    --detect_adapter_for_pe \
    --cut_front --cut_tail \
    --cut_window_size 4 --cut_mean_quality "$QUAL" \
    --length_required "$LEN" \
    -h "$FP_HTML" -j "$FP_JSON" \
    > "$FP_LOG" 2>&1
fi

run_megahit_continue () {
  rm -rf "$TMP_DIR"
  mkdir -p "$TMP_DIR"
  megahit --continue -o "$ASM_DIR" -t "$CPUS" --tmp-dir "$TMP_DIR" >>"$MH_LOG" 2>&1
}

run_megahit_fresh () {
  rm -rf "$ASM_DIR" "$TMP_DIR"
  mkdir -p "$TMP_DIR"
  megahit \
    -1 "$OUT_R1" -2 "$OUT_R2" \
    -o "$ASM_DIR" \
    -t "$CPUS" \
    --min-contig-len "$MIN_CONTIG_LEN" \
    --tmp-dir "$TMP_DIR" >>"$MH_LOG" 2>&1
}

{
  echo
  echo "========== MEGAHIT PIPELINE =========="
  echo "time   : $(date -Is)"
  echo "sample : $SAMPLE_NAME"
  echo "ASM_DIR: $ASM_DIR"
  echo "======================================"
} >> "$MH_LOG"

if [[ -d "$ASM_DIR" && ! -s "$ASM_DIR/final.contigs.fa" ]]; then
  echo "[RUN ] megahit --continue" | tee -a "$MH_LOG"
  if run_megahit_continue; then
    if [[ -s "$ASM_DIR/final.contigs.fa" ]]; then
      echo "[OK  ] continue produced final" | tee -a "$MH_LOG"
    else
      echo "[WARN] continue done but final missing -> fresh" | tee -a "$MH_LOG"
      run_megahit_fresh || { echo "[WARN] fresh failed -> retry fresh once" | tee -a "$MH_LOG"; run_megahit_fresh; }
    fi
  else
    echo "[WARN] continue failed -> fresh" | tee -a "$MH_LOG"
    run_megahit_fresh || { echo "[WARN] fresh failed -> retry fresh once" | tee -a "$MH_LOG"; run_megahit_fresh; }
  fi
else
  echo "[RUN ] megahit fresh" | tee -a "$MH_LOG"
  run_megahit_fresh || { echo "[WARN] fresh failed -> retry fresh once" | tee -a "$MH_LOG"; run_megahit_fresh; }
fi

if [[ ! -s "$ASM_DIR/final.contigs.fa" ]]; then
  echo "ERROR: final.contigs.fa missing after megahit" | tee -a "$MH_LOG"
  exit 2
fi

cp -f "$ASM_DIR/final.contigs.fa" "$FINAL_FA"
cp -f "$MH_LOG" "$FINAL_SAMPLE/" || true
cp -f "$FP_HTML" "$FP_JSON" "$FP_LOG" "$FINAL_SAMPLE/" || true
( cd "$FINAL_SAMPLE" && md5sum final.contigs.fa > final.contigs.fa.md5 )

echo "[DONE] $SAMPLE_NAME"
