#!/usr/bin/env bash
# Stage 07 -- DNABERT predictions from a local checkout (optional route).
# Required: DATA_DIRS MODEL_ROOTS DNABERT_EXAMPLES
# Optional: KMER MAXLEN_SPLICE MAXLEN_TFBS CONDA_ENV BATCH_SIZE
source "${PGVIS_BASH_LIB:-$(dirname "${BASH_SOURCE[0]}")}/lib/common.sh"
pgvis_require_var DATA_DIRS MODEL_ROOTS DNABERT_EXAMPLES
KMER="${KMER:-6}"; BATCH_SIZE="${BATCH_SIZE:-256}"
MAXLEN_SPLICE="${MAXLEN_SPLICE:-80}"; MAXLEN_TFBS="${MAXLEN_TFBS:-300}"
[ -n "${CONDA_ENV:-}" ] && { source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate "$CONDA_ENV"; }

# Fine-tuned models are saved as nested checkpoint-N / logstepN directories;
# take the highest-numbered one.
find_ckpt() { find "$1" -type d \( -name 'checkpoint-*' -o -name 'logstep*' \) 2>/dev/null | sort -V | tail -n1; }

cd "$DNABERT_EXAMPLES" || pgvis_die "cannot cd to $DNABERT_EXAMPLES"

read -r -a DATA_ARR  <<< "$DATA_DIRS"
read -r -a ROOT_ARR  <<< "$MODEL_ROOTS"

for data_dir in "${DATA_ARR[@]}"; do
  feature=$(basename "$data_dir")
  [ -f "$data_dir/dev.tsv" ] || { pgvis_log "no dev.tsv in $data_dir"; continue; }
  npy="$data_dir/preds/pred_results.npy"
  pgvis_done "$npy" && continue

  model=""
  for root in "${ROOT_ARR[@]}"; do
    [ -d "$root/$feature" ] && { model=$(find_ckpt "$root/$feature"); break; }
  done
  [ -n "$model" ] || { pgvis_log "no fine-tuned model found for '$feature' -- skipping"; continue; }

  case "$feature" in
    acceptor|donor) maxlen="$MAXLEN_SPLICE" ;;
    *)              maxlen="$MAXLEN_TFBS" ;;
  esac

  mkdir -p "$data_dir/preds"
  pgvis_log "predicting $feature (maxlen $maxlen) with $model"
  python run_finetune.py \
    --model_type dna --tokenizer_name "dna${KMER}" \
    --model_name_or_path "$model" --task_name dnaprom --do_predict \
    --data_dir "$data_dir" --max_seq_length "$maxlen" \
    --per_gpu_pred_batch_size="$BATCH_SIZE" \
    --output_dir "$model" --predict_dir "$data_dir/preds" \
    --n_process "${THREADS:-8}" \
    || pgvis_log "prediction FAILED for $feature"
done
pgvis_log "stage 07 complete"
