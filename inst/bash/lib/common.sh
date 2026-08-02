# PGViS shared shell helpers. Sourced by every stage.
# Nothing here is site-specific; configuration arrives as environment
# variables from the block PGViS::pgvis_render_script() prepends.

pgvis_log()  { printf '[PGViS %s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }
pgvis_die()  { pgvis_log "ERROR: $*"; exit 1; }

pgvis_need() {
  for t in "$@"; do
    command -v "$t" >/dev/null 2>&1 || pgvis_die "required tool not on PATH: $t"
  done
}

pgvis_require_var() {
  for v in "$@"; do
    [ -n "${!v:-}" ] || pgvis_die "required variable $v is unset"
  done
}

# Skip work that is already done, so every stage is restartable and a job
# killed by a walltime limit can simply be resubmitted.
# Set PGVIS_FORCE=1 to redo everything.
pgvis_done() {
  local out="$1"
  if [ "${PGVIS_FORCE:-0}" = "1" ]; then return 1; fi
  if [ -s "$out" ]; then pgvis_log "skip (exists, non-empty): $out"; return 0; fi
  if [ -f "$out" ]; then pgvis_log "removing empty placeholder: $out"; rm -f "$out"; fi
  return 1
}

# Write to a temporary file and move into place only on success, so an
# interrupted job never leaves a truncated output that a later run would skip.
pgvis_atomic() {
  local out="$1"; shift
  local tmp; tmp="$(mktemp "${out}.XXXXXX")"
  if "$@" > "$tmp"; then mv "$tmp" "$out"; else rm -f "$tmp"; return 1; fi
}

# Resolve this task's item from a list, for SLURM array jobs. Falls back to
# looping over everything when not running under an array.
pgvis_array_slice() {
  local -n _items=$1
  if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
    local i=$(( SLURM_ARRAY_TASK_ID - ${PGVIS_ARRAY_BASE:-1} ))
    if [ "$i" -lt 0 ] || [ "$i" -ge "${#_items[@]}" ]; then
      pgvis_log "array slot $SLURM_ARRAY_TASK_ID has no work; exiting"
      exit 0
    fi
    printf '%s\n' "${_items[$i]}"
  else
    printf '%s\n' "${_items[@]}"
  fi
}
