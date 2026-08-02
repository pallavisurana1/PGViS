# ---------------------------------------------------------------------------
# Model scoring bridge.
#
# PGViS does not run neural networks itself. It prepares model input, hands
# off to DNABERT / DeepVRegulome, and reads the predictions back. There are
# three supported routes, in decreasing order of convenience:
#
#   1. `pgvis_predict_python()`   - the `deepvregulome` PyPI package, via
#                                   reticulate. Handles model download from
#                                   Hugging Face itself.
#   2. `pgvis_predict_shell()`    - a local DNABERT
#                                   `run_finetune.py --do_predict` loop,
#                                   shipped as inst/bash/07_predict.sh.
#   3. bring your own             - write a table with the right columns and
#                                   pass it to `pgvis_variant_metrics()`.
#
# All three converge on the same tidy prediction table, so the rest of the
# package does not care which you used. Every function here is OPTIONAL:
# nothing in the metric or scoring layers requires Python.
# ---------------------------------------------------------------------------

#' Hugging Face repository holding the fine-tuned models
#'
#' The fine-tuned DNABERT models used by PGViS are published at
#' <https://huggingface.co/duttaprat/DeepVRegulome>: donor and acceptor
#' splice-site classifiers plus several hundred TFBS/RBP/histone-mark models,
#' one per ENCODE target.
#'
#' @return The repository ID, as a string.
#' @export
pgvis_model_repo <- function() "duttaprat/DeepVRegulome"

#' Check whether the optional Python scoring stack is available
#'
#' @param python Optional path to a Python executable or conda environment.
#' @return `TRUE` if `reticulate` can find a Python with `deepvregulome`
#'   importable, otherwise `FALSE` with an informative message.
#' @export
pgvis_python_available <- function(python = NULL) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    .msg("reticulate is not installed; Python scoring unavailable.")
    return(FALSE)
  }
  if (!is.null(python)) {
    try(reticulate::use_python(python, required = TRUE), silent = TRUE)
  }
  ok <- reticulate::py_module_available("deepvregulome")
  if (!ok) {
    .msg("Python found but 'deepvregulome' is not importable. Install with:\n",
         "    pip install deepvregulome")
  }
  ok
}

#' Score variants with the DeepVRegulome Python package
#'
#' Optional. Calls the `deepvregulome` package through `reticulate`,
#' which pulls the fine-tuned models from Hugging Face on first use and
#' caches them locally.
#'
#' @param input_dirs Directories written by [pgvis_write_model_input()]. The
#'   directory name is taken as the model/feature name.
#' @param outdir Where to write per-feature prediction files.
#' @param python Optional Python executable or conda environment path.
#' @param repo Hugging Face repository; see [pgvis_model_repo()].
#' @param cache_dir Local model cache. Defaults to the `deepvregulome` default.
#' @param batch_size Prediction batch size.
#' @param device `"cuda"`, `"cpu"`, or `NULL` to autodetect.
#' @param overwrite Re-run features that already have predictions. The default
#'   skip-if-present behaviour makes a run restartable after a job times out.
#' @return A tidy prediction table; see [pgvis_read_predictions()].
#' @export
pgvis_predict_python <- function(input_dirs, outdir,
                                 python = NULL,
                                 repo = pgvis_model_repo(),
                                 cache_dir = NULL,
                                 batch_size = 256L,
                                 device = NULL,
                                 overwrite = FALSE) {
  if (!pgvis_python_available(python)) {
    stop("The Python scoring backend is not available. Either\n",
         "  * pip install deepvregulome (and install.packages('reticulate')), or\n",
         "  * use pgvis_predict_shell() with a local DNABERT checkout, or\n",
         "  * supply your own predictions to pgvis_variant_metrics().",
         call. = FALSE)
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  script <- system.file("python", "run_deepvregulome.py", package = "PGViS")
  runner <- reticulate::import_from_path(
    "run_deepvregulome", path = dirname(script), convert = TRUE)

  runner$predict_dirs(
    input_dirs = as.list(input_dirs),
    outdir     = outdir,
    repo       = repo,
    cache_dir  = cache_dir,
    batch_size = as.integer(batch_size),
    device     = device,
    overwrite  = overwrite
  )
  pgvis_read_predictions(outdir)
}

#' Score variants with a local DNABERT checkout
#'
#' Optional. Renders and runs `inst/bash/07_predict.sh`, which is the
#' rendered from a template: it locates the deepest `checkpoint-*` /
#' `logstep*` directory under each model root, runs
#' `run_finetune.py --do_predict`, and skips features whose
#' `preds/pred_results.npy` already exists.
#'
#' @param input_dirs Directories from [pgvis_write_model_input()].
#' @param model_roots One or more directories containing fine-tuned models,
#'   one subdirectory per feature.
#' @param dnabert_examples Path to the DNABERT `examples/` directory.
#' @param outdir Prediction output directory.
#' @param kmer K-mer size (must match the tokenizer the models were trained with).
#' @param max_seq_length Named vector giving the maximum sequence length per
#'   feature type. Splice models use 80, TFBS models 300.
#' @param conda_env Optional conda environment to activate.
#' @param batch_size Per-GPU batch size.
#' @param cfg Optional [pgvis_config()]; used for scheduler settings.
#' @param submit Run now (`FALSE`) or emit a submittable script (`TRUE`).
#' @return Path to the rendered script, invisibly.
#' @export
pgvis_predict_shell <- function(input_dirs, model_roots, dnabert_examples,
                                outdir,
                                kmer = 6L,
                                max_seq_length = c(splice = 80L, tfbs = 300L),
                                conda_env = NULL,
                                batch_size = 256L,
                                cfg = NULL,
                                submit = FALSE) {
  vars <- list(
    DATA_DIRS        = paste(shQuote(input_dirs), collapse = " "),
    MODEL_ROOTS      = paste(shQuote(model_roots), collapse = " "),
    DNABERT_EXAMPLES = dnabert_examples,
    OUT_DIR          = outdir,
    KMER             = kmer,
    MAXLEN_SPLICE    = max_seq_length[["splice"]],
    MAXLEN_TFBS      = max_seq_length[["tfbs"]],
    CONDA_ENV        = conda_env %||% "",
    BATCH_SIZE       = batch_size
  )
  pgvis_render_script("07_predict.sh", vars, cfg = cfg, submit = submit)
}

#' Read model predictions into a tidy table
#'
#' Understands two layouts:
#'
#' * `<dir>/<feature>/preds/pred_results.npy` alongside `variants.tsv` --
#'   the raw DNABERT output. The `.npy` is a flat array of probabilities in
#'   reference/alternate pairs, reshaped to two columns and bound back to the
#'   variant rows in order. Row counts must agree exactly; a mismatch means
#'   the input was regenerated after the prediction and is an error, not a
#'   warning.
#' * `<dir>/<feature>.tsv` or `<dir>/<feature>/predictions.tsv` with explicit
#'   `p_ref` / `p_alt` columns -- what `pgvis_predict_python()` writes.
#'
#' @param path Prediction directory.
#' @param features Optional subset of features to read.
#' @return A table with `feature`, `func_site`, `var_site`, `p_ref`, `p_alt`.
#' @export
pgvis_read_predictions <- function(path, features = NULL) {
  .need("data.table", "reading predictions")
  dirs <- list.dirs(path, recursive = FALSE, full.names = TRUE)
  if (!is.null(features)) dirs <- dirs[basename(dirs) %in% features]
  if (!length(dirs)) stop("No prediction directories found under ", path, call. = FALSE)

  out <- lapply(dirs, function(d) {
    feature <- basename(d)
    tsv <- file.path(d, "predictions.tsv")
    if (file.exists(tsv)) {
      df <- data.table::fread(tsv, data.table = FALSE, showProgress = FALSE)
      df$feature <- feature
      return(.pgvis_finalise_preds(df, feature))
    }
    npy  <- file.path(d, "preds", "pred_results.npy")
    meta <- file.path(d, "variants.tsv")
    if (!file.exists(npy) || !file.exists(meta)) {
      .msg("Skipping ", feature, ": no predictions found.")
      return(NULL)
    }
    .need("reticulate", "reading .npy prediction arrays")
    np <- reticulate::import("numpy", convert = TRUE)
    arr <- matrix(as.numeric(np$load(npy)), ncol = 2L, byrow = TRUE)
    df <- data.table::fread(meta, data.table = FALSE, showProgress = FALSE)
    if (nrow(arr) != nrow(df)) {
      stop(sprintf(
        "%s: %d prediction pairs but %d input variants. The model input was regenerated after prediction; re-run the model for this feature.",
        feature, nrow(arr), nrow(df)), call. = FALSE)
    }
    df$p_ref <- arr[, 1L]
    df$p_alt <- arr[, 2L]
    .pgvis_finalise_preds(df, feature)
  })

  res <- data.table::rbindlist(Filter(Negate(is.null), out),
                               use.names = TRUE, fill = TRUE)
  .msg(sprintf("Read %d predictions across %d feature(s).",
               nrow(res), length(unique(res$feature))))
  res[]
}

.pgvis_finalise_preds <- function(df, feature) {
  df$feature <- feature
  keep <- intersect(c("feature", "func_site", "var_site", "p_ref", "p_alt"),
                    names(df))
  miss <- setdiff(c("func_site", "var_site", "p_ref", "p_alt"), keep)
  if (length(miss)) {
    stop(sprintf("Feature '%s': prediction table is missing %s.",
                 feature, paste(miss, collapse = ", ")), call. = FALSE)
  }
  unique(df[, keep, drop = FALSE])
}
