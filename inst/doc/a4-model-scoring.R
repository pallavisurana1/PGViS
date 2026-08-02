## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(PGViS)

## ----repo---------------------------------------------------------------------
pgvis_model_repo()

## ----python-check, eval = FALSE-----------------------------------------------
#  # Is the backend available? (needs reticulate + a Python with deepvregulome)
#  pgvis_python_available()

## ----python-run, eval = FALSE-------------------------------------------------
#  preds <- pgvis_predict_python(
#    input_dirs = dirs,            # from pgvis_write_model_input()
#    outdir     = cfg$dirs$preds,
#    python     = cfg$python       # NULL = resolve at call time
#  )

## ----shell-run, eval = FALSE--------------------------------------------------
#  pgvis_predict_shell(
#    input_dirs       = dirs,
#    model_roots      = c("/models/splice", "/models/tfbs"),
#    dnabert_examples = "/opt/DNABERT/examples",
#    outdir           = cfg$dirs$preds,
#    submit           = TRUE
#  )

## ----byo----------------------------------------------------------------------
preds <- pgvis_example("predictions")
head(preds)

## ----read, eval = FALSE-------------------------------------------------------
#  preds <- pgvis_read_predictions(cfg$dirs$preds)

## ----qc-----------------------------------------------------------------------
qc <- pgvis_qc(pgvis_example("predictions"), verbose = FALSE)
qc$degenerate_predictions
qc$median_prob_shift

