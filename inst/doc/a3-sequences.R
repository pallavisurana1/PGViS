## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(PGViS)

## ----mutate-------------------------------------------------------------------
mutate_sequence("A", "G", "TTTATTT", pos_start = 4, pos_end = 4)   # exact hit
mutate_sequence("A", "G", "TTTATTT", pos_start = 5, pos_end = 5)   # nominal pos off by one
mutate_sequence("CCCC", "G", "TTTATTT", pos_start = 4, pos_end = 7) # allele absent -> NA

## ----revcomp------------------------------------------------------------------
revcomp("ACCGGT")
revcomp(revcomp("ACGTTGCA"))    # round-trips

## ----build, eval = FALSE------------------------------------------------------
#  variants <- pgvis_example("variants")
#  seqs <- pgvis_build_sequences(
#    variants,
#    max_allele_len = 10,       # drop long indels the models weren't trained on
#    drop_failures  = TRUE      # remove rows whose alt sequence couldn't be built
#  )
#  seqs[, c("feature", "ref_kmer", "alt_kmer")]

## ----kmer---------------------------------------------------------------------
seq2kmer("ACGTACGTAC", k = 6)

## ----write, eval = FALSE------------------------------------------------------
#  dirs <- pgvis_write_model_input(seqs, cfg$dirs$seq)
#  dirs
#  #> [PGViS] Wrote model input for 4 feature(s) to .../07_seq

