"""Optional Python bridge from PGViS to the DeepVRegulome model collection.

PGViS calls this through reticulate. It is deliberately thin: PGViS owns
sequence construction, metric computation and scoring; this file owns only
"turn k-mer pairs into probabilities".

Requirements (none of which PGViS itself needs):

    pip install deepvregulome

Fine-tuned models are pulled from Hugging Face on first use and cached:
https://huggingface.co/duttaprat/DeepVRegulome

If the installed `deepvregulome` API differs from what is tried below, either
upgrade the package or fall back to `PGViS::pgvis_predict_shell()`, which
drives a local DNABERT checkout directly.
"""

from __future__ import annotations

import csv
import os
from typing import Iterable, Optional, Sequence

DEFAULT_REPO = "duttaprat/DeepVRegulome"


def _load_pairs(data_dir: str):
    """Read dev.tsv + variants.tsv written by PGViS::pgvis_write_model_input().

    dev.tsv alternates reference (label 0) and alternate (label 1) rows, so
    row 2*i of dev.tsv corresponds to row i of variants.tsv.
    """
    dev = os.path.join(data_dir, "dev.tsv")
    meta = os.path.join(data_dir, "variants.tsv")
    for p in (dev, meta):
        if not os.path.exists(p):
            raise FileNotFoundError(f"{p} not found; run pgvis_write_model_input() first")

    with open(dev, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    seqs = [r["Sequence"] for r in rows]

    with open(meta, newline="") as fh:
        variants = list(csv.DictReader(fh, delimiter="\t"))

    if len(seqs) != 2 * len(variants):
        raise ValueError(
            f"{data_dir}: dev.tsv has {len(seqs)} sequences but variants.tsv has "
            f"{len(variants)} rows (expected {2 * len(variants)}). "
            "The two files are out of sync; regenerate both."
        )
    return seqs, variants


def _get_predictor(feature: str, repo: str, cache_dir: Optional[str],
                   device: Optional[str]):
    """Resolve a callable that maps a list of k-mer strings to probabilities."""
    import deepvregulome as dvr  # noqa: F401  (import error is informative)

    kwargs = {}
    if cache_dir:
        kwargs["cache_dir"] = cache_dir
    if device:
        kwargs["device"] = device

    for attr in ("load_model", "from_pretrained", "Predictor", "DeepVRegulome"):
        fn = getattr(dvr, attr, None)
        if fn is None:
            continue
        try:
            return fn(feature, repo=repo, **kwargs)
        except TypeError:
            try:
                return fn(f"{repo}/{feature}", **kwargs)
            except Exception:
                continue
        except Exception:
            continue

    raise RuntimeError(
        f"Could not construct a DeepVRegulome predictor for '{feature}'. "
        "The installed deepvregulome API may have changed; use "
        "PGViS::pgvis_predict_shell() instead, or open an issue."
    )


def _probabilities(model, seqs: Sequence[str], batch_size: int):
    """Return P(positive class) for each sequence."""
    for attr in ("predict_proba", "predict", "__call__"):
        fn = getattr(model, attr, None)
        if fn is None:
            continue
        try:
            out = fn(list(seqs), batch_size=batch_size)
        except TypeError:
            out = fn(list(seqs))
        try:
            import numpy as np
            arr = np.asarray(out, dtype=float)
        except Exception:
            arr = out
        if hasattr(arr, "ndim") and arr.ndim == 2 and arr.shape[1] == 2:
            return arr[:, 1].tolist()
        return list(map(float, arr))
    raise RuntimeError("Predictor exposes no predict/predict_proba method.")


def predict_dirs(input_dirs: Iterable[str],
                 outdir: str,
                 repo: str = DEFAULT_REPO,
                 cache_dir: Optional[str] = None,
                 batch_size: int = 256,
                 device: Optional[str] = None,
                 overwrite: bool = False) -> list:
    """Score every PGViS model-input directory and write predictions.tsv.

    One output file per feature, with columns: func_site, var_site, p_ref,
    p_alt. Directories that already have a predictions.tsv are skipped unless
    ``overwrite`` is set, so a job killed partway through can be resubmitted.
    """
    os.makedirs(outdir, exist_ok=True)
    written = []

    for data_dir in input_dirs:
        feature = os.path.basename(os.path.normpath(data_dir))
        dest_dir = os.path.join(outdir, feature)
        dest = os.path.join(dest_dir, "predictions.tsv")

        if os.path.exists(dest) and not overwrite:
            print(f"[PGViS] skip {feature}: predictions.tsv exists")
            written.append(dest)
            continue

        try:
            seqs, variants = _load_pairs(data_dir)
            model = _get_predictor(feature, repo, cache_dir, device)
            probs = _probabilities(model, seqs, batch_size)
        except Exception as exc:  # keep going; one bad feature is not fatal
            print(f"[PGViS] FAILED {feature}: {exc}")
            continue

        if len(probs) != len(seqs):
            print(f"[PGViS] FAILED {feature}: got {len(probs)} probabilities "
                  f"for {len(seqs)} sequences")
            continue

        os.makedirs(dest_dir, exist_ok=True)
        tmp = dest + ".tmp"
        with open(tmp, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["func_site", "var_site", "p_ref", "p_alt"])
            for i, v in enumerate(variants):
                w.writerow([v.get("func_site", ""), v.get("var_site", ""),
                            probs[2 * i], probs[2 * i + 1]])
        os.replace(tmp, dest)
        print(f"[PGViS] {feature}: {len(variants)} variants -> {dest}")
        written.append(dest)

    return written


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description="Score PGViS model inputs with DeepVRegulome")
    ap.add_argument("input_dirs", nargs="+")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--repo", default=DEFAULT_REPO)
    ap.add_argument("--cache-dir", default=None)
    ap.add_argument("--batch-size", type=int, default=256)
    ap.add_argument("--device", default=None)
    ap.add_argument("--overwrite", action="store_true")
    a = ap.parse_args()
    predict_dirs(a.input_dirs, a.outdir, a.repo, a.cache_dir,
                 a.batch_size, a.device, a.overwrite)
