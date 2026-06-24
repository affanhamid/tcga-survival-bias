#!/usr/bin/env bash
# Robustness sampling runs — top-variance panels (#2) and p0 sweep (#3).
# These REUSE scripts/run_cohort.py entirely (no new fitting code). They are NOT launched
# automatically — inspect the fit counts below, then run a section yourself as an overnight batch:
#
#     bash scripts/robustness_runs.sh panel     # top-var panels (24 fits)
#     bash scripts/robustness_runs.sh p0        # p0 sweep        (12 new fits)
#     bash scripts/robustness_runs.sh plan      # just print the plan, run nothing (default)
#
# Thread pinning (OPENBLAS/OMP/MKL/NUMBA = 1) makes each of the 4 chains single-threaded so a
# single process uses ~4 cores cleanly; run two sections/cohort-groups in parallel to fill 8
# cores without the BLAS thread oversubscription (load >> ncores) seen otherwise.
# Runs are resumable: any (cohort,endpoint,panel,p0) already "ok" is skipped.
set -euo pipefail
cd "$(dirname "$0")/.."
PY=${PY:-.venv/bin/python}; [ -x "$PY" ] || PY=venv/bin/python; [ -x "$PY" ] || PY=python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=4

# ---- #2 top-variance panel robustness --------------------------------------------------
# 6 cohorts spanning the ICC range: ESCA (site-confounded), KIRC (genomic star), LUSC, LUAD,
# + STAD & COAD (mid). PRIMARY endpoint only (PFI). Panels topvar:500 and topvar:2000.
# Goal: show the ~0 gene-hit result and ICC_site pattern aren't artifacts of the hallmark panel.
PANEL_COHORTS="ESCA KIRC LUSC LUAD STAD COAD"
# fits = cohorts(6) x endpoints(1) x panels(2) x models(2) = 24
panel_run() {
  "$PY" scripts/run_cohort.py --cohorts $PANEL_COHORTS --endpoints PFI \
        --panels topvar:500 topvar:2000
}

# ---- #3 p0 (shrinkage) sensitivity -----------------------------------------------------
# 3 cohorts (ESCA, KIRC, LUSC), PFI only, hallmark panel, p0 in {10,50,100}.
# run_cohort.py takes ONE --p0 per call -> loop. p0=10 PFI is already done -> auto-skipped,
# so only p0=50 and p0=100 actually fit. Goal: ICC_site dominance + null-gene robust to shrinkage.
P0_COHORTS="ESCA KIRC LUSC"
# fits = cohorts(3) x endpoints(1) x panels(1) x models(2) x p0(3) = 18 total; 6 already done -> 12 new
p0_run() {
  for p0 in 10 50 100; do
    "$PY" scripts/run_cohort.py --cohorts $P0_COHORTS --endpoints PFI --p0 "$p0"
  done
}

plan() {
  echo "PLAN (nothing launched):"
  echo "  #2 top-var panels : cohorts={$PANEL_COHORTS} x PFI x {topvar:500,topvar:2000} x {A,B}"
  echo "                      = 6 x 1 x 2 x 2 = 24 fits"
  echo "  #3 p0 sweep       : cohorts={$P0_COHORTS} x PFI x hallmark x {A,B} x p0{10,50,100}"
  echo "                      = 3 x 1 x 1 x 2 x 3 = 18 fits (6 already done -> 12 NEW)"
  echo "  thread pinning    : OPENBLAS/OMP/MKL=1 (run 2 groups in parallel to fill 8 cores)"
  echo "  launch with       : bash scripts/robustness_runs.sh {panel|p0}"
}

case "${1:-plan}" in
  panel) panel_run ;;
  p0)    p0_run ;;
  plan)  plan ;;
  *) echo "usage: $0 {panel|p0|plan}"; exit 1 ;;
esac
