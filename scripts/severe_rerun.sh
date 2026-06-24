#!/usr/bin/env bash
# Re-run ONLY the convergence-SEVERE cohort/endpoint pairs flagged by triage_convergence.py,
# at higher target_accept + more tune/draws to fix the funnel mixing. run_cohort.py fits BOTH
# models per (cohort,endpoint), so the 10 severe traces = 9 (cohort,endpoint) pairs = 18 fits.
# --force overwrites the existing (bad) hallmark/p10 traces + their results.csv rows.
#
# Chained: waits for any in-flight --combat run to finish first (avoids 8-core oversubscription),
# then runs two balanced groups in parallel, thread-pinned so 2 procs x 4 chains = 8 cores cleanly.
set -uo pipefail
cd /root/tcga-survival-bias
PY=.venv/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=4
COMMON="--force --panels hallmark --p0 10 --target-accept 0.999 --tune 2000 --draws 2000 --chains 4"

echo "[severe] $(date) waiting for any --combat run to finish ..."
while pgrep -f "run_cohort.py --combat" >/dev/null 2>&1; do sleep 60; done
echo "[severe] $(date) combat done — starting severe re-runs"

# Group A (heaviest first): BRCA both endpoints, then GBM/OS, COAD/OS
groupA() {
  $PY scripts/run_cohort.py --cohorts BRCA --endpoints OS PFI $COMMON
  $PY scripts/run_cohort.py --cohorts GBM  --endpoints OS     $COMMON
  $PY scripts/run_cohort.py --cohorts COAD --endpoints OS     $COMMON
}
# Group B: LGG/OS, OV/PFI, BLCA/PFI, SARC/PFI, PAAD/PFI
groupB() {
  $PY scripts/run_cohort.py --cohorts LGG  --endpoints OS  $COMMON
  $PY scripts/run_cohort.py --cohorts OV   --endpoints PFI $COMMON
  $PY scripts/run_cohort.py --cohorts BLCA --endpoints PFI $COMMON
  $PY scripts/run_cohort.py --cohorts SARC --endpoints PFI $COMMON
  $PY scripts/run_cohort.py --cohorts PAAD --endpoints PFI $COMMON
}

groupA > /root/severeA.log 2>&1 &
A=$!
groupB > /root/severeB.log 2>&1 &
B=$!
wait $A $B
echo "[severe] $(date) ALL severe re-runs complete"
echo "[severe] NOTE: re-fit traces have NO log_likelihood group yet — re-run add_loglik.py +"
echo "[severe]       model_compare.py on these cohorts to refresh the LOO comparison."
