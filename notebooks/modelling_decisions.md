# Modelling decisions & gotchas (03_Modelling)

Hard-won fixes for the Weibull AFT + horseshoe survival model. These apply
identically to Model A (naive) and Model B (with site intercepts).

## The model

Weibull accelerated-failure-time (AFT) with a regularized horseshoe prior over
gene coefficients:

- `t_i ~ Weibull(alpha, lambda_i)`, scale `lambda_i = exp(mu + beta · x_i)`
- Right-censored patients contribute the log survival function via `pm.Potential`
- Horseshoe prior shrinks the many irrelevant genes toward 0

## Four fixes required to get clean sampling

A pilot (top-100 genes by variance, 495 LUSC patients) went from **400/400
divergences** to **0 divergences, r_hat ≈ 1.00, ess_bulk ~400–500** after all of:

### 1. macOS Apple-Silicon EOFError → `mp_ctx="spawn"`
`pm.sample(cores>1)` crashed with `EOFError` (worker process died). On macOS ARM
PyMC defaults the multiprocessing start method to `fork`, which is unsafe with
Accelerate/BLAS and fails dataset-size-dependently. Fix: `pm.sample(..., mp_ctx="spawn")`.
(Irrelevant if using `nuts_sampler="nutpie"`, which avoids Python multiprocessing.)

### 2. Weibull scale sign: `exp(eta)`, NOT `exp(-eta)`
The scale must be `lambda = exp(mu + beta·x)`. An earlier `exp(-eta)` made the
scale ≈ `exp(-6) ≈ 0.0025` while survival times are in the hundreds–thousands of
days → likelihood in a catastrophic region → exploding gradients → 100% divergence.
Censored term must match: `log S(t) = -(t/lambda)^alpha`.

### 3. Non-centered + REGULARIZED horseshoe
- Non-centered: `z ~ Normal(0,1); beta = z * tau * lam` (avoids Neal's funnel).
- Regularized (Piironen & Vehtari 2017): add a slab `c2 ~ InverseGamma` and use
  `lam_tilde = sqrt(c2 * lam^2 / (c2 + tau^2 * lam^2))`, `beta = z * tau * lam_tilde`.
  Plain HalfCauchy `lam` has heavy tails → some `beta_j` blow up → `eta` explodes
  the Weibull. The slab caps coefficients near `slab_scale` (~2 on the log-time scale).
- `sigma_tau = (p0/(p-p0)) * (sigma_est/sqrt(n))` (Piironen & Vehtari) is correct as-is.

### 4. Drop `OS_time <= 0`
LUSC has one censored patient with `OS_time == 0`. `log(0) = -inf` in the censored
term → NaN gradient → every transition diverges (chains frozen at init; the telltale
sign is identical degenerate diagnostics across runs: `ess≈2.08`, `r_hat≈9.9e15`).
`build_model_a` now filters `df["OS_time"] > 0` internally (drops 1 → 494 patients).
The parquet itself is left untouched.

## Numerical note
Use the stable censored form `log_sf = -exp(alpha * (log t - eta))` rather than
`-(t * exp(-eta))**alpha` to avoid intermediate overflow — but only valid once `t > 0`
(fix #4 must be in place first).

## Sampling / performance
- Pilot config: `draws=200, tune=500, chains=2, target_accept=0.95`.
- If divergences linger, bump `target_accept` to 0.99.
- Scaling 100 → 13,859 genes is much slower. Use `nuts_sampler="nutpie"`
  (`pip install nutpie numba`) — same NUTS algorithm, Rust backend, ~2–10× faster
  and better mass-matrix adaptation. Or keep a variance / MSigDB-hallmark gene subset.

## Environment
- Python 3.12 (`./venv`). pymc 6.0.1, pytensor 3.0.5, arviz 1.2.0 — these require
  3.12+, so a Python 3.11 downgrade is NOT viable (arviz 1.x is 3.12-only).
