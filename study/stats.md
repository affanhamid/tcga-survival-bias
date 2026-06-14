
### 1.2 Survival Analysis (deepening what you know)

**What you already know:**
Right-censoring, Cox partial likelihood, C-index, Weibull model, IBS.

**What you need to add:**

**Parametric survival families and when to use them**
- Exponential: constant hazard over time (memoryless). Rarely realistic.
- Weibull: hazard increases or decreases monotonically. Most flexible
  parametric choice. What you already use.
- Log-normal: hazard rises then falls. Good for some cancers where
  early mortality is high but survivors stabilise.
- Log-logistic: similar to log-normal but heavier tails.
- How to choose: LOO-CV via ArviZ, or WAIC. You will use this for
  model comparison in the paper.

**Accelerated failure time (AFT) models**
- Alternative to Cox. Instead of modelling the hazard, you model
  log(survival time) directly as a linear function of covariates.
- More interpretable: a coefficient of 0.3 means "this gene multiplies
  expected survival time by e^0.3 = 1.35"
- Your Weibull model IS an AFT model. This is the frame you use
  in the paper.

**What to read:**
- Klein & Moeschberger (2003). "Survival Analysis: Techniques for Censored
  and Truncated Data." Chapter 12 (parametric models). Available as PDF.
- The PyMC survival analysis examples:
  https://www.pymc.io/projects/examples/en/latest/survival_analysis/

---

### 1.3 Bayesian Methods (deepening what you know)

**What you already know:**
MCMC, NUTS sampler, posterior distributions, credible intervals, 
hierarchical models, random intercepts, ICC.

**What you need to add:**

**Horseshoe prior**
- A continuous shrinkage prior for high-dimensional sparse regression
- Each coefficient β_j has its own local shrinkage parameter λ_j and
  shares a global shrinkage parameter τ
- The "horseshoe" shape: most coefficients get shrunk near zero (noise),
  but a few large signals pass through nearly unshrunken
- Better than Lasso (L1) because it doesn't uniformly penalise large effects
- The regularised horseshoe (Piironen & Vehtari 2017) is the version
  you use — it adds a slab component for numerical stability

**What to read:**
- Piironen J & Vehtari A. (2017). "Sparsity information and regularization
  in the horseshoe and other shrinkage priors." Electronic Journal of Statistics.
  Read sections 1-3. This is the paper behind the horseshoe you'll implement.
- Austin Rochford's PyMC implementation:
  https://austinrochford.com/posts/2021-05-29-horseshoe-pymc3.html
  (Code examples you'll adapt directly.)

**Model comparison in Bayesian framework**
- LOO-CV (Leave-One-Out Cross Validation): computed from posterior samples
  using ArviZ's az.loo(). Compares models without refitting.
- WAIC (Widely Applicable Information Criterion): similar to LOO-CV,
  also from ArviZ. Use both.
- ELPD (Expected Log Predictive Density): the quantity LOO-CV estimates.
  Higher = better model.
- You use these to compare: Weibull vs log-normal vs log-logistic,
  and naive model vs site-adjusted model.

**What to read:**
- Vehtari A, Gelman A, Gabry J. (2017). "Practical Bayesian model
  evaluation using leave-one-out cross-validation and WAIC."
  Statistics and Computing. Read sections 1-4.
- ArviZ documentation on model comparison:
  https://python.arviz.org/en/stable/api/generated/arviz.compare.html

**Convergence diagnostics**
- R-hat < 1.01: chains have mixed well
- ESS (Effective Sample Size) > 400 per parameter: enough samples
- Divergences: ideally zero, max ~0.1% of draws. More = model
  misspecification or need for reparameterisation.
- You already know these from your project. The paper requires
  reporting them formally.
