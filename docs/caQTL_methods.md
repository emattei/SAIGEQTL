# Chromatin accessibility QTL mapping with SAIGE-QTL: statistical methods

Notation follows the Methods of the SAIGE-QTL preprint (Zhou et al., medRxiv 2024, v1). Code references are to commit `d89c238a`. Usage is described in [caQTL_usage.md](caQTL_usage.md).

SAIGE-QTL tests variants for association with a per-cell count, modelling cells nested within donors, without pseudobulking. Steps 1–3 make no assumption that the count belongs to a gene, so a peak count can be the phenotype and a window around the peak the cis region. No code changes are needed. This document states the model for accessibility counts, what each step computes for a peak, and when the tests are calibrated.

## Quantifying chromatin accessibility in single cells

For a peak and cell $j$ of individual $i$ ($i = 1,\dots,n$; $j = 1,\dots,m_i$; $N = \sum_i m_i$), let $y_{ij}$ be the number of Tn5 fragments with at least one insertion in the peak. A fragment with both insertions in the peak is counted once. This is paired-insertion counting (PIC; Miao & Kim 2024), the default `counting_strategy` of SnapATAC2's `make_peak_matrix` (version 2.10.0). Counting insertions instead counts such fragments twice, which inflates the number of even counts (Miao & Kim 2024; Martens et al. 2024).

A diploid cell has two copies of each locus, so $y_{ij}$ is mostly 0, 1 or 2. Typically 95–99% of the entries of a peak-by-cell matrix are zero. Under the Poisson model below,

$$\Pr(y_{ij} > 0) = 1 - e^{-\mu_{ij}}, \qquad \log\lbrace -\log(1 - \Pr(y_{ij} > 0))\rbrace = \log\mu_{ij}.$$

A model for binarized accessibility with a complementary log–log link therefore has the same linear predictor as the Poisson log-linear model. For the small $\mu_{ij}$ of most peaks the two models nearly coincide, so the counts can be modelled directly without binarizing.

## SAIGE-QTL for caQTLs

For cell $j$ of individual $i$,

$$y_{ij} \mid b_i \sim \mathrm{Poisson}(\mu_{ij}), \qquad \eta_{ij} = \log \mu_{ij} = o_{ij} + X_{di}\alpha_d + X_{cij}\alpha_c + G_i\beta + b_i.$$

- $o_{ij} = \log s_{ij}$ is a known offset, where $s_{ij}$ is a cell size factor such as the total in-peak PIC count of the cell. The model in the SAIGE-QTL paper has no offset. Here it enters through `--offsetCol`, which is added to $\eta$ unchanged and must therefore already be on the log scale.
- $X_{di}$ holds $p_d$ donor-level covariates: sex, age, batch indicators, genotype PCs and donor-level accessibility PCs.
- $X_{cij}$ holds $p_c$ cell-level covariates plus the intercept, for example TSS enrichment and leading spectral components.
- $G_i \in \lbrace 0, 1, 2\rbrace$ is the allele count.
- $b_i$ is the donor random effect.

In matrix form, $\eta = o + X_c\alpha_c + ZX_d\alpha_d + ZG\beta + Zb$, where $Z$ is the $N \times n$ incidence matrix of cells to donors. For unrelated donors $b \sim N(0, \tau_I I_n)$; for related donors $b \sim N(0, \tau_I I_n + \tau_G\Psi)$, with $\Psi$ the genetic relationship matrix (GRM). Write $X = [X_c, ZX_d]$ for the full fixed-effect design.

$\beta$ is the per-allele log fold change in expected accessibility per unit of $s_{ij}$. $\tau_I$ is the variance of donor-level log accessibility not explained by covariates. Given $b_i$, $\mathrm{Var}(y_{ij} \mid b_i) = \mu_{ij}$: the model has no cell-level overdispersion parameter.

## Estimating the variance component and other model parameters (Step 1)

Step 1 fits the null model ($\beta = 0$) by penalized quasi-likelihood (PQL; Breslow & Clayton 1993) with average-information REML (AI-REML; Gilmour et al. 1995). Given the current fit, the working vector and weights are

$$\tilde y = \hat\eta - o + \frac{y - \hat\mu}{\hat\mu}, \qquad \hat W = \mathrm{diag}(\hat\mu_{ij}),$$

(`R/SAIGE_fitGLMM_fast_func_multiV.R:6`). They define the working linear mixed model $\tilde y = X\alpha + Zb + \epsilon$ with

$$\Sigma = \tau_0\hat W^{-1} + \tau_I ZZ^\top \;\left(+\,\tau_G Z\Psi Z^\top\right).$$

For the Poisson model, $\tau_0 = 1$ is fixed and $\tau_I$ starts at 0.1 (`R/SAIGE_fitGLMM_fast_multiV.R:2281-2290`). Each iteration does two things:
1. It updates $\alpha$ and the predicted $b$ for the current $\tau$.
2. It takes one AI-REML step,

$$\tau \leftarrow \tau + \mathcal{I}^{-1}U, \qquad U_k = \tilde y^\top PA_kP\tilde y - \mathrm{tr}(PA_k), \qquad \mathcal{I}_{kl} = \tilde y^\top PA_kPA_lP\tilde y,$$

with $P = \Sigma^{-1} - \Sigma^{-1}X(X^\top\Sigma^{-1}X)^{-1}X^\top\Sigma^{-1}$ and $A_I = ZZ^\top$ (`src/Main.cpp:4871-4913`, `5214-5438`).

**Solving with $\Sigma$.** Without a GRM, $\Sigma$ is block diagonal over donors, $\Sigma_i = \hat W_i^{-1} + \tau_I\mathbf{1}\mathbf{1}^\top$. The Sherman–Morrison–Woodbury identity gives

$$\Sigma_i^{-1} = \hat W_i - \frac{\tau_I}{1 + \tau_I\hat\mu_{i\cdot}}\,\hat W_i\mathbf{1}\mathbf{1}^\top\hat W_i, \qquad \hat\mu_{i\cdot} = \sum_j\hat\mu_{ij},$$

so each product with $\Sigma^{-1}$ costs $O(N)$ (`src/Main.cpp:4629-4758`). With a GRM, the solves use preconditioned conjugate gradients (PCG).

**Traces.** $\mathrm{tr}(PA_k)$ is estimated with Hutchinson's randomized estimator, $R^{-1}\sum_{r=1}^{R}u_r^\top PA_ku_r$, where the $u_r$ are Rademacher vectors. $R$ starts at 30 (`--nrun`) and grows by 10 until the coefficient of variation of the estimate falls below `--traceCVcutoff` (0.0025; `src/Main.cpp:5440`).

**Convergence and the boundary.** With $\epsilon$ = `--tol`, iteration stops when

$$\max_k \frac{|\tau_k - \tau_k^{\mathrm{old}}|}{|\tau_k| + |\tau_k^{\mathrm{old}}| + \epsilon} < \epsilon$$

(`R/SAIGE_fitGLMM_fast_multiV.R:2454`). The same $\epsilon$ is also a floor on the variance component:
- after every update, any $\tau_k < \epsilon$ is set to 0 (`src/Main.cpp:4940-4961`);
- an iteration that sets $\tau_I = 0$ ends the fit with $\hat\tau_I = 0$ (`R/SAIGE_fitGLMM_fast_multiV.R:2450`).

The floor therefore applies to intermediate iterates as well as to the final estimate. In the two-component model, step-halving also keeps $\tau_I$ inside $(10^{-6}, 1 - 10^{-6})$ (`src/Main.cpp:4948-4958`), so the attainable range of $\hat\tau_I$ is $\lbrace 0\rbrace \cup [\epsilon, 1)$.

**Why the floor matters for accessibility.** Under the model, the donor sum of residuals $r_i = \sum_j(y_{ij} - \mu_{ij})$ has variance of about $\mu_{i\cdot}(1 + \tau_I\mu_{i\cdot})$. If $\hat\tau_I = 0$ while $\tau_I > 0$, the step 2 variance below drops the factor $1 + \tau_I\mu_{i\cdot}$, and the test statistic is inflated by about that factor. Deep donors (many cells, accessible peaks) can have large $\mu_{i\cdot}$, so even a small $\tau_I$ matters. $\epsilon$ must therefore satisfy $\epsilon\max_i\hat\mu_{i\cdot} \ll 1$, and the default `--tol=0.02` violates it once donors have more than a few expected counts ($\epsilon\mu_{i\cdot} = 0.2$ at $\mu_{i\cdot} = 10$). The condition is necessary but not sufficient: because the floor acts on the iterates, it can zero a true $\tau_I$ several times larger than $\epsilon$. Results from the deep-donor simulation ($\mu_{i\cdot} \approx 810$; 20 peaks per row, one genotype permutation):

| True $\tau_I$ | $\epsilon$ | Peaks with $\hat\tau_I = 0$ | $\lambda_{GC}$ | Type I error, $\alpha = 10^{-3}$ |
|---|---|---|---|---|
| $5\times10^{-4}$ | $10^{-3}$ | 20/20 | 1.41 | 0.0047 |
| $5\times10^{-4}$ | $10^{-4}$ | 6/20 | 1.11 | 0.0016 |
| $5\times10^{-4}$ | $10^{-5}$ | 0/20 | 1.05 | 0.00092 |
| $2\times10^{-3}$ | $10^{-3}$ | 10/20 | 1.58 | 0.019 |
| $2\times10^{-3}$ | $10^{-4}$ or $10^{-5}$ | 0/20 | 1.06 | 0.00083 |

- With $\epsilon = 10^{-3}$ and $\tau_I = 5\times10^{-4}$, $\lambda_{GC} = 1.41$ matches the predicted $1 + \tau_I\mu_{i\cdot} = 1.40$.
- $\epsilon = 10^{-4}$ still zeroed peaks whose true $\tau_I$ is five times larger.
- $\epsilon = 10^{-6}$ gave the same results as $10^{-5}$.

We use $\epsilon = 10^{-5}$, the value in the package's regression tests. It is validated up to $\mu_{i\cdot} \approx 810$; donors with larger expected counts may need a smaller $\epsilon$. The cost is convergence: the relative criterion at $10^{-5}$ is not always met. In the main simulation, 159 of 2,400 peaks (6.6%) failed the full model; the fallback below refitted 98 of them and 61 failed (Tables 1–3).

**Fallback.** If the fit fails or does not converge with `--isCovariateOffset=FALSE`, the step 1 script refits with `--isCovariateOffset=TRUE` (`extdata/step1_fitNULLGLMM_qtl.R:544-603`). The covariate effects are then fixed at their Poisson GLM estimates, fitted without the random effect, and moved into the offset; only the intercept and $\tau$ are estimated (`R/SAIGE_fitGLMM_fast_multiV.R:822-830`). Results from such peaks are reported separately below.

## Estimating the variance ratios

Let $g$ be the vector of donor genotypes, so the cell-level genotype is $Zg$. As in the paper, step 1 takes 30 random variants with minor allele count (MAC) ≥ 20 from `--plinkFile`, adding 10 at a time until the coefficient of variation of the ratios is below `--ratioCVcutoff` (0.001; computed on the `null_noXadj` ratios, `R/SAIGE_fitGLMM_fast_multiV.R:2115`). For each variant it computes two quantities (`R/SAIGE_fitGLMM_fast_multiV.R:1850-1945`):
- the variance of the score under the fitted null model, $v_1 = g^\top Z^\top\hat P Zg$;
- the variance used in step 2, $v_2 = \tilde g_d^\top Z^\top\hat\Sigma^{-1}Z\tilde g_d$, where

$$\tilde g_d = g - X_d(X_d^\top\hat MX_d)^{-1}X_d^\top\hat Mg, \qquad \hat M = \mathrm{diag}(\hat\mu_{i\cdot}),$$

adjusts the genotype for the donor-level covariates only.

The mean of $v_1/v_2$ over the variants, $\hat r$, is the `sparse` row of `<prefix>.varianceRatio.txt`. With the default `--is_fastTest=FALSE`, step 2 uses it for every variant (`src/Main.cpp:602-613`). The other rows are ratios to Poisson variances without $\tau$ and are used only by the fast test.

Without a GRM, $Z^\top\hat\Sigma^{-1}Z$ is diagonal:

$$\left(Z^\top\hat\Sigma^{-1}Z\right)_{ii} = \hat\mu_{i\cdot} - \frac{\hat\tau_I\hat\mu_{i\cdot}^2}{1 + \hat\tau_I\hat\mu_{i\cdot}} = \frac{\hat\mu_{i\cdot}}{1 + \hat\tau_I\hat\mu_{i\cdot}}.$$

$v_1$ and $v_2$ differ only in the projection: $v_1$ projects out all of $X$ in the $\hat\Sigma^{-1}$ metric, while $v_2$ projects out $X_d$ with weights $\hat M$. $\hat r$ is therefore close to 1 when cell-level covariates are nearly uncorrelated with donor identity ($\hat r = 0.9986$ in the simulations below).

## Score test for single variants (Step 2)

Step 2 works with donor sums of the null-model residuals and fitted means (`R/SAIGE_Test_main.R:494-496`):

$$r_i = \sum_j(y_{ij} - \hat\mu_{ij}), \qquad \hat\mu_{i\cdot} = \sum_j\hat\mu_{ij}.$$

The score statistic and its variance (no GRM; `src/SAIGE_test.cpp:238-264`) are

$$T = \sum_i\tilde g_{d,i}\,r_i, \qquad \widehat{\mathrm{Var}}(T) = \hat r\sum_i\tilde g_{d,i}^2\,\frac{\hat\mu_{i\cdot}}{1 + \hat\tau_I\hat\mu_{i\cdot}},$$

and $T^2/\widehat{\mathrm{Var}}(T)$ is compared with $\chi^2_1$. On a fitted model, both expressions matched the step 2 output (`Tstat` and `var`) to a relative difference below $5\times10^{-6}$.

At convergence the PQL estimating equations give $y - \hat\mu = \hat P\tilde y$, so $X^\top(y - \hat\mu) = 0$ and in particular $X_d^\top r = 0$. Then $T = g^\top Z^\top(y - \hat\mu)$, the statistic of the paper. In practice this holds up to the convergence of step 1: on the same model, the unadjusted $g^\top r$ differed from `Tstat` by up to $5\times10^{-4}$ relative. The identity does not hold after the `--isCovariateOffset=TRUE` fallback, where $T$ is the score of the covariate-adjusted genotype.

The reported effect size and standard error (`src/SAIGE_test.cpp:313-314`) are

$$\hat\beta = T/\widehat{\mathrm{Var}}(T), \qquad \mathrm{SE}(\hat\beta) = \widehat{\mathrm{Var}}(T)^{-1/2}.$$

$\hat\beta$ is the one-step score estimate of the per-allele log fold change and is accurate for small effects.

**Saddlepoint approximation.** When $|T|/\widehat{\mathrm{Var}}(T)^{1/2}$ exceeds `--SPAcutoff` (2), the p-value is recomputed by a saddlepoint approximation (SPA) (`src/SAIGE_test.cpp:693-718`). First the statistic is rescaled to the Poisson variance $v_0 = \sum_i\tilde g_{d,i}^2\hat\mu_{i\cdot}$:

$$q = T\left(v_0/\widehat{\mathrm{Var}}(T)\right)^{1/2}.$$

Its two tails are then evaluated with the cumulant generating function of $\sum_i\tilde g_{d,i}(Y_i - \hat\mu_{i\cdot})$, where the $Y_i \sim \mathrm{Poisson}(\hat\mu_{i\cdot})$ are independent:

$$K(t) = \sum_i\hat\mu_{i\cdot}\left(e^{\tilde g_{d,i}t} - \tilde g_{d,i}t - 1\right).$$

The tail probability is $1 - \Phi(w + w^{-1}\log(v/w))$, where $\hat\zeta$ solves $K'(\zeta) = q$, $w = \mathrm{sign}(\hat\zeta)\lbrace 2(\hat\zeta q - K(\hat\zeta))\rbrace^{1/2}$ and $v = \hat\zeta K''(\hat\zeta)^{1/2}$ (`src/SPA_survival.cpp:20-45`, `127-150`). The paper defers its SPA to a Supplementary Note that is not part of the preprint, so this description is read from the code. The SPA corrects the skewness of $T$ when few donors have counts. $\hat\beta$ and SE are not recomputed after SPA.

## Information content of sparse accessibility counts

**Invariance to aggregation.** Suppose the only cell-level covariate is the intercept. Then, given $b_i$, the likelihood of donor $i$'s cells factorizes as

$$\prod_j\frac{\mu_{ij}^{y_{ij}}e^{-\mu_{ij}}}{y_{ij}!} = \underbrace{\frac{y_{i\cdot}!}{\prod_j y_{ij}!}\prod_j\left(\frac{s_{ij}}{s_{i\cdot}}\right)^{y_{ij}}}_{\text{free of }\alpha,\,\beta,\,b_i}\times\frac{\mu_{i\cdot}^{y_{i\cdot}}e^{-\mu_{i\cdot}}}{y_{i\cdot}!}.$$

Here $\mu_{i\cdot} = s_{i\cdot}\exp(\alpha_0 + X_{di}\alpha_d + G_i\beta + b_i)$, and dots denote sums over cells. The cell-level model is thus the donor-level Poisson GLMM for the pseudobulk count $y_{i\cdot}$ with offset $\log s_{i\cdot}$. The same holds for the PQL working model. When $X$ has only donor-level columns, cells enter the fit, $T$ and $\widehat{\mathrm{Var}}(T)$ only through $y_{i\cdot}$ and $\hat\mu_{i\cdot}$.

We checked this in SAIGE-QTL's own fit. We fitted 50 cells per donor, and separately the same counts summed into two pseudo-cells per donor; step 1 fails with a single row per donor. With the package's Hutchinson estimator, the two fits differed by its Monte Carlo error (up to 1.2% in $\hat\tau_I$). With the trace computed exactly from the block structure of $\Sigma$ (a modification not in this commit), $\hat\tau_I$ agreed to six significant digits and step 2 statistics to $10^{-5}$. Relative to pseudobulk, the single-cell model therefore adds two things:
- cell-level covariates $X_c$;
- the count likelihood itself, in place of a linear model on normalized log counts.

**Per-donor information.** From $\widehat{\mathrm{Var}}(T)$, donor $i$ contributes information

$$\mathcal{I}_i = \frac{\mu_{i\cdot}}{1 + \tau_I\mu_{i\cdot}}$$

about $\beta$. The model has two regimes:
- **Count-limited ($\tau_I\mu_{i\cdot} \ll 1$).** $\mathcal{I}_i \approx \mu_{i\cdot}$, so more cells or deeper sequencing per donor increase power.
- **Donor-limited ($\tau_I\mu_{i\cdot} \gg 1$).** $\mathcal{I}_i \to 1/\tau_I$, so extra cells per donor add little and only more donors help.

The crossover is at $\mu_{i\cdot} = 1/\tau_I$. Both quantities are available from step 1 for every peak. A linear model on log-normalized pseudobulk instead weights donors equally after normalization and needs a pseudocount for zero totals.

## Calculating peak-level p-values for identifying caPeaks (Step 3)

For a peak spanning [start, end], the variants in the cis window [start − w, end + w] are tested (`--rangestoIncludeFile`). Larger windows add tests and dilute the peak-level combination; we use w = 10 kb. As in the paper, the peak-level p-value combines the $M$ variant p-values with the Cauchy combination (ACAT; Liu et al. 2019) with equal weights:

$$T_{\mathrm{ACAT}} = \frac{1}{M}\sum_{m=1}^{M}\tan\lbrace(0.5 - p_m)\pi\rbrace, \qquad p_{\mathrm{ACAT}} \approx \frac{1}{2} - \frac{\arctan T_{\mathrm{ACAT}}}{\pi}$$

(`get_CCT_pvalue`, `R/CCT_modified.R`). If any $p_m = 1$, the code returns $\min(1, M\min_m p_m)$ instead. Peak-level p-values are then adjusted across peaks with the Benjamini–Hochberg procedure or Storey's q-value, and the variant with the smallest p-value is reported as the lead variant.

## Data simulation for type I error evaluation

The simulation follows the paper's design but extends sparsity to ATAC levels.

**Genotypes.** We used the first 50,000 variants of `extdata/input/n.indep_100_n.cell_1_full` (100 unrelated donors); 3,854 of them pass MAC ≥ 20 and are tested.

**Null peaks.** Counts were drawn as

$$y_{ij} \sim \mathrm{Poisson}\lbrace\exp(\alpha_0 + o_{ij} + 0.1X_{1i} - 0.1X_{2i} + 0.2X_{3ij} + b_i)\rbrace,$$

with:
- $o_{ij} \sim N(0, 0.6^2)$;
- donor covariates $X_1 \sim N(0, 1)$ and $X_2 \sim \mathrm{Bernoulli}(0.5)$;
- a cell covariate $X_3 \sim N(0, 1)$;
- $b_i \sim N(0, \tau_I)$;
- $\alpha_0$ chosen per peak to give the target fraction of zeros.

**Scenarios.** The grid was:
- zeros: 80%, 95%, 99% or 99.5%;
- $\tau_I \in \lbrace 0, 0.1, 0.5\rbrace$;
- 50 or 300 cells per donor;
- 100 peaks per scenario.

A separate deep-donor scenario has 50% zeros, 1,000 cells per donor, $\tau_I \in \lbrace 5\times10^{-4}, 2\times10^{-3}\rbrace$ and 20 peaks per value.

**Analysis.** Donor labels in the genotype file were permuted three times (once for the deep-donor scenario). For each scenario we report:
- $\lambda_{GC}$;
- the type I error at $\alpha = 0.05$, $10^{-2}$, $10^{-3}$ and $10^{-4}$, with and without SPA, with 95% confidence intervals from a bootstrap over peaks (conditional on the permutations);
- the ACAT type I error for blocks of 50 adjacent variants.

Peaks where step 1 fell back are reported separately.

All results use $\epsilon = 10^{-5}$.

**Table 1.** Type I error by scenario. The last two rows are the deep-donor scenario (20 peaks, one permutation).

| Zeros | τ_I | Cells per donor | Peaks tested | Step 1 failed | λ_GC | α = 0.05 (95% CI) | α = 10⁻² | α = 10⁻³ (95% CI) | α = 10⁻⁴ | ACAT-50, α = 0.05 |
|---|---|---|---|---|---|---|---|---|---|---|
| 80% | 0 | 50 | 100 | 0 | 0.97 | 0.045 (0.043–0.047) | 0.0081 | 0.00063 (0.00057–0.00070) | 4.6e-05 | 0.032 |
| 80% | 0 | 300 | 100 | 0 | 0.97 | 0.046 (0.043–0.048) | 0.0083 | 0.00068 (0.00061–0.00076) | 5.6e-05 | 0.034 |
| 80% | 0.1 | 50 | 100 | 0 | 1.01 | 0.049 (0.049–0.049) | 0.0090 | 0.00069 (0.00064–0.00073) | 5.5e-05 | 0.035 |
| 80% | 0.1 | 300 | 100 | 0 | 1.05 | 0.053 (0.053–0.054) | 0.0104 | 0.00097 (0.00090–0.00104) | 7.5e-05 | 0.049 |
| 80% | 0.5 | 50 | 100 | 0 | 1.03 | 0.050 (0.049–0.050) | 0.0092 | 0.00066 (0.00061–0.00070) | 4.2e-05 | 0.034 |
| 80% | 0.5 | 300 | 100 | 0 | 1.03 | 0.050 (0.050–0.051) | 0.0091 | 0.00065 (0.00061–0.00069) | 4.7e-05 | 0.033 |
| 95% | 0 | 50 | 100 | 0 | 0.96 | 0.044 (0.042–0.045) | 0.0076 | 0.00057 (0.00050–0.00064) | 3.7e-05 | 0.028 |
| 95% | 0 | 300 | 100 | 0 | 0.97 | 0.045 (0.043–0.047) | 0.0083 | 0.00065 (0.00057–0.00072) | 4.1e-05 | 0.033 |
| 95% | 0.1 | 50 | 100 | 0 | 1.02 | 0.050 (0.050–0.051) | 0.0094 | 0.00074 (0.00070–0.00079) | 4.8e-05 | 0.038 |
| 95% | 0.1 | 300 | 100 | 0 | 1.04 | 0.053 (0.052–0.053) | 0.0102 | 0.00083 (0.00077–0.00089) | 6.7e-05 | 0.043 |
| 95% | 0.5 | 50 | 100 | 0 | 1.02 | 0.050 (0.049–0.050) | 0.0090 | 0.00064 (0.00060–0.00069) | 3.1e-05 | 0.032 |
| 95% | 0.5 | 300 | 100 | 0 | 1.04 | 0.051 (0.050–0.051) | 0.0094 | 0.00075 (0.00070–0.00080) | 5.2e-05 | 0.038 |
| 99% | 0 | 50 | 97 | 3% | 0.95 | 0.043 (0.041–0.045) | 0.0071 | 0.00043 (0.00038–0.00049) | 1.9e-05 | 0.022 |
| 99% | 0 | 300 | 100 | 0 | 0.98 | 0.045 (0.044–0.047) | 0.0079 | 0.00066 (0.00058–0.00073) | 4.0e-05 | 0.032 |
| 99% | 0.1 | 50 | 100 | 0 | 0.98 | 0.045 (0.043–0.046) | 0.0073 | 0.00047 (0.00041–0.00052) | 2.5e-05 | 0.023 |
| 99% | 0.1 | 300 | 100 | 0 | 1.05 | 0.053 (0.052–0.053) | 0.0103 | 0.00082 (0.00077–0.00087) | 6.4e-05 | 0.044 |
| 99% | 0.5 | 50 | 100 | 0 | 1.01 | 0.047 (0.046–0.048) | 0.0077 | 0.00048 (0.00043–0.00053) | 2.6e-05 | 0.024 |
| 99% | 0.5 | 300 | 100 | 0 | 1.05 | 0.053 (0.053–0.054) | 0.0101 | 0.00088 (0.00083–0.00093) | 3.9e-05 | 0.044 |
| 99.5% | 0 | 50 | 84 | 16% | 0.95 | 0.041 (0.039–0.043) | 0.0061 | 0.00030 (0.00026–0.00034) | 1.0e-05 | 0.016 |
| 99.5% | 0 | 300 | 100 | 0 | 0.94 | 0.043 (0.040–0.045) | 0.0074 | 0.00054 (0.00047–0.00061) | 3.7e-05 | 0.028 |
| 99.5% | 0.1 | 50 | 91 | 9% | 0.97 | 0.043 (0.041–0.045) | 0.0066 | 0.00038 (0.00033–0.00044) | 1.2e-05 | 0.019 |
| 99.5% | 0.1 | 300 | 100 | 0 | 1.02 | 0.051 (0.049–0.052) | 0.0095 | 0.00077 (0.00070–0.00084) | 4.3e-05 | 0.039 |
| 99.5% | 0.5 | 50 | 67 | 33% | 1.03 | 0.047 (0.046–0.048) | 0.0073 | 0.00039 (0.00032–0.00046) | 1.8e-05 | 0.019 |
| 99.5% | 0.5 | 300 | 100 | 0 | 1.04 | 0.052 (0.051–0.052) | 0.0094 | 0.00072 (0.00067–0.00077) | 3.5e-05 | 0.037 |
| 50% | 0.0005 | 1000 | 20 | 0 | 1.05 | 0.052 (0.050–0.054) | 0.0096 | 0.00092 (0.00069–0.00114) | 7.8e-05 | 0.042 |
| 50% | 0.002 | 1000 | 20 | 0 | 1.06 | 0.055 (0.054–0.056) | 0.0109 | 0.00083 (0.00067–0.00099) | 2.6e-05 | 0.044 |

**Table 2.** Type I error by step 1 outcome, pooled over the 24 main scenarios.

| Step 1 fit | Peaks | λ_GC range over scenarios | α = 0.05 | α = 10⁻² | α = 10⁻³ | α = 10⁻⁴ |
|---|---|---|---|---|---|---|
| full model | 2241 | 0.94–1.05 | 0.048 | 0.0085 | 0.00064 | 4.1e-05 |
| covariate-offset fallback | 98 | 1.01–1.08 | 0.052 | 0.0095 | 0.00074 | 3.8e-05 |

**Table 3.** Step 1 outcome by the number of donors with nonzero counts, pooled over the 24 main scenarios.

| Donors with nonzero counts (of 100) | Peaks | Step 1 failed | Fallback attempted |
|---|---|---|---|
| ≤ 15 | 11 | 36% | 36% |
| 15–20 | 112 | 33% | 43% |
| 20–25 | 129 | 12% | 29% |
| 25–30 | 62 | 3% | 18% |
| 30–40 | 200 | 0.5% | 8% |
| 40–60 | 92 | 1% | 8% |
| > 60 | 1794 | 0% | 2% |

**Results.**
- **Tail ($\alpha \le 10^{-3}$).** No scenario was inflated: every scenario was at or below the nominal level.
- **$\alpha = 0.05$.** Scenarios with $\tau_I > 0$ and 300 cells per donor, and the deep-donor scenarios, were slightly liberal: up to 0.055, with confidence intervals excluding 0.05, and $\lambda_{GC}$ up to 1.06.
- **Conservative scenarios.** Scenarios with $\tau_I = 0$, and sparse scenarios with 50 cells per donor, were conservative, more so in the tail (0.00030 at $\alpha = 10^{-3}$ for 99.5% zeros, $\tau_I = 0$, 50 cells).
- **ACAT.** ACAT over 50-variant blocks stayed at or below nominal (0.016–0.049), following the variant-level tail.
- **Fallback peaks.** These were slightly liberal at $\alpha = 0.05$ but not in the tail (Table 2).
- **Step 1 failures.** Failures concentrated in peaks with fewer than 30 donors having nonzero counts (Table 3). They are failures of the convergence criterion at $\epsilon = 10^{-5}$: with $\epsilon = 10^{-3}$ the same 2,400 peaks gave one fallback, no failures and $\lambda_{GC}$ 0.94–1.05, but $10^{-3}$ is not safe for deep donors (Step 1).

The confidence intervals come from a bootstrap over peaks. They are conditional on the genotype permutations shared by all peaks, so they understate the uncertainty.

## Assumptions and limitations

- **No cell-level overdispersion.** If counts are overdispersed within donors, $\mathrm{Var}(T)$ is underestimated. A diagnostic is the Pearson statistic $\sum_{ij}(y_{ij} - \hat\mu_{ij})^2/\hat\mu_{ij}$ divided by $N - p$, which should be near 1.
- **Bounds on $\hat\tau_I$.** The floor at `--tol` and the cap below 1 both underestimate $\mathrm{Var}(T)$ when the true $\tau_I$ lies outside $[\epsilon, 1)$; see Step 1.
- **PQL bias.** PQL tends to underestimate variance components for sparse counts (Breslow & Lin 1995). Without cell-level covariates, the fit depends on the data only through donor totals (see the information-content section), so the bias is governed by $\mu_{i\cdot}$ rather than by per-cell counts.
- **Few donors with counts.** Peaks with few donors having nonzero counts more often fail the step 1 convergence criterion at $\epsilon = 10^{-5}$ (Table 3). They carry little information about $\beta$, so dropping them before step 1 saves compute; calibration does not require it.
- **SE before SPA.** SE is computed before the SPA, so for SPA-adjusted variants $\hat\beta/\mathrm{SE}$ does not reproduce the reported p-value.
- **Poisson only.** `--traitType=count_nb` has no SPA branch in step 2.
- **Mapping bias.** Reads carrying the alternative allele can map less well; WASP-style filtering (van de Geijn et al. 2015) removes this bias.
- **Peak set.** Peaks must be called on all donors jointly. Peaks called within subsets of donors can create associations through the peak definition itself.

## Notation

| Symbol | Meaning |
|---|---|
| $y_{ij}$ | PIC count of the peak in cell $j$ of donor $i$ |
| $s_{ij}$, $o_{ij}$ | cell size factor and offset $o_{ij} = \log s_{ij}$ |
| $X_{di}$, $X_{cij}$ | donor-level and cell-level covariates; $X = [X_c, ZX_d]$ |
| $Z$ | $N \times n$ incidence matrix of cells to donors |
| $G_i$, $g$ | allele count of donor $i$; vector over donors |
| $b_i$, $\tau_I$, $\tau_G$, $\Psi$ | donor random effect, its variance components and the GRM |
| $\hat\mu_{ij}$, $\hat\mu_{i\cdot}$ | fitted null mean of a cell; sum over the donor's cells |
| $r_i$ | donor sum of null-model residuals |
| $\hat W$, $\Sigma$, $P$ | PQL weights, working covariance and projection |
| $\tilde g_d$ | genotype adjusted for donor-level covariates |
| $\hat r$ | variance ratio (`sparse` row) |
| $\epsilon$ | `--tol`: convergence threshold and floor on $\tau$ |

## References

- Breslow NE, Clayton DG. Approximate inference in generalized linear mixed models. *J Am Stat Assoc* 88, 9–25 (1993).
- Breslow NE, Lin X. Bias correction in generalised linear mixed models with a single component of dispersion. *Biometrika* 82, 81–91 (1995).
- Gilmour AR, Thompson R, Cullis BR. Average information REML: an efficient algorithm for variance parameter estimation in linear mixed models. *Biometrics* 51, 1440–1450 (1995).
- Hutchinson MF. A stochastic estimator of the trace of the influence matrix for Laplacian smoothing splines. *Commun Stat Simul Comput* 19, 433–450 (1990).
- Liu Y, et al. ACAT: a fast and powerful p value combination method for rare-variant analysis in sequencing studies. *Am J Hum Genet* 104, 410–421 (2019).
- Martens LD, et al. Modeling fragment counts improves single-cell ATAC-seq analysis. *Nat Methods* 21, 28–31 (2024).
- Miao Z, Kim J. Uniform quantification of single-nucleus ATAC-seq data with Paired-Insertion Counting (PIC) and a model-based insertion rate estimator. *Nat Methods* 21, 32–36 (2024).
- Storey JD. A direct approach to false discovery rates. *J R Stat Soc B* 64, 479–498 (2002).
- van de Geijn B, et al. WASP: allele-specific software for robust molecular quantitative trait locus discovery. *Nat Methods* 12, 1061–1063 (2015).
- Zhang K, et al. A fast, scalable and versatile tool for analysis of single-cell omics data. *Nat Methods* 21, 217–227 (2024).
- Zhou W, et al. Efficiently controlling for case-control imbalance and sample relatedness in large-scale genetic association studies. *Nat Genet* 50, 1335–1341 (2018).
- Zhou W, et al. Efficient and accurate mixed model association tool for single-cell eQTL analysis. medRxiv (2024); *Nat Genet* (2026), doi:10.1038/s41588-026-02798-2.
