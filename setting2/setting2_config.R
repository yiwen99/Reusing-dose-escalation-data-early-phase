# =====================================================================
# Setting 2: shared configuration and data-generating process
# =====================================================================
# Two-component Stage-1 mixture with pure component-ii in Stage 2. Sourced by
# every Setting-2 method module (naive, BOIN-MEM, BARD) and by run_setting2.R.
#
# component i  = Stage-1 base curves   (BOIN-MEM Table S3)
# component ii = Stage-2 target curves (BOIN-MEM Table S4); true OBD lives here.
#
# Each method computes its own BOIN boundaries and its own per-cohort sampling
# and summarization.
# =====================================================================

suppressPackageStartupMessages(library(BOIN))

# ---------------------------------------------------------------------
# 1. Scenario rate tables (15 scenarios x 5 doses)
# ---------------------------------------------------------------------
S1_tox <- rbind(
  c(.15,.20,.22,.37,.50), c(.05,.06,.20,.22,.45), c(.05,.07,.09,.30,.45),
  c(.05,.20,.25,.28,.36), c(.09,.12,.15,.20,.28), c(.11,.12,.24,.25,.45),
  c(.05,.07,.09,.11,.13), c(.05,.12,.18,.20,.21), c(.41,.50,.54,.62,.71),
  c(.05,.20,.30,.34,.35), c(.06,.10,.24,.40,.41), c(.17,.22,.35,.41,.50),
  c(.15,.21,.22,.28,.29), c(.07,.08,.14,.26,.50), c(.06,.06,.07,.09,.15)
)
S1_eff <- rbind(
  c(.10,.45,.44,.44,.46), c(.05,.10,.45,.65,.67), c(.10,.25,.46,.47,.48),
  c(.39,.40,.40,.46,.56), c(.45,.45,.45,.45,.45), c(.23,.45,.45,.48,.35),
  c(.20,.24,.28,.46,.50), c(.05,.06,.09,.11,.14), c(.25,.34,.37,.38,.38),
  c(.45,.50,.35,.28,.20), c(.12,.14,.35,.50,.60), c(.28,.39,.51,.49,.36),
  c(.25,.29,.50,.39,.40), c(.34,.39,.52,.67,.51), c(.24,.26,.32,.41,.58)
)
S2_tox <- rbind(
  c(.15,.20,.22,.37,.50), c(.05,.06,.20,.22,.45), c(.05,.07,.09,.30,.45),
  c(.05,.20,.25,.28,.36), c(.09,.12,.15,.20,.28), c(.11,.12,.24,.25,.45),
  c(.05,.07,.09,.11,.13), c(.05,.12,.18,.20,.21), c(.41,.50,.54,.62,.71),
  c(.05,.20,.30,.34,.35), c(.23,.26,.28,.46,.50), c(.26,.35,.42,.46,.52),
  c(.15,.22,.34,.41,.42), c(.08,.12,.26,.48,.66), c(.08,.09,.18,.27,.47)
)
S2_eff <- rbind(
  c(.31,.48,.57,.57,.56), c(.10,.41,.56,.71,.72), c(.26,.51,.51,.52,.52),
  c(.40,.45,.57,.66,.67), c(.61,.61,.61,.61,.61), c(.41,.63,.57,.61,.43),
  c(.56,.61,.63,.67,.75), c(.21,.22,.37,.51,.57), c(.43,.51,.53,.54,.54),
  c(.51,.67,.63,.53,.41), c(.41,.58,.67,.69,.70), c(.48,.49,.53,.46,.38),
  c(.55,.55,.56,.44,.41), c(.53,.67,.68,.51,.28), c(.47,.49,.62,.63,.69)
)

# True OBD (component-ii / Stage-2 rates), from BOIN-MEM Table S4.
trueOBD <- c(3, 4, 2, 4, 1, 2, 5, 5, NA, 2, 3, 1, 1, 2, 3)

# ---------------------------------------------------------------------
# 2. Shared design constants
# ---------------------------------------------------------------------
J        <- 5        # dose levels
cohort1  <- 3        # Stage-1 cohort size (naive/BARD call this `cohort`)
N1       <- 30       # Stage-1 max sample size
phiT     <- 0.30     # target toxicity (BOIN)
phiE     <- 0.25     # efficacy bar (BOIN-MEM; aligned with Setting 1)
alpha    <- 0.10     # one-sided level (naive + BARD primary selector)
cutoff.eli <- 0.95   # BOIN overdose elimination cutoff
bP1      <- 1.7      # prognostic slope on CP1 (BARD-derived)
bP2      <- -1.5     # prognostic slope on CP2
bP3      <- 1.3      # prognostic slope on CP3
rho      <- 0.5      # P(C3 copies C2); proxy quality for the mis pool
lam      <- 0.3      # poolable fraction P(tau = ii); held fixed
lambda_grid <- c(0.3) # naive file sweeps this; boinmem/bard use `lam` directly
psi      <- 0.2      # Gumbel toxicity-efficacy association (0 = independent)

# ---- accrual / assessment calendar (trial duration) ----
n.per.month <- 3          # accrual rate, patients per month (BARD paper)
dlt.window  <- 1          # DLT assessment window, months (BARD paper)
eff.window  <- 3          # efficacy assessment window, months
accrual     <- "poisson"  # "poisson" or "uniform"

NSIM <- 5000         # replicates per (scenario x lambda)
seed <- 2026         # each method re-seeds from this, so run order is irrelevant

# ---------------------------------------------------------------------
# 3. Truth / reporting (component-ii curve)
# ---------------------------------------------------------------------
# Oracle MTD on the component-ii toxicity curve: highest dose with tox <= phiT.
trueMTD <- sapply(1:15, function(s) {
  ok <- which(S2_tox[s, ] <= phiT); if (length(ok) == 0) NA_integer_ else max(ok)
})
util <- function(pT, pE) 0.3 - 0.3 * pT + 0.7 * pE - 0.2 * pT * pE

# OBD-vs-MTD family label (informational, not used in any metric).
family <- ifelse(is.na(trueOBD), "null",
          ifelse(trueOBD == trueMTD, "OBD=MTD",
          ifelse(trueOBD <  trueMTD, "OBD<MTD", "OBD>MTD")))

# ---------------------------------------------------------------------
# 4. Efficacy intercept calibration (both components)
# ---------------------------------------------------------------------
# logit pE(d, tau, CP1, CP2, CP3) = beta0[s,d] + bP1*CP1 + bP2*CP2 + bP3*CP3.
# component i -> S1_eff, component ii -> S2_eff. The unweighted 8-cell mean equals
# the population marginal because each covariate is Bern(1/2).
cp_cells <- expand.grid(c1 = c(0, 1), c2 = c(0, 1), c3 = c(0, 1))
cp_lp    <- bP1 * cp_cells$c1 + bP2 * cp_cells$c2 + bP3 * cp_cells$c3   # 8 values

calib_intercept <- function(target) {
  f <- function(b) mean(plogis(b + cp_lp)) - target
  uniroot(f, c(-30, 30))$root
}
beta0_i  <- matrix(NA_real_, nrow = 15, ncol = J)   # component i  (Stage-1 curve)
beta0_ii <- matrix(NA_real_, nrow = 15, ncol = J)   # component ii (Stage-2 curve)
for (s in 1:15) for (d in 1:J) {
  beta0_i[s, d]  <- calib_intercept(S1_eff[s, d])
  beta0_ii[s, d] <- calib_intercept(S2_eff[s, d])
}

# ---- phantom MTD-1 dose (used only when Stage 1 returns mtd == 1) ----
# See the Setting-1 config for the rationale. Stage 2 is pure component ii and
# Stage 1 never enrolls at the phantom, so only the component-ii curve is needed
# and there is no component-i phantom.
ph_tox  <- 0.8 * S2_tox[, 1]
ph_eff  <- 0.8 * S2_eff[, 1]
beta0_ph <- sapply(seq_len(15), function(s) calib_intercept(ph_eff[s]))

# ---------------------------------------------------------------------
# 5. Gumbel bivariate-binary draw (joint toxicity, efficacy)
# ---------------------------------------------------------------------
draw_gumbel <- function(pT, pE, psi) {
  g   <- (exp(psi) - 1) / (exp(psi) + 1)
  del <- pT * (1 - pT) * pE * (1 - pE) * g
  p00 <- (1 - pT) * (1 - pE) + del
  p01 <- (1 - pT) * pE       - del
  p10 <- pT       * (1 - pE) - del
  u   <- runif(length(pT))
  cat <- ifelse(u < p00, 1L,
         ifelse(u < p00 + p01, 2L,
         ifelse(u < p00 + p01 + p10, 3L, 4L)))
  list(tox = c(0L, 0L, 1L, 1L)[cat], eff = c(0L, 1L, 0L, 1L)[cat])
}

# Observed per-patient utility on the Gumbel cells (BOIN-MEM scores).
util_obs <- function(tox, eff) 0.3 * (1 - tox) * (1 - eff) + 1 * (1 - tox) * eff + 0.5 * tox * eff

# ---------------------------------------------------------------------
# 6. Small analysis helpers (shared)
# ---------------------------------------------------------------------
p_binom_greater <- function(y, n, p0) pbinom(y - 1, n, p0, lower.tail = FALSE)

pick_lowest_active <- function(cand, active) {
  if (!any(active)) return(0L)
  cand[which(active)[1]]
}

# ---------------------------------------------------------------------
# Phantom-aware accessors. PHANTOM is a dose CODE, not a dose INDEX: subscripting
# a length-J vector with 0 silently drops the element, so a phantom code that
# leaks into `x[cand]` would corrupt results without erroring. Every
# phantom-capable lookup goes through these.
# ---------------------------------------------------------------------
PHANTOM <- 0L

tox_at <- function(s, d) if (d == PHANTOM) ph_tox[s]   else S2_tox[s, d]
b0_at  <- function(s, d) if (d == PHANTOM) beta0_ph[s] else beta0_ii[s, d]
# Stage-1 tally at a candidate: the phantom always has none.
n1_at  <- function(v, d) if (d == PHANTOM) 0L else v[d]

arrival_gap <- function() {
  if (accrual == "poisson") rexp(1, rate = n.per.month) else runif(1, 0, 2 / n.per.month)
}

# =====================================================================
# 7. Shared Stage-1 escalation (naive + BOIN-MEM)
# ---------------------------------------------------------------------
# Both families run the same Stage-1 BOIN escalation through
# run_stage1_shared(). The routine is toxicity-driven only, so the dose path is
# identical; it also accumulates every Stage-1 quantity either family needs
# downstream (efficacy successes for pooling/combined-p, observed utility for
# MEM, and both matched pools), under both naming conventions.
#
# Common random numbers: the drivers reseed with stage1_seed(s, r) immediately
# before each replicate, so for a given (scenario, replicate) both families draw
# identical Stage-1 patients; the families then differ only in the OBD-selection
# rule. Stage 2 flows on from the shared stream and diverges by design. BARD is
# not part of this sharing (its BF-BOIN Stage 1 with backfill is a different
# design) and keeps its own sequential stream.
#
# The per-dose cap s1 and the boundary table are defined once here for both.
# =====================================================================
s1 <- 12                # per-dose Stage-1 cap, shared by naive and BOIN-MEM

# BOIN escalation boundaries (shared; naive and BOIN-MEM used identical calls)
gb_shared <- get.boundary(target = phiT, ncohort = ceiling(N1 / cohort1),
                          cohortsize = cohort1, cutoff.eli = cutoff.eli)
btab_sh   <- gb_shared$full_boundary_tab
rn_sh     <- rownames(btab_sh)
b_esc   <- suppressWarnings(as.numeric(btab_sh[grep("^Escalate",   rn_sh), ]))
b_deesc <- suppressWarnings(as.numeric(btab_sh[grep("^Deescalate", rn_sh), ]))
b_elim  <- suppressWarnings(as.numeric(btab_sh[grep("^Eliminate",  rn_sh), ]))
lam_d   <- gb_shared$lambda_d      # BOIN de-escalation rate, used by BOIN-MEM's select_R

# Deterministic per-(scenario, replicate) seed. The multiplier exceeds any
# realistic NSIM, so scenarios never collide.
stage1_seed <- function(s, r) as.integer(seed + 100003L * s + r)

# ---------------------------------------------------------------------
# run_stage1_shared(s, lam): one replicate of the shared Stage-1 escalation.
# Returns a superset list so both modules read their own field names unchanged:
#   escalation : n_pat/n1, n_tox/yT, elim, mtd, stage1_dur
#   efficacy   : y_all/yE          utility : u1
#   clean pool : n_ii_c/n1_cl, y_ii_c, u1_cl
#   mis   pool : n_ii_m/n1_ms, y_ii_m, u1_ms
# ---------------------------------------------------------------------
run_stage1_shared <- function(s, lam) {
  sql <- sqrt(lam)                 # P(C = 1) so that P(tau = ii) = lam

  n_pat <- integer(J); n_tox <- integer(J)
  y_all <- integer(J); u_all <- numeric(J)
  n_cl  <- integer(J); y_cl  <- integer(J); u_cl <- numeric(J)   # clean pool (C1 & C2)
  n_ms  <- integer(J); y_ms  <- integer(J); u_ms <- numeric(J)   # mis   pool (C1 & C3)
  elim  <- integer(J)

  d <- 1L; enrolled <- 0L
  clock <- 0; matur <- 0

  while (enrolled < N1 && max(n_pat) < s1) {     # shared total cap N1 AND per-dose cap s1
    t_arr <- numeric(cohort1)
    for (j in seq_len(cohort1)) {
      if (j > 1L) clock <- clock + arrival_gap()
      t_arr[j] <- clock
    }

    c1   <- rbinom(cohort1, 1, sql)
    c2   <- rbinom(cohort1, 1, sql)
    copy <- rbinom(cohort1, 1, rho)
    c3   <- ifelse(copy == 1L, c2, rbinom(cohort1, 1, sql))
    is_ii     <- (c1 == 1 & c2 == 1)
    ret_clean <- is_ii
    ret_mis   <- (c1 == 1 & c3 == 1)

    cp1 <- rbinom(cohort1, 1, 0.5)
    cp2 <- rbinom(cohort1, 1, 0.5)
    cp3 <- rbinom(cohort1, 1, 0.5)

    b0 <- ifelse(is_ii, beta0_ii[s, d], beta0_i[s, d])
    pE <- plogis(b0 + bP1 * cp1 + bP2 * cp2 + bP3 * cp3)
    pT <- ifelse(is_ii, S2_tox[s, d], S1_tox[s, d])
    oc <- draw_gumbel(pT, pE, psi)
    uu <- util_obs(oc$tox, oc$eff)

    dlt_time <- ifelse(oc$tox == 1L, runif(cohort1, 0, dlt.window), dlt.window)
    matur    <- max(t_arr + dlt_time)
    repeat { clock <- clock + arrival_gap(); if (clock > matur) break }

    n_pat[d] <- n_pat[d] + cohort1; n_tox[d] <- n_tox[d] + sum(oc$tox)
    y_all[d] <- y_all[d] + sum(oc$eff); u_all[d] <- u_all[d] + sum(uu)
    n_cl[d]  <- n_cl[d]  + sum(ret_clean)
    y_cl[d]  <- y_cl[d]  + sum(oc$eff[ret_clean]); u_cl[d] <- u_cl[d] + sum(uu[ret_clean])
    n_ms[d]  <- n_ms[d]  + sum(ret_mis)
    y_ms[d]  <- y_ms[d]  + sum(oc$eff[ret_mis]);   u_ms[d] <- u_ms[d] + sum(uu[ret_mis])
    enrolled <- enrolled + cohort1

    nn <- n_pat[d]; yy <- n_tox[d]
    eb <- b_elim[nn]; if (!is.na(eb) && yy >= eb) elim[d:J] <- 1L
    if (elim[d] == 1L) { if (d == 1L) break; d <- d - 1L }
    else if (yy <= b_esc[nn])   { if (d < J && elim[d + 1L] == 0L) d <- d + 1L }
    else if (yy >= b_deesc[nn]) { if (d > 1L) d <- d - 1L }
  }

  mtd <- select.mtd(target = phiT, npts = n_pat, ntox = n_tox, cutoff.eli = cutoff.eli)$MTD
  if (mtd == 99) mtd <- NA_integer_
  stage1_dur <- matur + max(dlt.window, eff.window)

  list(
    n_pat = n_pat, n_tox = n_tox, elim = elim, mtd = mtd, stage1_dur = stage1_dur,
    # naive field names
    y_all = y_all, n_ii_c = n_cl, y_ii_c = y_cl, n_ii_m = n_ms, y_ii_m = y_ms,
    # BOIN-MEM field names (same values)
    n1 = n_pat, yT = n_tox, yE = y_all, u1 = u_all,
    n1_cl = n_cl, u1_cl = u_cl, n1_ms = n_ms, u1_ms = u_ms
  )
}
