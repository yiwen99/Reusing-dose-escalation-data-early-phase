# =====================================================================
# Setting 1: shared configuration and data-generating process
# =====================================================================
# A single common dose-outcome curve is shared by Stage 1 and Stage 2, so all
# pooling is unbiased. Sourced by every Setting-1 method module (naive,
# BOIN-MEM, BARD) and by run_setting1.R.
#
# The common curve is the Stage-II target curve of the Setting-2 DGP (BOIN-MEM
# Simulation 2, Table S4): S_tox / S_eff below are the Setting-2 S2_tox / S2_eff
# tables and trueOBD is the corresponding OBD vector.
#
# Holds only what the three method families share: rate tables, truth, efficacy
# calibration, the Gumbel joint draw, the observed utility, the accrual calendar
# and small analysis helpers. Each method computes its own BOIN boundaries and
# its own per-cohort sampling and summarization.
# =====================================================================

suppressPackageStartupMessages(library(BOIN))

# ---------------------------------------------------------------------
# 1. Common rate table (Stage-II target curve, Table S4)
# ---------------------------------------------------------------------
S_tox <- rbind(
  c(.15,.20,.22,.37,.50),  # 1
  c(.05,.06,.20,.22,.45),  # 2
  c(.05,.07,.09,.30,.45),  # 3
  c(.05,.20,.25,.28,.36),  # 4
  c(.09,.12,.15,.20,.28),  # 5
  c(.11,.12,.24,.25,.45),  # 6
  c(.05,.07,.09,.11,.13),  # 7
  c(.05,.12,.18,.20,.21),  # 8
  c(.41,.50,.54,.62,.71),  # 9
  c(.05,.20,.30,.34,.35),  # 10
  c(.23,.26,.28,.46,.50),  # 11
  c(.26,.35,.42,.46,.52),  # 12
  c(.15,.22,.34,.41,.42),  # 13
  c(.08,.12,.26,.48,.66),  # 14
  c(.08,.09,.18,.27,.47)   # 15
)

S_eff <- rbind(
  c(.31,.48,.57,.57,.56),  # 1
  c(.10,.41,.56,.71,.72),  # 2
  c(.26,.51,.51,.52,.52),  # 3
  c(.40,.45,.57,.66,.67),  # 4
  c(.61,.61,.61,.61,.61),  # 5
  c(.41,.63,.57,.61,.43),  # 6
  c(.56,.61,.63,.67,.75),  # 7
  c(.21,.22,.37,.51,.57),  # 8
  c(.43,.51,.53,.54,.54),  # 9
  c(.51,.67,.63,.53,.41),  # 10
  c(.41,.58,.67,.69,.70),  # 11
  c(.48,.49,.53,.46,.38),  # 12
  c(.55,.55,.56,.44,.41),  # 13
  c(.53,.67,.68,.51,.28),  # 14
  c(.47,.49,.62,.63,.69)   # 15
)

# ---------------------------------------------------------------------
# 2. Shared design constants
# ---------------------------------------------------------------------
J        <- 5        # dose levels
cohort1  <- 3        # Stage-1 cohort size (naive/BARD call this `cohort`)
N1       <- 30       # Stage-1 max sample size
phiT     <- 0.30     # target toxicity (BOIN)
phiE     <- 0.25     # efficacy bar (BOIN-MEM Simulation)
alpha    <- 0.10     # one-sided level (naive + BARD primary selector)
cutoff.eli <- 0.95   # BOIN overdose elimination cutoff
bP1      <- 1.7      # prognostic slope on CP1 (BARD-derived)
bP2      <- -1.5     # prognostic slope on CP2
bP3      <- 1.3      # prognostic slope on CP3
psi      <- 0.2      # Gumbel toxicity-efficacy association (0 = independent)

# ---- poolable fraction for the naive partial-pooling arms ----
# Probability that a Stage-1 patient is tagged poolable (the `h` tag in
# run_stage1_shared below). Read only by the naive family, via n_pool / y_pool.
# The tag is an independent per-patient Bernoulli draw, so the pooled subset is
# a random fraction, not a fixed-size sample. The naive module's arm labels
# carry this fraction (pool_0.3, combined_p_0.3) and must be updated with it.
pool_frac <- 0.30

# ---- accrual / assessment calendar (trial duration) ----
# Shared by every method so the duration columns are comparable. Any change here
# propagates to all three families through this single definition.
n.per.month <- 3          # accrual rate, patients per month (BARD paper)
dlt.window  <- 1          # DLT assessment window, months (BARD paper)
eff.window  <- 3          # efficacy assessment window, months
accrual     <- "poisson"  # "poisson" or "uniform"

NSIM <- 5000         # replicates per scenario
seed <- 2026         # each method re-seeds from this, so run order is irrelevant

# ---------------------------------------------------------------------
# 3. True MTD, true OBD, OBD-vs-MTD family
# ---------------------------------------------------------------------
# True MTD: highest dose with true tox <= phiT (NA if none).
trueMTD <- sapply(1:15, function(s) {
  ok <- which(S_tox[s, ] <= phiT); if (length(ok) == 0) NA_integer_ else max(ok)
})

# True OBD on the common curve, from BOIN-MEM Table S4; cross-checked below
# against the utility-argmax recompute.
trueOBD <- c(3, 4, 2, 4, 1, 2, 5, 5, NA, 2, 3, 1, 1, 2, 3)

# BOIN-MEM utility (psi00=.3, psi01=1, psi10=0, psi11=.5):
#   U(d) = 0.3 - 0.3*pT + 0.7*pE - 0.2*pT*pE
util <- function(pT, pE) 0.3 - 0.3 * pT + 0.7 * pE - 0.2 * pT * pE

# Cross-check: recompute the utility-argmax OBD from the rates and warn if it
# disagrees with the hardcoded vector.
obd_from_rates <- sapply(1:15, function(s) {
  adm <- which(S_tox[s, ] <= phiT & S_eff[s, ] > phiE)
  if (length(adm) == 0) return(NA_integer_)
  adm[which.max(util(S_tox[s, adm], S_eff[s, adm]))]   # ties -> lowest dose
})
if (!isTRUE(all.equal(obd_from_rates, trueOBD))) {
  cat("WARNING: hardcoded trueOBD disagrees with utility-argmax recompute:\n")
  print(rbind(hardcoded = trueOBD, recomputed = obd_from_rates))
}

# OBD-vs-MTD family label (informational, not used in any metric). OBD>MTD can
# never fire: the OBD is chosen only among doses with tox <= phiT.
family <- ifelse(is.na(trueOBD), "null",
          ifelse(is.na(trueMTD), "null",
          ifelse(trueOBD == trueMTD, "OBD=MTD",
          ifelse(trueOBD <  trueMTD, "OBD<MTD", "OBD>MTD"))))

# ---------------------------------------------------------------------
# 4. Efficacy intercept calibration (single common curve)
# ---------------------------------------------------------------------
# logit pE(d, CP1, CP2, CP3) = beta0[s,d] + bP1*CP1 + bP2*CP2 + bP3*CP3.
# Solve beta0 so the 8-cell covariate average (CP ~ Bern(1/2)) matches S_eff.
cp_cells <- expand.grid(c1 = c(0, 1), c2 = c(0, 1), c3 = c(0, 1))
cp_lp    <- bP1 * cp_cells$c1 + bP2 * cp_cells$c2 + bP3 * cp_cells$c3   # 8 values

calib_intercept <- function(target) {
  f <- function(b) mean(plogis(b + cp_lp)) - target
  uniroot(f, c(-30, 30))$root
}
beta0 <- matrix(NA_real_, nrow = 15, ncol = J)
for (s in 1:15) for (d in 1:J) beta0[s, d] <- calib_intercept(S_eff[s, d])

# ---- phantom MTD-1 dose (used only when Stage 1 returns mtd == 1) ----
# A hypothetical dose one level below dose 1, at 80 percent of dose 1's true
# toxicity and efficacy. It exists only in Stage 2 and only in the naive and BARD
# families; Stage 1 never enrolls there, so it has no Stage-1 data under any
# pooling rule. Selecting it is scored as a failure to select an OBD.
# The efficacy target is a marginal rate, so it needs its own calibrated
# intercept; scaling beta0 would not scale the marginal rate.
ph_tox  <- 0.8 * S_tox[, 1]
ph_eff  <- 0.8 * S_eff[, 1]
beta0_ph <- sapply(seq_len(15), function(s) calib_intercept(ph_eff[s]))

# ---------------------------------------------------------------------
# 5. Gumbel bivariate-binary draw (joint toxicity, efficacy)
# ---------------------------------------------------------------------
# pT, pE are vectors (one entry per patient). Margins are exact in psi; the
# association enters only the joint cells.
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
# One-sided exact binomial p-value for H0: pE <= phiE vs H1: pE > phiE.
# n = 0 returns 1 (empty sample: no evidence).
p_binom_greater <- function(y, n, p0) pbinom(y - 1, n, p0, lower.tail = FALSE)

# Pick the LOWEST active candidate (both active -> MTD-1; none -> 0). cand sorted.
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

tox_at <- function(s, d) if (d == PHANTOM) ph_tox[s]  else S_tox[s, d]
b0_at  <- function(s, d) if (d == PHANTOM) beta0_ph[s] else beta0[s, d]
# Stage-1 tally at a candidate: the phantom always has none.
n1_at  <- function(v, d) if (d == PHANTOM) 0L else v[d]

# One inter-arrival gap (months); both accrual models have mean 1/n.per.month.
arrival_gap <- function() {
  if (accrual == "poisson") rexp(1, rate = n.per.month) else runif(1, 0, 2 / n.per.month)
}

# =====================================================================
# 7. Shared Stage-1 escalation (naive + BOIN-MEM)
# ---------------------------------------------------------------------
# Both families run the same Stage-1 BOIN escalation via run_stage1_shared(s).
# The routine also carries the naive random poolable subset (the h tag, drawn at
# rate pool_frac) and the observed utility, so both families read what they need
# with no remapping.
#
# Common random numbers: the drivers reseed with stage1_seed(s, r) immediately
# before each replicate, so for a given (scenario, replicate) both families draw
# identical Stage-1 patients. The h tag is part of the shared draw sequence, so
# the streams stay in lock-step. Stage 2 flows on and diverges by design. BARD is
# not part of this sharing and keeps its own sequential stream.
#
# The per-dose cap s1 and the boundary table are defined once here for both.
# =====================================================================
s1 <- 12                # per-dose Stage-1 cap, shared by naive and BOIN-MEM

gb_shared <- get.boundary(target = phiT, ncohort = ceiling(N1 / cohort1),
                          cohortsize = cohort1, cutoff.eli = cutoff.eli)
btab_sh   <- gb_shared$full_boundary_tab
rn_sh     <- rownames(btab_sh)
b_esc   <- suppressWarnings(as.numeric(btab_sh[grep("^Escalate",   rn_sh), ]))
b_deesc <- suppressWarnings(as.numeric(btab_sh[grep("^Deescalate", rn_sh), ]))
b_elim  <- suppressWarnings(as.numeric(btab_sh[grep("^Eliminate",  rn_sh), ]))
lam_d   <- gb_shared$lambda_d      # BOIN de-escalation rate, used by BOIN-MEM's select_R

stage1_seed <- function(s, r) as.integer(seed + 100003L * s + r)

# ---------------------------------------------------------------------
# run_stage1_shared(s): one replicate of the shared Stage-1 escalation.
# Returns a superset list so both modules read their own field names unchanged:
#   escalation : n_pat/n1, n_tox/yT, elim, mtd, stage1_dur
#   naive      : n_all, y_all, n_pool, y_pool     BOIN-MEM: yE (=y_all), u1
# ---------------------------------------------------------------------
run_stage1_shared <- function(s) {
  n_pat  <- integer(J); n_tox <- integer(J)
  y_all  <- integer(J); u_all <- numeric(J)
  n_pool <- integer(J); y_pool <- integer(J)          # random poolable subset
  elim   <- integer(J)

  d <- 1L; enrolled <- 0L
  clock <- 0; matur <- 0

  while (enrolled < N1 && max(n_pat) < s1) {     # shared total cap N1 AND per-dose cap s1
    t_arr <- numeric(cohort1)
    for (j in seq_len(cohort1)) {
      if (j > 1L) clock <- clock + arrival_gap()
      t_arr[j] <- clock
    }

    cp1 <- rbinom(cohort1, 1, 0.5)
    cp2 <- rbinom(cohort1, 1, 0.5)
    cp3 <- rbinom(cohort1, 1, 0.5)
    h   <- rbinom(cohort1, 1, pool_frac)           # poolability tag

    pE <- plogis(beta0[s, d] + bP1 * cp1 + bP2 * cp2 + bP3 * cp3)
    pT <- rep(S_tox[s, d], cohort1)
    oc <- draw_gumbel(pT, pE, psi)
    uu <- util_obs(oc$tox, oc$eff)

    dlt_time <- ifelse(oc$tox == 1L, runif(cohort1, 0, dlt.window), dlt.window)
    matur    <- max(t_arr + dlt_time)
    repeat { clock <- clock + arrival_gap(); if (clock > matur) break }

    n_pat[d]  <- n_pat[d]  + cohort1; n_tox[d] <- n_tox[d] + sum(oc$tox)
    y_all[d]  <- y_all[d]  + sum(oc$eff); u_all[d] <- u_all[d] + sum(uu)
    n_pool[d] <- n_pool[d] + sum(h)
    y_pool[d] <- y_pool[d] + sum(oc$eff[h == 1L])
    enrolled  <- enrolled  + cohort1

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
    n_all = n_pat, y_all = y_all, n_pool = n_pool, y_pool = y_pool,
    # aliases kept for back-compatibility (same values)
    n_half = n_pool, y_half = y_pool,
    # BOIN-MEM field names (same values)
    n1 = n_pat, yT = n_tox, yE = y_all, u1 = u_all
  )
}
