# =====================================================================
# Setting 1: BOIN-MEM design, following the reference package
# =====================================================================
# The Stage-2 logic follows the BOIN-MEM reference implementation
# (BOINMEM-functions.R, BOINMEM-implementation.R) wherever the package and the
# paper text differ. The points of divergence are:
#   - reselection fires after a fixed number of cohorts (sint / cohort2), not
#     when a dose first reaches sint cumulatively; hence the two-block enroll
#     structure below;
#   - the Eq (11) OBD sampler is the package weighted average of independent
#     Beta draws (maxu_MEM), not a draw from the mixture itself;
#   - the untried-dose valve is enabled at reselection as well as at Stage 1,
#     with the package's extra guards;
#   - the top-JR set is built by sequential argmax (pick, remove, recompute),
#     matching the select_dose_BM p2i loop, rather than one ranking of Pj;
#   - the per-cohort toxicity check is candidate-only within an enrollment
#     block, with allocation computed once at block entry;
#   - elimination is re-derived on every call rather than held in a persistent
#     eliminated[] vector. This equals a permanent flag in the common case (a
#     flagged dose stops accruing, so its counts freeze) and is more lenient
#     only where a lower dose un-flags as it accrues more data, which raises the
#     ceiling and re-admits the dose above it. The paper specifies permanent
#     elimination, so this is a paper-vs-package difference.
#
# The package's table() zero-cell behaviour is not replicated; tabulate() is
# used instead.
#
# Requires setting1_config.R sourced first (S_tox, S_eff, beta0, draw_gumbel,
# util, util_obs, arrival_gap, trueOBD, trueMTD, family, shared constants).
# =====================================================================

# ---- method-specific constants ----
# s1, the per-dose Stage-1 cap, is defined in setting1_config.R section 7.
cohort2  <- 4          # Stage-2 cohort size
N2       <- 40         # total Stage-2
s2       <- 40         # per-dose Stage-2 cap = N2 (paper convention)
sint     <- s2 / 2     # reselection trigger
JR       <- 2
C_T      <- 0.95
C_E      <- 0.90
a_prior  <- 1
b_prior  <- 1
M_rank   <- 10000L     # MC draws, matching the package
M_obd    <- 10000L     # MC draws for the final Eq (11) OBD dominance
RUN_REFERENCE <- FALSE # kept for schema parity

# The package OBD sampler is the omega-weighted average of independent Beta
# draws (maxu_MEM), not the Eq (11) mixture.
SAMPLER_PKG <- "reference"

VARIANTS <- data.frame(
  name  = c("boinmem1", "boinmem2", "boinmem3"),
  w     = c(0.5,        0.5,        0.0),
  gamma = c(0.0,      -Inf,         0.0),
  stringsAsFactors = FALSE
)
SAMPLERS <- SAMPLER_PKG

# ---------------------------------------------------------------------
# Component-specific cohort draw (Setting 1: single common curve)
# ---------------------------------------------------------------------
draw_cohort <- function(s, d, n) {
  cp1 <- rbinom(n, 1, 0.5); cp2 <- rbinom(n, 1, 0.5); cp3 <- rbinom(n, 1, 0.5)
  pE  <- plogis(beta0[s, d] + bP1 * cp1 + bP2 * cp2 + bP3 * cp3)
  oc  <- draw_gumbel(rep(S_tox[s, d], n), pE, psi)
  list(nt = sum(oc$tox), ne = sum(oc$eff), u = sum(util_obs(oc$tox, oc$eff)),
       tox = oc$tox)
}

# ---------------------------------------------------------------------
# BOIN boundaries (b_esc/b_deesc/b_elim) and lam_d are SHARED: defined in
# setting1_config.R section 7 (used by run_stage1_shared and by select_R).
# ---------------------------------------------------------------------

# =====================================================================
# MEM machinery
# =====================================================================
# The exclusivity test uses the nearest other candidate on each side separately
# (lo_rival below, hi_rival above), which is equivalent to the supplement S1
# all-candidate rule and to the package's ceiling((gap-1)/2) construction, and
# identical to the single-rival form at JR = 2.
build_models <- function(j1, lo_rival, hi_rival, J) {
  allow_lo <- function(k) if (is.na(lo_rival)) TRUE else abs(j1 - k) <= abs(lo_rival - k)
  allow_hi <- function(k) if (is.na(hi_rival)) TRUE else abs(j1 - k) <= abs(hi_rival - k)
  md <- 0L; k <- j1 - 1L; while (k >= 1L && allow_lo(k)) { md <- md + 1L; k <- k - 1L }
  mu <- 0L; k <- j1 + 1L; while (k <= J  && allow_hi(k)) { mu <- mu + 1L; k <- k + 1L }
  rows <- list()
  for (a in 0:md) for (b in 0:mu) { om <- integer(J); om[(j1 - a):(j1 + b)] <- 1L
    rows[[length(rows) + 1L]] <- om }
  rbind(integer(J), do.call(rbind, rows))
}
candidate_stats <- function(om, n1, u1, u2j, n2j) {
  a_l <- a_prior + as.numeric(om %*% u1)        + u2j
  b_l <- b_prior + as.numeric(om %*% (n1 - u1)) + (n2j - u2j)
  standalone <- lbeta(a_prior + u1, b_prior + (n1 - u1)) - lbeta(a_prior, b_prior)
  log_ml <- lbeta(a_l, b_l) - lbeta(a_prior, b_prior) + as.numeric((1 - om) %*% standalone)
  list(a = a_l, b = b_l, log_ml = log_ml)
}
lse <- function(x) { m <- max(x); if (!is.finite(m)) return(m); m + log(sum(exp(x - m))) }
model_weights <- function(stats, om, w, gamma) {
  L1 <- length(stats$log_ml); sizes <- rowSums(om)
  log_pre <- numeric(L1); log_pre[1] <- log(w)
  s2v <- sizes[-1]
  lg <- if (is.infinite(gamma) && gamma < 0) ifelse(s2v == 1, 0, -Inf) else gamma * log(s2v)
  log_pre[-1] <- log1p(-w) + (lg - lse(lg))
  lp <- stats$log_ml + log_pre; e <- exp(lp - max(lp)); e / sum(e)
}
draw_util <- function(a, b, wts, M, sampler) {
  if (sampler == "mixture") { idx <- sample.int(length(wts), M, replace = TRUE, prob = wts)
    rbeta(M, a[idx], b[idx])
  } else { mat <- vapply(seq_along(wts), function(l) rbeta(M, a[l], b[l]), numeric(M))
    as.numeric(mat %*% wts) }
}
# JR = 2 checks: dose 3 (upper rival 4) -> L = 3; dose 4 (lower rival 3) -> L = 2.
stopifnot(nrow(build_models(3, NA, 4, 5)) - 1L == 3L,
          nrow(build_models(4, 3, NA, 5)) - 1L == 2L)

# ---------------------------------------------------------------------
# Gates (paper direction) and simple-Beta dominance
# ---------------------------------------------------------------------
tox_admissible <- function(yT, n1) {
  ok <- rep(FALSE, J)
  for (j in 1:J) if (n1[j] > 0) {
    below <- 1:j
    pass  <- all((1 - pbeta(phiT, 1 + yT[below], 1 + n1[below] - yT[below])) < C_T)
    ok[j] <- pass
  }
  ok
}
eff_admissible <- function(yE, n1) (n1 > 0) & ((1 - pbeta(phiE, 1 + yE, 1 + n1 - yE)) > (1 - C_E))

simple_dominance <- function(doses, a, b, M) {
  if (length(doses) == 1) return(setNames(1, doses))
  X <- sapply(doses, function(j) rbeta(M, a[j], b[j]))   # M x length(doses)
  win <- max.col(X, ties.method = "first")
  p <- tabulate(win, nbins = length(doses)) / M
  setNames(p, doses)
}

# Sequential argmax ranking, matching select_dose_BM's p2i loop: the winner is
# removed and the dominance probability recomputed over the remaining set for
# each of the JR picks.
seq_top <- function(A, a, b, M) {
  R <- integer(0); pool <- A
  for (p in 1:JR) {
    if (length(pool) == 0) break
    if (length(pool) == 1L) { R <- c(R, pool); break }
    P <- simple_dominance(pool, a, b, M)
    win <- as.integer(names(P))[which.max(P)]
    R <- c(R, win); pool <- setdiff(pool, win)
  }
  R
}

# Untried-dose valve, with the package guards: it fires only if |R| < JR, the
# highest tried dose k < J, no dose is flagged toxic by Eq (14), and
# y_T1k / n1k < lam_d. The package's act guard (no BOIN-eliminated dose) is
# omitted because run_stage1_shared does not expose that flag and it is
# redundant given that k is the highest tried dose.
select_R <- function(A, a, b, n1, yT, M, allow_untried = TRUE, tox_flag_any = FALSE) {
  if (length(A) == 0) return(integer(0))
  R <- seq_top(A, a, b, M)                                # Step 4: sequential dominance
  if (allow_untried && length(R) < JR && !tox_flag_any) {
    tried <- which(n1 > 0); k <- if (length(tried)) max(tried) else 0L
    if (k >= 1L && k < J && (yT[k] / n1[k]) < lam_d && !((k + 1L) %in% R)) R <- c(R, k + 1L)
  }
  sort(unique(R))
}

# ---------------------------------------------------------------------
# One replicate: shared Stage 1, package Stage 2, then per-variant MEM OBD.
# Stage 1 is run_stage1_shared(s) from the config. The package Stage-2 enroll
# block (enroll_Stage2_BM) is an inner closure inside one_rep, so it captures
# n1, u1, yT1 and yE1 directly.
# ---------------------------------------------------------------------
NARM <- nrow(VARIANTS)
arm_names <- VARIANTS$name

one_rep <- function(s) {
  st1 <- run_stage1_shared(s)
  n1 <- st1$n1; yT1 <- st1$yT; yE1 <- st1$yE; u1 <- st1$u1
  expo <- n1
  na <- rep(NA_real_, NARM)
  stage1_dur <- st1$stage1_dur
  clock  <- stage1_dur
  matur2 <- stage1_dur; eff_all <- stage1_dur

  # ---- admissible set A and initial R (Eq 5 simple posterior) ----
  tox1 <- tox_admissible(yT1, n1)
  A0 <- which(tox1 & eff_admissible(yE1, n1))
  if (length(A0) == 0) {
    dc <- if (!any(tox1)) 2L else 1L
    return(list(sel = rep(0L, NARM), door = rep(dc, NARM),
        est = na, est_util = na, eom = na, expo = expo,
        dur1 = stage1_dur, dur2 = 0, dur2_roll = 0, dur = stage1_dur,
        s2_cohorts = 0, reached_s2 = FALSE))
  }

  a5 <- a_prior + u1;  b5 <- b_prior + (n1 - u1)           # Eq (5) params, all doses
  tox_flag1 <- any((1:J %in% which(!tox1)) & (n1 > 0))     # valve guard
  R  <- select_R(A0, a5, b5, n1, yT1, M_rank, allow_untried = TRUE, tox_flag_any = tox_flag1)
  if (length(R) == 0) return(list(sel = rep(0L, NARM), door = rep(1L, NARM),
      est = na, est_util = na, eom = na, expo = expo,
      dur1 = stage1_dur, dur2 = 0, dur2_roll = 0, dur = stage1_dur,
        s2_cohorts = 0, reached_s2 = FALSE))

  # ---- Stage-2 accumulators ----
  yT2 <- integer(J); yE2 <- integer(J); u2 <- numeric(J); n2 <- integer(J)

  # Replicates enroll_Stage2_BM. Allocation is computed once from the cumulative
  # pooled utility posteriors of the candidates R and held fixed across the
  # block; toxicity is checked per candidate that received patients,
  # deactivating that candidate and all higher ones; efficacy is deactivated
  # once at the end. Returns updated accumulators and the surviving-candidate
  # flags act2, aligned to R.
  enroll_block <- function(R, ncoh, n2, yT2, yE2, u2, expo, s2enr, clock, eff_all) {
    act2 <- rep(TRUE, length(R))
    aU <- a_prior + (u1 + u2);  bU <- b_prior + (n1 - u1) + (n2 - u2)   # pooled utility
    alloc <- simple_dominance(R, aU, bU, M_rank)
    names(alloc) <- as.character(R)
    if (ncoh >= 1L) for (ci in 1:ncoh) {
      if (!any(act2)) break
      av <- alloc; av[!act2] <- 0
      if (sum(av) == 0) break
      currch <- as.integer(rmultinom(1, cohort2, prob = av))
      tox_vec <- integer(0)
      for (k in seq_along(R)) {
        nk <- currch[k]; if (nk == 0L) next
        dk <- R[k]
        ct <- draw_cohort(s, dk, nk)
        n2[dk] <- n2[dk] + nk; yT2[dk] <- yT2[dk] + ct$nt
        yE2[dk] <- yE2[dk] + ct$ne; u2[dk] <- u2[dk] + ct$u
        expo[dk] <- expo[dk] + nk
        tox_vec <- c(tox_vec, ct$tox)
        pt <- 1 - pbeta(phiT, 1 + yT1[dk] + yT2[dk],
                        1 + (n1[dk] - yT1[dk]) + (n2[dk] - yT2[dk]))
        if (pt > C_T) act2[k:length(R)] <- FALSE      # deactivate this + higher candidates
      }
      s2enr <- s2enr + cohort2
      # toxicity-gated enrollment clock
      t_arr2 <- numeric(cohort2)
      for (j in seq_len(cohort2)) { clock <- clock + arrival_gap(); t_arr2[j] <- clock }
      if (length(tox_vec) > 0) {
        dlt_time2 <- ifelse(tox_vec == 1L, runif(length(tox_vec), 0, dlt.window), dlt.window)
        clock <- max(c(t_arr2[seq_along(tox_vec)] + dlt_time2, t_arr2))
      } else clock <- max(t_arr2)
      eff_all <- max(eff_all, max(t_arr2) + eff.window)
      if (any((n1 + n2)[R] >= s2)) break             # per-dose cap (cumulative)
    }
    # efficacy deactivation, once at block end
    for (k in seq_along(R)) {
      dk <- R[k]
      pe <- pbeta(phiE, 1 + yE1[dk] + yE2[dk],
                  1 + (n1[dk] - yE1[dk]) + (n2[dk] - yE2[dk]))
      if (pe > C_E) act2[k] <- FALSE
    }
    list(n2 = n2, yT2 = yT2, yE2 = yE2, u2 = u2, expo = expo,
         s2enr = s2enr, clock = clock, eff_all = eff_all, act2 = act2)
  }

  # ---- Stage-2 sequence: enroll(sint) -> Step 3 reselect -> enroll(rest) ----
  # Reselection is structural: it happens once, after a fixed sint / cohort2
  # cohorts. Enrollment always follows reselection, so a reselected dose cannot
  # enter Eq (11) with n2 = 0.
  ncoh_pre  <- as.integer(round(sint / cohort2))          # cohorts before reselection
  ncoh_post <- as.integer(round(N2 / cohort2)) - ncoh_pre # cohorts after (no recycling)
  s2enr <- 0L

  # Steps 1-2: enroll the first block on the Stage-1 R
  b1 <- enroll_block(R, ncoh_pre, n2, yT2, yE2, u2, expo, s2enr, clock, eff_all)
  n2 <- b1$n2; yT2 <- b1$yT2; yE2 <- b1$yE2; u2 <- b1$u2; expo <- b1$expo
  s2enr <- b1$s2enr; clock <- b1$clock; eff_all <- b1$eff_all

  # Step 3: single reselection from the admissible set, with the valve enabled.
  # Uses the same select_R as Stage-1 Steps 4 and 5, matching the package calling
  # select_dose_BM again. Efficacy must mature at this one evaluation point.
  clock <- max(clock, eff_all)
  a7 <- a_prior + u1 + u2;  b7 <- b_prior + (n1 - u1) + (n2 - u2)
  tox_all <- tox_admissible(yT1 + yT2, n1 + n2)
  Apool <- which(tox_all & eff_admissible(yE1 + yE2, n1 + n2))
  tox_flag_any <- any((1:J %in% which(!tox_all)) & ((n1 + n2) > 0))
  Rn <- select_R(Apool, a7, b7, n1 + n2, yT1 + yT2, M_rank,
                 allow_untried = TRUE, tox_flag_any = tox_flag_any)
  if (length(Rn) > 0) R <- Rn

  # Step 4: enroll the remaining block on the reselected R
  if (length(R) > 0 && ncoh_post > 0L) {
    b2 <- enroll_block(R, ncoh_post, n2, yT2, yE2, u2, expo, s2enr, clock, eff_all)
    n2 <- b2$n2; yT2 <- b2$yT2; yE2 <- b2$yE2; u2 <- b2$u2; expo <- b2$expo
    s2enr <- b2$s2enr; clock <- b2$clock; eff_all <- b2$eff_all
    act2_final <- b2$act2
  } else {
    act2_final <- if (length(R) > 0) rep(TRUE, length(R)) else logical(0)
  }

  # final OBD (Step 5) needs all Stage-2 efficacy, so the reported total duration
  # includes the last cohorts' efficacy maturation (a single trailing efficacy wait).
  matur2 <- max(clock, eff_all)

  # ---- final candidate set: the act2 flags surviving Step 4 ----
  # Step 5 chooses the OBD among the candidates that survived Stage-2 toxicity
  # and efficacy deactivation, not a fresh admissibility re-derivation.
  tox_f <- tox_admissible(yT1 + yT2, n1 + n2)
  dcf <- if (!any(tox_f)) 2L else 1L
  Rf <- if (length(R) > 0) R[act2_final] else integer(0)
  if (length(Rf) == 0) return(list(sel = rep(0L, NARM), door = rep(dcf, NARM),
      est = na, est_util = na, eom = na, expo = expo,
      dur1 = stage1_dur, dur2 = matur2 - stage1_dur,
      dur2_roll = if (s2enr > 0) s2enr / n.per.month + eff.window else 0,
      dur = matur2, s2_cohorts = s2enr / cohort2, reached_s2 = s2enr > 0))

  # ---- final OBD per variant: Eq (11) via the weighted-average sampler ----
  saved <- get(".Random.seed", envir = .GlobalEnv)

  nc <- length(Rf)
  oms <- vector("list", nc); stats <- vector("list", nc)
  for (i in 1:nc) {
    lo_rival <- if (i == 1L) NA_integer_ else Rf[i - 1L]     # per-side rivals
    hi_rival <- if (i == nc) NA_integer_ else Rf[i + 1L]
    oms[[i]]   <- build_models(Rf[i], lo_rival, hi_rival, J)
    stats[[i]] <- candidate_stats(oms[[i]], n1, u1, u2[Rf[i]], n2[Rf[i]])
  }

  sel <- integer(NARM); door <- integer(NARM)
  est <- rep(NA_real_, NARM); est_u <- rep(NA_real_, NARM); eom <- rep(NA_real_, NARM)
  arm <- 0L
  for (vi in 1:nrow(VARIANTS)) {
    w <- VARIANTS$w[vi]; gm <- VARIANTS$gamma[vi]
    wts <- lapply(1:nc, function(i) model_weights(stats[[i]], oms[[i]], w, gm))
    pmean <- sapply(1:nc, function(i) sum(wts[[i]] * stats[[i]]$a /
                                          (stats[[i]]$a + stats[[i]]$b)))
    eomega <- sapply(1:nc, function(i) 1 - wts[[i]][1])
    for (samp in SAMPLERS) {
      arm <- arm + 1L
      if (nc == 1L) { pick <- 1L } else {
        draws <- lapply(1:nc, function(i) draw_util(stats[[i]]$a, stats[[i]]$b,
                                                    wts[[i]], M_obd, samp))
        X <- do.call(cbind, draws)
        win <- max.col(X, ties.method = "first")
        Pj <- tabulate(win, nbins = nc) / M_obd
        pick <- which.max(Pj)
      }
      dd <- Rf[pick]
      sel[arm] <- dd; door[arm] <- 0L
      est[arm]   <- if (n2[dd] > 0) yE2[dd] / n2[dd] else NA_real_   # Stage-2 efficacy readout
      est_u[arm] <- pmean[pick]
      eom[arm]   <- eomega[pick]
    }
  }
  assign(".Random.seed", saved, envir = .GlobalEnv)
  list(sel = sel, door = door, est = est, est_util = est_u, eom = eom, expo = expo,
      dur1 = stage1_dur, dur2 = matur2 - stage1_dur,
      dur2_roll = if (s2enr > 0) s2enr / n.per.month + eff.window else 0,
      dur = matur2, s2_cohorts = s2enr / cohort2, reached_s2 = s2enr > 0)
}

# ---------------------------------------------------------------------
# Driver: emits the shared output tables plus a diagnostics file
# ---------------------------------------------------------------------
run_boinmem <- function(cfg) {
  set.seed(seed)
  full_rows <- list(); si <- 0
  alloc_rows <- list(); ai <- 0
  diag_rows <- list()

  for (s in 1:15) {
    ob <- trueOBD[s]; mt <- trueMTD[s]
    tox_over <- which(S_tox[s, ] > phiT)

    sel  <- matrix(NA_integer_, NSIM, NARM); door <- matrix(NA_integer_, NSIM, NARM)
    est  <- matrix(NA_real_, NSIM, NARM);    estu <- matrix(NA_real_, NSIM, NARM)
    eom  <- matrix(NA_real_, NSIM, NARM);    expo <- matrix(0L, NSIM, J)
    durv  <- rep(NA_real_, NSIM)
    dur1v <- rep(NA_real_, NSIM)
    dur2v <- rep(NA_real_, NSIM)
    dur2rv <- rep(NA_real_, NSIM)
    s2cv  <- rep(NA_real_, NSIM)
    s2v   <- logical(NSIM)

    for (r in 1:NSIM) {
      set.seed(stage1_seed(s, r))   # CRN: identical Stage 1 across families
      o <- one_rep(s)
      sel[r, ] <- o$sel; door[r, ] <- o$door; est[r, ] <- o$est
      estu[r, ] <- o$est_util; eom[r, ] <- o$eom; expo[r, ] <- o$expo
      durv[r] <- o$dur; dur1v[r] <- o$dur1; dur2v[r] <- o$dur2
      dur2rv[r] <- o$dur2_roll; s2cv[r] <- o$s2_cohorts; s2v[r] <- o$reached_s2
    }

    # ---- allocation, shared across the three variants (same trial) ----
    Nm    <- mean(rowSums(expo))
    belowN <- if (is.na(ob)) NA else mean(rowSums(expo[, seq_len(J) < ob, drop = FALSE]))
    atN    <- if (is.na(ob)) NA else mean(expo[, ob])
    overN  <- if (is.na(ob)) NA else mean(rowSums(expo[, seq_len(J) > ob, drop = FALSE]))
    abvN   <- if (length(tox_over) == 0) 0 else mean(rowSums(expo[, tox_over, drop = FALSE]))
    ai <- ai + 1
    alloc_rows[[ai]] <- data.frame(
      Scn = s, family = family[s], trueOBD = ob, trueMTD = mt, N_used = Nm,
      below_n = belowN, below_pct = 100 * belowN / Nm,
      at_n = atN,       at_pct = 100 * atN / Nm,
      over_n = overN,   over_pct = 100 * overN / Nm,
      aboveMTD_n = abvN, aboveMTD_pct = 100 * abvN / Nm,
      PCS_MTD = NA_real_,          # BOIN-MEM selects an OBD, not an MTD
      duration        = mean(durv),
      dur_stage1      = mean(dur1v),
      dur_stage2      = mean(dur2v),      # sequential: the design as implemented
      dur_stage2_roll = mean(dur2rv),     # one-shot counterfactual, accounting only
      dur_med         = median(durv),
      dur_q90         = unname(quantile(durv, 0.90)),
      dur_ifS2        = if (any(s2v)) mean(durv[s2v]) else NA,
      pct_noS2        = 100 * mean(!s2v),
      s2_cohorts      = mean(s2cv),
      stringsAsFactors = FALSE)

    # ---- per-variant selection metrics ----
    for (a in 1:NARM) {
      sa <- sel[, a]; da <- door[, a]; ea <- est[, a]; eu <- estu[, a]; em <- eom[, a]
      selected <- !is.na(sa) & sa > 0
      PCS   <- if (is.na(ob)) NA else 100 * mean(sa == ob, na.rm = TRUE)
      Under <- if (is.na(ob)) NA else 100 * mean(selected & sa < ob)
      Over  <- if (is.na(ob)) NA else 100 * mean(selected & sa > ob)
      noOBD <- 100 * mean(!selected)
      noMTD <- 100 * mean(da == 2)          # no tolerable dose
      noEff <- 100 * mean(da == 1)          # tolerable dose(s) but none efficacious
      tox_flag <- logical(NSIM); tox_flag[selected] <- S_tox[s, sa[selected]] > phiT
      ToxSel  <- 100 * mean(tox_flag); any_sel <- 100 * mean(selected)
      d1 <- 100*mean(sa==1); d2 <- 100*mean(sa==2); d3 <- 100*mean(sa==3)
      d4 <- 100*mean(sa==4); d5 <- 100*mean(sa==5)
      if (any(selected)) {
        Est <- mean(ea[selected], na.rm = TRUE); True <- mean(S_eff[s, sa[selected]])
        EstU <- mean(eu[selected]); TrueU <- mean(util(S_tox[s, sa[selected]], S_eff[s, sa[selected]]))
        Eom <- mean(em[selected]); Eom_sd <- sd(em[selected])
      } else { Est <- NA; True <- NA; EstU <- NA; TrueU <- NA; Eom <- NA; Eom_sd <- NA }
      si <- si + 1
      full_rows[[si]] <- data.frame(Scn = s, family = family[s], method = arm_names[a],
        trueOBD = ob, trueMTD = mt, PCS = PCS, UnderOBD = Under, OverOBD = Over,
        noOBD = noOBD, noMTD = noMTD, noEff = noEff, ToxSel = ToxSel, select_any = any_sel,
        Est = Est, True = True, Boost = Est - True,
        sel_d1 = d1, sel_d2 = d2, sel_d3 = d3, sel_d4 = d4, sel_d5 = d5,
        stringsAsFactors = FALSE)
      diag_rows[[si]] <- data.frame(Scn = s, method = arm_names[a],
        Est_util = EstU, True_util = TrueU, Boost_util = EstU - TrueU,
        Eomega_jj = Eom, Eomega_jj_sd = Eom_sd, stringsAsFactors = FALSE)
    }
    cat(sprintf("  [boinmem] scenario %2d done\n", s))
  }

  full  <- do.call(rbind, full_rows)
  alloc <- do.call(rbind, alloc_rows)
  diagn <- do.call(rbind, diag_rows)

  sel_tab <- full[, c("Scn","family","trueOBD","trueMTD","method",
                 "PCS","UnderOBD","OverOBD","noOBD","noMTD","noEff","ToxSel")]
  dist_tab <- data.frame(Scn = full$Scn, Method = full$method,
                   d1 = full$sel_d1, d2 = full$sel_d2, d3 = full$sel_d3,
                   d4 = full$sel_d4, d5 = full$sel_d5, noOBD = full$noOBD,
                   stringsAsFactors = FALSE)

  od <- cfg$outdir
  write.csv(full,  file.path(od, "setting1_boinmem_results_full.csv"),      row.names = FALSE)
  write.csv(sel_tab,  file.path(od, "setting1_boinmem_selection.csv"),   row.names = FALSE)
  write.csv(dist_tab, file.path(od, "setting1_boinmem_seldist.csv"),     row.names = FALSE)
  write.csv(alloc,    file.path(od, "setting1_boinmem_allocation.csv"),  row.names = FALSE)
  write.csv(diagn,    file.path(od, "setting1_boinmem_diagnostics.csv"), row.names = FALSE)

  cat("[boinmem] wrote setting1_boinmem_{selection,seldist,allocation}.csv, results_full, diagnostics to", od, "\n")

  invisible(list(full = full, sel = sel_tab, dist = dist_tab, alloc = alloc, diagn = diagn))
}
