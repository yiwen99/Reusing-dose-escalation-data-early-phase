# =====================================================================
# Setting 1: BARD (backfill + adaptive randomization), class "bard"
# =====================================================================
# Arms (3): bard_car3 (minimization balances CP1, CP2, CP3), bard_car2 (CP3
# masked), and bfboin_sr (1:1 simple randomization).
#
# bfboin_sr is not the published BOIN-SR comparator. It keeps BF-BOIN Stage 1
# with backfilling, the same n_cap, the same retained Stage-1 records, the same
# Stage-1-plus-Stage-2 combination at analysis, the same
# N2star = N2 - (retained Stage-1) sample size, and the same OBD selectors. The
# only difference from BARD is the Stage-2 assignment rule: an independent fair
# coin per patient instead of conditional Pocock-Simon minimization at gamma_car.
# The published BOIN-SR uses plain BOIN with no backfill and does not combine
# Stage-1 data, so it differs from BARD on three axes at once and cannot isolate
# the randomization rule. Simple randomization here is an independent fair coin,
# not permuted blocks and not forced 1:1 balance, so arm sizes are unequal and
# alloc_imb is larger than for the CAR arms by construction.
#
# Phantom MTD-1 dose: when Stage 1 returns mtd == 1 the candidate set becomes
# c(PHANTOM, 1), so there are always two candidates to randomize between. The
# phantom is a decoy at 80 percent of dose 1's toxicity and efficacy with no
# Stage-1 records; picking it is a selection failure, coded sel = -1, door = 3,
# reported as phantomSel, with its Stage-2 patients in phantom_n / phantom_pct.
#
# Requires setting1_config.R sourced first (S_tox, S_eff, beta0, trueOBD,
# trueMTD, family, phiT, phiE, alpha, cutoff.eli, bP1..bP3, psi, calendar, NSIM,
# seed).
#
# BARD keeps its own Gumbel helpers (gumbel_cell / outcome_at_u) because Stage-2
# common random numbers require the patient's uniform to be split out of the
# draw, and its own BF-BOIN boundary table, which is longer than the
# naive/BOIN-MEM table because backfill lets a dose accumulate up to n_cap
# patients.
# =====================================================================

# ---------------------------------------------------------------------
# method-specific constants
# ---------------------------------------------------------------------
ncohort  <- 10         # -> N_esc = ncohort * cohort1 = 30 escalation patients
n_stop   <- 9          # n.earlystop: stop at current dose with "stay"
n_cap    <- 12         # backfill cap per dose (escalation + backfill counted)
p.saf    <- 0.6 * phiT # 0.18
p.tox    <- 1.4 * phiT # 0.42

N2       <- 40         # Stage-2 target total at the two doses (inclusive of Stage 1)
n2_single <- 20        # per-dose target when only one candidate exists. The
                       # phantom makes length(cand) == 2 whenever mtd is not NA,
                       # so this branch is unreachable as configured; it is kept
                       # so the file still behaves if the phantom is disabled.
gamma_car <- 0.95      # minimization assignment probability (Pocock-Simon)

CAR_ENGINE <- "minirand"   # "minirand" (package) or "closedform" (paired closed form)
if (CAR_ENGINE == "minirand") suppressPackageStartupMessages(library(Minirand))

w0_alloc <- 0          # allocation-balance weight (Minirand faithful default 0)
delta_ni <- 0.05       # BARD efficacy-rate noninferiority margin (secondary)

GATES <- "study"                       # "paper" (BARD) or "study" (BOIN-MEM convention)
if (GATES == "paper") { C_T <- 0.90; C_E <- 0.95 } else { C_T <- 0.95; C_E <- 0.90 }
a_prior <- 1; b_prior <- 1     # Beta / Dirichlet flat prior
M_util  <- 4000L               # MC draws for the quasi-Bernoulli dominance selector

# ---------------------------------------------------------------------
# Gumbel cell selection with the patient's uniform supplied (CRN), plus draws
# ---------------------------------------------------------------------
gumbel_cell <- function(pT, pE, psi, u) {
  g <- (exp(psi) - 1) / (exp(psi) + 1); del <- pT*(1-pT)*pE*(1-pE)*g
  p00 <- (1-pT)*(1-pE)+del; p01 <- (1-pT)*pE-del; p10 <- pT*(1-pE)-del
  cat <- ifelse(u < p00,1L, ifelse(u < p00+p01,2L, ifelse(u < p00+p01+p10,3L,4L)))
  list(tox = c(0L,0L,1L,1L)[cat], eff = c(0L,1L,0L,1L)[cat])
}
draw_gumbel <- function(pT, pE, psi) gumbel_cell(pT, pE, psi, runif(length(pT)))
util_obs <- function(tox, eff) 0.3*(1-tox)*(1-eff) + 1*(1-tox)*eff + 0.5*tox*eff

# ---------------------------------------------------------------------
# Single-patient draws (common curve; identical mechanism in both stages)
# ---------------------------------------------------------------------
draw_cov <- function() {
  list(cp1 = rbinom(1,1,0.5), cp2 = rbinom(1,1,0.5), cp3 = rbinom(1,1,0.5))
}
outcome_at <- function(s, d, cov) {
  pE <- plogis(beta0[s,d] + bP1*cov$cp1 + bP2*cov$cp2 + bP3*cov$cp3)
  oc <- draw_gumbel(S_tox[s,d], pE, psi)
  list(tox = oc$tox, eff = oc$eff, util = util_obs(oc$tox, oc$eff))
}
# Stage-2 outcome draw, routed through the phantom-aware accessors because d may
# be PHANTOM. outcome_at above is Stage-1 only, where the phantom cannot occur,
# so it keeps its direct beta0 / S_tox lookups.
outcome_at_u <- function(s, d, cp1, cp2, cp3, u) {
  pE <- plogis(b0_at(s, d) + bP1*cp1 + bP2*cp2 + bP3*cp3)
  oc <- gumbel_cell(tox_at(s, d), pE, psi, u)
  list(tox = oc$tox, eff = oc$eff, util = util_obs(oc$tox, oc$eff))
}
draw_pat_s1 <- function(s, d) {
  cov <- draw_cov()
  oc  <- outcome_at(s, d, cov)
  dlt.time <- if (oc$tox == 1L) runif(1, 0, dlt.window) else dlt.window
  list(tox = oc$tox, eff = oc$eff, util = oc$util,
       cp1 = cov$cp1, cp2 = cov$cp2, cp3 = cov$cp3, dlt.time = dlt.time)
}

# ---------------------------------------------------------------------
# BOIN boundaries (extended so the table covers up to N_esc + J*n_cap)
# ---------------------------------------------------------------------
gb <- get.boundary(target = phiT,
                   ncohort = ncohort + J * ceiling(n_cap / cohort1),
                   cohortsize = cohort1,
                   n.earlystop = ncohort * cohort1 + J * n_cap,
                   p.saf = p.saf, p.tox = p.tox, cutoff.eli = cutoff.eli)$full_boundary_tab
b_esc  <- suppressWarnings(as.numeric(gb[2, ]))   # escalate   if y <= b_esc[n]
b_dees <- suppressWarnings(as.numeric(gb[3, ]))   # deescalate if y >= b_dees[n]
b_elim <- suppressWarnings(as.numeric(gb[4, ]))   # eliminate  if y >= b_elim[n]

# ---------------------------------------------------------------------
# Stage 1: BF-BOIN with backfill on the common curve
# ---------------------------------------------------------------------
run_stage1_bard <- function(s) {
  maxpat <- ncohort * cohort1 + J * n_cap + cohort1   # generous preallocation
  rec_dose <- integer(maxpat); rec_tox <- integer(maxpat); rec_eff <- integer(maxpat)
  rec_util <- numeric(maxpat)
  rec_cp1 <- integer(maxpat); rec_cp2 <- integer(maxpat); rec_cp3 <- integer(maxpat)
  nrec <- 0L
  add_rec <- function(d, p) {
    nrec <<- nrec + 1L
    rec_dose[nrec] <<- d; rec_tox[nrec] <<- p$tox; rec_eff[nrec] <<- p$eff
    rec_util[nrec] <<- p$util
    rec_cp1[nrec] <<- p$cp1; rec_cp2[nrec] <<- p$cp2; rec_cp3[nrec] <<- p$cp3
  }

  n_esc <- integer(J); y_esc <- integer(J)                    # escalation only
  bf_cal <- matrix(NA_real_, nrow = n_cap + 1L, ncol = J)     # backfill maturation calendar
  bf_evt <- matrix(0L,       nrow = n_cap + 1L, ncol = J)     # backfill DLT events
  bf_cnt <- integer(J)                                        # backfilled enrolled per dose

  first.resp <- rep(Inf, J)                                   # first response calendar per dose
  elim <- integer(J)
  d <- 1L
  clock <- 0
  n.cohort <- 1L
  pat.i <- 1L

  cal.cohort <- rep(NA_real_, cohort1)   # current escalation cohort maturation calendar
  evt.cohort <- integer(cohort1)         # current escalation cohort DLT events

  matured_backfill <- function(upto_d, now) {
    nb <- colSums(bf_cal[, 1:upto_d, drop = FALSE] < now, na.rm = TRUE)
    yb <- colSums(bf_evt[, 1:upto_d, drop = FALSE] == 1L &
                  bf_cal[, 1:upto_d, drop = FALSE] < now, na.rm = TRUE)
    list(nb = nb, yb = yb)
  }

  repeat {
    if (n.cohort > ncohort) break

    arr <- if (accrual == "poisson") rexp(1, rate = n.per.month) else runif(1, 0, 2/n.per.month)
    clock <- if (pat.i == 1L) 0 else clock + arr
    pat.i <- pat.i + 1L

    if (any(is.na(cal.cohort))) {
      pos <- max(c(0L, which(!is.na(cal.cohort)))) + 1L
      p <- draw_pat_s1(s, d)
      cal.cohort[pos] <- clock + p$dlt.time
      evt.cohort[pos] <- p$tox
      if (p$eff == 1L) first.resp[d] <- min(first.resp[d], clock + eff.window)
      add_rec(d, p)

    } else if (max(cal.cohort) < clock) {
      n_esc[d] <- n_esc[d] + cohort1
      y_esc[d] <- y_esc[d] + sum(evt.cohort)

      mb <- matured_backfill(d, clock)
      n.total <- mb$nb + n_esc[1:d]
      y.total <- mb$yb + y_esc[1:d]

      decision <- ifelse(y.total <= b_esc[pmax(1, n.total)], "escalate",
                  ifelse(y.total >= b_elim[pmax(1, n.total)], "eliminate",
                  ifelse(y.total >= b_dees[pmax(1, n.total)], "deescalate", "stay")))

      if (d == 1L) {
        if (decision[1] == "escalate") {
          if (elim[2] == 0L) d <- 2L else { if (n.total[1] >= n_stop) break }
        } else if (decision[1] == "stay") {
          if (n.total[1] >= n_stop) break
        } else if (decision[1] == "eliminate") {
          elim[1:J] <- 1L; break                    # dose 1 too toxic -> stop, no MTD
        }
        else if (n.total[1] >= n_stop) break

      } else if (decision[d] == "escalate" && all(decision[1:(d-1)] == "escalate")) {
        if (d != J && elim[d+1] == 0L) d <- d + 1L else { if (n.total[d] >= n_stop) break }
      } else if (decision[d] == "stay" && all(decision[1:(d-1)] %in% c("stay","escalate"))) {
        if (n.total[d] >= n_stop) break
      } else if (decision[d] == "deescalate" && all(decision[1:(d-1)] %in% c("stay","escalate"))) {
        d <- d - 1L
      } else if (decision[d] == "eliminate" && all(decision[1:(d-1)] %in% c("stay","escalate"))) {
        elim[d:J] <- 1L; d <- d - 1L
      } else if (decision[d] == "escalate") {
        b.star <- max(which(decision[1:(d-1)] != "escalate"))
        y.pooled <- sum(y.total[b.star:d]); n.pooled <- sum(n.total[b.star:d])
        if (y.pooled <= b_esc[n.pooled]) {
          if (d != J && elim[d+1] == 0L) d <- d + 1L else { if (n.total[d] >= n_stop) break }
        } else if (y.pooled >= b_dees[n.pooled]) {
          d.safe <- cumsum(y.total[b.star:d]) < b_dees[cumsum(n.total[b.star:d])]
          if (!any(d.safe)) {
            d <- b.star - 1L
            if (d == 0L) {
              if (y.total[1] >= b_elim[n.total[1]]) { elim[1:J] <- 1L; break } else d <- 1L
            }
          } else if (b.star - 1L + max(which(d.safe)) == d) {
            if (n.total[d] >= n_stop) break
          } else d <- b.star - 1L + max(which(d.safe))
        } else { if (n.total[d] >= n_stop) break }

      } else if (decision[d] %in% c("deescalate","eliminate","stay")) {
        if (decision[d] == "eliminate") elim[d:J] <- 1L
        b.star <- max(which(decision[1:(d-1)] %in% c("deescalate","eliminate")))
        y.pooled <- sum(y.total[b.star:d]); n.pooled <- sum(n.total[b.star:d])
        if (y.pooled <= b_esc[n.pooled]) {
          if (d != J && elim[d+1] == 0L) d <- d + 1L else { if (n.total[d] >= n_stop) break }
        } else if (y.pooled >= b_dees[n.pooled]) {
          d.safe <- cumsum(y.total[b.star:d]) < b_dees[cumsum(n.total[b.star:d])]
          if (!any(d.safe)) {
            d <- b.star - 1L
            if (d == 0L) {
              if (y.total[1] >= b_elim[n.total[1]]) { elim[1:J] <- 1L; break } else d <- 1L
            }
          } else if (b.star - 1L + max(which(d.safe)) == d) {
            if (n.total[d] >= n_stop) break
          } else d <- b.star - 1L + max(which(d.safe))
        } else { if (n.total[d] >= n_stop) break }
      }

      n.cohort <- n.cohort + 1L
      if (n.cohort > ncohort) break

      cal.cohort <- rep(NA_real_, cohort1); evt.cohort <- integer(cohort1)
      p <- draw_pat_s1(s, d)
      cal.cohort[1] <- clock + p$dlt.time
      evt.cohort[1] <- p$tox
      if (p$eff == 1L) first.resp[d] <- min(first.resp[d], clock + eff.window)
      add_rec(d, p)

    } else if (d > 1L && clock > min(first.resp[1:(d-1)])) {
      mb <- matured_backfill(d - 1L, clock)          # matured backfill n,y at doses 1:(d-1)
      yb_1 <- mb$yb; nb_1 <- mb$nb
      y.tot <- yb_1 + y_esc[1:(d-1)]                 # backfill + escalation at 1:(d-1)
      n.tot <- nb_1 + n_esc[1:(d-1)]
      nd_open <- sum(cal.cohort < clock, na.rm = TRUE)
      yd_open <- sum(evt.cohort == 1L & cal.cohort < clock, na.rm = TRUE)
      mid_y <- c(yb_1[-1], yd_open); mid_n <- c(nb_1[-1], nd_open)
      y.tot.plus <- y.tot + mid_y + y_esc[2:d]
      n.tot.plus <- n.tot + mid_n + n_esc[2:d]

      cond.1 <- cummax((first.resp[1:(d-1)] < clock) * 1L) == 1L   # activity, monotone
      cond.2 <- y.tot      <  b_dees[pmax(1, n.tot)]
      cond.3 <- y.tot.plus <  b_dees[pmax(1, n.tot.plus)]
      cond.4 <- (bf_cnt[1:(d-1)] + n_esc[1:(d-1)]) <= n_cap

      safe <- cond.2 | cond.3
      if (any(!safe)) safe[min(which(!safe)):(d-1)] <- FALSE
      open <- cond.1 & cond.4 & safe

      if (any(open)) {
        d.bf <- max(which(open))                      # highest open dose
        p <- draw_pat_s1(s, d.bf)
        slot <- bf_cnt[d.bf] + 1L
        if (slot <= n_cap + 1L) {        # cond.4 uses <=, so n_cap+1 is reachable
          bf_cal[slot, d.bf] <- clock + p$dlt.time
          bf_evt[slot, d.bf] <- p$tox
          bf_cnt[d.bf] <- slot
          if (p$eff == 1L) first.resp[d.bf] <- min(first.resp[d.bf], clock + eff.window)
          add_rec(d.bf, p)
        }
      }
    }
  }

  now <- if (all(is.na(cal.cohort))) clock else max(cal.cohort, na.rm = TRUE)
  nb <- colSums(bf_cal < now, na.rm = TRUE)
  yb <- colSums(bf_evt == 1L & bf_cal < now, na.rm = TRUE)
  n_pat <- n_esc + nb
  y_tox <- y_esc + yb

  mtd <- if (elim[1] == 1L) NA_integer_ else
         select.mtd(target = phiT, npts = n_pat, ntox = y_tox, cutoff.eli = cutoff.eli)$MTD
  if (!is.na(mtd) && mtd == 99) mtd <- NA_integer_

  stage1_dur <- now + max(dlt.window, eff.window)

  recs <- data.frame(
    dose = rec_dose[seq_len(nrec)], tox = rec_tox[seq_len(nrec)],
    eff = rec_eff[seq_len(nrec)], util = rec_util[seq_len(nrec)],
    cp1 = rec_cp1[seq_len(nrec)], cp2 = rec_cp2[seq_len(nrec)], cp3 = rec_cp3[seq_len(nrec)],
    stringsAsFactors = FALSE)

  list(n_pat = n_pat, y_tox = y_tox, mtd = mtd, recs = recs,
       B = sum(bf_cnt), stage1_dur = stage1_dur)
}

# ---------------------------------------------------------------------
# Covariate-adaptive randomization index (Minirand-equivalent)
# ---------------------------------------------------------------------
minirand_index <- function(x, m1, n_arm, car_idx, covwt, w0) {
  D <- numeric(2)
  for (i in 1:2) {
    v <- 0
    for (k in car_idx) {
      mk <- if (x[k] == 1L) m1[, k] else (n_arm - m1[, k])   # MATCHING-level counts
      mk[i] <- mk[i] + 1
      v <- v + covwt[k] * abs(mk[1] - mk[2])
    }
    nn <- n_arm; nn[i] <- nn[i] + 1
    D[i] <- v + w0 * abs(nn[1] - nn[2])
  }
  D
}

# ---------------------------------------------------------------------
# Stage 2: covariate-adaptive randomization on the candidate doses
# ---------------------------------------------------------------------
# `engine` selects the assignment rule: "car" is conditional Pocock-Simon
# minimization over the covariates named by car_idx, "sr" is 1:1 simple
# randomization (independent fair coin per patient), the BF-BOIN-SR arm.
run_stage2_car <- function(s, cand, ret, car_idx, stream, N2star, engine = "car") {
  d_high <- if (length(cand) == 2) cand[2] else NA_integer_
  nc <- length(cand)

  n_an  <- integer(nc)                        # analysis patients per candidate
  yE_an <- integer(nc)                        # efficacy successes
  cell  <- matrix(0L, nrow = nc, ncol = 4)    # utility cells (tox,eff): 00,01,10,11
  cp1_1 <- integer(nc); cp2_1 <- integer(nc); cp3_1 <- integer(nc)   # covariate 1-counts

  cix_of <- function(tox, eff) 1L + eff + 2L * tox   # 00->1, 01->2, 10->3, 11->4

  for (ci in seq_len(nc)) {
    idx <- which(ret$dose == cand[ci])
    if (length(idx)) {
      n_an[ci]  <- length(idx)
      yE_an[ci] <- sum(ret$eff[idx])
      cp1_1[ci] <- sum(ret$cp1[idx]); cp2_1[ci] <- sum(ret$cp2[idx]); cp3_1[ci] <- sum(ret$cp3[idx])
      for (k in idx) { cx <- cix_of(ret$tox[k], ret$eff[k]); cell[ci, cx] <- cell[ci, cx] + 1L }
    }
  }

  n1_seed <- n_an

  # car_idx is integer(0) on the SR path, where 1/length(car_idx) would be Inf.
  # covwt is used only by the closed-form CAR branch, so it is built only there.
  covwt <- if (engine == "car") rep(1/length(car_idx), 3) else rep(NA_real_, 3)
  s2_dur_pats <- 0L

  n1_rows <- nrow(ret)
  if (engine == "car" && CAR_ENGINE == "minirand" && !is.na(d_high) && N2star > 0L) {
    covmat_all <- rbind(as.matrix(ret[, c("cp1", "cp2", "cp3")]),
                        as.matrix(stream[, c("cp1", "cp2", "cp3")]))
    res_vec <- c(match(ret$dose, cand), rep(100L, N2star))
  }

  if (!is.na(d_high) && N2star > 0L) {
    for (q in seq_len(N2star)) {
      x1 <- stream$cp1[q]; x2 <- stream$cp2[q]; x3 <- stream$cp3[q]
      if (engine == "sr") {
        # Independent fair coin, no covariate or allocation balancing. It consumes
        # u_coin from the shared stream, so this arm sees the same patients and
        # the same outcome uniforms as the CAR arms, and calls no RNG of its own.
        target_ci <- if (stream$u_coin[q] < 0.5) 1L else 2L
      } else if (CAR_ENGINE == "minirand") {
        j <- n1_rows + q
        if (j == 1L) {
          target_ci <- sample.int(2, 1)          # nothing accrued yet; ratio 1:1
        } else {
          target_ci <- Minirand(covmat = covmat_all[, car_idx, drop = FALSE], j = j,
                                covwt = rep(1/length(car_idx), length(car_idx)),
                                ratio = c(1, 1), ntrt = 2, trtseq = c(1L, 2L),
                                method = "Range", result = res_vec, p = gamma_car)
        }
        res_vec[j] <- target_ci
      } else {
        D <- minirand_index(c(x1, x2, x3), rbind(c(cp1_1[1], cp2_1[1], cp3_1[1]),
                                                 c(cp1_1[2], cp2_1[2], cp3_1[2])),
                            n_an, car_idx, covwt, w0_alloc)
        if (D[1] == D[2]) {
          target_ci <- if (stream$u_tie[q] < 0.5) 1L else 2L
        } else {
          amin <- which.min(D)
          target_ci <- if (stream$u_coin[q] < gamma_car) amin else (3L - amin)
        }
      }
      po <- outcome_at_u(s, cand[target_ci], x1, x2, x3, stream$u_out[q])
      n_an[target_ci]  <- n_an[target_ci]  + 1L
      yE_an[target_ci] <- yE_an[target_ci] + po$eff
      cx <- cix_of(po$tox, po$eff); cell[target_ci, cx] <- cell[target_ci, cx] + 1L
      if (x1 == 1L) cp1_1[target_ci] <- cp1_1[target_ci] + 1L
      if (x2 == 1L) cp2_1[target_ci] <- cp2_1[target_ci] + 1L
      if (x3 == 1L) cp3_1[target_ci] <- cp3_1[target_ci] + 1L
    }
    s2_dur_pats <- N2star

  } else if (is.na(d_high) && N2star > 0L) {
    for (q in seq_len(N2star)) {
      po <- outcome_at_u(s, cand[1], stream$cp1[q], stream$cp2[q], stream$cp3[q],
                         stream$u_out[q])
      n_an[1]  <- n_an[1]  + 1L
      yE_an[1] <- yE_an[1] + po$eff
      cx <- cix_of(po$tox, po$eff); cell[1, cx] <- cell[1, cx] + 1L
      if (stream$cp1[q] == 1L) cp1_1[1] <- cp1_1[1] + 1L
      if (stream$cp2[q] == 1L) cp2_1[1] <- cp2_1[1] + 1L
      if (stream$cp3[q] == 1L) cp3_1[1] <- cp3_1[1] + 1L
    }
    s2_dur_pats <- N2star
  }

  imb <- rep(NA_real_, 3)
  if (!is.na(d_high) && n_an[1] > 0 && n_an[2] > 0) {
    imb[1] <- 100 * abs(cp1_1[1]/n_an[1] - cp1_1[2]/n_an[2])
    imb[2] <- 100 * abs(cp2_1[1]/n_an[1] - cp2_1[2]/n_an[2])
    imb[3] <- 100 * abs(cp3_1[1]/n_an[1] - cp3_1[2]/n_an[2])
  }

  s2_dur <- if (s2_dur_pats > 0) s2_dur_pats / n.per.month + eff.window else 0

  list(cand = cand, n_an = n_an, yE_an = yE_an, cell = cell,
       n2_alloc = n_an - n1_seed,          # fresh Stage-2 patients per candidate dose
       imb = imb, N2star = N2star, s2_dur = s2_dur,
       alloc_imb = if (!is.na(d_high)) abs(n_an[1] - n_an[2]) else NA_real_)
}

# ---- OBD selection on the combined analysis set ----------------------------
p_binom_greater <- function(y, n, p0) pbinom(y - 1, n, p0, lower.tail = FALSE)

u_cells <- c(0.3, 1.0, 0.0, 0.5)     # cells 00,01,10,11 (index = 1 + eff + 2*tox)

tox_gate <- function(yT, n) {
  yy <- yT; nn <- n
  if (length(n) == 2 && all(n > 0)) {
    if (yT[1]/n[1] > yT[2]/n[2]) { yy <- rep(sum(yT), 2); nn <- rep(sum(n), 2) }
  }
  (n > 0) & ((1 - pbeta(phiT, 1 + yy, 1 + nn - yy)) < C_T)
}
eff_gate <- function(yE, n) (n > 0) & ((1 - pbeta(phiE, 1 + yE, 1 + n - yE)) > (1 - C_E))

select_bard <- function(s, st2) {
  cand <- st2$cand; nc <- length(cand)
  n_an <- st2$n_an; yE_an <- st2$yE_an

  yT_an <- st2$cell[, 3] + st2$cell[, 4]

  tg <- tox_gate(yT_an, n_an); eg <- eff_gate(yE_an, n_an)
  adm <- which(tg & eg)

  # Every selector below works in candidate-index space (1..nc, 0 = nothing
  # selected) and is mapped to a reported dose code only at the end. As a dose
  # code PHANTOM == 0L is indistinguishable from "nothing selected", but index 1
  # and index 0 are not, and `adm` (a set of indices) never meets a dose number.

  # ---- PRIMARY: BARD utility, Dirichlet-multinomial posterior mean ----
  i_primary <- 0L
  if (length(adm) >= 1) {
    umean <- numeric(nc)
    for (ci in seq_len(nc)) {
      a4 <- a_prior/4 + st2$cell[ci, ]                 # Dirichlet(a/4,...) + cell counts
      umean[ci] <- sum((a4 / sum(a4)) * u_cells)
    }
    ua <- umean; ua[-adm] <- -Inf
    i_primary <- which.max(ua)                         # ties -> lower dose (which.max)
  }

  # ---- SECONDARY: lowest-active exact binomial (the naive-arm rule) ----
  active <- (p_binom_greater(yE_an, n_an, phiE) <= alpha) & (n_an > 0)
  i_binom <- if (any(active)) which(active)[1] else 0L

  # ---- SECONDARY: efficacy-rate noninferiority (BARD native) ----
  # The admissibility test is on the candidate INDEX, not the dose number: `adm`
  # is a set of indices, and the two coincide only when a dose number happens to
  # equal its position in cand.
  i_ni <- 0L
  if (length(adm) >= 1) {
    if (nc == 1) i_ni <- if (1L %in% adm) 1L else 0L
    else {
      pl <- if (n_an[1] > 0) yE_an[1]/n_an[1] else 0
      ph <- if (n_an[2] > 0) yE_an[2]/n_an[2] else 0
      pick <- if ((pl - ph) >= -delta_ni) 1L else 2L
      if (!(pick %in% adm)) pick <- 3L - pick
      i_ni <- if (pick %in% adm) pick else 0L
    }
  }

  # ---- SECONDARY: quasi-Bernoulli utility + dominance probability (BOIN-MEM rule) ----
  i_utilqb <- 0L
  if (length(adm) >= 1) {
    usum <- as.numeric(st2$cell %*% u_cells)
    if (nc == 1) {
      i_utilqb <- if (1L %in% adm) 1L else 0L
    } else {
      dr <- sapply(seq_len(nc), function(ci)
        rbeta(M_util, a_prior + usum[ci], b_prior + (n_an[ci] - usum[ci])))
      Pdom <- numeric(nc)
      for (ci in seq_len(nc)) Pdom[ci] <- mean(dr[, ci] >= apply(dr[, -ci, drop = FALSE], 1, max))
      Pa <- Pdom; Pa[-adm] <- -Inf
      i_utilqb <- which.max(Pa)
    }
  }

  # Efficacy estimate at the primary pick, taken on the index, so no reported
  # selection code (which may be -1) ever reaches a subscript or match().
  est <- NA_real_
  if (i_primary > 0) est <- yE_an[i_primary] / n_an[i_primary]

  # Map candidate index to reported selection code. Index 0 means "no admissible
  # candidate" and stays 0L; a phantom pick becomes -1L so the two stay
  # distinguishable downstream.
  code_pick <- function(i) {
    if (i == 0L) return(0L)
    d <- cand[i]
    if (d == PHANTOM) -1L else d
  }
  sel_primary <- code_pick(i_primary)
  sel_binom   <- code_pick(i_binom)
  sel_ni      <- code_pick(i_ni)
  sel_utilqb  <- code_pick(i_utilqb)

  door_primary <- if (sel_primary > 0) 0L else if (sel_primary == -1L) 3L else 1L

  list(sel_primary = sel_primary, sel_binom = sel_binom, sel_ni = sel_ni,
       sel_utilqb = sel_utilqb, est = est, door = door_primary)
}

# ---------------------------------------------------------------------
# Arm grid and one replicate
# ---------------------------------------------------------------------
ARMS <- data.frame(
  name   = c("bard_car3", "bard_car2", "bfboin_sr"),
  engine = c("car",       "car",       "sr"),
  car    = c("car3",      "car2",      NA_character_),
  stringsAsFactors = FALSE)
NARM <- nrow(ARMS)
# car_of is called with NA on the SR path; it returns integer(0) there, and
# car_idx is unused because the SR branch balances nothing.
car_of <- function(car) if (is.na(car)) integer(0) else if (car == "car3") 1:3 else 1:2

one_rep <- function(s) {
  st1 <- run_stage1_bard(s)
  mtd <- st1$mtd
  expo1 <- as.integer(table(factor(st1$recs$dose, levels = 1:J)))   # Stage-1 exposure per dose

  sel_p <- integer(NARM); sel_b <- integer(NARM)
  sel_n <- integer(NARM); sel_q <- integer(NARM)
  door  <- integer(NARM); est <- rep(NA_real_, NARM)
  imb   <- matrix(NA_real_, NARM, 3); N2s <- rep(NA_real_, NARM)
  dur   <- rep(NA_real_, NARM); allocimb <- rep(NA_real_, NARM)
  n_an_low <- rep(NA_real_, NARM); n_an_high <- rep(NA_real_, NARM)
  totalN <- rep(NA_real_, NARM)
  expo_arm <- matrix(0L, NARM, J)
  expo_ph  <- integer(NARM)          # phantom Stage-2 patients, per arm

  if (is.na(mtd)) {
    for (a in 1:NARM) expo_arm[a, ] <- expo1
    return(list(sel_p = rep(0L, NARM), sel_b = rep(0L, NARM),
                sel_n = rep(0L, NARM), sel_q = rep(0L, NARM),
                door = rep(2L, NARM), est = est, imb = imb, N2s = rep(0, NARM),
                dur = rep(st1$stage1_dur, NARM), allocimb = allocimb,
                n_an_low = n_an_low, n_an_high = n_an_high,
                totalN = rep(sum(expo1), NARM), mtd = mtd,
                expo_arm = expo_arm, expo_ph = expo_ph, expo1 = expo1, B = st1$B))
  }

  cand <- if (mtd == 1L) c(PHANTOM, 1L) else sort(intersect(c(mtd - 1L, mtd), 1:J))
  # Stage-1 records at the candidate doses. The phantom has none, since no
  # Stage-1 patient was treated below dose 1, so `ret` contains no rows for it
  # and its per-candidate tallies stay at zero.
  ret <- st1$recs[st1$recs$dose %in% cand[cand != PHANTOM], , drop = FALSE]

  n1_cand <- sapply(cand, function(d) if (d == PHANTOM) 0L else sum(ret$dose == d))
  # With the phantom, length(cand) is always 2 whenever mtd is not NA, so the
  # n2_single branch is unreachable as configured; it is kept so the file still
  # behaves if the phantom is disabled.
  N2star <- if (length(cand) == 2) max(0L, N2 - sum(n1_cand))
            else max(0L, n2_single - n1_cand[1])
  stream <- data.frame(
    cp1 = rbinom(N2star, 1, 0.5), cp2 = rbinom(N2star, 1, 0.5), cp3 = rbinom(N2star, 1, 0.5),
    u_out = runif(N2star), u_coin = runif(N2star), u_tie = runif(N2star))

  for (a in 1:NARM) {
    st2 <- run_stage2_car(s, cand, ret, car_of(ARMS$car[a]), stream, N2star,
                          engine = ARMS$engine[a])
    sd  <- select_bard(s, st2)

    sel_p[a] <- sd$sel_primary; sel_b[a] <- sd$sel_binom
    sel_n[a] <- sd$sel_ni;      sel_q[a] <- sd$sel_utilqb
    door[a]  <- sd$door; est[a] <- sd$est
    imb[a, ] <- st2$imb; N2s[a] <- st2$N2star
    allocimb[a] <- st2$alloc_imb
    n_an_low[a]  <- st2$n_an[1]
    n_an_high[a] <- if (length(cand) == 2) st2$n_an[2] else NA_real_
    totalN[a] <- sum(expo1) + st2$N2star     # Stage-1 enrolled (all doses) + fresh Stage-2
    dur[a] <- st1$stage1_dur + st2$s2_dur
    expo_arm[a, ] <- expo1
    # A phantom code cannot subscript the length-J exposure row, so the real and
    # phantom candidates are split and the phantom is tallied separately.
    real <- cand != PHANTOM
    expo_arm[a, cand[real]] <- expo_arm[a, cand[real]] + st2$n2_alloc[real]
    expo_ph[a] <- if (any(!real)) st2$n2_alloc[!real] else 0L
  }

  list(sel_p = sel_p, sel_b = sel_b, sel_n = sel_n, sel_q = sel_q, door = door, est = est,
       imb = imb, N2s = N2s, dur = dur, allocimb = allocimb,
       n_an_low = n_an_low, n_an_high = n_an_high, totalN = totalN,
       mtd = mtd, expo_arm = expo_arm, expo_ph = expo_ph, expo1 = expo1, B = st1$B)
}

# ---------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------
run_bard <- function(cfg) {
  set.seed(seed)
  full_rows <- list(); si <- 0
  alloc_rows <- list(); ai <- 0

  for (s in 1:15) {
    ob <- trueOBD[s]; mt <- trueMTD[s]
    tox_over <- which(S_tox[s, ] > phiT)

    sp <- matrix(NA_integer_, NSIM, NARM); sb <- matrix(NA_integer_, NSIM, NARM)
    sn <- matrix(NA_integer_, NSIM, NARM); sq <- matrix(NA_integer_, NSIM, NARM)
    dr <- matrix(NA_integer_, NSIM, NARM)
    es <- matrix(NA_real_, NSIM, NARM)
    i1 <- array(NA_real_, c(NSIM, NARM, 3))
    n2s <- matrix(NA_real_, NSIM, NARM); dur <- matrix(NA_real_, NSIM, NARM)
    aim <- matrix(NA_real_, NSIM, NARM); tN <- matrix(NA_real_, NSIM, NARM)
    nlo <- matrix(NA_real_, NSIM, NARM); nhi <- matrix(NA_real_, NSIM, NARM)
    expo <- array(0L, c(NSIM, NARM, J))          # per-arm per-dose exposure (Stage 1 + Stage 2)
    expo_ph_mat <- matrix(0L, NSIM, NARM)        # per-arm phantom Stage-2 patients
    expo1_mat <- matrix(0L, NSIM, J)             # Stage-1-only exposure (shared across arms)
    mtdv <- rep(NA_integer_, NSIM); Bv <- rep(NA_real_, NSIM)

    for (r in 1:NSIM) {
      o <- one_rep(s)
      sp[r, ] <- o$sel_p; sb[r, ] <- o$sel_b; sn[r, ] <- o$sel_n; sq[r, ] <- o$sel_q
      dr[r, ] <- o$door
      es[r, ] <- o$est; i1[r, , ] <- o$imb; n2s[r, ] <- o$N2s; dur[r, ] <- o$dur
      aim[r, ] <- o$allocimb; tN[r, ] <- o$totalN; nlo[r, ] <- o$n_an_low; nhi[r, ] <- o$n_an_high
      expo[r, , ] <- o$expo_arm; expo_ph_mat[r, ] <- o$expo_ph; expo1_mat[r, ] <- o$expo1
      mtdv[r] <- o$mtd; Bv[r] <- o$B
    }

    # ---- allocation, per arm ----
    # Phantom Stage-2 patients are folded into N_used and below_n, since the
    # phantom sits below dose 1 and so below any true OBD, and are broken out as
    # phantom_n / phantom_pct. at_n, over_n and aboveMTD_n are untouched.
    for (a in 1:NARM) {
      ea_    <- expo[, a, , drop = FALSE]; dim(ea_) <- c(NSIM, J)
      phv    <- expo_ph_mat[, a]
      Nm     <- mean(rowSums(ea_) + phv)
      belowN <- if (is.na(ob)) NA else
                mean(rowSums(ea_[, seq_len(J) < ob, drop = FALSE]) + phv)
      atN    <- if (is.na(ob)) NA else mean(ea_[, ob])
      overN  <- if (is.na(ob)) NA else mean(rowSums(ea_[, seq_len(J) > ob, drop = FALSE]))
      abvN   <- if (length(tox_over) == 0) 0 else mean(rowSums(ea_[, tox_over, drop = FALSE]))
      ai <- ai + 1
      alloc_rows[[ai]] <- data.frame(
        Scn = s, family = family[s], method = ARMS$name[a],
        trueOBD = ob, trueMTD = mt,
        N_used = Nm,                            # Stage 1 + Stage 2, all doses
        phantom_n = mean(phv), phantom_pct = 100 * mean(phv) / Nm,
        Nesc_stage1 = mean(rowSums(expo1_mat)), # Stage-1 only (escalation + backfill)
        below_n = belowN, at_n = atN, over_n = overN, aboveMTD_n = abvN,
        B_backfill = mean(Bv),
        PCS_MTD = if (is.na(mt)) NA else 100 * mean(!is.na(mtdv) & mtdv == mt),
        duration = mean(dur[, a], na.rm = TRUE),
        stringsAsFactors = FALSE)
    }

    # ---- per-arm selection + BARD-specific metrics ----
    for (a in 1:NARM) {
      sa <- sp[, a]; da <- dr[, a]; ea <- es[, a]
      selected <- !is.na(sa) & sa > 0
      PCS   <- if (is.na(ob)) NA else 100 * mean(sa == ob, na.rm = TRUE)
      Under <- if (is.na(ob)) NA else 100 * mean(selected & sa < ob)
      Over  <- if (is.na(ob)) NA else 100 * mean(selected & sa > ob)
      # sel = -1 (phantom pick) falls out of `selected`, so PCS / UnderOBD /
      # OverOBD / ToxSel / sel_d1..d5 are unaffected and it counts toward noOBD.
      # noOBD therefore exceeds noMTD + noEff by phantomSel.
      noOBD <- 100 * mean(!selected)
      noMTD <- 100 * mean(da == 2)
      noEff <- 100 * mean(da == 1)
      phantomSel <- 100 * mean(sa == -1, na.rm = TRUE)
      tox_flag <- logical(NSIM); tox_flag[selected] <- S_tox[s, sa[selected]] > phiT
      ToxSel <- 100 * mean(tox_flag); any_sel <- 100 * mean(selected)
      PCS_binom  <- if (is.na(ob)) NA else 100 * mean(sb[, a] == ob, na.rm = TRUE)
      PCS_ni     <- if (is.na(ob)) NA else 100 * mean(sn[, a] == ob, na.rm = TRUE)
      PCS_utilqb <- if (is.na(ob)) NA else 100 * mean(sq[, a] == ob, na.rm = TRUE)
      # <= 0 rather than == 0, so a phantom pick (-1) counts as "no real dose
      # selected", matching how noOBD treats the primary selector.
      noOBD_binom  <- 100 * mean(sb[, a] <= 0)
      noOBD_utilqb <- 100 * mean(sq[, a] <= 0)
      if (any(selected)) { Est <- mean(ea[selected], na.rm = TRUE); True <- mean(S_eff[s, sa[selected]]) }
      else { Est <- NA; True <- NA }
      d1 <- 100*mean(sa==1,na.rm=TRUE); d2 <- 100*mean(sa==2,na.rm=TRUE)
      d3 <- 100*mean(sa==3,na.rm=TRUE); d4 <- 100*mean(sa==4,na.rm=TRUE); d5 <- 100*mean(sa==5,na.rm=TRUE)

      si <- si + 1
      full_rows[[si]] <- data.frame(
        Scn = s, family = family[s],
        method = ARMS$name[a], engine = ARMS$engine[a], car = ARMS$car[a],
        trueOBD = ob, trueMTD = mt,
        PCS = PCS, UnderOBD = Under, OverOBD = Over,
        noOBD = noOBD, noMTD = noMTD, noEff = noEff, phantomSel = phantomSel,
        ToxSel = ToxSel, select_any = any_sel,
        PCS_binom = PCS_binom, PCS_ni = PCS_ni, PCS_utilqb = PCS_utilqb,
        noOBD_binom = noOBD_binom, noOBD_utilqb = noOBD_utilqb,
        Est = Est, True = True, Boost = Est - True,
        totalN = mean(tN[, a]), N2star = mean(n2s[, a], na.rm = TRUE),
        n_an_low = mean(nlo[, a], na.rm = TRUE), n_an_high = mean(nhi[, a], na.rm = TRUE),
        alloc_imb = mean(aim[, a], na.rm = TRUE),
        imb_CP1 = mean(i1[, a, 1], na.rm = TRUE),
        imb_CP2 = mean(i1[, a, 2], na.rm = TRUE),
        imb_CP3 = mean(i1[, a, 3], na.rm = TRUE),
        duration = mean(dur[, a], na.rm = TRUE),
        sel_d1 = d1, sel_d2 = d2, sel_d3 = d3, sel_d4 = d4, sel_d5 = d5,
        stringsAsFactors = FALSE)
    }
    cat(sprintf("  [bard] scenario %2d done\n", s))
  }

  full  <- do.call(rbind, full_rows)
  alloc <- do.call(rbind, alloc_rows)

  rnd <- function(x, k) round(x, k)

  sel_tab <- data.frame(
    Scn = full$Scn, family = full$family, method = full$method,
    trueOBD = full$trueOBD, trueMTD = full$trueMTD,
    PCS = rnd(full$PCS,1), UnderOBD = rnd(full$UnderOBD,1), OverOBD = rnd(full$OverOBD,1),
    noOBD = rnd(full$noOBD,1), noMTD = rnd(full$noMTD,1), noEff = rnd(full$noEff,1),
    phantomSel = rnd(full$phantomSel,1),
    ToxSel = rnd(full$ToxSel,1),
    PCS_binom = rnd(full$PCS_binom,1), PCS_ni = rnd(full$PCS_ni,1),
    PCS_utilqb = rnd(full$PCS_utilqb,1),
    noOBD_binom = rnd(full$noOBD_binom,1), noOBD_utilqb = rnd(full$noOBD_utilqb,1),
    stringsAsFactors = FALSE)

  dist_tab <- data.frame(
    Scn = full$Scn, method = full$method,
    d1 = rnd(full$sel_d1,1), d2 = rnd(full$sel_d2,1), d3 = rnd(full$sel_d3,1),
    d4 = rnd(full$sel_d4,1), d5 = rnd(full$sel_d5,1), noOBD = rnd(full$noOBD,1),
    stringsAsFactors = FALSE)

  alloc_tab <- data.frame(
    Scn = alloc$Scn, family = alloc$family, method = alloc$method,
    trueOBD = alloc$trueOBD, trueMTD = alloc$trueMTD,
    N_used = rnd(alloc$N_used,1),
    phantom_n = rnd(alloc$phantom_n,1), phantom_pct = rnd(alloc$phantom_pct,1),
    Nesc_stage1 = rnd(alloc$Nesc_stage1,1),
    below_n = rnd(alloc$below_n,1), at_n = rnd(alloc$at_n,1),
    over_n = rnd(alloc$over_n,1), aboveMTD_n = rnd(alloc$aboveMTD_n,1),
    B_backfill = rnd(alloc$B_backfill,1), PCS_MTD = rnd(alloc$PCS_MTD,1),
    duration = rnd(alloc$duration,1),
    stringsAsFactors = FALSE)

  metrics_tab <- data.frame(
    Scn = full$Scn, method = full$method,
    totalN = rnd(full$totalN,1), N2star = rnd(full$N2star,1),
    n_an_low = rnd(full$n_an_low,1), n_an_high = rnd(full$n_an_high,1),
    alloc_imb = rnd(full$alloc_imb,2),
    imb_CP1 = rnd(full$imb_CP1,2), imb_CP2 = rnd(full$imb_CP2,2), imb_CP3 = rnd(full$imb_CP3,2),
    duration = rnd(full$duration,1),
    stringsAsFactors = FALSE)

  od <- cfg$outdir
  write.csv(full, file.path(od, "setting1_bard_results_full.csv"),        row.names = FALSE)
  write.csv(sel_tab,     file.path(od, "setting1_bard_selection.csv"),  row.names = FALSE)
  write.csv(dist_tab,    file.path(od, "setting1_bard_seldist.csv"),    row.names = FALSE)
  write.csv(alloc_tab,   file.path(od, "setting1_bard_allocation.csv"), row.names = FALSE)
  write.csv(metrics_tab, file.path(od, "setting1_bard_metrics.csv"),   row.names = FALSE)

  cat("[bard] wrote setting1_bard_{selection,seldist,allocation,metrics}.csv and results_full to", od, "\n")

  invisible(list(full = full, sel = sel_tab, dist = dist_tab, alloc = alloc_tab, metrics = metrics_tab))
}
