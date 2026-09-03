# =====================================================================
# Setting 1: naive / pooling family (method class "naive")
# =====================================================================
# Arms (5): pool_all, stage2only, pool_0.3, combined_p_all, combined_p_0.3.
#
# Requires setting1_config.R to be sourced first (S_tox, S_eff, beta0,
# draw_gumbel, util_obs, arrival_gap, p_binom_greater, pick_lowest_active,
# trueOBD, trueMTD, family, the shared constants and calendar, the phantom block
# and the PHANTOM / tox_at / b0_at / n1_at accessors).
#
# Phantom MTD-1 dose: when Stage 1 returns mtd == 1 the candidate set becomes
# c(PHANTOM, 1) instead of a lone dose 1, so every replicate has two candidates.
# The phantom is a decoy at 80 percent of dose 1's toxicity and efficacy; picking
# it is a selection failure, coded sel = -1, door = 3, and reported as
# phantomSel. Its Stage-2 patients are counted in N_used and below_n and broken
# out as phantom_n / phantom_pct.
# =====================================================================

# ---- method-specific constants ----
cohort <- cohort1     # Stage-1 cohort size (config exposes it as cohort1)
n2     <- 20          # Stage-2 patients per candidate dose (fixed expansion)

# Stage 1 is the shared routine run_stage1_shared(s) in setting1_config.R
# section 7, together with the boundary vectors, the per-dose cap s1 and the
# per-replicate seed stage1_seed().

# pick_lowest_active returns 0L both for "nothing active" and for a phantom pick
# (PHANTOM == 0L), so the phantom case is remapped to -1L to keep the two
# distinguishable downstream.
resolve_pick <- function(cand, active) {
  if (!any(active)) return(0L)
  d <- cand[which(active)[1]]
  if (d == PHANTOM) -1L else d
}

# ---------------------------------------------------------------------
# One replicate: shared trial, then the five arms' tests
# ---------------------------------------------------------------------
# Arm order: 1 pool_all 2 stage2only 3 pool_0.3 4 combined_p_all 5 combined_p_0.3
# door: 0 selected | 1 MTD found, no candidate active | 2 no MTD | 3 phantom picked
one_rep <- function(s) {
  st1  <- run_stage1_shared(s)
  mtd  <- st1$mtd
  expo <- st1$n_pat

  if (is.na(mtd)) {
    return(list(sel = rep(0L, 5), door = rep(2L, 5),
                est = rep(NA_real_, 5), mtd = mtd, expo = expo, expo_ph = 0L,
                dur1 = st1$stage1_dur, dur2 = 0, dur = st1$stage1_dur,
                reached_s2 = FALSE))
  }

  # When Stage 1 returns mtd == 1 there is no real dose below it, so the phantom
  # is inserted. PHANTOM is a dose code, never an index.
  cand <- if (mtd == 1L) c(PHANTOM, 1L) else sort(intersect(c(mtd - 1L, mtd), 1:J))
  is_ph <- cand == PHANTOM

  # Phantom patients count toward N_used but cannot be recorded in the length-J
  # per-dose vector, so they are tallied separately and folded in by the driver.
  expo[cand[!is_ph]] <- expo[cand[!is_ph]] + n2
  expo_ph <- if (any(is_ph)) n2 else 0L

  dur2 <- (n2 * length(cand)) / n.per.month + eff.window

  # ---- Stage 2: n2 fresh patients at each candidate (common curve) ----
  y2 <- integer(length(cand))
  for (ci in seq_along(cand)) {
    d   <- cand[ci]
    cp1 <- rbinom(n2, 1, 0.5); cp2 <- rbinom(n2, 1, 0.5); cp3 <- rbinom(n2, 1, 0.5)
    pE  <- plogis(b0_at(s, d) + bP1 * cp1 + bP2 * cp2 + bP3 * cp3)
    oc  <- draw_gumbel(rep(tox_at(s, d), n2), pE, psi)
    y2[ci] <- sum(oc$eff)
  }

  # ---- borrowed Stage-1 at candidate doses: all and the random poolable subset ----
  # Routed through n1_at: the phantom has no Stage-1 patients under any rule.
  if (is.null(st1$n_pool) || is.null(st1$y_pool))
    stop("run_stage1_shared() returned no n_pool / y_pool: setting1_config.R is out of sync")
  nA <- sapply(cand, function(d) n1_at(st1$n_all,  d))
  yA <- sapply(cand, function(d) n1_at(st1$y_all,  d))
  nH <- sapply(cand, function(d) n1_at(st1$n_pool, d))
  yH <- sapply(cand, function(d) n1_at(st1$y_pool, d))

  # ---- (a) pool all ----
  p_a   <- p_binom_greater(yA + y2, nA + n2, phiE)
  sel_a <- resolve_pick(cand, p_a <= alpha)

  # ---- (b) Stage-2 only ----
  p_b   <- p_binom_greater(y2, n2, phiE)
  sel_b <- resolve_pick(cand, p_b <= alpha)

  # ---- (c) pool the random poolable subset ----
  p_c   <- p_binom_greater(yH + y2, nH + n2, phiE)
  sel_c <- resolve_pick(cand, p_c <= alpha)

  # ---- (d) combined-p, all (Fisher, df = 4; active if p < alpha) ----
  p_s2   <- p_binom_greater(y2, n2, phiE)
  p_s1A  <- p_binom_greater(yA, nA, phiE)
  Xd     <- -2 * (log(p_s1A) + log(p_s2))
  p_d    <- pchisq(Xd, df = 4, lower.tail = FALSE)
  sel_d  <- resolve_pick(cand, p_d < alpha)

  # ---- (e) combined-p, poolable subset ----
  p_s1H  <- p_binom_greater(yH, nH, phiE)
  Xe     <- -2 * (log(p_s1H) + log(p_s2))
  p_e    <- pchisq(Xe, df = 4, lower.tail = FALSE)
  sel_e  <- resolve_pick(cand, p_e < alpha)

  # ---- winner's-boost efficacy estimate at the selected dose ----
  # Guarded on > 0 rather than != 0, so a phantom pick (-1) reports est = NA
  # instead of reaching match().
  est_a <- if (sel_a > 0) { i <- match(sel_a, cand); (yA[i] + y2[i]) / (nA[i] + n2) } else NA
  est_b <- if (sel_b > 0) y2[match(sel_b, cand)] / n2 else NA
  est_c <- if (sel_c > 0) { i <- match(sel_c, cand); (yH[i] + y2[i]) / (nH[i] + n2) } else NA
  est_d <- if (sel_d > 0) y2[match(sel_d, cand)] / n2 else NA   # combined-p: Stage-2 readout
  est_e <- if (sel_e > 0) y2[match(sel_e, cand)] / n2 else NA

  sel  <- c(sel_a, sel_b, sel_c, sel_d, sel_e)
  door <- ifelse(sel > 0, 0L, ifelse(sel == -1L, 3L, 1L))
  list(sel = sel, door = door,
       est = c(est_a, est_b, est_c, est_d, est_e),
       mtd = mtd, expo = expo, expo_ph = expo_ph,
       dur1 = st1$stage1_dur, dur2 = dur2, dur = st1$stage1_dur + dur2,
       reached_s2 = TRUE)
}

# ---------------------------------------------------------------------
# Driver: loop scenarios x replicates, accumulate, summarize, write
# ---------------------------------------------------------------------
run_naive <- function(cfg) {
  set.seed(seed)
  # Arm labels, as written to the `method` column of every CSV. If pool_frac in
  # setting1_config.R changes, the two fraction-bearing labels must change too.
  arms    <- c("pool_all", "stage2only", "pool_0.3", "combined_p_all", "combined_p_0.3")
  NARM    <- length(arms)

  sel_rows   <- list(); si <- 0
  alloc_rows <- list(); ai <- 0

  for (s in 1:15) {
    ob <- trueOBD[s]; mt <- trueMTD[s]
    tox_over_doses <- which(S_tox[s, ] > phiT)

    sel  <- matrix(NA_integer_, NSIM, NARM)
    door <- matrix(NA_integer_, NSIM, NARM)
    est  <- matrix(NA_real_,    NSIM, NARM)
    mtdv <- rep(NA_integer_, NSIM)
    expo <- matrix(0L, NSIM, J)
    expo_ph_vec <- integer(NSIM)      # phantom Stage-2 patients (outside the J doses)
    durv  <- rep(NA_real_, NSIM)
    dur1v <- rep(NA_real_, NSIM)
    dur2v <- rep(NA_real_, NSIM)
    s2v   <- logical(NSIM)

    for (r in 1:NSIM) {
      set.seed(stage1_seed(s, r))   # CRN: identical Stage 1 across families
      o <- one_rep(s)
      sel[r, ]  <- o$sel
      door[r, ] <- o$door
      est[r, ]  <- o$est
      mtdv[r]   <- o$mtd
      expo[r, ] <- o$expo
      expo_ph_vec[r] <- o$expo_ph
      durv[r]   <- o$dur
      dur1v[r]  <- o$dur1
      dur2v[r]  <- o$dur2
      s2v[r]    <- o$reached_s2
    }

    # ---- shared allocation and escalation quality ----
    # Phantom patients are folded into N_used and below_n, so N_used means every
    # patient enrolled in the simulated trial; they are also broken out in
    # phantom_n / phantom_pct. The phantom sits below dose 1, hence below any
    # true OBD, so at_obd, over and aboveMTD are untouched.
    N_used   <- rowSums(expo) + expo_ph_vec
    below    <- if (is.na(ob)) rep(NA, NSIM) else
                rowSums(expo[, seq_len(J) <  ob, drop = FALSE]) + expo_ph_vec
    at_obd   <- if (is.na(ob)) rep(NA, NSIM) else expo[, ob]
    over     <- if (is.na(ob)) rep(NA, NSIM) else rowSums(expo[, seq_len(J) >  ob, drop = FALSE])
    aboveMTD <- if (length(tox_over_doses) == 0) rep(0, NSIM) else
                rowSums(expo[, tox_over_doses, drop = FALSE])
    Nm <- mean(N_used)

    ai <- ai + 1
    alloc_rows[[ai]] <- data.frame(
      Scn = s, family = family[s], trueOBD = ob, trueMTD = mt,
      N_used   = Nm,
      phantom_n = mean(expo_ph_vec), phantom_pct = 100 * mean(expo_ph_vec) / Nm,
      below_n  = mean(below),   below_pct = 100 * mean(below)   / Nm,
      at_n     = mean(at_obd),  at_pct    = 100 * mean(at_obd)  / Nm,
      over_n   = mean(over),    over_pct  = 100 * mean(over)    / Nm,
      aboveMTD_n = mean(aboveMTD), aboveMTD_pct = 100 * mean(aboveMTD) / Nm,
      PCS_MTD  = if (is.na(mt)) NA else 100 * mean(!is.na(mtdv) & mtdv == mt),
      duration   = mean(durv),
      dur_stage1 = mean(dur1v),
      dur_stage2 = mean(dur2v),
      dur_med    = median(durv),
      dur_q90    = unname(quantile(durv, 0.90)),
      dur_ifS2   = if (any(s2v)) mean(durv[s2v]) else NA,
      pct_noS2   = 100 * mean(!s2v),
      stringsAsFactors = FALSE
    )

    # ---- per-arm selection and efficacy ----
    for (a in 1:NARM) {
      sa <- sel[, a]; da <- door[, a]; ea <- est[, a]
      selected <- !is.na(sa) & sa > 0

      # sel = -1 (phantom pick) falls out of `selected`, so PCS / UnderOBD /
      # OverOBD / ToxSel / sel_d1..d5 are unaffected and it counts toward noOBD.
      # noOBD therefore exceeds noMTD + noEff by phantomSel (door 3).
      PCS   <- if (is.na(ob)) NA else 100 * mean(sa == ob, na.rm = TRUE)
      Under <- if (is.na(ob)) NA else 100 * mean(selected & sa < ob)
      Over  <- if (is.na(ob)) NA else 100 * mean(selected & sa > ob)
      noOBD <- 100 * mean(!selected)
      noMTD <- 100 * mean(da == 2, na.rm = TRUE)
      noEff <- 100 * mean(da == 1, na.rm = TRUE)
      phantomSel <- 100 * mean(sa == -1, na.rm = TRUE)

      tox_flag <- logical(NSIM)
      tox_flag[selected] <- S_tox[s, sa[selected]] > phiT
      ToxSel <- 100 * mean(tox_flag)

      any_sel <- 100 * mean(selected)

      if (any(selected)) {
        Est  <- mean(ea[selected])
        True <- mean(S_eff[s, sa[selected]])
      } else { Est <- NA; True <- NA }

      si <- si + 1
      sel_rows[[si]] <- data.frame(
        Scn = s, family = family[s], method = arms[a],
        trueOBD = ob, trueMTD = mt,
        PCS = PCS, UnderOBD = Under, OverOBD = Over,
        noOBD = noOBD, noMTD = noMTD, noEff = noEff, phantomSel = phantomSel,
        ToxSel = ToxSel,
        select_any = any_sel,
        Est = Est, True = True, Boost = Est - True,
        sel_d1 = 100 * mean(sa == 1, na.rm = TRUE), sel_d2 = 100 * mean(sa == 2, na.rm = TRUE),
        sel_d3 = 100 * mean(sa == 3, na.rm = TRUE), sel_d4 = 100 * mean(sa == 4, na.rm = TRUE),
        sel_d5 = 100 * mean(sa == 5, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
    cat(sprintf("  [naive] scenario %2d done\n", s))
  }

  res   <- do.call(rbind, sel_rows)
  alloc <- do.call(rbind, alloc_rows)

  # ---- output tables ----
  # Same vector as `arms`, so display order cannot drift from arm order.
  disp_methods <- arms
  res$method <- factor(res$method, levels = disp_methods)
  disp <- res[order(res$Scn, res$method), ]

  rnd <- function(x, k) round(x, k)

  sel_tab <- data.frame(
    Scn = disp$Scn, family = disp$family,
    trueOBD = disp$trueOBD, trueMTD = disp$trueMTD, method = as.character(disp$method),
    PCS = rnd(disp$PCS, 1), UnderOBD = rnd(disp$UnderOBD, 1), OverOBD = rnd(disp$OverOBD, 1),
    noOBD = rnd(disp$noOBD, 1), noMTD = rnd(disp$noMTD, 1), noEff = rnd(disp$noEff, 1),
    phantomSel = rnd(disp$phantomSel, 1),
    ToxSel = rnd(disp$ToxSel, 1),
    stringsAsFactors = FALSE
  )

  dist_tab <- data.frame(
    Scn = disp$Scn, Method = as.character(disp$method),
    d1 = rnd(disp$sel_d1, 1), d2 = rnd(disp$sel_d2, 1), d3 = rnd(disp$sel_d3, 1),
    d4 = rnd(disp$sel_d4, 1), d5 = rnd(disp$sel_d5, 1),
    noOBD = rnd(disp$noOBD, 1),
    stringsAsFactors = FALSE
  )

  alloc_tab <- data.frame(
    Scn = alloc$Scn, family = alloc$family, trueOBD = alloc$trueOBD, trueMTD = alloc$trueMTD,
    N_used = rnd(alloc$N_used, 1),
    phantom_n = rnd(alloc$phantom_n, 1), phantom_pct = rnd(alloc$phantom_pct, 1),
    below_n = rnd(alloc$below_n, 1), below_pct = rnd(alloc$below_pct, 1),
    at_n = rnd(alloc$at_n, 1),       at_pct = rnd(alloc$at_pct, 1),
    over_n = rnd(alloc$over_n, 1),   over_pct = rnd(alloc$over_pct, 1),
    aboveMTD_n = rnd(alloc$aboveMTD_n, 1), aboveMTD_pct = rnd(alloc$aboveMTD_pct, 1),
    PCS_MTD = rnd(alloc$PCS_MTD, 1),
    duration = rnd(alloc$duration, 1),
    dur_stage1 = rnd(alloc$dur_stage1, 1), dur_stage2 = rnd(alloc$dur_stage2, 1),
    dur_med = rnd(alloc$dur_med, 1), dur_q90 = rnd(alloc$dur_q90, 1),
    dur_ifS2 = rnd(alloc$dur_ifS2, 1), pct_noS2 = rnd(alloc$pct_noS2, 1),
    stringsAsFactors = FALSE
  )

  od <- cfg$outdir
  write.csv(res, file.path(od, "setting1_naive_results_full.csv"),     row.names = FALSE)
  write.csv(sel_tab,   file.path(od, "setting1_naive_selection.csv"),  row.names = FALSE)
  write.csv(dist_tab,  file.path(od, "setting1_naive_seldist.csv"),    row.names = FALSE)
  write.csv(alloc_tab, file.path(od, "setting1_naive_allocation.csv"), row.names = FALSE)

  cat("[naive] wrote setting1_naive_{selection,seldist,allocation}.csv and results_full to", od, "\n")

  invisible(list(res = res, sel = sel_tab, dist = dist_tab, alloc = alloc_tab))
}
