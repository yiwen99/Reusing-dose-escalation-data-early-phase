# =====================================================================
# Setting 2: naive / partial-pooling family (method class "naive")
# =====================================================================
# Arms (7): a_stage2only, b_partialpool_clean, b_partialpool_mis,
#           c_combinedp_clean, c_combinedp_mis, pool_all_exact,
#           pool_all_combinedp.
#
# Requires setting2_config.R sourced first (S1/S2 tox/eff, beta0_i/beta0_ii,
# draw_gumbel, p_binom_greater, pick_lowest_active, arrival_gap, trueOBD,
# trueMTD, family, rho, lambda_grid, the shared constants and calendar, the
# phantom block and the PHANTOM / tox_at / b0_at / n1_at accessors).
#
# Phantom MTD-1 dose: when Stage 1 returns mtd == 1 the candidate set becomes
# c(PHANTOM, 1) instead of a lone dose 1, so every replicate has two candidates.
# The phantom is a decoy at 80 percent of dose 1's component-ii toxicity and
# efficacy; picking it is a selection failure, coded sel = -1, door = 3, and
# reported as phantomSel. It has no Stage-1 patients under any matching rule,
# since no Stage-1 patient was treated below dose 1. Its Stage-2 patients are
# counted in N_used and below_n and broken out as phantom_n / phantom_pct.
# =====================================================================

# ---- method-specific constants ----
cohort <- cohort1     # Stage-1 cohort size (config exposes it as cohort1)
n2     <- 20          # Stage-2 patients per candidate dose (fixed expansion)

# Stage 1 is the shared routine run_stage1_shared(s, lam) in setting2_config.R
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
# One replicate: shared trial, then all seven arms' tests
# ---------------------------------------------------------------------
# Arm order: 1 a_stage2only 2 b_partialpool_clean 3 b_partialpool_mis
#            4 c_combinedp_clean [NA at lambda 0] 5 c_combinedp_mis [NA at lambda 0]
#            6 pool_all_exact 7 pool_all_combinedp
# door: 0 selected | 1 MTD found, no candidate active | 2 no MTD | 3 phantom picked
one_rep <- function(s, lam) {
  st1  <- run_stage1_shared(s, lam)
  mtd  <- st1$mtd
  expo <- st1$n_pat                                # Stage-1 exposure per dose

  if (is.na(mtd)) {
    cc <- if (lam > 0) 0 else NA
    dd <- if (lam > 0) 2 else NA
    sel  <- c(0, 0, 0, cc, cc, 0, 0)
    door <- c(2, 2, 2, dd, dd, 2, 2)
    est  <- rep(NA_real_, 7)
    return(list(sel = sel, door = door, est = est, mtd = mtd, expo = expo,
                expo_ph = 0L,
                dur1 = st1$stage1_dur, dur2 = 0, dur = st1$stage1_dur,
                reached_s2 = FALSE))
  }

  # When Stage 1 returns mtd == 1 there is no real dose below it, so the phantom
  # is inserted. PHANTOM is a dose code, never an index.
  cand <- if (mtd == 1L) c(PHANTOM, 1L) else sort(intersect(c(mtd - 1L, mtd), 1:J))
  is_ph <- cand == PHANTOM

  # Phantom patients count toward N_used but cannot be recorded in the length-J
  # per-dose vector, so they are tallied separately and folded in by the driver.
  expo[cand[!is_ph]] <- expo[cand[!is_ph]] + n2   # Stage-2 expansion, n2 per candidate
  expo_ph <- if (any(is_ph)) n2 else 0L

  dur2 <- (n2 * length(cand)) / n.per.month + eff.window

  # ---- Stage 2: n2 fresh pure-ii patients at each candidate (joint draw) ----
  y2 <- integer(length(cand))
  for (ci in seq_along(cand)) {
    d   <- cand[ci]
    cp1 <- rbinom(n2, 1, 0.5)
    cp2 <- rbinom(n2, 1, 0.5)
    cp3 <- rbinom(n2, 1, 0.5)
    pE  <- plogis(b0_at(s, d) + bP1 * cp1 + bP2 * cp2 + bP3 * cp3)
    oc  <- draw_gumbel(rep(tox_at(s, d), n2), pE, psi)
    y2[ci] <- sum(oc$eff)
  }

  # ---- retained Stage-1 (ii) at the candidate doses: clean and mis ----
  # Routed through n1_at. The phantom contributes zero to every pool: no Stage-1
  # patient was treated there, so none can be retained under any matching rule.
  n1c <- sapply(cand, function(d) n1_at(st1$n_ii_c, d))
  y1c <- sapply(cand, function(d) n1_at(st1$y_ii_c, d))
  n1m <- sapply(cand, function(d) n1_at(st1$n_ii_m, d))
  y1m <- sapply(cand, function(d) n1_at(st1$y_ii_m, d))
  # ---- ALL Stage-1 (full mixture) at the candidate doses (~1-lambda contam) ----
  nAll <- sapply(cand, function(d) n1_at(st1$n_pat, d))
  yAll <- sapply(cand, function(d) n1_at(st1$y_all, d))

  # ---- (a) Stage-2 only ----
  p_a   <- p_binom_greater(y2, n2, phiE)
  sel_a <- resolve_pick(cand, p_a <= alpha)

  # ---- (b) partial pooling, clean ----
  p_bc   <- p_binom_greater(y2 + y1c, n2 + n1c, phiE)
  sel_bc <- resolve_pick(cand, p_bc <= alpha)

  # ---- (b') partial pooling, mis-specified ----
  p_bm   <- p_binom_greater(y2 + y1m, n2 + n1m, phiE)
  sel_bm <- resolve_pick(cand, p_bm <= alpha)

  # ---- (d) pool-all, EXACT: pool the full Stage-1 mixture into the Stage-2 test ----
  p_pa   <- p_binom_greater(y2 + yAll, n2 + nAll, phiE)
  sel_pa <- resolve_pick(cand, p_pa <= alpha)

  p_s2 <- p_a   # Stage-2-only p-value; reused by every Fisher combination below

  # ---- (c/c') Combined-p, clean and mis: Fisher combine (df = 4) ----
  # The lam == 0 branch returns NA for these two arms.
  if (lam > 0) {
    p_s1c  <- p_binom_greater(y1c, n1c, phiE)
    Xc     <- -2 * (log(p_s1c) + log(p_s2))
    p_cc   <- pchisq(Xc, df = 4, lower.tail = FALSE)
    sel_cc <- resolve_pick(cand, p_cc < alpha)

    p_s1m  <- p_binom_greater(y1m, n1m, phiE)
    Xm     <- -2 * (log(p_s1m) + log(p_s2))
    p_cm   <- pchisq(Xm, df = 4, lower.tail = FALSE)
    sel_cm <- resolve_pick(cand, p_cm < alpha)
  } else {
    sel_cc <- NA_integer_; sel_cm <- NA_integer_
  }

  # ---- (e) Combined-p, ALL: Fisher on the full Stage-1 mixture p and Stage-2 p ----
  p_s1all <- p_binom_greater(yAll, nAll, phiE)
  Xall    <- -2 * (log(p_s1all) + log(p_s2))
  p_ca    <- pchisq(Xall, df = 4, lower.tail = FALSE)
  sel_ca  <- resolve_pick(cand, p_ca < alpha)

  # ---- efficacy estimate at the selected dose (winner's-boost readout) ----
  est_a  <- if (sel_a  > 0) y2[match(sel_a, cand)] / n2 else NA
  est_bc <- if (sel_bc > 0) { i <- match(sel_bc, cand); (y2[i] + y1c[i]) / (n2 + n1c[i]) } else NA
  est_bm <- if (sel_bm > 0) { i <- match(sel_bm, cand); (y2[i] + y1m[i]) / (n2 + n1m[i]) } else NA
  est_cc <- if (!is.na(sel_cc) && sel_cc > 0) y2[match(sel_cc, cand)] / n2 else NA
  est_cm <- if (!is.na(sel_cm) && sel_cm > 0) y2[match(sel_cm, cand)] / n2 else NA
  est_pa <- if (sel_pa > 0) { i <- match(sel_pa, cand); (y2[i] + yAll[i]) / (n2 + nAll[i]) } else NA
  est_ca <- if (sel_ca > 0) y2[match(sel_ca, cand)] / n2 else NA

  door_of <- function(x) if (x > 0) 0 else if (x == -1) 3 else 1
  door_a  <- door_of(sel_a)
  door_bc <- door_of(sel_bc)
  door_bm <- door_of(sel_bm)
  door_cc <- if (is.na(sel_cc)) NA else door_of(sel_cc)
  door_cm <- if (is.na(sel_cm)) NA else door_of(sel_cm)
  door_pa <- door_of(sel_pa)
  door_ca <- door_of(sel_ca)

  list(
    sel  = c(sel_a, sel_bc, sel_bm, sel_cc, sel_cm, sel_pa, sel_ca),
    door = c(door_a, door_bc, door_bm, door_cc, door_cm, door_pa, door_ca),
    est  = c(est_a, est_bc, est_bm, est_cc, est_cm, est_pa, est_ca),
    mtd  = mtd,
    expo = expo, expo_ph = expo_ph,
    dur1 = st1$stage1_dur, dur2 = dur2, dur = st1$stage1_dur + dur2,
    reached_s2 = TRUE
  )
}

# ---------------------------------------------------------------------
# Driver: loop scenarios x lambda x replicates, accumulate, summarize, write
# ---------------------------------------------------------------------
run_naive <- function(cfg) {
  set.seed(seed)
  arms <- c("a_stage2only", "b_partialpool_clean", "b_partialpool_mis",
            "c_combinedp_clean", "c_combinedp_mis",
            "pool_all_exact", "pool_all_combinedp")
  NARM <- length(arms)
  comb_arms <- c(4L, 5L)    # clean/mis combined-p arms: undefined at lambda 0

  method_label <- function(arm, lam) {
    if (arm == "a_stage2only") {
      "stage2only"
    } else if (arm == "b_partialpool_clean") {
      sprintf("partialpool_clean_%.1f", lam)
    } else if (arm == "b_partialpool_mis") {
      sprintf("partialpool_mis_%.1f", lam)
    } else if (arm == "c_combinedp_clean") {
      sprintf("combined_p_clean_%.1f", lam)
    } else if (arm == "c_combinedp_mis") {
      sprintf("combined_p_mis_%.1f", lam)
    } else if (arm == "pool_all_exact") {
      sprintf("pool_all_%.1f", lam)
    } else {
      sprintf("combined_p_all_%.1f", lam)
    }
  }

  sel_rows   <- list(); si <- 0     # per (scn, lambda, arm)
  alloc_rows <- list(); ai <- 0     # per (scn, lambda), shared across arms

  for (s in 1:15) {
    ob <- trueOBD[s]; mt <- trueMTD[s]
    tox_over_doses <- which(S2_tox[s, ] > phiT)      # doses above the MTD (overdosing)

    for (lam in lambda_grid) {

      sel  <- matrix(NA_integer_, NSIM, NARM)
      door <- matrix(NA_integer_, NSIM, NARM)
      est  <- matrix(NA_real_,    NSIM, NARM)
      mtdv <- rep(NA_integer_, NSIM)
      expo <- matrix(0L, NSIM, J)
      expo_ph_vec <- integer(NSIM)    # phantom Stage-2 patients (outside the J doses)
      durv  <- rep(NA_real_, NSIM)
      dur1v <- rep(NA_real_, NSIM)
      dur2v <- rep(NA_real_, NSIM)
      s2v   <- logical(NSIM)

      for (r in 1:NSIM) {
        set.seed(stage1_seed(s, r))   # CRN: identical Stage 1 across families
        o <- one_rep(s, lam)
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
      # Phantom patients are folded into N_used and below_n, so N_used means
      # every patient enrolled in the simulated trial; they are also broken out
      # in phantom_n / phantom_pct. The phantom sits below dose 1, hence below
      # any true OBD, so at_obd, over and aboveMTD are untouched.
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
        Scn = s, trueOBD = ob, trueMTD = mt, lambda = lam,
        regime  = if (lam == 0) "stage2only" else sprintf("pooled_%.1f", lam),
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
        if (a %in% comb_arms && lam == 0) next     # combined-p undefined with no Stage-1 pool
        sa <- sel[, a]; da <- door[, a]; ea <- est[, a]
        selected <- !is.na(sa) & sa > 0

        # sel = -1 (phantom pick) falls out of `selected`, so PCS / UnderOBD /
        # OverOBD / ToxSel / sel_d1..d5 are unaffected and it counts toward
        # noOBD. noOBD therefore exceeds noMTD + noEff by phantomSel (door 3).
        PCS   <- if (is.na(ob)) NA else 100 * mean(sa == ob, na.rm = TRUE)
        Under <- if (is.na(ob)) NA else 100 * mean(selected & sa < ob)
        Over  <- if (is.na(ob)) NA else 100 * mean(selected & sa > ob)
        noOBD <- 100 * mean(!selected)
        noMTD <- 100 * mean(da == 2, na.rm = TRUE)
        noEff <- 100 * mean(da == 1, na.rm = TRUE)
        phantomSel <- 100 * mean(sa == -1, na.rm = TRUE)

        tox_flag <- logical(NSIM)
        tox_flag[selected] <- S2_tox[s, sa[selected]] > phiT
        ToxSel <- 100 * mean(tox_flag)

        any_sel <- 100 * mean(selected)

        if (any(selected)) {
          Est  <- mean(ea[selected])
          True <- mean(S2_eff[s, sa[selected]])
        } else { Est <- NA; True <- NA }

        si <- si + 1
        sel_rows[[si]] <- data.frame(
          Scn = s, family = family[s],
          lambda = lam, arm = arms[a], method = method_label(arms[a], lam),
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
    }
    cat(sprintf("  [naive] scenario %2d done\n", s))
  }

  res   <- do.call(rbind, sel_rows)
  alloc <- do.call(rbind, alloc_rows)

  # ---- output tables (headline methods, at the display lambda) ----
  # Built from lambda_grid with the same %.1f format method_label() uses, so the
  # two cannot drift apart. A lambda with two significant decimals would collide
  # under %.1f in both places.
  LAMBDA_DISP  <- lambda_grid[1]
  disp_methods <- c("stage2only",
                    sprintf("partialpool_clean_%.1f", LAMBDA_DISP),
                    sprintf("partialpool_mis_%.1f",   LAMBDA_DISP),
                    sprintf("pool_all_%.1f",          LAMBDA_DISP),
                    sprintf("combined_p_clean_%.1f",  LAMBDA_DISP),
                    sprintf("combined_p_mis_%.1f",    LAMBDA_DISP),
                    sprintf("combined_p_all_%.1f",    LAMBDA_DISP))
  res$method   <- factor(res$method, levels = c(disp_methods,
                          setdiff(unique(res$method), disp_methods)))
  disp <- res[res$method %in% disp_methods, ]
  disp <- disp[order(disp$Scn, disp$method), ]

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
    Scn = alloc$Scn, trueOBD = alloc$trueOBD, trueMTD = alloc$trueMTD, regime = alloc$regime,
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
  write.csv(res, file.path(od, "setting2_naive_results_full.csv"),     row.names = FALSE)
  write.csv(sel_tab,   file.path(od, "setting2_naive_selection.csv"),  row.names = FALSE)
  write.csv(dist_tab,  file.path(od, "setting2_naive_seldist.csv"),    row.names = FALSE)
  write.csv(alloc_tab, file.path(od, "setting2_naive_allocation.csv"), row.names = FALSE)

  cat("[naive] wrote setting2_naive_{selection,seldist,allocation}.csv and results_full to", od, "\n")

  invisible(list(res = res, sel = sel_tab, dist = dist_tab, alloc = alloc_tab))
}
