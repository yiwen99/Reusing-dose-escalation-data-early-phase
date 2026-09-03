# =====================================================================
# Setting 2: BOIN-MEM fast-read outcome variants
# =====================================================================
# Five timing arms, otherwise identical:
#
#   standard                    cohort 4;  no shortened window
#   fast efficacy read          cohort 4;  eff.window shortened
#   fast toxicity read          cohort 4;  dlt.window shortened
#   fast toxicity read          cohort 10; dlt.window shortened
#   fast toxicity + efficacy    cohort 10; both windows shortened
#
# Only the named window(s) change; every other element of the duration model is
# as in the standard run, and all fast arms use the same FAST_WINDOW.
#
# The toxicity gate paces Stage 2 between cohorts, so it has the obvious claim on
# trial duration, but the trial ends only once the last patient's efficacy has
# also matured, and eff.window is 3 months against dlt.window's 1. The fifth arm
# shortens both, so its difference from the fourth is exactly the residual
# efficacy tail at the end of Stage 2.
#
# The shortened windows apply to Stage 2 only. Stage 1 keeps the original
# assessment windows, so the question asked is what a faster Stage-2 readout
# would buy given a conventionally run Stage 1 shared with the other method
# families. This is also the only sound implementation, because Stage 1's
# inter-cohort stagger is
#
#   dlt_time <- ifelse(oc$tox == 1L, runif(cohort1, 0, dlt.window), dlt.window)
#   matur    <- max(t_arr + dlt_time)
#   repeat { clock <- clock + arrival_gap(); if (clock > matur) break }
#
# and that repeat loop consumes a variable number of random draws determined by
# matur, which depends on dlt.window. Shortening the window globally would make
# the loop run a different number of times, shifting the whole random stream and
# changing the Stage-1 dose-escalation path. The windows are therefore never
# assigned in the global environment or edited in setting2_config.R; they are
# assigned into the module's private environment after sys.source.
#
# Consequently the selection results must not move across the three cohort-4
# arms, nor between the two cohort-10 arms. In Stage 2 neither window feeds any
# decision rule: tox_admissible, eff_admissible, select_R, simple_dominance and
# the sint trigger all work on patient counts, and the windows enter only the
# timing arithmetic. runif(cohort2, 0, dlt.window) consumes cohort2 uniforms
# whatever the upper bound, and the arrival loop is a fixed seq_len(cohort2), so
# Stage 2 is RNG-neutral to both windows. Only the duration columns may differ
# within a cohort group.
#
# The cohort size is the one thing that is not RNG-neutral: it sets the rmultinom
# draw size, the length of the Stage-2 arrival loop and the runif count for
# dlt_time2, so the Stage-2 stream diverges and selection results legitimately
# move between the cohort-4 and cohort-10 groups. The matching standard cohort-10
# reference is not re-run here; it exists in
# results_setting2_boinmem_sweep/cohort10 at the same NSIM and seed.
#
# The standard arm is run rather than reused from results_setting2/, so all five
# bars come from one invocation at one NSIM. It is expected to reproduce the
# existing cohort-4 results exactly.
#
# Run from the setting folder (setting2/), as run_setting2.R assumes.
# =====================================================================

# ---------------------------------------------------------------------
CODE_DIR    <- "."                                      # folder holding the config + module
CONFIG_FILE <- "setting2_config.R"
MODULE      <- "s2_method_boinmem_following_package.R"
OUTROOT     <- "results_setting2_boinmem_fastread_v2"   # per-arm subdirectories

FAST_WINDOW <- 0.25    # months, every fast-read arm. The standard windows are
                       # dlt.window = 1 and eff.window = 3 months. A value of
                       # zero would model an instantaneous assay, making the fast
                       # arms a different design rather than a faster one.

NSIM_OVERRIDE <- NA    # leave NA to use NSIM from setting2_config.R
# ---------------------------------------------------------------------

source(file.path(CODE_DIR, CONFIG_FILE))   # into globalenv: shared DGP + truth

if (!is.na(NSIM_OVERRIDE)) { NSIM <<- as.integer(NSIM_OVERRIDE); cat("NSIM overridden to", NSIM, "\n") }

need_files <- c(CONFIG_FILE, MODULE)
miss <- need_files[!file.exists(file.path(CODE_DIR, need_files))]
if (length(miss))
  stop("file(s) not found under CODE_DIR = '", CODE_DIR, "': ",
       paste(miss, collapse = ", "),
       "\n  working directory is: ", getwd())

if (is.na(FAST_WINDOW) || FAST_WINDOW <= 0)
  stop("FAST_WINDOW must be a positive number of months; a zero or negative ",
       "window models an instantaneous assay, which is a different design.")
if (FAST_WINDOW >= min(dlt.window, eff.window))
  stop("FAST_WINDOW = ", FAST_WINDOW, " is not shorter than the standard windows ",
       "(dlt.window = ", dlt.window, ", eff.window = ", eff.window, ").")

# Output isolation: the protected trees below hold results that this runner
# reuses or that the plotting scripts read, and are never written to here.
protected <- c("results_setting2", "results_setting2_boinmem_fastread",
               "results_setting2_boinmem_sweep")
if (normalizePath(OUTROOT, mustWork = FALSE) %in% normalizePath(protected, mustWork = FALSE))
  stop("OUTROOT must not point at an existing result tree (", 
       paste(protected, collapse = ", "), "); those hold results that must not ",
       "be overwritten.")

if (!dir.exists(OUTROOT)) dir.create(OUTROOT, recursive = TRUE)

# ---- the five timing arms, as a registry ----
# NA leaves the module reading the module or global value for that field.
TIMING <- list(
  standard      = list(dir = "standard",      dlt = NA,          eff = NA,          cohort = NA,
                       label = "standard, cohort 4"),
  fast_eff      = list(dir = "fast_eff",      dlt = NA,          eff = FAST_WINDOW, cohort = NA,
                       label = "fast efficacy read, cohort 4"),
  fast_tox      = list(dir = "fast_tox",      dlt = FAST_WINDOW, eff = NA,          cohort = NA,
                       label = "fast toxicity read, cohort 4"),
  fast_tox_c10  = list(dir = "fast_tox_c10",  dlt = FAST_WINDOW, eff = NA,          cohort = 10L,
                       label = "fast toxicity read, cohort 10"),
  fast_both_c10 = list(dir = "fast_both_c10", dlt = FAST_WINDOW, eff = FAST_WINDOW, cohort = 10L,
                       label = "fast toxicity + efficacy read, cohort 10"))

# ---- read N2 and the module's default cohort2 from the module itself. The
# probe environment is discarded; every arm gets a fresh one.
probe   <- new.env(parent = globalenv())
sys.source(file.path(CODE_DIR, MODULE), envir = probe)
N2_mod  <- as.integer(get("N2",      envir = probe))
C2_mod  <- as.integer(get("cohort2", envir = probe))   # the module writes it as a double
rm(probe)

# Divisibility guard, as in run_setting2_boinmem_sweep.R: the Stage-2 loop breaks
# on a test at the top of the loop and does not truncate its final cohort, so a
# cohort size that does not divide N2 overshoots the target sample size.
req_cohorts <- vapply(TIMING, function(x) if (is.na(x$cohort)) C2_mod else as.integer(x$cohort), integer(1))
bad <- unique(req_cohorts[N2_mod %% req_cohorts != 0])
if (length(bad))
  stop("Stage-2 cohort size must divide N2 = ", N2_mod,
       "; the Stage-2 loop does not truncate its final cohort, so a non-divisor ",
       "overshoots the target sample size: ", paste(bad, collapse = ", "))

sel_parts   <- list()
alloc_parts <- list()

for (tn in names(TIMING)) {
  tm <- TIMING[[tn]]
  od <- file.path(OUTROOT, tm$dir)
  if (!dir.exists(od)) dir.create(od, recursive = TRUE)

  cat("\nSetting 2 | BOIN-MEM + matched (5 variants) | timing arm:", tm$label, "\n")

  # Fresh environment per arm: a reused one would keep a previous window
  # assignment alive if a later arm failed.
  env <- new.env(parent = globalenv())
  sys.source(file.path(CODE_DIR, MODULE), envir = env)

  # Assigning into `env` reconfigures the module only. sys.source gives every
  # function defined in the module an enclosing environment of `env`, so the
  # Stage-2 timing arithmetic in one_rep finds the shortened window before the
  # global one. run_stage1_shared is defined in setting2_config.R, sourced into
  # the global environment, so it keeps reading the unmodified global windows and
  # Stage 1 stays identical across arms. Nothing is ever assigned to dlt.window
  # or eff.window outside `envir = env`.
  if (!is.na(tm$dlt)) assign("dlt.window", tm$dlt, envir = env)   # Stage 2 only
  if (!is.na(tm$eff)) assign("eff.window", tm$eff, envir = env)   # Stage 2 only

  # Same mechanism: overrides the module's own cohort2 for this arm only.
  # run_stage1_shared uses cohort1 and is unaffected either way.
  if (!is.na(tm$cohort)) assign("cohort2", as.integer(tm$cohort), envir = env)
  k <- get("cohort2", envir = env)

  cat(sprintf("  Stage-2 windows: dlt.window = %.3f, eff.window = %.3f  (global: %.3f, %.3f)\n",
              get("dlt.window", envir = env), get("eff.window", envir = env),
              dlt.window, eff.window))
  cat(sprintf("  Stage-2 cohort size: %d  (module default: %d)\n", k, C2_mod))

  cfg <- list(outdir = od, seed = seed)
  res <- env$run_boinmem(cfg)          # invisible: full, sel, dist, alloc, diagn

  sl <- res$sel;   sl$timing <- tn; sl$cohort2 <- k
  al <- res$alloc; al$timing <- tn; al$cohort2 <- k
  sel_parts[[tn]]   <- sl
  alloc_parts[[tn]] <- al
}

sel_all   <- do.call(rbind, sel_parts)
alloc_all <- do.call(rbind, alloc_parts)

# The selection file exists so the equality check can be run on it, not because
# it will be plotted.
write.csv(sel_all,   file.path(OUTROOT, "setting2_boinmem_fastread_selection.csv"),  row.names = FALSE)
write.csv(alloc_all, file.path(OUTROOT, "setting2_boinmem_fastread_allocation.csv"), row.names = FALSE)

meta <- data.frame(
  timing = names(TIMING),
  label  = vapply(TIMING, function(x) x$label, character(1)),
  cohort2 = req_cohorts,
  dlt_window_stage2 = vapply(TIMING, function(x) if (is.na(x$dlt)) dlt.window else x$dlt, numeric(1)),
  eff_window_stage2 = vapply(TIMING, function(x) if (is.na(x$eff)) eff.window else x$eff, numeric(1)),
  dlt_window_stage1 = dlt.window,
  eff_window_stage1 = eff.window,
  FAST_WINDOW = FAST_WINDOW,
  NSIM = NSIM,
  outdir = file.path(OUTROOT, vapply(TIMING, function(x) x$dir, character(1))),
  stringsAsFactors = FALSE)
write.csv(meta, file.path(OUTROOT, "setting2_boinmem_fastread_meta.csv"), row.names = FALSE)

cat(sprintf("\nBOIN-MEM fast-read study complete: %s  ->  %s\n",
            paste(names(TIMING), collapse = ", "), OUTROOT))

# Mean duration per arm. The difference between the two cohort-10 arms is the
# residual Stage-2 efficacy tail.
d <- tapply(alloc_all$duration, alloc_all$timing, mean)
cat("\nmean trial duration (months), by arm:\n")
for (tn in names(TIMING)) cat(sprintf("  %-14s %6.2f   (%s)\n", tn, d[[tn]], TIMING[[tn]]$label))
cat(sprintf("\n  fast_both_c10 minus fast_tox_c10 = %.2f months\n",
            d[["fast_both_c10"]] - d[["fast_tox_c10"]]))

# =====================================================================
# Equality check, to be run by hand after the run.
#
# Within a cohort group cohort2 is fixed, so the Stage-2 random stream is fixed
# and the windows are RNG-neutral: the three cohort-4 arms must agree on every
# non-duration column, and so must the two cohort-10 arms. Across the groups
# nothing is claimed, because cohort2 moves the stream. The standard arm must
# also reproduce the stored cohort-4 results in results_setting2/ column for
# column, duration included. Compare with all.equal at a small tolerance rather
# than identical, since the tables go through a write.csv / read.csv round trip.
# =====================================================================
