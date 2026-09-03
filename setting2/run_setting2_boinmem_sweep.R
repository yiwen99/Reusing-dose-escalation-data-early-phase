# =====================================================================
# Setting 2: BOIN-MEM Stage-2 cohort size sweep (cohorts 8 and 10)
# =====================================================================
# Runs the BOIN-MEM module at Stage-2 cohort size 8 and 10, then stacks those
# results with the existing cohort-4 results. Nothing outside BOIN-MEM is run.
#
# A separate run per cohort size is required because cohort2 is not a variant
# axis. The five VARIANTS differ only in w, gamma and source, all applied to one
# finished trial after the fact, so they stay paired. cohort2 instead sets the
# rmultinom draw size, the length of the t_arr2 arrival loop and the runif count
# for dlt_time2, so the random stream diverges at the first Stage-2 cohort and
# the Stage-2 loop runs a different number of times, which also moves when the
# sint reselection trigger fires.
#
# The module file is not edited. The cohort size is changed from outside by
# assigning into the module's private environment after sys.source: sys.source
# gives every function in the module an enclosing environment of `env`, so a name
# assigned into `env` after sourcing is found before the global environment.
#
# Cohort 4 is not re-run. The driver calls set.seed(stage1_seed(s, r))
# immediately before every one_rep, so Stage 1 is identical for a given (s, r)
# across separate runs and only Stage 2 diverges. This holds only if NSIM
# matches; the metadata CSV written below records what was used.
#
# Run from the setting folder (setting2/), as run_setting2.R assumes.
# =====================================================================

# ---------------------------------------------------------------------
CODE_DIR    <- "."                                # folder holding the config + module
CONFIG_FILE <- "setting2_config.R"
MODULE      <- "s2_method_boinmem_following_package.R"
OUTROOT     <- "results_setting2_boinmem_sweep"   # per-cohort subdirectories
COHORTS     <- c(8L, 10L)     # cohort 4 is not re-run; it is read from BASE_DIR
BASE_DIR    <- "results_setting2"                 # existing cohort-4 results
BASE_COHORT <- 4L

# Leave NA to use NSIM from setting2_config.R. A sweep run at one NSIM cannot be
# compared against cohort-4 results produced at another.
NSIM_OVERRIDE <- NA
# ---------------------------------------------------------------------

source(file.path(CODE_DIR, CONFIG_FILE))   # into globalenv: shared DGP + truth

if (!is.na(NSIM_OVERRIDE)) { NSIM <<- as.integer(NSIM_OVERRIDE); cat("NSIM overridden to", NSIM, "\n") }

# ---- fail loudly and early on paths -------------------------------------
need_files <- c(CONFIG_FILE, MODULE)
miss <- need_files[!file.exists(file.path(CODE_DIR, need_files))]
if (length(miss))
  stop("file(s) not found under CODE_DIR = '", CODE_DIR, "': ",
       paste(miss, collapse = ", "),
       "\n  working directory is: ", getwd())

base_sel_file   <- file.path(BASE_DIR, "setting2_boinmem_selection.csv")
base_alloc_file <- file.path(BASE_DIR, "setting2_boinmem_allocation.csv")
base_miss <- c(base_sel_file, base_alloc_file)[!file.exists(c(base_sel_file, base_alloc_file))]
if (length(base_miss))
  stop("existing cohort-", BASE_COHORT, " result file(s) not found: ",
       paste(base_miss, collapse = ", "),
       "\n  BASE_DIR = '", BASE_DIR, "', working directory is: ", getwd(),
       "\n  the combined sweep files must carry all three cohort sizes;",
       " a two-cohort figure is not what this task produces.")

# ---- output isolation ---------------------------------------------------
# BASE_DIR holds the cohort-4 results this script reuses and is never written to.
if (normalizePath(OUTROOT, mustWork = FALSE) == normalizePath(BASE_DIR, mustWork = FALSE))
  stop("OUTROOT must not point at BASE_DIR ('", BASE_DIR,
       "'), which holds the existing cohort-4 results and must not be overwritten.")
if (BASE_COHORT %in% COHORTS)
  stop("COHORTS must not contain BASE_COHORT = ", BASE_COHORT,
       "; that cohort is read from BASE_DIR, not re-run, and re-running it here",
       " would put duplicate rows into the combined files.")

# ---- read N2 from the module --------------------------------------------
# N2 is defined in the module, not in the config, so it is read from a throwaway
# sourced environment. Every cohort below gets its own fresh environment.
probe  <- new.env(parent = globalenv())
sys.source(file.path(CODE_DIR, MODULE), envir = probe)
N2_mod <- get("N2", envir = probe)
rm(probe)

# ---- divisibility guard -------------------------------------------------
# The Stage-2 loop breaks on a test at the top of the loop and does not truncate
# its final cohort, so a cohort size that does not divide N2 overshoots the
# target sample size and is not comparable.
all_cohorts <- c(BASE_COHORT, COHORTS)
if (any(N2_mod %% all_cohorts != 0))
  stop("cohort size must divide N2 = ", N2_mod,
       "; the Stage-2 loop does not truncate its final cohort, so a ",
       "non-divisor overshoots the target sample size: ",
       paste(all_cohorts[N2_mod %% all_cohorts != 0], collapse = ", "))

if (!dir.exists(OUTROOT)) dir.create(OUTROOT, recursive = TRUE)

# ---- run one configuration per cohort size ------------------------------
sel_parts   <- list()
alloc_parts <- list()
meta_rows   <- list()

for (k in COHORTS) {
  od <- file.path(OUTROOT, sprintf("cohort%02d", k))
  if (!dir.exists(od)) dir.create(od, recursive = TRUE)

  cat("\nSetting 2 | BOIN-MEM + matched (5 variants) | Stage-2 cohort size", k, "\n")

  # Fresh environment per cohort: reusing one would leave the previous
  # assignment in place if a later cohort failed.
  env <- new.env(parent = globalenv())
  sys.source(file.path(CODE_DIR, MODULE), envir = env)

  # Overrides the module's own `cohort2` for this environment only. The module
  # file is untouched, and run_stage1_shared, defined in the config and so
  # enclosing the global environment, is unaffected.
  assign("cohort2", k, envir = env)

  cfg <- list(outdir = od, seed = seed)
  res <- env$run_boinmem(cfg)          # invisible: full, sel, dist, alloc, diagn

  sl <- res$sel;     sl$cohort2 <- k
  al <- res$alloc;   al$cohort2 <- k
  sel_parts[[as.character(k)]]   <- sl
  alloc_parts[[as.character(k)]] <- al
  meta_rows[[as.character(k)]] <- data.frame(
    cohort2 = k, NSIM = NSIM, outdir = od,
    note = "run by run_setting2_boinmem_sweep.R",
    stringsAsFactors = FALSE)
}

# ---- stack in the existing cohort-4 results -----------------------------
base_sel <- read.csv(base_sel_file,   stringsAsFactors = FALSE)
base_al  <- read.csv(base_alloc_file, stringsAsFactors = FALSE)
base_sel$cohort2 <- BASE_COHORT
base_al$cohort2  <- BASE_COHORT

# Schema drift between the stored cohort-4 files and what the module writes now
# would silently misalign the stacked frame, so refuse rather than rbind blindly.
new_sel_names <- names(sel_parts[[1]])
new_al_names  <- names(alloc_parts[[1]])
if (!identical(sort(names(base_sel)), sort(new_sel_names)))
  stop("column mismatch between ", base_sel_file, " and the new selection tables.",
       "\n  stored: ", paste(sort(names(base_sel)), collapse = ", "),
       "\n  new   : ", paste(sort(new_sel_names), collapse = ", "),
       "\n  the stored cohort-4 file was written by a different version of the module;",
       " re-export it or drop it from the sweep.")
if (!identical(sort(names(base_al)), sort(new_al_names)))
  stop("column mismatch between ", base_alloc_file, " and the new allocation tables.",
       "\n  stored: ", paste(sort(names(base_al)), collapse = ", "),
       "\n  new   : ", paste(sort(new_al_names), collapse = ", "))

sel_all   <- do.call(rbind, c(list(base_sel[, new_sel_names]), sel_parts))
alloc_all <- do.call(rbind, c(list(base_al[, new_al_names]),  alloc_parts))

# The cohort-4 metadata row. If no meta file sits beside the stored results,
# NSIM is NA here and the plotting script refuses to draw the cohort-4 point
# until the value is confirmed by hand.
base_meta_file <- file.path(BASE_DIR, "setting2_boinmem_meta.csv")
base_nsim <- NA_integer_
base_note <- paste("pre-dates the metadata convention: NSIM not recorded.",
                   "Confirm it by hand in the plotting script.")
if (file.exists(base_meta_file)) {
  bm <- read.csv(base_meta_file, stringsAsFactors = FALSE)
  base_nsim <- as.integer(bm$NSIM[1])
  base_note <- paste("NSIM read from", base_meta_file)
}
meta_rows[["base"]] <- data.frame(
  cohort2 = BASE_COHORT, NSIM = base_nsim, outdir = BASE_DIR,
  note = base_note, stringsAsFactors = FALSE)

meta_all <- do.call(rbind, meta_rows[c("base", as.character(COHORTS))])
meta_all <- meta_all[order(meta_all$cohort2), ]

write.csv(sel_all,   file.path(OUTROOT, "setting2_boinmem_cohort_sweep_selection.csv"),  row.names = FALSE)
write.csv(alloc_all, file.path(OUTROOT, "setting2_boinmem_cohort_sweep_allocation.csv"), row.names = FALSE)
write.csv(meta_all,  file.path(OUTROOT, "setting2_boinmem_cohort_sweep_meta.csv"),       row.names = FALSE)

cat(sprintf("\nBOIN-MEM cohort sweep complete: cohorts %s re-run, cohort %d reused from %s  ->  %s\n",
            paste(COHORTS, collapse = ", "), BASE_COHORT, BASE_DIR, OUTROOT))

if (is.na(base_nsim))
  cat("\nNOTE: the cohort-4 metadata row carries NSIM = NA. Confirm the NSIM those\n",
      "     results used and set BASE_NSIM_CONFIRMED in plots/boinmem_cohort_sweep_plot.R\n",
      "     before drawing the figure.\n", sep = "")
