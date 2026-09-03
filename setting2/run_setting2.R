# =====================================================================
# Setting 2: runner
# =====================================================================
# Sources the shared config into the global environment, then sources each
# requested method module into its own private environment and runs it through
# its run_<class>(cfg) entry point.
#
# The private environments matter because every method module defines top-level
# helpers under the same names (run_stage1, one_rep, NARM, draw_gumbel, ...).
# Sourcing them into one environment would let the last-sourced module overwrite
# the earlier ones. Each environment has the global environment as parent, so the
# shared config stays visible, and each run_<class> closes over its own
# environment.
#
# Each method re-seeds from the shared `seed`, so run order and the RUN subset
# never perturb another method's numbers.
#
# To add a method: drop a module defining run_<newclass>(cfg) into this folder,
# add one row to `specs`, and list its name in RUN.
# =====================================================================

# ---------------------------------------------------------------------
# Run from the project root (the folder holding setting1/, setting2/, plots/).
# ---------------------------------------------------------------------
CODE_DIR    <- "."                 # folder holding the config + modules
CONFIG_FILE <- "setting2_config.R"
OUTDIR      <- "results_setting2"        # CSV output directory (created if missing)
RUN         <- c("naive", "boinmem", "bard")

# Leave NA to use NSIM from setting2_config.R.
NSIM_OVERRIDE <- NA
# ---------------------------------------------------------------------

source(file.path(CODE_DIR, CONFIG_FILE))   # into globalenv: shared DGP + truth

if (!is.na(NSIM_OVERRIDE)) { NSIM <<- as.integer(NSIM_OVERRIDE); cat("NSIM overridden to", NSIM, "\n") }

if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)
cfg <- list(outdir = OUTDIR, seed = seed)

# method registry: class -> module file, run-function name, label
specs <- list(
  naive   = list(file = "s2_method_naive.R",   fn = "run_naive",   label = "naive / partial-pooling (7 arms)"),
  boinmem = list(file = "s2_method_boinmem.R", fn = "run_boinmem", label = "BOIN-MEM + matched (5 variants)"),
  bard    = list(file = "s2_method_bard.R",    fn = "run_bard",
                 label = "BARD (6 arms: match x CAR) + BF-BOIN-SR (3 arms: match)")
)

# Fail loudly and early if the config or a module is not where the block says.
need_files <- c(CONFIG_FILE, vapply(specs, function(x) x$file, character(1)))
miss <- need_files[!file.exists(file.path(CODE_DIR, need_files))]
if (length(miss))
  stop("module(s) not found under CODE_DIR = '", CODE_DIR, "': ",
       paste(miss, collapse = ", "),
       "\n  working directory is: ", getwd())

unknown <- setdiff(RUN, names(specs))
if (length(unknown)) stop("unknown method(s) in RUN: ", paste(unknown, collapse = ", "))

for (m in RUN) {
  cat("\nSetting 2 |", specs[[m]]$label, "\n")
  env <- new.env(parent = globalenv())                       # private per-module scope
  sys.source(file.path(CODE_DIR, specs[[m]]$file), envir = env)
  env[[ specs[[m]]$fn ]](cfg)
}
cat(sprintf("\nSetting 2 complete: %s  ->  %s\n", paste(RUN, collapse = ", "), OUTDIR))
