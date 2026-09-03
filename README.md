# Simulation code

Code for the two-stage OBD-selection simulation study. Setting 1 and Setting 2
are separate and self-contained, each in its own folder.

## Layout

```
setting1/
  setting1_config.R              shared setup and data-generating process
  s1_method_naive.R              naive pooling / combined-p
  s1_method_boinmem.R            BOIN-MEM
  s1_method_bard.R               BARD + BF-BOIN-SR
  run_setting1.R                 runs all three methods

setting2/
  setting2_config.R              shared setup and data-generating process
  s2_method_naive.R              naive pooling / combined-p
  s2_method_boinmem.R            BOIN-MEM (+ matched variants)
  s2_method_bard.R               BARD + BF-BOIN-SR (+ matched variants)
  run_setting2.R                 runs all three methods
  run_setting2_boinmem_sweep.R   Stage-2 cohort size lever
  run_setting2_boinmem_fastread.R   assessment window lever
```

**The config file** has the settings the three methods share: the scenario rate
tables, the true MTD and OBD, the accrual calendar, and the Stage-1 escalation.
It is also where the main knobs live, including `NSIM`, `phiT`, `phiE` and `N1`.
Change something here and all three methods use the new value.

**The method files** are one script per method, each running all of that
method's variants. Settings that belong to a single method, such as its Stage-2
sample size or its arm list, sit at the top of its own file rather than in the
config.

**The wrapper** `run_settingX.R` loads the config, runs each method, and saves
the results.

## Running it

You need R with the `BOIN` and `Minirand` packages.

Run from inside the setting folder:

```r
setwd("setting1")
source("run_setting1.R")
```

Same for Setting 2. Results are saved to `results_setting1/` or
`results_setting2/`, next to the code.

To try a quick version first, set `NSIM_OVERRIDE` near the top of the wrapper to
something small like 50. The full study uses `NSIM = 5000`.

## Output

Each method saves these files to the results folder:

| file | contents |
|---|---|
| `..._selection.csv` | PCS, under/over selection, no-OBD breakdown |
| `..._seldist.csv` | selection percentage by dose |
| `..._allocation.csv` | patient allocation, trial duration |
| `..._results_full.csv` | everything, unrounded |

Setting 1 BOIN-MEM, for example, gives `setting1_boinmem_selection.csv` and so
on. BOIN-MEM also saves `..._diagnostics.csv` and BARD `..._metrics.csv`.

## The two duration levers

The manuscript looks at two ways of shortening a BOIN-MEM trial, and each has
its own script. `run_setting2_boinmem_sweep.R` varies the Stage-2 cohort size;
`run_setting2_boinmem_fastread.R` shortens the toxicity and efficacy
assessment windows. Both run BOIN-MEM only and save to their own folders. Run
them the same way:

```r
setwd("setting2")
source("run_setting2_boinmem_sweep.R")
```

The cohort sweep reuses the existing cohort-4 results from `results_setting2/`,
so run `run_setting2.R` first.
