source("est.R")
cal <- calibrate_complete()                          # full run, prints progress
cal <- calibrate_complete(save_dir = "results/tables")  # also writes calibration_limits.csv and calibration.rds
cal$L; cal$alpha; cal$UCL_T2; cal$check