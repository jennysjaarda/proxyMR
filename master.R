#master script

library(workflowr)
wflow_start("proxyMR")

renv::init()
tar_renv()

install.packages("remotes")
remotes::install_github("ropensci/targets")
remotes::install_github("ropensci/tarchetypes")

library(targets)

#manually create R/ folder
#create 'R/function.R' file
#add '_drake.R' file

wflow_use_github("jennysjaarda")

## On terminal

#git clone https://github.com/jennysjaarda/proxyMR.git

# only do this once

Neale_SGG_dir_org <- paste0(Neale_summary_dir,"/Neale_SGG_directory.csv") ## not tracking this as input because it changes overtime
Neale_SGG_dir_cp <- paste0("data/Neale_SGG_directory_", format(Sys.Date(), "%d_%m_%Y"), ".csv")
file.copy(Neale_SGG_dir_org, Neale_SGG_dir_cp)

dir.create("output/figures")
dir.create("output/tables")


stats1 <- 2
stats2 <- 2
while ( (length(stats1) > 1 | length(stats2) > 1)){
  stats1 <- suppressWarnings(system(paste("squeue -n", "process_Neale"), intern = TRUE))
  stats2 <- suppressWarnings(system(paste("squeue -n", "clump_Neale_IVs"), intern = TRUE))
  print("Still running...")
  Sys.sleep(60)
}
