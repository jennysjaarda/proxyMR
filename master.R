#master script

library(workflowr)
wflow_start("proxyMR")

renv::init()

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
Neale_SGG_dir_cp <- file.copy(Neale_SGG_dir_org, Neale_SGG_dir_cp)
