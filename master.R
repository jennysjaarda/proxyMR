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


