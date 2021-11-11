#Working scripts for package development.

#check status of projects and directories
usethis::proj_sitrep()


#test package as is
devtools::load_all()


#This has been added to .Rbuildignore using usethis::use_build_ignore()
usethis::use_build_ignore("./working_dev_script.R")


#the usethis package is very useful for most required package dev
#library(usethis)


#See instructions for creating documentation with roxegon: https://r-pkgs.org/man.html
#uses command: devtools::document()


#Require a package in this package
usethis::use_package("dplyr")
usethis::use_package("tibble")


#It's handy to have raw data that is compiled to a .RData format for formal use
