#Working scripts for package development.

#check status of projects and directories
usethis::proj_sitrep()


#test package as is
devtools::load_all() #Ctrl + Shift + L


#This has been added to .Rbuildignore using usethis::use_build_ignore()
usethis::use_build_ignore("./working_dev_script.R")


#the usethis package is very useful for most required package dev
#library(usethis)


#See instructions for creating documentation with roxegon: https://r-pkgs.org/man.html
#uses command: devtools::document()


#Require a package in this package
usethis::use_package("dplyr")
usethis::use_package("tibble")
usethis::use_package("crayon")
usethis::use_package("gprofiler2")
usethis::use_package("pbapply")

#uses .data from rlang
usethis::use_import_from("rlang", ".data")


#It's handy to have raw data that is compiled to a .RData format for formal use



#added use of pipe internally with
usethis::use_pipe()




########################################################################## #
# Documenting functions
# Template

#' Title
#'
#' \code{words_to_look_like_code} short description
#'
#' long description
#'
#' @param data dataframe.
#'
#' @export


########################################################################### #






#notes
two <- coldata %>% select(rowname) %>% magrittr::extract2("rowname")
