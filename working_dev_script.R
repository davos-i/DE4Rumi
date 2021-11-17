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
usethis::use_package("DESeq2")
usethis::use_package("rlang")
usethis::use_package("IHW")
usethis::use_package("S4Vectors", type = "Depends") #IHW relies on names that aren't defined properly in their package
usethis::use_package("SummarizedExperiment")
usethis::use_package("stringr")

usethis::use_package("BiocGenerics")
usethis::use_package("data.table")
usethis::use_package("genefilter")
usethis::use_package("ggplot2")
usethis::use_package("ggrepel")
usethis::use_package("pheatmap")
usethis::use_package("purrr")
usethis::use_package("tidyr")
usethis::use_package("tidyselect")
usethis::use_package("graphics")
usethis::use_package("grDevices")

#uses .data from rlang
usethis::use_import_from("rlang", ".data")
#usethis::use_import_from("tidyselect", "where") #fail
#have to put utils::globalVariables("where") inside one of the .R files
usethis::use_import_from("S4Vectors", "mcols")
usethis::use_import_from("S4Vectors", "mcols<-")
usethis::use_import_from("S4Vectors", "metadata")
usethis::use_import_from("S4Vectors", "metadata<-")

#usethis::use_import_from("rlang", "ensym")


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
#' @return returns a value of...
#'
#' @export


########################################################################### #






#notes
two <- coldata %>% select(rowname) %>% magrittr::extract2("rowname")
