## code to prepare `DATASET` dataset goes here

#import data about columns, that is each column should be metadata about the sample/animal. Particularly include all treatment info and animal ID in columns.
coldata <-  data.table::fread(file = "data-raw/coldata_new_names_20201105.csv")

coldata <- dplyr::mutate(coldata,
                         sample_names = paste0(Region, "_", ID, "__", Diet))

#test to throw error
#coldata <-  data.table::fread(file = "./coldata_test.csv") %>%
# dplyr::mutate(sample_names = paste0(Region, "_", ID, "__", Diet))


cts_all <- data.table::fread(file= "data-raw/cts_renamed_all_Sheep.txt")



usethis::use_data(coldata, cts_all, overwrite = TRUE)
