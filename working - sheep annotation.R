#Biomart

biomaRt::listEnsembl()
ensembl <- biomaRt::useEnsembl(biomart = "genes")
biomaRt::listDatasets(ensembl) #this shows the list of species available

ensembl <- biomaRt::useDataset("oaries_gene_ensembl", mart = ensembl, verbose = TRUE)

biomaRt::listAttributes(ensembl)
biomaRt::listFilters(ensembl)
#clear cache
biomaRt::biomartCacheClear()


test <-
  biomaRt::getBM(attributes = c('external_gene_name'),
      filters = 'ensembl_gene_id',
      values = test_annot$gene_ensembl,
      mart = ensembl,
      verbose = TRUE)


test
ens <- test_annot$gene_name %>% stringr::str_detect(pattern = "ENSOARG")

non_ens <- test_annot$gene_name[!ens]

length(non_ens %>% unique())



unofficial_names <-
  data.table::fread(file="data/Clarke_et_al_gene_annotations.csv")

colnames(unofficial_names)[1] <- "gene_ensembl"


#reference!!
link_to_ref <- unofficial_names %>%
  dplyr::select(reference) %>% .[1]



test_annot_2 <- test_annot %>%
  dplyr::left_join(unofficial_names)

unofficial_names_2 <-
  data.table::fread(file="data/Clarke_et_al_longer_descriptions.csv")

colnames(unofficial_names_2)[1] <- "gene_ensembl"

test_annot_3 <-
  test_annot_2 %>%
  dplyr::left_join(unofficial_names_2)
