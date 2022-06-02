## FTICR-MS DRAKE PLAN
## USE THIS SCRIPT/PLAN TO PROCESS, ANALYZE, AND GRAPH FTICR-MS DATA
## KAIZAD F. PATEL
## OCT-06-2020

##############################
##############################

## SOURCE THE FUNCTION FILES FIRST, AND THEN RUN THE DRAKE PLAN
## DON'T RUN THE PROCESSING PLAN MULTIPLE TIMES. ONCE IS ENOUGH.

##############################
##############################


# 0. load packages --------------------------------------------------------
library(drake)
library(tidyverse)
library(PNWColors)
library(soilpalettes)
library(readxl)

# 1. SET input file paths -------------------------------
#COREKEY = "Ficus_data/Ficus_sample_key.txt"
REPORT = "1-data/fticr/xtra_51537_FICUS_Met_21T_Report.xlsx"
REPORT_12T = "1-data/fticr/xtra_51537_FICUS_Met_12T_Report.xlsx"

COREKEY = "1-data/Ficus_sample_key.csv"
REPORT_LIPID = "1-data/fticr/xtra_51537_FICUS_Lip_12T_Report.xlsx"

## SET the treatment variables
TREATMENTS = quos(Incubation_Temp, A_priori_CUE)

# 2. source the functions --------------------------------------------------------
source("2-code/fticrrr/a-functions_processing.R")
source("2-code/fticrrr/b-functions_relabund.R")
source("2-code/fticrrr/c-functions_vankrevelen.R")
source("2-code/fticrrr/d-functions_statistics.R")

# 3. load drake plans -----------------------------------------------------
fticr_processing_plan = drake_plan(
  # 1. METABOLITES 21T ----
  ## a. PROCESSING 
  datareport = read_excel(file_in(REPORT)),

  corekey = read.csv(file_in(COREKEY)),

  fticr_meta = make_fticr_meta(datareport)$meta2,
  fticr_data_longform = make_fticr_data(datareport, 
                                        corekey, 
                                        TREATMENTS)$data_long_key_repfiltered,
  fticr_data_trt = make_fticr_data(datareport, corekey, TREATMENTS)$data_long_trt,
  
  ## b. RELATIVE ABUNDANCE  
  relabund_cores = fticr_data_longform %>% 
    compute_relabund_cores(fticr_meta, TREATMENTS),
  
  gg_relabund_bar = relabund_cores %>% plot_relabund(TREATMENTS),
  
  ## create relabund table
  
  ## c. VAN KREVELEN PLOTS
  gg_vankrevelen_domains = plot_vankrevelen_domains(fticr_meta),
  gg_vankrevelens = plot_vankrevelens(fticr_data_trt, fticr_meta),
  gg_vankrevelen_unique = plot_vk_unique(fticr_data_trt, fticr_meta),

  ## d. STATISTICS
  ## PERMANOVA
  fticr_permanova = compute_permanova(relabund_cores),

  ## PCA
  gg_pca = compute_fticr_pca(relabund_cores), 

  # e. OUTPUT FILES
  fticr_meta %>% write.csv("1-data/processed/fticr/fticr_meta_21T.csv", row.names = FALSE),
  fticr_data_trt %>% write.csv("1-data/processed/fticr/fticr_data_by_treatment_21T.csv", row.names = FALSE),
  fticr_data_longform %>% write.csv("1-data/processed/fticr/fticr_data_longform_21T.csv", row.names = FALSE), 
  
  
  #
  # 2. LIPIDOMICS 12T ----
  
  ## a. PROCESSING
  datareport_lipid = read_excel(file_in(REPORT_LIPID)),
  
  lipid_meta = make_fticr_meta(datareport_lipid)$meta2,
  lipid_data_longform = make_lipid_data(datareport_lipid, 
                                        corekey, 
                                        TREATMENTS)$data_long_key_repfiltered,
  lipid_data_trt = make_lipid_data(datareport_lipid, corekey, TREATMENTS)$data_long_trt,
  
  ## b. RELATIVE ABUNDANCE
  relabund_cores_lipid = lipid_data_longform %>% 
    compute_relabund_cores(lipid_meta, TREATMENTS),
  
  gg_relabund_bar_lipid = relabund_cores_lipid %>% plot_relabund(TREATMENTS),
  
  ## create relabund table
  
  ## c. VAN KREVELEN PLOTS
  gg_vankrevelen_domains_lipid = plot_vankrevelen_domains(lipid_meta),
  gg_vankrevelens_lipid = plot_vankrevelens(lipid_data_trt, fticr_meta),
  gg_vankrevelen_unique_lipid = plot_vk_unique(lipid_data_trt, fticr_meta),
  
  ## d. STATISTICS
  ## PERMANOVA
  lipid_permanova = compute_permanova(relabund_cores_lipid),
  
  ## PCA
  #gg_pca = compute_fticr_pca(relabund_cores_lipid), 
  
  ## e. OUTPUT FILES
  lipid_meta %>% write.csv("1-data/processed/fticr/lipid_meta_12T.csv", row.names = FALSE),
  lipid_data_trt %>% write.csv("1-data/processed/fticr/lipid_data_by_treatment_12T.csv", row.names = FALSE),
  lipid_data_longform %>% write.csv("1-data/processed/fticr/lipid_data_longform_12T.csv", row.names = FALSE), 
  
  #
  # 2. METABOLITES 12T ----
  
  ## a. PROCESSING
  datareport_12T = read_excel(file_in(REPORT_12T)),
  
  meta_12T = make_fticr_meta(datareport_12T)$meta2,
  data_longform_12T = make_fticr_data_12T(datareport_12T, 
                                        corekey, 
                                        TREATMENTS)$data_long_key_repfiltered,
  data_trt_12T = make_fticr_data_12T(datareport_12T, corekey, TREATMENTS)$data_long_trt,
  
  ## b. RELATIVE ABUNDANCE
  relabund_cores_12T = data_longform_12T %>% 
    compute_relabund_cores(meta_12T, TREATMENTS),
  
  gg_relabund_bar_12T = relabund_cores_12T %>% plot_relabund(TREATMENTS),
  
  ## create relabund table
  
  ## c. VAN KREVELEN PLOTS
  gg_vankrevelen_domains_12T = plot_vankrevelen_domains(meta_12T),
  gg_vankrevelens_12T = plot_vankrevelens(data_trt_12T, meta_12T),
  gg_vankrevelen_unique_12T = plot_vk_unique(data_trt_12T, meta_12T),
  
  ## d. STATISTICS
  ## PERMANOVA
  permanova_12T = compute_permanova(relabund_cores_12T),
  
  ## PCA
  gg_pca_12T = compute_fticr_pca(relabund_cores_12T), 
  
  ## e. OUTPUT FILES
  meta_12T %>% write.csv("1-data/processed/fticr/metab_meta_12T.csv", row.names = FALSE),
  data_trt_12T %>% write.csv("1-data/processed/fticr/metab_data_by_treatment_12T.csv", row.names = FALSE),
  data_longform_12T %>% write.csv("1-data/processed/fticr/metab_data_longform_12T.csv", row.names = FALSE), 
  
  #
  # REPORT ----
  outputreport = rmarkdown::render(
    knitr_in("reports/report_fticr.Rmd"),
    output_format = rmarkdown::github_document())
)


# 4. make plans -------------------------------------------------------------------------
make(fticr_processing_plan, lock_cache = F)


