# FTICRRR: fticr results in R
# Kaizad F. Patel
# October 2020

################################################## #

# `functions_processing.R`

################################################## #

## this script will load functions for:
## (a) processing FTICR reports generated in Formularity
## -- (a1) filtering peaks 
## -- (a2) computing indices, element composition, class assignment for metadata 
## -- (a3) cleaning the data file and creating a longform version 
## -- (a4) converting peak intensities into presence/absence
## (b) computing relative abundance by compound class, for each core

# INSTRUCTIONS:
## source this file in the `fticr_drake_plan.R` file, do not run the script here.
## This script can (generally) be used as is for most data that follow this format. No modifications needed 

## 20-Nov-2020 edit: `data_presence2` file is unique to this experiment because of the shitty sample nomenclature.


################################################## #
################################################## #

# 1. PROCESSING FUNCTIONS -------------------------------------------------
## LEVEL I FUNCTIONS -------------------------------------------------------
## for metadata file
apply_filter_report = function(report){
  report %>% 
    # filter appropriate mass range
    filter(Mass>200 & Mass<900) %>% 
    # remove isotopes
    filter(C13==0) %>% 
    # remove peaks without C assignment
    filter(C>0)
}
compute_indices = function(dat){
  dat %>% 
    dplyr::select(Mass, C:P) %>% 
    dplyr::mutate(AImod = round((1+C-(0.5*O)-S-(0.5*(N+P+H)))/(C-(0.5*O)-S-N-P),4),
                  NOSC =  round(4-(((4*C)+H-(3*N)-(2*O)-(2*S))/C),4),
                  HC = round(H/C,2),
                  OC = round(O/C,2),
                  DBE_AI = 1+C-O-S-0.5*(N+P+H),
                  DBE =  1 + ((2*C-H + N + P))/2,
                  DBE_C = round(DBE_AI/C,4)) %>% 
    dplyr::select(-c(C:P))
}
compute_mol_formula = function(dat){
  dat %>% 
    dplyr::select(Mass, C:P) %>% 
    dplyr::mutate(formula_c = if_else(C>0,paste0("C",C),as.character(NA)),
                  formula_h = if_else(H>0,paste0("H",H),as.character(NA)),
                  formula_o = if_else(O>0,paste0("O",O),as.character(NA)),
                  formula_n = if_else(N>0,paste0("N",N),as.character(NA)),
                  formula_s = if_else(S>0,paste0("S",S),as.character(NA)),
                  formula_p = if_else(P>0,paste0("P",P),as.character(NA)),
                  formula = paste0(formula_c,formula_h, formula_o, formula_n, formula_s, formula_p),
                  formula = str_replace_all(formula,"NA","")) %>% 
    dplyr::select(Mass, formula)
}
assign_class_seidel = function(meta_clean, meta_indices){
  meta_clean %>%
    left_join(meta_indices, by = "Mass") %>% 
    mutate(Class = case_when(AImod>0.66 ~ "condensed aromatic",
                             AImod<=0.66 & AImod > 0.50 ~ "aromatic",
                             AImod <= 0.50 & HC < 1.5 ~ "unsaturated/lignin",
                             HC >= 1.5 ~ "aliphatic"),
           Class = replace_na(Class, "other"),
           Class_detailed = case_when(AImod>0.66 ~ "condensed aromatic",
                                      AImod<=0.66 & AImod > 0.50 ~ "aromatic",
                                      AImod <= 0.50 & HC < 1.5 ~ "unsaturated/lignin",
                                      HC >= 2.0 & OC >= 0.9 ~ "carbohydrate",
                                      HC >= 2.0 & OC < 0.9 ~ "lipid",
                                      HC < 2.0 & HC >= 1.5 & N==0 ~ "aliphatic",
                                      HC < 2.0 & HC >= 1.5 & N > 0 ~ "aliphatic+N")) %>% 
    dplyr::select(Mass, Class, Class_detailed)
}

## for data file
compute_presence = function(dat){
  dat %>% 
    pivot_longer(-("Mass"), values_to = "presence", names_to = "CoreID") %>% 
    # convert intensities to presence==1/absence==0  
    dplyr::mutate(presence = if_else(presence>0,1,0)) %>% 
    # keep only peaks present
    filter(presence>0)
}
apply_replication_filter = function(data_long_key, TREATMENTS){
  max_replicates = 
    data_long_key %>% 
    ungroup() %>% 
    group_by(!!!TREATMENTS) %>% 
    distinct(CoreID) %>% 
    dplyr::summarise(reps = n())
  
  
  # second, join the `reps` file to the long_key file
  # and then use the replication filter  
  data_long_key %>% 
    group_by(formula, !!!TREATMENTS) %>% 
    dplyr::mutate(n = n()) %>% 
    left_join(max_replicates) %>% 
    ungroup() %>% 
    mutate(keep = n >= (2/3)*reps) %>% 
    filter(keep) %>% 
    dplyr::select(-keep, -reps)
  
}

## LEVEL II FUNCTIONS ------------------------------------------------------

make_fticr_meta = function(report){
  fticr_report = (apply_filter_report(report))
  
  meta_clean = 
    fticr_report %>% 
    # select only the relevant columns for the formula assignments
    dplyr::select(Mass, C, H, O, N, S, P, El_comp)
  
  meta_indices = compute_indices(meta_clean)
  meta_formula = compute_mol_formula(meta_clean)
  meta_class = assign_class_seidel(meta_clean, meta_indices)
  
  # output
  meta2 = meta_formula %>% 
    left_join(meta_class, by = "Mass") %>% 
    left_join(meta_indices, by = "Mass") %>% dplyr::select(-Mass) %>% distinct(.)
  
  
  list(meta2 = meta2,
       meta_formula = meta_formula)
}
make_fticr_data = function(report, corekey, TREATMENTS){
  fticr_report = (apply_filter_report(report))
  mass_to_formula = make_fticr_meta(report)$meta_formula
  
  data_columns = fticr_report %>% dplyr::select(Mass, 
                                                starts_with(c("15C", "25C", "45C")))
  
  data_presence = compute_presence(data_columns) %>% 
    left_join(mass_to_formula, by = "Mass") %>% 
    dplyr::select(formula, CoreID, presence) 
    
  
  #  data_presence2 = 
  #    data_presence %>% 
  #    separate(CoreID, sep = "_Alder", into = c("ID", "random1")) %>% 
  #    separate(ID, sep = "_", into = c("Fans", "random2", "DOC_ID")) %>% 
  #    dplyr::select(formula, DOC_ID, presence) %>% 
  #    mutate(DOC_ID = str_replace(DOC_ID, "DOC", "DOC-"))
  
  corekey2 = 
    corekey %>% 
    mutate(CoreID = str_remove(Sample, "TPC_")) %>% 
    dplyr::select(CoreID, Incubation_Temp, A_priori_CUE)
  
  
  data_long_key = data_presence %>% left_join(corekey2, by = "CoreID")
  
  data_long_key_repfiltered = apply_replication_filter(data_long_key, TREATMENTS)
  
  data_long_trt = data_long_key_repfiltered %>% 
    distinct(formula, !!!TREATMENTS)
  
  list(data_long_trt = data_long_trt,
       data_long_key_repfiltered = data_long_key_repfiltered)
  
}

make_lipid_data = function(report, corekey, TREATMENTS){
  fticr_report = (apply_filter_report(report))
  mass_to_formula = make_fticr_meta(report)$meta_formula
  
  data_columns = fticr_report %>% dplyr::select(Mass, 
                                                starts_with(c("15C", "25C", "45C")))
  
  data_presence = compute_presence(data_columns) %>% 
    left_join(mass_to_formula, by = "Mass") %>% 
    dplyr::select(formula, CoreID, presence) %>% 
    separate(CoreID, sep = "_Lip_", into = c("CoreID", "rep")) %>% 
    group_by(formula, CoreID) %>% 
    mutate(n = n()) %>% 
    filter(n >= (2/3)*8) %>% 
    dplyr::select(-n) %>% 
    distinct(formula, CoreID, presence)

  
  corekey2 = 
    corekey %>% 
    rename(CoreID = Lipid_ID) %>% 
    dplyr::select(CoreID, Sample, Incubation_Temp, A_priori_CUE)
  
  
  data_long_key = data_presence %>% left_join(corekey2, by = "CoreID")
  
  data_long_key_repfiltered = apply_replication_filter(data_long_key, TREATMENTS)
  
  data_long_trt = data_long_key_repfiltered %>% 
    distinct(formula, !!!TREATMENTS)
  
  list(data_long_trt = data_long_trt,
       data_long_key_repfiltered = data_long_key_repfiltered)
  
}

make_fticr_data_12T = function(report, corekey, TREATMENTS){
  fticr_report = (apply_filter_report(report))
  mass_to_formula = make_fticr_meta(report)$meta_formula
  
  data_columns = fticr_report %>% dplyr::select(Mass, 
                                                starts_with(c("15C", "25C", "45C")))
  
  data_presence = compute_presence(data_columns) %>% 
    left_join(mass_to_formula, by = "Mass") %>% 
    dplyr::select(formula, CoreID, presence) %>% 
    separate(CoreID, sep = "_Met_", into = c("CoreID", "rep")) %>% 
    group_by(formula, CoreID) %>% 
    mutate(n = n()) %>% 
    filter(n >= (2/3)*8) %>% 
    dplyr::select(-n) %>% 
    distinct(formula, CoreID, presence)
  
  
  corekey2 = 
    corekey %>% 
    rename(CoreID = Lipid_ID) %>% 
    dplyr::select(CoreID, Sample, Incubation_Temp, A_priori_CUE)
  
  
  data_long_key = data_presence %>% left_join(corekey2, by = "CoreID")
  
  data_long_key_repfiltered = apply_replication_filter(data_long_key, TREATMENTS)
  
  data_long_trt = data_long_key_repfiltered %>% 
    distinct(formula, !!!TREATMENTS)
  
  list(data_long_trt = data_long_trt,
       data_long_key_repfiltered = data_long_key_repfiltered)
  
}




# 2. RELATIVE ABUNDANCE COMPUTE FUNCTIONS -------------------------------------------------
compute_relabund_cores = function(fticr_data_longform, fticr_meta, TREATMENTS){
  fticr_data_longform %>% 
    # add the Class column to the data
    left_join(dplyr::select(fticr_meta, formula, Class), by = "formula") %>% 
    # calculate abundance of each Class as the sum of all counts
    group_by(CoreID, Class, !!!TREATMENTS) %>%
    dplyr::summarise(abund = sum(presence)) %>%
    ungroup %>% 
    # create a new column for total counts per core assignment
    # and then calculate relative abundance  
    group_by(CoreID, !!!TREATMENTS) %>% 
    dplyr::mutate(total = sum(abund),
                  relabund  = round((abund/total)*100,2))
}



# 3. MISC FUNCTIONS -------------------------------------------------------

