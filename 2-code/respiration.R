
library(tidyverse)
library(readxl)   
library(data.table)


# load files --------------------------------------------------------------
# the data are spread across multiple tabs (wtf)
## first, import all of the CO2 data ----

# function to import all tabs
read_excel_allsheets <- function(filename, tibble = TRUE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) 
    readxl::read_excel(filename, sheet = X, skip = 1) %>% 
      mutate_all(as.character) %>% 
      mutate(source = X))
  names(x) <- sheets
  x
}

resp = read_excel_allsheets("1-data/respiration copy.xlsx")
# turn into dataframe
resp_df = rbindlist(resp, use.names = FALSE)
resp_df2 = 
  resp_df %>% 
  rename(time_initial = `Initial time`,
         time_final = `end time`,
         minutes_incubated = `Minutes incubated`,
         CO2_ppm = CO2) %>% 
  dplyr::select(Sample, source, minutes_incubated, CO2_ppm)
#

## next, import all background data ----

# the bg data are on row 1, and in col3, col4, or col5. not consistent (wtf)
# so, extract col3:5 of row1 and then clean

read_excel_allsheets_bg <- function(filename, tibble = TRUE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) 
    readxl::read_excel(filename, sheet = X, col_names = FALSE)[1,3:5] %>% 
      mutate(source = X))
  names(x) <- sheets
  x
}

bg = read_excel_allsheets_bg("1-data/respiration copy.xlsx")

bg_df = rbindlist(bg)
bg_df2 = 
  bg_df %>% 
  mutate(bg_combined = paste(`...3`, `...4`, `...5`)) %>% 
  separate(bg_combined, sep = " =", into = c("ignore", "background_ppm")) %>% 
  mutate(background_ppm = str_remove_all(background_ppm, " NA")) %>% 
  mutate(background_ppm = as.numeric(background_ppm)) %>% 
  dplyr::select(source, background_ppm)

#
## finally, combine CO2 with bg ----
respiration = 
  resp_df2 %>% 
  left_join(bg_df2, by = "source")

respiration_processed = 
  respiration %>% 
  mutate(source = str_remove(source, "Day "),
         source = as.numeric(source),
         CO2_ppm = as.numeric(CO2_ppm),
         minutes_incubated = as.numeric(minutes_incubated),
         CO2_bg_corrected = CO2_ppm - background_ppm,
         CO2_min = CO2_bg_corrected/minutes_incubated) %>% 
  rename(day = source) %>% 
  mutate(Sample2 = Sample) %>% 
  separate(Sample2, sep = "_", into = c("vial", "temp", "T"))


respiration_processed %>% 
  ggplot(aes(x = day, y = CO2_min,
             color = vial))+
  #geom_point()+
  geom_line()+
  facet_grid(temp ~ .)


### KFP 2022-04-07. check bg CO2 values!