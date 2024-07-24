here::i_am("data/01_data-cleaning.R")
library(tidyverse)
library(here)

# clean the stoichiometry resource data from specific papers

## cruz-rivera2000

cruz_rivera2000stoic <- read.table(
  here("data/cruz-rivera2000_stoic.txt"),
  sep = ",", header = TRUE) %>%
  mutate(CNmol = (perC/perN)*(14.007/12.011)) %>%
  select(sourceID, dietID, CNmol)

## demott1998
demott1998stoic <- read.table(
  here("data/demott1998_stoic.txt"),
  sep = ",", header = TRUE) %>%
  pivot_wider(id_cols = c(sourceID, dietID), names_from = stoich_var, values_from = stoich) %>%
  select(sourceID, dietID, CPmol = 'cp')

## ferrao-filho2005

ferrao2005stoic <- read.table(
  here("data/ferrao-filho2005_stoic.txt"),
  sep = ",", header = TRUE) %>%
  pivot_wider(id_cols = c(sourceID, dietID), names_from = stoich_var, values_from = stoich) %>%
  select(sourceID, dietID, CPmol = 'cp', CNmol = 'cn')

## fink2006

fink2006stoic <- read.table(
  here("data/fink2006_stoic.txt"),
  sep = ',', header = TRUE) %>%
  pivot_wider(id_cols = c(sourceID, dietID), names_from = stoich_var, values_from = stoich) %>%
  select(sourceID, dietID, CPmol = 'cp', CNmol = 'cn', NPmol = 'np')

## van geest2007

vangeest2007stoic <- read.table(
  here("data/vangeest2007_stoic.txt"),
  sep = ',', header = TRUE) %>%
  pivot_wider(id_cols = c(sourceID, dietID), names_from = stoich_var, values_from = stoich) %>%
  select(sourceID, dietID, CPmol = 'cp', CNmol = 'cn')

##bind source stoic

source_stoic <- list(cruz_rivera2000stoic,
                     demott1998stoic,
                     ferrao2005stoic,
                     fink2006stoic,
                     vangeest2007stoic) %>%
  map_df(I) %>%
  mutate(across(contains('mol'), ~round(.x, 2)))

## URABE
urabe_cla <- read.csv(here("data/URABE_data_CLAEglm.csv"), header = TRUE) %>%
  mutate(sourceID = 'urabe2017', ingestion_var = 'mL_h-1_ind-1') %>%
  select(sourceID, PC, assimilation_eff = 'AE_C', ingestion_var, ingestion = 'CR_C') %>%
  mutate(PC = PC/1000,
         CPmol = (1/PC)*(30.974/12.011)) %>%
  select(-PC)

urabe_growth <- read.csv(here("data/URABE_data_Growth_s.csv"), header = TRUE) %>%
  mutate(sourceID = 'urabe2017', growth_var = 'd-1', ingestion_var = 'ugC_ind-1_d-1', n = 5) %>%
  select(sourceID, growth_var, growth ='Growth.rate', PC_ratio, NC_ratio, ingestion_var, ingestion = 'IR_C', n) %>%
  mutate(CPmol = 1/(PC_ratio/1000)*(30.974/12.011),
         CNmol = 1/NC_ratio*(14.007/12.011)) %>%
  select(-PC_ratio,-NC_ratio)

# merge the different data frames together
## assimilation
assimilation <- read.table(
  here("data/assimilation.txt"),
  sep = ',', header = TRUE, strip.white = TRUE) %>%
  mutate(dietID = ifelse(is.na(stoich_var), gsub("dietID=(\\w{1,}.)\\;.*", "\\1", notes),NA)) %>%
  merge(source_stoic %>% pivot_longer(cols = contains('mol'),names_to = 'stoich_var', values_to = 'stoich'),
            by = c('sourceID', 'dietID')) %>%
  select(-stoich_var.x, -stoich.x) %>%
  select(sourceID, dietID, stoich_var = 'stoich_var.y', stoich = 'stoich.y', assimilation_eff, notes) %>%
  mutate(assimilation_eff = case_when(assimilation_eff > 2 ~assimilation_eff/100,
                                      between(assimilation_eff,1,2) ~ 1,
                                      .default = assimilation_eff)) %>%
  filter(!grepl("element=P",notes, ignore.case = TRUE)) %>%
  mutate(assimilation_eff = round(assimilation_eff,3))

## growth-eff

growth_eff <- read.table(
  here("data/growth-eff.txt"),
  sep = ',', header = TRUE, strip.white = TRUE) %>%
  mutate(dietID = ifelse(is.na(stoich_var), gsub("dietID=(\\w{1,}.)\\;.*", "\\1", notes),NA)) %>%
  merge(source_stoic %>% pivot_longer(cols = contains('mol'),names_to = 'stoich_var', values_to = 'stoich'),
        by = c('sourceID', 'dietID')) %>%
  select(-stoich_var.x, -stoich.x) %>%
  select(sourceID, dietID, stoich_var = 'stoich_var.y', stoich = 'stoich.y', eff_var, eff, notes) %>%
  filter(!is.na(stoich),
         !grepl("element=P",notes, ignore.case = TRUE)) %>%
  mutate(eff = round(eff,3))

## growth

growth <- read.table(
  here("data/growth.txt"),
  sep = ',', header = TRUE, strip.white = TRUE) %>%
  mutate(dietID = case_when(is.na(stoich_var) & grepl("dietID=(\\w{1,}.)\\;.*", notes) ~ gsub("dietID=(\\w{1,}.)\\;.*", "\\1", notes),
                            is.na(stoich_var) & grepl("dietID=(\\d\\:\\d)\\;.*", notes) ~ gsub("dietID=(\\d\\:\\d)\\;.*", "\\1", notes),
                            is.na(stoich_var) & grepl("dietID=(.P.N)\\;.*", notes) ~ gsub("dietID=(.P.N)\\;.*", "\\1", notes),
                            .default = NA_character_)) %>%
  merge(source_stoic %>% pivot_longer(cols = contains('mol'),names_to = 'stoich_var', values_to = 'stoich'),
        by = c('sourceID', 'dietID')) %>%
  select(-stoich_var.x, -stoich.x) %>%
  select(sourceID, dietID, stoich_var = 'stoich_var.y', stoich = 'stoich.y', growth_var, growth, notes) %>%
  filter(!is.na(stoich),
         !grepl("element=P",notes, ignore.case = TRUE)) %>%
  mutate(growth = ifelse(growth < 0, 0.001, growth),
         growth = round(growth,3))

## ingestion

ingestion <- read.table(
  here("data/ingestion.txt"),
  sep = ",", header = TRUE, strip.white = TRUE) %>%
  mutate(dietID = case_when(is.na(stoich_var) & grepl("dietID=(\\w{1,}.)\\;.*", notes) ~ gsub("dietID=(\\w{1,}.)\\;.*", "\\1", notes),
                            is.na(stoich_var) & grepl("dietID=(\\d\\:\\d)\\;.*", notes) ~ gsub("dietID=(\\d\\:\\d)\\;.*", "\\1", notes),
                            is.na(stoich_var) & grepl("dietID=(.P.N)\\;.*", notes) ~ gsub("dietID=(.P.N)\\;.*", "\\1", notes),
                            .default = NA_character_)) %>%
  merge(source_stoic %>% pivot_longer(cols = contains('mol'),names_to = 'stoich_var', values_to = 'stoich'),
        by = c('sourceID', 'dietID')) %>%
  select(-stoich_var.x, -stoich.x) %>%
  select(sourceID, dietID, stoich_var = 'stoich_var.y', stoich = 'stoich.y', ingestion_var, ingestion, notes) %>%
  filter(!is.na(stoich),
         !grepl("element=P",notes, ignore.case = TRUE)) %>%
  mutate(ingestion = round(ingestion,3))
