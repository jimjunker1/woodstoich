library(here)
i_am("code/get-amp-pars.R")
source(here("code/get-pars-functions.R"))
library(rvest)
library(tidyverse)
## load the parameter file
load(file = here("data/allStat.Rda"))
# get all species available
allDEB.species <- unlist(labels(allStat$allStat))
allDEB.species <- allDEB.species[1:(length(allDEB.species) -2)]

# get par names for isolating groups
par.names <- unlist(labels(allStat$allStat[[1]]))
head(par.names,10)

# I created a small rough function to subset the dataset
# based on some level e.g., class, order, species, phylum, etc.
# quick example with porifera
z = allStat_select(variable = 'phylum', string = 'Chordata')

# checking that the names are all still there
labels(z[[279]])
z[[279]][24]
## interesting variables
vars = c("species", # species
         "family", #family
         "order",
         "class",
         "phylum",
         "T.typical",
         "J.X.Am", # mass-based max ingestion
         "p.Xm", # energy-based max ingestion
         "J.E.Am", # mass-based max assimilation
         "p.Am", # energy-based max assimilation
         "K", # type II function half saturation
         "F.m", # type II specific searching rate
         "r.B", # specific growth rate. "r.B" = von Bertalanffy growth rates
         "M.V", # volume specific structural mass
         "kap", # fraction of reserve sent to somatic
         "kap.G", # growth efficiency
         # , # fraction of rejected flux reserve from SU
         "m.Em", # mass-based max reserve density
         "E.m", # energy based max reserve density
          # , # max somatic structural body Volume
         "del.V", # max somatic structural body biomass # this is a guess based on fraction of max weight that is structure
         "y.V.E", # yield of structure produced from reserve
         "y.E.X", # yield of reserve produced from food
         "k.M", # somatic maintenance coefficient
         "k.J", # development maintenance rate coefficient
         "v"#, # energy conductance
         # , # binding probability
         )
## Chordates
chordates_allStat = allStat_select(variable = 'phylum', string = 'Chordata')

# select just the variables of interest
debugonce(allStat_get_pars)
chordates_varsList = allStat_get_pars(variables = vars, list = chordates_allStat)

chordates_df = do.call("rbind", chordates_varsList) %>%
  data.frame %>%
  setNames(nm = vars)


#load population
# load(file = here("data/popStat.Rda"))
