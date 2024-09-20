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
## this is a rough function to subset the allStat list by variable.
## I don't like it in general but it gets the job done and doesn't completely
## wipe the attributes
## it takes `variable` order, family, species, etc.
## and `string` ""
allStat_select = function(variable = NULL, string = NULL,...){
  matches = lapply(allStat$allStat, function(x){
    labs = labels(x)[[1]]
    variable_pos = which(labs == variable)
    unlist(x)[variable_pos] == string
  })
  match_list = allStat$allStat[unlist(matches)]

  return(match_list)
}

# debugonce(allStat_select)
# quick example with porifera
z = allStat_select(variable = 'phylum', string = 'Porifera')

# checking that the names are all still there
labels(z[[1]])


#load population
# load(file = here("data/popStat.Rda"))
