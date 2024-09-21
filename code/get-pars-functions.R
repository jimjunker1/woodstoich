getDEB.species <- function() {
  require(rvest)
  library(rvest)
  url <- "https://www.bio.vu.nl/thb/deb/deblab/add_my_pet/species_list.html"
  d1 <- read_html(url)

  phylum <- d1 %>% html_nodes("td:nth-child(1)") %>% html_text()
  class <- d1 %>% html_nodes("td:nth-child(2)") %>% html_text()
  order <- d1 %>% html_nodes("td:nth-child(3)") %>% html_text()
  family <- d1 %>% html_nodes("td:nth-child(4)") %>% html_text()
  species <- d1 %>% html_nodes("td:nth-child(5)") %>% html_text()
  common <- d1 %>% html_nodes("td:nth-child(6)") %>% html_text()
  type <- d1 %>% html_nodes("td:nth-child(7)") %>% html_text()
  mre <- d1 %>% html_nodes("td:nth-child(8)") %>% html_text()
  smre <- d1 %>% html_nodes("td:nth-child(9)") %>% html_text()
  complete <- d1 %>% html_nodes("td:nth-child(10)") %>% html_text()
  all.species <- as.data.frame(cbind(phylum, class, order,
                                     family, species, common, type, mre, smre, complete),
                               stringsAsFactors = FALSE)
  all.species$species <- gsub(" ", "_", all.species$species)
  all.species$mre <- as.numeric(mre)
  all.species$smre <- as.numeric(smre)
  all.species$complete <- as.numeric(complete)
  return(all.species)
}

getDEB.pars <- function(species) {
  require(rvest)
  library(rvest)
  baseurl <- "https://www.bio.vu.nl/thb/deb/deblab/add_my_pet/entries_web/"
  d1 <- read_html(paste0(baseurl, species, "/", species, "_par.html"))
  symbol1 <- d1 %>% html_nodes("td:nth-child(1)") %>% html_text()

  value1 <- d1 %>% html_nodes("td:nth-child(2)") %>% html_text()

  units1 <- d1 %>% html_nodes("td:nth-child(3)") %>% html_text()

  description1 <- d1 %>% html_nodes("td:nth-child(4)") %>%
    html_text()

  extra1 <- d1 %>% html_nodes("td:nth-child(5)") %>% html_text()

  extra2 <- d1 %>% html_nodes("td:nth-child(6)") %>% html_text()
  end <- which(symbol1 == "T_ref")
  symbol <- symbol1[1:end]
  value <- value1[1:end]
  units <- units1[1:end]
  description <- description1[1:end]

  pars <- as.data.frame(cbind(symbol, value, units, description))
  pars$symbol <- as.character(symbol)
  pars$value <- as.numeric(value)
  pars$units <- as.character(units)
  pars$description <- as.character(description)

  chempot <- c(value1[end + 1], units1[end + 1], description1[end +
                                                                1], extra1[1])
  dens <- c(value1[end + 2], units1[end + 2], description1[end +
                                                             2], extra1[2])
  org.C <- c(units1[end + 3], description1[end + 3], extra1[3],
             extra2[1])
  org.H <- c(value1[end + 4], units1[end + 4], description1[end +
                                                              4], extra1[4])
  org.O <- c(value1[end + 5], units1[end + 5], description1[end +
                                                              5], extra1[5])
  org.N <- c(value1[end + 6], units1[end + 6], description1[end +
                                                              6], extra1[6])
  min.C <- c(units1[end + 7], description1[end + 7], extra1[7],
             extra2[2])
  min.H <- c(value1[end + 8], units1[end + 8], description1[end +
                                                              8], extra1[8])
  min.O <- c(value1[end + 9], units1[end + 9], description1[end +
                                                              9], extra1[9])
  min.N <- c(value1[end + 10], units1[end + 10], description1[end +
                                                                10], extra1[10])

  organics <- rbind(org.C, org.H, org.O, org.N)
  minerals <- rbind(min.C, min.H, min.O, min.N)
  colnames(organics) <- c("X", "V", "E", "P")
  colnames(minerals) <- c("CO2", "H2O", "O2", "N-waste")
  rownames(organics) <- c("C", "H", "O", "N")
  rownames(minerals) <- c("C", "H", "O", "N")
  class(chempot) <- "numeric"
  class(dens) <- "numeric"
  class(organics) <- "numeric"
  class(minerals) <- "numeric"

  return(list(pars = pars, chempot = chempot, dens = dens,
              organics = organics, minerals = minerals))
}

getDEB.implied <- function(species) {
  require(rvest)
  library(rvest)
  baseurl <- "https://www.bio.vu.nl/thb/deb/deblab/add_my_pet/entries_web/"
  d1 <- read_html(paste0(baseurl, species, "/", species, "_stat.html"))
  symbol <- d1 %>% html_nodes("td:nth-child(1)") %>% html_text()

  value <- d1 %>% html_nodes("td:nth-child(2)") %>% html_text()

  units <- d1 %>% html_nodes("td:nth-child(3)") %>% html_text()

  description <- d1 %>% html_nodes("td:nth-child(4)") %>% html_text()

  final <- as.data.frame(cbind(symbol, value, units, description))
  final$symbol <- as.character(symbol)
  final$value <- as.numeric(value)
  final$units <- as.character(units)
  final$description <- as.character(description)
  return(final)
}

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

##
allStat_get_pars = function(variables = NULL, list = NULL,...){
  library(purrr)
  # library(rlist)
  cleanList = map(list, \(x){
    pList = x
    # pList = pluck(x,1)
    var_vec = which(unlist(labels(pList)) %in% variables)
    pList[var_vec]
  })
  return(cleanList)
}

