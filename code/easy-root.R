library(rootSolve)
library(tidyverse)
#
# model <- function(x){
# # x[1] = rdot; x[2] = JECc
#   Bs = 1
#   kg = 1
#   MV = 0.007
#   kM = 0.25
#   YVEc = 0.99
#   YVEn = 0.99
#   kappa = 0.58
#   sigma = 0  # <- simpilfied for this case
#   JAEcmax = 0.6
#   JAEnmax = 0.8
#   MEnmax = 0.13
#   MEcmax = 0.385
#   MEc = 0.1
#   MEn = 0.04
# # set out the function of fluxes
# F1 <- 1/MV*(kg/(1+(YVEc*(kappa*x[2]-(1+sigma)*kM*MV/YVEc)^(-1)+(YVEc*(kappa*x[2]-(1+sigma)*kM*MV/YVEc))))^(-1)) - x[1]
#
# F2 <- MEc*(JAEcmax/MEcmax)-MV*x[1] - x[2]
# q = c(F1 = F1, F2 = F2)
# return(q)}
# rdotSt = 0.05
# JEcGSt = 0.008
# JEnGSt = 0.004
# (ss <- multiroot(f = model, start = c(rdotSt, JEcGSt), maxiter = 1e4))
#
#
# #
#   model <- function(x){
#     #x[1] = rdot, x[2] = JECc, x[3] = JECn
#     kg = 1
#     MV = 1
#     kM = 0.25
#     YVEc = 0.5
#     # YSEn = 0.01
#     sigma = 0  # <- simpilfied for this case
#     JAEcmax = 0.6
#     JAEnmax = 0.8
#     MEnmax = 0.13
#     MEcmax = 0.385
#     MEc = 0.1
#     MEn = 0.04
#   #Quadratic equation
#  F1 <- 1/MV*(kg/(1+(YVEc*(kappa*x[2])))) - x[1]
#  F2 <- YSEc*JSEc*(JSEc+JSEn) - ((1+sigma)*kM*MV)*(Bs*JSEc+JSEc+JSEn) - x[2]
#  F3 <- -JSEc*(JSEc + JSEn)*((1+sigma)*kM*MV) - x[3]
#  q = c(F1 = F1, F2 = F2, F3 = F3)
#  return(q)}
#   rdot = 0.008
#   JEcGst = 300
#   JEnGst = 100
#   (ss <- multiroot(f = model, start = c(rdot, JEcGst, JEnGst), maxiter = 1e4))
#
#
# rootDf = rootList %>% map(~.x %>% pluck('root') %>% setNames(., nm = c('rdot','JECc','JECn'))) %>% bind_rows %>%
#   bind_cols(parsDf,.)
#
# nrow(rootDf %>% filter(all(rdot > 0 & JECc > 0 & JECn > 0 )))


### simplified root


nboots = 5
parsDf = expand.grid(
            kg = seq(0.1,1, length.out = nboots),
            MV = seq(1,5,length.out = nboots),
            kM = seq(0.1,0.9, length.out = nboots),
            YSEc = 0.8,
            YVEc = 0.8,
            kappa = seq(0.1,0.99, length.out = nboots),
            sigma = 0,
            JAEcmax = 0.002,#0.5,
            JAEnmax = 0.005,#0.9,
            MEcmax = seq(0.1,1e1, length.out = nboots),
            MEnmax = seq(0.1, 1e1, length.out = nboots),
            kEn = 0.5,
            MEc = 1000,
            MEn = 200
            )
parsList = split(parsDf, seq(nrow(parsDf)))

rootList= purrr::map(parsList, \(y){
  model <- function(x){
    #x[1] = rdot, x[2] = JECc, x[3] = JECn
    kg = y[['kg']]
    MV = y[['MV']]
    kM = y[['kM']]
    YSEc = y[['YSEc']]
    kappa = y[['kappa']]
    sigma = y[['sigma']]  # <- simpilfied for this case
    JAEcmax = y[['JAEcmax']]
    JAEnmax = y[['JAEnmax']]
    MEnmax = y[['MEnmax']]
    MEcmax = y[['MEcmax']]
    kEn = y[['kEn']]
    MEc = y[['MEc']]
    MEn = y[['MEn']]
    YVEc = y[['YVEc']]
    v = 0.4
    MV_bracket = 0.1
    JS = (1+sigma)*kM*MV/YSEc
    #Quadratic equation
    F1 <- 1/MV*(kg/(1+(YVEc*(kappa*x[2] - JS))^(-1) + (kappa*x[3])^(-1) + (YVEc*(kappa*x[2] - JS) + kappa*x[3])^(-1))) - x[1]
    F2 <- (MEc*(v*(MV/MV_bracket)^(-0.33)) - x[1]) - x[2]
    F3 <- (MEn*(v*(MV/MV_bracket)^(-0.33)) - x[1]) - x[3]
    q = c(F1 = F1, F2 = F2, F3 = F3)
    return(q)}
rdot = 1
JECc = 5
JECn = 3
  (ss <- multiroot(f = model, start = c(rdot, JECc, JECn), maxiter = 1e2))
})

rootDf = rootList %>% map(~.x %>% pluck('root') %>% setNames(., nm = c('rdot','JECc','JECn'))) %>% bind_rows %>%
  bind_cols(parsDf,.)

nrow(rootDf %>% filter(all(rdot > 0 & JECc > 0 & JECn > 0 )))


