
if(!require(rootSolve)){
  install.packages("rootSolve")
}
library(rootSolve)
calc_roots_somatic = function(y){
  model <- function(x){
    k_G = y[1]
    M_V = y[2]
    M_V_SqBrack = y[3]
    k_M = y[4]
    kappa = y[5]
    sigma = y[6]
    Mdens_Ec_max = y[7]
    Mdens_En_max = y[8]
    M_Ec = y[9]
    M_En = y[10]
    nu = y[11]
    rho_S_N = y[12]
    rho_S_V = y[13]
    Y_VS = y[14]
    Y_EcV = y[15]
    Y_EnV = y[16]
    J_S = (1+sigma)*k_M*M_V
    # Functions for root finding
    F1 <- (M_Ec*(nu*(M_V/M_V_SqBrack)^(-1/3)) - x[7]) - x[1]
    F2 <- (M_En*(nu*(M_V/M_V_SqBrack)^(-1/3)) - x[7]) - x[2]
    F3 <- (x[6] * ((kappa * x[1])/(((kappa * x[1]) + x[6]))))- x[3]
    F4 <- (x[6] * ((kappa * x[2] * x[6])/(((kappa * x[1]) + x[6]) * ((kappa * x[1]) + (kappa * x[2])  + (rho_S_N * x[6])))))- x[4]
    F5 <- (x[6] * ((kappa * x[2] * x[6])/(((kappa * x[1]) + x[6]) * ((kappa * x[1]) + (kappa * x[2])  + (rho_S_N * x[6])) * ((kappa * x[1]) + (kappa * x[2])  + x[5] + (rho_S_V * x[6])))))- (x[5] * Y_VS)
    F6 <- x[3] + x[4] + (x[5] * Y_VS) - J_S
    F7 <- (1/M_V) * (((((k_G)^-1) + (((x[1] - x[3])/Y_EcV)^-1) + (((x[2] - x[4])/Y_EnV)^-1) + (((x[1] - x[3])/Y_EcV)+((x[2] - x[4])/Y_EnV))^-1)^-1) - x[5])
    q = c(F1 = F1, F2 = F2, F3 = F3, F4 = F4, F5 = F5, F6 = F6, F7 = F7)
    return(q)}

  J_EcC = 0.5
  J_EnC = 0.4
  J_EcS = 0.5
  J_EnS = 0.4
  J_VS = 0.1
  k_S = 0.01
  r_dot = 0

  (ss <- multiroot(f = model, start = c(J_EcC, J_EnC, J_EcS, J_EnS, J_VS, k_S, r_dot), maxiter = 1e2))
  ss.root = ss$root

  return(ss.root)
}
