
if(!require(rootSolve)){
  install.packages("rootSolve")
}
library(rootSolve)
calc_roots_devrepro = function(y){
  model <- function(x){
    kappa = y[1]
    J_EcC = y[2]
    J_EnC = y[3]
    rho_D_N = y[4]
    rho_D_HR = y[5]
    Y_HRD = y[6]
    Y_EcR = y[7]
    Y_EnR = y[8]
    k_J = y[9]
    M_R = y[10]
    J_D = k_J*M_R
    # Functions for root finding
    F1 <- (x[4] * (((1-kappa) * J_EcC)/((((1-kappa) *J_EcC) + x[4]))))- x[1]
    F2 <- (x[4] * (((1-kappa) * J_EnC * x[4])/((((1-kappa) * J_EcC) + x[4]) * (((1-kappa) * J_EcC) + ((1-kappa) * J_EnC)  + (rho_D_N * x[4])))))- x[2]
    F3 <- (x[4] * (((1-kappa) * J_EnC * x[4])/((((1-kappa) * J_EcC) + x[4]) * (((1-kappa) * J_EcC) + ((1-kappa) * J_EnC)  + (rho_D_N * x[4])) * (((1-kappa) * J_EcC) + ((1-kappa) * J_EnC)  + x[5] + (rho_D_HR * x[4])))))- (x[3] * Y_HRD)
    F4 <- x[1] + x[2] + (x[3] * Y_HRD) - J_D
    q = c(F1 = F1, F2 = F2, F3 = F3, F4 = F4)
    return(q)}
  

  J_EcD = 0.5
  J_EnD = 0.4
  J_HRD = 0.1
  k_D = 0.01

  (ss <- multiroot(f = model, start = c(J_EcD, J_EnD, J_HRD, k_D), maxiter = 1e2))
  ss.root = ss$root

  return(ss.root)
}