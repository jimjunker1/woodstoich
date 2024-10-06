if(!require(rootSolve)){
  install.packages("rootSolve")
}
library(rootSolve)
calc_roots = function(y){
  model <- function(x){
    #x[1] = rdot, x[2] = JECc, x[3] = JECn
    k_G = y[1]
    M_V = y[2]
    M_V_SqBrack = y[3]
    k_M = y[4]
    Y_SEc = y[5]
    Y_VEc = y[6]
    kappa = y[7]
    sigma = y[8]  # <- simpilfied for this case
    Mdens_Ec_max = y[9]
    Mdens_En_max = y[10]
    M_Ec = y[11]
    M_En = y[12]
    nu = y[13]
    JS = (1+sigma)*k_M*M_V/Y_SEc
    #Quadratic equation
    F1 <- 1/M_V*(k_G/(1+(Y_VEc*(kappa*x[2] - JS))^(-1) + (kappa*x[3])^(-1) + (Y_VEc*(kappa*x[2] - JS) + kappa*x[3])^(-1))) - x[1]
    F2 <- (M_Ec*(nu*(M_V/M_V_SqBrack)^(-0.33)) - x[1]) - x[2]
    F3 <- (M_En*(nu*(M_V/M_V_SqBrack)^(-0.33)) - x[1]) - x[3]
    q = c(F1 = F1, F2 = F2, F3 = F3)
    return(q)}
  rdot = 1
  JECc = 5
  JECn = 3
  (ss <- multiroot(f = model, start = c(rdot, JECc, JECn), maxiter = 1e2))
  ss.root = ss$root
  #ss.root[1]rdot, [2]JECc,[3]JECn

  return(ss.root)
}
