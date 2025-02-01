if(!require(rootSolve)){
  install.packages("rootSolve")
}
library(rootSolve)

flux_calculations <- function(y){
  
  # Initial state conditions
  M_V = y[1]
  M_H = y[2]
  M_Ec = y[3]
  M_En = y[4]
  
  # Somatic flux parameters
  k_G = y[5]
  M_V_SqBrack = y[6]
  k_M = y[7]
  kappa = y[8]
  sigma = y[9]
  nu = y[10]
  rho_S_N = y[11]
  rho_S_V = y[12]
  Y_EcV = y[13]
  Y_EnV = y[14]
  
  # Maturity/repro flux parameters
  k_R = y[15]
  rho_D_N = y[16]
  rho_D_H = y[17]
  Y_EcR = y[18]
  Y_EnR = y[19]
  k_J = y[20]
  
  J_S = (1 + sigma)*(k_M * M_V)
  J_VC = J_S
  
  J_D = (k_J * M_H)
  J_HC = J_D
  
  solve_somaticterms <- function(x) c(
    F1 = (M_Ec * (nu * ((M_V / M_V_SqBrack)^(-1/3))) - x[5]) - x[1],
    F2 = (M_En * (nu * ((M_V / M_V_SqBrack)^(-1/3))) - x[5]) - x[2],
    F3 = ((kappa * x[1] * x[3]) / (kappa * x[1] + x[3])) +
      ((kappa * x[2] * x[3] * x[3] * rho_S_N) / ((kappa * x[1] + x[3]) * (kappa * x[1] + kappa * x[2] + rho_S_N * x[3]))) +
      (x[3] * rho_S_V * ((J_VC * x[3] * (kappa * x[1] + rho_S_N * x[3]))/ ((kappa * x[1] + x[3]) * (kappa * x[1] + kappa * x[2] + rho_S_N * x[3]) * (kappa * x[1] + kappa * x[2] + J_VC + rho_S_V * x[3])))) -
      J_S,
    F4 = (
      ((1/k_G) + 
         ((((kappa * x[1]) - 
              ((kappa * x[1] * x[3]) / (kappa * x[1] + x[3]))
         )/Y_EcV)^-1)+  
         ((((kappa * x[2]) - 
              ((kappa * x[2] * x[3] * x[3] * rho_S_N) / ((kappa * x[1] + x[3]) * (kappa * x[1] + kappa * x[2] + rho_S_N * x[3])))
         )/Y_EnV)^-1) - (
           (((kappa * x[1]) - 
               ((kappa * x[1] * x[3]) / (kappa * x[1] + x[3]))
           )/Y_EcV) + (((kappa * x[2]) - 
                          ((kappa * x[2] * x[3] * x[3] * rho_S_N) / ((kappa * x[1] + x[3]) * (kappa * x[1] + kappa * x[2] + rho_S_N * x[3])))
           )/Y_EnV))^-1)^-1)- 
      x[4],
    F5 = ((1 / M_V) * (x[4] -
                         (x[3] * rho_S_V * ((J_VC * x[3] * (kappa * x[1] + rho_S_N * x[3]))/ ((kappa * x[1] + x[3]) * (kappa * x[1] + kappa * x[2] + rho_S_N * x[3]) * (kappa * x[1] + kappa * x[2] + J_VC + rho_S_V * x[3]))))
    )) - x[5]
  )
  
  # Set a reasonable initial guess with the correct length (5 elements for 5 equations)
  initial_guesses_somatic <- c(1, 1, 0.1, 1, 1)
  
  # Solve the system
  somatic_roots <- multiroot(f = solve_somaticterms, start = initial_guesses_somatic)
  
  # Calculate the rest of the somatic fluxes
  J_EcC = somatic_roots$root[1]
  J_EnC = somatic_roots$root[2]
  k_S = somatic_roots$root[3]
  J_GV = somatic_roots$root[4]
  r_dot = somatic_roots$root[5]
  J_EcS = ((kappa * J_EcC * k_S) / (kappa * J_EcC + k_S))
  J_EnS = ((kappa * J_EnC * k_S * k_S * rho_S_N) / ((kappa * J_EcC + k_S) * (kappa * J_EcC + kappa * J_EnC + rho_S_N * k_S)))
  J_VS = (k_S * rho_S_V * ((J_VC * k_S * (kappa * J_EcC + rho_S_N * k_S))/ ((kappa * J_EcC + k_S) * (kappa * J_EcC + kappa * J_EnC + rho_S_N * k_S) * (kappa * J_EcC + kappa * J_EnC + J_VC + rho_S_V * k_S))))
  J_EcG = (kappa * J_EcC) - J_EcS
  J_EnG = (kappa * J_EnC) - J_EnS
  J_EcG_plus = J_GV * Y_EcV
  J_EnG_plus = J_GV * Y_EnV
  J_EcG_neg = J_EcG - J_EcG_plus
  J_EnG_neg = J_EnG - J_EnG_plus
  
  somatic_fluxes = c(
    J_EcC,
    J_EnC,
    J_GV,
    r_dot,
    J_EcS,
    J_EnS,
    J_VS,
    J_EcG,
    J_EnG,
    J_EcG_plus,
    J_EnG_plus,
    J_EcG_neg,
    J_EnG_neg
  )
  
  
  #Solve for maturity maintenance
  solve_devterms <- function(x) c(
    F1 = (((1-kappa) * J_EcC * x[1]) / ((1-kappa) * J_EcC + x[1])) +
      (((1-kappa) * J_EnC * x[1] * x[1] * rho_D_N) / (((1-kappa) * J_EcC + x[1]) * ((1-kappa) * J_EcC + (1-kappa) * J_EnC + rho_D_N * x[1]))) +
      (x[1] * rho_D_H * ((J_HC * x[1] * ((1-kappa) * J_EcC + rho_D_N * x[1]))/ (((1-kappa) * J_EcC + x[1]) * ((1-kappa) * J_EcC + (1-kappa) * J_EnC + rho_D_N * x[1]) * ((1-kappa) * J_EcC + (1-kappa) * J_EnC + J_HC + rho_D_H * x[1])))) -
      J_D
  )
  
  # Set a reasonable initial guess with the correct length
  initial_guesses_dev <- c(0)
  
  # Solve the system
  dev_root <- multiroot(f = solve_devterms, start = initial_guesses_dev)
  k_D = dev_root$root
  
  #Solve the rest of the fluxes
  J_EcD = (((1-kappa) * J_EcC * k_D) / ((1-kappa) * J_EcC + k_D))
  J_EnD = (((1-kappa) * J_EnC * k_D * k_D * rho_D_N) / (((1-kappa) * J_EcC + k_D) * ((1-kappa) * J_EcC + (1-kappa) * J_EnC + rho_D_N * k_D)))
  J_HD = (k_D* rho_D_H * ((J_HC * k_D * ((1-kappa) * J_EcC + rho_D_N * k_D))/ (((1-kappa) * J_EcC + k_D) * ((1-kappa) * J_EcC + (1-kappa) * J_EnC + rho_D_N * k_D) * ((1-kappa) * J_EcC + (1-kappa) * J_EnC + J_HC + rho_D_H * k_D))))
  J_EcR = ((1-kappa) * J_EcC) - J_EcD 
  J_EnR = ((1-kappa) * J_EnC) - J_EnD 
  J_RR = ((1 / k_R) + (J_EcR/ Y_EcR)^(-1) + (J_EnR/ Y_EnR)^(-1) - ((J_EcR/ Y_EcR)  +  (J_EnR/ Y_EnR))^(-1))^(-1)
  J_EcR_plus = J_RR * Y_EcR
  J_EnR_plus = J_RR * Y_EnR
  J_EcR_neg = J_EcR - J_EcR_plus
  J_EnR_neg = J_EnR - J_EnR_plus
  
  maturity_fluxes = c(
    J_EcD,
    J_EnD,
    J_HD,
    J_EcR,
    J_EnR,
    J_RR,
    J_EcR_plus,
    J_EnR_plus,
    J_EcR_neg,
    J_EnR_neg
  )
  
  ## Combine all fluxes
  fluxes = c(somatic_fluxes,maturity_fluxes)
  
  return(fluxes)
}