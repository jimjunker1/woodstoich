library(rootSolve)
library(tidyverse)

# Define terms
q_V_C = 4.5
q_V_N = 1
q_Ec_C = 1
q_Ec_N = 0
q_En_C = 4
q_En_N = 1
n_V_En = 1
n_V = 1
n_V_Ec = ((q_V_C * n_V) - (q_En_C * n_V_En) / q_Ec_C)
Y_EcV = ((n_V_Ec * q_Ec_C) / (n_V * q_V_C))
Y_EnV = ((n_V_En * q_En_C) / (n_V * q_V_C))
Y_EcF = 0.6
Y_EnF = 0.8
M_V_SqBrack = 0.01034
M_V = 0.3
nu = 0.03993
i_M_CrlBrack = 0.00461 
J_AEc_CrlBrack = i_M_CrlBrack * Y_EcF
J_AEn_CrlBrack = i_M_CrlBrack * Y_EnF
kappaG_Ec = 0.8
kappaG_En = 0.8
Mdens_Ec_M = J_AEc_CrlBrack/((1 - kappaG_Ec)*M_V_SqBrack*nu)
Mdens_En_M = J_AEn_CrlBrack/((1 - kappaG_En)*M_V_SqBrack*nu)
Mdens_Ec = Mdens_Ec_M
Mdens_En = Mdens_En_M
M_Ec = Mdens_Ec * M_V
M_En = Mdens_En * M_V
k_M = 0.54772
J_S = (M_V * k_M)
J_VC = J_S
kappa = 0.84233
rho_S_N = 0.5
rho_S_V = 0.05
k_G = 1

# Define the function to find roots
solve_somaticterms <- function(x) c(
  F1 <- ((kappa *x[7] * x[4]) / (kappa * x[7] + x[4])) - x[1],
  F2 <- ((kappa * x[8] * x[4] * x[4] * rho_S_N) / ((kappa * x[7] + x[4]) * (kappa *x[7] + kappa * x[8] + rho_S_N * x[4]))) - x[2],
  F3 <- (x[4] * rho_S_V * ((J_VC * x[4] * (kappa * x[7] + rho_S_N * x[4]))/ ((kappa * x[7] + x[4]) * (kappa * x[7] + kappa * x[8] + rho_S_N * x[4]) * (kappa * x[7] + kappa * x[8] + J_VC + rho_S_V * x[4])))) - x[3],
  F4 <- x[1] + x[2] + x[3] - J_S,
  F5 <- (((1/k_G) + ((((kappa * x[7]) - x[1])/Y_EcV)^-1)+  ((((kappa * x[8]) - x[2])/Y_EnV)^-1) - ((((kappa * x[7]) - x[1])/Y_EcV) + (((kappa * x[8]) - x[2])/Y_EnV))^-1)^-1) -x[5],
  F6 <- ((1 / M_V) * (x[5] -x[3])) - x[6],
  F7 <- (M_Ec * (nu * ((M_V / M_V_SqBrack)^(-1/3))) - x[6]) - x[7],
  F8 <- (M_En * (nu * ((M_V / M_V_SqBrack)^(-1/3))) - x[6]) - x[8]
  )
# Use multiroot from rootSolve with an initial guess for k_S
# initial_guesses = c(1,1,0.1,1,0)  # Set a reasonable initial guess for k_S
initial_guesses = c(1,1,0,0.1,1,0,1,1)  # Set a reasonable initial guess for k_S
somatic_roots <- multiroot(f = solve_somaticterms, start = initial_guesses)

somatic_terms = c(
J_EcS = somatic_roots$root[1],
J_EnS = somatic_roots$root[2],
J_VS = somatic_roots$root[3],
k_S = somatic_roots$root[4],
J_GV = somatic_roots$root[5],
r_dot = somatic_roots$root[6],
J_EcC = somatic_roots$root[7],
J_EnC = somatic_roots$root[8]
)

print(somatic_terms)