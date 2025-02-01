# if(!require(rootSolve)){
#   install.packages("rootSolve")
# }
library(rootSolve)
library(tidyverse)

# Define parameters

### Stoich
q_V_C = 7
q_V_N = 1
q_R_C = q_V_C
q_R_N = 1
q_Ec_C = 1
q_Ec_N = 0
q_En_C = 4
q_En_N = 1

### Yields
n_V_En = 1
n_V = 1
n_V_Ec = n_V * ((q_V_C  - (q_En_C * n_V_En)) / q_Ec_C)
Y_EcV = ((n_V_Ec * q_Ec_C) / (n_V * q_V_C))
Y_EnV = ((n_V_En * q_En_C) / (n_V * q_V_C))
n_R_En = 1
n_R = 1
n_R_Ec = n_V * ((q_R_C  - (q_En_C * n_R_En)) / q_Ec_C)
Y_EcR = ((n_R_Ec * q_Ec_C) / (n_R * q_R_C))
Y_EnR = ((n_R_En * q_En_C) / (n_R * q_R_C))

### Brown Hare Parms/estimates
M_V_SqBrack = 0.0125523	
nu = 	0.08996
mdens_M =	6.183382
mdens_M_Ec = (mdens_M * 0.6)
mdens_M_En = (mdens_M * 0.4)
kappa = 0.90653
k_M =	0.06313964   
k_J = 0.002
M_Hp = 	11.76545


# Define initial states
M_V = 0.29332253728153584
V = (M_V / M_V_SqBrack)
# M_Ec = mdens_Ec * M_V
M_Ec = 1.2621961871466234
# M_En = mdens_En * M_V
M_En = 0.9747746533685965
M_H = 0.07048299802863989
  
### Guess/arbitrary multireserve parms
rho_S_N = 0.1
rho_S_V = 0.001
k_G = 1
rho_D_N = 0.1
rho_D_H = 0.001
k_R = 1
sigma = 0

# Define pre-determined flux values
J_S = (1 + sigma)*(M_V * k_M)
J_VC = J_S
J_D = (M_H * k_J)
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
# print(somatic_roots$root)

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
J_EcG_pos = J_GV * Y_EcV
J_EnG_pos = J_GV * Y_EnV
J_EcG_neg = J_EcG - J_EcG_pos
J_EnG_neg = J_EnG - J_EnG_pos

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
  J_EcG_pos,
  J_EnG_pos,
  J_EcG_neg,
  J_EnG_neg
  )

print(somatic_fluxes)

# Check that the somatic maintenance fluxes are actually solved correctly
calculated_sommaintenance_fluxes = c(check_som =J_EcS + J_EnS + J_VS,
                                     J_S = J_S)
calculated_sommaintenance_fluxes["difference"] = (calculated_sommaintenance_fluxes[2] - calculated_sommaintenance_fluxes[1])
print(calculated_sommaintenance_fluxes)


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

J_EcD = (((1-kappa) * J_EcC * k_D) / ((1-kappa) * J_EcC + k_D))
J_EnD = (((1-kappa) * J_EnC * k_D * k_D * rho_D_N) / (((1-kappa) * J_EcC + k_D) * ((1-kappa) * J_EcC + (1-kappa) * J_EnC + rho_D_N * k_D)))
J_HD = (k_D* rho_D_H * ((J_HC * k_D * ((1-kappa) * J_EcC + rho_D_N * k_D))/ (((1-kappa) * J_EcC + k_D) * ((1-kappa) * J_EcC + (1-kappa) * J_EnC + rho_D_N * k_D) * ((1-kappa) * J_EcC + (1-kappa) * J_EnC + J_HC + rho_D_H * k_D))))
J_EcR = ((1-kappa) * J_EcC) - J_EcD 
J_EnR = ((1-kappa) * J_EnC) - J_EnD 
J_RR = ((1 / k_R) + (J_EcR/ Y_EcR)^(-1) + (J_EnR/ Y_EnR)^(-1) - ((J_EcR/ Y_EcR)  +  (J_EnR/ Y_EnR))^(-1))^(-1)
J_EcR_pos = J_RR * Y_EcR
J_EnR_pos = J_RR * Y_EnR
J_EcR_neg = J_EcR - J_EcR_pos
J_EnR_neg = J_EnR - J_EnR_pos

maturity_fluxes = c(
  J_EcD,
  J_EnD,
  J_HD,
  J_EcR,
  J_EnR,
  J_RR,
  J_EcR_pos,
  J_EnR_pos,
  J_EcR_neg,
  J_EnR_neg
)

print(maturity_fluxes)

# Check that the maturity maintenance fluxes are actually solved correctly
calculated_devmaintenance_fluxes = c(check_dev =J_EcD + J_EnD + J_HD,
                                     J_D = J_D)
calculated_devmaintenance_fluxes["difference"] = (calculated_devmaintenance_fluxes[2] - calculated_devmaintenance_fluxes[1])
print(calculated_devmaintenance_fluxes)


## Combine all fluxes
fluxes = c(somatic_fluxes,maturity_fluxes)

named_fluxes = c(
  J_EcC=J_EcC,
  J_EnC=J_EnC,
  J_GV=J_GV,
  r_dot=r_dot,
  J_EcS=J_EcS,
  J_EnS=J_EnS,
  J_VS=J_VS,
  J_EcG=J_EcG,
  J_EnG=J_EnG,
  J_EcG_pos=J_EcG_pos,
  J_EnG_pos=J_EnG_pos,
  J_EcG_neg=J_EcG_neg,
  J_EnG_neg=J_EnG_neg,
  J_EcD=J_EcD,
  J_EnD=J_EnD,
  J_HD=J_HD,
  J_EcR=J_EcR,
  J_EnR=J_EnR,
  J_RR=J_RR,
  J_EcR_pos=J_EcR_pos,
  J_EnR_pos=J_EnR_pos,
  J_EcR_neg=J_EcR_neg,
  J_EnR_neg =J_EnR_neg 
)

print(named_fluxes)

# Below: if dM_Vdt < 0, the organism shrinks. If dM_Hdt < 0, it rejuvenates
dM_Vdt = J_GV - J_VS
dM_Hdt = J_RR - J_HD
print(c(dM_Vdt,dM_Hdt))