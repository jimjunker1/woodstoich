if(!require(rootSolve)){
  install.packages("rootSolve")
  #install.packages("tidyverse")
  #install.packages("rjson")
}
library(rootSolve)
#library(tidyverse)
#library(rjson)
    
  DEB_ODE_calc <- function(netlogo_input){
        
    Food = netlogo_input[1]
    M_V = netlogo_input[2]
    V = netlogo_input[3]
    M_Ec = netlogo_input[4]
    M_En = netlogo_input[5]
    M_H = netlogo_input[6]
    
    q_F_C = netlogo_input[7]
    q_V_C = netlogo_input[8]
    q_R_C = netlogo_input[9]
    q_En_C = netlogo_input[10]
    q_F_N = netlogo_input[11]
    q_V_N = netlogo_input[12]
    q_R_N = netlogo_input[13]
    q_En_N = netlogo_input[14]
    
    J_S = netlogo_input[15]
    J_D = netlogo_input[16]
    
    kappa_Ec = netlogo_input[17]
    kappa_En = netlogo_input[18]
    kappaG_Ec = netlogo_input[19]
    kappaG_En = netlogo_input[20]
    kappaR_Ec = netlogo_input[21]
    kappaR_En = netlogo_input[22]

    F_h = netlogo_input[23]
    J_XA_M_CrlBrack = netlogo_input[24]
    nu = netlogo_input[25]
    mdens_M_Ec = netlogo_input[26]
    mdens_M_En = netlogo_input[27]
    
    Y_EcF = netlogo_input[28]
    Y_EnF = netlogo_input[29]
    rho_S_N = netlogo_input[30]
    rho_S_V = netlogo_input[31]
    Y_EcV = netlogo_input[32]
    Y_EnV = netlogo_input[33]
    rho_D_N = netlogo_input[34]
    rho_D_V = netlogo_input[35]
    Y_EcR = netlogo_input[36]
    Y_EnR = netlogo_input[37]
    M_Hp = netlogo_input[38]
    
    J_EcC_rootguess = netlogo_input[39]
    J_EnC_rootguess = netlogo_input[40]
    kS_rootguess = netlogo_input[41]
    kD_rootguess = netlogo_input[42]
    rdot_rootguess = netlogo_input[43]
    
    # Calculate roots from somatic flux equations: 
    # x[1] = J_EcC
    # x[2] = J_EnC
    # x[3] = k_S
    # x[4] = k_D
    # x[5] = rdot
    
    solve_rootfcns <- function(x){
          ## Perform convenience calculations
          L = V^(1/3)
          J_VC_S = J_S
          J_VC_D = J_D
          
          ## Somatic maintenance equations
          ### Use the root solving function above to determine the SU disassociation rate, k_S:
          
          ### Denominator terms for calculating somatic maintenance terms
          S_A = (kappa_Ec * x[1] + x[3])
          S_B = (kappa_Ec * x[1] + kappa_En * x[2] + rho_S_N * x[3])
          S_C = (kappa_Ec * x[1] + kappa_En * x[2] + J_VC_S + rho_S_V * x[3])
          
          ### Flux calculations
          J_EcS = x[3] * ((kappa_Ec * x[1]) / S_A)
          J_EnS = rho_S_N * x[3] * ((kappa_En * x[2] * x[3]) / (S_A * S_B))
          J_VS = rho_S_V * x[3] * ((J_VC_S * x[3] * (kappa_Ec * x[1] + rho_S_N * x[3]))/ (S_A * S_B * S_C))
          
          ## Calculate somatic growth SU fluxes
          J_EcG = (kappa_Ec * x[1] - J_EcS)
          J_EnG = (kappa_En * x[2] - J_EnS)
          J_GV = (
            (
              ((J_EcG/Y_EcV)^-1) +  
                ((J_EnG/Y_EnV)^-1) - (
                  (J_EcG/Y_EcV) + (J_EnG/Y_EnV))^(-1)
            )
            ^-1)
          
          ## Denominator terms for calculating somatic maintenance terms
          D_A = ((1 - kappa_Ec) * x[1] + x[4])
          D_B = ((1 - kappa_Ec) * x[1] + (1 - kappa_En) * x[2] + rho_D_N * x[4])
          D_C = ((1 - kappa_Ec) * x[1] + (1 - kappa_En) * x[2] + J_VC_D + rho_D_V * x[4])
          
          ## Flux calculations
          J_EcD = x[4] * (((1 - kappa_Ec) * x[1]) / D_A)
          J_EnD = rho_D_N * x[4] *(((1 - kappa_En) * x[2] * x[4]) / (D_A * D_B))
          J_VD = rho_D_V * x[4] * ((J_VC_D * x[4]* ((1 - kappa_Ec) * x[1] + rho_D_N * x[4]))/ (D_A * D_B * D_C))
          
          # EQNS for rootsolver
          F1 = (M_Ec * ((nu / L) - x[5])) - x[1]
          F2 = (M_En * ((nu / L) - x[5])) - x[2]
          F3 = (J_EcS + J_EnS + J_VS) - J_S
          F4 = (J_EcD + J_EnD + J_VD) - J_D
          F5 = ((1 / M_V) * (J_GV - J_VS - J_VD)) - x[5]
          root_fcns = c(F1 = F1, F2 = F2, F3 = F3, F4 = F4, F5 = F5)
          
          return(root_fcns)
        }

    # Perform convenience calculations
    L = V^(1/3)
    J_VC_S = J_S
    J_VC_D = J_D
    
    # Flux calculations
    ## Ingestion flux
    if (Food > 0) {
      f = (Food / (Food + F_h))
      J_FA_max = f * J_XA_M_CrlBrack * L^2
      
      if (J_FA_max > Food) {
        J_FA = Food
      } else {
        J_FA = J_FA_max
      }
    } else {
      J_FA = 0
    }

    ## Assimilation
    J_AEn = (q_F_N / q_F_C) *(q_En_C / q_En_N)* Y_EnF * J_FA
    J_AEc = (Y_EcF * J_FA) - J_AEn
    
    ## Waste
    J_Wc = ((1 - Y_EcF) * J_FA)
    J_Wn = (q_F_N / q_F_C) * ((1 - Y_EnF) * J_FA)
    
    ## Root solve for reserve and structure mobilization
    previous_roots = as.numeric(c(J_EcC_rootguess, J_EnC_rootguess, kS_rootguess, kD_rootguess, rdot_rootguess))
    ODEeqn_roots <- multiroot(f = solve_rootfcns, start = previous_roots)
    
    ## Mobilization fluxes
    J_EcC = (ODEeqn_roots$root)[1]
    J_EnC = (ODEeqn_roots$root)[2]
    
    ## Growth rate
    r_dot = (ODEeqn_roots$root)[5]
    
    ## Somatic maintenance
    k_S_root = (ODEeqn_roots$root)[3]
    J_VC_S = J_S
    
    ### Denominator terms for calculating somatic maintenance terms
    S_A = (kappa_Ec * J_EcC + k_S_root)
    S_B = (kappa_Ec * J_EcC + kappa_En * J_EnC + rho_S_N * k_S_root)
    S_C = (kappa_Ec * J_EcC + kappa_En * J_EnC + J_VC_S + rho_S_V * k_S_root)
    
    ### Flux calculations
    J_EcS = k_S_root * ((kappa_Ec * J_EcC) / S_A)
    J_EnS = rho_S_N * k_S_root * ((kappa_En * J_EnC * k_S_root) / (S_A * S_B))
    J_VS = rho_S_V * k_S_root * ((J_VC_S * k_S_root * (kappa_Ec * J_EcC + rho_S_N * k_S_root))/ (S_A * S_B * S_C))

    ## Calculate somatic growth SU fluxes
    J_EcG = (kappa_Ec * J_EcC - J_EcS)
    J_EnG = (kappa_En * J_EnC - J_EnS)
    J_GV = (
      (
        ((J_EcG/Y_EcV)^-1) +
          ((J_EnG/Y_EnV)^-1) - (
            (J_EcG/Y_EcV) + (J_EnG/Y_EnV))^(-1)
      )
      ^-1)
 
    J_EcG_plus = (J_GV * Y_EcV)
    J_EnG_plus = (J_GV * Y_EnV)
    J_EcG_neg = (J_EcG - J_EcG_plus)
    J_EnG_neg = (J_EnG - J_EnG_plus)
    
    ## Calculate maturity maintenance
    k_D_root = (ODEeqn_roots$root)[4]
    
    ## Denominator terms for calculating somatic maintenance terms
    D_A = ((1 - kappa_Ec) * J_EcC + k_D_root)
    D_B = ((1 - kappa_Ec) * J_EcC + (1 - kappa_En) * J_EnC + rho_D_N * k_D_root)
    D_C = ((1 - kappa_Ec) * J_EcC + (1 - kappa_En) * J_EnC + J_VC_D + rho_D_V * k_D_root)
    
    ## Flux calculations
    J_EcD = k_D_root * (((1 - kappa_Ec) * J_EcC) / D_A)
    J_EnD = rho_D_N * k_D_root *(((1 - kappa_En) * J_EnC * k_D_root) / (D_A * D_B))
    J_VD = rho_D_V * k_D_root * ((J_VC_D * k_D_root* ((1 - kappa_Ec) * J_EcC + rho_D_N * k_D_root))/ (D_A * D_B * D_C))
    
    ## Reproduction/Development fluxes
    J_EcR = ((1 - kappa_Ec) * J_EcC - J_EcD)
    J_EnR = ((1 - kappa_En) * J_EnC - J_EnD)
    J_RR = (
      (
        ((J_EcR/Y_EcR)^-1) +
          ((J_EnR/Y_EnR)^-1) - (
            (J_EcR/Y_EcR) + (J_EnR/Y_EnR))^(-1)
      )
      ^-1)
    J_EcR_plus = (J_RR * Y_EcR)
    J_EnR_plus = (J_RR * Y_EnR)
    J_EcR_neg = (J_EcR - J_EcR_plus)
    J_EnR_neg = (J_EnR - J_EnR_plus)
    
    ## Excretion fluxes
    J_GREc = (kappaG_Ec * J_EcG_neg) + (kappaR_Ec * J_EcR_neg)
    J_GREn = (kappaG_En * J_EnG_neg) + (kappaR_En * J_EnR_neg)
    J_GRXEc = ((1 - kappaG_Ec) * J_EcG_neg) + ((1 - kappaR_Ec) * J_EcR_neg)
    J_GRXEn = ((1 - kappaG_En) * J_EnG_neg) + ((1 - kappaR_En) * J_EnR_neg)
    
    # ODE: CALCULATE CHANGE IN STATE VARIABLES
    
    dFooddt = -J_FA  # total consumption of food from environment
    
    ## Somatic growth
    dM_Vdt = J_GV - J_VS - J_VD
    
    ## Reserves dynamics
    ### Carbon reserve
    dM_Ecdt = J_AEc + J_GREc - J_EcC
    
    ### Nitrogen reserve
    dM_Endt = J_AEn + J_GREn - J_EnC
    
    
    ## Maturity
    dM_Hdt = ifelse(
      M_H < M_Hp,
      J_RR,
      0
    )
    
    ## Reproduction
    dM_Rdt = ifelse(
      M_H < M_Hp,
      0,
      J_RR
    )
    
    ## Excreted nutrients
    
    dX_Cdt = ifelse(M_H < M_Hp,
                    (J_EcS + J_EcD + J_EnS + J_EnD + J_GRXEc + J_GRXEn + J_RR),
                    (J_EcS + J_EcD + J_EnS + J_EnD + J_GRXEc + J_GRXEn)
    )
    
    ### Nitrogen
    
    dX_Ndt = ifelse(M_H < M_Hp,
                    (((q_En_N / q_En_C) * (J_EnS + J_EnD + J_GRXEn)) + ((q_V_N / q_V_C) * (J_VS + J_VD)) + (q_R_N / q_R_C) * J_RR),
                    (((q_En_N / q_En_C) * (J_EnS + J_EnD + J_GRXEn)) + ((q_V_N / q_V_C) * (J_VS + J_VD)))
    )
    
    ## Feces (Unassimilated "Waste")
    ### Carbon
    dW_Cdt = J_Wc
    
    ### Nitrogen
    dW_Ndt = J_Wn
    
    # Store ODE results
    r_calcs_output <- c(dFooddt, dM_Vdt, dM_Ecdt, dM_Endt, dM_Hdt, dM_Rdt, dX_Cdt, dX_Ndt, dW_Cdt, dW_Ndt)
    
    return(c(r_calcs_output, ODEeqn_roots$root))
  }