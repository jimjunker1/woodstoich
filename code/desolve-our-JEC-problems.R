k_m = 0.1
k_c = 0.2
beta_s = 1.5
YSEc = 0.8
YSEn = 0.7
kv_c = 0.1
kv_n = 0.1
Yvc = 0.8
Yvn = 0.7
k_g = 0.5
beta_Ec = 1
beta_En = 1
tolerance = 0.001
M_V = 10
Ec = 2
En = 2

J_S = k_m * M_V

J_Ec_guess = J_S/seq(1, 5, length.out = 50)
# while(abs(J_S_calc - J_S) > tolerance){
  J_Ec = J_Ec_guess
  J_En = J_S - J_Ec
 if(J_En < 0){stop("J_En outside range")}

  theta_Ec = J_Ec/(k_c + J_Ec)
  theta_En = (J_Ec *beta_s * k_c)/((k_c + J_Ec) * (J_Ec + J_En +  beta_s * k_c))
  k_s_temp = J_S / (YSEc * theta_Ec + beta_s * YSEn * theta_En)

  J_S_calc = k_s_temp * (YSEc * theta_Ec + YSEn * theta_En)

  if(J_S_calc - J_S < 0 )

  if (!(f(a) < 0) && (f(b) > 0)) {
    stop('signs of f(a) and f(b) differ')
  } else if ((f(a) > 0) && (f(b) < 0)) {
    stop('signs of f(a) and f(b) differ')
  }

  for (i in 1:n) {
    c <- (a + b) / 2 # Calculate midpoint

    # If the function equals 0 at the midpoint or the midpoint is below the desired tolerance, stop the
    # function and return the root.
    if ((f(c) == 0) || ((b - a) / 2) < tol) {
      return(c)
    }

    # If another iteration is required,
    # check the signs of the function at the points c and a and reassign
    # a or b accordingly as the midpoint to be used in the next iteration.
    ifelse(sign(f(c)) == sign(f(a)),
           a <- c,
           b <- c)
  }

}

theta_Ec = J_Ec/(k_c + J_Ec)
theta_En = (J_Ec *beta_s * k_c)/((k_c + J_Ec) * (J_Ec + J_En +  beta_s * k_c))
k_s = J_S / (YSEc * theta_Ec + beta_s * YSEn * theta_En)
J_EcG = Ec - J_Ec
J_EnG = En - J_En


