if(!require(rootSolve)){
  install.packages("rootSolve")
}
library(rootSolve)
calc_multi_roots = function(y){
model <- function(x){
  #x[1] = ks, x[2] = J_SEc, x[3] = J_SEn
  sigma = y[[1]]
  M_V = y[[2]]
  k_M = y[[3]]
  Y_SEc = y[[4]]
  B_S = y[[5]]
  Y_SEn = y[[6]]
  ## maintenance cost
  JS = (1+sigma)*k_M*M_V
  # functions to solve for
  A <- B_S*(Y_SEc * x[2] + Y_SEn * x[3] - JS)
  B <- Y_SEc*x[2]*(x[2] + x[3]) - JS * (B_S*x[2] + x[2] + x[3])
  C <- -x[2]*(x[2]+x[3])*JS

  ## find the discriminate
  discriminant = ((B^2) - 4*A*C)

  # Compute both solutions for k_s
  # k_s1 <- (-B + sqrt(discriminant)) / (2 * A)
  # k_s2 <- (-B - sqrt(discriminant)) / (2 * A)
  # ks_disc <- ifelse(k_s1 > 0, k_s1, k_s2)

  F1 <- (-B + sqrt(discriminant)) / (2 * A) - x[1]
  F2 <- x[1] * Y_SEc * (x[1] + x[2])/(x[2]) + B_S * x[1] * Y_SEn * (x[2]*B_S*x[1])/((x[1] + x[2])*(x[2]+x[3]+B_S*x[1])) - JS
  F3 <- ((B_S^2)*(x[1]^2)*Y_SEn*x[2])/(JS - x[1] * Y_SEc * (x[1])) - x[2] - B_S * x[1] - x[3]
  q = c(F1 = F1, F2 = F2, F3 = F3)
  return(q)}
ks = 0.1
J_SEc = 10
J_SEn = 5
(ss <- multiroot(f = model, start = c(ks, J_SEc, J_SEn), maxiter = 1e3))
ss.root = ss$root
# ss.root[1]k_s, [2]J_SEc, [3]J_SEn
return(ss.root)
}
