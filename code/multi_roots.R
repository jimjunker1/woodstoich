library(rootSolve)
library(tidyverse)
nboots = 5
parsDf = expand.grid(
  sigma = 3,
  M_V = seq(10,100,length.out = nboots),
  k_M = seq(0.1,0.9, length.out = nboots),
  Y_SEc = seq(0.1, 0.9, length.out = nboots),
  # k_S = 0.5,
  B_S = 0.3,
  Y_SEn = seq(0.1, 0.5, length.out = nboots)
)

parsList = split(parsDf, seq(nrow(parsDf)))
rootList= purrr::map(parsList, \(y){
# calc_multi_roots = function(){
  model <- function(x){
    #x[1] = ks, x[2] = J_SEc, x[3] = J_SEn
    sigma = y[[1]]
    M_V = y[[2]]
    k_M = y[[3]]
    Y_SEc = y[[4]]
    # k_S = y[[5]]
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
    k_s1 <- (-B + sqrt(discriminant)) / (2 * A)
    k_s2 <- (-B - sqrt(discriminant)) / (2 * A)
    ks_disc <- ifelse(k_s1 > 0, k_s1, k_s2)

    F1 <- ((-B + sqrt((B^2) - 4*A*C))/2*A) - x[1]
    F2 <- x[1] * Y_SEc * (x[1] + x[2])/(x[2]) + B_S * x[1] * Y_SEn * (x[2]*B_S*x[1])/((x[1] + x[2])*(x[2]+x[3]+B_S*x[1])) - JS
    F3 <- ((B_S^2)*(x[1]^2)*Y_SEn*x[2])/(JS - x[1] * Y_SEc * (x[1])) - x[2] - B_S * x[1] - x[3]
    q = c(F1 = F1, F2 = F2, F3 = F3)
    return(q)}
  ks = 0.1
  J_SEc = 10
  J_SEn = 5
  (ss <- multiroot(f = model, start = c(ks, J_SEc, J_SEn), maxiter = 1e2))
  ss.root = ss$root
  return(ss.root)
# }
  })

rootDf = rootList %>% map(~.x %>% setNames(., nm = c('ks','J_SEc','J_SEn'))) %>% bind_rows %>%
  bind_cols(parsDf,.)

rootDf %>% filter(ks > 0 & J_SEc > 0 & J_SEn > 0 ) %>% View()

rootDf %>%
  filter(ks > 0 & J_SEc > 0 & J_SEn > 0 ) %>%
  ggplot(aes(x = ks, y = J_SEn))+
  geom_point()

rootDf %>%
  filter(ks > 0 & J_SEc > 0 & J_SEn > 0 ) %>%
  ggplot(aes(x = ks, y = J_SEc))+
  geom_point()

###
### calculate discriminant from quadratic formula
parsList = split(parsDf, seq(nrow(parsDf)))
discList= purrr::map(parsList, \(y){
  # calc_multi_roots = function(){
  model <- function(x){
    #x[1] = ks, x[2] = J_SEc, x[3] = J_SEn
    sigma = y[[1]]
    M_V = y[[2]]
    k_M = y[[3]]
    Y_SEc = y[[4]]
    k_S = y[[5]]
    B_S = y[[6]]
    Y_SEn = y[[7]]
    ## maintenance cost
    JS = (1+sigma)*k_M*M_V
    # functions to solve for
    A <- B_S*(Y_SEc * x[2] + Y_SEn * x[3] - JS)
    B <- Y_SEc*x[2]*(x[2] + x[3]) - JS * (B_S*x[2] + x[2] + x[3])
    C <- -x[2]*(x[2]+x[3])*JS

    ## find the discriminate
    discriminant = ((B^2) - 4*A*C)

    # Compute both solutions for k_s
    k_s1 <- (-B + sqrt(discriminant)) / (2 * A)
    k_s2 <- (-B - sqrt(discriminant)) / (2 * A)

    F1 <- max(c(((-B + sqrt((B^2) - 4*A*C))/2*A), ((-B - sqrt((B^2) - 4*A*C))/2*A))) - x[1]
    F2 <- x[1] * Y_SEc * (x[1] + x[2])/(x[2]) + B_S * x[1] * Y_SEn * (x[2]*B_S*x[1])/((x[1] + x[2])*(x[2]+x[3]+B_S*x[1])) - JS
    F3 <- ((B_S^2)*(x[1]^2)*Y_SEn*x[2])/(JS - x[1] * Y_SEc * (k_S)) - x[2] - B_S * x[1] - x[3]
    q = c(F1 = F1, F2 = F2, F3 = F3)
    return(q)
    }
  ks = 0.1
  J_SEc = 10
  J_SEn = 5
  (ss <- multiroot(f = model, start = c(ks, J_SEc, J_SEn), maxiter = 1e2))
  ss.root = ss$root
  discmodel <- function(x){
    #x[1] = ks, x[2] = J_SEc, x[3] = J_SEn
    sigma = y[[1]]
    M_V = y[[2]]
    k_M = y[[3]]
    Y_SEc = y[[4]]
    k_S = y[[5]]
    B_S = y[[6]]
    Y_SEn = y[[7]]
    ## maintenance cost
    JS = (1+sigma)*k_M*M_V
    # functions to solve for
    A <- B_S*(Y_SEc * x[2] + Y_SEn * x[3] - JS)
    B <- Y_SEc*x[2]*(x[2] + x[3]) - JS * (B_S*x[2] + x[2] + x[3])
    C <- -x[2]*(x[2]+x[3])*JS

    ## find the discriminate
    discriminant = ((B^2) - 4*A*C)

    # Compute both solutions for k_s
    k_s1 <- (-B + sqrt(discriminant)) / (2 * A)
    k_s2 <- (-B - sqrt(discriminant)) / (2 * A)
    ifelse(k_s1 > 0, ks_disc = k_s1, ks_disc = k_s2)

    F1 <-  ks_disc - x[1]
    F2 <- x[1] * Y_SEc * (x[1] + x[2])/(x[2]) + B_S * x[1] * Y_SEn * (x[2]*B_S*x[1])/((x[1] + x[2])*(x[2]+x[3]+B_S*x[1])) - JS
    F3 <- ((B_S^2)*(x[1]^2)*Y_SEn*x[2])/(JS - x[1] * Y_SEc * (k_S)) - x[2] - B_S * x[1] - x[3]
    q = c(F1 = F1, F2 = F2, F3 = F3)

    return(ks.list = c(
      'ks1' = k_s1,
      'ks2' = k_s2
      ))}

  return(list(
    "ss.root" = ss.root,
    "ks" = discmodel$ks.list
              ))

  # }
})




debugonce(calc_multi_roots)
calc_multi_roots(y)
rootList= purrr::map(parsList, \(y){

# calc_multi_roots = function(y) {

  model <- function(x) {
    # x[1] = rdot, x[2] = JECc, x[3] = JECn
    sigma = y[[1]]
    M_V = y[[2]]
    k_M = y[[3]]
    Y_SEc = y[[4]]
    k_S = y[[5]]
    B_S = y[[6]]
    Y_SEn = y[[7]]

    # Maintenance cost
    JS = (1 + sigma) * k_M * M_V

    # Functions to solve for
    A <- B_S * (Y_SEc * x[2] + Y_SEn * x[3] - JS)
    B <- Y_SEc * x[2] * (x[2] + x[3]) - JS * (B_S * x[2] + x[2] + x[3])
    C <- -x[2] * (x[2] + x[3]) * JS  # Note: Changed from 'c' to 'C' to avoid name conflicts

    # Be careful about complex roots and division by zero
    discriminant <- B^2 - 4 * A * C
    if (discriminant < 0) {
      return(c(F1 = NA, F2 = NA, F3 = NA))  # Return NA if discriminant is negative (no real roots)
    }

    F1 <- (-B + sqrt(discriminant)) / (2 * A) - x[1]

    # Second equation
    F2 <- x[1] * Y_SEc * (x[1] + x[2]) / x[2] + B_S * x[1] * Y_SEn * (x[2] * B_S * x[1]) /
      ((x[1] + x[2]) * (x[2] + x[3] + B_S * x[1])) - JS

    # Third equation
    F3 <- (B_S^2 * (x[1]^2) * Y_SEn * x[2]) / (JS - x[1] * Y_SEc * k_S) - x[2] - B_S * x[1] - x[3]

    # Return the system of equations as a numeric vector
    return(c(F1 = F1, F2 = F2, F3 = F3))
  }

  # Initial guesses for the roots
  ks = 0.1
  J_SEc = 10
  J_SEn = 5

  # Solve using multiroot from rootSolve
  ss <- tryCatch(
    multiroot(f = model, start = c(ks, J_SEc, J_SEn), maxiter = 1e2),
    error = function(e) {
      message("Error in root finding: ", e)
      return(NULL)
    }
  )

  if (!is.null(ss)) {
    ss.root <- ss$root
    # ss.root[1] = rdot, ss.root[2] = JECc, ss.root[3] = JECn
    return(ss.root)
  } else {
    return(NA)  # If root solving fails, return NA
  }
})



#####
# Load the necessary package
library(rootSolve)

# Define the system of equations for J_Ec, J_En, and k_s
model <- function(x, params) {
  # x[1] = J_Ec, x[2] = J_En, x[3] = k_s
  J_Ec <- x[1]
  J_En <- x[2]
  k_s <- x[3]

  # Extract the parameters
  J_S <- params$J_S
  Y_Ec <- params$Y_Ec
  Y_En <- params$Y_En
  B_S <- params$B_S

  # First equation for J_Ec, J_En, and k_s
  f1 <- J_S - k_s * Y_Ec * (k_s + J_Ec) / J_Ec -
    B_S * k_s * Y_En * (J_Ec * B_S * k_s) / ((k_s + J_Ec) * (J_Ec + J_En + B_S * k_s))

  # Second equation for J_En
  f2 <- J_En - (B_S^2 * k_s^2 * Y_En * J_Ec) /
    (J_S - k_s * Y_Ec * (k_s + J_Ec) / J_Ec) + J_Ec + B_S * k_s

  # Return the system of equations as a numeric vector
  return(c(f1 = f1, f2 = f2))
}

# Set the parameters
params <- list(
  J_S = 1.0,     # Example maintenance cost
  Y_Ec = 0.8,    # Yield coefficient for Ec
  Y_En = 0.6,    # Yield coefficient for En
  B_S = 0.5      # Biomass coefficient
)

# Initial guesses for J_Ec, J_En, and k_s
start_values <- c(J_Ec = 0.1, J_En = 0.1, k_s = 0.1)  # Replace with reasonable guesses

# Solve using multiroot
solution <- multiroot(f = model, start = start_values, parms = params)

# Extract the roots
J_Ec_solution <- solution$root[1]
J_En_solution <- solution$root[2]
k_s_solution <- solution$root[3]

# Output the results
cat("Solutions:\n")
cat("J_Ec =", J_Ec_solution, "\n")
cat("J_En =", J_En_solution, "\n")
cat("k_s =", k_s_solution, "\n")
####
