#script to model the temperature dependence of growth
#combined with resource competition

#Patrick and Joey
#Feb 15, 2019

#general temp dependence####
Ea <- 0.65 #activation energy
k <- 8.617e-03 #boltzmann constant
Td <- function(temp) { #temperature dependence function
  exp(-Ea/(k*(temp+273.15)))
  }

#intrinsic rate of growth####
e <- matrix(c(0.6, 0.4),1,2) #efficiency of species 1 and 2
a <- matrix(c(0.55, 0.55),1,2) #consumption rate of species 1 and 2
R0 <- 1:200 #supply rate of resource 1 and 2
d <- matrix(c(0.1, 0.1), 1, 2)
r <- function(temp) {rep(e, each = length(temp)) * (Td(temp) %*% a)} #Eq 9 from Fronhofer et al. 2018 
#with R0 removed because r should be independent of resource supply, and deaths removed
matplot(r(1:40), type = 'l', xlab = "*C", ylab = "r")

k_2 = c(50,50) #half saturation resource of growth 
temp <- 40
dn_dt <- (R0 %*% r(temp))/(R0 + k_2) - rep((d * Td(temp)),each = length(R0))
R_star <- (k_2 * (d * Td(temp)))/(r(temp) - d* Td(temp))
matplot(dn_dt, type = 'l', lty = 1)
abline(v = R_star, lty = 2, col = c(1,2))
abline(a = 0, b = 0, lty = 2)
