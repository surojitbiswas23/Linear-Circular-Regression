rm(list=ls())

# Load required libraries
library(CircStats)
library(numDeriv)
library(cubature)
library(parallel)
library(foreach)
library(doParallel)
library(circular)


# Set up parallel backend to use available CPU cores
numCores <- detectCores() - 1
cl <- makeCluster(numCores)
registerDoParallel(cl)

# Parameters
itrn=3000
nn = 200
par = mat.or.vec(nn, 5)
tx <- rnorm(nn, mean = 0, sd = 1)
#tx=rcauchy(nn, location = 0, scale = 1)
#plot(tx)
ss <- length(tx)
x <<- sort(tx)

# Define the function for optimization
optimization_fn <- function(q) {
  ##############################################################
  

  #th <<- rvonmises(nn,0,kappa = 1)#rwrappednormal(n, mu = 0, rho = NULL, sd = 1)
  
  #ty=mat.or.vec(ss,1)
  
  
  beta1=complex(real=-0.7,imaginary=2.4)
  beta0=complex(real=cos(pi/6),imaginary=sin(pi/6))
  tyy=Arg(beta0 * (x - beta1) / (x- Conj(beta1))) %% (2 * pi)
  
  kappa=10#1#0.5
  ty=(tyy +(rvonmises(ss,mu=0,kappa)%% (2 * pi)))%% (2 * pi)
  
  # rho=0.3
  # ty=(tyy +(rwrappedcauchy(ss,mu=0,rho)%% (2 * pi)))%% (2 * pi)
  
  
  th_prd =numeric(0)
  
  
  
  r=1
  R=1
  area_com<- function(tt1,tt2,pp1, pp2) {
    
    ph= (pp2-pp1)
    th= ((tt2-tt1)+((r/R)*(sin(tt2)-sin(tt1))))
    return(ph*th/(r/R))
    
  }
  
  
  torus.area_com<-function(phi1,theta1,phi2,theta2){
    
    I1=abs(area_com(theta1, theta2,phi1,phi2))
    I2=abs(area_com(theta1, theta2,2*pi-phi1, phi2))
    I3=abs(area_com(2*pi-theta1, theta2, phi1, phi2))
    I4=abs(area_com(2*pi-theta1, theta2, 2*pi-phi1, phi2))
    I=abs(c(I1,I2,I3,I4))/(4*pi^2*r*R)
    
    min(abs(I[which(I>0)]))
    
  }
  
  
  ########## variance calculation for Zero centered #######  ### Remark: variance is radious invariant
  vartorus<-function(agl){
    #vtorus<-torus.area(0,0,agl,agl)
    vtorus<-torus.area_com(0,0,agl,agl)
  }
  
  fn=function(d){
    beta0e=complex(real=cos(d[1]),imaginary=sin(d[1]))
    beta1e<-complex(real=d[2],imaginary=d[3])
    
    th_prd <- Arg(beta0e * (x - beta1e) / (x - Conj(beta1e))) %% (2 * pi)
    
    df_th <- (th_prd - ty) %% (2 * pi)
    
    # Assuming vartorus is defined elsewhere
    ar <- mean(unlist(lapply(df_th, vartorus)))
    return(ar)
  }
  
  
  

  # Optimization with constraints
  di=c(pi/3,.1,1)
  op=optim(di,fn,lower=c(-pi,-10,0),
           upper=c(pi,10,10),method="L-BFGS-B",hessian = TRUE)
  
  if (op$convergence == 0) {
    e <- eigen(op$hessian)$value
    e2 <- e[e > 0]
    if (length(e2) == length(e)) {
      return(op$par)
    }
  }
  
}
# Measure the computation time using system.time()
start_time <- Sys.time()  # Start time tracking

# Use parallel foreach loop to perform the optimization in parallel
results <- foreach(q = 1:itrn, .combine = rbind, .packages = c("CircStats",'circular', "numDeriv", "cubature")) %dopar% {
  optimization_fn(q)
}

# End time tracking
end_time <- Sys.time()  # End time tracking
computation_time <- end_time - start_time  # Calculate total time taken

# Close the parallel cluster
stopCluster(cl)

# Remove rows with NA values
par <- na.omit(results)

# Calculate estimates and standard errors
b1 = mean(par[, 2])
b2 = mean(par[, 3])
th0 = circ.mean(par[, 1])

se_b1 = sd(par[, 2]) / sqrt(nrow(par))
se_b2 = sd(par[, 3]) / sqrt(nrow(par))

se_th0 = sqrt(var.circular(par[, 1]))

# Output the estimates and their standard errors
cat('Estimated parameters (b1, b2, theta0):\n')
cat(b1, b2, th0, "\n")

cat('Standard errors for (b1, b2, theta0):\n')
cat(se_b1, se_b2,  se_th0, "\n")

# Output the total computation time
cat('Total computation time:', computation_time, '\n')

# 
# Estimated parameters (b1, b2, theta0):
#   1.497631 0.4999023 -2.172562e-05 
# Standard errors for (b1, b2, theta0):
#   0.001116118 0.001140089 0.06922943 
# Total computation time: 13.50721 for kappa=1