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
itrn=1000   #number of iteration
nn = 100   #number of data points
par = mat.or.vec(nn, 5)
tx <- rnorm(nn, mean = 0, sd = 1)
#tx=rcauchy(nn, location = 0, scale = 1)
#plot(tx)
ss <- length(tx)
x <<- sort(tx)

##############################################################
#################### Define the function for optimization ################
optimization_fn <- function(q) {

  beta1=complex(real=-0.7,imaginary=2.4)
  beta0=complex(real=cos(pi/6),imaginary=sin(pi/6))
  tyy=Arg(beta0 * (x - beta1) / (x- Conj(beta1))) %% (2 * pi)
  
  kappa=10#1#0.5 #change kappa accordingly 
  ty=(tyy +(rvonmises(ss,mu=0,kappa)%% (2 * pi)))%% (2 * pi) #von Mises error.
  
  # rho=0.3 #change rho accordingly 
  # ty=(tyy +(rwrappedcauchy(ss,mu=0,rho)%% (2 * pi)))%% (2 * pi)  #wrapped Cauchy error.
  
  
  th_prd =numeric(0)
  
  
  
  
  fn=function(d){
    beta0e=complex(real=cos(d[1]),imaginary=sin(d[1]))
    beta1e<-complex(real=d[2],imaginary=d[3])
    
    th_prd <- Arg(beta0e * (x - beta1e) / (x - Conj(beta1e))) %% (2 * pi)
    
    df_th <- (th_prd - ty) %% (2 * pi)
    
    # calling sq.angle is defined elsewhere
    ar <- mean(unlist(lapply(df_th, sq.angle)))
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

