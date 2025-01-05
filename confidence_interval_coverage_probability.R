
# Measure the computation time using system.time()
start_time <- date()  # Start time tracking
# Load required libraries
library(CircStats)
library(numDeriv)
library(cubature)
library(parallel)
library(foreach)
library(doParallel)
library(circular)

# Set up parallel backend to use available CPU cores
numCores <- detectCores() 
cl <- makeCluster(numCores)
registerDoParallel(cl)

# Parameters
itrn=1000
nn = 500
boot=500
par = mat.or.vec(nn, 5)
tx <- rnorm(nn,  0, 1)
#plot(tx)
n <- length(tx)
x <<- sort(tx)
xnew=3

# Define the function for optimization
beta1=complex(real=1.1,imaginary=.7)
beta0=complex(real=cos(0),imaginary=sin(0))
ty_org=Arg(beta0 * (xnew - beta1) / (xnew- Conj(beta1))) %% (2 * pi) 
rho=3


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


cover_ci<-function(q)
{
  
  tyy=Arg(beta0 * (x - beta1) / (x- Conj(beta1))) %% (2 * pi)
  
  th_prd =numeric(0)
  library(circular)
  ty= (tyy +(rvonmises(n,mu=0,rho)%% (2 * pi)))%% (2 * pi)
  
  
  
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
  di=c(0,.1,0.2)
  op=optim(di,fn,lower=c(-pi,-10,0),
           upper=c(pi,10,10),method="L-BFGS-B",hessian = TRUE)
  
  p<- if (op$convergence == 0 && all(eigen(op$hessian)$values > 0)) {
    op$par  # If optimization converged and all eigenvalues are positive, use the parameters
  } else {
    c(0, 0, 0)  # Otherwise, use the fallback vector
  }
  
  #cat('estimated parameters vcalues=', p,'\n')
  
  beta0f=complex(real=cos(p[1]),imaginary=sin(p[1]))
  beta1f<-complex(real=p[2],imaginary=p[3])
  yhat0= Arg(beta0f * (x - beta1f) / (x - Conj(beta1f))) %% (2 * pi)
  
  er0=(ty-yhat0)%% (2 * pi)
  
  #plot(as.numeric(er0))
  
  
  tynew=Arg(beta0f * (xnew - beta1f) / (xnew- Conj(beta1f))) %% (2 * pi)
  
 # mat=matrix(0,nrow = itrn,ncol = 2)
  covg <- numeric(itrn)
  #covgp <- numeric(itrn)
  #tynew<- numeric(itrn)
 # intvmat<-array(0,dim = c(itrn,4))
  yhat_boot=numeric(0)
  
  
  for (j in 1:boot) {
    
    
    th_prd1 =numeric(0)
    ty= (yhat0 +sample(er0,replace = T))%% (2 * pi)
    
    
    
    fn=function(d){
      beta0e=complex(real=cos(d[1]),imaginary=sin(d[1]))
      beta1e<-complex(real=d[2],imaginary=d[3])
      
      th_prd1 <- Arg(beta0e * (x - beta1e) / (x - Conj(beta1e))) %% (2 * pi)
      
      df_th <- (th_prd1 - ty) %% (2 * pi)
      
      # Assuming vartorus is defined elsewhere
      ar <- mean(unlist(lapply(df_th, vartorus)))
      return(ar)
    }
    
    
    
    
    # Optimization with constraints
    di=c(0,.1,0.2)
    op_boot=optim(di,fn,lower=c(-pi,-10,0),
             upper=c(pi,10,10),method="L-BFGS-B",hessian = TRUE)
    
    p1 <- if (op_boot$convergence == 0 && all(eigen(op_boot$hessian)$values > 0)) {
      op_boot$par  # If optimization converged and all eigenvalues are positive, use the parameters
    } else {
      c(0, 0, 0)  # Otherwise, use the fallback vector
    }
    
    #cat('estimated parameters vcalues=', p,'\n')
    
    beta0fb=complex(real=cos(p1[1]),imaginary=sin(p1[1]))
    beta1fb<-complex(real=p1[2],imaginary=p1[3])
    yhat_boot[j]= Arg(beta0fb * (xnew - beta1fb) / (xnew- Conj(beta1fb))) %% (2 * pi)
    
  }  ####end of for loop of bootstarp
  
    cnf=(quantile.circular(yhat_boot,probs = c(0.025,0.975)))%% (2 * pi)
    cnf_min=min(cnf)
    cnf_max=max(cnf)#-2*pi
    
    agldata<-array(0,dim = c(2,4))
    agldata[1,]<-c("L","R","cm_s","mu")
    
    
    cm_s<-tynew
    mu=ty_org
    L<-cnf_min#min(mat[i,1], mat[i,2])
    R<-cnf_max#max(mat[i,1], mat[i,2])
    
    agldata[2,]<-as.numeric(c(L,R,cm_s,mu))
    ss<-agldata[1,order(agldata[2,])]
    covg<-(ss[2]=="L")*(ss[3]=="R")+(ss[1]=="L")*
      (ss[4]=="R")+(ss[1]=="L")*(ss[2]=="R")+(ss[3]=="L")*(ss[4]=="R")
    # 
    # covgp[j]<-mean(covg)
    # intvmat[i,]<-c(mat[i,1],mat[i,2],cm_s,mu)
   # print(covgp[j])
    #cat(i, tynew[i],mu,"\n")
    return(c(p[1],p[2],p[3],L,R,cm_s,mu,covg))
  }








# Use parallel foreach loop to perform the optimization in parallel
results_ci <- foreach(q = 1:itrn, .combine = rbind, .packages = c("CircStats", "numDeriv", "cubature")) %dopar% {
  cover_ci(q)
}

# End time tracking
end_time <- date()  # End time tracking
#computation_time <- end_time - start_time  # Calculate total time taken

# Close the parallel cluster
stopCluster(cl)

# Output the total computation time
cat('start time:', start_time, '\n')
cat('end  time:', end_time, '\n')

#plot(covgp, type='l', ylim = c(0.9,1.01))
#abline(h=0.95)
cat("COVERGE PROBABILITY=",mean(results_ci[,8]),"\n")
cat("does not converges for", itrn*mean((results_ci[,1]==0)*(results_ci[,2]==0)*(results_ci[,3]==0)),'times of total iteration')

res_matrix= array(0, dim = c(itrn,8))#matrix(0,nrow = itrn,ncol = 8)# 
for (i in 1:8) {
  
  res_matrix[,i]=results_ci[,i]
}

write.csv(x = res_matrix,file = 'data matrix of coverage probability_von.csv')

cov_pro=read.csv('data matrix of coverage probability_von.csv')
nzt=sum(cov_pro$V6==0)
sum((cov_pro$V8))/(1000-nzt)














