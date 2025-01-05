
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


vartorus<-function(agl){
 
  vtorus<-torus.area_com(0,0,agl,agl)
}

zero_slide<-function(x){
  dd=as.numeric(x)
  dd1=(dd-pi)%%(2*pi)-pi
}


data<-read.csv("daily_btc_stats_2017_with_avg.csv")
lf=1# 1st jan 2017
rt=214#31st july 2017
ht=data$high_time_radians 


  
 abab=(data$high/ data$low)^(-1)*(data$close-data$open)   
  
 

library(CircStats)
library(numDeriv)

tx<-cl_low


tyy<<-ht[lf:rt]
xx=tx[lf:rt]
x=xx

#pdf(file="histogram_x_bit.pdf")
#postscript(file="true_data_plot_bit.eps")
hist(x,breaks = 50, main='')
#dev.off()

ty<<-tyy

# pdf(file="rose_plot_max_bit_for_regression.pdf")
# # #postscript(file="rose_plot_max_etherum.eps")
# # Plot 4: Circular plot of ty
# plot(as.circular(ty), stack = TRUE, 
#      main = "", 
#      zero = 0,  # Align zero with the positive real axis
#      clockwise = TRUE, 
#      col = "blue", pch = 20)
#  dev.off()


ss<-length(x)


th_prd =numeric(0)

fn=function(d){
  beta0e=complex(real=cos(d[1]),imaginary=sin(d[1]))
  beta1e<-complex(real=d[2],imaginary=d[3])
  
  th_prd <- Arg(beta0e * (x - beta1e) / (x - Conj(beta1e))) %% (2 * pi)
  
  df_th <- (th_prd - ty) %% (2 * pi)
  
  ar <- mean(unlist(lapply(df_th, vartorus)))
  return(ar)
}

it=10
dm=array(data = 0,dim = c(it,9))  
for (j in 1:it) {
  
srtdata=matrix(0,nrow = length(ty), ncol = 2)


conv = 100
while (conv > 0) {
  di <- c(runif(1, 0 * pi, 2 * pi), runif(1, -10, 10), runif(1, 0, 10))
  op = optim(di, fn, lower = c(0, -10, -10),
             upper = c(2 * pi, 10, 10), method = "L-BFGS-B", hessian = TRUE)
  
  # Check if the Hessian has positive eigenvalues to ensure it's invertible
  e <- eigen(op$hessian)$value
  e2 <- e[e > 0]
  conv <- op$convergence
  
  #print(c(length(e2), length(e)))
  
  if (length(e2) < length(e)) {
    conv = 0  # Exit if there are non-positive eigenvalues
  } else {
    # Calculate the standard error if Hessian is positive definite
    hessian_inv <- solve(op$hessian)
    standard_errors <- sqrt(diag(hessian_inv))
  }
}



zf<-mat.or.vec(ss,1)
beta1f<-complex(real=op$par[2],imaginary=op$par[3])
beta0f <- complex(real = cos(op$par[1]), imaginary = sin(op$par[1]))

yhat= Arg(beta0f * (x - beta1f) / (x - Conj(beta1f))) %% (2 * pi)



er=(ty-yhat)%% (2 * pi)



qp=qqplot(ty,yhat,main=paste(j))
sum1=0
for (i in 1:length(ty)) {

  pd=abs(qp$x[i]+qp$y[i])/sqrt(2)

  sum1=sum1+pd
}

cat('At iteration',j,'estimated parameters=',op$par,
    'correlation true vs pred=', cor(ty,yhat),'\n')

cat('standard error=',standard_errors,'\n')


dm[j,]=c(op$par[1],op$par[2],op$par[3],
         standard_errors[1], standard_errors[2], standard_errors[3],sum(standard_errors^2),cor(ty,yhat),sum1) ###storing par, se, corr
# 


wt=watson.test(er,alpha = 0.05,dist = 'vonmises')
print(wt)
par(mfrow=c(2,2))
#pdf(file="true_data_plot_bit.pdf")
#postscript(file="true_data_plot_bit.eps")
plot((ty), ylim = c(0*pi, 2*pi), yaxt = "n", ylab ='Observed',xlab='n')
# Customize the y-axis ticks
axis(2, at = c( 0, pi,2*pi), 
     labels = expression(0, pi,2*pi))
#dev.off()



#pdf(file="predicted_data_plot_bit.pdf")
#postscript(file="predicted_data_plot_bit.eps")
plot(yhat, ylim = c(0*pi, 2*pi), yaxt = "n", ylab ='Predicted',xlab='n')
# Customize the y-axis ticks
axis(2, at = c( 0, pi,2*pi), 
     labels = expression(0, pi,2*pi))
#dev.off()




#pdf(file="error_plot_bit.pdf")
#postscript(file="error_plot_bit.eps")
plot(zero_slide(as.numeric(er)), ylim = c(-pi, pi), yaxt = "n", ylab ='Residual',xlab='n')
# Customize the y-axis ticks
axis(2, at = c( -pi, 0,pi), 
     labels = expression(-pi, 0,pi))
abline(h=0,col=6,lwd=2)
legend("topright", legend = 'Resudual=0', 
       col = 6, lwd = 2)
#dev.off()

#pdf(file="QQ_plot_dataset_bit.pdf")
#postscript(file="QQ_plot_dataset_bit.eps")
qqplot((ty),(yhat),ylim = c(0, 2*pi),
       xlim = c(0, 2*pi),
       yaxt = "n",  # Suppress default y-axis
       xaxt = "n",  # Suppress default x-axis
       ylab ='Predicted',
       xlab='Observed')  # X-axis label with pi symbol)
abline(a=0,b=1,col=6,lwd=2)
# Add custom y-axis ticks with pi symbols
axis(2, at = c( 0, pi,2*pi), 
     labels = expression(0, pi,2*pi))

# Add custom x-axis ticks with pi symbols
axis(1, at = c( 0, pi,2*pi), 
     labels = expression(0, pi,2*pi))
abline(a=0,b=1,col=6,lwd=2)
#dev.off()


}


par(mfrow=c(1,1))
cat('parameters and error:', dm[which.min(dm[,7]),])

beta1fp<-complex(real= dm[which.min(dm[,7]),2],imaginary= dm[which.min(dm[,7]),3])
beta0fp<- complex(real = cos( dm[which.min(dm[,7]),1]), imaginary = sin( dm[which.min(dm[,7]),1]))

typ= Arg(beta0fp * (x - beta1fp) / (x - Conj(beta1fp))) %% (2 * pi)

par(mfrow=c(2,2))
#pdf(file="true_data_plot_bit.pdf")
#postscript(file="true_data_plot_bit.eps")
plot((ty), ylim = c(0*pi, 2*pi), yaxt = "n", ylab ='Observed',xlab='n')
# Customize the y-axis ticks
axis(2, at = c( 0, pi,2*pi), 
     labels = expression(0, pi,2*pi))
#dev.off()



#pdf(file="predicted_data_plot_bit.pdf")
#postscript(file="predicted_data_plot_bit.eps")
plot(typ, ylim = c(0*pi, 2*pi), yaxt = "n", ylab ='Predicted',xlab='n')
# Customize the y-axis ticks
axis(2, at = c( 0, pi,2*pi), 
     labels = expression(0, pi,2*pi))
#dev.off()



#pdf(file="error_plot_bit.pdf")
#postscript(file="error_plot_bit.eps")
plot(zero_slide(as.numeric(er)), ylim = c(-pi, pi), yaxt = "n", ylab ='Residual',xlab='n')
# Customize the y-axis ticks
axis(2, at = c( -pi, 0,pi), 
     labels = expression(-pi, 0,pi))
abline(h=0,col=6,lwd=2)

legend("topleft", legend = 'Resudual=0', 
       col = 6, lwd = 2)
#dev.off()

#pdf(file="QQ_plot_dataset_bit.pdf")
#postscript(file="QQ_plot_dataset_bit.eps")
qqplot((ty),(typ),ylim = c(0, 2*pi),
       xlim = c(0, 2*pi),
       yaxt = "n",  # Suppress default y-axis
       xaxt = "n",  # Suppress default x-axis
       ylab ='Predicted',
       xlab='Observed')  # X-axis label with pi symbol)
abline(a=0,b=1,col=6,lwd=2)
# Add custom y-axis ticks with pi symbols
axis(2, at = c( 0, pi,2*pi), 
     labels = expression(0, pi,2*pi))

# Add custom x-axis ticks with pi symbols
axis(1, at = c( 0, pi,2*pi), 
     labels = expression(0, pi,2*pi))
abline(a=0,b=1,col=6,lwd=2)
#dev.off()



