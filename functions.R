
r=1 ###radious of vertical circel 
R=1 ###radious of horizontal circel

## area calculation of curved torus
area_com<- function(tt1,tt2,pp1, pp2) {
  
  ph= (pp2-pp1)
  th= ((tt2-tt1)+((r/R)*(sin(tt2)-sin(tt1))))
  return(ph*th/(r/R))
  
}

## minimum among all four different area.
torus.area_com<-function(phi1,theta1,phi2,theta2){
  
  I1=abs(area_com(theta1, theta2,phi1,phi2))
  I2=abs(area_com(theta1, theta2,2*pi-phi1, phi2))
  I3=abs(area_com(2*pi-theta1, theta2, phi1, phi2))
  I4=abs(area_com(2*pi-theta1, theta2, 2*pi-phi1, phi2))
  I=abs(c(I1,I2,I3,I4))/(4*pi^2*r*R)
  
  min(abs(I[which(I>0)]))
  
}


########## square of an agle calculation
sq.angle<-function(agl){
  vtorus<-torus.area_com(0,0,agl,agl)
}


##################for plotting -pi to pi #########################
zero_slide<-function(x){
  dd=as.numeric(x)
  dd1=(dd-pi)%%(2*pi)-pi
}
