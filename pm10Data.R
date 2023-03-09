## Supplementary R code for the paper 
## "A general framework for circular local likelihood regression" 

# Workflow for obtaining right panels of Figure 1
# and statistical analysis related to the pm10 data

################
library(plotrix)
################

###########################################
source("Functions_implementing_methods.R")
###########################################


dat<-read.table("pm10.txt")			# Read data

dir<-circular(dat$direction,units="degrees")   # Wind direction in degrees
X1<-conversion.circular(dir,units="radians")   # Wind direction in radians
X2<-dat$speed					     # Wind speed
Y<-dat$pm10						     # Pm10 concentration


nt<-500									# Number of points where the target function is to be estimated
xseq<-circular(seq(0,2*pi,length=nt),rotation="clock")	# Points where the target function is to be estimated

startv<-c(log(mean(Y)),0,0,0)		# Starting values for the estimation process


lowerk<-0.3		# Minimum value of kappa for the concentration parameter selection search
upperk<-100		# Maximum value of kappa for the concentration parameter selection search
lowerk_ast<-0.3	# Minimum value of kappa star for the concentration parameter selection search
upperk_ast<-15	# Maximum value of kappa star for the concentration parameter selection search



# Selection of the smoothing parameter (Refined procedure)

kappa_ast_p3<-optimize(function(kappa)gamma_CRSC_int(xseq,X1,Y,kappa,startv),lower=lowerk_ast,upper=upperk_ast)$minimum
kappa_ast_p3<-kappa_ast_p3*(192/167)^(2/9)   # Pilot kappa
kappa_opt<-optimize(function(kappa)gamma_MSE_int(xseq,X1,Y,kappa,kappa_ast_p3,startv),lower=lowerk,upper=upperk)$minimum # Refined value of kappa


# Estimation
modREFINED<-gamma_p1(xseq,X1,Y,kappa_opt,startv)[[1]]


# Bias and variance
sesgo<-bias_gamma_p3(xseq,X1,Y,kappa=kappa_opt,kappa_ast=kappa_ast_p3,startv)
varianza<-variance_gamma_p1(xseq,X1,Y,kappa=kappa_opt,kappa_ast_p3,startv)


## Kernel weighted bias 
sesgo_kernel<-matrix(nrow=nt,ncol=nt)
for (i in 1:nt){
	sesgo_kernel[i,]<-sesgo*dvonmises(xseq,xseq[i],kappa=kappa_opt)
}
sesgok_final<-numeric(nt)
for (i in 1:nt){
	sesgok_final[i]<-sintegral(xseq,sesgo_kernel[i,])$int
}


## Kernel weighted variance
var_kernel<-matrix(nrow=nt,ncol=nt)
for (i in 1:nt){
	var_kernel[i,]<-varianza*dvonmises(xseq,xseq[i],kappa=kappa_opt)
}
vark_final<-numeric(nt)
for (i in 1:nt){
	vark_final[i]<-sintegral(xseq,var_kernel[i,])$int
}



######################################################
# 	PLANAR REPRESENTATION (Figure 3: bottom right)
######################################################

plot(as.numeric(X1),Y,pch=1,xlim=c(0,2*pi),xlab="Wind direction",ylab="pm10",xaxt="n")
axis(1,c(0,pi/2,pi,3*pi/2,2*pi),c(expression(0),expression(pi/2),expression(pi),expression(3*pi/2),expression(2*pi)))
points(as.numeric(xseq),exp(modREFINED),type="l",col=2,lwd=3)

# Confidence intervals
int_sup<-modREFINED-sesgok_final+qnorm(1-0.05/2)*sqrt(vark_final)		
int_inf<-modREFINED-sesgok_final-qnorm(1-0.05/2)*sqrt(vark_final)

points(as.numeric(xseq),exp(int_sup),type="l",col="lightblue",lwd=3)
points(as.numeric(xseq),exp(int_inf),type="l",col="lightblue",lwd=3)

for (i in 1:nt){
	segments(xseq[i],exp(int_inf[i]),xseq[i],exp(int_sup[i]),lwd=2,col="lightblue")
}

points(as.numeric(xseq),exp(modREFINED),type="l",col=2,lwd=3)
points(as.numeric(X1),Y,pch=1)


### Zoomed version (the one in Figure 1, middle)

plot(as.numeric(X1),Y,pch=1,xlim=c(0,2*pi),xlab="Wind direction",ylab="pm10",xaxt="n",ylim=c(0,50))
axis(1,c(0,pi/2,pi,3*pi/2,2*pi),c(expression(0),expression(pi/2),expression(pi),expression(3*pi/2),expression(2*pi)))
points(as.numeric(xseq),exp(modREFINED),type="l",col=2,lwd=3)

int_sup<-modREFINED-sesgok_final+qnorm(1-0.05/2)*sqrt(vark_final)
int_inf<-modREFINED-sesgok_final-qnorm(1-0.05/2)*sqrt(vark_final)

points(as.numeric(xseq),exp(int_sup),type="l",col="lightblue",lwd=3)
points(as.numeric(xseq),exp(int_inf),type="l",col="lightblue",lwd=3)

for (i in 1:nt){
	segments(xseq[i],exp(int_inf[i]),xseq[i],exp(int_sup[i]),lwd=2,col="lightblue")
}

points(as.numeric(xseq),exp(modREFINED),type="l",col=2,lwd=3)
abline(v=conversion.circular(circular(243,units="degrees"),units="radians"),col="grey30",lwd=2,lty=2)  # location of the factory



######################################################
# 	RADIAL REPRESENTATION (Figure 3: left)
######################################################

radial.plot(as.numeric(Y),X1,rp.type = "s",radial.lim=c(-20,115),point.symbols = 1,cex=1,labels=c("N","NE","E","SE","S","SW","W","NW"),label.pos=c(0,pi/4,pi/2,3*pi/4,pi,5*pi/4,3*pi/2,7*pi/4))
radial.plot(exp(modREFINED),circular(xseq),rp.type="l",add=TRUE,lwd=3,line.col=2,radial.lim=c(-20,115))

# location of the factory
radial.plot(120,rep(conversion.circular(circular(243,units="degrees",modulo="2pi",rotation="clock"),units="radians",modulo="2pi",rotation="clock"),50),rp.type="s",add=T,point.symbols=8,lty=2,point.col=4,radial.lim=c(-20,115))
radial.plot(120,rep(conversion.circular(circular(243,units="degrees",modulo="2pi",rotation="clock"),units="radians",modulo="2pi",rotation="clock"),50),rp.type="s",add=T,point.symbols=16,lty=2,point.col=4,radial.lim=c(-20,115))
text(x=-82,y = -130, labels = "FACTORY", font=6,cex=0.75)	 # location of the factory



######################################################
# 	PARTIALLY LINEAR MODEL
######################################################


# Pm10 ~ alpha0 + alpha1*WSpeed + m(WDir)


startv<-c(log(mean(Y)),0,0,0)	# starting values for the estimation


mod_lin<-gamma_1(X2,Y,startv)	# estimation of the linear part
alpha0<-mod_lin[1]
alpha1<-mod_lin[2]

# Refined procedure for the selection of the smoothing parameter (TAKES A LONG TIME!, FINAL RESULT IS kappa_opt=5.711)
lowerk<-0.3
upperk<-15
lowerk_ast<-0.3
upperk_ast<-5
xseq<-seq(0,2*pi,length=100)
kappa_ast_p3<-optimize(function(kappa)gamma_CRSC_int(xseq,X1,x0=X2,Y,kappa,alpha0,alpha1,startv),lower=lowerk_ast,upper=upperk_ast)$minimum
kappa_ast_p3<-kappa_ast_p3*(192/167)^(2/9)
# Selection of optimal kapppa 
kappa_opt<-optimize(function(kappa)MSE_int(xseq,X1,x0=X2,Y,kappa,kappa_ast_p3,alpha0,alpha1,startv),lower=lowerk,upper=upperk)$minimum



# Estimation

xseq<-seq(0,2*pi,length=500)
xseq_theta<-c(xseq,X1)
k<-length(xseq)
L2<-gamma_2(xseq_theta,X1,X2,Y,kappa_opt,alpha0,alpha1,startv)
L3<-gamma_3(X2,Y,(L2[[1]])[(k+1):(length(xseq_theta))],startv)
L2<-L2[[1]]


# Representation

xi<-seq(1,max(X2),length=200)
final<-matrix(nrow=length(xseq),ncol=200)
for (i in 1:200){

	x0<-xi[i]
	final[,i]<-exp(L3[1]+L3[2]*rep(x0,length(xseq))+L2[1:length(xseq)])	# estimation in a grid of values of x and theta

}

final_f<-as.numeric(final)
f1<-min(final_f)
f2<-max(final_f)
f2-f1

micol<-colorRampPalette(c('lightgrey', 'indianred3'))(15)	# color scale

radial.plot(1,1,rp.type = "s",radial.labels=c("",0,5,10,15,20,25),radial.lim=c(-2,22.5),labels=c("N","NE","E","SE","S","SW","W","NW"),label.pos=c(0,pi/4,pi/2,3*pi/4,pi,5*pi/4,3*pi/2,7*pi/4),point.symbols=16,point.col="white",grid.col="black")


xi<-seq(1,max(X2),length=200)
for (i in 1:200){
	for(j in 1:length(xseq)){
		x0<-xi[i]
		estimacion<-exp(L3[1]+L3[2]*x0+L2[j])
		est_red<-round(estimacion,0)
		radial.plot(x0,xseq[j],rp.type="s",add=T,point.symbols=16,lty=1,point.col=micol[est_red-12],radial.lim=c(-2,22.5),grid.col="black")
	}
}


radial.plot.labels(lengths=seq(0,25,length=6),radial.pos=rep(0,6),units="radians",radial.lim=c(-2,25),
  start=0,clockwise=FALSE,c("0","5","10","15","20","25"),adj=NULL,pos=NULL,boxed.labels=FALSE)


# Legend
legend_image <- as.raster(matrix(micol, ncol=1))
rasterImage(legend_image, 27, 27, 29,10)

# Factory's location
text(x=31, y = seq(11,27,l=5), labels = seq(13,27,length=5),cex=0.75)
text(x=28,y = 29, labels = "pm10", font=6)
ence<-conversion.circular(circular(243,units="degrees",modulo="2pi",rotation="clock"),units="radians",modulo="2pi",rotation="clock")
radial.plot(25,as.numeric(ence),rp.type="s",add=T,point.symbols=8,lty=1,point.col=4,radial.lim=c(-2,22.5),grid.col="black",cex=1)
radial.plot(25,as.numeric(ence),rp.type="s",add=T,point.symbols=16,lty=1,point.col=4,radial.lim=c(-2,22.5),grid.col="black",cex=1)
text(x=-16,y = -25, labels = "FACTORY", font=6,cex=0.75)





