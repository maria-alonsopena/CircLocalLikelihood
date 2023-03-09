## Supplementary R code for the paper 
## "A general framework for circular local likelihood regression" 

# Workflow for obtaining left panels of Figure 3
# and statistical analysis related to the Human Resonance Data

################
library(plotrix)
################

###########################################
source("Functions_implementing_methods.R")
###########################################


dat<-read.table("HumanMotorResonance.txt")

X<-dat$Angular.position
Y<-dat$Reflex.Amplitude

nt<-300
xseq<-seq(0,2*pi,length=nt)

lowerk<-0.3  # Minimum value possible for kappa
upperk<-100  # Maximum value possible for kappa

lowerk_ast<-0.3   	# Minimum value possible for kappa star
upperk_ast<-20	 	# Maximum value possible for kappa star


# Selection of the smoothing parameter

kappa_ast_p3<-optimize(function(kappa)norm_CRSC_int(xseq,X,Y,kappa),lower=lowerk_ast,upper=upperk_ast)$minimum
kappa_ast_p3<-kappa_ast_p3*(192/167)^(2/9)
kappa_opt<-optimize(function(kappa)norm_MSE_int(xseq,X,Y,kappa,kappa_ast_p3),lower=lowerk,upper=upperk)$minimum


# Estimation
modREFINED<-(norm_p1(xseq,X,Y,kappa_opt))[[1]]


# Bias and variance
sesgo<-bias_norm_p3(xseq,X,Y,kappa=kappa_opt,kappa_ast=kappa_ast_p3)
varianza<-variance_norm_p3(xseq,X,Y,kappa=kappa_opt,kappa_ast=kappa_ast_p3)


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
# 	PLANAR REPRESENTATION (Figure 3: bottom left)
######################################################

plot(X,Y,pch=16,xlab=expression(Theta),ylab="Amplitude",xaxt="n" )
axis(1,c(0,pi/2,pi,3*pi/2,2*pi),c(expression(0),expression(pi/2),expression(pi),expression(3*pi/2),expression(2*pi)))
points(xseq,modREFINED,type="l",col=2,lwd=3)

int_sup<-modREFINED-sesgok_final+qnorm(1-0.05/2)*sqrt(vark_final)
int_inf<-modREFINED-sesgok_final-qnorm(1-0.05/2)*sqrt(vark_final)

points(xseq,int_sup,type="l",col="lightblue",lwd=3)
points(xseq,int_inf,type="l",col="lightblue",lwd=3)

for (i in 1:nt){
	segments(xseq[i],int_inf[i],xseq[i],int_sup[i],lwd=2,col="lightblue")
}

points(xseq,modREFINED,type="l",col=2,lwd=3)
points(X,Y,pch=16,col="grey50")



######################################################
# 	RADIAL REPRESENTATION (Figure 3: top left)
######################################################

radial.plot(as.numeric(Y),X,rp.type = "s",radial.lim=c(-150,152),point.symbols = 16,cex=1,labels=c(expression(0),expression(pi/4),expression(pi/2),expression(3*pi/4),expression(pi),expression(5*pi/4),expression(3*pi/2),expression(7*pi/4)),label.pos=c(0,pi/4,pi/2,3*pi/4,pi,5*pi/4,3*pi/2,7*pi/4))
radial.plot(modREFINED,circular(xseq),rp.type="l",add=TRUE,lwd=3,line.col=2,radial.lim=c(-150,152))

radial.plot(int_sup,circular(xseq),rp.type="l",lwd=2,line.col="lightblue",add=TRUE,radial.lim=c(-150,152))
radial.plot(int_inf,circular(xseq),rp.type="l",lwd=2,line.col="lightblue",add=TRUE,radial.lim=c(-150,152))


for (i in 1:nt){
	xpos0 <- cos(xseq[i]) * int_sup[i]
	ypos0 <- sin(xseq[i]) * int_sup[i]
	xpos1 <- cos(xseq[i]) * int_inf[i]
	ypos1 <- sin(xseq[i]) * int_inf[i] 

	segments(xpos0+150*cos(xseq[i]),ypos0+150*sin(xseq[i]),xpos1+150*cos(xseq[i]),ypos1+150*sin(xseq[i]),lwd=4,col="lightblue")
}
radial.plot(modREFINED,circular(xseq),rp.type="l",add=TRUE,lwd=3,line.col=2,radial.lim=c(-150,152))
radial.plot(Y,X,add=TRUE,rp.type = "s",point.symbols = 16,cex=1,radial.lim=c(-150,152),labels=c(expression(0),expression(pi/4),expression(pi/2),expression(3*pi/4),expression(pi),expression(5*pi/4),expression(3*pi/2),expression(7*pi/4)))


