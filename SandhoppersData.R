## Supplementary R code for the paper 
## "A general framework for circular local likelihood regression" 

# Workflow for obtaining center panels of Figure 3
# and statistical analysis related to the sandhoppers data

################
library(plotrix)
library(HDiR)
################

###########################################
source("Functions_implementing_methods.R")
###########################################


data(sandhoppers)

X<-sandhoppers$angle[sandhoppers$species!="ND"]
Y<-sandhoppers$species[sandhoppers$species!="ND"]
Y<-as.numeric(gsub("salt", 0, gsub("brito", 1, Y)))

nt<-500
xseq<-seq(0,2*pi,length=nt)

startv<-c(log(mean(Y)/(1-mean(Y))),0,0,0)

lowerk<-0.3
upperk<-50

lowerk_ast<-0.3
upperk_ast<-10

kappa_ast_seq<-seq(lowerk_ast,upperk_ast,length=30)
kappa_seq<-seq(lowerk,upperk,length=100)

# Selection of the smoothing parameter

# Refined procedure (p=3-->p=1)  (takes a long time! it gives kappa_opt=7.83)
k_ast_p3_res<-sapply(1:length(kappa_ast_seq),function(i){logistic_CRSC_int(xseq,X,Y,kappa_ast_seq[i],startv)})
kappa_ast_p3<-kappa_ast_seq[which.min(k_ast_p3_res)]		
kappa_ast_p3<-kappa_ast_p3*(192/167)^(2/9)
k_res<-sapply(1:length(kappa_seq),function(i){logistic_MSE_int(xseq,X,Y,kappa_seq[i],kappa_ast_p3,startv)})
kappa_opt<-kappa_seq[which.min(k_res)]	



# Estimation
modREFINED<-logistic_p1(xseq,X,Y,kappa_opt,startv)[[1]]


# Bias and variance
sesgo<-bias_logistic_p3(xseq,X,Y,kappa=kappa_opt,kappa_ast=kappa_ast_p3,startv)
varianza<-variance_logistic_p1(xseq,X,Y,kappa=kappa_opt,startv)


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
# 	PLANAR REPRESENTATION (Figure 3: bottom middle)
######################################################

plot(as.numeric(X),Y,xlab=expression(Theta),ylab="Probability of Talorchestia brito",xaxt="n",pch=16)
axis(1,c(0,pi/2,pi,3*pi/2,2*pi),c(expression(0),expression(pi/2),expression(pi),expression(3*pi/2),expression(2*pi)))
points(xseq,exp(modREFINED)/(1+exp(modREFINED)),type="l",col=2,lwd=3)


int_sup<-modREFINED-sesgok_final+qnorm(1-0.05/2)*sqrt(vark_final)
int_inf<-modREFINED-sesgok_final-qnorm(1-0.05/2)*sqrt(vark_final)

points(xseq,exp(int_sup)/(1+exp(int_sup)),type="l",col="lightblue",lwd=3)
points(xseq,exp(int_inf)/(1+exp(int_inf)),type="l",col="lightblue",lwd=3)


for (i in 1:nt){
	segments(xseq[i],exp(int_sup[i])/(1+exp(int_sup[i])),xseq[i],exp(int_inf[i])/(1+exp(int_inf[i])),lwd=2,col="lightblue")
}

points(xseq,exp(modREFINED)/(1+exp(modREFINED)),type="l",col=2,lwd=3)


######################################################
# 	RADIAL REPRESENTATION (Figure 3: top middle)
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



radial.plot(0.1,0.1,rp.type = "s",point.col = "white",point.symbols = 1,cex=0.1,radial.lim=c(0,1),labels=c(expression(0),expression(pi/4),expression(pi/2),expression(3*pi/4),expression(pi),expression(5*pi/4),expression(3*pi/2),expression(7*pi/4)),label.pos=c(0,pi/4,pi/2,3*pi/4,pi,5*pi/4,3*pi/2,7*pi/4))
radial.plot(exp(modREFINED)/(1+exp(modREFINED)),circular(xseq),rp.type="l",add=TRUE,lwd=3,line.col=2,radial.lim=c(0,1))

radial.plot(exp(int_sup)/(1+exp(int_sup)),circular(xseq),rp.type="l",lwd=2,line.col="lightblue",add=TRUE,radial.lim=c(0,1))
radial.plot(exp(int_inf)/(1+exp(int_inf)),circular(xseq),rp.type="l",lwd=2,line.col="lightblue",add=TRUE,radial.lim=c(0,1))


for (i in 1:nt){
	xpos0 <- cos(xseq[i]) * exp(int_sup[i])/(1+exp(int_sup[i]))
	ypos0 <- sin(xseq[i]) * exp(int_sup[i])/(1+exp(int_sup[i]))
	xpos1 <- cos(xseq[i]) * exp(int_inf[i])/(1+exp(int_inf[i]))
	ypos1 <- sin(xseq[i]) * exp(int_inf[i] )/(1+exp(int_inf[i]))

	segments(xpos0+00*cos(xseq[i]),ypos0+00*sin(xseq[i]),xpos1+00*cos(xseq[i]),ypos1+00*sin(xseq[i]),lwd=4,col="lightblue")
}
radial.plot(exp(modREFINED)/(1+exp(modREFINED)),circular(xseq),rp.type="l",add=TRUE,lwd=3,line.col=2,radial.lim=c(0,1))



