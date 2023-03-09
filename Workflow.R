## Supplementary R code for the paper 
## "A general framework for circular local likelihood regression" 

# Workflow for obtaining (part of) Table 2 in the paper.
# This workflow leads to the results in Table 2 for 
# sample sizes of n=70. 
# For the remaining sample sizes, one should just change the values 
# of n and the maximum value for the search of the smoothing parameter.

# Running times are very large, except for the normal models.
# We do not recommend running code for models B1, B2, P1, P2, G1, G2,
# since they involve an iterative estimation.
# Results for models N1 and N2 can be otained in reasonable times.



#################
library(NPCirc)
library(Bolstad2)
#################

##########################################
source("Functions_implementing_methods.R")
##########################################




###############################################################
###############################################################

## SIMULATION STUDY


#### Model N1 ####

n<-70  # Sample size
xseq<-seq(0,2*pi,length=100)	# grid where the regression funciton will be evaluated

lowerk<-0.3  # Minimum value possible for kappa
upperk<-100  # Maximum value possible for kappa

lowerk_ast<-0.3   	# Minimum value possible for kappa star
upperk_ast<-20	 	# Maximum value possible for kappa star


m<-sin(2*xseq)*cos(xseq)    		# Regression function evaluated at xseq
m2_int<-sintegral(xseq,m^2)$int	# Integral of m^2

B<-500  # Number of replicates
resultsPILOT<-numeric(B)     	# The approx. ISE with the pilot method will be stored here
resultsREFINED<-numeric(B)		# The approx. ISE with the refined method will be stored here
resultsCV<-numeric(B)			# The approx. ISE with the cross-validation method will be stored here
kappaPILOT<-numeric(B)			# The selected kappas with the pilot method will be stored here
kappaREFINED<-numeric(B)		# The selected kappas with the refined method will be stored here
kappaCV<-numeric(B)			# The selected kappas with the cross-validation method will be stored here

for(b in 1:B){
  
  set.seed(888*b)					# Seed for reproducing the same results
  x<-runif(n,0,2*pi)				# Simulating values of the covariate
  y<-sin(2*x)*cos(x)+rnorm(n,sd=0.35)	# simulating values of the response)
  
  # Selection of the pilot kappa (p=1)
  kappa_ast<-optimize(function(kappa)norm_CRSC_p1_int(xseq,x,y,kappa),lower=lowerk,upper=upperk)$minimum
  kappaPILOT[b]<-kappa_ast
  
  
  # Selection with the refined procedure (pilot with p=3, estimation with p=1)
  kappa_ast_p3<-optimize(function(kappa)norm_CRSC_int(xseq,x,y,kappa),lower=lowerk_ast,upper=upperk_ast)$minimum
  kappa_ast_p3<-kappa_ast_p3*(192/167)^(2/9)
  kappa_opt<-optimize(function(kappa)norm_MSE_int(xseq,x,y,kappa,kappa_ast_p3),lower=lowerk,upper=upperk)$minimum
  kappaREFINED[b]<-kappa_opt
  
  # Selection by cross-validation
  kappa_cv<-optimize(function(kappa)lscv(x,y,kappa),lower=lowerk,upper=upperk)$minimum
  kappaCV[b]<-kappa_cv
  
  
 
  modPILOT<-(norm_p1(xseq,x,y,kappa_ast))[[1]] 		# Estimation with pilot kappa
  modREFINED<-(norm_p1(xseq,x,y,kappa_opt))[[1]]		# Estimation with refined kappa
  modCV<-(norm_p1(xseq,x,y,kappa_cv))[[1]]			# Estimation with cross-validation kappa
  
  aaPILOT<-(modPILOT-m)^2					# L2 distance between true and estimated functions (pilot)
  aaREFINED<-(modREFINED-m)^2				# L2 distance between true and estimated functions (refined)
  aaCV<-(modCV-m)^2						# L2 distance between true and estimated functions (cross-validation)
  
  resultsPILOT[b]<-sintegral(xseq,aaPILOT)$int/m2_int		# Approx. ISE (pilot)
  resultsREFINED[b]<-sintegral(xseq,aaREFINED)$int/m2_int		# Approx. ISE (refined)
  resultsCV[b]<-sintegral(xseq,aaCV)$int/m2_int			# Approx. ISE (cross-validation)
  
  print(c(b,resultsPILOT[b],resultsREFINED[b],resultsCV[b]))

}

kappas_N1<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		# Concentration parameters selected with each method  on each data replication
ISE_N1<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication

colMeans(ISE_N1)		# Monte Carlo averages of the approx. ISE (Table 2, row 1, columns 1 to 3)


#### Model N2 ####

n<-70  # Sample size
xseq<-seq(0,2*pi,length=100)	# grid where the regression funciton will be evaluated

lowerk<-0.3  # Minimum value possible for kappa
upperk<-100  # Maximum value possible for kappa

lowerk_ast<-0.3   	# Minimum value possible for kappa star
upperk_ast<-20		# Maximum value possible for kappa star


m<-1.75*cos(1*xseq-pi)*sin(1*xseq)+cos(xseq)   	# Regression function evaluated at xseq
m2_int<-sintegral(xseq,m^2)$int				# Integral of m^2

B<-500  # Number of replicates
resultsPILOT<-numeric(B)     	# The approx. ISE with the pilot method will be stored here
resultsREFINED<-numeric(B)		# The approx. ISE with the refined method will be stored here
resultsCV<-numeric(B)			# The approx. ISE with the cross-validation method will be stored here
kappaPILOT<-numeric(B)			# The selected kappas with the pilot method will be stored here
kappaREFINED<-numeric(B)		# The selected kappas with the refined method will be stored here
kappaCV<-numeric(B)			# The selected kappas with the cross-validation method will be stored here

for(b in 1:B){
  
  set.seed(8887*b)
  x<-runif(n,0,2*pi)
  y<-1.75*cos(1*x-pi)*sin(1*x)+cos(x)+rnorm(n,sd=0.5)
  
  # Selection of the pilot kappa (p=1)
  kappa_ast<-optimize(function(kappa)norm_CRSC_p1_int(xseq,x,y,kappa),lower=lowerk,upper=upperk)$minimum
  kappaPILOT[b]<-kappa_ast
  
  
  # Selection with the refined procedure (pilot with p=3, estimation with p=1)
  kappa_ast_p3<-optimize(function(kappa)norm_CRSC_int(xseq,x,y,kappa),lower=lowerk_ast,upper=upperk_ast)$minimum
  kappa_ast_p3<-kappa_ast_p3*(192/167)^(2/9)
  kappa_opt<-optimize(function(kappa)norm_MSE_int(xseq,x,y,kappa,kappa_ast_p3),lower=lowerk,upper=upperk)$minimum
  kappaREFINED[b]<-kappa_opt
  
  # Selection by cross-validation
  kappa_cv<-optimize(function(kappa)lscv(x,y,kappa),lower=lowerk,upper=upperk)$minimum
  kappaCV[b]<-kappa_cv
  
 
  modPILOT<-(norm_p1(xseq,x,y,kappa_ast))[[1]] 		# Estimation with pilot kappa
  modREFINED<-(norm_p1(xseq,x,y,kappa_opt))[[1]]		# Estimation with refined kappa
  modCV<-(norm_p1(xseq,x,y,kappa_cv))[[1]]			# Estimation with cross-validation kappa
  
  aaPILOT<-(modPILOT-m)^2					# L2 distance between true and estimated functions (pilot)
  aaREFINED<-(modREFINED-m)^2				# L2 distance between true and estimated functions (refined)
  aaCV<-(modCV-m)^2						# L2 distance between true and estimated functions (cross-validation)
  
  resultsPILOT[b]<-sintegral(xseq,aaPILOT)$int/m2_int		# Approx. ISE (pilot)
  resultsREFINED[b]<-sintegral(xseq,aaREFINED)$int/m2_int		# Approx. ISE (refined)
  resultsCV[b]<-sintegral(xseq,aaCV)$int/m2_int			# Approx. ISE (cross-validation)
  
  print(c(b,resultsPILOT[b],resultsREFINED[b],resultsCV[b]))

}

kappas_N2<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		# Concentration parameters selected with each method  on each data replication
ISE_N2<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication

colMeans(ISE_N2)		# Monte Carlo averages of the approx. ISE (Table 2, row 1, columns 4 to 6)



#### Model B1 ####


n<-70  # Sample size
xseq<-seq(0,2*pi,length=100)	# grid where the regression funciton will be evaluated

lowerk<-0.3  # Minimum value possible for kappa
upperk<-45  # Maximum value possible for kappa

lowerk_ast<-0.3   # Minimum value possible for kappa star
upperk_ast<-10	# Maximum value possible for kappa star

# In the logistic case, for some values of kappa, the estimator may not exist, so we perform the minimization 
# over a grid of valeus of kappa, in order to be able to discard those kappas for which the estimator 
# does not exist

kappa_ast_seq<-seq(lowerk_ast,upperk_ast,length=50)   	# Grid for kappa star
kappa_seq<-seq(lowerk,upperk,length=200)			 	# Grid for kappa	

m1<-2*sin(xseq)*cos(2*xseq)
m<-exp(m1)/(1+exp(m1))					# Regression function evaluated at xseq
m2_int<-sintegral(xseq,m^2)$int			# Integral of m^2


B<-500  # Number of replicates
resultsPILOT<-numeric(B)     	# The approx. ISE with the pilot method will be stored here
resultsREFINED<-numeric(B)		# The approx. ISE with the refined method will be stored here
resultsCV<-numeric(B)			# The approx. ISE with the cross-validation method will be stored here
kappaPILOT<-numeric(B)			# The selected kappas with the pilot method will be stored here
kappaREFINED<-numeric(B)		# The selected kappas with the refined method will be stored here
kappaCV<-numeric(B)			# The selected kappas with the cross-validation method will be stored here

for(b in 1:B){
  
  set.seed(888*b)	
  x<-runif(n,0,2*pi)
  logit<-2*sin(x)*cos(2*x)
  y<-rbinom(n,size=1, prob=exp(logit)/(1+exp(logit)))

  startv<-c(log(mean(y)/(1-mean(y))),0,0,0)  # Starting values

  
  	# Pilot kappa (p=1)
	k_ast_res<-sapply(1:length(kappa_seq),function(i){logistic_CRSC_p1_int(xseq,x,y,kappa_seq[i],startv[1:2])})
	if(all(is.na(k_ast_res))){kappaPILOT[b]<-NA}else{
		kappa_ast<-kappa_seq[which.min(k_ast_res)]	
		kappaPILOT[b]<-kappa_ast
	}
	

	# Refined procedure (p=3-->p=1)
	k_ast_p3_res<-sapply(1:length(kappa_ast_seq),function(i){logistic_CRSC_int(xseq,x,y,kappa_ast_seq[i],startv)})
	if(all(is.na(k_ast_p3_res))){kappa_ast_p3<-NA}else{
		kappa_ast_p3<-kappa_ast_seq[which.min(k_ast_p3_res)]		
		kappa_ast_p3<-kappa_ast_p3*(192/167)^(2/9)
	}
	if(!is.na(kappa_ast_p3)){
		k_res<-sapply(1:length(kappa_seq),function(i){logistic_MSE_int(xseq,x,y,kappa_seq[i],kappa_ast_p3,startv)})
		if(all(is.na(k_res))){kappaREFINED[b]<-NA}else{
			kappa_opt<-kappa_seq[which.min(k_res)]	
			kappaREFINED[b]<-kappa_opt
		}
	}else{kappaREFINED[b]<-NA}


	# cross-validation
	k_cv<-sapply(1:length(kappa_seq),function(i){logistic_loglik_cv(x,y,kappa_seq[i],startv)})
	kappa_cv<-kappa_seq[which.max(k_cv)]	
	kappaCV[b]<-kappa_cv
  

 	if(!is.na(kappaPILOT[b])){
		modPILOT<-(logistic_p1(xseq,x,y,kappaPILOT[b],startv))[[1]]
		aaPILOT<-(exp(modPILOT)/(1+exp(modPILOT))-m)^2
		resultsPILOT[b]<-sintegral(xseq,aaPILOT)$int/m2_int
	}else{
		resultsPILOT[b]<-NA
	}

	if(!is.na(kappaREFINED[b])){
		modREFINED<-(logistic_p1(xseq,x,y,kappaREFINED[b],startv))[[1]]
		aaREFINED<-(exp(modREFINED)/(1+exp(modREFINED))-m)^2
		resultsREFINED[b]<-sintegral(xseq,aaREFINED)$int/m2_int
	}else{
		resultsREFINED[b]<-NA
	}

	modCV<-(logistic_p1(xseq,x,y,kappaCV[b],startv))[[1]]
	aaCV<-(exp(modCV)/(1+exp(modCV))-m)^2
	resultsCV[b]<-sintegral(xseq,aaCV)$int/m2_int


  print(c(b,resultsPILOT[b],resultsREFINED[b],resultsCV[b]))

}

kappas_B1<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		# Concentration parameters selected with each method  on each data replication
ISE_B1<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication

colMeans(ISE_B1)			# Monte Carlo averages of the approx. ISE (Table 2, row 6, columns 1 to 3)


#### Model B2 ####


n<-70  # Sample size
xseq<-seq(0,2*pi,length=100)	# grid where the regression funciton will be evaluated

lowerk<-0.3  # Minimum value possible for kappa
upperk<-45  # Maximum value possible for kappa

lowerk_ast<-0.3   # Minimum value possible for kappa star
upperk_ast<-10	# Maximum value possible for kappa star

# In the logistic case, for some values of kappa, the estimator may not exist, so we perform the minimization 
# over a grid of valeus of kappa, in order to be able to discard those kappas for which the estimator 
# does not exist

kappa_ast_seq<-seq(lowerk_ast,upperk_ast,length=50)   	# Grid for kappa star
kappa_seq<-seq(lowerk,upperk,length=200)			 	# Grid for kappa	

m1<-log(1.6+1.5*sin(xseq)+exp(cos(xseq))/10)
m<-exp(m1)/(1+exp(m1))					# Regression function evaluated at xseq
m2_int<-sintegral(xseq,m^2)$int			# Integral of m^2


B<-500  # Number of replicates
resultsPILOT<-numeric(B)     	# The approx. ISE with the pilot method will be stored here
resultsREFINED<-numeric(B)		# The approx. ISE with the refined method will be stored here
resultsCV<-numeric(B)			# The approx. ISE with the cross-validation method will be stored here
kappaPILOT<-numeric(B)			# The selected kappas with the pilot method will be stored here
kappaREFINED<-numeric(B)		# The selected kappas with the refined method will be stored here
kappaCV<-numeric(B)			# The selected kappas with the cross-validation method will be stored here

for(b in 1:B){
  
  set.seed(888*b)	
  x<-runif(n,0,2*pi)
  logit<-log(1.6+1.5*sin(x)+exp(cos(x))/10)
  y<-rbinom(n,size=1, prob=exp(logit)/(1+exp(logit)))

  startv<-c(log(mean(y)/(1-mean(y))),0,0,0)  # Starting values

  
  	# Pilot kappa (p=1)
	k_ast_res<-sapply(1:length(kappa_seq),function(i){logistic_CRSC_p1_int(xseq,x,y,kappa_seq[i],startv[1:2])})
	if(all(is.na(k_ast_res))){kappaPILOT[b]<-NA}else{
		kappa_ast<-kappa_seq[which.min(k_ast_res)]	
		kappaPILOT[b]<-kappa_ast
	}
	

	# Refined procedure (p=3-->p=1)
	k_ast_p3_res<-sapply(1:length(kappa_ast_seq),function(i){logistic_CRSC_int(xseq,x,y,kappa_ast_seq[i],startv)})
	if(all(is.na(k_ast_p3_res))){kappa_ast_p3<-NA}else{
		kappa_ast_p3<-kappa_ast_seq[which.min(k_ast_p3_res)]		
		kappa_ast_p3<-kappa_ast_p3*(192/167)^(2/9)
	}
	if(!is.na(kappa_ast_p3)){
		k_res<-sapply(1:length(kappa_seq),function(i){logistic_MSE_int(xseq,x,y,kappa_seq[i],kappa_ast_p3,startv)})
		if(all(is.na(k_res))){kappaREFINED[b]<-NA}else{
			kappa_opt<-kappa_seq[which.min(k_res)]	
			kappaREFINED[b]<-kappa_opt
		}
	}else{kappaREFINED[b]<-NA}


	# cross-validation
	k_cv<-sapply(1:length(kappa_seq),function(i){logistic_loglik_cv(x,y,kappa_seq[i],startv)})
	kappa_cv<-kappa_seq[which.max(k_cv)]	
	kappaCV[b]<-kappa_cv
  

 	if(!is.na(kappaPILOT[b])){
		modPILOT<-(logistic_p1(xseq,x,y,kappaPILOT[b],startv))[[1]]
		aaPILOT<-(exp(modPILOT)/(1+exp(modPILOT))-m)^2
		resultsPILOT[b]<-sintegral(xseq,aaPILOT)$int/m2_int
	}else{
		resultsPILOT[b]<-NA
	}

	if(!is.na(kappaREFINED[b])){
		modREFINED<-(logistic_p1(xseq,x,y,kappaREFINED[b],startv))[[1]]
		aaREFINED<-(exp(modREFINED)/(1+exp(modREFINED))-m)^2
		resultsREFINED[b]<-sintegral(xseq,aaREFINED)$int/m2_int
	}else{
		resultsREFINED[b]<-NA
	}

	modCV<-(logistic_p1(xseq,x,y,kappaCV[b],startv))[[1]]
	aaCV<-(exp(modCV)/(1+exp(modCV))-m)^2
	resultsCV[b]<-sintegral(xseq,aaCV)$int/m2_int


  print(c(b,resultsPILOT[b],resultsREFINED[b],resultsCV[b]))

}

kappas_B2<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		# Concentration parameters selected with each method  on each data replication
ISE_B2<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication


colMeans(ISE_B2)		# Monte Carlo averages of the approx. ISE (Table 2, row 6, columns 4 to 6)




#### Model P1 ####

n<-70  # Sample size
xseq<-seq(0,2*pi,length=100)	# grid where the regression funciton will be evaluated

lowerk<-0.3  # Minimum value possible for kappa
upperk<-50  # Maximum value possible for kappa

lowerk_ast<-0.3   	# Minimum value possible for kappa star
upperk_ast<-15	 	# Maximum value possible for kappa star


m<-log(5+exp(1.5*sin(2*xseq-3)))		# Regression function evaluated at xseq
m2_int<-sintegral(xseq,m^2)$int		# Integral of m^2

B<-500  # Number of replicates
resultsPILOT<-numeric(B)     	# The approx. ISE with the pilot method will be stored here
resultsREFINED<-numeric(B)		# The approx. ISE with the refined method will be stored here
resultsCV<-numeric(B)			# The approx. ISE with the cross-validation method will be stored here
kappaPILOT<-numeric(B)			# The selected kappas with the pilot method will be stored here
kappaREFINED<-numeric(B)		# The selected kappas with the refined method will be stored here
kappaCV<-numeric(B)			# The selected kappas with the cross-validation method will be stored here

for(b in 1:B){
  
	set.seed(888*b)
	x<-runif(n,0,2*pi)
	meanpois<-5+exp(1.5*sin(2*x-3))
	y<-rpois(n,lambda=meanpois)

	startv<-c(log(mean(y)),0,0,0) 	# Starting values

  
  # Selection of the pilot kappa (p=1)
  kappa_ast<-optimize(function(kappa)pois_CRSC_p1_int(xseq,x,y,kappa,startv[1:2]),lower=lowerk,upper=upperk)$minimum
  kappaPILOT[b]<-kappa_ast

 
  # Selection with the refined procedure (pilot with p=3, estimation with p=1)
  kappa_ast_p3<-optimize(function(kappa)pois_CRSC_int(xseq,x,y,kappa,startv),lower=lowerk_ast,upper=upperk_ast)$minimum
  kappa_ast_p3<-kappa_ast_p3*(192/167)^(2/9)
  kappa_opt<-optimize(function(kappa)pois_MSE_int(xseq,x,y,kappa,kappa_ast_p3,startv),lower=lowerk,upper=upperk)$minimum
  kappaREFINED[b]<-kappa_opt
  
  # Selection by cross-validation
  kappa_cv<-optimize(function(kappa)pois_loglik_cv(x,y,kappa,startv),lower=lowerk,upper=upperk,maximum=TRUE)$maximum
  kappaCV[b]<-kappa_cv


  
  
 
  modPILOT<-(pois_p1(xseq,x,y,kappa_ast,startv))[[1]] 		# Estimation with pilot kappa
  modREFINED<-(pois_p1(xseq,x,y,kappa_opt,startv))[[1]]		# Estimation with refined kappa
  modCV<-(pois_p1(xseq,x,y,kappa_cv,startv))[[1]]			# Estimation with cross-validation kappa
  
  aaPILOT<-(modPILOT-m)^2					# L2 distance between true and estimated functions (pilot)
  aaREFINED<-(modREFINED-m)^2				# L2 distance between true and estimated functions (refined)
  aaCV<-(modCV-m)^2						# L2 distance between true and estimated functions (cross-validation)
  
  resultsPILOT[b]<-sintegral(xseq,aaPILOT)$int/m2_int		# Approx. ISE (pilot)
  resultsREFINED[b]<-sintegral(xseq,aaREFINED)$int/m2_int		# Approx. ISE (refined)
  resultsCV[b]<-sintegral(xseq,aaCV)$int/m2_int			# Approx. ISE (cross-validation)
  
  print(c(b,resultsPILOT[b],resultsREFINED[b],resultsCV[b]))

}

kappas_P1<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		# Concentration parameters selected with each method  on each data replication
ISE_P1<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication

colMeans(ISE_P1)		# Monte Carlo averages of the approx. ISE (Table 2, row 11, columns 1 to 3)




#### Model P2 ####

n<-70  # Sample size
xseq<-seq(0,2*pi,length=100)	# grid where the regression funciton will be evaluated


lowerk<-0.3  # Minimum value possible for kappa
upperk<-200  # Maximum value possible for kappa

lowerk_ast<-0.3   	# Minimum value possible for kappa star
upperk_ast<-40	 	# Maximum value possible for kappa star


m<-log(40+20*sin(3*xseq))		# Regression function evaluated at xseq
m2_int<-sintegral(xseq,m^2)$int		# Integral of m^2

B<-500  # Number of replicates
resultsPILOT<-numeric(B)     	# The approx. ISE with the pilot method will be stored here
resultsREFINED<-numeric(B)		# The approx. ISE with the refined method will be stored here
resultsCV<-numeric(B)			# The approx. ISE with the cross-validation method will be stored here
kappaPILOT<-numeric(B)			# The selected kappas with the pilot method will be stored here
kappaREFINED<-numeric(B)		# The selected kappas with the refined method will be stored here
kappaCV<-numeric(B)			# The selected kappas with the cross-validation method will be stored here

for(b in 1:B){
  
	set.seed(8887*b)
	x<-runif(n,0,2*pi)
	meanpois<-40+20*sin(3*x)
	y<-rpois(n,lambda=meanpois)

	startv<-c(log(mean(y)),0,0,0) 	# Starting values

  
  # Selection of the pilot kappa (p=1)
  kappa_ast<-optimize(function(kappa)pois_CRSC_p1_int(xseq,x,y,kappa,startv[1:2]),lower=lowerk,upper=upperk)$minimum
  kappaPILOT[b]<-kappa_ast

 
  # Selection with the refined procedure (pilot with p=3, estimation with p=1)
  kappa_ast_p3<-optimize(function(kappa)pois_CRSC_int(xseq,x,y,kappa,startv),lower=lowerk_ast,upper=upperk_ast)$minimum
  kappa_ast_p3<-kappa_ast_p3*(192/167)^(2/9)
  kappa_opt<-optimize(function(kappa)pois_MSE_int(xseq,x,y,kappa,kappa_ast_p3,startv),lower=lowerk,upper=upperk)$minimum
  kappaREFINED[b]<-kappa_opt
  
  # Selection by cross-validation
  kappa_cv<-optimize(function(kappa)pois_loglik_cv(x,y,kappa,startv),lower=lowerk,upper=upperk,maximum=TRUE)$maximum
  kappaCV[b]<-kappa_cv

  
 
  modPILOT<-(pois_p1(xseq,x,y,kappa_ast,startv))[[1]] 		# Estimation with pilot kappa
  modREFINED<-(pois_p1(xseq,x,y,kappa_opt,startv))[[1]]		# Estimation with refined kappa
  modCV<-(pois_p1(xseq,x,y,kappa_cv,startv))[[1]]			# Estimation with cross-validation kappa
  
  aaPILOT<-(modPILOT-m)^2					# L2 distance between true and estimated functions (pilot)
  aaREFINED<-(modREFINED-m)^2				# L2 distance between true and estimated functions (refined)
  aaCV<-(modCV-m)^2						# L2 distance between true and estimated functions (cross-validation)
  
  resultsPILOT[b]<-sintegral(xseq,aaPILOT)$int/m2_int		# Approx. ISE (pilot)
  resultsREFINED[b]<-sintegral(xseq,aaREFINED)$int/m2_int		# Approx. ISE (refined)
  resultsCV[b]<-sintegral(xseq,aaCV)$int/m2_int			# Approx. ISE (cross-validation)
  
  print(c(b,resultsPILOT[b],resultsREFINED[b],resultsCV[b]))

}

kappas_P2<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		# Concentration parameters selected with each method  on each data replication
ISE_P2<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication

colMeans(ISE_P2)		# Monte Carlo averages of the approx. ISE (Table 2, row 11, columns 1 to 3)





#### Model G1 ####

n<-70  # Sample size
xseq<-seq(0,2*pi,length=100)	# grid where the regression funciton will be evaluated
nu<-1

lowerk<-0.3  # Minimum value possible for kappa
upperk<-40  # Maximum value possible for kappa

lowerk_ast<-0.3   	# Minimum value possible for kappa star
upperk_ast<-5	 	# Maximum value possible for kappa star


m<-log(4+4*sin(2*xseq)*cos(xseq)) 		# Regression function evaluated at xseq
m2_int<-sintegral(xseq,m^2)$int		# Integral of m^2

B<-500  # Number of replicates
resultsPILOT<-numeric(B)     	# The approx. ISE with the pilot method will be stored here
resultsREFINED<-numeric(B)		# The approx. ISE with the refined method will be stored here
resultsCV<-numeric(B)			# The approx. ISE with the cross-validation method will be stored here
kappaPILOT<-numeric(B)			# The selected kappas with the pilot method will be stored here
kappaREFINED<-numeric(B)		# The selected kappas with the refined method will be stored here
kappaCV<-numeric(B)			# The selected kappas with the cross-validation method will be stored here

for(b in 1:B){
  
  set.seed(888*b)
  x<-runif(n,0,2*pi)
  meangamma <- 4+4*sin(2*x)*cos(x)
  y<-rgamma(n,shape=nu,rate=nu/meangamma)

  startv<-c(log(mean(y)),0,0,0) 	# Starting values

  
  # Selection of the pilot kappa (p=1)
  kappa_ast<-optimize(function(kappa)gamma_CRSC_p1_int(xseq,x,y,kappa,startv[1:2]),lower=lowerk,upper=upperk)$minimum
  kappaPILOT[b]<-kappa_ast


 
  # Selection with the refined procedure (pilot with p=3, estimation with p=1)
  kappa_ast_p3<-optimize(function(kappa)gamma_CRSC_int(xseq,x,y,kappa,startv),lower=lowerk_ast,upper=upperk_ast)$minimum
  kappa_ast_p3<-kappa_ast_p3*(192/167)^(2/9)
  kappa_opt<-optimize(function(kappa)gamma_MSE_int(xseq,x,y,kappa,kappa_ast_p3,startv),lower=lowerk,upper=upperk)$minimum
  kappaREFINED[b]<-kappa_opt

  
  # Selection by cross-validation
  kappa_cv<-optimize(function(kappa)gamma_loglik_cv(x,y,kappa,startv),lower=lowerk,upper=upperk)$minimum
  kappaCV[b]<-kappa_cv

  
 
  modPILOT<-(gamma_p1(xseq,x,y,kappa_ast,startv))[[1]] 		# Estimation with pilot kappa
  modREFINED<-(gamma_p1(xseq,x,y,kappa_opt,startv))[[1]]		# Estimation with refined kappa
  modCV<-(gamma_p1(xseq,x,y,kappa_cv,startv))[[1]]			# Estimation with cross-validation kappa
  
  aaPILOT<-(modPILOT-m)^2					# L2 distance between true and estimated functions (pilot)
  aaREFINED<-(modREFINED-m)^2				# L2 distance between true and estimated functions (refined)
  aaCV<-(modCV-m)^2						# L2 distance between true and estimated functions (cross-validation)
  
  resultsPILOT[b]<-sintegral(xseq,aaPILOT)$int/m2_int		# Approx. ISE (pilot)
  resultsREFINED[b]<-sintegral(xseq,aaREFINED)$int/m2_int		# Approx. ISE (refined)
  resultsCV[b]<-sintegral(xseq,aaCV)$int/m2_int			# Approx. ISE (cross-validation)
  
  print(c(b,resultsPILOT[b],resultsREFINED[b],resultsCV[b]))

}

kappas_G1<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		# Concentration parameters selected with each method  on each data replication
ISE_G1<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication

colMeans(ISE_G1)		# Monte Carlo averages of the approx. ISE (Table 2, row 16, columns 1 to 3)





#### Model G2 ####

n<-70  # Sample size
xseq<-seq(0,2*pi,length=100)	# grid where the regression funciton will be evaluated
nu<-2

lowerk<-0.3  # Minimum value possible for kappa
upperk<-100  # Maximum value possible for kappa

lowerk_ast<-0.3   	# Minimum value possible for kappa star
upperk_ast<-12	 	# Maximum value possible for kappa star


m<-log(5+2*cos(3+1.5*sin(xseq))) 		# Regression function evaluated at xseq
m2_int<-sintegral(xseq,m^2)$int		# Integral of m^2

B<-500  # Number of replicates
resultsPILOT<-numeric(B)     	# The approx. ISE with the pilot method will be stored here
resultsREFINED<-numeric(B)		# The approx. ISE with the refined method will be stored here
resultsCV<-numeric(B)			# The approx. ISE with the cross-validation method will be stored here
kappaPILOT<-numeric(B)			# The selected kappas with the pilot method will be stored here
kappaREFINED<-numeric(B)		# The selected kappas with the refined method will be stored here
kappaCV<-numeric(B)			# The selected kappas with the cross-validation method will be stored here

for(b in 1:B){
  
  set.seed(888*b)
  x<-runif(n,0,2*pi)
  meangamma <- 5+2*cos(3+1.5*sin(x))
  y<-rgamma(n,shape=nu,rate=nu/meangamma)

  startv<-c(log(mean(y)),0,0,0) 	# Starting values

  
  # Selection of the pilot kappa (p=1)
  kappa_ast<-optimize(function(kappa)gamma_CRSC_p1_int(xseq,x,y,kappa,startv[1:2]),lower=lowerk,upper=upperk)$minimum
  kappaPILOT[b]<-kappa_ast

 
  # Selection with the refined procedure (pilot with p=3, estimation with p=1)
  kappa_ast_p3<-optimize(function(kappa)gamma_CRSC_int(xseq,x,y,kappa,startv),lower=lowerk_ast,upper=upperk_ast)$minimum
  kappa_ast_p3<-kappa_ast_p3*(192/167)^(2/9)
  kappa_opt<-optimize(function(kappa)gamma_MSE_int(xseq,x,y,kappa,kappa_ast_p3,startv),lower=lowerk,upper=upperk)$minimum
  kappaREFINED[b]<-kappa_opt

  
  # Selection by cross-validation
  kappa_cv<-optimize(function(kappa)gamma_loglik_cv(x,y,kappa,startv),lower=lowerk,upper=upperk)$minimum
  kappaCV[b]<-kappa_cv

  
  
 
  modPILOT<-(gamma_p1(xseq,x,y,kappa_ast,startv))[[1]] 		# Estimation with pilot kappa
  modREFINED<-(gamma_p1(xseq,x,y,kappa_opt,startv))[[1]]		# Estimation with refined kappa
  modCV<-(gamma_p1(xseq,x,y,kappa_cv,startv))[[1]]			# Estimation with cross-validation kappa
  
  aaPILOT<-(modPILOT-m)^2					# L2 distance between true and estimated functions (pilot)
  aaREFINED<-(modREFINED-m)^2				# L2 distance between true and estimated functions (refined)
  aaCV<-(modCV-m)^2						# L2 distance between true and estimated functions (cross-validation)
  
  resultsPILOT[b]<-sintegral(xseq,aaPILOT)$int/m2_int		# Approx. ISE (pilot)
  resultsREFINED[b]<-sintegral(xseq,aaREFINED)$int/m2_int		# Approx. ISE (refined)
  resultsCV[b]<-sintegral(xseq,aaCV)$int/m2_int			# Approx. ISE (cross-validation)
  
  print(c(b,resultsPILOT[b],resultsREFINED[b],resultsCV[b]))

}

kappas_G2<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		# Concentration parameters selected with each method  on each data replication
ISE_G2<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication

colMeans(ISE_G2)		# Monte Carlo averages of the approx. ISE (Table 2, row 16, columns 4 to 6)



