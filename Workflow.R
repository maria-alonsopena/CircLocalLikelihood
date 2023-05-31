## Supplementary R code for the paper 
## "A general framework for circular local likelihood regression" 

# Workflow for obtaining Table 2 and Figure 2 of the paper and the figures 
# S1-S16 in the Supplementary Material.


# This workflow consists of three parts:
# 1. Simulation study
# 2. Figure 2 of the paper
# 3. Figures in the Supplementary Material (S1-S16)


# PART 1: The first part of this workflow leads to the results in Table 2 for 
# all sample sizes: n=70, n=100, n=250, n=500 and n=1500. 

# IMPORTANT:
# Running times are very large, except for the normal models,
# because of the data-driven selection of the smoothing parameters
# with the three different methods.
# We do not recommend running code for Poisson, Bernoulli and 
# gamma models in a regular machine, since they involve an 
# iterative estimation.
# Results for models N1 and N2 can be obtained in reasonable times
# for small sample sizes, although running times may also be large
# for the largest sample sizes (n=500, n=1500).
# For this first part, the authors run the code in a supercomputer
# consisting of 357 nodes and powered by 714 Intel Xeon Ice Lake 
# 8352Y processors.


# PART 2: The second part of this workflow leads to the eight plots in
# Figure 2 of the paper regarding simulation results. 
# All of these plots refer to sample sizes n=250. Therefore, for 
# obtaining these plots it is only necessary to run the code of PART 1 
# for n=250.


# PART 3: The third part of this workflow leads to the figures in
# the supplementary material regarding simulation results (S1-S16). 
# In order to obtain all the panels in all the figures, one must run  
# the code of the whole simulation study (PART 1). Individual panels 
# and figures for a specific model/sample size can be obtained by 
# only running the code of PART 1 corresponding to such model and 
# sample size.


#################
library(NPCirc)
library(Bolstad2)
#################



# Load the code for all the methods

##########################################
source("Functions_implementing_methods.R")
##########################################




###############################################################
###############################################################

##        PART 1:          SIMULATION STUDY                  ##

###############################################################
###############################################################

##################
#### Model N1 ####
####  n = 70  ####  
##################

n<-70  # Sample size
xseq<-seq(0,2*pi,length=100)	# grid where the regression function will be evaluated

lowerk<-0.3  # Minimum value possible for kappa
upperk<-100  # Maximum value possible for kappa

lowerk_ast<-0.3   	# Minimum value possible for kappa star
upperk_ast<-20	 	  # Maximum value possible for kappa star


m<-sin(2*xseq)*cos(xseq)    		# Regression function evaluated at xseq
m2_int<-sintegral(xseq,m^2)$int	# Integral of m^2

B<-500  # Number of replicates
resultsPILOT<-numeric(B)     	# The approx. ISE with the pilot method will be stored here
resultsREFINED<-numeric(B)		# The approx. ISE with the refined method will be stored here
resultsCV<-numeric(B)			    # The approx. ISE with the cross-validation method will be stored here
kappaPILOT<-numeric(B)			  # The selected kappas with the pilot method will be stored here
kappaREFINED<-numeric(B)		  # The selected kappas with the refined method will be stored here
kappaCV<-numeric(B)		      	# The selected kappas with the cross-validation method will be stored here

for(b in 1:B){
  
  set.seed(888*b)					  # Seed for reproducing the same results
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
  
  aaPILOT<-(modPILOT-m)^2					  # L2 distance between true and estimated functions (pilot)
  aaREFINED<-(modREFINED-m)^2				# L2 distance between true and estimated functions (refined)
  aaCV<-(modCV-m)^2					      	# L2 distance between true and estimated functions (cross-validation)
  
  resultsPILOT[b]<-sintegral(xseq,aaPILOT)$int/m2_int		    # Approx. ISE (pilot)
  resultsREFINED[b]<-sintegral(xseq,aaREFINED)$int/m2_int		# Approx. ISE (refined)
  resultsCV[b]<-sintegral(xseq,aaCV)$int/m2_int		        	# Approx. ISE (cross-validation)
  
  print(paste("Replication", b, "out of", B))

}

kappas_N1_n70<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		  # Concentration parameters selected with each method  on each data replication
ISE_N1_n70<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication

MISE_N1_n70<-colMeans(ISE_N1_n70,na.rm=TRUE)		# Monte Carlo averages of the approx. ISE (Table 2, row 1, columns 1 to 3)

##################
#### Model N1 ####
####  n = 100 ####
##################

n<-100  # Sample size
xseq<-seq(0,2*pi,length=100)	# grid where the regression function will be evaluated

lowerk<-0.3  # Minimum value possible for kappa
upperk<-100  # Maximum value possible for kappa

lowerk_ast<-0.3   	# Minimum value possible for kappa star
upperk_ast<-20	 	  # Maximum value possible for kappa star


m<-sin(2*xseq)*cos(xseq)    		# Regression function evaluated at xseq
m2_int<-sintegral(xseq,m^2)$int	# Integral of m^2

B<-500  # Number of replicates
resultsPILOT<-numeric(B)     	# The approx. ISE with the pilot method will be stored here
resultsREFINED<-numeric(B)		# The approx. ISE with the refined method will be stored here
resultsCV<-numeric(B)			    # The approx. ISE with the cross-validation method will be stored here
kappaPILOT<-numeric(B)			  # The selected kappas with the pilot method will be stored here
kappaREFINED<-numeric(B)		  # The selected kappas with the refined method will be stored here
kappaCV<-numeric(B)		      	# The selected kappas with the cross-validation method will be stored here

for(b in 1:B){
  
  set.seed(888*b)					  # Seed for reproducing the same results
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
  
  aaPILOT<-(modPILOT-m)^2					  # L2 distance between true and estimated functions (pilot)
  aaREFINED<-(modREFINED-m)^2				# L2 distance between true and estimated functions (refined)
  aaCV<-(modCV-m)^2					      	# L2 distance between true and estimated functions (cross-validation)
  
  resultsPILOT[b]<-sintegral(xseq,aaPILOT)$int/m2_int		    # Approx. ISE (pilot)
  resultsREFINED[b]<-sintegral(xseq,aaREFINED)$int/m2_int		# Approx. ISE (refined)
  resultsCV[b]<-sintegral(xseq,aaCV)$int/m2_int		        	# Approx. ISE (cross-validation)
  
  print(paste("Replication", b, "out of", B))
  
}

kappas_N1_n100<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		  # Concentration parameters selected with each method  on each data replication
ISE_N1_n100<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication

MISE_N1_n100<-colMeans(ISE_N1_n100,na.rm=TRUE)		# Monte Carlo averages of the approx. ISE (Table 2, row 2, columns 1 to 3)


##################
#### Model N1 ####
#### n = 250  ####
##################

n<-250  # Sample size
xseq<-seq(0,2*pi,length=100)	# grid where the regression function will be evaluated

lowerk<-0.3  # Minimum value possible for kappa
upperk<-150  # Maximum value possible for kappa

lowerk_ast<-0.3   	# Minimum value possible for kappa star
upperk_ast<-30	 	  # Maximum value possible for kappa star


m<-sin(2*xseq)*cos(xseq)    		# Regression function evaluated at xseq
m2_int<-sintegral(xseq,m^2)$int	# Integral of m^2

B<-500  # Number of replicates
resultsPILOT<-numeric(B)     	# The approx. ISE with the pilot method will be stored here
resultsREFINED<-numeric(B)		# The approx. ISE with the refined method will be stored here
resultsCV<-numeric(B)			    # The approx. ISE with the cross-validation method will be stored here
kappaPILOT<-numeric(B)			  # The selected kappas with the pilot method will be stored here
kappaREFINED<-numeric(B)		  # The selected kappas with the refined method will be stored here
kappaCV<-numeric(B)		      	# The selected kappas with the cross-validation method will be stored here

for(b in 1:B){
  
  set.seed(888*b)					  # Seed for reproducing the same results
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
  
  aaPILOT<-(modPILOT-m)^2					  # L2 distance between true and estimated functions (pilot)
  aaREFINED<-(modREFINED-m)^2				# L2 distance between true and estimated functions (refined)
  aaCV<-(modCV-m)^2					      	# L2 distance between true and estimated functions (cross-validation)
  
  resultsPILOT[b]<-sintegral(xseq,aaPILOT)$int/m2_int		    # Approx. ISE (pilot)
  resultsREFINED[b]<-sintegral(xseq,aaREFINED)$int/m2_int		# Approx. ISE (refined)
  resultsCV[b]<-sintegral(xseq,aaCV)$int/m2_int		        	# Approx. ISE (cross-validation)
  
  print(paste("Replication", b, "out of", B))
  
}

kappas_N1_n250<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		  # Concentration parameters selected with each method  on each data replication
ISE_N1_n250<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication

MISE_N1_n250<-colMeans(ISE_N1_n250,na.rm=TRUE)		# Monte Carlo averages of the approx. ISE (Table 2, row 3, columns 1 to 3)



##################
#### Model N1 ####
#### n = 500  ####
##################

n<-500  # Sample size
xseq<-seq(0,2*pi,length=100)	# grid where the regression function will be evaluated

lowerk<-0.3  # Minimum value possible for kappa
upperk<-150  # Maximum value possible for kappa

lowerk_ast<-0.3   	# Minimum value possible for kappa star
upperk_ast<-40	 	  # Maximum value possible for kappa star


m<-sin(2*xseq)*cos(xseq)    		# Regression function evaluated at xseq
m2_int<-sintegral(xseq,m^2)$int	# Integral of m^2

B<-500  # Number of replicates
resultsPILOT<-numeric(B)     	# The approx. ISE with the pilot method will be stored here
resultsREFINED<-numeric(B)		# The approx. ISE with the refined method will be stored here
resultsCV<-numeric(B)			    # The approx. ISE with the cross-validation method will be stored here
kappaPILOT<-numeric(B)			  # The selected kappas with the pilot method will be stored here
kappaREFINED<-numeric(B)		  # The selected kappas with the refined method will be stored here
kappaCV<-numeric(B)		      	# The selected kappas with the cross-validation method will be stored here

for(b in 1:B){
  
  set.seed(888*b)					  # Seed for reproducing the same results
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
  
  aaPILOT<-(modPILOT-m)^2					  # L2 distance between true and estimated functions (pilot)
  aaREFINED<-(modREFINED-m)^2				# L2 distance between true and estimated functions (refined)
  aaCV<-(modCV-m)^2					      	# L2 distance between true and estimated functions (cross-validation)
  
  resultsPILOT[b]<-sintegral(xseq,aaPILOT)$int/m2_int		    # Approx. ISE (pilot)
  resultsREFINED[b]<-sintegral(xseq,aaREFINED)$int/m2_int		# Approx. ISE (refined)
  resultsCV[b]<-sintegral(xseq,aaCV)$int/m2_int		        	# Approx. ISE (cross-validation)
  
  print(paste("Replication", b, "out of", B))
  
}

kappas_N1_n500<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		  # Concentration parameters selected with each method  on each data replication
ISE_N1_n500<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication

MISE_N1_n500<-colMeans(ISE_N1_n500,na.rm=TRUE)		# Monte Carlo averages of the approx. ISE (Table 2, row 4, columns 1 to 3)


##################
#### Model N1 ####
#### n = 1500 ####
##################

n<-1500  # Sample size
xseq<-seq(0,2*pi,length=100)	# grid where the regression function will be evaluated

lowerk<-0.3  # Minimum value possible for kappa
upperk<-200  # Maximum value possible for kappa

lowerk_ast<-0.3   	# Minimum value possible for kappa star
upperk_ast<-40	 	  # Maximum value possible for kappa star


m<-sin(2*xseq)*cos(xseq)    		# Regression function evaluated at xseq
m2_int<-sintegral(xseq,m^2)$int	# Integral of m^2

B<-500  # Number of replicates
resultsPILOT<-numeric(B)     	# The approx. ISE with the pilot method will be stored here
resultsREFINED<-numeric(B)		# The approx. ISE with the refined method will be stored here
resultsCV<-numeric(B)			    # The approx. ISE with the cross-validation method will be stored here
kappaPILOT<-numeric(B)			  # The selected kappas with the pilot method will be stored here
kappaREFINED<-numeric(B)		  # The selected kappas with the refined method will be stored here
kappaCV<-numeric(B)		      	# The selected kappas with the cross-validation method will be stored here

for(b in 1:B){
  
  set.seed(888*b)					  # Seed for reproducing the same results
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
  
  aaPILOT<-(modPILOT-m)^2					  # L2 distance between true and estimated functions (pilot)
  aaREFINED<-(modREFINED-m)^2				# L2 distance between true and estimated functions (refined)
  aaCV<-(modCV-m)^2					      	# L2 distance between true and estimated functions (cross-validation)
  
  resultsPILOT[b]<-sintegral(xseq,aaPILOT)$int/m2_int		    # Approx. ISE (pilot)
  resultsREFINED[b]<-sintegral(xseq,aaREFINED)$int/m2_int		# Approx. ISE (refined)
  resultsCV[b]<-sintegral(xseq,aaCV)$int/m2_int		        	# Approx. ISE (cross-validation)
  
  print(paste("Replication", b, "out of", B))
  
}

kappas_N1_n1500<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		  # Concentration parameters selected with each method  on each data replication
ISE_N1_n1500<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication

MISE_N1_n1500<-colMeans(ISE_N1_n1500,na.rm=TRUE)		# Monte Carlo averages of the approx. ISE (Table 2, row 5, columns 1 to 3)


##################
#### Model N2 ####
####  n = 70  ####
##################

n<-70  # Sample size
xseq<-seq(0,2*pi,length=100)	# grid where the regression function will be evaluated

lowerk<-0.3  # Minimum value possible for kappa
upperk<-100  # Maximum value possible for kappa

lowerk_ast<-0.3   # Minimum value possible for kappa star
upperk_ast<-20		# Maximum value possible for kappa star


m<-1.75*cos(1*xseq-pi)*sin(1*xseq)+cos(xseq)   	# Regression function evaluated at xseq
m2_int<-sintegral(xseq,m^2)$int				# Integral of m^2

B<-500  # Number of replicates
resultsPILOT<-numeric(B)     	# The approx. ISE with the pilot method will be stored here
resultsREFINED<-numeric(B)		# The approx. ISE with the refined method will be stored here
resultsCV<-numeric(B)			    # The approx. ISE with the cross-validation method will be stored here
kappaPILOT<-numeric(B)			  # The selected kappas with the pilot method will be stored here
kappaREFINED<-numeric(B)		  # The selected kappas with the refined method will be stored here
kappaCV<-numeric(B)		      	# The selected kappas with the cross-validation method will be stored here

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
  aaREFINED<-(modREFINED-m)^2			# L2 distance between true and estimated functions (refined)
  aaCV<-(modCV-m)^2						    # L2 distance between true and estimated functions (cross-validation)
  
  resultsPILOT[b]<-sintegral(xseq,aaPILOT)$int/m2_int		    # Approx. ISE (pilot)
  resultsREFINED[b]<-sintegral(xseq,aaREFINED)$int/m2_int		# Approx. ISE (refined)
  resultsCV[b]<-sintegral(xseq,aaCV)$int/m2_int			        # Approx. ISE (cross-validation)
  
  print(paste("Replication", b, "out of", B))

}

kappas_N2_n70<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		  # Concentration parameters selected with each method  on each data replication
ISE_N2_n70<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication

MISE_N2_n70<-colMeans(ISE_N2_n70,na.rm=TRUE)		# Monte Carlo averages of the approx. ISE (Table 2, row 1, columns 4 to 6)


##################
#### Model N2 ####
####  n = 100 #### 
##################

n<-100  # Sample size
xseq<-seq(0,2*pi,length=100)	# grid where the regression function will be evaluated

lowerk<-0.3  # Minimum value possible for kappa
upperk<-150  # Maximum value possible for kappa

lowerk_ast<-0.3   # Minimum value possible for kappa star
upperk_ast<-20		# Maximum value possible for kappa star


m<-1.75*cos(1*xseq-pi)*sin(1*xseq)+cos(xseq)   	# Regression function evaluated at xseq
m2_int<-sintegral(xseq,m^2)$int				# Integral of m^2

B<-500  # Number of replicates
resultsPILOT<-numeric(B)     	# The approx. ISE with the pilot method will be stored here
resultsREFINED<-numeric(B)		# The approx. ISE with the refined method will be stored here
resultsCV<-numeric(B)			    # The approx. ISE with the cross-validation method will be stored here
kappaPILOT<-numeric(B)			  # The selected kappas with the pilot method will be stored here
kappaREFINED<-numeric(B)		  # The selected kappas with the refined method will be stored here
kappaCV<-numeric(B)		      	# The selected kappas with the cross-validation method will be stored here

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
  aaREFINED<-(modREFINED-m)^2			# L2 distance between true and estimated functions (refined)
  aaCV<-(modCV-m)^2						    # L2 distance between true and estimated functions (cross-validation)
  
  resultsPILOT[b]<-sintegral(xseq,aaPILOT)$int/m2_int		    # Approx. ISE (pilot)
  resultsREFINED[b]<-sintegral(xseq,aaREFINED)$int/m2_int		# Approx. ISE (refined)
  resultsCV[b]<-sintegral(xseq,aaCV)$int/m2_int			        # Approx. ISE (cross-validation)
  
  print(paste("Replication", b, "out of", B))
  
}

kappas_N2_n100<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		  # Concentration parameters selected with each method  on each data replication
ISE_N2_n100<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication

MISE_N2_n100<-colMeans(ISE_N2_n100,na.rm=TRUE)		# Monte Carlo averages of the approx. ISE (Table 2, row 2, columns 4 to 6)


##################
#### Model N2 ####
#### n = 250  ####
##################

n<-250  # Sample size
xseq<-seq(0,2*pi,length=100)	# grid where the regression function will be evaluated

lowerk<-0.3  # Minimum value possible for kappa
upperk<-200  # Maximum value possible for kappa

lowerk_ast<-0.3   # Minimum value possible for kappa star
upperk_ast<-20		# Maximum value possible for kappa star


m<-1.75*cos(1*xseq-pi)*sin(1*xseq)+cos(xseq)   	# Regression function evaluated at xseq
m2_int<-sintegral(xseq,m^2)$int				# Integral of m^2

B<-500  # Number of replicates
resultsPILOT<-numeric(B)     	# The approx. ISE with the pilot method will be stored here
resultsREFINED<-numeric(B)		# The approx. ISE with the refined method will be stored here
resultsCV<-numeric(B)			    # The approx. ISE with the cross-validation method will be stored here
kappaPILOT<-numeric(B)			  # The selected kappas with the pilot method will be stored here
kappaREFINED<-numeric(B)		  # The selected kappas with the refined method will be stored here
kappaCV<-numeric(B)		      	# The selected kappas with the cross-validation method will be stored here

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
  aaREFINED<-(modREFINED-m)^2			# L2 distance between true and estimated functions (refined)
  aaCV<-(modCV-m)^2						    # L2 distance between true and estimated functions (cross-validation)
  
  resultsPILOT[b]<-sintegral(xseq,aaPILOT)$int/m2_int		    # Approx. ISE (pilot)
  resultsREFINED[b]<-sintegral(xseq,aaREFINED)$int/m2_int		# Approx. ISE (refined)
  resultsCV[b]<-sintegral(xseq,aaCV)$int/m2_int			        # Approx. ISE (cross-validation)
  
  print(paste("Replication", b, "out of", B))
  
}

kappas_N2_n250<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		  # Concentration parameters selected with each method  on each data replication
ISE_N2_n250<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication

MISE_N2_n250<-colMeans(ISE_N2_n250,na.rm=TRUE)		# Monte Carlo averages of the approx. ISE (Table 2, row 3, columns 4 to 6)



##################
#### Model N2 ####
#### n = 500  ####
##################

n<-500  # Sample size
xseq<-seq(0,2*pi,length=100)	# grid where the regression function will be evaluated

lowerk<-0.3  # Minimum value possible for kappa
upperk<-200  # Maximum value possible for kappa

lowerk_ast<-0.3   # Minimum value possible for kappa star
upperk_ast<-20		# Maximum value possible for kappa star


m<-1.75*cos(1*xseq-pi)*sin(1*xseq)+cos(xseq)   	# Regression function evaluated at xseq
m2_int<-sintegral(xseq,m^2)$int				# Integral of m^2

B<-500  # Number of replicates
resultsPILOT<-numeric(B)     	# The approx. ISE with the pilot method will be stored here
resultsREFINED<-numeric(B)		# The approx. ISE with the refined method will be stored here
resultsCV<-numeric(B)			    # The approx. ISE with the cross-validation method will be stored here
kappaPILOT<-numeric(B)			  # The selected kappas with the pilot method will be stored here
kappaREFINED<-numeric(B)		  # The selected kappas with the refined method will be stored here
kappaCV<-numeric(B)		      	# The selected kappas with the cross-validation method will be stored here

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
  aaREFINED<-(modREFINED-m)^2			# L2 distance between true and estimated functions (refined)
  aaCV<-(modCV-m)^2						    # L2 distance between true and estimated functions (cross-validation)
  
  resultsPILOT[b]<-sintegral(xseq,aaPILOT)$int/m2_int		    # Approx. ISE (pilot)
  resultsREFINED[b]<-sintegral(xseq,aaREFINED)$int/m2_int		# Approx. ISE (refined)
  resultsCV[b]<-sintegral(xseq,aaCV)$int/m2_int			        # Approx. ISE (cross-validation)
  
  print(paste("Replication", b, "out of", B))
  
}

kappas_N2_n500<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		  # Concentration parameters selected with each method  on each data replication
ISE_N2_n500<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication

MISE_N2_n500<-colMeans(ISE_N2_n500,na.rm=TRUE)		# Monte Carlo averages of the approx. ISE (Table 2, row 4, columns 4 to 6)


##################
#### Model N2 ####
#### n = 1500 ####
##################

n<-1500  # Sample size
xseq<-seq(0,2*pi,length=100)	# grid where the regression function will be evaluated

lowerk<-0.3  # Minimum value possible for kappa
upperk<-200  # Maximum value possible for kappa

lowerk_ast<-0.3   # Minimum value possible for kappa star
upperk_ast<-40		# Maximum value possible for kappa star


m<-1.75*cos(1*xseq-pi)*sin(1*xseq)+cos(xseq)   	# Regression function evaluated at xseq
m2_int<-sintegral(xseq,m^2)$int				# Integral of m^2

B<-500  # Number of replicates
resultsPILOT<-numeric(B)     	# The approx. ISE with the pilot method will be stored here
resultsREFINED<-numeric(B)		# The approx. ISE with the refined method will be stored here
resultsCV<-numeric(B)			    # The approx. ISE with the cross-validation method will be stored here
kappaPILOT<-numeric(B)			  # The selected kappas with the pilot method will be stored here
kappaREFINED<-numeric(B)		  # The selected kappas with the refined method will be stored here
kappaCV<-numeric(B)		      	# The selected kappas with the cross-validation method will be stored here

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
  aaREFINED<-(modREFINED-m)^2			# L2 distance between true and estimated functions (refined)
  aaCV<-(modCV-m)^2						    # L2 distance between true and estimated functions (cross-validation)
  
  resultsPILOT[b]<-sintegral(xseq,aaPILOT)$int/m2_int		    # Approx. ISE (pilot)
  resultsREFINED[b]<-sintegral(xseq,aaREFINED)$int/m2_int		# Approx. ISE (refined)
  resultsCV[b]<-sintegral(xseq,aaCV)$int/m2_int			        # Approx. ISE (cross-validation)
  
  print(paste("Replication", b, "out of", B))
  
}

kappas_N2_n1500<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		  # Concentration parameters selected with each method  on each data replication
ISE_N2_n1500<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication

MISE_N2_n1500<-colMeans(ISE_N2_n1500,na.rm=TRUE)		# Monte Carlo averages of the approx. ISE (Table 2, row 5, columns 4 to 6)


##################
#### Model P1 ####
####  n = 70  ####
##################

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
  
  print(paste("Replication", b, "out of", B))

}

kappas_P1_n70<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		# Concentration parameters selected with each method  on each data replication
ISE_P1_n70<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication

MISE_P1_n70<-colMeans(ISE_P1_n70,na.rm=TRUE)		# Monte Carlo averages of the approx. ISE (Table 2, row 11, columns 1 to 3)


##################
#### Model P1 ####
#### n = 100  ####
##################

n<-100  # Sample size
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
  
  print(paste("Replication", b, "out of", B))
  
}

kappas_P1_n100<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		# Concentration parameters selected with each method  on each data replication
ISE_P1_n100<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication

MISE_P1_n100<-colMeans(ISE_P1_n100,na.rm=TRUE)		# Monte Carlo averages of the approx. ISE (Table 2, row 12, columns 1 to 3)


##################
#### Model P1 ####
####  n = 250 ####
##################

n<-250  # Sample size
xseq<-seq(0,2*pi,length=100)	# grid where the regression funciton will be evaluated

lowerk<-0.3  # Minimum value possible for kappa
upperk<-100  # Maximum value possible for kappa

lowerk_ast<-0.3   	# Minimum value possible for kappa star
upperk_ast<-25	 	# Maximum value possible for kappa star


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
  
  print(paste("Replication", b, "out of", B))
  
}

kappas_P1_n250<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		# Concentration parameters selected with each method  on each data replication
ISE_P1_n250<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication

MISE_P1_n250<-colMeans(ISE_P1_n250,na.rm=TRUE)		# Monte Carlo averages of the approx. ISE (Table 2, row 13, columns 1 to 3)


##################
#### Model P1 ####
#### n = 500  ####
##################

n<-500  # Sample size
xseq<-seq(0,2*pi,length=100)	# grid where the regression funciton will be evaluated

lowerk<-0.3  # Minimum value possible for kappa
upperk<-150  # Maximum value possible for kappa

lowerk_ast<-0.3   	# Minimum value possible for kappa star
upperk_ast<-45	 	# Maximum value possible for kappa star


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
  
  print(paste("Replication", b, "out of", B))
  
}

kappas_P1_n500<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		# Concentration parameters selected with each method  on each data replication
ISE_P1_n500<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication

MISE_P1_n500<-colMeans(ISE_P1_n500,na.rm=TRUE)		# Monte Carlo averages of the approx. ISE (Table 2, row 14, columns 1 to 3)



##################
#### Model P1 ####
#### n = 1500 ####
##################

n<-1500  # Sample size
xseq<-seq(0,2*pi,length=100)	# grid where the regression funciton will be evaluated

lowerk<-0.3  # Minimum value possible for kappa
upperk<-200  # Maximum value possible for kappa

lowerk_ast<-0.3   	# Minimum value possible for kappa star
upperk_ast<-55	 	# Maximum value possible for kappa star


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
  
  print(paste("Replication", b, "out of", B))
  
}

kappas_P1_n1500<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		# Concentration parameters selected with each method  on each data replication
ISE_P1_n1500<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication

MISE_P1_n1500<-colMeans(ISE_P1_n1500,na.rm=TRUE)		# Monte Carlo averages of the approx. ISE (Table 2, row 15, columns 1 to 3)


##################
#### Model P2 ####
####  n = 70  ####
##################

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
  
  print(paste("Replication", b, "out of", B))

}

kappas_P2_n70<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		# Concentration parameters selected with each method  on each data replication
ISE_P2_n70<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication

MISE_P2_n70<-colMeans(ISE_P2_n70,na.rm=TRUE)		# Monte Carlo averages of the approx. ISE (Table 2, row 11, columns 4 to 6)

##################
#### Model P2 ####
#### n = 100  ####
##################

n<-100  # Sample size
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
  
  print(paste("Replication", b, "out of", B))
  
}

kappas_P2_n100<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		# Concentration parameters selected with each method  on each data replication
ISE_P2_n100<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication

MISE_P2_n100<-colMeans(ISE_P2_n100,na.rm=TRUE)		# Monte Carlo averages of the approx. ISE (Table 2, row 12, columns 4 to 6)


##################
#### Model P2 ####
#### n = 250  ####
##################

n<-250  # Sample size
xseq<-seq(0,2*pi,length=100)	# grid where the regression funciton will be evaluated


lowerk<-0.3  # Minimum value possible for kappa
upperk<-200  # Maximum value possible for kappa

lowerk_ast<-0.3   	# Minimum value possible for kappa star
upperk_ast<-45	 	# Maximum value possible for kappa star


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
  
  print(paste("Replication", b, "out of", B))
  
}

kappas_P2_n250<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		# Concentration parameters selected with each method  on each data replication
ISE_P2_n250<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication

MISE_P2_n250<-colMeans(ISE_P2_n250,na.rm=TRUE)		# Monte Carlo averages of the approx. ISE (Table 2, row 13, columns 4 to 6)


##################
#### Model P2 ####
###  n = 500  ####
##################

n<-500  # Sample size
xseq<-seq(0,2*pi,length=100)	# grid where the regression funciton will be evaluated


lowerk<-0.3  # Minimum value possible for kappa
upperk<-200  # Maximum value possible for kappa

lowerk_ast<-0.3   	# Minimum value possible for kappa star
upperk_ast<-55	 	# Maximum value possible for kappa star


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
  
  print(paste("Replication", b, "out of", B))
  
}

kappas_P2_n500<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		# Concentration parameters selected with each method  on each data replication
ISE_P2_n500<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication

MISE_P2_n500<-colMeans(ISE_P2_n500,na.rm=TRUE)		# Monte Carlo averages of the approx. ISE (Table 2, row 14, columns 4 to 6)


##################
#### Model P2 ####
###  n = 1500 ####
##################

n<-1500  # Sample size
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
  
  
  kappa_ast_seq<-seq(0.5,40,length=100)
  kappa_seq<-seq(60,200,length=200)
  # Selection of the pilot kappa (p=1)
  k_ast_res<-sapply(1:length(kappa_seq),function(i){CRSC_p1_int(xseq,x,y,kappa_seq[i],startv[1:2])})
  if(all(is.na(k_ast_res))){kappaPILOT[b]<-NA}else{
    kappa_ast<-kappa_seq[which.min(k_ast_res)]	
    kappaPILOT[b]<-kappa_ast
  }
  if(!is.na(kappaPILOT[b])){
    modPILOT<-(pois_p1(xseq,x,y,kappa_ast,startv))[[1]]	
    aaPILOT<-(modPILOT-m)^2
    resultsPILOT[b]<-sintegral(xseq,aaPILOT)$int/m2_int
  }else{
    resultsPILOT[b]<-NA
  }
  # Selection with the refined procedure (pilot with p=3, estimation with p=1)
  k_ast_p3_res<-sapply(1:length(kappa_ast_seq),function(i){pois_CRSC_int(xseq,x,y,kappa_ast_seq[i],startv)})
  if(all(is.na(k_ast_p3_res))){kappa_ast_p3<-NA}else{
    kappa_ast_p3<-kappa_ast_seq[which.min(k_ast_p3_res)]		
    kappa_ast_p3<-kappa_ast_p3*(192/167)^(2/9)
  }
  if(!is.na(kappa_ast_p3)){
    k_res<-sapply(1:length(kappa_seq),function(i){MSE_int(xseq,x,y,kappa_seq[i],kappa_ast_p3,startv)})
    if(all(is.na(k_res))){kappaREFINED[b]<-NA}else{
      kappa_opt<-kappa_seq[which.min(k_res)]	
      kappaREFINED[b]<-kappa_opt
    }
  }else{kappaREFINED[b]<-NA}
  if(!is.na(kappaREFINED[b])){
    modREFINED<-(pois_p1(xseq,x,y,kappa_opt,startv))[[1]]
    aaREFINED<-(modREFINED-m)^2
    resultsREFINED[b]<-sintegral(xseq,aaREFINED)$int/m2_int
  }else{
    resultsREFINED[b]<-NA
  }
  # Selection by cross-validation
  k_cv<-sapply(1:length(kappa_seq),function(i){loglik_cv(x,y,kappa_seq[i],startv)})
  if(all(is.na(k_cv))){kappaCV[b]<-NA}else{
    kappa_cv<-kappa_seq[which.max(k_cv)]	
    kappaCV[b]<-kappa_cv
  }
  if(!is.na(kappaCV[b])){
    modCV<-(pois_p1(xseq,x,y,kappa_cv,startv))[[1]]
    aaCV<-(modCV-m)^2
    resultsCV[b]<-sintegral(xseq,aaCV)$int/m2_int
  }else{
    resultsCV[b]<-NA
  }
  print(paste("Replication", b, "out of", B))
  
}

kappas_P2_n1500<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		# Concentration parameters selected with each method  on each data replication
ISE_P2_n1500<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication

MISE_P2_n1500<-colMeans(ISE_P2_n1500,na.rm=TRUE)		# Monte Carlo averages of the approx. ISE (Table 2, row 15, columns 4 to 6)


##################
#### Model B1 ####
###   n = 70  ####
##################

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
  if(!is.na(kappaCV[b])){
    modCV<-(logistic_p1(xseq,x,y,kappaCV[b],startv))[[1]]
    aaCV<-(exp(modCV)/(1+exp(modCV))-m)^2
    resultsCV[b]<-sintegral(xseq,aaCV)$int/m2_int	
  }else{
    resultsCV[b]<-NA
  }
  
  print(paste("Replication", b, "out of", B))
  
}

kappas_B1_n70<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		# Concentration parameters selected with each method  on each data replication
ISE_B1_n70<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication

MISE_B1_n70<-colMeans(ISE_B1_n70,na.rm=TRUE)			# Monte Carlo averages of the approx. ISE (Table 2, row 6, columns 1 to 3)


##################
#### Model B1 ####
###  n = 100  ####
##################

n<-100  # Sample size
xseq<-seq(0,2*pi,length=100)	# grid where the regression function will be evaluated

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
  if(!is.na(kappaCV[b])){
    modCV<-(logistic_p1(xseq,x,y,kappaCV[b],startv))[[1]]
    aaCV<-(exp(modCV)/(1+exp(modCV))-m)^2
    resultsCV[b]<-sintegral(xseq,aaCV)$int/m2_int	
  }else{
    resultsCV[b]<-NA
  }
  
  print(paste("Replication", b, "out of", B))
  
}

kappas_B1_n100<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		# Concentration parameters selected with each method  on each data replication
ISE_B1_n100<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication

MISE_B1_n100<-colMeans(ISE_B1_n100,na.rm=TRUE)			# Monte Carlo averages of the approx. ISE (Table 2, row 7, columns 1 to 3)


##################
#### Model B1 ####
###  n = 250  ####
##################

n<-250  # Sample size
xseq<-seq(0,2*pi,length=100)	# grid where the regression function will be evaluated

lowerk<-0.3  # Minimum value possible for kappa
upperk<-70  # Maximum value possible for kappa

lowerk_ast<-0.3   # Minimum value possible for kappa star
upperk_ast<-15	# Maximum value possible for kappa star

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
  if(!is.na(kappaCV[b])){
    modCV<-(logistic_p1(xseq,x,y,kappaCV[b],startv))[[1]]
    aaCV<-(exp(modCV)/(1+exp(modCV))-m)^2
    resultsCV[b]<-sintegral(xseq,aaCV)$int/m2_int	
  }else{
    resultsCV[b]<-NA
  }
  
  print(paste("Replication", b, "out of", B))
  
}

kappas_B1_n250<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		# Concentration parameters selected with each method  on each data replication
ISE_B1_n250<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication

MISE_B1_n250<-colMeans(ISE_B1_n250,na.rm=TRUE)			# Monte Carlo averages of the approx. ISE (Table 2, row 8, columns 1 to 3)



##################
#### Model B1 ####
###  n = 500  ####
##################

n<-500  # Sample size
xseq<-seq(0,2*pi,length=100)	# grid where the regression function will be evaluated

lowerk<-0.3  # Minimum value possible for kappa
upperk<-100  # Maximum value possible for kappa

lowerk_ast<-0.3   # Minimum value possible for kappa star
upperk_ast<-18	# Maximum value possible for kappa star

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
  if(!is.na(kappaCV[b])){
    modCV<-(logistic_p1(xseq,x,y,kappaCV[b],startv))[[1]]
    aaCV<-(exp(modCV)/(1+exp(modCV))-m)^2
    resultsCV[b]<-sintegral(xseq,aaCV)$int/m2_int	
  }else{
    resultsCV[b]<-NA
  }
  
  print(paste("Replication", b, "out of", B))
  
}

kappas_B1_n500<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		# Concentration parameters selected with each method  on each data replication
ISE_B1_n500<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication

MISE_B1_n500<-colMeans(ISE_B1_n500,na.rm=TRUE)			# Monte Carlo averages of the approx. ISE (Table 2, row 9, columns 1 to 3)


##################
#### Model B1 ####
###  n = 1500 ####
##################

n<-1500  # Sample size
xseq<-seq(0,2*pi,length=100)	# grid where the regression function will be evaluated

lowerk<-0.3  # Minimum value possible for kappa
upperk<-200  # Maximum value possible for kappa

lowerk_ast<-0.3   # Minimum value possible for kappa star
upperk_ast<-30	# Maximum value possible for kappa star

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
  if(!is.na(kappaCV[b])){
    modCV<-(logistic_p1(xseq,x,y,kappaCV[b],startv))[[1]]
    aaCV<-(exp(modCV)/(1+exp(modCV))-m)^2
    resultsCV[b]<-sintegral(xseq,aaCV)$int/m2_int	
  }else{
    resultsCV[b]<-NA
  }
  
  print(paste("Replication", b, "out of", B))
  
}

kappas_B1_n1500<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		# Concentration parameters selected with each method  on each data replication
ISE_B1_n1500<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication

MISE_B1_n1500<-colMeans(ISE_B1_n1500,na.rm=TRUE)			# Monte Carlo averages of the approx. ISE (Table 2, row 10, columns 1 to 3)


##################
#### Model B2 ####
###   n = 70  ####
##################

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
  if(!is.na(kappaCV[b])){
    modCV<-(logistic_p1(xseq,x,y,kappaCV[b],startv))[[1]]
    aaCV<-(exp(modCV)/(1+exp(modCV))-m)^2
    resultsCV[b]<-sintegral(xseq,aaCV)$int/m2_int	
  }else{
    resultsCV[b]<-NA
  }
  
  
  print(paste("Replication", b, "out of", B))
  
}

kappas_B2_n70<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		# Concentration parameters selected with each method  on each data replication
ISE_B2_n70<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication


MISE_B2_n70<-colMeans(ISE_B2_n70,na.rm=TRUE)		# Monte Carlo averages of the approx. ISE (Table 2, row 6, columns 4 to 6)


##################
#### Model B2 ####
###   n = 100 ####
##################

n<-100 # Sample size
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
  if(!is.na(kappaCV[b])){
    modCV<-(logistic_p1(xseq,x,y,kappaCV[b],startv))[[1]]
    aaCV<-(exp(modCV)/(1+exp(modCV))-m)^2
    resultsCV[b]<-sintegral(xseq,aaCV)$int/m2_int	
  }else{
    resultsCV[b]<-NA
  }
  
  
  print(paste("Replication", b, "out of", B))
  
}

kappas_B2_n100<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		# Concentration parameters selected with each method  on each data replication
ISE_B2_n100<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication


MISE_B2_n100<-colMeans(ISE_B2_n100,na.rm=TRUE)		# Monte Carlo averages of the approx. ISE (Table 2, row 7, columns 4 to 6)


##################
#### Model B2 ####
###  n = 250  ####
##################

n<-250 # Sample size
xseq<-seq(0,2*pi,length=100)	# grid where the regression funciton will be evaluated

lowerk<-0.3  # Minimum value possible for kappa
upperk<-70  # Maximum value possible for kappa

lowerk_ast<-0.3   # Minimum value possible for kappa star
upperk_ast<-15	# Maximum value possible for kappa star

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
  if(!is.na(kappaCV[b])){
    modCV<-(logistic_p1(xseq,x,y,kappaCV[b],startv))[[1]]
    aaCV<-(exp(modCV)/(1+exp(modCV))-m)^2
    resultsCV[b]<-sintegral(xseq,aaCV)$int/m2_int	
  }else{
    resultsCV[b]<-NA
  }
  
  
  print(paste("Replication", b, "out of", B))
  
}

kappas_B2_n250<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		# Concentration parameters selected with each method  on each data replication
ISE_B2_n250<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication


MISE_B2_n250<-colMeans(ISE_B2_n250,na.rm=TRUE)		# Monte Carlo averages of the approx. ISE (Table 2, row 8, columns 4 to 6)


##################
#### Model B2 ####
###  n = 500  ####
##################

n<-500 # Sample size
xseq<-seq(0,2*pi,length=100)	# grid where the regression funciton will be evaluated

lowerk<-0.3  # Minimum value possible for kappa
upperk<-100  # Maximum value possible for kappa

lowerk_ast<-0.3   # Minimum value possible for kappa star
upperk_ast<-18	# Maximum value possible for kappa star

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
  if(!is.na(kappaCV[b])){
    modCV<-(logistic_p1(xseq,x,y,kappaCV[b],startv))[[1]]
    aaCV<-(exp(modCV)/(1+exp(modCV))-m)^2
    resultsCV[b]<-sintegral(xseq,aaCV)$int/m2_int	
  }else{
    resultsCV[b]<-NA
  }
  
  
  print(paste("Replication", b, "out of", B))
  
}

kappas_B2_n500<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		# Concentration parameters selected with each method  on each data replication
ISE_B2_n500<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication


MISE_B2_n500<-colMeans(ISE_B2_n500,na.rm=TRUE)		# Monte Carlo averages of the approx. ISE (Table 2, row 9, columns 4 to 6)


##################
#### Model B2 ####
###  n = 1500 ####
##################

n<-1500 # Sample size
xseq<-seq(0,2*pi,length=100)	# grid where the regression funciton will be evaluated

lowerk<-0.3  # Minimum value possible for kappa
upperk<-200  # Maximum value possible for kappa

lowerk_ast<-0.3   # Minimum value possible for kappa star
upperk_ast<-30	# Maximum value possible for kappa star

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
  if(!is.na(kappaCV[b])){
    modCV<-(logistic_p1(xseq,x,y,kappaCV[b],startv))[[1]]
    aaCV<-(exp(modCV)/(1+exp(modCV))-m)^2
    resultsCV[b]<-sintegral(xseq,aaCV)$int/m2_int	
  }else{
    resultsCV[b]<-NA
  }
  
  
  print(paste("Replication", b, "out of", B))
  
}

kappas_B2_n1500<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		# Concentration parameters selected with each method  on each data replication
ISE_B2_n1500<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication


MISE_B2_n1500<-colMeans(ISE_B2_n1500,na.rm=TRUE)		# Monte Carlo averages of the approx. ISE (Table 2, row 10, columns 4 to 6)


##################
#### Model G1 ####
###   n = 70  ####
##################

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

  if(!is.na(kappaPILOT[b])){
    modPILOT<-(gamma_p1(xseq,x,y,kappa_ast,startv))[[1]] 
    aaPILOT<-(modPILOT-m)^2	
    resultsPILOT[b]<-sintegral(xseq,aaPILOT)$int/m2_int	
  }else{
    resultsPILOT[b]<-NA
  }
  if(!is.na(kappaREFINED[b])){
    modREFINED<-(gamma_p1(xseq,x,y,kappa_opt,startv))[[1]]	
    aaREFINED<-(modREFINED-m)^2
    resultsREFINED[b]<-sintegral(xseq,aaREFINED)$int/m2_int
  }else{
    resultsREFINED[b]<-NA
  }
  if(!is.na(kappaCV[b])){
    modCV<-(gamma_p1(xseq,x,y,kappa_cv,startv))[[1]]		
    aaCV<-(modCV-m)^2			
    resultsCV[b]<-sintegral(xseq,aaCV)$int/m2_int		
  }else{
    resultsCV[b]<-NA
  }
  
  print(paste("Replication", b, "out of", B))
  
}

kappas_G1_n70<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		# Concentration parameters selected with each method  on each data replication
ISE_G1_n70<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication

MISE_G1_n70<-colMeans(ISE_G1_n70,na.rm=TRUE)		# Monte Carlo averages of the approx. ISE (Table 2, row 16, columns 1 to 3)

##################
#### Model G1 ####
###  n = 100  ####
##################

n<-100  # Sample size
xseq<-seq(0,2*pi,length=100)	# grid where the regression funciton will be evaluated
nu<-1

lowerk<-0.3  # Minimum value possible for kappa
upperk<-100  # Maximum value possible for kappa

lowerk_ast<-0.3   	# Minimum value possible for kappa star
upperk_ast<-12	 	# Maximum value possible for kappa star


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
  
  if(!is.na(kappaPILOT[b])){
    modPILOT<-(gamma_p1(xseq,x,y,kappa_ast,startv))[[1]] 
    aaPILOT<-(modPILOT-m)^2	
    resultsPILOT[b]<-sintegral(xseq,aaPILOT)$int/m2_int	
  }else{
    resultsPILOT[b]<-NA
  }
  if(!is.na(kappaREFINED[b])){
    modREFINED<-(gamma_p1(xseq,x,y,kappa_opt,startv))[[1]]	
    aaREFINED<-(modREFINED-m)^2
    resultsREFINED[b]<-sintegral(xseq,aaREFINED)$int/m2_int
  }else{
    resultsREFINED[b]<-NA
  }
  if(!is.na(kappaCV[b])){
    modCV<-(gamma_p1(xseq,x,y,kappa_cv,startv))[[1]]		
    aaCV<-(modCV-m)^2			
    resultsCV[b]<-sintegral(xseq,aaCV)$int/m2_int		
  }else{
    resultsCV[b]<-NA
  }
  
  print(paste("Replication", b, "out of", B))
  
}

kappas_G1_n100<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		# Concentration parameters selected with each method  on each data replication
ISE_G1_n100<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication

MISE_G1_n100<-colMeans(ISE_G1_n100,na.rm=TRUE)		# Monte Carlo averages of the approx. ISE (Table 2, row 17, columns 1 to 3)


##################
#### Model G1 ####
###  n = 250  ####
##################

n<-250  # Sample size
xseq<-seq(0,2*pi,length=100)	# grid where the regression funciton will be evaluated
nu<-1

lowerk<-0.3  # Minimum value possible for kappa
upperk<-200  # Maximum value possible for kappa

lowerk_ast<-0.3   	# Minimum value possible for kappa star
upperk_ast<-20	 	# Maximum value possible for kappa star


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
  
  if(!is.na(kappaPILOT[b])){
    modPILOT<-(gamma_p1(xseq,x,y,kappa_ast,startv))[[1]] 
    aaPILOT<-(modPILOT-m)^2	
    resultsPILOT[b]<-sintegral(xseq,aaPILOT)$int/m2_int	
  }else{
    resultsPILOT[b]<-NA
  }
  if(!is.na(kappaREFINED[b])){
    modREFINED<-(gamma_p1(xseq,x,y,kappa_opt,startv))[[1]]	
    aaREFINED<-(modREFINED-m)^2
    resultsREFINED[b]<-sintegral(xseq,aaREFINED)$int/m2_int
  }else{
    resultsREFINED[b]<-NA
  }
  if(!is.na(kappaCV[b])){
    modCV<-(gamma_p1(xseq,x,y,kappa_cv,startv))[[1]]		
    aaCV<-(modCV-m)^2			
    resultsCV[b]<-sintegral(xseq,aaCV)$int/m2_int		
  }else{
    resultsCV[b]<-NA
  }
  
  print(paste("Replication", b, "out of", B))

}

kappas_G1_n250<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		# Concentration parameters selected with each method  on each data replication
ISE_G1_n250<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication

MISE_G1_n250<-colMeans(ISE_G1_n250,na.rm=TRUE)		# Monte Carlo averages of the approx. ISE (Table 2, row 18, columns 1 to 3)



##################
#### Model G1 ####
###  n = 500  ####
##################

n<-500  # Sample size
xseq<-seq(0,2*pi,length=100)	# grid where the regression funciton will be evaluated
nu<-1

lowerk<-0.3  # Minimum value possible for kappa
upperk<-200  # Maximum value possible for kappa

lowerk_ast<-0.3   	# Minimum value possible for kappa star
upperk_ast<-20	 	# Maximum value possible for kappa star


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
  
  if(!is.na(kappaPILOT[b])){
    modPILOT<-(gamma_p1(xseq,x,y,kappa_ast,startv))[[1]] 
    aaPILOT<-(modPILOT-m)^2	
    resultsPILOT[b]<-sintegral(xseq,aaPILOT)$int/m2_int	
  }else{
    resultsPILOT[b]<-NA
  }
  if(!is.na(kappaREFINED[b])){
    modREFINED<-(gamma_p1(xseq,x,y,kappa_opt,startv))[[1]]	
    aaREFINED<-(modREFINED-m)^2
    resultsREFINED[b]<-sintegral(xseq,aaREFINED)$int/m2_int
  }else{
    resultsREFINED[b]<-NA
  }
  if(!is.na(kappaCV[b])){
    modCV<-(gamma_p1(xseq,x,y,kappa_cv,startv))[[1]]		
    aaCV<-(modCV-m)^2			
    resultsCV[b]<-sintegral(xseq,aaCV)$int/m2_int		
  }else{
    resultsCV[b]<-NA
  }
  
  print(paste("Replication", b, "out of", B))

}

kappas_G1_n500<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		# Concentration parameters selected with each method  on each data replication
ISE_G1_n500<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication

MISE_G1_n500<-colMeans(ISE_G1_n500,na.rm=TRUE)		# Monte Carlo averages of the approx. ISE (Table 2, row 19, columns 1 to 3)


##################
#### Model G1 ####
###  n = 1500 ####
##################

n<-1500  # Sample size
xseq<-seq(0,2*pi,length=100)	# grid where the regression funciton will be evaluated
nu<-1

lowerk<-0.3  # Minimum value possible for kappa
upperk<-250  # Maximum value possible for kappa

lowerk_ast<-0.3   	# Minimum value possible for kappa star
upperk_ast<-40	 	# Maximum value possible for kappa star


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
  
  if(!is.na(kappaPILOT[b])){
    modPILOT<-(gamma_p1(xseq,x,y,kappa_ast,startv))[[1]] 
    aaPILOT<-(modPILOT-m)^2	
    resultsPILOT[b]<-sintegral(xseq,aaPILOT)$int/m2_int	
  }else{
    resultsPILOT[b]<-NA
  }
  if(!is.na(kappaREFINED[b])){
    modREFINED<-(gamma_p1(xseq,x,y,kappa_opt,startv))[[1]]	
    aaREFINED<-(modREFINED-m)^2
    resultsREFINED[b]<-sintegral(xseq,aaREFINED)$int/m2_int
  }else{
    resultsREFINED[b]<-NA
  }
  if(!is.na(kappaCV[b])){
    modCV<-(gamma_p1(xseq,x,y,kappa_cv,startv))[[1]]		
    aaCV<-(modCV-m)^2			
    resultsCV[b]<-sintegral(xseq,aaCV)$int/m2_int		
  }else{
    resultsCV[b]<-NA
  }
  
  print(paste("Replication", b, "out of", B))

}

kappas_G1_n1500<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		# Concentration parameters selected with each method  on each data replication
ISE_G1_n1500<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication

MISE_G1_n1500<-colMeans(ISE_G1_n1500,na.rm=TRUE)		# Monte Carlo averages of the approx. ISE (Table 2, row 20, columns 1 to 3)


##################
#### Model G2 ####
###   n = 70  ####
##################

n<-70  # Sample size
xseq<-seq(0,2*pi,length=100)	# grid where the regression funciton will be evaluated
nu<-2

lowerk<-0.3  # Minimum value possible for kappa
upperk<-100  # Maximum value possible for kappa

lowerk_ast<-0.3   	# Minimum value possible for kappa star
upperk_ast<-12	 	  # Maximum value possible for kappa star


m<-log(5+2*cos(3+1.5*sin(xseq))) 		# Regression function evaluated at xseq
m2_int<-sintegral(xseq,m^2)$int		# Integral of m^2

B<-500  # Number of replicates
resultsPILOT<-numeric(B)     	# The approx. ISE with the pilot method will be stored here
resultsREFINED<-numeric(B)		# The approx. ISE with the refined method will be stored here
resultsCV<-numeric(B)			    # The approx. ISE with the cross-validation method will be stored here
kappaPILOT<-numeric(B)			  # The selected kappas with the pilot method will be stored here
kappaREFINED<-numeric(B)		  # The selected kappas with the refined method will be stored here
kappaCV<-numeric(B)			      # The selected kappas with the cross-validation method will be stored here

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
  if(!is.na(kappaPILOT[b])){
    modPILOT<-(gamma_p1(xseq,x,y,kappa_ast,startv))[[1]] 
    aaPILOT<-(modPILOT-m)^2	
    resultsPILOT[b]<-sintegral(xseq,aaPILOT)$int/m2_int	
  }else{
    resultsPILOT[b]<-NA
  }
  if(!is.na(kappaREFINED[b])){
    modREFINED<-(gamma_p1(xseq,x,y,kappa_opt,startv))[[1]]	
    aaREFINED<-(modREFINED-m)^2
    resultsREFINED[b]<-sintegral(xseq,aaREFINED)$int/m2_int
  }else{
    resultsREFINED[b]<-NA
  }
  if(!is.na(kappaCV[b])){
    modCV<-(gamma_p1(xseq,x,y,kappa_cv,startv))[[1]]		
    aaCV<-(modCV-m)^2			
    resultsCV[b]<-sintegral(xseq,aaCV)$int/m2_int		
  }else{
    resultsCV[b]<-NA
  }
  
  print(paste("Replication", b, "out of", B))

}

kappas_G2_n70<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		# Concentration parameters selected with each method  on each data replication
ISE_G2_n70<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication

MISE_G2_n70<-colMeans(ISE_G2_n70,na.rm=TRUE)		# Monte Carlo averages of the approx. ISE (Table 2, row 16, columns 4 to 6)


##################
#### Model G2 ####
###  n = 100  ####
##################

n<-100  # Sample size
xseq<-seq(0,2*pi,length=100)	# grid where the regression funciton will be evaluated
nu<-2

lowerk<-0.3  # Minimum value possible for kappa
upperk<-100  # Maximum value possible for kappa

lowerk_ast<-0.3   	# Minimum value possible for kappa star
upperk_ast<-12	 	  # Maximum value possible for kappa star


m<-log(5+2*cos(3+1.5*sin(xseq))) 		# Regression function evaluated at xseq
m2_int<-sintegral(xseq,m^2)$int		# Integral of m^2

B<-500  # Number of replicates
resultsPILOT<-numeric(B)     	# The approx. ISE with the pilot method will be stored here
resultsREFINED<-numeric(B)		# The approx. ISE with the refined method will be stored here
resultsCV<-numeric(B)			    # The approx. ISE with the cross-validation method will be stored here
kappaPILOT<-numeric(B)			  # The selected kappas with the pilot method will be stored here
kappaREFINED<-numeric(B)		  # The selected kappas with the refined method will be stored here
kappaCV<-numeric(B)			      # The selected kappas with the cross-validation method will be stored here

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
  if(!is.na(kappaPILOT[b])){
    modPILOT<-(gamma_p1(xseq,x,y,kappa_ast,startv))[[1]] 
    aaPILOT<-(modPILOT-m)^2	
    resultsPILOT[b]<-sintegral(xseq,aaPILOT)$int/m2_int	
  }else{
    resultsPILOT[b]<-NA
  }
  if(!is.na(kappaREFINED[b])){
    modREFINED<-(gamma_p1(xseq,x,y,kappa_opt,startv))[[1]]	
    aaREFINED<-(modREFINED-m)^2
    resultsREFINED[b]<-sintegral(xseq,aaREFINED)$int/m2_int
  }else{
    resultsREFINED[b]<-NA
  }
  if(!is.na(kappaCV[b])){
    modCV<-(gamma_p1(xseq,x,y,kappa_cv,startv))[[1]]		
    aaCV<-(modCV-m)^2			
    resultsCV[b]<-sintegral(xseq,aaCV)$int/m2_int		
  }else{
    resultsCV[b]<-NA
  }
  
  print(paste("Replication", b, "out of", B))
  
}

kappas_G2_n100<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		# Concentration parameters selected with each method  on each data replication
ISE_G2_n100<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication

MISE_G2_n100<-colMeans(ISE_G2_n100,na.rm=TRUE)		# Monte Carlo averages of the approx. ISE (Table 2, row 17, columns 4 to 6)



##################
#### Model G2 ####
###  n = 250  ####
##################

n<-250  # Sample size
xseq<-seq(0,2*pi,length=100)	# grid where the regression funciton will be evaluated
nu<-2

lowerk<-0.3  # Minimum value possible for kappa
upperk<-100  # Maximum value possible for kappa

lowerk_ast<-0.3   	# Minimum value possible for kappa star
upperk_ast<-12	 	  # Maximum value possible for kappa star


m<-log(5+2*cos(3+1.5*sin(xseq))) 		# Regression function evaluated at xseq
m2_int<-sintegral(xseq,m^2)$int		# Integral of m^2

B<-500  # Number of replicates
resultsPILOT<-numeric(B)     	# The approx. ISE with the pilot method will be stored here
resultsREFINED<-numeric(B)		# The approx. ISE with the refined method will be stored here
resultsCV<-numeric(B)			    # The approx. ISE with the cross-validation method will be stored here
kappaPILOT<-numeric(B)			  # The selected kappas with the pilot method will be stored here
kappaREFINED<-numeric(B)		  # The selected kappas with the refined method will be stored here
kappaCV<-numeric(B)			      # The selected kappas with the cross-validation method will be stored here

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
  if(!is.na(kappaPILOT[b])){
    modPILOT<-(gamma_p1(xseq,x,y,kappa_ast,startv))[[1]] 
    aaPILOT<-(modPILOT-m)^2	
    resultsPILOT[b]<-sintegral(xseq,aaPILOT)$int/m2_int	
  }else{
    resultsPILOT[b]<-NA
  }
  if(!is.na(kappaREFINED[b])){
    modREFINED<-(gamma_p1(xseq,x,y,kappa_opt,startv))[[1]]	
    aaREFINED<-(modREFINED-m)^2
    resultsREFINED[b]<-sintegral(xseq,aaREFINED)$int/m2_int
  }else{
    resultsREFINED[b]<-NA
  }
  if(!is.na(kappaCV[b])){
    modCV<-(gamma_p1(xseq,x,y,kappa_cv,startv))[[1]]		
    aaCV<-(modCV-m)^2			
    resultsCV[b]<-sintegral(xseq,aaCV)$int/m2_int		
  }else{
    resultsCV[b]<-NA
  }
  
  print(paste("Replication", b, "out of", B))
  
}

kappas_G2_n250<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		# Concentration parameters selected with each method  on each data replication
ISE_G2_n250<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication

MISE_G2_n250<-colMeans(ISE_G2_n250,na.rm=TRUE)		# Monte Carlo averages of the approx. ISE (Table 2, row 18, columns 4 to 6)


##################
#### Model G2 ####
###  n = 500  ####
##################

n<-500  # Sample size
xseq<-seq(0,2*pi,length=100)	# grid where the regression funciton will be evaluated
nu<-2

lowerk<-0.3  # Minimum value possible for kappa
upperk<-150  # Maximum value possible for kappa

lowerk_ast<-0.3   	# Minimum value possible for kappa star
upperk_ast<-16	 	  # Maximum value possible for kappa star


m<-log(5+2*cos(3+1.5*sin(xseq))) 		# Regression function evaluated at xseq
m2_int<-sintegral(xseq,m^2)$int		# Integral of m^2

B<-500  # Number of replicates
resultsPILOT<-numeric(B)     	# The approx. ISE with the pilot method will be stored here
resultsREFINED<-numeric(B)		# The approx. ISE with the refined method will be stored here
resultsCV<-numeric(B)			    # The approx. ISE with the cross-validation method will be stored here
kappaPILOT<-numeric(B)			  # The selected kappas with the pilot method will be stored here
kappaREFINED<-numeric(B)		  # The selected kappas with the refined method will be stored here
kappaCV<-numeric(B)			      # The selected kappas with the cross-validation method will be stored here

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
  if(!is.na(kappaPILOT[b])){
    modPILOT<-(gamma_p1(xseq,x,y,kappa_ast,startv))[[1]] 
    aaPILOT<-(modPILOT-m)^2	
    resultsPILOT[b]<-sintegral(xseq,aaPILOT)$int/m2_int	
  }else{
    resultsPILOT[b]<-NA
  }
  if(!is.na(kappaREFINED[b])){
    modREFINED<-(gamma_p1(xseq,x,y,kappa_opt,startv))[[1]]	
    aaREFINED<-(modREFINED-m)^2
    resultsREFINED[b]<-sintegral(xseq,aaREFINED)$int/m2_int
  }else{
    resultsREFINED[b]<-NA
  }
  if(!is.na(kappaCV[b])){
    modCV<-(gamma_p1(xseq,x,y,kappa_cv,startv))[[1]]		
    aaCV<-(modCV-m)^2			
    resultsCV[b]<-sintegral(xseq,aaCV)$int/m2_int		
  }else{
    resultsCV[b]<-NA
  }
  
  print(paste("Replication", b, "out of", B))
  
}

kappas_G2_n500<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		# Concentration parameters selected with each method  on each data replication
ISE_G2_n500<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication

MISE_G2_n500<-colMeans(ISE_G2_n500,na.rm=TRUE)		# Monte Carlo averages of the approx. ISE (Table 2, row 19, columns 4 to 6)


##################
#### Model G2 ####
###  n = 1500 ####
##################

n<-1500  # Sample size
xseq<-seq(0,2*pi,length=100)	# grid where the regression funciton will be evaluated
nu<-2

lowerk<-0.3  # Minimum value possible for kappa
upperk<-150  # Maximum value possible for kappa

lowerk_ast<-0.3   	# Minimum value possible for kappa star
upperk_ast<-16	 	  # Maximum value possible for kappa star


m<-log(5+2*cos(3+1.5*sin(xseq))) 		# Regression function evaluated at xseq
m2_int<-sintegral(xseq,m^2)$int		# Integral of m^2

B<-500  # Number of replicates
resultsPILOT<-numeric(B)     	# The approx. ISE with the pilot method will be stored here
resultsREFINED<-numeric(B)		# The approx. ISE with the refined method will be stored here
resultsCV<-numeric(B)			    # The approx. ISE with the cross-validation method will be stored here
kappaPILOT<-numeric(B)			  # The selected kappas with the pilot method will be stored here
kappaREFINED<-numeric(B)		  # The selected kappas with the refined method will be stored here
kappaCV<-numeric(B)			      # The selected kappas with the cross-validation method will be stored here

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
  if(!is.na(kappaPILOT[b])){
    modPILOT<-(gamma_p1(xseq,x,y,kappa_ast,startv))[[1]] 
    aaPILOT<-(modPILOT-m)^2	
    resultsPILOT[b]<-sintegral(xseq,aaPILOT)$int/m2_int	
  }else{
    resultsPILOT[b]<-NA
  }
  if(!is.na(kappaREFINED[b])){
    modREFINED<-(gamma_p1(xseq,x,y,kappa_opt,startv))[[1]]	
    aaREFINED<-(modREFINED-m)^2
    resultsREFINED[b]<-sintegral(xseq,aaREFINED)$int/m2_int
  }else{
    resultsREFINED[b]<-NA
  }
  if(!is.na(kappaCV[b])){
    modCV<-(gamma_p1(xseq,x,y,kappa_cv,startv))[[1]]		
    aaCV<-(modCV-m)^2			
    resultsCV[b]<-sintegral(xseq,aaCV)$int/m2_int		
  }else{
    resultsCV[b]<-NA
  }
  
  print(paste("Replication", b, "out of", B))
  
}

kappas_G2_n1500<-cbind(kappaPILOT,kappaREFINED,kappaCV) 		# Concentration parameters selected with each method  on each data replication
ISE_G2_n1500<-cbind(resultsPILOT,resultsREFINED,resultsCV)		# Approx. ISE with each method on each data replication

MISE_G2_n1500<-colMeans(ISE_G2_n1500,na.rm=TRUE)		# Monte Carlo averages of the approx. ISE (Table 2, row 20, columns 4 to 6)



####################################################################
##              TABLE WITH SIMULATION RESULTS                     ##  
####################################################################

MISE_N1_n100<-c(1,2,3)
MISE_N2_n100<-c(4,5,6)

C1<-rbind(MISE_N1_n70,MISE_N1_n100,MISE_N1_n250,MISE_N1_n500,MISE_N1_n1500,
          MISE_B1_n70,MISE_B1_n100,MISE_B1_n250,MISE_B1_n500,MISE_B1_n1500,
          MISE_P1_n70,MISE_P1_n100,MISE_P1_n250,MISE_P1_n500,MISE_P1_n1500,
          MISE_G1_n70,MISE_G1_n100,MISE_G1_n250,MISE_G1_n500,MISE_G1_n1500)
C2<-rbind(MISE_N2_n70,MISE_N2_n100,MISE_N2_n250,MISE_N2_n500,MISE_N2_n1500,
          MISE_B2_n70,MISE_B2_n100,MISE_B2_n250,MISE_B2_n500,MISE_B2_n1500,
          MISE_P2_n70,MISE_P2_n100,MISE_P2_n250,MISE_P2_n500,MISE_P2_n1500,
          MISE_G2_n70,MISE_G2_n100,MISE_G2_n250,MISE_G2_n500,MISE_G2_n1500)

TABLE2<-rbind(C1,C2)
TABLE2



####################################################################
####################################################################
#                                                                  #
#  FIGURE 2      (Representative estimators n=250)                 #
#                                                                  #
####################################################################
####################################################################


library(plotrix)
library(circular)

################
### FIGURE 2 (a)
################
iN1n250<-MISE_N1_n250
kN1n250<-kappas_N1_n250
qq<-quantile(iN1n250[,2],prob=c(0.05,0.5,0.95)) # Finding the representatives

# Representative quantile 0.05
qq[1]
ind<-which(iN1n100[,2]>0.03073&iN1n100[,2]<0.03078) 
# Regression function
xseq<-seq(0,2*pi,length=300)
m<-sin(2*xseq)*cos(xseq) 
# We replicate the data again
n<-250
b<-ind
set.seed(888*b)
x<-rcircularuniform(n)
y<-sin(2*x)*cos(x)+rnorm(n,sd=0.35)
modREFINED1<-kern.reg.circ.lin(x,y,t=xseq,bw=kN1n250[ind,2]) # estimation representative 0.05

# Representative quantile 0.5
qq[2]
ind2<-which(iN1n250[,2]>0.03027&iN1n250[,2]<0.03030)
b<-ind2
set.seed(888*b)
x<-rcircularuniform(n)
y<-sin(2*x)*cos(x)+rnorm(n,sd=0.35)
modREFINED2<-kern.reg.circ.lin(x,y,t=xseq,bw=kN1n250[ind2,2]) # estimation representative 0.5

# Representative quantile 0.95
qq[3]
ind3<-which(iN1n250[,2]>0.04985&iN1n250[,2]<0.04993)
b<-ind3
set.seed(888*b)
x<-rcircularuniform(n)
y<-sin(2*x)*cos(x)+rnorm(n,sd=0.35)
modREFINED3<-kern.reg.circ.lin(x,y,t=xseq,bw=kN1n250[ind3,2]) # estimation representative 0.95

# radial representation (FIGURE 2 (a))
radial.plot(as.numeric(m),circular(xseq),rp.type="l",lwd=3,radial.lim=c(-1.15,1),labels=c(expression(0),expression(pi/4),expression(pi/2),expression(3*pi/4),expression(pi),expression(5*pi/4),expression(3*pi/2),expression(7*pi/4)),label.pos=c(0,pi/4,pi/2,3*pi/4,pi,5*pi/4,3*pi/2,7*pi/4))
radial.plot(modREFINED1$y,modREFINED1$x,rp.type="l",add=T,lwd=3,lty=2,line.col=2,radial.lim=c(-1.15,1))
radial.plot(modREFINED2$y,modREFINED2$x,rp.type="l",add=T,lwd=3,lty=3,line.col=4,radial.lim=c(-1.15,1))
radial.plot(modREFINED3$y,modREFINED3$x,rp.type="l",add=T,lwd=3,lty=4,line.col=3,radial.lim=c(-1.15,1))


################
### FIGURE 2 (b)
################

iN2n250<-MISE_N2_n250
kN2n250<-kappas_N2_n250

qq<-quantile(iN2n250[,2],prob=c(0.05,0.5,0.95))

# Representative quantile 0.05
qq[1]
ind<-which(iN2n250[,2]>0.00619&iN2n250[,2]<0.0062)
xseq<-seq(0,2*pi,length=300)
m<-1.75*cos(1*xseq-pi)*sin(1*xseq)+cos(xseq)
n<-250
b<-ind
set.seed(8887*b)
x<-rcircularuniform(n)
y<-1.75*cos(1*x-pi)*sin(1*x)+cos(x)+rnorm(n,sd=0.5)
modREFINED1<-kern.reg.circ.lin(x,y,t=xseq,bw=kN2n250[ind,2])

# Representative quantile 0.5
qq[2]
ind2<-which(iN2n250[,2]>0.01283&iN2n250[,2]<0.01285)
b<-ind2
set.seed(8887*b)
x<-rcircularuniform(n)
y<-1.75*cos(1*x-pi)*sin(1*x)+cos(x)+rnorm(n,sd=0.5)
modREFINED2<-kern.reg.circ.lin(x,y,t=xseq,bw=kN2n250[ind2,2])


# Representative quantile 0.95
qq[3]
ind3<-which(iN2n250[,2]>0.02401&iN2n250[,2]<0.02406)
b<-ind3
set.seed(8887*b)
x<-rcircularuniform(n)
y<-1.75*cos(1*x-pi)*sin(1*x)+cos(x)+rnorm(n,sd=0.5)
modREFINED3<-kern.reg.circ.lin(x,y,t=xseq,bw=kN2n250[ind3,2])

# radial representation
radial.plot(as.numeric(m),circular(xseq),rp.type="l",lwd=3,radial.lim=c(-1.7,1.7),labels=c(expression(0),expression(pi/4),expression(pi/2),expression(3*pi/4),expression(pi),expression(5*pi/4),expression(3*pi/2),expression(7*pi/4)),label.pos=c(0,pi/4,pi/2,3*pi/4,pi,5*pi/4,3*pi/2,7*pi/4))
radial.plot(modREFINED1$y,modREFINED1$x,rp.type="l",add=T,lwd=3,lty=2,line.col=2,radial.lim=c(-1.7,1.7))
radial.plot(modREFINED2$y,modREFINED2$x,rp.type="l",add=T,lwd=3,lty=3,line.col=4,radial.lim=c(-1.7,1.7))
radial.plot(modREFINED3$y,modREFINED3$x,rp.type="l",add=T,lwd=3,lty=4,line.col=3,radial.lim=c(-1.7,1.7))


################
### FIGURE 2 (c)
################

iB1n250<-MISE_B1_n250
kB1n250<-kappas_B1_n250

qq<-quantile(iB1n250[,2],prob=c(0.05,0.5,0.95))


# Representative quantile 0.05
qq[1]
ind<-which(iB1n250[,2]>0.01105&iB1n250[,2]<0.0111)
xseq<-seq(0,2*pi,length=300)
m1<-2*sin(xseq)*cos(2*xseq)
m<-exp(m1)/(1+exp(m1))
n<-250
b<-ind
set.seed(888*b)
x<-rcircularuniform(n)
logit<-2*sin(x)*cos(2*x)
y<-rbinom(n,size=1, prob=exp(logit)/(1+exp(logit)))
startv<-c(log(mean(y)/(1-mean(y))),0,0,0)
modREFINED1<-(logistic_p1(xseq,x,y,kB1n250[ind,2],startv))[[1]]

# Representative quantile 0.5
qq[2]
ind2<-which(iB1n250[,2]>0.0278&iB1n250[,2]<0.0279)
b<-ind2
set.seed(888*b)
x<-rcircularuniform(n)
logit<-2*sin(x)*cos(2*x)
y<-rbinom(n,size=1, prob=exp(logit)/(1+exp(logit)))
startv<-c(log(mean(y)/(1-mean(y))),0,0,0)
modREFINED2<-(logistic_p1(xseq,x,y,kB1n250[ind2,2],startv))[[1]]

# Representative quantile 0.95
qq[3]
ind3<-which(iB1n250[,2]>0.053665&iB1n250[,2]<0.05368)
b<-ind3
set.seed(888*(b+1))
x<-rcircularuniform(n)
logit<-2*sin(x)*cos(2*x)
y<-rbinom(n,size=1, prob=exp(logit)/(1+exp(logit)))
startv<-c(log(mean(y)/(1-mean(y))),0,0,0)
modREFINED3<-(logistic_p1(xseq,x,y,kB1n250[ind3,2],startv))[[1]]

# radial representation
radial.plot(as.numeric(m),circular(xseq),rp.type="l",lwd=3,radial.lim=c(0,1),labels=c(expression(0),expression(pi/4),expression(pi/2),expression(3*pi/4),expression(pi),expression(5*pi/4),expression(3*pi/2),expression(7*pi/4)),label.pos=c(0,pi/4,pi/2,3*pi/4,pi,5*pi/4,3*pi/2,7*pi/4))
radial.plot(exp(modREFINED1)/(1+exp(modREFINED1)),xseq,rp.type="l",add=T,lwd=3,lty=2,line.col=2,radial.lim=c(0,1))
radial.plot(exp(modREFINED2)/(1+exp(modREFINED2)),xseq,rp.type="l",add=T,lwd=3,lty=3,line.col=4,radial.lim=c(0,1))
radial.plot(exp(modREFINED3)/(1+exp(modREFINED3)),xseq,rp.type="l",add=T,lwd=3,lty=4,line.col=3,radial.lim=c(0,1))


################
### FIGURE 2 (d)
################

iB2n250<-MISE_B2_n250
kB2n250<-kappas_B2_n250

qq<-quantile(iB2n250[,2],prob=c(0.05,0.5,0.95))

# Representative quantile 0.05
qq[1]
ind<-which(iB2n250[,2]>0.00383&iB2n250[,2]<0.00384)
xseq<-seq(0,2*pi,length=300)
m1<-log(1.6+1.5*sin(xseq)+exp(cos(xseq))/10)
m<-exp(m1)/(1+exp(m1))
n<-250
b<-ind
set.seed(888*b)
x<-rcircularuniform(n)
logit<-log(1.6+1.5*sin(x)+exp(cos(x))/10)
y<-rbinom(n,size=1, prob=exp(logit)/(1+exp(logit)))
startv<-c(log(mean(y)/(1-mean(y))),0,0,0)
modREFINED1<-(logistic_p1(xseq,x,y,kB2n250[ind,2],startv))[[1]]

# Representative quantile 0.5
qq[2]
ind2<-which(iB2n250[,2]>0.011215&iB2n250[,2]<0.011218)
b<-ind2
set.seed(888*b)
x<-rcircularuniform(n)
logit<-log(1.6+1.5*sin(x)+exp(cos(x))/10)
y<-rbinom(n,size=1, prob=exp(logit)/(1+exp(logit)))
startv<-c(log(mean(y)/(1-mean(y))),0,0,0)
modREFINED2<-(logistic_p1(xseq,x,y,kB2n250[ind2,2],startv))[[1]]

# Representative quantile 0.95
qq[3]
ind3<-which(iB2n250[,2]>0.0283&iB2n250[,2]<0.0284)
b<-ind3
set.seed(888*b)
x<-rcircularuniform(n)
logit<-log(1.6+1.5*sin(x)+exp(cos(x))/10)
y<-rbinom(n,size=1, prob=exp(logit)/(1+exp(logit)))
startv<-c(log(mean(y)/(1-mean(y))),0,0,0)
modREFINED3<-(logistic_p1(xseq,x,y,kB2n250[ind3,2],startv))[[1]]

# radial representation
radial.plot(as.numeric(m),circular(xseq),rp.type = "l",lwd=3,radial.lim=c(0,1),labels=c(expression(0),expression(pi/4),expression(pi/2),expression(3*pi/4),expression(pi),expression(5*pi/4),expression(3*pi/2),expression(7*pi/4)),label.pos=c(0,pi/4,pi/2,3*pi/4,pi,5*pi/4,3*pi/2,7*pi/4))
radial.plot(exp(modREFINED1)/(1+exp(modREFINED1)),xseq,rp.type="l",add=T,lwd=3,lty=2,line.col=2,radial.lim=c(0,1))
radial.plot(exp(modREFINED2)/(1+exp(modREFINED2)),xseq,rp.type="l",add=T,lwd=3,lty=3,line.col=4,radial.lim=c(0,1))
radial.plot(exp(modREFINED3)/(1+exp(modREFINED3)),xseq,rp.type="l",add=T,lwd=3,lty=4,line.col=3,radial.lim=c(0,1))


################
### FIGURE 2 (e)
################

iP1n250<-MISE_P1_n250
kP1n250<-kappas_P1_n250

qq<-quantile(iP1n250[,2],prob=c(0.05,0.5,0.95))

# Representative quantile 0.05
qq[1]
ind<-which(iP1n250[,2]>0.0007542&iP1n250[,2]<0.0007547)
xseq<-seq(0,2*pi,length=300)
m<-log(5+exp(1.5*sin(2*xseq-3)))
n<-250
b<-ind
set.seed(888*b)
x<-rcircularuniform(n)
meanpois<-5+exp(1.5*sin(2*x-3))
y<-rpois(n,lambda=meanpois)
startv<-c(log(mean(y)),0,0,0)
modREFINED1<-(pois_p1(xseq,x,y,kP1n250[ind,2],startv))[[1]]

# Representative quantile 0.5
qq[2]
ind2<-which(iP1n250[,2]>0.001528&iP1n250[,2]<0.001531)
b<-ind2
set.seed(888*b)
x<-rcircularuniform(n)
meanpois<-5+exp(1.5*sin(2*x-3))
y<-rpois(n,lambda=meanpois)
startv<-c(log(mean(y)),0,0,0)
modREFINED2<-(pois_p1(xseq,x,y,kP1n250[ind2,2],startv))[[1]]

# Representative quantile 0.95
qq[3]
ind3<-which(iP1n250[,2]>0.00289&iP1n250[,2]<0.0029)
b<-ind3
set.seed(888*b)
x<-rcircularuniform(n)
meanpois<-5+exp(1.5*sin(2*x-3))
y<-rpois(n,lambda=meanpois)
startv<-c(log(mean(y)),0,0,0)
modREFINED3<-(pois_p1(xseq,x,y,kP1n250[ind3,2],startv))[[1]]

# radial representation
radial.plot(as.numeric(m),circular(xseq),rp.type = "l",lwd=3,radial.lim=c(1.5,2.3),labels=c(expression(0),expression(pi/4),expression(pi/2),expression(3*pi/4),expression(pi),expression(5*pi/4),expression(3*pi/2),expression(7*pi/4)),label.pos=c(0,pi/4,pi/2,3*pi/4,pi,5*pi/4,3*pi/2,7*pi/4))
radial.plot(modREFINED1,xseq,rp.type="l",add=T,lwd=3,lty=2,line.col=2,radial.lim=c(1.5,2.3))
radial.plot(modREFINED2,xseq,rp.type="l",add=T,lwd=3,lty=3,line.col=4,radial.lim=c(1.5,2.3))
radial.plot(modREFINED3,xseq,rp.type="l",add=T,lwd=3,lty=4,line.col=3,radial.lim=c(1.5,2.3))


################
### FIGURE 2 (f)
################

iP2n250<-MISE_P2_n250
kP2n250<-kappas_P2_n250

qq<-quantile(iP2n250[,2],prob=c(0.05,0.5,0.95))

# Representative quantile 0.05
qq[1]
ind<-which(iP2n250[,2]>0.0001075&iP2n250[,2]<0.0001079)
xseq<-seq(0,2*pi,length=300)
m<-log(40+20*sin(3*xseq))
n<-250
b<-ind
set.seed(8887*b)
x<-rcircularuniform(n)
meanpois<-40+20*sin(3*x)
y<-rpois(n,lambda=meanpois)
startv<-c(log(mean(y)),0,0,0)
modREFINED1<-(pois_p1(xseq,x,y,kP2n250[ind,2],startv))[[1]]

# Representative quantile 0.5
qq[1]
ind2<-which(iP2n250[,2]>0.000187765&iP2n250[,2]<0.00018787)
b<-ind2
set.seed(8887*b)
x<-rcircularuniform(n)
meanpois<-40+20*sin(3*x)
y<-rpois(n,lambda=meanpois)
startv<-c(log(mean(y)),0,0,0)
modREFINED2<-(pois_p1(xseq,x,y,kP2n250[ind2,2],startv))[[1]]

# Representative quantile 0.95
qq[3]
ind3<-which(iP2n250[,2]>0.0003305&iP2n250[,2]<0.000331)
b<-ind3
set.seed(8887*b)
x<-rcircularuniform(n)
meanpois<-40+20*sin(3*x)
y<-rpois(n,lambda=meanpois)
startv<-c(log(mean(y)),0,0,0)
modREFINED3<-(pois_p1(xseq,x,y,kP2n250[ind3,2],startv))[[1]]

# radial representation
radial.plot(as.numeric(m),circular(xseq),rp.type = "l",lwd=3,radial.lim=c(2.9,4.3),labels=c(expression(0),expression(pi/4),expression(pi/2),expression(3*pi/4),expression(pi),expression(5*pi/4),expression(3*pi/2),expression(7*pi/4)),label.pos=c(0,pi/4,pi/2,3*pi/4,pi,5*pi/4,3*pi/2,7*pi/4))
radial.plot(modREFINED1,xseq,rp.type="l",add=T,lwd=3,lty=2,line.col=2,radial.lim=c(2.9,4.3))
radial.plot(modREFINED2,xseq,rp.type="l",add=T,lwd=3,lty=3,line.col=4,radial.lim=c(2.9,4.3))
radial.plot(modREFINED3,xseq,rp.type="l",add=T,lwd=3,lty=4,line.col=3,radial.lim=c(2.9,4.3))



################
### FIGURE 2 (g)
################

iG1n250<-MISE_G1_n250
kG1n250<-kappas_G1_n250

qq<-quantile(iG1n250[,2],prob=c(0.05,0.5,0.95))

# Representative quantile 0.05
qq[1]
ind<-which(iG1n250[,2]>0.01278&iG1n250[,2]<0.0128)
xseq<-seq(0,2*pi,length=300)
m<-log(4+4*sin(2*xseq)*cos(xseq))
n<-250
b<-ind
nu<-1
set.seed(888*b)
x<-rcircularuniform(n)
meangamma <- 4+4*sin(2*x)*cos(x)
y<-rgamma(n,shape=nu,rate=nu/meangamma)
startv<-c(log(mean(y)),0,0,0)
modREFINED1<-(gamma_p1(xseq,x,y,kG1n250[ind,2],startv))[[1]]

# Representative quantile 0.5
qq[2]
ind2<-which(iG1n250[,2]> 0.0263&iG1n250[,2]< 0.0264)
b<-ind2
nu<-1
set.seed(888*b)
x<-rcircularuniform(n)
meangamma <- 4+4*sin(2*x)*cos(x)
y<-rgamma(n,shape=nu,rate=nu/meangamma)
startv<-c(log(mean(y)),0,0,0)
modREFINED2<-(gamma_p1(xseq,x,y,kG1n250[ind2,2],startv))[[1]]

# Representative quantile 0.95
qq[3]
ind3<-which(iG1n250[,2]>0.04665&iG1n250[,2]<0.04675)
b<-ind3
nu<-1
set.seed(888*b)
x<-rcircularuniform(n)
meangamma <- 4+4*sin(2*x)*cos(x)
y<-rgamma(n,shape=nu,rate=nu/meangamma)
startv<-c(log(mean(y)),0,0,0)
modREFINED3<-(gamma_p1(xseq,x,y,kG1n250[ind3,2],startv))[[1]]

# radial representation
radial.plot(as.numeric(m),circular(xseq),rp.type = "l",lwd=3,radial.lim=c(-0.5,2.3),labels=c(expression(0),expression(pi/4),expression(pi/2),expression(3*pi/4),expression(pi),expression(5*pi/4),expression(3*pi/2),expression(7*pi/4)),label.pos=c(0,pi/4,pi/2,3*pi/4,pi,5*pi/4,3*pi/2,7*pi/4))
radial.plot(modREFINED1,xseq,rp.type="l",add=T,lwd=3,lty=2,line.col=2,radial.lim=c(-0.5,2.3))
radial.plot(modREFINED2,xseq,rp.type="l",add=T,lwd=3,lty=3,line.col=4,radial.lim=c(-0.5,2.3))
radial.plot(modREFINED3,xseq,rp.type="l",add=T,lwd=3,lty=4,line.col=3,radial.lim=c(-0.5,2.3))


################
### FIGURE 2 (h)
################

iG2n250<-MISE_G2_n250
kG2n250<-kappas_G2_n250

qq<-quantile(iG2n250[,2],prob=c(0.05,0.5,0.95))

# Representative quantile 0.05
qq[1]
ind<-which(iG2n250[,2]>0.00208&iG2n250[,2]<0.00209)
xseq<-seq(0,2*pi,length=300)
m<-log(5+2*cos(3+1.5*sin(xseq)))
n<-250
nu<-2
b<-ind
set.seed(888*b)
x<-rcircularuniform(n)
meangamma <- 5+2*cos(3+1.5*sin(x))
y<-rgamma(n,shape=nu,rate=nu/meangamma)
startv<-c(log(mean(y)),0,0,0)
modREFINED1<-(gamma_p1(xseq,x,y,kG2n250[ind,2],startv))[[1]]

# Representative quantile 0.5
qq[2]
ind2<-which(iG2n250[,2]>0.006035&iG2n250[,2]<0.00605)
b<-ind2
set.seed(888*b)
x<-rcircularuniform(n)
meangamma <- 5+2*cos(3+1.5*sin(x))
y<-rgamma(n,shape=nu,rate=nu/meangamma)
startv<-c(log(mean(y)),0,0,0)
modREFINED2<-(gamma_p1(xseq,x,y,kG2n250[ind2,2],startv))[[1]]

# Representative quantile 0.95
qq[3]
ind3<-which(iG2n250[,2]>0.01363&iG2n250[,2]<0.01365)
b<-ind3
set.seed(888*b)
x<-rcircularuniform(n)
meangamma <- 5+2*cos(3+1.5*sin(x))
y<-rgamma(n,shape=nu,rate=nu/meangamma)
startv<-c(log(mean(y)),0,0,0)
modREFINED3<-(gamma_p1(xseq,x,y,kG2n250[ind3,2],startv))[[1]]

# radial representation
radial.plot(as.numeric(m),circular(xseq),rp.type = "l",lwd=3,radial.lim=c(1,2),labels=c(expression(0),expression(pi/4),expression(pi/2),expression(3*pi/4),expression(pi),expression(5*pi/4),expression(3*pi/2),expression(7*pi/4)),label.pos=c(0,pi/4,pi/2,3*pi/4,pi,5*pi/4,3*pi/2,7*pi/4))
radial.plot(modREFINED1,xseq,rp.type="l",add=T,lwd=3,lty=2,line.col=2,radial.lim=c(1,2))
radial.plot(modREFINED2,xseq,rp.type="l",add=T,lwd=3,lty=3,line.col=4,radial.lim=c(1,2))
radial.plot(modREFINED3,xseq,rp.type="l",add=T,lwd=3,lty=4,line.col=3,radial.lim=c(1,2))


####################################################################
####################################################################
#                                                                  #
#  FIGURE S1      (Supplementary Material)                         #
#                                                                  #
####################################################################
####################################################################

################
### FIGURE S1 (a)
################

kN1n70<-kappas_N1_n70

xseq<-seq(0,2*pi,length=1000)
m1<- -5*cos(xseq)*sin(2*xseq)-4*sin(xseq)*cos(2*xseq) 
kappa_opt0<-(2*sqrt(pi)*(mean(m1^2))*70/(2*pi*0.35^2))^(2/5)
kappa_opt0

plot(density(kN1n70[,2],bw="SJ"),lwd=3,col=2,xlim=c(0,110),ylim=c(0,0.16),main="",xlab=expression(kappa))
points(density(kN1n70[,3],bw="SJ"),lwd=3,col=4,type="l",lty=2)
points(density(kN1n70[,1],bw="SJ"),lwd=3,col="forestgreen",type="l",lty=3)
abline(v=kappa_opt0,lwd=2,col="darkgrey")


################
### FIGURE S1 (b)
################

kN1n100<-kappas_N1_n100

xseq<-seq(0,2*pi,length=1000)
m1<- -5*cos(xseq)*sin(2*xseq)-4*sin(xseq)*cos(2*xseq) 
kappa_opt1<-(2*sqrt(pi)*(mean(m1^2))*100/(2*pi*0.35^2))^(2/5)
kappa_opt1

plot(density(kN1n100[,2],bw="SJ"),lwd=3,col=2,xlim=c(0,110),ylim=c(0,0.16),main="",xlab=expression(kappa))
points(density(kN1n100[,3],bw="SJ"),lwd=3,col=4,type="l",lty=2)
points(density(kN1n100[,1],bw="SJ"),lwd=3,col="forestgreen",type="l",lty=3)
abline(v=kappa_opt1,lwd=2,col="darkgrey")

################
### FIGURE S1 (c)
################

kN1n250<-kappas_N1_n250

xseq<-seq(0,2*pi,length=1000)
m1<- -5*cos(xseq)*sin(2*xseq)-4*sin(xseq)*cos(2*xseq) 
kappa_opt2<-(2*sqrt(pi)*mean(m1^2)*250/(2*pi*0.35^2))^(2/5)
kappa_opt2

plot(density(kN1n250[,2],bw="SJ"),lwd=3,col=2,xlim=c(0,110),ylim=c(0,0.16),main="",xlab=expression(kappa))
points(density(kN1n250[,3],bw="SJ"),lwd=3,col=4,type="l",lty=2)
points(density(kN1n250[,1],bw="SJ"),lwd=3,col="forestgreen",type="l",lty=3)
abline(v=kappa_opt2,lwd=2,col="darkgrey")


################
### FIGURE S1 (d)
################

kN1n500<-kappas_N1_n500

xseq<-seq(0,2*pi,length=1000)
m1<- -5*cos(xseq)*sin(2*xseq)-4*sin(xseq)*cos(2*xseq) 
kappa_opt3<-(2*sqrt(pi)*mean(m1^2)*500/(2*pi*0.35^2))^(2/5)
kappa_opt3

plot(density(kN1n500[,2],bw="SJ"),lwd=3,col=2,xlim=c(0,110),ylim=c(0,0.16),main="",xlab=expression(kappa))
points(density(kN1n500[,3],bw="SJ"),lwd=3,col=4,type="l",lty=2)
points(density(kN1n500[,1],bw="SJ"),lwd=3,col="forestgreen",type="l",lty=3)
abline(v=kappa_opt3,lwd=2,col="darkgrey")


################
### FIGURE S1 (e)
################

kN1n1500<-kappas_N1_n1500

xseq<-seq(0,2*pi,length=1000)
m1<- -5*cos(xseq)*sin(2*xseq)-4*sin(xseq)*cos(2*xseq) 
kappa_opt4<-(2*sqrt(pi)*mean(m1^2)*1500/(2*pi*0.35^2))^(2/5)
kappa_opt4

plot(density(kN1n1500[,2],bw="SJ"),lwd=3,col=2,xlim=c(0,110),ylim=c(0,0.16),main="",xlab=expression(kappa))
points(density(kN1n1500[,3],bw="SJ"),lwd=3,col=4,type="l",lty=2)
points(density(kN1n1500[,1],bw="SJ"),lwd=3,col="forestgreen",type="l",lty=3)
abline(v=kappa_opt4,lwd=2,col="darkgrey")


####################################################################
####################################################################
#                                                                  #
#  FIGURE S2      (Supplementary Material)                         #
#                                                                  #
####################################################################
####################################################################


################
### FIGURE S2 (a)
################

kN2n70<-kappas_N2_n70

m2<- 3.5*sin(xseq)-cos(xseq) 
kappa_opt0<-(2*sqrt(pi)*mean(m2^2)*70/(2*pi*0.5^2))^(2/5)
kappa_opt0

plot(density(kN2n70[,2],bw="SJ"),lwd=3,col=2,xlim=c(0,80),ylim=c(0,0.3),main="",xlab=expression(kappa))
points(density(kN2n70[,3],bw="SJ"),lwd=3,col=4,type="l",lty=2)
points(density(kN2n70[,1],bw="SJ"),lwd=3,col="forestgreen",type="l",lty=3)
abline(v=kappa_opt1,lwd=2,col="darkgrey")


################
### FIGURE S2 (b)
################

kN2n100<-kappas_N2_n100

m2<- 3.5*sin(xseq)-cos(xseq) 
kappa_opt1<-(2*sqrt(pi)*mean(m2^2)*100/(2*pi*0.5^2))^(2/5)
kappa_opt1

plot(density(kN2n100[,2],bw="SJ"),lwd=3,col=2,xlim=c(0,80),ylim=c(0,0.3),main="",xlab=expression(kappa))
points(density(kN2n100[,3],bw="SJ"),lwd=3,col=4,type="l",lty=2)
points(density(kN2n100[,1],bw="SJ"),lwd=3,col="forestgreen",type="l",lty=3)
abline(v=kappa_opt1,lwd=2,col="darkgrey")


################
### FIGURE S2 (c)
################

kN2n250<-kappas_N2_n250

m2<- 3.5*sin(xseq)-cos(xseq) 
kappa_opt2<-(2*sqrt(pi)*mean(m2^2)*250/(2*pi*0.5^2))^(2/5)
kappa_opt2

plot(density(kN2n250[,2],bw="SJ"),lwd=3,col=2,xlim=c(0,80),ylim=c(0,0.3),main="",xlab=expression(kappa))
points(density(kN2n250[,3],bw="SJ"),lwd=3,col=4,type="l",lty=2)
points(density(kN2n250[,1],bw="SJ"),lwd=3,col="forestgreen",type="l",lty=3)
abline(v=kappa_opt2,lwd=2,col="darkgrey")


################
### FIGURE S2 (d)
################

kN2n500<-kappas_N2_n500

m2<- 3.5*sin(xseq)-cos(xseq) 
kappa_opt3<-(2*sqrt(pi)*mean(m2^2)*500/(2*pi*0.5^2))^(2/5)
kappa_opt3

plot(density(kN2n500[,2],bw="SJ"),lwd=3,col=2,xlim=c(0,80),ylim=c(0,0.3),main="",xlab=expression(kappa))
points(density(kN2n500[,3],bw="SJ"),lwd=3,col=4,type="l",lty=2)
points(density(kN2n500[,1],bw="SJ"),lwd=3,col="forestgreen",type="l",lty=3)
abline(v=kappa_opt3,lwd=2,col="darkgrey")


################
### FIGURE S2 (e)
################

kN2n1500<-kappas_N2_n1500

m2<- 3.5*sin(xseq)-cos(xseq) 

kappa_opt4<-(2*sqrt(pi)*mean(m2^2)*1500/(2*pi*0.5^2))^(2/5)
kappa_opt4

plot(density(kN2n1500[,2],bw="SJ"),lwd=3,col=2,xlim=c(0,80),ylim=c(0,0.3),main="",xlab=expression(kappa))
points(density(kN2n1500[,3],bw="SJ"),lwd=3,col=4,type="l",lty=2)
points(density(kN2n1500[,1],bw="SJ"),lwd=3,col="forestgreen",type="l",lty=3)
abline(v=kappa_opt4,lwd=2,col="darkgrey")


####################################################################
####################################################################
#                                                                  #
#  FIGURE S3      (Supplementary Material)                         #
#                                                                  #
####################################################################
####################################################################


################
### FIGURE S3 (a)
################

iN1n70<-MISE_N1_n70

dfs70 <- stack(iN1n70)
boxplot(values~ind,data=dfs70, main="",
        xlab="", ylab="",col="lightblue",names=c("CRSC","Refined","CV"),ylim=c(0,0.3))


################
### FIGURE S3 (b)
################

iN1n100<-MISE_N1_n100

dfs100 <- stack(iN1n100)
boxplot(values~ind,data=dfs100, main="",
        xlab="", ylab="",col="lightblue",names=c("CRSC","Refined","CV"),ylim=c(0,0.3))


################
### FIGURE S3 (c)
################

iN1n250<-MISE_N1_n250

dfs250 <- stack(iN1n250)
boxplot(values~ind,data=dfs250, main="",
        xlab="", ylab="",col="lightblue",names=c("CRSC","Refined","CV"))


################
### FIGURE S3 (d)
################

iN1n500<-MISE_N1_n500

dfs500 <- stack(iN1n500)
boxplot(values~ind,data=dfs500, main="",
        xlab="", ylab="",col="lightblue",names=c("CRSC","Refined","CV"))

################
### FIGURE S3 (e)
################

iN1n1500<-MISE_N1_n1500

dfs1500 <- stack(iN1n1500)
boxplot(values~ind,data=dfs1500, main="",
        xlab="", ylab="",col="lightblue",names=c("CRSC","Refined","CV"))


####################################################################
####################################################################
#                                                                  #
#  FIGURE S4      (Supplementary Material)                         #
#                                                                  #
####################################################################
####################################################################

################
### FIGURE S4 (a)
################

iN2n70<-MISE_N2_n70

dfs70 <- stack(iN2n70)
boxplot(values~ind,data=dfs70, main="",
        xlab="", ylab="",col="lightblue",names=c("CRSC","Refined","CV"))


################
### FIGURE S4 (b)
################

iN2n100<-MISE_N2_n100

dfs100 <- stack(iN2n100)
boxplot(values~ind,data=dfs100, main="",
        xlab="", ylab="",col="lightblue",names=c("CRSC","Refined","CV"))


################
### FIGURE S4 (c)
################

iN2n250<-MISE_N2_n250

dfs250 <- stack(iN2n250)
boxplot(values~ind,data=dfs250, main="",
        xlab="", ylab="",col="lightblue",names=c("CRSC","Refined","CV"))


################
### FIGURE S4 (d)
################

iN2n500<-MISE_N2_n500

dfs500 <- stack(iN2n500)
boxplot(values~ind,data=dfs500, main="",
        xlab="", ylab="",col="lightblue",names=c("CRSC","Refined","CV"))


################
### FIGURE S4 (e)
################

iN2n1500<-MISE_N2_n1500

dfs1500 <- stack(iN2n1500)
boxplot(values~ind,data=dfs1500, main="",
        xlab="", ylab="",col="lightblue",names=c("CRSC","Refined","CV"))


####################################################################
####################################################################
#                                                                  #
#  FIGURE S5      (Supplementary Material)                         #
#                                                                  #
####################################################################
####################################################################

################
### FIGURE S5 (a)
################

kB1n70<-kappas_B1_n70

plot(density(kB1n70[-which(is.na(kB1n70[,2])),2]),lwd=3,col=2,xlim=c(0,60),ylim=c(0,0.22),main="",xlab=expression(kappa))
points(density(kB1n70[,3]),lwd=3,col=4,type="l",lty=2)
points(density(kB1n70[,1]),lwd=3,col="forestgreen",type="l",lty=3)

################
### FIGURE S5 (b)
################

kB1n100<-kappas_B1_n100

plot(density(kB1n100[-which(is.na(kB1n100[,2])),2]),lwd=3,col=2,xlim=c(0,60),ylim=c(0,0.22),main="",xlab=expression(kappa))
points(density(kB1n100[,3]),lwd=3,col=4,type="l",lty=2)
points(density(kB1n100[,1]),lwd=3,col="forestgreen",type="l",lty=3)

################
### FIGURE S5 (c)
################

kB1n250<-kappas_B1_n250

plot(density(kB1n250[,2],bw="SJ"),lwd=3,col=2,xlim=c(0,60),ylim=c(0,0.22),main="",xlab=expression(kappa))
points(density(kB1n250[,3],bw="SJ"),lwd=3,col=4,type="l",lty=2)
points(density(kB1n250[,1],bw="SJ"),lwd=3,col="forestgreen",type="l",lty=3)

################
### FIGURE S5 (d)
################

kB1n500<-kappas_B1_n500

plot(density(kB1n500[,2],bw="SJ"),lwd=3,col=2,xlim=c(0,60),ylim=c(0,0.22),main="",xlab=expression(kappa))
points(density(kB1n500[,3],bw="SJ"),lwd=3,col=4,type="l",lty=2)
points(density(kB1n500[,1],bw="SJ"),lwd=3,col="forestgreen",type="l",lty=3)

################
### FIGURE S5 (e)
################

kB1n1500<-kappas_B1_n1500

plot(density(kB1n1500[,2]),lwd=3,col=2,xlim=c(0,60),ylim=c(0,0.22),main="",xlab=expression(kappa))
points(density(kB1n1500[,3]),lwd=3,col=4,type="l",lty=2)
points(density(kB1n1500[,1]),lwd=3,col="forestgreen",type="l",lty=3)


####################################################################
####################################################################
#                                                                  #
#  FIGURE S6      (Supplementary Material)                         #
#                                                                  #
####################################################################
####################################################################

################
### FIGURE S6 (a)
################

kB2n70<-kappas_B2_n70

plot(density(kB2n70[-which(is.na(kB2n70[,2])),2]),lwd=3,col=2,xlim=c(0,20),ylim=c(0,0.55),main="",xlab=expression(kappa))
points(density(kB2n70[,3]),lwd=3,col=4,type="l",lty=2)
points(density(kB2n70[,1]),lwd=3,col="forestgreen",type="l",lty=3)

################
### FIGURE S6 (b)
################

kB2n100<-kappas_B2_n100

plot(density(kB2n100[-which(is.na(kB2n100[,2])),2]),lwd=3,col=2,xlim=c(0,20),ylim=c(0,0.55),main="",xlab=expression(kappa))
points(density(kB2n100[,3]),lwd=3,col=4,type="l",lty=2)
points(density(kB2n100[,1]),lwd=3,col="forestgreen",type="l",lty=3)

################
### FIGURE S6 (c)
################

kB2n250<-kappas_B2_n250

plot(density(kB2n250[,2]),lwd=3,col=2,xlim=c(0,20),main="",xlab=expression(kappa))
points(density(kB2n250[,3]),lwd=3,col=4,type="l",lty=2)
points(density(kB2n250[,1]),lwd=3,col="forestgreen",type="l",lty=3)

################
### FIGURE S6 (d)
################

kB2n500<-kappas_B2_n500

plot(density(kB2n500[,2]),lwd=3,col=2,xlim=c(0,20),main="",xlab=expression(kappa))
points(density(kB2n500[,3]),lwd=3,col=4,type="l",lty=2)
points(density(kB2n500[,1]),lwd=3,col="forestgreen",type="l",lty=3)

################
### FIGURE S6 (e)
################

kB2n1500<-kappas_B2_n1500

plot(density(kB2n1500[,2],bw=0.45),lwd=3,col=2,xlim=c(0,20),main="",xlab=expression(kappa))
points(density(kB2n1500[,3]),lwd=3,col=4,type="l",lty=2)
points(density(kB2n1500[,1]),lwd=3,col="forestgreen",type="l",lty=3)


####################################################################
####################################################################
#                                                                  #
#  FIGURE S7      (Supplementary Material)                         #
#                                                                  #
####################################################################
####################################################################

################
### FIGURE S7 (a)
################

iB1n70<-MISE_B1_n70

dfs70 <- stack(iB1n70)
boxplot(values~ind,data=dfs70, main="",
        xlab="", ylab="",col="lightblue",names=c("ECRSC","Refined","CV"))

################
### FIGURE S7 (b)
################

iB1n100<-MISE_B1_n100

dfs100 <- stack(iB1n100)
boxplot(values~ind,data=dfs100, main="",
        xlab="", ylab="",col="lightblue",names=c("ECRSC","Refined","CV"))

################
### FIGURE S7 (c)
################

iB1n250<-MISE_B1_n250

dfs250 <- stack(iB1n250)
boxplot(values~ind,data=dfs250, main="",
        xlab="", ylab="",col="lightblue",names=c("ECRSC","Refined","CV"))

################
### FIGURE S7 (d)
################

iB1n500<-MISE_B1_n500

dfs500 <- stack(iB1n500)
boxplot(values~ind,data=dfs500, main="",
        xlab="", ylab="",col="lightblue",names=c("ECRSC","Refined","CV"))

################
### FIGURE S7 (e)
################

iB1n1500<-MISE_B1_n1500

dfs1500 <- stack(iB1n1500)
boxplot(values~ind,data=dfs1500, main="",
        xlab="", ylab="",col="lightblue",names=c("ECRSC","Refined","CV"))


####################################################################
####################################################################
#                                                                  #
#  FIGURE S8      (Supplementary Material)                         #
#                                                                  #
####################################################################
####################################################################

################
### FIGURE S8 (a)
################

iB2n70<-MISE_B2_n70

dfs70 <- stack(iB2n70)
boxplot(values~ind,data=dfs70, main="",
        xlab="", ylab="",col="lightblue",names=c("ECRSC","Refined","CV"))

################
### FIGURE S8 (b)
################

iB2n100<-MISE_B2_n100

dfs100 <- stack(iB2n100)
boxplot(values~ind,data=dfs100, main="",
        xlab="", ylab="",col="lightblue",names=c("ECRSC","Refined","CV"))

################
### FIGURE S8 (c)
################

iB2n250<-MISE_B2_n250

dfs250 <- stack(iB2n250)
boxplot(values~ind,data=dfs250, main="",
        xlab="", ylab="",col="lightblue",names=c("ECRSC","Refined","CV"))

################
### FIGURE S8 (d)
################

iB2n500<-MISE_B2_n500

dfs500 <- stack(iB2n500)
boxplot(values~ind,data=dfs500, main="",
        xlab="", ylab="",col="lightblue",names=c("ECRSC","Refined","CV"))

################
### FIGURE S8 (e)
################

iB2n1500<-MISE_B2_n1500

dfs1500 <- stack(iB2n1500)
boxplot(values~ind,data=dfs1500, main="",
        xlab="", ylab="",col="lightblue",names=c("ECRSC","Refined","CV"))



####################################################################
####################################################################
#                                                                  #
#  FIGURE S9      (Supplementary Material)                         #
#                                                                  #
####################################################################
####################################################################

################
### FIGURE S9 (a)
################

kP1n70<-kappas_P1_n70

plot(density(kP1n70[,2],bw="SJ"),lwd=3,col=2,xlim=c(0,55),ylim=c(0,0.25),main="",xlab=expression(kappa))
points(density(kP1n70[,3],bw="SJ"),lwd=3,col=4,type="l",lty=2)
points(density(kP1n70[,1],bw="SJ"),lwd=3,col="forestgreen",type="l",lty=3)

################
### FIGURE S9 (b)
################

kP1n100<-kappas_P1_n100

plot(density(kP1n100[,2],bw="SJ"),lwd=3,col=2,xlim=c(0,55),ylim=c(0,0.25),main="",xlab=expression(kappa))
points(density(kP1n100[,3],bw="SJ"),lwd=3,col=4,type="l",lty=2)
points(density(kP1n100[,1],bw="SJ"),lwd=3,col="forestgreen",type="l",lty=3)

################
### FIGURE S9 (c)
################

kP1n250<-kappas_P1_n250

plot(density(kP1n250[,2],bw="SJ"),lwd=3,col=2,xlim=c(0,55),ylim=c(0,0.25),main="",xlab=expression(kappa))
points(density(kP1n250[,3],bw="SJ"),lwd=3,col=4,type="l",lty=2)
points(density(kP1n250[,1],bw="SJ"),lwd=3,col="forestgreen",type="l",lty=3)

################
### FIGURE S9 (d)
################

kP1n500<-kappas_P1_n500

plot(density(kP1n500[,2],bw="SJ"),lwd=3,col=2,xlim=c(0,55),ylim=c(0,0.25),main="",xlab=expression(kappa))
points(density(kP1n500[,3],bw="SJ"),lwd=3,col=4,type="l",lty=2)
points(density(kP1n500[,1],bw="SJ"),lwd=3,col="forestgreen",type="l",lty=3)

################
### FIGURE S9 (d)
################

kP1n1500<-kappas_P1_n1500

plot(density(kP1n1500[,2],bw="SJ"),lwd=3,col=2,xlim=c(0,55),ylim=c(0,0.25),main="",xlab=expression(kappa))
points(density(kP1n1500[,3],bw="SJ"),lwd=3,col=4,type="l",lty=2)
points(density(kP1n1500[,1],bw="SJ"),lwd=3,col="forestgreen",type="l",lty=3)


####################################################################
####################################################################
#                                                                  #
#  FIGURE S10      (Supplementary Material)                        #
#                                                                  #
####################################################################
####################################################################

#################
### FIGURE S10 (a)
#################

kP2n70<-kappas_P2_n70

plot(density(kP2n70[,2],bw="SJ"),lwd=3,col=2,xlim=c(0,200),ylim=c(0,0.1),main="",xlab=expression(kappa))
points(density(kP2n70[,3],bw="SJ"),lwd=3,col=4,type="l",lty=2)
points(density(kP2n70[,1],bw="SJ"),lwd=3,col="forestgreen",type="l",lty=3)

#################
### FIGURE S10 (b)
#################

kP2n100<-kappas_P2_n100

plot(density(kP2n100[,2],bw="SJ"),lwd=3,col=2,xlim=c(0,200),ylim=c(0,0.1),main="",xlab=expression(kappa))
points(density(kP2n100[,3],bw="SJ"),lwd=3,col=4,type="l",lty=2)
points(density(kP2n100[,1],bw="SJ"),lwd=3,col="forestgreen",type="l",lty=3)

#################
### FIGURE S10 (c)
#################

kP2n250<-kappas_P2_n250

plot(density(kP2n250[,2],bw="SJ"),lwd=3,col=2,xlim=c(0,200),ylim=c(0,0.1),main="",xlab=expression(kappa))
points(density(kP2n250[,3],bw="SJ"),lwd=3,col=4,type="l",lty=2)
points(density(kP2n250[,1],bw="SJ"),lwd=3,col="forestgreen",type="l",lty=3)

#################
### FIGURE S10 (d)
#################

kP2n500<-kappas_P2_n500

plot(density(kP2n500[,2],bw="SJ"),lwd=3,col=2,xlim=c(0,200),ylim=c(0,0.1),main="",xlab=expression(kappa))
points(density(kP2n500[,3],bw="SJ"),lwd=3,col=4,type="l",lty=2)
points(density(kP2n500[,1],bw="SJ"),lwd=3,col="forestgreen",type="l",lty=3)

#################
### FIGURE S10 (e)
#################

kP2n1500<-kappas_P2_n1500

plot(density(kP2n1500[,2],bw=4),lwd=3,col=2,xlim=c(0,200),ylim=c(0,0.1),main="",xlab=expression(kappa))
points(density(kP2n1500[,3],bw="SJ"),lwd=3,col=4,type="l",lty=2)
points(density(kP2n1500[,1],bw="SJ"),lwd=3,col="forestgreen",type="l",lty=3)


####################################################################
####################################################################
#                                                                  #
#  FIGURE S11      (Supplementary Material)                        #
#                                                                  #
####################################################################
####################################################################

#################
### FIGURE S11 (a)
#################

iP1n70<-MISE_P1_n70

dfs70 <- stack(iP1n70)
boxplot(values~ind,data=dfs70, main="",
        xlab="", ylab="",col="lightblue",names=c("CRSC","Refined","CV"),ylim=c(0,0.025))

#################
### FIGURE S11 (b)
#################

iP1n100<-MISE_P1_n100

dfs100 <- stack(iP1n100)
boxplot(values~ind,data=dfs100, main="",
        xlab="", ylab="",col="lightblue",names=c("CRSC","Refined","CV"),ylim=c(0,0.025))

#################
### FIGURE S11 (c)
#################

iP1n250<-MISE_P1_n250

dfs250 <- stack(iP1n250)
boxplot(values~ind,data=dfs250, main="",
        xlab="", ylab="",col="lightblue",names=c("CRSC","Refined","CV"))

#################
### FIGURE S11 (d)
#################

iP1n500<-MISE_P1_n500

dfs500 <- stack(iP1n500)
boxplot(values~ind,data=dfs500, main="",
        xlab="", ylab="",col="lightblue",names=c("CRSC","Refined","CV"))

#################
### FIGURE S11 (e)
#################

iP1n15000<-MISE_P1_n15000

dfs1500 <- stack(iP1n1500)
boxplot(values~ind,data=dfs1500, main="",
        xlab="", ylab="",col="lightblue",names=c("ECRSC","Refined","CV"))


####################################################################
####################################################################
#                                                                  #
#  FIGURE S12      (Supplementary Material)                        #
#                                                                  #
####################################################################
####################################################################

#################
### FIGURE S12 (a)
#################

iP2n70<-MISE_P2_n70

dfs70 <- stack(iP2n70)
boxplot(values~ind,data=dfs100, main="",
        xlab="", ylab="",col="lightblue",names=c("CRSC","Refined","CV"))

#################
### FIGURE S12 (b)
#################

iP2n100<-MISE_P2_n100

dfs100 <- stack(iP2n100)
boxplot(values~ind,data=dfs100, main="",
        xlab="", ylab="",col="lightblue",names=c("CRSC","Refined","CV"))

#################
### FIGURE S12 (c)
#################

iP2n250<-MISE_P2_n250

dfs250 <- stack(iP2n250)
boxplot(values~ind,data=dfs250, main="",
        xlab="", ylab="",col="lightblue",names=c("CRSC","Refined","CV"))

#################
### FIGURE S12 (d)
#################

iP2n500<-MISE_P2_n500

dfs500 <- stack(iP2n500)
boxplot(values~ind,data=dfs500, main="",
        xlab="", ylab="",col="lightblue",names=c("CRSC","Refined","CV"))

#################
### FIGURE S12 (e)
#################

iP2n1500<-MISE_P2_n1500

dfs1500 <- stack(iP2n1500)
boxplot(values~ind,data=dfs1500, main="",
        xlab="", ylab="",col="lightblue",names=c("ECRSC","Refined","CV"))


####################################################################
####################################################################
#                                                                  #
#  FIGURE S13      (Supplementary Material)                        #
#                                                                  #
####################################################################
####################################################################

#################
### FIGURE S13 (a)
#################

kG1n70<-kappas_G1_n70

plot(density(kG1n70[-which(is.na(kG1n70[,2])),2],bw="SJ"),lwd=3,col=2,xlim=c(0,130),ylim=c(0,0.11),main="",xlab=expression(kappa))
points(density(kG1n70[,3],bw="SJ"),lwd=3,col=4,type="l",lty=2)
points(density(kG1n70[-which(is.na(kG1n70[,1])),1],bw="SJ"),lwd=3,col="forestgreen",type="l",lty=3)


#################
### FIGURE S13 (b)
#################

kG1n100<-kappas_G1_n100

plot(density(kG1n100[,2],bw="SJ"),lwd=3,col=2,xlim=c(0,130),ylim=c(0,0.11),main="",xlab=expression(kappa))
points(density(kG1n100[,3],bw="SJ"),lwd=3,col=4,type="l",lty=2)
points(density(kG1n100[,1],bw="SJ"),lwd=3,col="forestgreen",type="l",lty=3)

#################
### FIGURE S13 (c)
#################

kG1n250<-kappas_G1_n250

plot(density(kG1n250[,2],bw="SJ"),lwd=3,col=2,xlim=c(0,130),ylim=c(0,0.11),main="",xlab=expression(kappa))
points(density(kG1n250[,3],bw="SJ"),lwd=3,col=4,type="l",lty=2)
points(density(kG1n250[,1],bw="SJ"),lwd=3,col="forestgreen",type="l",lty=3)


#################
### FIGURE S13 (d)
#################

kG1n500<-kappas_G1_n500

plot(density(kG1n500[,2],bw="SJ"),lwd=3,col=2,xlim=c(0,130),ylim=c(0,0.11),main="",xlab=expression(kappa))
points(density(kG1n500[,3],bw="SJ"),lwd=3,col=4,type="l",lty=2)
points(density(kG1n500[,1],bw="SJ"),lwd=3,col="forestgreen",type="l",lty=3)


#################
### FIGURE S13 (e)
#################

kG1n1500<-kappas_G1_n1500

plot(density(kG1n1500[,2],bw="SJ"),lwd=3,col=2,xlim=c(0,130),ylim=c(0,0.11),main="",xlab=expression(kappa))
points(density(kG1n1500[,3],bw="SJ"),lwd=3,col=4,type="l",lty=2)
points(density(kG1n1500[,1],bw="SJ"),lwd=3,col="forestgreen",type="l",lty=3)


####################################################################
####################################################################
#                                                                  #
#  FIGURE S14      (Supplementary Material)                        #
#                                                                  #
####################################################################
####################################################################

#################
### FIGURE S14 (a)
#################

kG2n70<-kappas_G2_n70

plot(density(kG2n70[,2],bw="SJ"),lwd=3,col=2,xlim=c(0,40),ylim=c(0,0.35),main="",xlab=expression(kappa))
points(density(kG2n70[,3],bw="SJ"),lwd=3,col=4,type="l",lty=2)
points(density(kG2n70[,1],bw="SJ"),lwd=3,col="forestgreen",type="l",lty=3)

#################
### FIGURE S14 (b)
#################

kG2n100<-kappas_G2_n100

plot(density(kG2n100[,2],bw="SJ"),lwd=3,col=2,xlim=c(0,40),ylim=c(0,0.35),main="",xlab=expression(kappa))
points(density(kG2n100[,3],bw="SJ"),lwd=3,col=4,type="l",lty=2)
points(density(kG2n100[,1],bw="SJ"),lwd=3,col="forestgreen",type="l",lty=3)

#################
### FIGURE S14 (c)
#################

kG2n250<-kappas_G2_n250

plot(density(kG2n250[,2],bw="SJ"),lwd=3,col=2,xlim=c(0,40),ylim=c(0,0.35),main="",xlab=expression(kappa))
points(density(kG2n250[,3],bw="SJ"),lwd=3,col=4,type="l",lty=2)
points(density(kG2n250[,1],bw="SJ"),lwd=3,col="forestgreen",type="l",lty=3)

#################
### FIGURE S14 (d)
#################

kG2n500<-kappas_G2_n500

plot(density(kG2n500[,2],bw="SJ"),lwd=3,col=2,xlim=c(0,40),ylim=c(0,0.35),main="",xlab=expression(kappa))
points(density(kG2n500[,3],bw="SJ"),lwd=3,col=4,type="l",lty=2)
points(density(kG2n500[,1],bw="SJ"),lwd=3,col="forestgreen",type="l",lty=3)

#################
### FIGURE S14 (e)
#################

kG2n1500<-kappas_G2_n1500

plot(density(kG2n1500[,2],bw="SJ"),lwd=3,col=2,xlim=c(0,40),ylim=c(0,0.35),main="",xlab=expression(kappa))
points(density(kG2n1500[,3],bw="SJ"),lwd=3,col=4,type="l",lty=2)
points(density(kG2n1500[,1],bw="SJ"),lwd=3,col="forestgreen",type="l",lty=3)


####################################################################
####################################################################
#                                                                  #
#  FIGURE S15      (Supplementary Material)                        #
#                                                                  #
####################################################################
####################################################################

#################
### FIGURE S15 (a)
#################

iG1n70<-MISE_G1_n70

dfs70 <- stack(iG1n70)
boxplot(values~ind,data=dfs70, main="",
        xlab="", ylab="",col="lightblue",names=c("ECRSC","Refined","CV"),ylim=c(0,0.3))

#################
### FIGURE S15 (b)
#################

iG1n100<-MISE_G1_n100

dfs100 <- stack(iG1n100)
boxplot(values~ind,data=dfs100, main="",
        xlab="", ylab="",col="lightblue",names=c("ECRSC","Refined","CV"),ylim=c(0,0.3))

#################
### FIGURE S15 (c)
#################

iG1n250<-MISE_G1_n250

dfs250 <- stack(iG1n250)
boxplot(values~ind,data=dfs250, main="",
        xlab="", ylab="",col="lightblue",names=c("ECRSC","Refined","CV"))

#################
### FIGURE S15 (d)
#################

iG1n500<-MISE_G1_n500

dfs500 <- stack(iG1n500)
boxplot(values~ind,data=dfs500, main="",
        xlab="", ylab="",col="lightblue",names=c("ECRSC","Refined","CV"))


#################
### FIGURE S15 (e)
#################

iG1n1500<-MISE_G1_n1500

dfs1500 <- stack(iG1n1500)
boxplot(values~ind,data=dfs1500, main="",
        xlab="", ylab="",col="lightblue",names=c("ECRSC","Refined","CV"),ylim=c(0,0.018))


####################################################################
####################################################################
#                                                                  #
#  FIGURE S16      (Supplementary Material)                        #
#                                                                  #
####################################################################
####################################################################

#################
### FIGURE S16 (a)
#################

iG2n70<-MISE_G2_n70

dfs70 <- stack(iG2n70)
boxplot(values~ind,data=dfs70, main="",
        xlab="", ylab="",col="lightblue",names=c("ECRSC","Refined","CV"))

#################
### FIGURE S16 (b)
#################

iG2n100<-MISE_G2_n100

dfs100 <- stack(iG2n100)
boxplot(values~ind,data=dfs100, main="",
        xlab="", ylab="",col="lightblue",names=c("ECRSC","Refined","CV"))

#################
### FIGURE S16 (c)
#################

iG2n250<-MISE_G2_n250

dfs250 <- stack(iG2n250)
boxplot(values~ind,data=dfs250, main="",
        xlab="", ylab="",col="lightblue",names=c("ECRSC","Refined","CV"))

#################
### FIGURE S16 (d)
#################

iG2n500<-MISE_G2_n500

dfs500 <- stack(iG2n500)
boxplot(values~ind,data=dfs500, main="",
        xlab="", ylab="",col="lightblue",names=c("ECRSC","Refined","CV"))


#################
### FIGURE S16 (e)
#################

iG2n1500<-MISE_G2_n1500

dfs1500 <- stack(iG2n1500)
boxplot(values~ind,data=dfs1500, main="",
        xlab="", ylab="",col="lightblue",names=c("ECRSC","Refined","CV"))
