#Below we provide example R code to run the bivariate mixed model on oak and frass data.

#The analyses below require that the MCMCglmm package has been installed.
library(MCMCglmm)

#Read in the oak and frass data frame, which is downloadadble from http://dx.doi.org/10.7488/ds/2215.

oakfrass<-read.csv("...../oakandfrass.txt",sep="\t")
#there are 10073 rows of oak data and 695 of frass data. 
#The data types differ in that oak has a single estimate of first leaf, whereas for frass we have interval censored phenology data.


#################################################################################################################
###### PRIORS

#model priors. For the G random terms we specify paramater expanded priors.
a<-1000

pa_prior<-list(R=list(V=diag(2), nu=0.002), 
G=list(G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*a),
G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*a),
G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*a),
G1=list(V=diag(1), nu=1, alpha.mu=c(0), alpha.V=diag(1)*a)))

#################################################################################################################
###### MODEL

#In the bivariate response model below we treat both oak and frass phenology as interval censored ("cenguassian"). However for oak phenology.lower = phenology.upper.
#Specifying the us term for a random effect allows us to estimate a separate variance for each species and the covariance between species across levels of that random effect.

model<-MCMCglmm(cbind(phenology.lower, phenology.upper)~species-1+species: mean.centre.latitude, 
random=
~us(species):grid_50km+
us(species):grid_5km+
us(species):year+
us(at.level(species,"frass")): tray,
rcov=~idh(species):units,family=c("cengaussian"),nitt=400000,burnin=40000,thin=100,data= oakfrass,prior=pa_prior)


#Below we check for convergence, autocorrelation and adequate effective sample sizes of focal parameter estimates.

plot(model$Sol)
plot(model$VCV)

autocorr(model$Sol)
autocorr(model$VCV)

summary(model)

#Mean phenology of oak at the latitudinal midpoint (followed by credible intervals)
mean(model$Sol[,"speciesoak"])
HPDinterval(model$Sol[,"speciesoak"])

#Mean phenology of frass at the latitudinal midpoint (followed by credible intervals)
mean(model$Sol[,"speciesfrass"])
HPDinterval(model$Sol[,"speciesfrass"])

#Difference in timing at the latitudinal midpoint (followed by credible intervals)
mean(model$Sol[,"speciesfrass"]-model$Sol[,"speciesoak"])
HPDinterval(model$Sol[,"speciesfrass"]-model$Sol[,"speciesoak"])

#Latitudinal slope (and credible interval) for oak and frass.
#oak
mean(model$Sol[,"speciesoak:mean.centre.latitude"])
HPDinterval(model$Sol[,"speciesoak:mean.centre.latitude"])
#frass
mean(model$Sol[,"speciesfrass:mean.centre.latitude"])
HPDinterval(model$Sol[,"speciesfrass:mean.centre.latitude"])


#################################################################################################################
###### ESTIMATE MAJOR AXIS SLOPE
#Code adapted from the line.cis function in the smatr library - authors: Warton, D.I., J. Ormerod, & S. Taskinen

#Temporal major axis slope
	vary <-model$VCV[,"speciesfrass:speciesfrass.year"]
	covxy<-model$VCV[,"speciesoak:speciesfrass.year"]
	varx <-model$VCV[,"speciesoak:speciesoak.year"]
	fac <- vary - varx
	major.axis <- (fac + sqrt(fac^2 + 4 * covxy^2))/2/covxy
	
mean(major.axis)
HPDinterval(major.axis)


#################################################################################################################
###### ESTIMATE CORRELATION

#Temporal correlation
posterior.correlation<-c()
for(x in 1:length(model$VCV[,1])){
	mat<-matrix(nrow=2,ncol=2,
	c(model$VCV[x,"speciesoak:speciesoak.year"],
	model$VCV[x,"speciesoak:speciesfrass.year"],
	model$VCV[x,"speciesoak:speciesfrass.year"],
	model$VCV[x,"speciesfrass:speciesfrass.year"]))
	
	posterior.correlation [x]<-cov2cor(mat)[1,2]
}

mean(posterior.correlation)
HPDinterval(mcmc(posterior.correlation))

