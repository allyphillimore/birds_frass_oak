library("MCMCglmm")

library(mvtnorm)
#required for simulating bivariate distributions.

#the paramater expanded priors we will use
a<-1000
pa_prior<-list(R=list(V=diag(2), nu=0.002), 
G=list(G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*a),
G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*a),
G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*a)))


datastructure<-read.csv("datastructure.txt",sep="\t",header=T)
#The power analysis uses the structure of real data - i.e. the same replication of years, sites and latitudes.
#The datastructure file is in the same repository on Github


#con = consumer species
#res = resource species

#you will also need to specify yourfolder

power.analysis<-
function(iteration,con_spec,res_spec,con_intercept,res_intercept,con_latslope,res_latslope,con_spacevar,conres_spacecov,res_spacevar,
con_5kmvar,conres_5kmcov,res_5kmvar,con_tempvar,conres_tempcov,res_tempvar,con_residvar,res_residvar)
{
	#Parameter estimates required.
	#Intercepts
	#Slopes of phenology on latitude
	#Spatial VCV
	#5km grid cell VCV
	#Temporal VCV
	#Residual V

	corr_space<-conres_spacecov/(sqrt(con_spacevar)*sqrt(res_spacevar))
	corr_time<-conres_tempcov/(sqrt(con_tempvar)*sqrt(res_tempvar))

	#select the structure for the focal species
	consumer<-which(as.character(datastructure$species)==con_spec)
	resource<-which(as.character(datastructure$species)==res_spec)
	speciesdata<-datastructure[c(consumer,resource),]
	speciesdata$spec<-as.factor(c(rep("con",length(consumer)),rep("res",length(resource))))
	speciesdata$year<-as.factor(speciesdata$year)
	speciesdata<-droplevels(speciesdata)

	#Simulated effect sizes
	#latitudinal slope difference that is sufficient to generate 2 days difference per degree latitude

	#spatial slope - MA slope of 0.55 compared against slope of 1
	#spatial slope - correlation of 0.9

	#temporal slope - MA slope of 0.55 compared against slope of 1
	#temporal slope - correlation of 0.9

	store<-matrix(nrow=iteration,ncol=14,NA)

	maintitle<-c("latitudinal slope con","latitudinal slope res","slope diff","latslopediff P","space cor","space cor P","time cor","time cor P","MA space","Ma space P != 1","Ma space P != 0","Ma time", "Ma time p != 1","Ma time p != 0")

	for(iter in 1:iteration){

		spacemeans<-rmvnorm(n=length(levels(speciesdata$grid50)),sigma=matrix(nrow=2,ncol=2,c(con_spacevar, conres_spacecov,conres_spacecov,res_spacevar)),mean=c(0,0))
		rownames(spacemeans)<-levels(speciesdata$grid50)

		grid5meandevs<-rmvnorm(n=length(levels(speciesdata$grid5)),sigma=matrix(nrow=2,ncol=2,c(con_5kmvar, conres_5kmcov ,conres_5kmcov,res_5kmvar)),mean=c(0,0))
		rownames(grid5meandevs)<-levels(speciesdata$grid5)

		yearmeandevs<-rmvnorm(n=length(levels(speciesdata$year)),sigma=matrix(nrow=2,ncol=2,c(con_tempvar, conres_tempcov ,conres_tempcov,res_tempvar)),mean=c(0,0))
		rownames(yearmeandevs)<-levels(speciesdata$year)

		residdevs<-c(rnorm(length(consumer),mean=0,sqrt(con_residvar)),rnorm(length(resource),mean=0,sqrt(res_residvar)))

		consumerphenol<-con_intercept+ con_latslope*speciesdata$latcentre[speciesdata$spec=="con"]+spacemeans[speciesdata$grid50[speciesdata$spec=="con"],1]+ grid5meandevs[speciesdata$grid5[speciesdata$spec=="con"],1]+yearmeandevs[speciesdata$year[speciesdata$spec=="con"],1]

		resourcephenol<-res_intercept+ res_latslope*speciesdata$latcentre[speciesdata$spec=="res"]+spacemeans[speciesdata$grid50[speciesdata$spec=="res"],2]+ grid5meandevs[speciesdata$grid5[speciesdata$spec=="res"],2]+yearmeandevs[speciesdata$year[speciesdata$spec=="res"],2]

		speciesdata$phenol<-c(consumerphenol, resourcephenol)+ residdevs


#run the MCMCglmm model. 

		model<-MCMCglmm(phenol~spec-1+spec:latcentre, 
random=~us(spec):grid50+
us(spec):grid5+
us(spec):year,
rcov=~idh(spec):units,family=c("gaussian"),nitt=100000,burnin=20000,data=speciesdata,prior=pa_prior,pr=TRUE)


		#lat slopes and lat slope difference and the P value of the slope difference.
		latslopediff<-model$Sol[,"speccon:latcentre"]-model$Sol[,"specres:latcentre"]

		store[iter,1:4]<-c(mean(model$Sol[,"speccon:latcentre"]),mean(model$Sol[,"specres:latcentre"]),
		mean(latslopediff),2*min(c(length(which(latslopediff>0)),length(which(latslopediff <0))))/length(latslopediff))

		#correlation space (and signif)
		spacecorr<-model$VCV[,"specres:speccon.grid50"]/(sqrt(model$VCV[,"speccon:speccon.grid50"])*sqrt(model$VCV[,"specres:specres.grid50"]))
		spacecorrp<-2*min(c(length(which(spacecorr>0)),length(which(spacecorr<0))))/length(spacecorr)
		store[iter,5:6]<-c(mean(spacecorr),spacecorrp)

		#correlation time (and signif)
		timecorr<-model$VCV[,"specres:speccon.year"]/(sqrt(model$VCV[,"speccon:speccon.year"])*sqrt(model$VCV[,"specres:specres.year"]))
		timecorrp<-2*min(c(length(which(timecorr>0)),length(which(timecorr<0))))/length(timecorr)
		store[iter,7:8]<-c(mean(timecorr),timecorrp)


#Used to calculate the major axis slope
Macalc<-function(xvarcol,covcol,yvarcol){
	#based on code in the smatr package - https://cran.r-project.org/web/packages/smatr/smatr.pdf
	
	varx <-model$VCV[,xvarcol]
	covxy<-model$VCV[,covcol]
	vary <-model$VCV[,yvarcol]

		fac <- varx - vary
        b <- (fac + sqrt(fac^2 + 4 * covxy^2))/2/covxy
	return(b)}

###############



		#MA slope space (and signif)
		maspace<-Macalc(1,3,4)
		spacemap<-2*min(c(length(which(maspace>1)),length(which(maspace <1))))/length(maspace)
		spacemap0<-2*min(c(length(which(maspace>0)),length(which(maspace <0))))/length(maspace)
		store[iter,9:11]<-c(mean(maspace),spacemap,spacemap0)

		#MA slope time (and signif)
		matime<-Macalc(9,10,12)
		timemap<-2*min(c(length(which(matime>1)),length(which(matime <1))))/length(matime)
		timemap0<-2*min(c(length(which(matime>0)),length(which(matime <0))))/length(matime)
		store[iter,12:14]<-c(mean(matime),timemap,timemap0)

write.table(store,paste("yourfolder/",con_spec,".txt",sep=""),sep="\t",col.names= maintitle,row.names=F)

}
}


####################
#Used to test correlation and MA slope prior to run and identify parameters that will generate the desired effect sizes
#e.g., 
	varx <-41
	covxy<-63
	vary <-120

	fac <- varx - vary
        b <- (fac + sqrt(fac^2 + 4 * covxy^2))/2/covxy
corr<-covxy/(sqrt(varx)*sqrt(vary))
corr
b

#Here the paramaters are as generated by models with the exception of the spatial and temporal covariances and the consumer latitudinal slope, which we altered to generate effect sizes as described above.
#for the lat slope we simulate a difference of 2.

power.analysis(100,"piedf","frass",135.04 ,152.59,1.73,1.73+2,41/4,63/4,120/4,3.13,-0.37,5.75,41,63,120,44.68,33.33)

power.analysis(100,"frass","oak",144.26 ,116.65 ,3.01-2,3.01,41/4,63/4,120/4,7.23,1.48,24.82,41,63,120,33.35,54.8)

power.analysis(100,"greti","frass",118.96 ,148.73 ,1.93,1.93+2,41/4,63/4,120/4,17.4,-2.64,6.56,41,63,120,61.28,33.36)

power.analysis(100,"bluti","frass",118.3 ,147.85 ,1.67,1.67+2,41/4,63/4,120/4,10.18,-0.65,5.96,41,63,120,44.25,33.4)

##############
bluti<-read.csv("yourfolder/bluti.txt",sep="\t")

greti<-read.csv("yourfolder/greti.txt",sep="\t")

piedf<-read.csv("yourfolder/piedf.txt",sep="\t")

frass<-read.csv("yourfolder/frass.txt",sep="\t")


par(mfrow=c(3,2),mar=c(4,4,2,2),cex.lab=1.2)

barplot(c(length(intersect(which(bluti$latslopediff.P<=0.05),which(bluti$slope.diff<0)))/length(na.omit(bluti[,1])),
length(intersect(which(greti$latslopediff.P<=0.05),which(greti$slope.diff<0)))/length(na.omit(greti[,1])),
length(intersect(which(piedf$latslopediff.P<=0.05),which(piedf$slope.diff<0)))/length(na.omit(piedf[,1])),
length(intersect(which(frass$latslopediff.P<=0.05),which(frass$slope.diff<0)))/length(na.omit(frass[,1])))
,names.arg=c("bt","gt","pf","cat"),ylim=c(0,1),col=c("blue","yellow","black","dark green"),las=1,ylab="Power")

mtext(side=3,adj=0,"a) latitudinal slope difference ",cex=0.8,line=1)

plot(NULL,NULL)

barplot(c(length(intersect(which(bluti$space.cor.P <=0.05),which(bluti$space.cor >0)))/length(na.omit(bluti[,1])),
length(intersect(which(greti$space.cor.P <=0.05),which(greti$space.cor >0)))/length(na.omit(greti[,1])),
length(intersect(which(piedf$space.cor.P <=0.05),which(piedf$space.cor >0)))/length(na.omit(piedf[,1])),
length(intersect(which(frass$space.cor.P <=0.05),which(frass$space.cor >0)))/length(na.omit(frass[,1])))
,names.arg=c("bt","gt","pf","cat"),ylim=c(0,1),col=c("blue","yellow","black","dark green"),las=1,ylab="Power")

mtext(side=3,adj=0,"b) spatial correlation",cex=0.8,line=1)

barplot(c(length(intersect(which(bluti$Ma.space.P....1 <=0.05),which(bluti$MA.space <1)))/length(na.omit(bluti[,1])),
length(intersect(which(greti$Ma.space.P....1 <=0.05),which(greti$MA.space <1)))/length(na.omit(greti[,1])),
length(intersect(which(piedf$Ma.space.P....1 <=0.05),which(piedf$MA.space <1)))/length(na.omit(piedf[,1])),
length(intersect(which(frass$Ma.space.P....1 <=0.05),which(frass$MA.space < 1)))/length(na.omit(frass[,1])))
,names.arg=c("bt","gt","pf","cat"),ylim=c(0,1),col=c("blue","yellow","black","dark green"),las=1,ylab="Power")

mtext(side=3,adj=0,"c) spatial major axis",cex=0.8,line=1)


barplot(c(length(intersect(which(bluti$time.cor.P <=0.05),which(bluti$time.cor >0)))/length(na.omit(bluti[,1])),
length(intersect(which(greti$time.cor.P <=0.05),which(greti$time.cor >0)))/length(na.omit(greti[,1])),
length(intersect(which(piedf$time.cor.P <=0.05),which(piedf$time.cor >0)))/length(na.omit(piedf[,1])),
length(intersect(which(frass$time.cor.P <=0.05),which(frass$time.cor > 0)))/length(na.omit(frass[,1])))
,names.arg=c("bt","gt","pf","cat"),ylim=c(0,1),col=c("blue","yellow","black","dark green"),las=1,ylab="Power")

mtext(side=3,adj=0,"d) temporal correlation",cex=0.8,line=1)

barplot(c(length(intersect(which(bluti$Ma.time.p....1 <=0.05),which(bluti$Ma.time <1)))/length(na.omit(bluti[,1])),
length(intersect(which(greti$Ma.time.p....1 <=0.05),which(greti$Ma.time <1)))/length(na.omit(greti[,1])),
length(intersect(which(piedf$Ma.time.p....1 <=0.05),which(piedf$Ma.time <1)))/length(na.omit(piedf[,1])),
length(intersect(which(frass$Ma.time.p....1 <=0.05),which(frass$Ma.time < 1)))/length(na.omit(frass[,1])))
,names.arg=c("bt","gt","pf","cat"),ylim=c(0,1),col=c("blue","yellow","black","dark green"),las=1,ylab="Power")


mtext(side=3,adj=0,"e) temporal major axis",cex=0.8,line=1)


