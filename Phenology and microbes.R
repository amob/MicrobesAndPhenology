library(UpSetR)
library(MCMCglmm)

# analyze WOS records for AjB lit review on influence of microbes plant phenology
setwd("~/MicrobesAndPhenology/")

# wos_recs <- read.csv("~/saved records WOS friday sep 11 2020 - WOS records.tsv",
# 					sep="\t",stringsAsFactors=F,na.strings = c("NA",""))

cwos_recs <- read.csv("saved records WOS friday sep 11 2020 - Copy of WOS records.tsv",
					sep="\t",stringsAsFactors=F,na.strings = c("NA",""))


relevant <- cwos_recs$"relevantpt2"=="y"
PiM <- as.numeric(as.factor(cwos_recs$"Phenology....microbes...if.only.this..stop.."))==3 & !is.na(cwos_recs$"Phenology....microbes...if.only.this..stop..") # plant phenology influences microbes
PoM <- as.numeric(as.factor(cwos_recs$"phenology.of.micriobes...if.only.this..stop.."))==3 & !is.na(cwos_recs$"phenology.of.micriobes...if.only.this..stop..") #phenology of microbes
PiMPoM <- PiM | PoM
MiP <- !is.na(cwos_recs$What.phenological.trait.) # microbes influence plant phenology
MaltSonP <- cwos_recs$"Did.they.measure.selection...on.the.plant.only.."=="y" #phenology of microbes


minyear <- min(cwos_recs$Publication.Year,na.rm=T)
wosbyyear <- sapply(minyear:2020, function(yr) sum(cwos_recs$Publication.Year == yr,na.rm=T) )
wosbyyearR <- sapply(minyear:2020, function(yr) sum(cwos_recs[relevant,]$Publication.Year == yr,na.rm=T) )
wosbyyearPiM <- sapply(minyear:2020, function(yr) sum(cwos_recs[PiM,]$Publication.Year == yr,na.rm=T) )
wosbyyearPoM <- sapply(minyear:2020, function(yr) sum(cwos_recs[PoM,]$Publication.Year == yr,na.rm=T) )
wosbyyearPiMPoM <- sapply(minyear:2020, function(yr) sum(cwos_recs[PiMPoM,]$Publication.Year == yr,na.rm=T) )
wosbyyearMiP <- sapply(minyear:2020, function(yr) sum(cwos_recs[MiP,]$Publication.Year == yr,na.rm=T) )
wosbyyearSEL <- sapply(minyear:2020, function(yr) sum(cwos_recs[MaltSonP,]$Publication.Year == yr,na.rm=T) )

years <- c(minyear:2020)
pdf("WOS_AjB_micr_phen.pdf",height=3.5,width=3.5)
par(mar=c(3,3,1,1))
par(oma=c(0,0,0,0))
plot(wosbyyear~years,type="l",lty=2,xlab="",ylab="")
	lines(wosbyyearR~years,col=rgb(0,0,0))
# 	lines(wosbyyearPiM~years,col=rgb(0.75,0,0))
# 	lines(wosbyyearPoM~years,col=rgb(0,0,0.75))
	lines(wosbyyearPiMPoM~years,col=rgb(0,0,0.75))
	lines(wosbyyearMiP~years,col=rgb(0.75,0,0))
	legend(1990,78,c("all records","relevant","phen -> micr, or phen of micr","micr -> phen"),
			 lty=c(2,1,1,1,1),col=rgb( c(0,0,0.75,0),0,c(0,0,0,0.75) ) ,bty="n"   ) 
	mtext("WOS records",side=2,line=2)
	mtext("Year",side=1,line=2)
dev.off()


 pdf("WOS_AjB_micr_phen_r_sim.pdf",height=3.5,width=5.5)
par(mar=c(3,3,1,1))
par(oma=c(0,0,0,0))
plot(wosbyyearR~years,type="l",lty=2,xlab="",ylab="",ylim=c(0,25))
	lines(wosbyyearPiMPoM~years)
	lines(wosbyyearMiP~years, col=rgb(0.75,0,0))
	lines(wosbyyearSEL~years,col=rgb(0,0,0.75))
	legend(1987,26,c("plant microbes & phenology","microbe phenology (+ plant phenology -> microbes)","microbes affect phenology?","microbes alter selection on phenology?"),
			 lty=c(2,1,1 ),col=c(rgb(0,0,0),rgb(0,0,0),rgb(0.75,0,0),rgb(0,0,0.75)) ,bty="n"   ) 
	mtext("WOS records",side=2,line=2)
	mtext("Year",side=1,line=2)
dev.off()

table(relevant)
table(PiMPoM)
table(MiP)
table(cwos_recs$any.effect.in.study)

wos_rel <- cwos_recs[which(relevant),]

sim2and <- gsub("; ", wos_rel$microbetax.sim2,replacement="&")
simand <- gsub("; ", wos_rel$microbetax.sim,replacement="&")
loc2and <- gsub("; ", wos_rel$Microbelocation.2,replacement="&")
phen.and <- gsub("; ", wos_rel$phen.trait,replacement="&")

pdf("upsetrtax.pdf",width=4,height=5)
upset(fromExpression(table(sim2and)),nsets=15, order.by = "freq")
dev.off()
pdf("upsetrtaxcomp.pdf",width=4,height=5)
upset(fromExpression(table(simand)),nsets=15, order.by = "freq")
dev.off()
pdf("upsetrloc.pdf",width=4,height=5)
upset(fromExpression(table(loc2and)),nsets=10, order.by = "freq")
dev.off()
pdf("upsetrtrt.pdf",width=4,height=5)
upset(fromExpression(table(phen.and)), nsets=15, order.by = "freq")
dev.off()

# upset(fromExpression(table(phen.and)),set.metadata = cbind(phen.and), nsets=15, order.by = "freq")

effobs <- ifelse(wos_rel$any.effect.in.study=="y",1,0)

unqtype <- unique(unlist(sapply(sim2and, function(z) strsplit(z,"&"))))
unqloc <- unique(unlist(sapply(loc2and, function(z) strsplit(z,"&"))))
unqtrt <- unique(unlist(sapply(phen.and, function(z) strsplit(z,"&"))))

micrtypetab <- as.data.frame(sapply(unqtype[!is.na(unqtype)], function(mtype) as.numeric(grepl(as.character(mtype),sim2and))))
micrloctab <- as.data.frame(sapply(unqloc[!is.na(unqloc)], function(ltype) as.numeric(grepl(as.character(ltype),loc2and))))
traittab <- as.data.frame(sapply(unqtrt[!is.na(unqtrt)], function(ttype) as.numeric(grepl(as.character(ttype),phen.and))))

tapply(effobs,phen.and,mean,na.rm=T)
tapply(effobs,sim2and,mean,na.rm=T)
tapply(effobs,sim2and,sum,na.rm=T)

tapply(effobs,micrtypetab$bacteria,sum,na.rm=T)
tapply(effobs,micrtypetab$MF,mean,na.rm=T)


tapply(effobs,loc2and,mean,na.rm=T)
# trtexp <- as.factor(as.character(unlist(sapply(wos_rel$phen.trait, function(x) unlist(strsplit(x,split="; "))))))
fails <- 1-effobs

effcatdat <- data.frame(cbind(effobs, fails, micrtypetab,micrloctab,traittab))


modmtypeBact <- MCMCglmm(cbind(effobs,fails)~bacteria,family ="multinomial2" ,data=effcatdat,pr=T,verbose=F,nitt=100000,burnin=10000,thin=100)#n.s. more sig effs than avg rate without
modmtypeFungi <- MCMCglmm(cbind(effobs,fails)~fungi,family ="multinomial2" ,data=effcatdat,pr=T,verbose=F,nitt=100000,burnin=10000,thin=100)#n.s. more sig effects than avg rate without
modmtypeMF <- MCMCglmm(cbind(effobs,fails)~MF,family ="multinomial2" ,data=effcatdat,pr=T,verbose=F,nitt=100000,burnin=10000,thin=100) #n.s. fewer significant effects than average rate without
modmtypeMixed <- MCMCglmm(cbind(effobs,fails)~mixed,family ="multinomial2" ,data=effcatdat,pr=T,verbose=F,nitt=100000,burnin=10000,thin=100)#n.s. more sig effs than avg rate without
modmtypeVirus <- MCMCglmm(cbind(effobs,fails)~virus,family ="multinomial2" ,data=effcatdat,pr=T,verbose=F,nitt=100000,burnin=10000,thin=100) #too few obs. not to be trusted

modlocroot <- MCMCglmm(cbind(effobs,fails)~root,family ="multinomial2" ,data=effcatdat,pr=T,verbose=F,nitt=100000,burnin=10000,thin=100)#n.s. more sig effs than avg rate without
modlocshoot <- MCMCglmm(cbind(effobs,fails)~shoot,family ="multinomial2" ,data=effcatdat,pr=T,verbose=F,nitt=100000,burnin=10000,thin=100)#n.s. fewer sig effs than avg rate without
modlocreproductive <- MCMCglmm(cbind(effobs,fails)~reproductive,family ="multinomial2" ,data=effcatdat,pr=T,verbose=F,nitt=100000,burnin=10000,thin=100)#n.s. fewer sig effs than avg rate without
modlocseed <- MCMCglmm(cbind(effobs,fails)~seed,family ="multinomial2" ,data=effcatdat,pr=T,verbose=F,nitt=100000,burnin=10000,thin=100)#n.s. more sig effs than avg rate without

##model fit issues....effective samples is low in most
modtraitflwrT <- MCMCglmm(cbind(effobs,fails)~flowering.time,family ="multinomial2" ,data=effcatdat,pr=T,verbose=F,nitt=500000,burnin=10000,thin=100)#n.s. fewer sig effs than avg rate without
modtraitfbT <- MCMCglmm(cbind(effobs,fails)~floral.budset.time,family ="multinomial2" ,data=effcatdat,pr=T,verbose=F,nitt=500000,burnin=10000,thin=100)#n.s. more sig effs than avg rate without
modtraitgT <- MCMCglmm(cbind(effobs,fails)~germination.time,family ="multinomial2" ,data=effcatdat,pr=T,verbose=F,nitt=500000,burnin=10000,thin=100)#n.s. more sig effs than avg rate without
modtraitgP <- MCMCglmm(cbind(effobs,fails)~germination.prb,family ="multinomial2" ,data=effcatdat,pr=T,verbose=F,nitt=500000,burnin=10000,thin=100)#n.s. more sig effs than avg rate without
modtraitphyl <- MCMCglmm(cbind(effobs,fails)~phyllochron,family ="multinomial2" ,data=effcatdat,pr=T,verbose=F,nitt=500000,burnin=10000,thin=100)#"sig" more, but there are only 2 studies
modtraitsncT <- MCMCglmm(cbind(effobs,fails)~senescence.time,family ="multinomial2" ,data=effcatdat,pr=T,verbose=F,nitt=500000,burnin=10000,thin=100)#n.s. fewer sig effs than avg rate without
modtraitmatT <- MCMCglmm(cbind(effobs,fails)~maturation.time,family ="multinomial2" ,data=effcatdat,pr=T,verbose=F,nitt=500000,burnin=10000,thin=100)#n.s. more sig effs than avg rate without
modtraitfrtT <- MCMCglmm(cbind(effobs,fails)~fruiting.time,family ="multinomial2" ,data=effcatdat,pr=T,verbose=F,nitt=500000,burnin=10000,thin=100)#n.s. slightly fewer sig effs than avg rate without
modtraitflwrD <- MCMCglmm(cbind(effobs,fails)~flowering.duration,family ="multinomial2" ,data=effcatdat,pr=T,verbose=F,nitt=500000,burnin=10000,thin=100)#n.s. more sig effs than avg rate without
modtraitrdT <- MCMCglmm(cbind(effobs,fails)~root.development.timing,family ="multinomial2" ,data=effcatdat,pr=T,verbose=F,nitt=500000,burnin=10000,thin=100)#n.s. more sig effs than avg rate without
modtraitbudB <- MCMCglmm(cbind(effobs,fails)~budburst,family ="multinomial2" ,data=effcatdat,pr=T,verbose=F,nitt=500000,burnin=10000,thin=100)#n.s. more sig effs than avg rate without; close to sig, but there's only three studies here
