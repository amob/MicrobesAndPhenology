library(UpSetR) #for set plots
library(MCMCglmm) #for models
library(grDevices) #for color

HPDi <- function(vect,prob){
	int <- HPDinterval(as.mcmc(vect),prob=prob)
	return(int)
}

invlogistic <- function(x){exp(x)/(1+ exp(x))}
std.error <- function(dat, na.rm=TRUE) {sd(dat,na.rm=na.rm)/sqrt(length(dat))}#defaults to na.rm=T


################################################################################
####read in data
################################################################################


# analyze WOS records for AjB lit review on influence of microbes plant phenology
setwd("~/MicrobesAndPhenology/")

# wos_recs <- read.csv("~/saved records WOS friday sep 11 2020 - WOS records.tsv",
# 					sep="\t",stringsAsFactors=F,na.strings = c("NA",""))

cwos_recs <- read.csv("saved records WOS friday sep 11 2020 - Copy of WOS records.tsv",
					sep="\t",stringsAsFactors=F,na.strings = c("NA","")) #might consider a simplified version of this, so that column names are not as crazy
split_recs <- read.csv("trimmed RECORDS FOR SPLITTING.csv",
					sep="\t",stringsAsFactors=F,na.strings = c("NA","")) 
			#hand removed some columns, re-enforced the 5 options for direction of effect
			#also re-enforced options for mating, microbe nutrient-pathogen-other; Micrloc.2; state.relative.too, incl 1 "sterilized soils" to "without focal microbe"
			#also, there seem to be some extra rows still, perhaps those should be reomoved.
			#there are a few notes about changes in how records are scored, might want to make sure those are reflected in the previous tab as well?
			#Deleted a few effect directions to be blank because it seemed like they were a pooled version of split records? (e.g. still had two traits mixed --"flowering time; fruiting time" -- this might not be right)
split_recs_scr <- split_recs[!is.na(split_recs$direction.effect.split),]

################################################################################
####Analysis of all broadly relevant records			
################################################################################


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
wos_MiP <- cwos_recs[which(MiP),]

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

#dataframe for models
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




################################################################################
####Analysis of phenological timing records
################################################################################

#quick survey
table(paste(split_recs_scr$micrtax.sim2,split_recs_scr$direction.effect.split))
# 
#   bacteria delayed   bacteria earlier  bacteria expanded  bacteria narrowed      bacteria none 
#                  9                 68                 47                  5                127 
#      fungi delayed      fungi earlier     fungi expanded     fungi narrowed         fungi none 
#                 62                 20                  2                 17                 58 
#         MF delayed         MF earlier        MF expanded        MF narrowed            MF none 
#                 11                 72                  8                  2                 46 
#  MF; fungi earlier MF; fungi expanded     MF; fungi none      mixed delayed      mixed earlier 
#                  2                  3                  7                 38                 29 
#     mixed expanded     mixed narrowed         mixed none 
#                  8                 10                228 
table(paste(split_recs_scr$culture.sim,split_recs_scr$direction.effect.split))
#         culture delayed         culture earlier        culture expanded        culture narrowed 
#                      79                     147                      56                      12 
#            culture none culture; direct delayed            cutlure none          direct delayed 
#                     179                       4                       2                      27 
#          direct earlier         direct expanded         direct narrowed             direct none 
#                      33                       7                       9                     215 
#        observed earlier       observed expanded         removed delayed         removed earlier 
#                       1                       1                      10                      10 
#        removed expanded        removed narrowed            removed none 
#                       4                      13                      70 
table(paste(split_recs_scr$Strainvcomm.1,split_recs_scr$direction.effect.split))
#   community delayed   community earlier  community expanded  community narrowed 
#                  38                  35                  11                  12 
#      community none      strain delayed      strain earlier     strain expanded 
#                 240                  58                 134                  49 
#  strain mix delayed  strain mix earlier strain mix expanded     strain mix none 
#                  24                  22                   8                  21 
#     strain narrowed         strain none 
#                  22                 205 
table(paste(split_recs_scr$phen.trait.2,split_recs_scr$direction.effect.split))

split_recs_scr$culture.sim[split_recs_scr$culture.sim== "culture; direct"] <- "culture" # making a decision on Lu2018 
split_recs_scr <- split_recs_scr[split_recs_scr$micrtax.sim2!= "MF; fungi",] #subtract the few records in this funny category
#might consider keeping them for other analyses?
split_recs_scr$micrtax.sim2[split_recs_scr$micrtax.sim2=="fungi"] <- "otherfungi" #changing alphabetical position for convenience


numeric.phen.eff <- rep(NA, length=nrow(split_recs_scr))
numeric.phen.eff[split_recs_scr$direction.effect.split=="earlier"] <- -1
numeric.phen.eff[split_recs_scr$direction.effect.split=="delayed"] <- 1
numeric.phen.eff[split_recs_scr$direction.effect.split=="none"] <- 0
numeric.phen.eff[split_recs_scr$direction.effect.split=="expanded"] <- 1
numeric.phen.eff[split_recs_scr$direction.effect.split=="narrowed"] <- -1
split_recs_scr$numeric.phen.eff <- numeric.phen.eff
split_recs_scr$lifeform2 <- split_recs_scr$lifeform
split_recs_scr$lifeform2[  split_recs_scr$lifeform %in% c("annual in cultivation (perennnial)","perennial; annual (both)","both") ] <-  "both"
split_recs_scr$lifeform2[  split_recs_scr$lifeform =="perennial (biennial?)" ] <-  "perennial"
split_recs_scr$mating2 <- split_recs_scr$mating
split_recs_scr$mating2 [is.na(split_recs_scr$mating)] <- "unk"
split_recs_scr$mating2 [split_recs_scr$mating=="both"] <- "s+o"

split_recs_scr$micrloc3 <- split_recs_scr$Micrloc.2
split_recs_scr$micrloc3 [split_recs_scr$Micrloc.2%in% c("reproductive; shoot","root; shoot; reproductive","root; shoot; reproductive; seed")] <- "various incl reproductive"
split_recs_scr$micrloc3 [split_recs_scr$Micrloc.2%in% c("seed; shoot","root; shoot","root; seed")] <- "various not incl reproductive"
split_recs_scr$micrloc3 [split_recs_scr$Micrloc.2=="seed "] <- "seed"

##variables for binomial versions
split_recs_scr$isEarlyOrNrw <- ifelse(numeric.phen.eff=="-1",1,0)
split_recs_scr$notEarlyOrNrw <- 1-split_recs_scr$isEarlyOrNrw
split_recs_scr$isDelayOrExp <- ifelse(numeric.phen.eff=="1",1,0)
split_recs_scr$notDelayOrExp <- 1-split_recs_scr$isDelayOrExp
split_recs_scr$isSig <- ifelse(numeric.phen.eff=="1" | numeric.phen.eff=="-1",1,0)
split_recs_scr$notSig <- 1- split_recs_scr$isSig 

# numeric.phen.effpos <- split_recs_scr$numeric.phen.eff +2
#  split_recs_scr$numeric.phen.effpos <- numeric.phen.effpos

#inermediate microb taxa variable
split_recs_scr$micrtax.sim1.5 <- split_recs_scr$micrtax.sim1
split_recs_scr$micrtax.sim1.5[split_recs_scr$micrtax.sim1.5%in% c("AMF","EMF","EricoidMF","OrchidMF")] <- "MF"
split_recs_scr$micrtax.sim1.5[split_recs_scr$micrtax.sim1.5=="fungi"] <- "otherfungi"
split_recs_scr$micrtax.sim1.5[split_recs_scr$micrtax.sim1.5=="bacteria"] <- "bacteria,other"
split_recs_scr$micrtax.sim1.5[split_recs_scr$micrtax.sim1.5 %in% c("bacteria; rhizobia","rhizobia")] <- "rhizobia+" #to avoid additional category...

#simpler phenotypes variable
split_recs_scr$phen.trait.3 <- split_recs_scr$phen.trait.2
split_recs_scr$phen.trait.3[split_recs_scr$phen.trait.2 %in% c("maturation time","floral budset time","flowering time","peak flowering","flower senescence time","fruiting time")] <- "reproductive time"
split_recs_scr$phen.trait.3[split_recs_scr$phen.trait.2 %in% c("budburst","phyllochron")] <- "vegetative time"


flowering <- split_recs_scr[split_recs_scr$phen.trait.2=="flowering time",]
	flowering$dummyWOSrelrank <- flowering$WOS.relevance.rank
	flowering$dummyWOSrelrank[flowering$dummyWOSrelran%in%names(table(flowering$WOS.relevance.rank)[table(flowering$WOS.relevance.rank)<13])] <- "binned"
floweringMAT <- split_recs_scr[split_recs_scr$phen.trait.2=="flowering time" & split_recs_scr$mating2 != "unk" & split_recs_scr$mating2 != "s+o",]
	floweringMAT$dummyWOSrelrank <- floweringMAT$WOS.relevance.rank
	floweringMAT$dummyWOSrelrank[floweringMAT$dummyWOSrelran%in%names(table(floweringMAT$WOS.relevance.rank)[table(floweringMAT$WOS.relevance.rank)<13])] <- "binned"
germtime <- split_recs_scr[split_recs_scr$phen.trait.2=="germination time",]
	germtime$dummyWOSrelrank <- germtime$WOS.relevance.rank
	germtime$dummyWOSrelrank[germtime$dummyWOSrelran%in%names(table(germtime$WOS.relevance.rank)[table(germtime$WOS.relevance.rank)<13])] <- "binned"
germtimeMAT <- split_recs_scr[split_recs_scr$phen.trait.2=="germination time" & split_recs_scr$mating2 != "unk" & split_recs_scr$mating2 != "s+o",]
	germtimeMAT$dummyWOSrelrank <- germtimeMAT$WOS.relevance.rank
	germtimeMAT$dummyWOSrelrank[germtimeMAT$dummyWOSrelran%in%names(table(germtimeMAT$WOS.relevance.rank)[table(germtimeMAT$WOS.relevance.rank)<13])] <- "binned"
germprb <- split_recs_scr[split_recs_scr$phen.trait.2=="germination prb"& split_recs_scr$mating2 != "unk" & split_recs_scr$mating2 != "s+o",]
	germprb$dummyWOSrelrank <- germprb$WOS.relevance.rank
	germprb$dummyWOSrelrank[germprb$dummyWOSrelran%in%names(table(germprb$WOS.relevance.rank)[table(germprb$WOS.relevance.rank)<13])] <- "binned"
# 	germprb<- germprb[-which(germprb$culture.sim=="observed"),] #removing the sole record that is not manipulated
#no longer need to remove the observed record; as not checking how this relates to results 
germprbMAT <- split_recs_scr[split_recs_scr$phen.trait.2=="germination prb" & split_recs_scr$mating2 != "unk" & split_recs_scr$mating2 != "s+o",]
	germprbMAT$dummyWOSrelrank <- germprbMAT$WOS.relevance.rank
	germprbMAT$dummyWOSrelrank[germprbMAT$dummyWOSrelran%in%names(table(germprbMAT$WOS.relevance.rank)[table(germprbMAT$WOS.relevance.rank)<13])] <- "binned"
# 	germprbMAT<- germprbMAT[-which(germprbMAT$culture.sim=="observed"),] #removing the sole record that is not manipulated


priornr=list(R=list(V= 1, fix=1))
#I'm still not sure why/if R should be fixed, but JH has it fixed in ALL example priors, so I think so?
prior=list(R=list(V= 1, fix=1), G=list(G1=list(V=1, nu=0)))
priornuup=list(R=list(V= 1, fix=1), G=list(G1=list(V=1, nu=8)))#turning this degree of belief way up solved sampling issues for flowering time
	#estimates for intercept and random effects became more independent, because this shunts variation to be explained by fixed effect when there isn't information in the data for the random effect
prior2=list(R=list(V= 1, fix=1), G=list(G1=list(V=1, nu=0),G2=list(V=1, nu=0)))
priornuup2=list(R=list(V= 1, fix=1), G=list(G1=list(V=1, nu=8),G2=list(V=1, nu=8)))
#after fitting the flowering model with various priors for the case with random effects(prior 2 onwards, the latter two of which are meant for no-intercept models); and various levels of binning,
###I conclude that the prior doesn't seem to help with the correlation issues between the posterior fixed and random effects. These issues are exacerbated by two little or too much binning
####effective sampling for the random effect was higher with the most complex prior, but it reduced DIC and also effective sampling (though still acceptable) for cutpoint
####no-intercept models seem to make all problems worse.

######BINOMIALS
# flowerbiINT <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~1,random = ~dummyWOSrelrank,data=flowering,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)

flowerbitax <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~micrtax.sim2,random = ~dummyWOSrelrank,data=flowering,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
#a a b b 
flowerbiscom <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~Strainvcomm.1,random = ~dummyWOSrelrank,data=flowering,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
flowerbimat <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~mating2,random = ~dummyWOSrelrank,data=floweringMAT,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
flowerbilh <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~lifeform2,random = ~dummyWOSrelrank,data=flowering,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
flowerbiloc <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~micrloc3,random = ~dummyWOSrelrank,data=flowering,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
flowerDbitax <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~micrtax.sim2,random = ~dummyWOSrelrank,data=flowering,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
flowerDbiscom <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~Strainvcomm.1,random = ~dummyWOSrelrank,data=flowering,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
flowerDbimat <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~mating2,random = ~dummyWOSrelrank,data=floweringMAT,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
flowerDbilh <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~lifeform2,random = ~dummyWOSrelrank,data=flowering,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
flowerDbiloc <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~micrloc3,random = ~dummyWOSrelrank,data=flowering,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
earlyFmns <- lapply(c("micrtax.sim2","Strainvcomm.1","mating2","lifeform2","micrloc3"),
		 function(z) tapply(flowering$isEarlyOrNrw,flowering[,colnames(flowering)==z],mean) )
	earlyFmns[[3]] <- tapply(floweringMAT$isEarlyOrNrw,floweringMAT[,colnames(floweringMAT)=="mating2"],mean)
	earlyFmns[[6]] <- mean(floweringMAT$isEarlyOrNrw)
earlyFses <- lapply(c("micrtax.sim2","Strainvcomm.1","mating2","lifeform2","micrloc3"),
		 function(z) tapply(flowering$isEarlyOrNrw,flowering[,colnames(flowering)==z],std.error) )
	earlyFses[[3]] <- tapply(floweringMAT$isEarlyOrNrw,floweringMAT[,colnames(floweringMAT)=="mating2"],std.error)
	earlyFses[[6]] <- std.error(floweringMAT$isEarlyOrNrw)
delayFmns <- lapply(c("micrtax.sim2","Strainvcomm.1","mating2","lifeform2","micrloc3"),
		 function(z) tapply(flowering$isDelayOrExp,flowering[,colnames(flowering)==z],mean) )
	delayFmns[[3]] <- tapply(floweringMAT$isDelayOrExp,floweringMAT[,colnames(floweringMAT)=="mating2"],mean)
	delayFmns[[6]] <- mean(floweringMAT$isDelayOrExp)
delayFses <- lapply(c("micrtax.sim2","Strainvcomm.1","mating2","lifeform2","micrloc3"),
		 function(z) tapply(flowering$isDelayOrExp,flowering[,colnames(flowering)==z],std.error) )
	delayFses[[3]] <- tapply(floweringMAT$isDelayOrExp,floweringMAT[,colnames(floweringMAT)=="mating2"],std.error)
	delayFses[[6]] <- std.error(floweringMAT$isDelayOrExp)

#more likely to be earlier:bacteria vs fungi (simplified micr only); culture inoc vs removed exps; strains vs mixed; ; annuals vs perennials
germtimebitax <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~micrtax.sim2,random = ~dummyWOSrelrank,data=germtime,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
#b a c b
germtimebiscom <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~Strainvcomm.1,random = ~dummyWOSrelrank,data=germtime,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
germtimebimat <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~mating2,random = ~dummyWOSrelrank,data=germtimeMAT,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
germtimebilh <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~lifeform2,random = ~dummyWOSrelrank,data=germtime,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
germtimebiloc <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~micrloc3,random = ~dummyWOSrelrank,data=germtime,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
germtimeDbitax <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~micrtax.sim2,random = ~dummyWOSrelrank,data=germtime,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
germtimeDbiscom <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~Strainvcomm.1,random = ~dummyWOSrelrank,data=germtime,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
germtimeDbimat <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~mating2,random = ~dummyWOSrelrank,data=germtimeMAT,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
germtimeDbilh <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~lifeform2,random = ~dummyWOSrelrank,data=germtime,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
germtimeDbiloc <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~micrloc3,random = ~dummyWOSrelrank,data=germtime,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
earlyGmns <- lapply(c("micrtax.sim2","Strainvcomm.1","mating2","lifeform2","micrloc3"),
		 function(z) tapply(germtime$isEarlyOrNrw,germtime[,colnames(germtime)==z],mean) )
	earlyGmns[[3]] <- tapply(germtimeMAT$isEarlyOrNrw,germtimeMAT[,colnames(germtimeMAT)=="mating2"],mean)
	earlyGmns[[6]] <- mean(germtime$isEarlyOrNrw)
earlyGses <- lapply(c("micrtax.sim2","Strainvcomm.1","mating2","lifeform2","micrloc3"),
		 function(z) tapply(germtime$isEarlyOrNrw,germtime[,colnames(germtime)==z],std.error) )
	earlyGses[[3]] <- tapply(germtimeMAT$isEarlyOrNrw,germtimeMAT[,colnames(germtimeMAT)=="mating2"],std.error)
	earlyGses[[6]] <- std.error(germtime$isEarlyOrNrw)
delayGmns <- lapply(c("micrtax.sim2","Strainvcomm.1","mating2","lifeform2","micrloc3"),
		 function(z) tapply(germtime$isDelayOrExp,germtime[,colnames(germtime)==z],mean) )
	delayGmns[[3]] <- tapply(germtimeMAT$isDelayOrExp,germtimeMAT[,colnames(germtimeMAT)=="mating2"],mean)
	delayGmns[[6]] <- mean(germtime$isDelayOrExp)
delayGses <- lapply(c("micrtax.sim2","Strainvcomm.1","mating2","lifeform2","micrloc3"),
		 function(z) tapply(germtime$isDelayOrExp,germtime[,colnames(germtime)==z],std.error) )
	delayGses[[3]] <- tapply(germtimeMAT$isDelayOrExp,germtimeMAT[,colnames(germtimeMAT)=="mating2"],std.error)
	delayGses[[6]] <- std.error(germtime$isDelayOrExp)


germprbbitax <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~micrtax.sim2,random = ~dummyWOSrelrank,data=germprb,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
germprbbiscom <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~Strainvcomm.1,random = ~dummyWOSrelrank,data=germprb,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
germprbbimat <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~mating2,random = ~dummyWOSrelrank,data=germprb,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
germprbbilh <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~lifeform2,random = ~dummyWOSrelrank,data=germprb,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
germprbbiloc <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~micrloc3,random = ~dummyWOSrelrank,data=germprb,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
germprbDbitax <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~micrtax.sim2,random = ~dummyWOSrelrank,data=germprb,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
germprbDbiscom <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~Strainvcomm.1,random = ~dummyWOSrelrank,data=germprb,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
germprbDbimat <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~mating2,random = ~dummyWOSrelrank,data=germprb,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
germprbDbilh <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~lifeform2,random = ~dummyWOSrelrank,data=germprb,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
germprbDbiloc <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~micrloc3,random = ~dummyWOSrelrank,data=germprb,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
earlyGPmns <- lapply(c("micrtax.sim2","Strainvcomm.1","mating2","lifeform2","micrloc3"),
		 function(z) tapply(germprb$isEarlyOrNrw,germprb[,colnames(germprb)==z],mean) )
	earlyGPmns[[3]] <- tapply(germprbMAT$isEarlyOrNrw,germprbMAT[,colnames(germprbMAT)=="mating2"],mean)
	earlyGPmns[[6]] <- mean(germprb$isEarlyOrNrw)
earlyGPses <- lapply(c("micrtax.sim2","Strainvcomm.1","mating2","lifeform2","micrloc3"),
		 function(z) tapply(germprb$isEarlyOrNrw,germprb[,colnames(germprb)==z],std.error) )
	earlyGPses[[3]] <- tapply(germprbMAT$isEarlyOrNrw,germprbMAT[,colnames(germprbMAT)=="mating2"],std.error)
	earlyGPses[[6]] <- std.error(germprb$isEarlyOrNrw)
delayGPmns <- lapply(c("micrtax.sim2","Strainvcomm.1","mating2","lifeform2","micrloc3"),
		 function(z) tapply(germprb$isDelayOrExp,germprb[,colnames(germprb)==z],mean) )
	delayGPmns[[3]] <- tapply(germprbMAT$isDelayOrExp,germprbMAT[,colnames(germprbMAT)=="mating2"],mean)
	delayGPmns[[6]] <- mean(germprb$isDelayOrExp)
delayGPses <- lapply(c("micrtax.sim2","Strainvcomm.1","mating2","lifeform2","micrloc3"),
		 function(z) tapply(germprb$isDelayOrExp,germprb[,colnames(germprb)==z],std.error) )
	delayGPses[[3]] <- tapply(germprbMAT$isDelayOrExp,germprbMAT[,colnames(germprbMAT)=="mating2"],std.error)
	delayGPses[[6]] <- std.error(germprb$isDelayOrExp)


pdf("means_ses_and prelim fitted diffs from binom.pdf",height=8,width=5)
xvals <- c(seq(from=0,to=1,length.out=c(4+3+2+3+5+1) ))
# xlab <- c("B","MF","MX","OF","B","MF","MX","OF","R","cltr","dir","rm","str","strs","com","out","os","self","unk","ann","b","per")
xlab <- c("bacteria","myc fungi","mixed","other fungi","community","strain","strain mix","outcrosser","selfer","annual","both","perennial","root","seed","shoot","v. w/ reprod","v. w/o reprod","All")
par(mfrow=c(3,1))
par(mar=c(1,0,1,0))
par(oma=c(5,4,1,1))
plot(c(0,1)~I(c(0,1)),pch=NA,xlab="",ylab="",bty="n",xaxt="n",ylim=c(-0.2,1.2),yaxt="n")
	points(earlyFmns[[1]]~xvals[1:4]  ,pch=16:20); 	points(delayFmns[[1]]~I(xvals[1:4] +0.01 ),pch=16:20,col="gray")
	arrows(x0=xvals[1:4] , y0=earlyFmns[[1]] -earlyFses[[1]], y1 = earlyFmns[[1]] +earlyFses[[1]],length=0     )
	arrows(x0=xvals[1:4]+0.01 , y0=delayFmns[[1]] -delayFses[[1]], y1 = delayFmns[[1]] +delayFses[[1]],length=0  ,col="gray"   )
	abline(v= sum(xvals[4:5])/2,lty=3 )
	text(x=xvals[1:4] ,y=-0.2,c("a","a","ab","b"))
	text(x=xvals[1:4] ,y=1.2,c("b","ab","a","ab"),col="gray")
	points(earlyFmns[[2]]~xvals[5:7]  ,pch=16:20); 		points(delayFmns[[2]]~I(xvals[5:7] +0.01)  ,pch=16:20,col="gray")
	arrows(x0=xvals[5:7] , y0=earlyFmns[[2]] -earlyFses[[2]], y1 = earlyFmns[[2]] +earlyFses[[2]],length=0     )
	arrows(x0=xvals[5:7] + 0.01, y0=delayFmns[[2]] -delayFses[[2]], y1 = delayFmns[[2]] +delayFses[[2]],length=0 ,col="gray"    )
	text(x=xvals[5:7] ,y=-0.2,c("b","ab","a")) #
	text(x=xvals[5:7] ,y=1.2,c("a","b","b"),col="gray")#NOTE THIS PREDICTION LOOKS REVERSED
	abline(v= sum(xvals[7:8])/2,lty=3 )
	points(earlyFmns[[3]]~xvals[8:9]  ,pch=16:20); 	points(delayFmns[[3]]~I(xvals[8:9] +0.01) ,pch=16:20, col="gray")
	arrows(x0=xvals[8:9] , y0=earlyFmns[[3]] -earlyFses[[3]], y1 = earlyFmns[[3]] +earlyFses[[3]],length=0     )
	arrows(x0=xvals[8:9]+0.01 , y0=delayFmns[[3]] -delayFses[[3]], y1 = delayFmns[[3]] +delayFses[[3]],length=0 ,col="gray"    )
	text(x=xvals[8:9] ,y=-0.2,c(" "," "))
	text(x=xvals[8:9] ,y=1.2,c(" "," "),col="gray")
	abline(v= sum(xvals[9:10])/2,lty=3 )
	points(earlyFmns[[4]]~xvals[10:12]  ,pch=16:20); 	points(delayFmns[[4]]~I(xvals[10:12]+0.01)  ,pch=16:21,col="gray")
	arrows(x0=xvals[10:12] , y0=earlyFmns[[4]] -earlyFses[[4]], y1 = earlyFmns[[4]] +earlyFses[[4]],length=0     )
	arrows(x0=xvals[10:12]+0.01 , y0=delayFmns[[4]] -delayFses[[4]], y1 = delayFmns[[4]] +delayFses[[4]],length=0 ,col="gray"    )
	text(x=xvals[10:12] ,y=-0.2,c("a","ab","b"))
	text(x=xvals[10:12] ,y=1.2,c(" "," "," "),col="gray")
	abline(v= sum(xvals[12:13])/2,lty=3 )
	points(earlyFmns[[5]]~xvals[13:17]  ,pch=16:20); 	points(delayFmns[[5]]~I(xvals[13:17]+0.01)  ,pch=16:20,col="gray")
	arrows(x0=xvals[13:17] , y0=earlyFmns[[5]] -earlyFses[[5]], y1 = earlyFmns[[5]] +earlyFses[[5]],length=0     )
	arrows(x0=xvals[13:17] +0.01, y0=delayFmns[[5]] -delayFses[[5]], y1 = delayFmns[[5]] +delayFses[[5]],length=0 ,col="gray"    )
	text(x=xvals[13:17] ,y=-0.2,c("a","b","b","b","ab"))
	text(x=xvals[13:17] ,y=1.2,c("a","b","ab","ab","ab"),col="gray")
	abline(v= sum(xvals[17:18])/2,lty=3 )
	points(earlyFmns[[6]]~xvals[18]  ,pch=16:20); 	points(delayFmns[[6]]~I(xvals[18]+0.01)  ,pch=16:20,col="gray")
	arrows(x0=xvals[18] , y0=earlyFmns[[6]] -earlyFses[[6]], y1 = earlyFmns[[6]] +earlyFses[[6]],length=0     )
	arrows(x0=xvals[18] +0.01, y0=delayFmns[[6]] -delayFses[[6]], y1 = delayFmns[[6]] +delayFses[[6]],length=0 ,col="gray"    )
	mtext(side=3,"flowering time")
# 	mtext(side=2,"Prb of observing Effect",line=2)
	axis(side=2,at=c(0,0.25,0.5,0.75,1),labels=c(0," ",0.5, " ",1))
legend(xvals[12],y=1.2,c("Earlier","Delay"),fill=c("black","gray"),bty="n")
plot(c(0,1)~I(c(0,1)),pch=NA,xlab="",ylab="",bty="n",xaxt="n",ylim=c(-0.2,1.2),yaxt="n")
	points(earlyGmns[[1]]~xvals[1:4]  ,pch=16:20); 	points(delayGmns[[1]]~I(xvals[1:4] +0.01 ),pch=16:20,col="gray")
	arrows(x0=xvals[1:4] , y0=earlyGmns[[1]] -earlyGses[[1]], y1 = earlyGmns[[1]] +earlyGses[[1]],length=0     )
	arrows(x0=xvals[1:4]+0.01 , y0=delayGmns[[1]] -delayGses[[1]], y1 = delayGmns[[1]] +delayGses[[1]],length=0  ,col="gray"   )
	abline(v= sum(xvals[4:5])/2,lty=3 )
	text(x=xvals[1:4] ,y=-0.2,c("b","a","c","bc"))
	text(x=xvals[1:4] ,y=1.2,c("b","c","a","a"),col="gray")
	points(earlyGmns[[2]]~xvals[5:7]  ,pch=16:20); 		points(delayGmns[[2]]~I(xvals[5:7] +0.01)  ,pch=16:20,col="gray")
	arrows(x0=xvals[5:7] , y0=earlyGmns[[2]] -earlyGses[[2]], y1 = earlyGmns[[2]] +earlyGses[[2]],length=0     )
	arrows(x0=xvals[5:7] + 0.01, y0=delayGmns[[2]] -delayGses[[2]], y1 = delayGmns[[2]] +delayGses[[2]],length=0 ,col="gray"    )
	text(x=xvals[5:7] ,y=-0.2,c(" "," "," "))
	text(x=xvals[5:7] ,y=1.2,c("a","b"," b"),col="gray")
	abline(v= sum(xvals[7:8])/2,lty=3 )
	points(earlyGmns[[3]]~xvals[8:9]  ,pch=16:20); 	points(delayGmns[[3]]~I(xvals[8:9] +0.01) ,pch=16:20, col="gray")
	arrows(x0=xvals[8:9] , y0=earlyGmns[[3]] -earlyGses[[3]], y1 = earlyGmns[[3]] +earlyGses[[3]],length=0     )
	arrows(x0=xvals[8:9]+0.01 , y0=delayGmns[[3]] -delayGses[[3]], y1 = delayGmns[[3]] +delayGses[[3]],length=0 ,col="gray"    )
	text(x=xvals[8:9] ,y=-0.2,c(" "," "))
	text(x=xvals[8:9] ,y=1.2,c(" "," "),col="gray")
	abline(v= sum(xvals[9:10])/2,lty=3 )
	points(earlyGmns[[4]]~xvals[10:12]  ,pch=16:20); 	points(delayGmns[[4]]~I(xvals[10:12]+0.01)  ,pch=16:21,col="gray")
	arrows(x0=xvals[10:12] , y0=earlyGmns[[4]] -earlyGses[[4]], y1 = earlyGmns[[4]] +earlyGses[[4]],length=0     )
	arrows(x0=xvals[10:12]+0.01 , y0=delayGmns[[4]] -delayGses[[4]], y1 = delayGmns[[4]] +delayGses[[4]],length=0 ,col="gray"    )
	text(x=xvals[10:12] ,y=-0.2,c("b","ab","a"))
	text(x=xvals[10:12] ,y=1.2,c("a","ab","b"),col="gray")
	abline(v= sum(xvals[12:13])/2,lty=3 )
	points(earlyGmns[[5]]~xvals[13:17]  ,pch=16:20); 	points(delayGmns[[5]]~I(xvals[13:17]+0.01)  ,pch=16:20,col="gray")
	arrows(x0=xvals[13:17] , y0=earlyGmns[[5]] -earlyGses[[5]], y1 = earlyGmns[[5]] +earlyGses[[5]],length=0     )
	arrows(x0=xvals[13:17] +0.01, y0=delayGmns[[5]] -delayGses[[5]], y1 = delayGmns[[5]] +delayGses[[5]],length=0 ,col="gray"    )
	text(x=xvals[13:17] ,y=-0.2,c(" "," "," "," "," "))
	text(x=xvals[13:17] ,y=1.2,c("a","a","b","a","b"),col="gray")
	abline(v= sum(xvals[17:18])/2,lty=3 )
	points(earlyGmns[[6]]~xvals[18]  ,pch=16:20); 	points(delayGmns[[6]]~I(xvals[18]+0.01)  ,pch=16:20,col="gray")
	arrows(x0=xvals[18] , y0=earlyGmns[[6]] -earlyGses[[6]], y1 = earlyGmns[[6]] +earlyGses[[6]],length=0     )
	arrows(x0=xvals[18] +0.01, y0=delayGmns[[6]] -delayGses[[6]], y1 = delayGmns[[6]] +delayGses[[6]],length=0 ,col="gray"    )
	mtext(side=3,"germination time")
	mtext(side=2,"Prb of observing Effect",line=2)
	axis(side=2,at=c(0,0.25,0.5,0.75,1),labels=c(0," ",0.5, " ",1))
legend(xvals[12],y=1.2,c("Earlier","Delay"),fill=c("black","gray"),bty="n")
plot(c(0,1)~I(c(0,1)),pch=NA,xlab="",ylab="",bty="n",xaxt="n",ylim=c(-0.2,1.2),yaxt="n")
	points(earlyGPmns[[1]]~xvals[1:4]  ,pch=16:20); 	points(delayGPmns[[1]]~I(xvals[1:4] +0.01 ),pch=16:20,col="gray")
	arrows(x0=xvals[1:4] , y0=earlyGPmns[[1]] -earlyGPses[[1]], y1 = earlyGPmns[[1]] +earlyGPses[[1]],length=0     )
	arrows(x0=xvals[1:4]+0.01 , y0=delayGPmns[[1]] -delayGPses[[1]], y1 = delayGPmns[[1]] +delayGPses[[1]],length=0  ,col="gray"   )
	abline(v= sum(xvals[4:5])/2,lty=3 )
	text(x=xvals[1:4] ,y=-0.2,c("b","c","a","a"))
	text(x=xvals[1:4] ,y=1.2,c("b","a","bc","c"),col="gray")
	points(earlyGPmns[[2]]~xvals[5:7]  ,pch=16:20); 		points(delayGPmns[[2]]~I(xvals[5:7] +0.01)  ,pch=16:20,col="gray")
	arrows(x0=xvals[5:7] , y0=earlyGPmns[[2]] -earlyGPses[[2]], y1 = earlyGPmns[[2]] +earlyGPses[[2]],length=0     )
	arrows(x0=xvals[5:7] + 0.01, y0=delayGPmns[[2]] -delayGPses[[2]], y1 = delayGPmns[[2]] +delayGPses[[2]],length=0 ,col="gray"    )
	text(x=xvals[5:7] ,y=-0.2,c("a","a","b")) 
	text(x=xvals[5:7] ,y=1.2,c("b","b","a"),col="gray")
	abline(v= sum(xvals[7:8])/2,lty=3 )
	points(earlyGPmns[[3]]~xvals[8:9]  ,pch=16:20); 	points(delayGPmns[[3]]~I(xvals[8:9] +0.01) ,pch=16:20, col="gray")
	arrows(x0=xvals[8:9] , y0=earlyGPmns[[3]] -earlyGPses[[3]], y1 = earlyGPmns[[3]] +earlyGPses[[3]],length=0     )
	arrows(x0=xvals[8:9]+0.01 , y0=delayGPmns[[3]] -delayGPses[[3]], y1 = delayGPmns[[3]] +delayGPses[[3]],length=0 ,col="gray"    )
	text(x=xvals[8:9] ,y=-0.2,c(" "," "))
	text(x=xvals[8:9] ,y=1.2,c("b","a"),col="gray")
	abline(v= sum(xvals[9:10])/2,lty=3 )
	points(earlyGPmns[[4]]~xvals[10:12]  ,pch=16:20); 	points(delayGPmns[[4]]~I(xvals[10:12]+0.01)  ,pch=16:21,col="gray")
	arrows(x0=xvals[10:12] , y0=earlyGPmns[[4]] -earlyGPses[[4]], y1 = earlyGPmns[[4]] +earlyGPses[[4]],length=0     )
	arrows(x0=xvals[10:12]+0.01 , y0=delayGPmns[[4]] -delayGPses[[4]], y1 = delayGPmns[[4]] +delayGPses[[4]],length=0 ,col="gray"    )
	text(x=xvals[10:12] ,y=-0.2,c(" "," "," "))
	text(x=xvals[10:12] ,y=1.2,c("b","ab","a"),col="gray")#THIS ONE ALSO LOOKS REVERSED
	abline(v= sum(xvals[12:13])/2,lty=3 )
	points(earlyGPmns[[5]]~xvals[13:17]  ,pch=16:20); 	points(delayGPmns[[5]]~I(xvals[13:17]+0.01)  ,pch=16:20,col="gray")
	arrows(x0=xvals[13:17] , y0=earlyGPmns[[5]] -earlyGPses[[5]], y1 = earlyGPmns[[5]] +earlyGPses[[5]],length=0     )
	arrows(x0=xvals[13:17] +0.01, y0=delayGPmns[[5]] -delayGPses[[5]], y1 = delayGPmns[[5]] +delayGPses[[5]],length=0 ,col="gray"    )
	text(x=xvals[13:17] ,y=-0.2,c("a","a","b","a","a"))
	text(x=xvals[13:17] ,y=1.2,c("a","a","b","a","a"),col="gray")
	abline(v= sum(xvals[17:18])/2,lty=3 )
	points(earlyGPmns[[6]]~xvals[18]  ,pch=16:20); 	points(delayGPmns[[6]]~I(xvals[18]+0.01)  ,pch=16:20,col="gray")
	arrows(x0=xvals[18] , y0=earlyGPmns[[6]] -earlyGPses[[6]], y1 = earlyGPmns[[6]] +earlyGPses[[6]],length=0     )
	arrows(x0=xvals[18] +0.01, y0=delayGPmns[[6]] -delayGPses[[6]], y1 = delayGPmns[[6]] +delayGPses[[6]],length=0 ,col="gray"    )
	mtext(side=3,"germination probability")
# 	mtext(side=2,"Prb of observing Effect",line=2)
	axis(side=2,at=c(0,0.25,0.5,0.75,1),labels=c(0," ",0.5, " ",1))
legend(xvals[12],y=1.2,c("Narrowed","Expanded"),fill=c("black","gray"),bty="n")
axis(side=1,at=xvals, labels=xlab,las=2)
dev.off()




flower.tax <- MCMCglmm(numeric.phen.eff~micrtax.sim2,random = ~dummyWOSrelrank,data=flowering,family="ordinal",verbose=F,prior=priornuup,pl=T,pr=T,nitt=100000,thin=50,burnin=1000)
flower.scom <- MCMCglmm(numeric.phen.eff~Strainvcomm.1,random = ~dummyWOSrelrank,data=flowering,family="ordinal",verbose=F,prior=priornuup,pl=T,pr=T,nitt=100000,thin=50,burnin=1000)
flower.mat <- MCMCglmm(numeric.phen.eff~mating2,random = ~dummyWOSrelrank,data=flowering,family="ordinal",verbose=F,prior=priornuup,pl=T,pr=T,nitt=100000,thin=50,burnin=1000)
flower.lh <- MCMCglmm(numeric.phen.eff~lifeform2,random = ~dummyWOSrelrank,data=flowering,family="ordinal",verbose=F,prior=priornuup,pl=T,pr=T,nitt=100000,thin=50,burnin=1000)
#even with nuup, there is still some autocorr in the posterior for one of the random effects, and still some correlation bt fixed and random effects, but MUCH less, and effective sampls is SO much better.
# flower.taxnr <- MCMCglmm(numeric.phen.eff~micrtax.sim2-1,random = ~dummyWOSrelrank,data=flowering,family="ordinal",verbose=F,prior=prior.m2c.4,pl=T,pr=T,nitt=100000,thin=50,burnin=1000)

germprb.tax <- MCMCglmm(numeric.phen.eff~micrtax.sim2,random = ~dummyWOSrelrank,data=germprb,family="ordinal",verbose=F,prior=priornuup,pl=T,pr=T,nitt=100000,thin=50,burnin=1000)
#eff samples issue for MF
germprb.scom <- MCMCglmm(numeric.phen.eff~Strainvcomm.1,random = ~dummyWOSrelrank,data=germprb,family="ordinal",verbose=F,prior=priornuup,pl=T,pr=T,nitt=100000,thin=50,burnin=1000)
#eff sample issue for strainmix
#.tax best by DIC; BUT ROUNDS FLUCTUATE A LOT for .scom & .tax;
germprb.mat <- MCMCglmm(numeric.phen.eff~mating2,random = ~dummyWOSrelrank,data=germprb,family="ordinal",verbose=F,prior=priornuup,pl=T,pr=T,nitt=100000,thin=50,burnin=1000)
germprb.lh <- MCMCglmm(numeric.phen.eff~lifeform2,random = ~dummyWOSrelrank,data=germprb,family="ordinal",verbose=F,prior=priornuup,pl=T,pr=T,nitt=100000,thin=50,burnin=1000)

germtime.tax <- MCMCglmm(numeric.phen.eff~micrtax.sim2,random = ~dummyWOSrelrank,data=germtime,family="ordinal",verbose=F,prior=priornuup,pl=T,pr=T,nitt=100000,thin=50,burnin=1000)
#eff samples issue for MF
germtime.scom <- MCMCglmm(numeric.phen.eff~Strainvcomm.1,random = ~dummyWOSrelrank,data=germtime,family="ordinal",verbose=F,prior=priornuup,pl=T,pr=T,nitt=100000,thin=50,burnin=1000)
#.tax best by DIC; but rounds fluctuate somewhat for DIC
germtime.mat <- MCMCglmm(numeric.phen.eff~mating2,random = ~dummyWOSrelrank,data=germtime,family="ordinal",verbose=F,prior=priornuup,pl=T,pr=T,nitt=100000,thin=50,burnin=1000)
germtime.lh <- MCMCglmm(numeric.phen.eff~lifeform2,random = ~dummyWOSrelrank,data=germtime,family="ordinal",verbose=F,prior=priornuup,pl=T,pr=T,nitt=100000,thin=50,burnin=1000)


mnpredbins <- function(int,other,cp){
	intpred <- c(pnorm(-int), pnorm(cp-int)-pnorm(-int), 1-pnorm(cp-int))
	otherplusint <- other + int
	otherpred <- t(sapply(otherplusint, function(p) c(pnorm(-p), pnorm(cp-p)-pnorm(-p), 1-pnorm(cp-p) ) ) )
	allpred <- data.frame(rbind(intpred, otherpred))
	colnames(allpred) <- c("P(k1|X)","P(k2|X)","P(k3|X)")
	return(allpred)
}


wb <- colorRampPalette(c(rgb(1,1,1),rgb(0,0,1)))

mtypes  <- c("bacteria","MF","mixed","otherfungi")
resptypes <- c(-1,0,1)
wosdf <- unique(flowering$dummyWOSrelrank)
fxmxrxw <- lapply(mtypes, function(m)  t(sapply(wosdf, function(w) sapply(resptypes, function(r) sum(flowering$dummyWOSrelrank==w &flowering$micrtax.sim2==m & flowering$numeric.phen.eff==r) )  )))
flwr.m <- t(sapply(fxmxrxw, function(m) colMeans(m/rowSums(m),na.rm=T) ))
gtxmxrxw <- lapply(mtypes, function(m)  t(sapply(wosdf, function(w) sapply(resptypes, function(r) sum(germtime$dummyWOSrelrank==w &germtime$micrtax.sim2==m & germtime$numeric.phen.eff==r) )  )))
germT.m <- t(sapply(gtxmxrxw, function(m) colMeans(m/rowSums(m),na.rm=T) ))
oxmxrxw <- lapply(mtypes, function(m)  t(sapply(wosdf, function(w) sapply(resptypes, function(r) sum(othertime$dummyWOSrelrank==w &othertime$micrtax.sim2==m & othertime$numeric.phen.eff==r) )  )))
othT.m <- t(sapply(oxmxrxw, function(m) colMeans(m/rowSums(m),na.rm=T) ))
gpxmxrxw <- lapply(mtypes, function(m)  t(sapply(wosdf, function(w) sapply(resptypes, function(r) sum(germprb$dummyWOSrelrank==w &germprb$micrtax.sim2==m & germprb$numeric.phen.eff==r) )  )))
germP.m <- t(sapply(gpxmxrxw, function(m) colMeans(m/rowSums(m),na.rm=T) ))
ttypes <- sort(unique(alltime$phen.trait.3))
axmxrxw <- lapply(mtypes, function(m) lapply(ttypes, function(trt)  t(sapply(wosdf, function(w) sapply(resptypes, function(r) sum(alltime$dummyWOSrelrank==w &alltime$micrtax.sim2==m & alltime$numeric.phen.eff==r & alltime$phen.trait.3==trt) )  )) ))
allT.m1 <- lapply(axmxrxw, function(m) rbind(m[[1]],m[[2]],m[[3]],m[[4]]))#depends on knowing the length of each element, which is the # of trait types
#unlist(axmxrxw,recursive=F)
allT.m <- t(sapply(allT.m1,function(mt) colMeans(mt/rowSums(mt),na.rm=T) ) )
# apxmxrxw <- lapply(mtypes, function(m)  t(sapply(wosdf, function(w) sapply(resptypes, function(r) sum(alltime$dummyWOSrelrank==w &alltime$micrtax.sim2==m & alltime$numeric.phen.eff==r) )  )))
# allT.m <- t(sapply(apxmxrxw, function(m) colMeans(m/rowSums(m),na.rm=T) ))
rxmxrxw <- lapply(mtypes, function(m)  t(sapply(wosdf, function(w) sapply(resptypes, function(r) sum(reprotime$dummyWOSrelrank==w &reprotime$micrtax.sim2==m & reprotime$numeric.phen.eff==r) )  )))
rprd.m <- t(sapply(rxmxrxw, function(m) colMeans(m/rowSums(m),na.rm=T) ))
vxmxrxw <- lapply(mtypes, function(m)  t(sapply(wosdf, function(w) sapply(resptypes, function(r) sum(vegtime$dummyWOSrelrank==w &vegtime$micrtax.sim2==m & vegtime$numeric.phen.eff==r) )  )))
veg.m <- t(sapply(vxmxrxw, function(m) colMeans(m/rowSums(m),na.rm=T) ))


##differences based on models
flwr.mgrps <- matrix(rep(c("a","b"),each=2),nrow=4,ncol=1,byrow=T)
germT.mgrps <- matrix(c(rep(c("b","a","c","c"),each=1)),nrow=4,ncol=1,byrow=T)
othT.mgrps <- matrix(rep(c("c","b","c","a"),times=1),nrow=4,ncol=1,byrow=T)
allT.mgrps <- matrix(rep(c("a","b"),each=2),nrow=4,ncol=1,byrow=T)
germP.mgrps <- matrix(rep(c("b","c","a","a"),each=1),nrow=4,ncol=1,byrow=T)
rpd.mgrps <- matrix(c(rep(c("b","a","c","abc"),each=1)),nrow=4,ncol=1,byrow=T)
veg.mgrps <- matrix(c(rep(c("b","a","a"," "),each=1)),nrow=4,ncol=1,byrow=T)

#check similarity of model predictions and results
mnpredbins(summary(flower.tax)$solutions[1,1],summary(flower.tax)$solutions[-1,1],summary(flower.tax)$cutpoints[1])
#prediction medium ok
mnpredbins(summary(germtime.tax)$solutions[1,1],summary(germtime.tax)$solutions[-1,1],summary(germtime.tax)$cutpoints[1])
#predictions for germination time  medium ok
mnpredbins(summary(othertime.tax)$solutions[1,1],summary(othertime.tax)$solutions[-1,1],summary(othertime.tax)$cutpoints[1])
#predictions for othertime a bit off
mnpredbins(summary(alltime.tax2r)$solutions[1,1],summary(alltime.tax2r)$solutions[-1,1],summary(alltime.tax2r)$cutpoints[1])
#predictions medium ok
mnpredbins(summary(germprb.tax)$solutions[1,1],summary(germprb.tax)$solutions[-1,1],summary(germprb.tax)$cutpoints[1])
#a bit off / pretty good
mnpredbins(summary(rprd.tax)$solutions[1,1],summary(rprd.tax)$solutions[-1,1],summary(rprd.tax)$cutpoints[1])
mnpredbins(summary(veg.tax)$solutions[1,1],summary(veg.tax)$solutions[-1,1],summary(veg.tax)$cutpoints[1])

pdf("Phenology effects of microbe types.pdf",width=6,height=3)
par(mfrow=c(1,4))
par(oma=c(0,5,0,1))
par(mar=c(5,0,3,0))
image(t(flwr.m), col=wb(100),zlim=c(0,1),xaxt="n",yaxt="n",ylab="",xlab="")
	axis(side = 1, at = c(0,0.5,1),labels=c("earlier","none","delayed"),las=2)
	axis(side = 2, at = c(0,0.33,0.66,1),labels=mtypes,las=2)
	text(rep(c(0,0.5,1),each=4),rep( c(0,0.33,0.66,1),times=3), round(as.vector(flwr.m),digits=2))
# 	text(rep(c(0,0.5,1)+.1,each=4),rep( c(0,0.33,0.66,1)+.1,times=3), flwr.mgrps,cex=0.75)
	text(-0.2,rep( c(0,0.33,0.66,1)+.1, times=1), flwr.mgrps,cex=0.75)
	mtext(side=3,"Flowering time")
image(t(germT.m), col=wb(100),zlim=c(0,1),xaxt="n",yaxt="n",ylab="",xlab="")
	axis(side = 1, at = c(0,0.5,1),labels=c("earlier","none","delayed"),las=2)
	text(rep(c(0,0.5,1),each=4),rep( c(0,0.33,0.66,1),times=3), round(as.vector(germT.m),digits=2))
# 	text(rep(c(0,0.5,1)+.1,each=4),rep( c(0,0.33,0.66,1)+.1,times=3), germT.mgrps,cex=0.75)
	text(-0.2,rep( c(0,0.33,0.66,1)+.1, times=1), germT.mgrps,cex=0.75)
	mtext(side=3,"Germination time")
image(t(othT.m), col=wb(100),zlim=c(0,1),xaxt="n",yaxt="n",ylab="",xlab="")
	axis(side = 1, at = c(0,0.5,1),labels=c("earlier","none","delayed"),las=2)
	text(rep(c(0,0.5,1),each=4),rep( c(0,0.33,0.66,1),times=3), round(as.vector(othT.m),digits=2))
# 	text(rep(c(0,0.5,1)+.1,each=4),rep( c(0,0.33,0.66,1)+.1,times=3), otherT.mgrps,cex=0.75)
	text(-0.2,rep( c(0,0.33,0.66,1)+.1, times=1), othT.mgrps,cex=0.75)
	mtext(side=3,"Other times")
image(t(allT.m), col=wb(100),zlim=c(0,1),xaxt="n",yaxt="n",ylab="",xlab="")
	axis(side = 1, at = c(0,0.5,1),labels=c("earlier","none","delayed"),las=2)
	text(rep(c(0,0.5,1),each=4),rep( c(0,0.33,0.66,1),times=3), round(as.vector(allT.m),digits=2))
	text(-0.2,rep( c(0,0.33,0.66,1)+.1, times=1), allT.mgrps,cex=0.75)
	mtext(side=3,"All times")
# 	text(rep(c(0,0.5,1)+.1,each=4),rep( c(0,0.33,0.66,1)+.1,times=3), allT.mgrps,cex=0.75)
dev.off()

pdf("Germination Prb effects of microbe types.pdf",width=3,height=3)
par(oma=c(0,5,0,1))
par(mar=c(5,0,3,0))
image(t(germP.m), col=wb(100),zlim=c(0,1),xaxt="n",yaxt="n",ylab="",xlab="")
	axis(side = 1, at = c(0,0.5,1),labels=c("narrowed","none","expanded"),las=2)
	axis(side = 2, at = c(0,0.33,0.66,1),labels=mtypes,las=2)
	text(rep(c(0,0.5,1),each=4),rep( c(0,0.33,0.66,1),times=3), round(as.vector(germP.m),digits=2))
# 	text(rep(c(0,0.5,1)+.1,each=4),rep( c(0,0.33,0.66,1)+.1,times=3), germP.mgrps,cex=0.75)
	text(-0.2,rep( c(0,0.33,0.66,1)+.1, times=1), germP.mgrps,cex=0.75)
	mtext(side=3,"Germination Prb")
dev.off()


ltypes  <- c("annual","both","perennial") #leaves out observed, but obs is only 1 datapoint; might consider whole record?
fxlxrxw <- lapply(ltypes, function(e)  t(sapply(wosdf, function(w) sapply(resptypes, function(r) sum(flowering$dummyWOSrelrank==w &flowering$lifeform2==e & flowering$numeric.phen.eff==r) )  )))
flwr.l <- t(sapply(fxlxrxw, function(e) colMeans(e/rowSums(e),na.rm=T) ))
gpxlxrxw <- lapply(ltypes, function(e)  t(sapply(wosdf, function(w) sapply(resptypes, function(r) sum(germprb$dummyWOSrelrank==w &germprb$lifeform2==e & germprb$numeric.phen.eff==r) )  )))
germP.l <- t(sapply(gpxlxrxw, function(e) colMeans(e/rowSums(e),na.rm=T) ))
gtxlxrxw <- lapply(ltypes, function(e)  t(sapply(wosdf, function(w) sapply(resptypes, function(r) sum(germtime$dummyWOSrelrank==w &germtime$lifeform2==e & germtime$numeric.phen.eff==r) )  )))
germT.l <- t(sapply(gtxlxrxw, function(e) colMeans(e/rowSums(e),na.rm=T) ))
flwr.lgrps <- matrix(c("a","a","b"),nrow=3,ncol=1,byrow=T)
germT.lgrps <- matrix(c("b","a","a"),nrow=3,ncol=1,byrow=T)
germP.lgrps <- matrix(c("a","a","b"),nrow=3,ncol=1,byrow=T)
pdf("Phenology effects of lifehistory types.pdf",width=6,height=3) #general weirdness about pred vs data 
par(mfrow=c(1,3))
par(oma=c(0,5,0,1))
par(mar=c(5,0,3,0))
image(t(flwr.l), col=wb(100),zlim=c(0,1),xaxt="n",yaxt="n",ylab="",xlab="")
	axis(side = 1, at = c(0,0.5,1),labels=c("earlier","none","delayed"),las=2)
	axis(side = 2, at = c(0,0.5,1),labels=ltypes,las=2)
	text(rep(c(0,0.5,1),each=3),rep( c(0,0.5,1),times=3), round(as.vector(flwr.l),digits=2))
	text(-0.2,rep( c(0,0.5,1)+.1, times=1), flwr.lgrps,cex=0.75)
# 	text(1.1,1.1,"*" )
	mtext(side=3,"Flowering time")
image(t(germT.l), col=wb(100),zlim=c(0,1),xaxt="n",yaxt="n",ylab="",xlab="")
	axis(side = 1, at = c(0,0.5,1),labels=c("earlier","none","delayed"),las=2)
	text(rep(c(0,0.5,1),each=3),rep( c(0,0.5,1),times=3), round(as.vector(germT.l),digits=2))
	text(-0.2,rep( c(0,0.5,1)+.1, times=1), germT.lgrps,cex=0.75)
	mtext(side=3,"Germination time")
image(t(germP.l), col=wb(100),zlim=c(0,1),xaxt="n",yaxt="n",ylab="",xlab="")
	axis(side = 1, at = c(0,0.5,1),labels=c("narrowed","none","expanded"),las=2)
	text(rep(c(0,0.5,1),each=3),rep( c(0,0.5,1),times=3), round(as.vector(germP.l),digits=2))
 	text(-0.2,rep( c(0,0.5,1)+.1, times=1), germP.lgrps,cex=0.75)
	mtext(side=3,"Germ prb")
dev.off()


matypes  <- c("outcross","s+o","selfer","unk") #leaves out observed, but obs is only 1 datapoint; might consider whole record?
gpxmaxrxw <- lapply(matypes, function(e)  t(sapply(wosdf, function(w) sapply(resptypes, function(r) sum(germprb$dummyWOSrelrank==w &germprb$mating2==e & germprb$numeric.phen.eff==r) )  )))
germP.ma <- t(sapply(gpxmaxrxw, function(e) colMeans(e/rowSums(e),na.rm=T) ))
germP.magrps <- matrix(c("a","b","a"),nrow=3,ncol=1,byrow=T)
pdf("Germprb mating types.pdf",width=3,height=3) #general weirdness about pred vs data 
par(oma=c(0,5,0,1))
par(mar=c(5,0,3,0))
image(t(germP.ma), col=wb(100),zlim=c(0,1),xaxt="n",yaxt="n",ylab="",xlab="")
	axis(side = 1, at = c(0,0.5,1),labels=c("narrowed","none","expanded"),las=2)
	axis(side = 2, at = c(0,0.33,0.66,1),labels=matypes,las=2)
	text(rep(c(0,0.5,1),each=4),rep( c(0,0.33,0.66,1),times=3), round(as.vector(germP.ma),digits=2))
 	text(-0.2,rep( c(0,0.5,1)+.1, times=1), germP.magrps,cex=0.75)
	mtext(side=3,"Germ prb")
dev.off()

itypes  <- c("community","strain","strain mix") #leaves out observed, but obs is only 1 datapoint; might consider whole record?
resptypes <- c(-1,0,1)
wosdf <- unique(flowering$dummyWOSrelrank)
fxixrxw <- lapply(itypes, function(i)  t(sapply(wosdf, function(w) sapply(resptypes, function(r) sum(flowering$dummyWOSrelrank==w &flowering$Strainvcomm.1==i & flowering$numeric.phen.eff==r) )  )))
flwr.i <- t(sapply(fxixrxw, function(i) colMeans(i/rowSums(i),na.rm=T) ))
gtxixrxw <- lapply(itypes, function(i)  t(sapply(wosdf, function(w) sapply(resptypes, function(r) sum(germtime$dummyWOSrelrank==w &germtime$Strainvcomm.1==i & germtime$numeric.phen.eff==r) )  )))
germT.i <- t(sapply(gtxixrxw, function(i) colMeans(i/rowSums(i),na.rm=T) ))
oxixrxw <- lapply(itypes, function(i)  t(sapply(wosdf, function(w) sapply(resptypes, function(r) sum(othertime$dummyWOSrelrank==w &othertime$Strainvcomm.1==i & othertime$numeric.phen.eff==r) )  )))
othT.i <- t(sapply(oxixrxw, function(i) colMeans(i/rowSums(i),na.rm=T) ))
ttypes <- sort(unique(alltime$phen.trait.3))
axixrxw <- lapply(itypes, function(i) lapply(ttypes, function(trt)  t(sapply(wosdf, function(w) sapply(resptypes, function(r) sum(alltime$dummyWOSrelrank==w &alltime$Strainvcomm.1==i & alltime$numeric.phen.eff==r & alltime$phen.trait.3==trt) )  )) ))
allT.i1 <- lapply(axixrxw, function(i) rbind(i[[1]],i[[2]],i[[3]],i[[4]]))#depends on knowing the length of each element, which is the # of trait types
allT.i <- t(sapply(allT.i1,function(it) colMeans(it/rowSums(it),na.rm=T) ) )
gpxixrxw <- lapply(itypes, function(i)  t(sapply(wosdf, function(w) sapply(resptypes, function(r) sum(germprb$dummyWOSrelrank==w &germprb$Strainvcomm.1==i & germprb$numeric.phen.eff==r) )  )))
germP.i <- t(sapply(gpxixrxw, function(i) colMeans(i/rowSums(i),na.rm=T) ))
flwr.igrps <- matrix(rep(c("a","b","b"),each=1),nrow=3,ncol=1,byrow=T)
germT.igrps <- matrix(c(rep(c("a","a","a"),each=1)),nrow=3,ncol=1,byrow=T)
othT.igrps <- matrix(rep(c("b","a","ab"),times=1),nrow=3,ncol=1,byrow=T)
allT.igrps <- matrix(rep(c("a","b","b"),each=2),nrow=3,ncol=1,byrow=T)
germP.igrps <- matrix(rep(c("a","a","b"),each=1),nrow=3,ncol=1,byrow=T)
mnpredbins(summary(flower.scom)$solutions[1,1],summary(flower.scom)$solutions[-1,1],summary(flower.scom)$cutpoints[1])
##PREDICTS OPPOSITE DIRECTION....
mnpredbins(summary(germtime.scom)$solutions[1,1],summary(germtime.scom)$solutions[-1,1],summary(germtime.scom)$cutpoints[1])
mnpredbins(summary(othertime.scom)$solutions[1,1],summary(othertime.scom)$solutions[-1,1],summary(othertime.scom)$cutpoints[1])
mnpredbins(summary(alltime.scom2r)$solutions[1,1],summary(alltime.scom2r)$solutions[-1,1],summary(alltime.scom2r)$cutpoints[1])
mnpredbins(summary(germprb.scom)$solutions[1,1],summary(germprb.scom)$solutions[-1,1],summary(germprb.scom)$cutpoints[1])


pdf("Phenology effects of inoculation types.pdf",width=6,height=3) #again weirdness on some (flowering esp, alltime to lesser extent), bc of same problem as for strains
par(mfrow=c(1,4))
par(oma=c(0,5,0,1))
par(mar=c(5,0,3,0))
image(t(flwr.i), col=wb(100),zlim=c(0,1),xaxt="n",yaxt="n",ylab="",xlab="")
	axis(side = 1, at = c(0,0.5,1),labels=c("earlier","none","delayed"),las=2)
	axis(side = 2, at = c(0,0.5,1),labels=itypes,las=2)
	text(rep(c(0,0.5,1),each=3),rep( c(0,0.5,1),times=3), round(as.vector(flwr.i),digits=2))
	text(-0.2,rep( c(0,0.5,1)+.1, times=1), flwr.igrps,cex=0.75)
	mtext(side=3,"Flowering time")#these model results are  VERY odd
image(t(germT.i), col=wb(100),zlim=c(0,1),xaxt="n",yaxt="n",ylab="",xlab="")
	axis(side = 1, at = c(0,0.5,1),labels=c("earlier","none","delayed"),las=2)
	text(rep(c(0,0.5,1),each=3),rep( c(0,0.5,1),times=3), round(as.vector(germT.i),digits=2))
	text(-0.2,rep( c(0,0.5,1)+.1, times=1), germT.igrps,cex=0.75)
	mtext(side=3,"Germination time")
image(t(othT.i), col=wb(100),zlim=c(0,1),xaxt="n",yaxt="n",ylab="",xlab="")
	axis(side = 1, at = c(0,0.5,1),labels=c("earlier","none","delayed"),las=2)
	text(rep(c(0,0.5,1),each=3),rep( c(0,0.5,1),times=3), round(as.vector(othT.i),digits=2))
	text(-0.2,rep( c(0,0.5,1)+.1, times=1), othT.igrps,cex=0.75)
	mtext(side=3,"Other times")
image(t(allT.i), col=wb(100),zlim=c(0,1),xaxt="n",yaxt="n",ylab="",xlab="")
	axis(side = 1, at = c(0,0.5,1),labels=c("earlier","none","delayed"),las=2)
	text(rep(c(0,0.5,1),each=3),rep( c(0,0.5,1),times=3), round(as.vector(allT.i),digits=2))
	mtext(side=3,"All times")
	text(-0.2,rep( c(0,0.5,1)+.1, times=1), allT.igrps,cex=0.75)
dev.off()


pdf("Germination Prb effects of inoculation types.pdf",width=3,height=3)
par(oma=c(0,5,0,1))
par(mar=c(5,0,3,0))
image(t(germP.i), col=wb(100),zlim=c(0,1),xaxt="n",yaxt="n",ylab="",xlab="")
	axis(side = 1, at = c(0,0.5,1),labels=c("narrowed","none","expanded"),las=2)
	axis(side = 2, at = c(0,0.5,1),labels=itypes,las=2)
	text(rep(c(0,0.5,1),each=3),rep( c(0,0.5,1),times=3), round(as.vector(germP.i),digits=2))
	text(-0.2,rep( c(0,0.5,1)+.1, times=1), germP.igrps,cex=0.75)
	mtext(side=3,"Germination Prb")
dev.off()





####BINOMIAL 
# 	abline(v= sum(xvals[19:20])/2,lty=3 )
# 	points(earlyFmns[[6]]~xvals[20:22]  ,pch=16:20); 	points(delayFmns[[6]]~I(xvals[20:22]+ 0.01)  ,pch=16:20, col="gray")
# 	arrows(x0=xvals[20:22] , y0=earlyFmns[[6]] -earlyFses[[6]], y1 = earlyFmns[[6]] +earlyFses[[6]],length=0     )
# 	arrows(x0=xvals[20:22]+0.01 , y0=delayFmns[[6]] -delayFses[[6]], y1 = delayFmns[[6]] +delayFses[[6]],length=0 ,col="gray"    )
# 	text(x=xvals[20:22] ,y=-0.2,c("a"," ","b"))
# 	text(x=xvals[20:22] ,y=1.2,c(" "," "," "),col="gray")

##### RESPONSE TRAITS NO LONGER OF INTEReST
# othertime <- split_recs_scr[split_recs_scr$phen.trait.2%in% 
# 				c("budburst", "floral budset time", "flower senescence time", "fruiting time","maturation time" ,  "peak flowering" , "phyllochron","senescence time") ,]
# 	othertime$dummyWOSrelrank <- othertime$WOS.relevance.rank
# 	othertime$dummyWOSrelrank[othertime$dummyWOSrelran%in%names(table(othertime$WOS.relevance.rank)[table(othertime$WOS.relevance.rank)<13])] <- "binned"
# 	othertime<- othertime[-which(othertime$culture.sim=="observed"),]#removing the sole record that is not manipulated
# alltime <- split_recs_scr[split_recs_scr$phen.trait.2%in% 
# 				c("germination time","budburst", "floral budset time", "flowering time", "flower senescence time", "fruiting time","maturation time" ,  "peak flowering" , "phyllochron","senescence time") ,]
# 	alltime$dummyWOSrelrank <- alltime$WOS.relevance.rank
# 	alltime$dummyWOSrelrank[alltime$dummyWOSrelran%in%names(table(alltime$WOS.relevance.rank)[table(alltime$WOS.relevance.rank)<13])] <- "binned"
# 	alltime<- alltime[-which(alltime$culture.sim=="observed"),]##removing the sole record that is not manipulated
# reprotime <- split_recs_scr[split_recs_scr$phen.trait.3=="reproductive time",]
# 	reprotime$dummyWOSrelrank <- reprotime$WOS.relevance.rank
# 	reprotime$dummyWOSrelrank[reprotime$dummyWOSrelran%in%names(table(reprotime$WOS.relevance.rank)[table(reprotime$WOS.relevance.rank)<13])] <- "binned"
# vegtime <- split_recs_scr[split_recs_scr$phen.trait.3=="vegetative time",]
# 	vegtime$dummyWOSrelrank <- vegtime$WOS.relevance.rank
# 	vegtime$dummyWOSrelrank[vegtime$dummyWOSrelran%in%names(table(vegtime$WOS.relevance.rank)[table(vegtime$WOS.relevance.rank)<13])] <- "binned"
# 	vegtime<- vegtime[-which(vegtime$culture.sim=="observed"),]#removing the sole record that is not manipulated
# EarlyresF <- lapply(list(flowerbitax,flowerbitax15,flowerbicltr,flowerbiscom,flowerbimat,flowerbilh) , function(z) get.mnpred.binom(z,0.95))


#####CATEGORIES OR RESPONSE TRAITS NO LONGER OF INTEReST
# flowerDbitax15 <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~micrtax.sim1.5,random = ~dummyWOSrelrank,data=flowering,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
# flowerDbicltr <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~culture.sim,random = ~dummyWOSrelrank,data=flowering,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
# # EarlyresF <- lapply(list(flowerbitax,flowerbitax15,flowerbicltr,flowerbiscom,flowerbimat,flowerbilh) , function(z) get.mnpred.binom(z,0.95))
# flowerbitax15 <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~micrtax.sim1.5,random = ~dummyWOSrelrank,data=flowering,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
# flowerbicltr <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~culture.sim,random = ~dummyWOSrelrank,data=flowering,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
# germtimebitax15 <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~micrtax.sim1.5,random = ~dummyWOSrelrank,data=germtime,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
# germtimebicltr <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~culture.sim,random = ~dummyWOSrelrank,data=germtime,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
# germtimeDbitax15 <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~micrtax.sim1.5,random = ~dummyWOSrelrank,data=germtime,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
# germtimeDbicltr <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~culture.sim,random = ~dummyWOSrelrank,data=germtime,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
# 
# 
# alltimebitax <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~micrtax.sim2,random = ~dummyWOSrelrank + phen.trait.3,data=alltime,family="multinomial2",prior=priornuup2,verbose=F,pr=T,nitt=100000,thin=10,burnin=1000)
# alltimebitax15 <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~micrtax.sim1.5,random = ~dummyWOSrelrank + phen.trait.3,data=alltime,family="multinomial2",prior=priornuup2,verbose=F,pr=T,nitt=100000,thin=10,burnin=1000)
# alltimebicltr <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~culture.sim,random = ~dummyWOSrelrank + phen.trait.3,data=alltime,family="multinomial2",prior=priornuup2,verbose=F,pr=T,nitt=100000,thin=10,burnin=1000)
# alltimebiscom <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~Strainvcomm.1,random = ~dummyWOSrelrank + phen.trait.3,data=alltime,family="multinomial2",prior=priornuup2,verbose=F,pr=T,nitt=100000,thin=10,burnin=1000)
# alltimebimat <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~mating2,random = ~dummyWOSrelrank + phen.trait.3,data=alltime,family="multinomial2",prior=priornuup2,verbose=F,pr=T,nitt=100000,thin=10,burnin=1000)
# alltimebilh <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~lifeform2,random = ~dummyWOSrelrank + phen.trait.3,data=alltime,family="multinomial2",prior=priornuup2,verbose=F,pr=T,nitt=100000,thin=10,burnin=1000)
# #more likely to be earlier:bacteria vs fungi / mixed -- unless compl tax in which case MF and rhiz earlier; culture inoc vs removed exps; strains & strain mix vs mixed; ; no lh no mat
# EarlyresA <- lapply(list(alltimebitax,alltimebitax15,alltimebicltr,alltimebiscom,alltimebimat,alltimebilh) , function(z) get.mnpred.binom(z,0.95))
# earlyAmns <- lapply(c("micrtax.sim2","micrtax.sim1.5","culture.sim","Strainvcomm.1","mating2","lifeform2"),
# 		 function(z) tapply(alltime$isEarlyOrNrw,alltime[,colnames(alltime)==z],mean) )
# earlyAses <- lapply(c("micrtax.sim2","micrtax.sim1.5","culture.sim","Strainvcomm.1","mating2","lifeform2"),
# 		 function(z) tapply(alltime$isEarlyOrNrw,alltime[,colnames(alltime)==z],std.error) )


# get.mnpred.binom <- function(model,prb){
# 	Sol <- model$Sol
# 	ncats <- nrow(summary(model)$solutions)
# 	int <- c(HPDi(invlogistic(Sol[,1]),prob=prb)[1], mean(invlogistic(Sol[,1])), HPDi(invlogistic(Sol[,1]),prob=prb)[2])
# 	betacats <- 2:ncats
# 	betahpdi <- sapply(betacats, function(b) c( HPDi(invlogistic( Sol[,1]+Sol[,b] ),prob=prb)[1],  
# 												mean(invlogistic( Sol[,1]+Sol[,b] )), 
# 												HPDi(invlogistic( Sol[,1]+Sol[,b] ),prob=prb)[2]  ) )
# 	Expprb <- rbind(int, t(betahpdi))
# 	rownames(Expprb) <- colnames(Sol)[1:ncats]
# 	return(Expprb )
# }

# get.mnpred.binomni <- function(model,prb){
# 	Sol <- model$Sol
# 	betacats <- 1:nrow(summary(model)$solutions)
# 	Expprb <- t(sapply(betacats, function(b) c( HPDi(invlogistic( Sol[,b] ),prob=prb)[1],  
# 												mean(invlogistic( Sol[,b] )), 
# 												HPDi(invlogistic( Sol[,b] ),prob=prb)[2]  ) ) )
# 	rownames(Expprb) <- colnames(Sol)[1:nrow(summary(model)$solutions)]
# 	return(Expprb )
# } #tried this once, model and strange overlap of sig. diff predictions was the same.





### NO LONGER INTERESTED IN MICR15/ CULTURE groupings
# germprb.tax15 <- MCMCglmm(numeric.phen.eff~micrtax.sim1.5,random = ~dummyWOSrelrank,data=germprb,family="ordinal",verbose=F,prior=priornuup,pl=T,pr=T,nitt=100000,thin=50,burnin=1000)
# germprb.cltr <- MCMCglmm(numeric.phen.eff~culture.sim,random = ~dummyWOSrelrank,data=germprb,family="ordinal",verbose=F,prior=priornuup,pl=T,pr=T,nitt=100000,thin=50,burnin=1000)
# germtime.tax15 <- MCMCglmm(numeric.phen.eff~micrtax.sim1.5,random = ~dummyWOSrelrank,data=germtime,family="ordinal",verbose=F,prior=priornuup,pl=T,pr=T,nitt=100000,thin=50,burnin=1000)
# ##no diff bt bacts
# germtime.cltr <- MCMCglmm(numeric.phen.eff~culture.sim,random = ~dummyWOSrelrank,data=germtime,family="ordinal",verbose=F,prior=priornuup,pl=T,pr=T,nitt=100000,thin=50,burnin=1000)
###eff samp issue in wosrelrank rand effect


#### NO LONGER INTERESTED IN CLASSIFYING PHENOTYPES THIS WAY
# pdf("Phenology effects of microbe types - diff trait bins.pdf",width=6,height=3)
# par(mfrow=c(1,4))
# par(oma=c(0,5,0,1))
# par(mar=c(5,0,3,0))
# image(t(rprd.m), col=wb(100),zlim=c(0,1),xaxt="n",yaxt="n",ylab="",xlab="")
# 	axis(side = 1, at = c(0,0.5,1),labels=c("earlier","none","delayed"),las=2)
# 	axis(side = 2, at = c(0,0.33,0.66,1),labels=mtypes,las=2)
# 	text(rep(c(0,0.5,1),each=4),rep( c(0,0.33,0.66,1),times=3), round(as.vector(rprd.m),digits=2))
# # 	text(rep(c(0,0.5,1)+.1,each=4),rep( c(0,0.33,0.66,1)+.1,times=3), flwr.mgrps,cex=0.75)
# 	text(-0.2,rep( c(0,0.33,0.66,1)+.1, times=1), rpd.mgrps,cex=0.75)
# 	mtext(side=3,"Reproductive time")
# image(t(germT.m), col=wb(100),zlim=c(0,1),xaxt="n",yaxt="n",ylab="",xlab="")
# 	axis(side = 1, at = c(0,0.5,1),labels=c("earlier","none","delayed"),las=2)
# 	text(rep(c(0,0.5,1),each=4),rep( c(0,0.33,0.66,1),times=3), round(as.vector(germT.m),digits=2))
# # 	text(rep(c(0,0.5,1)+.1,each=4),rep( c(0,0.33,0.66,1)+.1,times=3), germT.mgrps,cex=0.75)
# 	text(-0.2,rep( c(0,0.33,0.66,1)+.1, times=1), germT.mgrps,cex=0.75)
# 	mtext(side=3,"Germination time")
# image(t(veg.m), col=wb(100),zlim=c(0,1),xaxt="n",yaxt="n",ylab="",xlab="")
# 	axis(side = 1, at = c(0,0.5,1),labels=c("earlier","none","delayed"),las=2)
# 	text(rep(c(0,0.5,1),each=4),rep( c(0,0.33,0.66,1),times=3), round(as.vector(veg.m),digits=2))
# # 	text(rep(c(0,0.5,1)+.1,each=4),rep( c(0,0.33,0.66,1)+.1,times=3), otherT.mgrps,cex=0.75)
# 	text(-0.2,rep( c(0,0.33,0.66,1)+.1, times=1), veg.mgrps,cex=0.75)
# 	mtext(side=3,"Vegetative times")
# image(t(allT.m), col=wb(100),zlim=c(0,1),xaxt="n",yaxt="n",ylab="",xlab="")
# 	axis(side = 1, at = c(0,0.5,1),labels=c("earlier","none","delayed"),las=2)
# 	text(rep(c(0,0.5,1),each=4),rep( c(0,0.33,0.66,1),times=3), round(as.vector(allT.m),digits=2))
# 	text(-0.2,rep( c(0,0.33,0.66,1)+.1, times=1), allT.mgrps,cex=0.75)
# 	mtext(side=3,"All times (+snsc)")
# # 	text(rep(c(0,0.5,1)+.1,each=4),rep( c(0,0.33,0.66,1)+.1,times=3), allT.mgrps,cex=0.75)
# dev.off()

####NO LONGER INTERESTED IN CLASSIFYING MICR THIS WAY
# mtypes  <- c("bacteria,other","MF","mixed","otherfungi", "rhizobia+")
# resptypes <- c(-1,0,1)
# wosdf <- unique(flowering$dummyWOSrelrank)
# fxmxrxw15 <- lapply(mtypes, function(m)  t(sapply(wosdf, function(w) sapply(resptypes, function(r) sum(flowering$dummyWOSrelrank==w &flowering$micrtax.sim1.5==m & flowering$numeric.phen.eff==r) )  )))
# flwr.m15 <- t(sapply(fxmxrxw15, function(m) colMeans(m/rowSums(m),na.rm=T) ))
# gtxmxrxw15 <- lapply(mtypes, function(m)  t(sapply(wosdf, function(w) sapply(resptypes, function(r) sum(germtime$dummyWOSrelrank==w &germtime$micrtax.sim1.5==m & germtime$numeric.phen.eff==r) )  )))
# germT.m15 <- t(sapply(gtxmxrxw15, function(m) colMeans(m/rowSums(m),na.rm=T) ))
# oxmxrxw15 <- lapply(mtypes, function(m)  t(sapply(wosdf, function(w) sapply(resptypes, function(r) sum(othertime$dummyWOSrelrank==w &othertime$micrtax.sim1.5==m & othertime$numeric.phen.eff==r) )  )))
# othT.m15 <- t(sapply(oxmxrxw15, function(m) colMeans(m/rowSums(m),na.rm=T) ))
# gpxmxrxw15 <- lapply(mtypes, function(m)  t(sapply(wosdf, function(w) sapply(resptypes, function(r) sum(germprb$dummyWOSrelrank==w &germprb$micrtax.sim1.5==m & germprb$numeric.phen.eff==r) )  )))
# germP.m15 <- t(sapply(gpxmxrxw15, function(m) colMeans(m/rowSums(m),na.rm=T) ))
# ttypes <- sort(unique(alltime$phen.trait.3))
# axmxrxw15 <- lapply(mtypes, function(m) lapply(ttypes, function(trt)  t(sapply(wosdf, function(w) sapply(resptypes, function(r) sum(alltime$dummyWOSrelrank==w &alltime$micrtax.sim1.5==m & alltime$numeric.phen.eff==r & alltime$phen.trait.3==trt) )  )) ))
# allT.m115 <- lapply(axmxrxw15, function(m) rbind(m[[1]],m[[2]],m[[3]],m[[4]]))#depends on knowing the length of each element, which is the # of trait types
# allT.m15 <- t(sapply(allT.m115,function(mt) colMeans(mt/rowSums(mt),na.rm=T) ) )
# ##differences based on models
# flwr.mgrps15 <- matrix(c("a","ab","b","b","ab"),nrow=5,ncol=1,byrow=T)
# germT.mgrps15 <- matrix(c("b","a","c","c","bc"),nrow=5,ncol=1,byrow=T)
# othT.mgrps15 <- matrix(c("c","b","bc","a","bc"),nrow=5,ncol=1,byrow=T)
# allT.mgrps15 <- matrix(c("c","ab","d","d","a"),nrow=5,ncol=1,byrow=T)
# germP.mgrps15 <- matrix(c("b","c","a","a"," "),nrow=5,ncol=1,byrow=T)
# #check similarity of model predictions and results
# mnpredbins(summary(flower.tax15)$solutions[1,1],summary(flower.tax15)$solutions[-1,1],summary(flower.tax15)$cutpoints[1])
# mnpredbins(summary(germtime.tax15)$solutions[1,1],summary(germtime.tax15)$solutions[-1,1],summary(germtime.tax15)$cutpoints[1])
# mnpredbins(summary(othertime.tax15)$solutions[1,1],summary(othertime.tax15)$solutions[-1,1],summary(othertime.tax15)$cutpoints[1])
# mnpredbins(summary(alltime.tax2r15)$solutions[1,1],summary(alltime.tax2r15)$solutions[-1,1],summary(alltime.tax2r15)$cutpoints[1])
# mnpredbins(summary(germprb.tax15)$solutions[1,1],summary(germprb.tax15)$solutions[-1,1],summary(germprb.tax15)$cutpoints[1])
# 
# pdf("Phenology effects of microbe types ammended.pdf",width=6,height=3)
# par(mfrow=c(1,4))
# par(oma=c(0,5,0,1))
# par(mar=c(5,0,3,0))
# image(t(flwr.m15), col=wb(100),zlim=c(0,1),xaxt="n",yaxt="n",ylab="",xlab="")
# 	axis(side = 1, at = c(0,0.5,1),labels=c("earlier","none","delayed"),las=2)
# 	axis(side = 2, at = c(0,0.25,0.5,0.75,1),labels=mtypes,las=2)
# 	text(rep(c(0,0.5,1),each=5),rep( c(0,0.25,0.5,0.75,1),times=3), round(as.vector(flwr.m15),digits=2))
# 	text(-0.2,rep( c(0,0.25,0.5,0.75,1)+.1, times=1), flwr.mgrps15,cex=0.75)
# 	mtext(side=3,"Flowering time")
# image(t(germT.m15), col=wb(100),zlim=c(0,1),xaxt="n",yaxt="n",ylab="",xlab="")
# 	axis(side = 1, at = c(0,0.5,1),labels=c("earlier","none","delayed"),las=2)
# 	text(rep(c(0,0.5,1),each=5),rep( c(0,0.25,0.5,0.75,1),times=3), round(as.vector(germT.m15),digits=2))
# 	text(-0.2,rep( c(0,0.25,0.5,0.75,1)+.1, times=1), germT.mgrps15,cex=0.75)
# 	mtext(side=3,"Germination time")
# image(t(othT.m15), col=wb(100),zlim=c(0,1),xaxt="n",yaxt="n",ylab="",xlab="")
# 	axis(side = 1, at = c(0,0.5,1),labels=c("earlier","none","delayed"),las=2)
# 	text(rep(c(0,0.5,1),each=5),rep( c(0,0.25,0.5,0.75,1),times=3), round(as.vector(othT.m15),digits=2))
# 	text(-0.2,rep( c(0,0.25,0.5,0.75,1)+.1, times=1), othT.mgrps15,cex=0.75)
# 	mtext(side=3,"Other times")
# image(t(allT.m15), col=wb(100),zlim=c(0,1),xaxt="n",yaxt="n",ylab="",xlab="")
# 	axis(side = 1, at = c(0,0.5,1),labels=c("earlier","none","delayed"),las=2)
# 	text(rep(c(0,0.5,1),each=5),rep( c(0,0.25,0.5,0.75,1),times=3), round(as.vector(allT.m15),digits=2))
# 	text(-0.2,rep( c(0,0.25,0.5,0.75,1)+.1, times=1), allT.mgrps15,cex=0.75)
# 	mtext(side=3,"All times")
# dev.off()
# pdf("Germination Prb effects of microbe types ammended.pdf",width=3,height=3)
# par(oma=c(0,5,0,1))
# par(mar=c(5,0,3,0))
# image(t(germP.m15), col=wb(100),zlim=c(0,1),xaxt="n",yaxt="n",ylab="",xlab="")
# 	axis(side = 1, at = c(0,0.5,1),labels=c("narrowed","none","expanded"),las=2)
# 	axis(side = 2, at = c(0,0.25,0.5,0.75,1),labels=mtypes,las=2)
# 	text(rep(c(0,0.5,1),each=5),rep( c(0,0.25,0.5,0.75,1),times=3), round(as.vector(germP.m15),digits=2))
# 	text(-0.2,rep( c(0,0.25,0.5,0.75,1)+.1, times=1), germP.mgrps15,cex=0.75)
# 	mtext(side=3,"Germination Prb")
# dev.off()



###ORDINAL PROBLEM AND no longer interested in these explanatory variables
# etypes  <- c("culture","direct","removed") #leaves out observed, but obs is only 1 datapoint; might consider whole record?
# resptypes <- c(-1,0,1)
# wosdf <- unique(flowering$dummyWOSrelrank)
# fxexrxw <- lapply(etypes, function(e)  t(sapply(wosdf, function(w) sapply(resptypes, function(r) sum(flowering$dummyWOSrelrank==w &flowering$culture.sim==e & flowering$numeric.phen.eff==r) )  )))
# flwr.e <- t(sapply(fxexrxw, function(e) colMeans(e/rowSums(e),na.rm=T) ))
# gtxexrxw <- lapply(etypes, function(e)  t(sapply(wosdf, function(w) sapply(resptypes, function(r) sum(germtime$dummyWOSrelrank==w &germtime$culture.sim==e & germtime$numeric.phen.eff==r) )  )))
# germT.e <- t(sapply(gtxexrxw, function(e) colMeans(e/rowSums(e),na.rm=T) ))
# oxexrxw <- lapply(etypes, function(e)  t(sapply(wosdf, function(w) sapply(resptypes, function(r) sum(othertime$dummyWOSrelrank==w &othertime$culture.sim==e & othertime$numeric.phen.eff==r) )  )))
# othT.e <- t(sapply(oxexrxw, function(e) colMeans(e/rowSums(e),na.rm=T) ))
# ttypes <- sort(unique(alltime$phen.trait.3))
# axexrxw <- lapply(etypes, function(e) lapply(ttypes, function(trt)  t(sapply(wosdf, function(w) sapply(resptypes, function(r) sum(alltime$dummyWOSrelrank==w &alltime$culture.sim==e & alltime$numeric.phen.eff==r & alltime$phen.trait.3==trt) )  )) ))
# allT.e1 <- lapply(axexrxw, function(e) rbind(e[[1]],e[[2]],e[[3]],e[[4]]))#depends on knowing the length of each element, which is the # of trait types
# allT.e <- t(sapply(allT.e1,function(et) colMeans(et/rowSums(et),na.rm=T) ) )
# gpxexrxw <- lapply(etypes, function(e)  t(sapply(wosdf, function(w) sapply(resptypes, function(r) sum(germprb$dummyWOSrelrank==w &germprb$culture.sim==e & germprb$numeric.phen.eff==r) )  )))
# germP.e <- t(sapply(gpxexrxw, function(e) colMeans(e/rowSums(e),na.rm=T) ))
# flwr.egrps <- matrix(rep(c("a","a","b"),each=1),nrow=3,ncol=1,byrow=T)
# germT.egrps <- matrix(c(rep(c("a","a","a"),each=1)),nrow=3,ncol=1,byrow=T)
# othT.egrps <- matrix(rep(c("a","a","b"),times=1),nrow=3,ncol=1,byrow=T)
# allT.egrps <- matrix(rep(c("a","a","b"),each=2),nrow=3,ncol=1,byrow=T)
# germP.egrps <- matrix(rep(c("b","a","a"),each=1),nrow=3,ncol=1,byrow=T)
# mnpredbins(summary(flower.cltr)$solutions[1,1],summary(flower.cltr)$solutions[-1,1],summary(flower.cltr)$cutpoints[1])
# #PREDICTIONS TERRIBLE. probably bc flwoering time with culture seems to be upside down parabola....which violates assumptions of model
# mnpredbins(summary(germtime.cltr)$solutions[1,1],summary(germtime.cltr)$solutions[-1,1],summary(germtime.cltr)$cutpoints[1])
# mnpredbins(summary(othertime.cltr)$solutions[1,1],summary(othertime.cltr)$solutions[-1,1],summary(othertime.cltr)$cutpoints[1])
# mnpredbins(summary(alltime.cltr2r)$solutions[1,1],summary(alltime.cltr2r)$solutions[-1,1],summary(alltime.cltr2r)$cutpoints[1])
# mnpredbins(summary(germprb.cltr)$solutions[1,1],summary(germprb.cltr)$solutions[-1,1],summary(germprb.cltr)$cutpoints[1])
# 
# 
# pdf("Phenology effects of experiment types.pdf",width=6,height=3) #general weirdness about pred vs data 
# par(mfrow=c(1,4))
# par(oma=c(0,5,0,1))
# par(mar=c(5,0,3,0))
# image(t(flwr.e), col=wb(100),zlim=c(0,1),xaxt="n",yaxt="n",ylab="",xlab="")
# 	axis(side = 1, at = c(0,0.5,1),labels=c("earlier","none","delayed"),las=2)
# 	axis(side = 2, at = c(0,0.5,1),labels=etypes,las=2)
# 	text(rep(c(0,0.5,1),each=3),rep( c(0,0.5,1),times=3), round(as.vector(flwr.e),digits=2))
# 	text(-0.2,rep( c(0,0.5,1)+.1, times=1), flwr.egrps,cex=0.75)
# # 	text(1.1,1.1,"*" )
# 	mtext(side=3,"Flowering time")
# image(t(germT.e), col=wb(100),zlim=c(0,1),xaxt="n",yaxt="n",ylab="",xlab="")
# 	axis(side = 1, at = c(0,0.5,1),labels=c("earlier","none","delayed"),las=2)
# 	text(rep(c(0,0.5,1),each=3),rep( c(0,0.5,1),times=3), round(as.vector(germT.e),digits=2))
# 	text(-0.2,rep( c(0,0.5,1)+.1, times=1), germT.egrps,cex=0.75)
# 	mtext(side=3,"Germination time")
# image(t(othT.e), col=wb(100),zlim=c(0,1),xaxt="n",yaxt="n",ylab="",xlab="")
# 	axis(side = 1, at = c(0,0.5,1),labels=c("earlier","none","delayed"),las=2)
# 	text(rep(c(0,0.5,1),each=3),rep( c(0,0.5,1),times=3), round(as.vector(othT.e),digits=2))
# 	text(-0.2,rep( c(0,0.5,1)+.1, times=1), othT.egrps,cex=0.75)
# 	mtext(side=3,"Other times")
# image(t(allT.e), col=wb(100),zlim=c(0,1),xaxt="n",yaxt="n",ylab="",xlab="")
# 	axis(side = 1, at = c(0,0.5,1),labels=c("earlier","none","delayed"),las=2)
# 	text(rep(c(0,0.5,1),each=3),rep( c(0,0.5,1),times=3), round(as.vector(allT.e),digits=2))
# 	mtext(side=3,"All times")
# 	text(-0.2,rep( c(0,0.5,1)+.1, times=1), allT.egrps,cex=0.75)
# dev.off()
# 
# pdf("Germination Prb effects of experiment types.pdf",width=3,height=3)#general weirdness about pred vs data / sig diffs
# par(oma=c(0,5,0,1))
# par(mar=c(5,0,3,0))
# image(t(germP.e), col=wb(100),zlim=c(0,1),xaxt="n",yaxt="n",ylab="",xlab="")
# 	axis(side = 1, at = c(0,0.5,1),labels=c("narrowed","none","expanded"),las=2)
# 	axis(side = 2, at = c(0,0.5,1),labels=etypes,las=2)
# 	text(rep(c(0,0.5,1),each=3),rep( c(0,0.5,1),times=3), round(as.vector(germP.e),digits=2))
# 	text(-0.2,rep( c(0,0.5,1)+.1, times=1), germP.egrps,cex=0.75)
# 	mtext(side=3,"Germination Prb")
# dev.off()




###ORDINAL PROBLEM AND no longer using these phenotype classes.
# othertime.tax <- MCMCglmm(numeric.phen.eff~micrtax.sim2,random = ~dummyWOSrelrank,data=othertime,family="ordinal",verbose=F,prior=priornuup,pl=T,pr=T,nitt=100000,thin=50,burnin=1000)
# # has issue described in handbook in estimating significance of mixed vs bact -- should probably be significant....
# othertime.tax15 <- MCMCglmm(numeric.phen.eff~micrtax.sim1.5,random = ~dummyWOSrelrank,data=othertime,family="ordinal",verbose=F,prior=priornuup,pl=T,pr=T,nitt=100000,thin=50,burnin=1000)
# #### only mixed bacteria differ
# othertime.cltr <- MCMCglmm(numeric.phen.eff~culture.sim,random = ~dummyWOSrelrank,data=othertime,family="ordinal",verbose=F,prior=priornuup,pl=T,pr=T,nitt=100000,thin=50,burnin=1000)
# othertime.scom <- MCMCglmm(numeric.phen.eff~Strainvcomm.1,random = ~dummyWOSrelrank,data=othertime,family="ordinal",verbose=F,prior=priornuup,pl=T,pr=T,nitt=100000,thin=50,burnin=1000)
# #.tax best by DIC; but rounds fluctuate somewhat for DIC
# othertime.mat <- MCMCglmm(numeric.phen.eff~mating2,random = ~dummyWOSrelrank,data=othertime,family="ordinal",verbose=F,prior=priornuup,pl=T,pr=T,nitt=100000,thin=50,burnin=1000)
# othertime.lh <- MCMCglmm(numeric.phen.eff~lifeform2,random = ~dummyWOSrelrank,data=othertime,family="ordinal",verbose=F,prior=priornuup,pl=T,pr=T,nitt=100000,thin=50,burnin=1000)
# 
# # alltime.chk1 <- MCMCglmm(numeric.phen.eff~1,data=alltime,family="ordinal",prior=priornr,pl=T,verbose=F,pr=T,nitt=100000,thin=10,burnin=1000)
# # alltime.chkw <- MCMCglmm(numeric.phen.eff~1,random = ~dummyWOSrelrank ,data=alltime,family="ordinal",prior=priornuup,pl=T,verbose=F,pr=T,nitt=100000,thin=10,burnin=1000)
# # alltime.chkp <- MCMCglmm(numeric.phen.eff~1,random = ~ phen.trait.3,data=alltime,family="ordinal",prior=priornuup,pl=T,verbose=F,pr=T,nitt=100000,thin=10,burnin=1000)
# # alltime.chkpw <- MCMCglmm(numeric.phen.eff~1,random = ~dummyWOSrelrank + phen.trait.3,data=alltime,family="ordinal",prior=priornuup2,pl=T,verbose=F,pr=T,nitt=100000,thin=10,burnin=1000)
# #adding phen trait improves fit (DIC) by only a little, keeping for now
# #using phen.trait.3 vs phen.trait.2 only improves things a very little, that nu parameter seems to be doing its job
# alltime.tax2r <- MCMCglmm(numeric.phen.eff~micrtax.sim2,random = ~dummyWOSrelrank + phen.trait.3,data=alltime,family="ordinal",prior=priornuup2,pl=T,verbose=F,pr=T,nitt=100000,thin=10,burnin=1000)
# alltime.tax2r15 <- MCMCglmm(numeric.phen.eff~micrtax.sim1.5,random = ~dummyWOSrelrank + phen.trait.3,data=alltime,family="ordinal",prior=priornuup2,pl=T,verbose=F,pr=T,nitt=10000,thin=10,burnin=1000)
# ##
# alltime.cltr2r <- MCMCglmm(numeric.phen.eff~culture.sim,random = ~dummyWOSrelrank + phen.trait.3,data=alltime,family="ordinal",prior=priornuup2,pl=T,verbose=F,pr=T,nitt=100000,thin=10,burnin=1000)
# alltime.scom2r <- MCMCglmm(numeric.phen.eff~Strainvcomm.1,random = ~dummyWOSrelrank + phen.trait.3,data=alltime,family="ordinal",prior=priornuup2,pl=T,verbose=F,pr=T,nitt=100000,thin=10,burnin=1000)
# #eff sampl issue in wolrelrank rand effect (.tax and .cltr have mild issue with this)
# #.tax best by DIC; no real fluctuations
# alltime.lh <- MCMCglmm(numeric.phen.eff~lifeform2,random = ~dummyWOSrelrank + phen.trait.3,data=alltime,family="ordinal",prior=priornuup2,pl=T,verbose=F,pr=T,nitt=100000,thin=10,burnin=1000)
# alltime.mat <- MCMCglmm(numeric.phen.eff~mating2,random = ~dummyWOSrelrank + phen.trait.3,data=alltime,family="ordinal",prior=priornuup2,pl=T,verbose=F,pr=T,nitt=100000,thin=10,burnin=1000)
# 
# rprd.tax <- MCMCglmm(numeric.phen.eff~micrtax.sim2,random = ~dummyWOSrelrank,data=reprotime,family="ordinal",verbose=F,prior=priornuup,pl=T,pr=T,nitt=100000,thin=10,burnin=1000)
# rprd.tax15 <- MCMCglmm(numeric.phen.eff~micrtax.sim1.5,random = ~dummyWOSrelrank,data=reprotime,family="ordinal",verbose=F,prior=priornuup,pl=T,pr=T,nitt=100000,thin=10,burnin=1000)
# #fits worse; only mixed bacteria differs
# rprd.cltr <- MCMCglmm(numeric.phen.eff~culture.sim,random = ~dummyWOSrelrank,data=reprotime,family="ordinal",verbose=F,prior=priornuup,pl=T,pr=T,nitt=100000,thin=10,burnin=1000)
# rprd.scom <- MCMCglmm(numeric.phen.eff~Strainvcomm.1,random = ~dummyWOSrelrank,data=reprotime,family="ordinal",verbose=F,prior=priornuup,pl=T,pr=T,nitt=100000,thin=10,burnin=1000)
# 
# rprd.lh <- MCMCglmm(numeric.phen.eff~mating2,random = ~dummyWOSrelrank,data=reprotime,family="ordinal",verbose=F,prior=priornuup,pl=T,pr=T,nitt=100000,thin=10,burnin=1000)
# rprd.mat <- MCMCglmm(numeric.phen.eff~lifeform2,random = ~dummyWOSrelrank,data=reprotime,family="ordinal",verbose=F,prior=priornuup,pl=T,pr=T,nitt=100000,thin=10,burnin=1000)
# veg.mat <- MCMCglmm(numeric.phen.eff~mating2,random = ~dummyWOSrelrank,data=vegtime,family="ordinal",verbose=F,prior=priornuup,pl=T,pr=T,nitt=100000,thin=10,burnin=1000)
# veg.lh <- MCMCglmm(numeric.phen.eff~lifeform2,random = ~dummyWOSrelrank,data=vegtime,family="ordinal",verbose=F,prior=priornuup,pl=T,pr=T,nitt=100000,thin=10,burnin=1000)
# 
# veg.tax <- MCMCglmm(numeric.phen.eff~micrtax.sim2,random = ~dummyWOSrelrank,data=vegtime,family="ordinal",verbose=F,prior=priornuup,pl=T,pr=T,nitt=100000,thin=10,burnin=1000)
# veg.cltr <- MCMCglmm(numeric.phen.eff~culture.sim,random = ~dummyWOSrelrank,data=vegtime,family="ordinal",verbose=F,prior=priornuup,pl=T,pr=T,nitt=100000,thin=10,burnin=1000)
# veg.scom <- MCMCglmm(numeric.phen.eff~Strainvcomm.1,random = ~dummyWOSrelrank,data=vegtime,family="ordinal",verbose=F,prior=priornuup,pl=T,pr=T,nitt=100000,thin=10,burnin=1000)
# veg.tax15 <- MCMCglmm(numeric.phen.eff~micrtax.sim1.5,random = ~dummyWOSrelrank,data=vegtime,family="ordinal",verbose=F,prior=priornuup,pl=T,pr=T,nitt=100000,thin=10,burnin=1000)
# 



##FROM JARROD HADFIELD ON RSIG ; 2010
#https://stat.ethz.ch/pipermail/r-sig-mixed-models/2010q2/003671.html
###how to find probabilities of falling into different categoires. l is the latent variable computed from the parameters as estimated in the model for fixed slopes (beta) and random slopes/ints (u)
# 
# l = Xb+Zu+e
# 
# The probabilities of falling into each of the four categories are:
# 
# 
# pnorm(-l)
# 
# pnorm(cp[1]-l)-pnorm(-l)
# 
#   pnorm(cp[2]-l)-pnorm(cp[1]-l)
# 
# 1-pnorm(cp[2]-l)

##so I infer that when k is 3 the probabilities of each k are:
# pnorm(-l)
# 
# pnorm(cp[1]-l)-pnorm(-l)
# # 
# 1-pnorm(cp[1]-l)







