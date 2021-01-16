library(UpSetR) #for set plots
library(MCMCglmm) #for models
library(grDevices) #for color

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

flowering <- split_recs_scr[split_recs_scr$phen.trait.2=="flowering time",]
flower.tax <- MCMCglmm(numeric.phen.eff~micrtax.sim2,random = ~WOS.relevance.rank,data=flowering,family="ordinal",verbose=F,pr=T,nitt=1000000,thin=10,burnin=1000)
flower.cltr <- MCMCglmm(numeric.phen.eff~culture.sim,random = ~WOS.relevance.rank,data=flowering,family="ordinal",verbose=F,pr=T,nitt=1000000,thin=10,burnin=1000)
flower.scom <- MCMCglmm(numeric.phen.eff~Strainvcomm.1,random = ~WOS.relevance.rank,data=flowering,family="ordinal",verbose=F,pr=T,nitt=1000000,thin=10,burnin=1000)

germprb <- split_recs_scr[split_recs_scr$phen.trait.2=="germination prb",]

germprb.tax <- MCMCglmm(numeric.phen.eff~micrtax.sim2,random = ~WOS.relevance.rank,data=germprb,family="ordinal",verbose=F,pr=T,nitt=1000000,thin=10,burnin=1000)
germprb.cltr <- MCMCglmm(numeric.phen.eff~culture.sim,random = ~WOS.relevance.rank,data=germprb,family="ordinal",verbose=F,pr=T,nitt=1000000,thin=10,burnin=1000)
germprb.scom <- MCMCglmm(numeric.phen.eff~Strainvcomm.1,random = ~WOS.relevance.rank,data=germprb,family="ordinal",verbose=F,pr=T,nitt=1000000,thin=10,burnin=1000)


germtime <- split_recs_scr[split_recs_scr$phen.trait.2=="germination time",]

germtime.tax <- MCMCglmm(numeric.phen.eff~micrtax.sim2,random = ~WOS.relevance.rank,data=germtime,family="ordinal",verbose=F,pr=T,nitt=1000000,thin=10,burnin=1000)
germtime.cltr <- MCMCglmm(numeric.phen.eff~culture.sim,random = ~WOS.relevance.rank,data=germtime,family="ordinal",verbose=F,pr=T,nitt=1000000,thin=10,burnin=1000)
germtime.scom <- MCMCglmm(numeric.phen.eff~Strainvcomm.1,random = ~WOS.relevance.rank,data=germtime,family="ordinal",verbose=F,pr=T,nitt=1000000,thin=10,burnin=1000)

othertime <- split_recs_scr[split_recs_scr$phen.trait.2%in% 
				c("budburst", "floral budset time", "flower senescence time", "fruiting time","maturation time" ,  "peak flowering" , "phyllochron","senescence time") ,]
othertime.tax <- MCMCglmm(numeric.phen.eff~micrtax.sim2,random = ~WOS.relevance.rank,data=othertime,family="ordinal",verbose=F,pr=T,nitt=1000000,thin=10,burnin=1000)
othertime.cltr <- MCMCglmm(numeric.phen.eff~culture.sim,random = ~WOS.relevance.rank,data=othertime,family="ordinal",verbose=F,pr=T,nitt=1000000,thin=10,burnin=1000)
othertime.scom <- MCMCglmm(numeric.phen.eff~Strainvcomm.1,random = ~WOS.relevance.rank,data=othertime,family="ordinal",verbose=F,pr=T,nitt=1000000,thin=10,burnin=1000)

alltime <- split_recs_scr[split_recs_scr$phen.trait.2%in% 
				c("germination time","budburst", "floral budset time", "flowering time", "flower senescence time", "fruiting time","maturation time" ,  "peak flowering" , "phyllochron","senescence time") ,]
# alltime.tax <- MCMCglmm(numeric.phen.eff~micrtax.sim2,random = ~WOS.relevance.rank + phen.trait.2,data=alltime,family="ordinal",verbose=F,pr=T,nitt=1000000,thin=10,burnin=1000)
alltime.tax2r <- MCMCglmm(numeric.phen.eff~micrtax.sim2,random = ~WOS.relevance.rank + phen.trait.2,data=alltime,family="ordinal",verbose=F,pr=T,nitt=1000000,thin=10,burnin=1000)
alltime.cltr2r <- MCMCglmm(numeric.phen.eff~culture.sim,random = ~WOS.relevance.rank + phen.trait.2,data=alltime,family="ordinal",verbose=F,pr=T,nitt=1000000,thin=10,burnin=1000)
alltime.scom2r <- MCMCglmm(numeric.phen.eff~Strainvcomm.1,random = ~WOS.relevance.rank + phen.trait.2,data=alltime,family="ordinal",verbose=F,pr=T,nitt=1000000,thin=10,burnin=1000)



wb <- colorRampPalette(c(rgb(1,1,1),rgb(0,0,1)))

mtypes  <- c("bacteria","MF","mixed","otherfungi")
resptypes <- c(-1,0,1)
flwr.m1 <- matrix(sapply(mtypes, function(m) sapply(resptypes, function(r) 
		sum(flowering$micrtax.sim2==m & flowering$numeric.phen.eff==r ))),
		ncol=3,byrow=T)
flwr.m <- flwr.m1/rowSums(flwr.m1)
germT.m1 <- matrix(sapply(mtypes, function(m) sapply(resptypes, function(r) 
		sum(germtime$micrtax.sim2==m & germtime$numeric.phen.eff==r ))),
		ncol=3,byrow=T)
germT.m <- germT.m1/rowSums(germT.m1)
othT.m1 <- matrix(sapply(mtypes, function(m) sapply(resptypes, function(r) 
		sum(othertime$micrtax.sim2==m & othertime$numeric.phen.eff==r ))),
		ncol=3,byrow=T)
othT.m <- othT.m1/rowSums(othT.m1)
allT.m1 <- matrix(sapply(mtypes, function(m) sapply(resptypes, function(r) 
		sum(alltime$micrtax.sim2==m & alltime$numeric.phen.eff==r ))),
		ncol=3,byrow=T)
allT.m <- allT.m1/rowSums(allT.m1)

germP.m1 <- matrix(sapply(mtypes, function(m) sapply(resptypes, function(r) 
		sum(germprb$micrtax.sim2==m & germprb$numeric.phen.eff==r ))),
		ncol=3,byrow=T)
germP.m <- germP.m1/rowSums(germP.m1)



pdf("Phenology effects of microbe types.pdf",width=6,height=3)
par(mfrow=c(1,4))
par(oma=c(0,5,0,1))
par(mar=c(5,0,3,0))
image(t(flwr.m), col=wb(100),zlim=c(0,1),xaxt="n",yaxt="n",ylab="",xlab="")
	axis(side = 1, at = c(0,0.5,1),labels=c("earlier","none","delayed"),las=2)
	axis(side = 2, at = c(0,0.33,0.66,1),labels=mtypes,las=2)
	text(rep(c(0,0.5,1),each=4),rep( c(0,0.33,0.66,1),times=3), as.vector(flwr.m1))
	text(c(1.1,1.1),c(0.76,1.1),c("*",".") )
	mtext(side=3,"Flowering time")
image(t(germT.m), col=wb(100),zlim=c(0,1),xaxt="n",yaxt="n",ylab="",xlab="")
	axis(side = 1, at = c(0,0.5,1),labels=c("earlier","none","delayed"),las=2)
	text(rep(c(0,0.5,1),each=4),rep( c(0,0.33,0.66,1),times=3), as.vector(germT.m1))
	text(c(0.1,1.1),c(0.43,0.76),c(".","*") )
	mtext(side=3,"Germination time")
image(t(othT.m), col=wb(100),zlim=c(0,1),xaxt="n",yaxt="n",ylab="",xlab="")
	axis(side = 1, at = c(0,0.5,1),labels=c("earlier","none","delayed"),las=2)
	text(rep(c(0,0.5,1),each=4),rep( c(0,0.33,0.66,1),times=3), as.vector(othT.m1))
	text(0.1,1.1,"*" )
	mtext(side=3,"Other times")
image(t(allT.m), col=wb(100),zlim=c(0,1),xaxt="n",yaxt="n",ylab="",xlab="")
	axis(side = 1, at = c(0,0.5,1),labels=c("earlier","none","delayed"),las=2)
	text(rep(c(0,0.5,1),each=4),rep( c(0,0.33,0.66,1),times=3), as.vector(allT.m1))
	mtext(side=3,"All times")
	text(c(0.1,1.1,1.1),c(0.1,0.66,1.1),c(".","***","*") )
dev.off()

pdf("Germination Prb effects of microbe types.pdf",width=3,height=3)
par(oma=c(0,5,0,1))
par(mar=c(5,0,3,0))
image(t(germP.m), col=wb(100),zlim=c(0,1),xaxt="n",yaxt="n",ylab="",xlab="")
	axis(side = 1, at = c(0,0.5,1),labels=c("narrowed","none","expanded"),las=2)
	axis(side = 2, at = c(0,0.33,0.66,1),labels=mtypes,las=2)
	text(rep(c(0,0.5,1),each=4),rep( c(0,0.33,0.66,1),times=3), as.vector(germP.m1))
	text(c(1.1,1.1,0.1,0.1),c(0.1,0.43,0.76,1.1),c("***","*","**","***") )
	mtext(side=3,"Germination Prb")
dev.off()



etypes  <- c("culture","direct","removed") #leaves out observed, but obs is only 1 datapoint; might consider whole record?
resptypes <- c(-1,0,1)
flwr.e1 <- matrix(sapply(etypes, function(e) sapply(resptypes, function(r) 
		sum(flowering$culture.sim==e & flowering$numeric.phen.eff==r ))),
		ncol=3,byrow=T)
flwr.e <- flwr.e1/rowSums(flwr.e1)
germT.e1 <- matrix(sapply(etypes, function(e) sapply(resptypes, function(r) 
		sum(germtime$culture.sim==e & germtime$numeric.phen.eff==r ))),
		ncol=3,byrow=T)
germT.e <- germT.e1/rowSums(germT.e1)
othT.e1 <- matrix(sapply(etypes, function(e) sapply(resptypes, function(r) 
		sum(othertime$culture.sim==e & othertime$numeric.phen.eff==r ))),
		ncol=3,byrow=T)
othT.e <- othT.e1/rowSums(othT.e1)
allT.e1 <- matrix(sapply(etypes, function(e) sapply(resptypes, function(r) 
		sum(alltime$culture.sim==e & alltime$numeric.phen.eff==r ))),
		ncol=3,byrow=T)
allT.e <- allT.e1/rowSums(allT.e1)

germP.e1 <- matrix(sapply(etypes, function(e) sapply(resptypes, function(r) 
		sum(germprb$culture.sim==e & germprb$numeric.phen.eff==r ))),
		ncol=3,byrow=T)
germP.e <- germP.e1/rowSums(germP.e1)

pdf("Phenology effects of experiment types.pdf",width=6,height=3)
par(mfrow=c(1,4))
par(oma=c(0,5,0,1))
par(mar=c(5,0,3,0))
image(t(flwr.e), col=wb(100),zlim=c(0,1),xaxt="n",yaxt="n",ylab="",xlab="")
	axis(side = 1, at = c(0,0.5,1),labels=c("earlier","none","delayed"),las=2)
	axis(side = 2, at = c(0,0.5,1),labels=etypes,las=2)
	text(rep(c(0,0.5,1),each=3),rep( c(0,0.5,1),times=3), as.vector(flwr.e1))
	text(1.1,1.1,"*" )
	mtext(side=3,"Flowering time")
image(t(germT.e), col=wb(100),zlim=c(0,1),xaxt="n",yaxt="n",ylab="",xlab="")
	axis(side = 1, at = c(0,0.5,1),labels=c("earlier","none","delayed"),las=2)
	text(rep(c(0,0.5,1),each=3),rep( c(0,0.5,1),times=3), as.vector(germT.e1))
	text(c(1.1,1.1),c(0.6,1.1),c("*","*") )
	mtext(side=3,"Germination time")
image(t(othT.e), col=wb(100),zlim=c(0,1),xaxt="n",yaxt="n",ylab="",xlab="")
	axis(side = 1, at = c(0,0.5,1),labels=c("earlier","none","delayed"),las=2)
	text(rep(c(0,0.5,1),each=3),rep( c(0,0.5,1),times=3), as.vector(othT.e1))
	text(1.1,1.1,"*" )
	mtext(side=3,"Other times")
image(t(allT.e), col=wb(100),zlim=c(0,1),xaxt="n",yaxt="n",ylab="",xlab="")
	axis(side = 1, at = c(0,0.5,1),labels=c("earlier","none","delayed"),las=2)
	text(rep(c(0,0.5,1),each=3),rep( c(0,0.5,1),times=3), as.vector(allT.e1))
	mtext(side=3,"All times")
 	text(c(1.1,1.1),c(0.6,1.1),c(".","*") )
dev.off()

pdf("Germination Prb effects of experiment types.pdf",width=3,height=3)
par(oma=c(0,5,0,1))
par(mar=c(5,0,3,0))
image(t(germP.e), col=wb(100),zlim=c(0,1),xaxt="n",yaxt="n",ylab="",xlab="")
	axis(side = 1, at = c(0,0.5,1),labels=c("narrowed","none","expanded"),las=2)
	axis(side = 2, at = c(0,0.5,1),labels=etypes,las=2)
	text(rep(c(0,0.5,1),each=3),rep( c(0,0.5,1),times=3), as.vector(germP.e1))
	text(c(1.1,0.1),c(0.1,0.6),c("***","*") )
	mtext(side=3,"Germination Prb")
dev.off()




itypes  <- c("community","strain","strain mix") #leaves out observed, but obs is only 1 datapoint; might consider whole record?
resptypes <- c(-1,0,1)
flwr.i1 <- matrix(sapply(itypes, function(i) sapply(resptypes, function(r) 
		sum(flowering$Strainvcomm.1==i & flowering$numeric.phen.eff==r ))),
		ncol=3,byrow=T)
flwr.i <- flwr.i1/rowSums(flwr.i1)
germT.i1 <- matrix(sapply(itypes, function(i) sapply(resptypes, function(r) 
		sum(germtime$Strainvcomm.1==i & germtime$numeric.phen.eff==r ))),
		ncol=3,byrow=T)
germT.i <- germT.i1/rowSums(germT.i1)
othT.i1 <- matrix(sapply(itypes, function(i) sapply(resptypes, function(r) 
		sum(othertime$Strainvcomm.1==i & othertime$numeric.phen.eff==r ))),
		ncol=3,byrow=T)
othT.i <- othT.i1/rowSums(othT.i1)
allT.i1 <- matrix(sapply(itypes, function(i) sapply(resptypes, function(r) 
		sum(alltime$Strainvcomm.1==i & alltime$numeric.phen.eff==r ))),
		ncol=3,byrow=T)
allT.i <- allT.i1/rowSums(allT.i1)

germP.i1 <- matrix(sapply(itypes, function(i) sapply(resptypes, function(r) 
		sum(germprb$Strainvcomm.1==i & germprb$numeric.phen.eff==r ))),
		ncol=3,byrow=T)
germP.i <- germP.i1/rowSums(germP.i1)


pdf("Phenology effects of inoculation types.pdf",width=6,height=3)
par(mfrow=c(1,4))
par(oma=c(0,5,0,1))
par(mar=c(5,0,3,0))
image(t(flwr.i), col=wb(100),zlim=c(0,1),xaxt="n",yaxt="n",ylab="",xlab="")
	axis(side = 1, at = c(0,0.5,1),labels=c("earlier","none","delayed"),las=2)
	axis(side = 2, at = c(0,0.5,1),labels=itypes,las=2)
	text(rep(c(0,0.5,1),each=3),rep( c(0,0.5,1),times=3), as.vector(flwr.i1))
	text(c(1.1,0.1,0.1),c(0.1,0.6,1.1),c(".","*","*") )
	mtext(side=3,"Flowering time")
image(t(germT.i), col=wb(100),zlim=c(0,1),xaxt="n",yaxt="n",ylab="",xlab="")
	axis(side = 1, at = c(0,0.5,1),labels=c("earlier","none","delayed"),las=2)
	text(rep(c(0,0.5,1),each=3),rep( c(0,0.5,1),times=3), as.vector(germT.i1))
	text(c(0.1),c(0.6),c("**") )
	mtext(side=3,"Germination time")
image(t(othT.i), col=wb(100),zlim=c(0,1),xaxt="n",yaxt="n",ylab="",xlab="")
	axis(side = 1, at = c(0,0.5,1),labels=c("earlier","none","delayed"),las=2)
	text(rep(c(0,0.5,1),each=3),rep( c(0,0.5,1),times=3), as.vector(othT.i1))
	text(c(0.1),c(0.6),c(".") )
	mtext(side=3,"Other times")
image(t(allT.i), col=wb(100),zlim=c(0,1),xaxt="n",yaxt="n",ylab="",xlab="")
	axis(side = 1, at = c(0,0.5,1),labels=c("earlier","none","delayed"),las=2)
	text(rep(c(0,0.5,1),each=3),rep( c(0,0.5,1),times=3), as.vector(allT.i1))
	mtext(side=3,"All times")
	text(c(0.1,0.1),c(0.6,1.1),c("***","*") )
dev.off()
#these model results are a bit odd
pdf("Germination Prb effects of inoculation types.pdf",width=3,height=3)
par(oma=c(0,5,0,1))
par(mar=c(5,0,3,0))
image(t(germP.i), col=wb(100),zlim=c(0,1),xaxt="n",yaxt="n",ylab="",xlab="")
	axis(side = 1, at = c(0,0.5,1),labels=c("narrowed","none","expanded"),las=2)
	axis(side = 2, at = c(0,0.5,1),labels=itypes,las=2)
	text(rep(c(0,0.5,1),each=3),rep( c(0,0.5,1),times=3), as.vector(germP.i1))
	text(c(1.1,1.1),c(0.6,1.1),c("**","*") )
	mtext(side=3,"Germination Prb")
dev.off()
