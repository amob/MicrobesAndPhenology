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
split_recs <- read.csv("trimmed RECORDS FOR SPLITTING v4.csv", #sep="\t",
					stringsAsFactors=F,na.strings = c("NA","")) 
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


##these questions are all covered later
# modmtypeBact <- MCMCglmm(cbind(effobs,fails)~bacteria,family ="multinomial2" ,data=effcatdat,pr=T,verbose=F,nitt=100000,burnin=10000,thin=100)#n.s. more sig effs than avg rate without
# modmtypeFungi <- MCMCglmm(cbind(effobs,fails)~fungi,family ="multinomial2" ,data=effcatdat,pr=T,verbose=F,nitt=100000,burnin=10000,thin=100)#n.s. more sig effects than avg rate without
# modmtypeMF <- MCMCglmm(cbind(effobs,fails)~MF,family ="multinomial2" ,data=effcatdat,pr=T,verbose=F,nitt=100000,burnin=10000,thin=100) #n.s. fewer significant effects than average rate without
# modmtypeMixed <- MCMCglmm(cbind(effobs,fails)~mixed,family ="multinomial2" ,data=effcatdat,pr=T,verbose=F,nitt=100000,burnin=10000,thin=100)#n.s. more sig effs than avg rate without
# modmtypeVirus <- MCMCglmm(cbind(effobs,fails)~virus,family ="multinomial2" ,data=effcatdat,pr=T,verbose=F,nitt=100000,burnin=10000,thin=100) #too few obs. not to be trusted
# 
# modlocroot <- MCMCglmm(cbind(effobs,fails)~root,family ="multinomial2" ,data=effcatdat,pr=T,verbose=F,nitt=100000,burnin=10000,thin=100)#n.s. more sig effs than avg rate without
# modlocshoot <- MCMCglmm(cbind(effobs,fails)~shoot,family ="multinomial2" ,data=effcatdat,pr=T,verbose=F,nitt=100000,burnin=10000,thin=100)#n.s. fewer sig effs than avg rate without
# modlocreproductive <- MCMCglmm(cbind(effobs,fails)~reproductive,family ="multinomial2" ,data=effcatdat,pr=T,verbose=F,nitt=100000,burnin=10000,thin=100)#n.s. fewer sig effs than avg rate without
# modlocseed <- MCMCglmm(cbind(effobs,fails)~seed,family ="multinomial2" ,data=effcatdat,pr=T,verbose=F,nitt=100000,burnin=10000,thin=100)#n.s. more sig effs than avg rate without
# 
# ##model fit issues....effective samples is low in most
# modtraitflwrT <- MCMCglmm(cbind(effobs,fails)~flowering.time,family ="multinomial2" ,data=effcatdat,pr=T,verbose=F,nitt=500000,burnin=10000,thin=100)#n.s. fewer sig effs than avg rate without
# modtraitfbT <- MCMCglmm(cbind(effobs,fails)~floral.budset.time,family ="multinomial2" ,data=effcatdat,pr=T,verbose=F,nitt=500000,burnin=10000,thin=100)#n.s. more sig effs than avg rate without
# modtraitgT <- MCMCglmm(cbind(effobs,fails)~germination.time,family ="multinomial2" ,data=effcatdat,pr=T,verbose=F,nitt=500000,burnin=10000,thin=100)#n.s. more sig effs than avg rate without
# modtraitgP <- MCMCglmm(cbind(effobs,fails)~germination.prb,family ="multinomial2" ,data=effcatdat,pr=T,verbose=F,nitt=500000,burnin=10000,thin=100)#n.s. more sig effs than avg rate without
# modtraitphyl <- MCMCglmm(cbind(effobs,fails)~phyllochron,family ="multinomial2" ,data=effcatdat,pr=T,verbose=F,nitt=500000,burnin=10000,thin=100)#"sig" more, but there are only 2 studies
# modtraitsncT <- MCMCglmm(cbind(effobs,fails)~senescence.time,family ="multinomial2" ,data=effcatdat,pr=T,verbose=F,nitt=500000,burnin=10000,thin=100)#n.s. fewer sig effs than avg rate without
# modtraitmatT <- MCMCglmm(cbind(effobs,fails)~maturation.time,family ="multinomial2" ,data=effcatdat,pr=T,verbose=F,nitt=500000,burnin=10000,thin=100)#n.s. more sig effs than avg rate without
# modtraitfrtT <- MCMCglmm(cbind(effobs,fails)~fruiting.time,family ="multinomial2" ,data=effcatdat,pr=T,verbose=F,nitt=500000,burnin=10000,thin=100)#n.s. slightly fewer sig effs than avg rate without
# modtraitflwrD <- MCMCglmm(cbind(effobs,fails)~flowering.duration,family ="multinomial2" ,data=effcatdat,pr=T,verbose=F,nitt=500000,burnin=10000,thin=100)#n.s. more sig effs than avg rate without
# modtraitrdT <- MCMCglmm(cbind(effobs,fails)~root.development.timing,family ="multinomial2" ,data=effcatdat,pr=T,verbose=F,nitt=500000,burnin=10000,thin=100)#n.s. more sig effs than avg rate without
# modtraitbudB <- MCMCglmm(cbind(effobs,fails)~budburst,family ="multinomial2" ,data=effcatdat,pr=T,verbose=F,nitt=500000,burnin=10000,thin=100)#n.s. more sig effs than avg rate without; close to sig, but there's only three studies here




################################################################################
####Analysis of phenological timing records
################################################################################

##consistency edits

split_recs_scr$MicrobeMPH [split_recs_scr$MicrobeMPH %in% c("other benefical","other benefial")] <- "other beneficial"
split_recs_scr$MicrobeMPH [split_recs_scr$MicrobeMPH=="pethogen"] <- "pathogen"
split_recs_scr$MicrobeMPH [split_recs_scr$MicrobeMPH=="nutrient"] <- "nutrients"
split_recs_scr$MicrobeMPH [split_recs_scr$MicrobeMPH=="hormones"] <- "phytohormones"

split_recs_scr$direction.effect.split[split_recs_scr$direction.effect.split== "later"] <- "delayed" 
split_recs_scr$direction.effect.split[split_recs_scr$direction.effect.split %in% c("decreased","reduced")] <- "narrowed" 
split_recs_scr$direction.effect.split[split_recs_scr$direction.effect.split %in% c("increased","extended")] <- "expanded" 
split_recs_scr$direction.effect.split[split_recs_scr$direction.effect.split== "none "] <- "none" 

split_recs_scr$culture.sim[split_recs_scr$culture.sim %in% c("culture; direct","cutlure")] <- "culture" # making a decision on Lu2018 
split_recs_scr$micrtax.sim2[split_recs_scr$micrtax.sim2=="fungi"] <- "otherfungi" #changing alphabetical position for convenience
split_recs_scr$micrtax.sim2[split_recs_scr$micrtax.sim2=="rhizobia"] <- "bacteria" #rhizobia not allowed in this category
length(unique(split_recs_scr$WOS.relevance.rank[split_recs_scr$micrtax.sim2=="MF"] ))
length(unique(split_recs_scr$WOS.relevance.rank[split_recs_scr$micrtax.sim2=="bacteria"] ))
length(unique(split_recs_scr$WOS.relevance.rank[split_recs_scr$micrtax.sim2=="otherfungi"] ))
length(unique(split_recs_scr$WOS.relevance.rank[split_recs_scr$micrtax.sim2=="mixed"] ))


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
length(unique(split_recs_scr$WOS.relevance.rank[split_recs_scr$lifeform2=="annual" | split_recs_scr$lifeform2=="both"]))
length(unique(split_recs_scr$WOS.relevance.rank[split_recs_scr$lifeform2=="perennial"]))
split_recs_scr$mating2 <- split_recs_scr$mating
split_recs_scr$mating2 [is.na(split_recs_scr$mating)] <- "unk"
split_recs_scr$mating2 [split_recs_scr$mating=="outcross?"] <- "unk"
split_recs_scr$mating2 [split_recs_scr$mating=="both"] <- "s+o"
split_recs_scr$mating2 [split_recs_scr$mating %in% c("selfie","selifer","selfier")] <- "selfer" #correcting spelling errors
split_recs_scr$phen.trait.2[split_recs_scr$phen.trait.2=='"budburst"'] <- "budburst" #correcting spelling errors



split_recs_scr$micrloc3 <- split_recs_scr$Micrloc.2
# split_recs_scr$micrloc3 [split_recs_scr$Micrloc.2%in% c("reproductive; shoot","root; shoot; reproductive","root; shoot; reproductive; seed")] <- "various incl reproductive"
# split_recs_scr$micrloc3 [split_recs_scr$Micrloc.2%in% c("seed; shoot","root; shoot","root; seed","root; shoot; seed")] <- "various not incl reproductive"
split_recs_scr$micrloc3 [split_recs_scr$Micrloc.2%in% c("reproductive; shoot","root; shoot; reproductive","root; shoot; reproductive; seed","seed; shoot","root; shoot","root; seed","root; shoot; seed")] <- "multiple"
split_recs_scr$micrloc3 [split_recs_scr$Micrloc.2=="seed "] <- "seed"
split_recs_scr$micrloc3 [split_recs_scr$Micrloc.2=="leaves"] <- "shoot"


##variables for binomial versions
split_recs_scr$isEarlyOrNrw <- ifelse(numeric.phen.eff=="-1",1,0)
split_recs_scr$notEarlyOrNrw <- 1-split_recs_scr$isEarlyOrNrw
split_recs_scr$isDelayOrExp <- ifelse(numeric.phen.eff=="1",1,0)
split_recs_scr$notDelayOrExp <- 1-split_recs_scr$isDelayOrExp
split_recs_scr$isSig <- ifelse(numeric.phen.eff=="1" | numeric.phen.eff=="-1",1,0)
split_recs_scr$notSig <- 1- split_recs_scr$isSig 

#remove singleton taxa type records from all for now
split_recs_scr <- split_recs_scr[-which(split_recs_scr$micrtax.sim2 %in% c("virus","MF; fungi")),]


##numbers for text
split_recs_scr$anyeffnum <- as.numeric(as.factor(split_recs_scr$anyeffect))-1
sum(   tapply(split_recs_scr$anyeffnum[split_recs_scr$micrtax.sim2=="MF"],split_recs_scr$WOS.relevance.rank[split_recs_scr$micrtax.sim2=="MF"],sum) ==0)#4
length(tapply(split_recs_scr$anyeffnum[split_recs_scr$micrtax.sim2=="MF"],split_recs_scr$WOS.relevance.rank[split_recs_scr$micrtax.sim2=="MF"],sum))#17  17-4 = 13
sum(   tapply(split_recs_scr$anyeffnum[split_recs_scr$micrtax.sim2=="bacteria"],split_recs_scr$WOS.relevance.rank[split_recs_scr$micrtax.sim2=="bacteria"],sum) ==0)#0
length(tapply(split_recs_scr$anyeffnum[split_recs_scr$micrtax.sim2=="bacteria"],split_recs_scr$WOS.relevance.rank[split_recs_scr$micrtax.sim2=="bacteria"],sum))#13
table(split_recs_scr$WOS.relevance.rank[split_recs_scr$micrloc3=="various incl reproductive"])
length(unique(split_recs_scr$WOS.relevance.rank[split_recs_scr$phen.trait.2=="fruiting time"]))
length(unique(split_recs_scr$WOS.relevance.rank[split_recs_scr$phen.trait.2=="senescence time"]))
length(unique(split_recs_scr$WOS.relevance.rank[split_recs_scr$phen.trait.2 %in% c("phyllochron","budburst") ]))
table(split_recs_scr$MicrobeMPH)


#rotate var levels
#
split_recs_scr$micrtaxmfa <- split_recs_scr$micrtax.sim2
split_recs_scr$micrtaxmfa[split_recs_scr$micrtaxmfa=="MF"] <- "0MF"
split_recs_scr$micrtaxmixa <- split_recs_scr$micrtax.sim2
split_recs_scr$micrtaxmixa[split_recs_scr$micrtaxmixa=="mixed"] <- "0mix"
#
split_recs_scr$strainvcommsa <- split_recs_scr$Strainvcomm.1
split_recs_scr$strainvcommsa[split_recs_scr$strainvcommsa=="strain"] <- "0strain"
#
split_recs_scr$lifeform2pa <- split_recs_scr$lifeform2
split_recs_scr$lifeform2pa[split_recs_scr$lifeform2pa == "perennial"] <- "0perennial"
#
split_recs_scr$mphpa <- split_recs_scr$MicrobeMPH
split_recs_scr$mphpa[split_recs_scr$mphpa == "pathogen"] <- "0path"
split_recs_scr$mphba <- split_recs_scr$MicrobeMPH
split_recs_scr$mphba[split_recs_scr$mphba == "other beneficial"] <- "0ob"
#
split_recs_scr$micrloc3seeda <- split_recs_scr$micrloc3
split_recs_scr$micrloc3seeda[split_recs_scr$micrloc3seeda=="seed"] <- "0seed"
split_recs_scr$micrloc3shoota <- split_recs_scr$micrloc3
split_recs_scr$micrloc3shoota[split_recs_scr$micrloc3shoota=="shoot"] <- "0shoot"
# split_recs_scr$micrloc3vira <- split_recs_scr$micrloc3
# split_recs_scr$micrloc3vira[split_recs_scr$micrloc3vira=="multiple"] <- "0vir"

split_recs_scr$phen.trait.2fba <- split_recs_scr$phen.trait.2
split_recs_scr$phen.trait.2fba[split_recs_scr$phen.trait.2=="floral budset time"] <- "0fbt"
split_recs_scr$phen.trait.2fla <- split_recs_scr$phen.trait.2
split_recs_scr$phen.trait.2fla[split_recs_scr$phen.trait.2=="flowering time"] <- "0flwr"
split_recs_scr$phen.trait.2fra <- split_recs_scr$phen.trait.2
split_recs_scr$phen.trait.2fra[split_recs_scr$phen.trait.2=="fruiting time"] <- "0frt"
split_recs_scr$phen.trait.2ga <- split_recs_scr$phen.trait.2
split_recs_scr$phen.trait.2ga[split_recs_scr$phen.trait.2=="germination time"] <- "0gt"
split_recs_scr$phen.trait.2fda <- split_recs_scr$phen.trait.2
split_recs_scr$phen.trait.2fda[split_recs_scr$phen.trait.2=="flowering duration"] <- "0fd"

wosbin <- function(dataframe, woscolname,nrecs){
	dataframe$dummyWOSrelrank <- dataframe[,woscolname]
	dataframe$dummyWOSrelrank[dataframe$dummyWOSrelran%in%names(table(dataframe$WOS.relevance.rank)[table(dataframe$WOS.relevance.rank)<nrecs])] <- "binned"
	return(dataframe)
}

flowering <- wosbin(split_recs_scr[split_recs_scr$phen.trait.2=="flowering time",],"WOS.relevance.rank",13)
	#consider removing: seed = 6  *** BUT NOTE THAT SOMETIMES THERE IS STILL ONLY ONE STUDY PER CATEGORY ***
	floweringLOC <- wosbin(split_recs_scr[split_recs_scr$phen.trait.2=="flowering time" & split_recs_scr$micrloc3 != "seed",],"WOS.relevance.rank",13)
	floweringMAT <- wosbin(split_recs_scr[split_recs_scr$phen.trait.2=="flowering time" & split_recs_scr$mating2 != "unk" & split_recs_scr$mating2 != "s+o",],"WOS.relevance.rank",13)
	#consider removing: phytoH = 6
	floweringMPH <- wosbin(split_recs_scr[split_recs_scr$phen.trait.2=="flowering time" & !split_recs_scr$MicrobeMPH%in%c( "unknown","phytohormones"),],"WOS.relevance.rank",13)
germtime <- wosbin(split_recs_scr[split_recs_scr$phen.trait.2=="germination time",],"WOS.relevance.rank",13)
	#consider removing: MF size = 1
	germtimeTAX <- wosbin(split_recs_scr[split_recs_scr$phen.trait.2=="germination time" & split_recs_scr$micrtax.sim2 !="MF",],"WOS.relevance.rank",13)
	#consider removing: shoot = 4 , var incl reprod = 8, var not incl reprod = 12
	germtimeLOC <- wosbin(split_recs_scr[split_recs_scr$phen.trait.2=="germination time" & split_recs_scr$micrloc3 != "shoot",],"WOS.relevance.rank",13)
	germtimeMAT <- wosbin(split_recs_scr[split_recs_scr$phen.trait.2=="germination time" & split_recs_scr$mating2 != "unk" & split_recs_scr$mating2 != "s+o",],"WOS.relevance.rank",13)
	#phytoH = 3 path = 5
	germtimeMPH <- wosbin(split_recs_scr[split_recs_scr$phen.trait.2=="germination time" & !split_recs_scr$MicrobeMPH %in%c("unknown","pathogen","phytohormones"),],"WOS.relevance.rank",13)
germprb <- wosbin(split_recs_scr[split_recs_scr$phen.trait.2=="germination prb",],"WOS.relevance.rank",13)
	#consider removing: MF, size = 1
	germprbTAX <- wosbin(split_recs_scr[split_recs_scr$phen.trait.2=="germination prb" & split_recs_scr$micrtax.sim2 !="MF",],"WOS.relevance.rank",13)
	germprbMAT <- wosbin(split_recs_scr[split_recs_scr$phen.trait.2=="germination prb" & split_recs_scr$mating2 != "unk" & split_recs_scr$mating2 != "s+o",],"WOS.relevance.rank",13)
	germprbMPH <- wosbin(split_recs_scr[split_recs_scr$phen.trait.2=="germination prb" & split_recs_scr$MicrobeMPH != "unknown",],"WOS.relevance.rank",13)
	#consider removing: strain mix = 7
	germprbSCM <- wosbin(split_recs_scr[split_recs_scr$phen.trait.2=="germination prb" & split_recs_scr$Strainvcomm.1 != "strain mix",],"WOS.relevance.rank",13)# allphen <- split_recs_scr[!is.na(split_recs_scr$phen.trait.2) & !(split_recs_scr$phen.trait.2%in% c("floral bud duration","flowering duration","germination prb")),]
	#consider removing: shoot = 4 
	germprbLOC <- wosbin(split_recs_scr[split_recs_scr$phen.trait.2=="germination prb" & split_recs_scr$micrloc3 != "shoot",],"WOS.relevance.rank",13)

#consider removing, mat time=4, senec time = 2, peak flower = 9, flower sen. time = 4; floral bud duration = 9
allphen <- split_recs_scr[!is.na(split_recs_scr$phen.trait.2) & (split_recs_scr$phen.trait.2%in% c("budburst","floral budset time","flowering time","phyllochron","fruiting time","germination time")),]
	allphen$dummyWOSrelrank <- allphen$WOS.relevance.rank
	allphen$dummyWOSrelrank[allphen$dummyWOSrelran%in%names(table(allphen$WOS.relevance.rank)[table(allphen$WOS.relevance.rank)<13])] <- "binned"
allduration <- split_recs_scr[!is.na(split_recs_scr$phen.trait.2) & (split_recs_scr$phen.trait.2%in% c("flowering duration","germination prb")),]
	allduration$dummyWOSrelrank <- allduration$WOS.relevance.rank
	allduration$dummyWOSrelrank[allduration$dummyWOSrelran%in%names(table(allduration$WOS.relevance.rank)[table(allduration$WOS.relevance.rank)<13])] <- "binned"


###for numbers in text:
flowering$anyeffnum <- as.numeric(as.factor(flowering$anyeffect))-1
sum(tapply(flowering$anyeffnum,flowering$WOS.relevance.rank,sum) == 0) #6
length(tapply(flowering$anyeffnum,flowering$WOS.relevance.rank,sum))#39
germtime$anyeffnum <- as.numeric(as.factor(germtime$anyeffect))-1
sum(tapply(germtime$anyeffnum,germtime$WOS.relevance.rank,sum) ==0)#1 -- note this record only appears in germtime, not germprb
length(tapply(germtime$anyeffnum,germtime$WOS.relevance.rank,sum))#19
germprb$anyeffnum <- as.numeric(as.factor(germprb$anyeffect))-1
sum(tapply(germprb$anyeffnum,germprb$WOS.relevance.rank,sum) ==0)#0
length(tapply(germprb$anyeffnum,germprb$WOS.relevance.rank,sum))#13
length(unique(c(germprb$WOS.relevance.rank, germtime$WOS.relevance.rank)))#overlap


#NOTE R MUST be fixed. J Hadfield states, "The residual variance is not  identifiable in the likelihood for binary data:I have tried to  explain this intuitively in Section 2.6 of the CourseNotes using a  tasteless example of hospital deaths."
priornr=list(R=list(V= 1, fix=1))
prior=list(R=list(V= 1, fix=1), G=list(G1=list(V=1, nu=0)))
priornuup=list(R=list(V= 1, fix=1), G=list(G1=list(V=1, nu=8)))
	#turning this degree of belief up to solve sampling issues for random effects when studies have few tests, or tests only in one category (in which case, this biases towards all the variance being due to the treatment, rather than study)
	#estimates for intercept and random effects became more independent, because this shunts variation to be explained by fixed effect when there isn't information in the data for the random effect
priornuup2=list(R=list(V= 1, fix=1), G=list(G1=list(V=1, nu=8),G2=list(V=1, nu=8)))
#after fitting the flowering model with various priors for the case with random effects(prior 2 onwards, the latter two of which are meant for no-intercept models); and various levels of binning,

######BINOMIALS
 durdiffs <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~phen.trait.2,random = ~dummyWOSrelrank,data=allduration,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
 durDdiffs <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~phen.trait.2,random = ~dummyWOSrelrank,data=allduration,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
###NO DIFFERENCES. so no figure.
 traitdiffs <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~phen.trait.2,random = ~dummyWOSrelrank,data=allphen,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
	 traitdiffs2 <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~phen.trait.2fba,random = ~dummyWOSrelrank,data=allphen,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
	 traitdiffs3 <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~phen.trait.2fla,random = ~dummyWOSrelrank,data=allphen,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
	 traitdiffs4 <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~phen.trait.2fra,random = ~dummyWOSrelrank,data=allphen,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
	 traitdiffs5 <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~phen.trait.2ga,random = ~dummyWOSrelrank,data=allphen,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
 traitDdiffs <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~phen.trait.2,random = ~dummyWOSrelrank,data=allphen,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
	 traitDdiffs2 <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~phen.trait.2fba,random = ~dummyWOSrelrank,data=allphen,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
	 traitDdiffs3 <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~phen.trait.2fla,random = ~dummyWOSrelrank,data=allphen,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
	 traitDdiffs4 <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~phen.trait.2fra,random = ~dummyWOSrelrank,data=allphen,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
	 traitDdiffs5 <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~phen.trait.2ga,random = ~dummyWOSrelrank,data=allphen,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
###NO DIFFERENCES. so no figure.


flowerbitax <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~micrtax.sim2,random = ~dummyWOSrelrank,data=flowering,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	flowerbitax2 <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~micrtaxmfa,random = ~dummyWOSrelrank,data=flowering,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	flowerbitax3 <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~micrtaxmixa,random = ~dummyWOSrelrank,data=flowering,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
flowerbiscom <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~Strainvcomm.1,random = ~dummyWOSrelrank,data=flowering,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	flowerbiscom2 <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~strainvcommsa,random = ~dummyWOSrelrank,data=flowering,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
flowerbimat <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~mating2,random = ~dummyWOSrelrank,data=floweringMAT,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
flowerbilh <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~lifeform2,random = ~dummyWOSrelrank,data=flowering,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	flowerbilh2 <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~lifeform2pa,random = ~dummyWOSrelrank,data=flowering,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
flowerbimph <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~MicrobeMPH,random = ~dummyWOSrelrank,data=floweringMPH,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	flowerbimph2 <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~mphba,random = ~dummyWOSrelrank,data=floweringMPH,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
flowerbiloc <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~micrloc3,random = ~dummyWOSrelrank,data=floweringLOC,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	flowerbiloc2 <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~micrloc3shoota,random = ~dummyWOSrelrank,data=floweringLOC,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
#
flowerDbitax <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~micrtax.sim2,random = ~dummyWOSrelrank,data=flowering,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	flowerDbitax2 <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~micrtaxmfa,random = ~dummyWOSrelrank,data=flowering,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	flowerDbitax3 <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~micrtaxmixa,random = ~dummyWOSrelrank,data=flowering,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
flowerDbiscom <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~Strainvcomm.1,random = ~dummyWOSrelrank,data=flowering,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	flowerDbiscom2 <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~strainvcommsa,random = ~dummyWOSrelrank,data=flowering,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
flowerDbimat <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~mating2,random = ~dummyWOSrelrank,data=floweringMAT,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
flowerDbilh <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~lifeform2,random = ~dummyWOSrelrank,data=flowering,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	flowerDbilh2 <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~lifeform2pa,random = ~dummyWOSrelrank,data=flowering,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
flowerDbimph <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~MicrobeMPH,random = ~dummyWOSrelrank,data=floweringMPH, family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	flowerDbimph2 <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~mphba,random = ~dummyWOSrelrank,data=floweringMPH, family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
flowerDbiloc <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~micrloc3,random = ~dummyWOSrelrank,data=floweringLOC,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	flowerDbiloc2 <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~micrloc3shoota,random = ~dummyWOSrelrank,data=floweringLOC,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)

germtimebitax <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~micrtax.sim2,random = ~dummyWOSrelrank,data=germtimeTAX,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	germtimebitax2 <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~micrtaxmixa,random = ~dummyWOSrelrank,data=germtimeTAX,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
germtimebiscom <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~Strainvcomm.1,random = ~dummyWOSrelrank,data=germtime,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	germtimebiscom2 <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~strainvcommsa,random = ~dummyWOSrelrank,data=germtime,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
germtimebimat <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~mating2,random = ~dummyWOSrelrank,data=germtimeMAT,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
germtimebimph <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~MicrobeMPH,random = ~dummyWOSrelrank,data=germtimeMPH,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
germtimebilh <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~lifeform2,random = ~dummyWOSrelrank,data=germtime,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	germtimebilh2 <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~lifeform2pa,random = ~dummyWOSrelrank,data=germtime,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
germtimebiloc <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~micrloc3,random = ~dummyWOSrelrank,data=germtimeLOC,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	germtimebiloc2 <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~micrloc3seeda,random = ~dummyWOSrelrank,data=germtimeLOC,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
#
germtimeDbitax <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~micrtax.sim2,random = ~dummyWOSrelrank,data=germtimeTAX,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	germtimeDbitax2 <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~micrtaxmixa,random = ~dummyWOSrelrank,data=germtimeTAX,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
germtimeDbiscom <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~Strainvcomm.1,random = ~dummyWOSrelrank,data=germtime,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	germtimeDbiscom2 <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~strainvcommsa,random = ~dummyWOSrelrank,data=germtime,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
germtimeDbimat <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~mating2,random = ~dummyWOSrelrank,data=germtimeMAT,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
germtimeDbimph <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~MicrobeMPH,random = ~dummyWOSrelrank,data=germtimeMPH,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
germtimeDbilh <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~lifeform2,random = ~dummyWOSrelrank,data=germtime,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	germtimeDbilh2 <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~lifeform2pa,random = ~dummyWOSrelrank,data=germtime,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
germtimeDbiloc <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~micrloc3,random = ~dummyWOSrelrank,data=germtimeLOC,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	germtimeDbiloc2 <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~micrloc3seeda,random = ~dummyWOSrelrank,data=germtimeLOC,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)

germprbbitax <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~micrtax.sim2,random = ~dummyWOSrelrank,data=germprbTAX,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	germprbbitax2 <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~micrtaxmixa,random = ~dummyWOSrelrank,data=germprbTAX,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
germprbbiscom <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~Strainvcomm.1,random = ~dummyWOSrelrank,data=germprbSCM,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
germprbbimat <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~mating2,random = ~dummyWOSrelrank,data=germprbMAT,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
germprbbimph <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~MicrobeMPH,random = ~dummyWOSrelrank,data=germprbMPH,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	germprbbimph2 <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~mphba,random = ~dummyWOSrelrank,data=germprbMPH,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
germprbbilh <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~lifeform2,random = ~dummyWOSrelrank,data=germprb,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	germprbbilh2 <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~lifeform2pa,random = ~dummyWOSrelrank,data=germprb,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
germprbbiloc <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~micrloc3,random = ~dummyWOSrelrank,data=germprbLOC,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	germprbbiloc2 <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~micrloc3seeda,random = ~dummyWOSrelrank,data=germprbLOC,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
#
germprbDbitax <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~micrtax.sim2,random = ~dummyWOSrelrank,data=germprbTAX,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	germprbDbitax2 <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~micrtaxmixa,random = ~dummyWOSrelrank,data=germprbTAX,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
germprbDbiscom <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~Strainvcomm.1,random = ~dummyWOSrelrank,data=germprbSCM,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
germprbDbimat <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~mating2,random = ~dummyWOSrelrank,data=germprbMAT,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
germprbDbimph <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~MicrobeMPH,random = ~dummyWOSrelrank,data=germprbMPH,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	germprbDbimph2 <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~mphba,random = ~dummyWOSrelrank,data=germprbMPH,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
germprbDbilh <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~lifeform2,random = ~dummyWOSrelrank,data=germprb,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	germprbDbilh2 <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~lifeform2pa,random = ~dummyWOSrelrank,data=germprb,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
germprbDbiloc <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~micrloc3,random = ~dummyWOSrelrank,data=germprbLOC,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	germprbDbiloc2 <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~micrloc3seeda,random = ~dummyWOSrelrank,data=germprbLOC,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)



weightmns2 <- function(dataframelist, cols,dfwithnames){
	typelist <- lapply(1:length(cols), function(x) sort(unique( dfwithnames[ ,cols[x] ] )))
# 	lengthlist <- lapply(typelist, length)
	wos <- lapply(dataframelist, function(dat) sort(unique(dat$dummyWOSrelrank) ) )#random effects in model are automatically sorted
	xEmn <- lapply(1:length(typelist), function(type) sapply(typelist[[type]], function(m)  sapply(wos[[type]], function(w)  
			mean(dataframelist[[type]]$isEarlyOrNrw[dataframelist[[type]]$dummyWOSrelrank==w & dataframelist[[type]][,which(colnames(dataframelist[[type]])==cols[type])]==m] )
			   )))
	xE.m <- unlist( lapply(1:length(xEmn), function(T)   colMeans(xEmn[[T]],na.rm=T)))
	xEse <- lapply(1:length(typelist), function(type) sapply(typelist[[type]], function(m)  sapply(wos[[type]], function(w)  
			std.error(dataframelist[[type]]$isEarlyOrNrw[dataframelist[[type]]$dummyWOSrelrank==w & dataframelist[[type]][,which(colnames(dataframelist[[type]])==cols[type])]==m] )  
			 )))
	xE.se <-unlist(lapply(1:length(xEse), function(type)  sapply(1:ncol(xEse[[type]]), function(m) mean(xEse[[type]][,m],na.rm=T))))
 	xDmn <-  lapply(1:length(typelist), function(type) sapply(typelist[[type]], function(m)  sapply(wos[[type]], function(w)  
			mean(dataframelist[[type]]$isDelayOrExp[dataframelist[[type]]$dummyWOSrelrank==w & dataframelist[[type]][,which(colnames(dataframelist[[type]])==cols[type])]==m] )
			   )))
	xDse <- lapply(1:length(typelist), function(type) sapply(typelist[[type]], function(m)  sapply(wos[[type]], function(w)  
			std.error(dataframelist[[type]]$isDelayOrExp[dataframelist[[type]]$dummyWOSrelrank==w & dataframelist[[type]][,which(colnames(dataframelist[[type]])==cols[type])]==m] )  
			 )))
	xD.m <- unlist(lapply(1:length(xDmn), function(type) colMeans(xDmn[[type]],na.rm=T)))
	xD.se <-unlist(lapply(1:length(xDse), function(type)  sapply(1:ncol(xDse[[type]]), function(m) mean(xDse[[type]][,m],na.rm=T))))
	return(list(xE.m,xE.se,xD.m,xD.se))
}

cols1 <- c("MicrobeMPH","micrtax.sim2","micrloc3") #"Strainvcomm.1","lifeform2",
flwrdlist <- list(flowering,flowering,flowering)
germtlist <- list(germtime,germtime,germtime)
germplist <- list(germprb,germprb,germprb)

flwrE.m <-   weightmns2(flwrdlist,cols=cols1,flowering)[[1]][-5] 
flwrE.se <-  weightmns2(flwrdlist,cols=cols1,flowering)[[2]][-5]
flwrD.m <-   weightmns2(flwrdlist,cols=cols1,flowering)[[3]][-5]
flwrD.se <-  weightmns2(flwrdlist,cols=cols1,flowering)[[4]][-5]
germTE.m <-  weightmns2(germtlist,cols=cols1,germtime)[[1]][-5] 
germTE.se <- weightmns2(germtlist,cols=cols1,germtime)[[2]][-5] 
germTD.m <-  weightmns2(germtlist,cols=cols1,germtime)[[3]][-5] 
germTD.se <- weightmns2(germtlist,cols=cols1,germtime)[[4]][-5] 
germPE.m <-  weightmns2(germplist,cols=cols1,germprb)[[1]] 
germPE.se <- weightmns2(germplist,cols=cols1,germprb)[[2]] 
germPD.m <-  weightmns2(germplist,cols=cols1,germprb)[[3]] 
germPD.se <- weightmns2(germplist,cols=cols1,germprb)[[4]] 
#subtracting unknown category in MPH, which is not useful. leaving in germP to fill gap from phytohormones

#fill in from model results
siggroupsFE <- c("a","b","ab", "-", #,"ab" #mph , at 0.9, pathogen marginally less than nutrients (is b only)
				"a","a","ab","b",#tax #at 0.9 mixed is marginally less than MF but not than bacteria
# 				"a","b","b","a","a")#, #loc
				"a","a","-","b")#, #loc
siggroupsFD <- c(" "," "," ","-", #"b", #mph
				"b","b","a","ab",#tax
				" "," ","-"," ")#, #loc
# 				"a","b","a","a","a")#, #loc
siggroupsGtE <- c( "b","a","-","-",#"bc","a", #mph
# 				"b","a","c","b",#tax  
				"a","-","b","a",#tax  
				"","","","-")#, #loc	
siggroupsGtD <- c(" "," ","-","-",#"a","c", #mph
# 				"b","c","a","a",#tax
				"b","-","a","a",#tax
				" "," "," ","-")#,#loc  
# 				"a","a","b","a","b")#,#loc 
siggroupsGpE <- c("b","b","a","-", #mph ****no phytohormone records***
				"b","-","a","a",#tax
				" "," "," ","-")#, #loc
# 				"a","a","b","a","a")#, #loc
siggroupsGpD <- c("b","a","b","-", #mph ****no phytohormone records***
				"a","-","ab","b",#tax
# 				"b","a","bc","c",#tax
				" "," "," ","-")#, #loc
# 				"a","a","b","a","a")#, #loc

pdf("means_ses_and prelim fitted diffs from binom slim.pdf",height=6,width=4)
xvals <- c(seq(from=0,to=1,length.out=c(4+4+4) ))
# xlab <- c("nutrients","other beneficial","pathogen","phytohormones","bacteria","myc fungi","mixed","other fungi","root","seed","shoot","mult. w/ reprod","mult. w/o reprod")
xlab <- c("nutrients","other beneficial","pathogen","phytohormones","bacteria","myc fungi","mixed","other fungi","multiple","root","seed","shoot")
par(mfrow=c(3,1))
par(mar=c(1,0,1,0))
par(oma=c(7,4,1,1))
plot(c(0,1)~I(c(0,1)),pch=NA,xlab="",ylab="",bty="n",xaxt="n",ylim=c(-0.3,1.2),yaxt="n")
	points(unlist(flwrE.m)~xvals,pch=16)
	points(unlist(flwrD.m)~I(xvals+0.01),pch=16,col="gray")
	arrows(x0=xvals , y0=flwrE.m -flwrE.se, y1 = flwrE.m +flwrE.se,length=0     )
	arrows(x0=xvals+0.01 , y0=flwrD.m -flwrD.se, y1 = flwrD.m +flwrD.se,length=0  ,col="gray"   )
	abline(v= c(sum(xvals[4:5])/2, sum(xvals[8:9])/2 ),lty=3 )
	text(xvals,y=1.2,siggroupsFE)
	text(xvals,y=-0.2,siggroupsFD,col="gray")
	mtext(side=3,"flowering time")
	axis(side=2,at=c(0,0.25,0.5,0.75,1),labels=c(0," ",0.5, " ",1))
legend(xvals[10],y=1.2,c("Earlier","Delay"),fill=c("black","gray"),bty="n")
plot(c(0,1)~I(c(0,1)),pch=NA,xlab="",ylab="",bty="n",xaxt="n",ylim=c(-0.3,1.2),yaxt="n")
	points(unlist(germTE.m)~xvals,pch=16)
	points(unlist(germTD.m)~I(xvals+0.01),pch=16,col="gray")
	arrows(x0=xvals , y0=germTE.m -germTE.se, y1 = germTE.m +germTE.se,length=0     )
	arrows(x0=xvals+0.01 , y0=germTD.m -germTD.se, y1 = germTD.m +germTD.se,length=0  ,col="gray"   )
	abline(v= c(sum(xvals[4:5])/2, sum(xvals[8:9])/2),lty=3 )
	text(xvals,y=1.2,siggroupsGtE)
	text(xvals,y=-0.2,siggroupsGtD,col="gray")
	mtext(side=3,"germination time")
	mtext(side=2,"Prb of observing Effect",line=2)
	axis(side=2,at=c(0,0.25,0.5,0.75,1),labels=c(0," ",0.5, " ",1))
legend(xvals[10],y=1.2,c("Earlier","Delay"),fill=c("black","gray"),bty="n")
plot(c(0,1)~I(c(0,1)),pch=NA,xlab="",ylab="",bty="n",xaxt="n",ylim=c(-0.3,1.2),yaxt="n")
	points(unlist(germPE.m)~xvals,pch=16)
	points(unlist(germPD.m)~I(xvals+0.01),pch=16,col="gray")
	arrows(x0=xvals, y0=germPE.m -germPE.se, y1 = germPE.m +germPE.se,length=0     )
	arrows(x0=xvals+0.01 , y0=germPD.m -germPD.se, y1 = germPD.m +germPD.se,length=0  ,col="gray"   )
	abline(v= c(sum(xvals[4:5])/2, sum(xvals[8:9])/2),lty=3 )
	text(xvals,y=1.2,siggroupsGpE)
	text(xvals,y=-0.2,siggroupsGpD,col="gray")
	mtext(side=3,"germination probability")
	axis(side=2,at=c(0,0.25,0.5,0.75,1),labels=c(0," ",0.5, " ",1))
legend(xvals[10],y=1.2,c("Narrowed","Expanded"),fill=c("black","gray"),bty="n")
axis(side=1,at=xvals, labels=xlab,las=2)
dev.off()


###Supplemental figure. (replaces variables for easier code, watchout if running piecewise)

cols1 <- c("Strainvcomm.1","lifeform2","mating2") #"Strainvcomm.1","lifeform2",
flwrdlist <- list(flowering,flowering,floweringMAT)
germtlist <- list(germtime,germtime,germtimeMAT)
germplist <- list(germprb,germprb,germprbMAT)

flwrE.m <-   weightmns2(flwrdlist,cols=cols1,flowering)[[1]][-c(8,10)] 
flwrE.se <-  weightmns2(flwrdlist,cols=cols1,flowering)[[2]][-c(8,10)] 
flwrD.m <-   weightmns2(flwrdlist,cols=cols1,flowering)[[3]][-c(8,10)] 
flwrD.se <-  weightmns2(flwrdlist,cols=cols1,flowering)[[4]][-c(8,10)] 
germTE.m <-  weightmns2(germtlist,cols=cols1,germtime)[[1]] [-9]
germTE.se <- weightmns2(germtlist,cols=cols1,germtime)[[2]] [-9]
germTD.m <-  weightmns2(germtlist,cols=cols1,germtime)[[3]] [-9]
germTD.se <- weightmns2(germtlist,cols=cols1,germtime)[[4]] [-9]
germPE.m <-  weightmns2(germplist,cols=cols1,germprb)[[1]] [-9]
germPE.se <- weightmns2(germplist,cols=cols1,germprb)[[2]] [-9]
germPD.m <-  weightmns2(germplist,cols=cols1,germprb)[[3]] [-9]
germPD.se <- weightmns2(germplist,cols=cols1,germprb)[[4]] [-9]

#fill in from model results
siggroupsFE <- c("","","",#scom
 				"a","ab","b",#lh  if 0.9, definitely,  at .95 sometimes there is no difference...: "a","ab","b"; doesn't look perfect
 				"","") #mat
siggroupsFD <- c("a","b","b",#scom;  looks backwards!
				"","","",#lh
				"","") #mat
siggroupsGtE <- c("","","",#scom
				"b","ab","a",#lh
				"","") #mat
siggroupsGtD <- c("a","b","ab",#scom # at 0.9, strain mix is b only, recall is comm, strain, strain mix
				"","","",#lh
				"","") #mat
siggroupsGpE <- c(" "," ","-",#scom
				"","","",#lh
				"","") #mat
siggroupsGpD <- c(" "," ","-",#scom ; recall is comm, strain, strain mix
				"ab","b","a",#lh
				"b","a") # mat looks backwards

pdf("means_ses_and prelim fitted diffs from binom supp cats.pdf",height=6,width=3)
xvals <- c(seq(from=0,to=1,length.out=c(3+3+2) ))
xlab <- c("community","strain","strain mix","annual","both","perennial","outcrosser","selfer")
par(mfrow=c(3,1))
par(mar=c(1,0,1,0))
par(oma=c(7,4,1,1))
plot(c(0,1)~I(c(0,1)),pch=NA,xlab="",ylab="",bty="n",xaxt="n",ylim=c(-0.3,1.2),yaxt="n")
	points(unlist(flwrE.m)~xvals,pch=16)
	points(unlist(flwrD.m)~I(xvals+0.01),pch=16,col="gray")
	arrows(x0=xvals , y0=flwrE.m -flwrE.se, y1 = flwrE.m +flwrE.se,length=0     )
	arrows(x0=xvals+0.01 , y0=flwrD.m -flwrD.se, y1 = flwrD.m +flwrD.se,length=0  ,col="gray"   )
	abline(v= c(sum(xvals[3:4])/2, sum(xvals[6:7])/2 ),lty=3 )
	text(xvals,y=1.2,siggroupsFE)
	text(xvals,y=-0.2,siggroupsFD,col="gray")
	mtext(side=3,"flowering time")
	axis(side=2,at=c(0,0.25,0.5,0.75,1),labels=c(0," ",0.5, " ",1))
legend(xvals[4],y=1.2,c("Earlier","Delay"),fill=c("black","gray"),bty="n")
plot(c(0,1)~I(c(0,1)),pch=NA,xlab="",ylab="",bty="n",xaxt="n",ylim=c(-0.3,1.2),yaxt="n")
	points(unlist(germTE.m)~xvals,pch=16)
	points(unlist(germTD.m)~I(xvals+0.01),pch=16,col="gray")
	arrows(x0=xvals , y0=germTE.m -germTE.se, y1 = germTE.m +germTE.se,length=0     )
	arrows(x0=xvals+0.01 , y0=germTD.m -germTD.se, y1 = germTD.m +germTD.se,length=0  ,col="gray"   )
	abline(v= c(sum(xvals[3:4])/2, sum(xvals[6:7])/2 ),lty=3 )
	text(xvals,y=1.2,siggroupsGtE)
	text(xvals,y=-0.2,siggroupsGtD,col="gray")
	mtext(side=3,"germination time")
	mtext(side=2,"Prb of observing Effect",line=2)
	axis(side=2,at=c(0,0.25,0.5,0.75,1),labels=c(0," ",0.5, " ",1))
legend(xvals[4],y=1.2,c("Earlier","Delay"),fill=c("black","gray"),bty="n")
plot(c(0,1)~I(c(0,1)),pch=NA,xlab="",ylab="",bty="n",xaxt="n",ylim=c(-0.3,1.2),yaxt="n")
	points(unlist(germPE.m)~xvals,pch=16)
	points(unlist(germPD.m)~I(xvals+0.01),pch=16,col="gray")
	arrows(x0=xvals, y0=germPE.m -germPE.se, y1 = germPE.m +germPE.se,length=0     )
	arrows(x0=xvals+0.01 , y0=germPD.m -germPD.se, y1 = germPD.m +germPD.se,length=0  ,col="gray"   )
	abline(v= c(sum(xvals[3:4])/2, sum(xvals[6:7])/2 ),lty=3 )
	text(xvals,y=1.2,siggroupsGpE)
	text(xvals,y=-0.2,siggroupsGpD,col="gray")
	mtext(side=3,"germination probability")
	axis(side=2,at=c(0,0.25,0.5,0.75,1),labels=c(0," ",0.5, " ",1))
legend(xvals[4],y=1.2,c("Narrowed","Expanded"),fill=c("black","gray"),bty="n")
axis(side=1,at=xvals, labels=xlab,las=2)
dev.off()


##get means, but adjusting for by study random effect category
### the goal is to extract and subtract from the raw means the estimated study effect
	#because the effects are estimated on the latent scale, we have to indirectly get the adjustment for the data on raw scale
		#we get the difference between the mean prediction at each response variable category at the mean of the random effect distribution (0) from each actual study bin's random effect
		#this table gives the amount we have estimated that our observed data is inflated/deflated due to the random effects
		#we then apply this adjustment to the raw means by study & treatment, and take the row means.
		#we have not (and cannot with this data) estimate effects of studies on the variance in the data, thus there is no adjustment to the SE
# weightmns <- function(dataframe, cols,modslistE,modslistD){
# 	typelist <- lapply(cols, function(x) sort(unique(dataframe[,x])))
# 	wos <- sort(unique(dataframe$dummyWOSrelrank))#random effects in model are automatically sorted
# 	xEmn <- lapply(1:length(typelist), function(type) sapply(typelist[[type]], function(m)  sapply(wos, function(w)  mean(dataframe$isEarlyOrNrw[dataframe$dummyWOSrelrank==w & dataframe[,which(colnames(dataframe)==cols[type])]==m] )   )))
# 	ltrts <- lapply(typelist,length)
# 	if(ltrts>2){  #bad code , but works because the only time I use this function when ltrts is less than two, there's only one model, so it is first.
# 		prdwrE <-lapply(1:length(typelist), function(T) sapply( (ltrts[[T]]+1):(ltrts[[T]]+length(wos)), function(x) 
# 				invlogistic(mean(modslistE[[T]]$Sol[,1])  + mean(modslistE[[T]]$Sol[,x])+c(0,colMeans(modslistE[[T]]$Sol[,2:ltrts[[T]]])) )    )   )
# 		pbaseE <-lapply(1:length(typelist), function(T) 
# 				invlogistic(mean(modslistE[[T]]$Sol[,1])  +c(0,colMeans(modslistE[[T]]$Sol[,2:ltrts[[T]]])) )       )
# 	} else { #assume at least two treatments!
# 		prdwrE <-lapply(1:length(typelist), function(T) sapply( (ltrts[[T]]+1):(ltrts[[T]]+length(wos)), function(x) 
# 				invlogistic(mean(modslistE[[T]]$Sol[,1])  + mean(modslistE[[T]]$Sol[,x])+c(0,mean(modslistE[[T]]$Sol[,2])) )    )   )
# 		pbaseE <-lapply(1:length(typelist), function(T) 
# 				invlogistic(mean(modslistE[[T]]$Sol[,1])  +c(0,mean(modslistE[[T]]$Sol[,2])) )       )
# 	}
# 	xE.m <- unlist( lapply(1:length(prdwrE), function(T)   rowMeans(t(xEmn[[T]]) - (prdwrE[[T]] - pbaseE[[T]]), na.rm=T) ) )
# 	xEse <- lapply(1:length(typelist), function(type) sapply(typelist[[type]], function(m)  sapply(wos, function(w)  std.error(dataframe$isEarlyOrNrw[dataframe$dummyWOSrelrank==w & dataframe[,which(colnames(dataframe)==cols[type])]==m] )   )))
# 	xE.se <-unlist(lapply(1:length(xEse), function(type)  sapply(1:ncol(xEse[[type]]), function(m) mean(xEse[[type]][,m],na.rm=T))))
# # 	xE.m <- unlist(lapply(1:length(xEmn), function(type) colMeans(xEmn[[type]],na.rm=T)))
# 	if(ltrts>2){
# 		prdwrD <-lapply(1:length(typelist), function(T) sapply( (ltrts[[T]]+1):(ltrts[[T]]+length(wos)), function(x) 
# 				invlogistic(mean(modslistD[[T]]$Sol[,1])  + mean(modslistD[[T]]$Sol[,x])+c(0,colMeans(modslistD[[T]]$Sol[,2:ltrts[[T]]])) )    )   )
# 		pbaseD <-lapply(1:length(typelist), function(T) 
# 				invlogistic(mean(modslistD[[T]]$Sol[,1])  +c(0,colMeans(modslistD[[T]]$Sol[,2:ltrts[[T]]])) )       )
# 	} else { #assume at least two treatments!
# 		prdwrD <-lapply(1:length(typelist), function(T) sapply( (ltrts[[T]]+1):(ltrts[[T]]+length(wos)), function(x) 
# 				invlogistic(mean(modslistD[[T]]$Sol[,1])  + mean(modslistD[[T]]$Sol[,x])+c(0,mean(modslistD[[T]]$Sol[,2])) )    )   )
# 		pbaseD <-lapply(1:length(typelist), function(T) 
# 				invlogistic(mean(modslistD[[T]]$Sol[,1])  +c(0,mean(modslistD[[T]]$Sol[,2])) )       )
# 	}
#  	xDmn <- lapply(1:length(typelist), function(type) sapply(typelist[[type]], function(m)  sapply(wos, function(w)  mean(dataframe$isDelayOrExp[dataframe$dummyWOSrelrank==w & dataframe[,which(colnames(dataframe)==cols[type])]==m] )   )))
# 	xD.m <- unlist( lapply(1:length(prdwrD), function(T)   rowMeans(t(xDmn[[T]]) - (prdwrD[[T]] - pbaseD[[T]]), na.rm=T) ) )
# 	xDse <- lapply(1:length(typelist), function(type) sapply(typelist[[type]], function(m)  sapply(wos, function(w)  std.error(dataframe$isDelayOrExp[dataframe$dummyWOSrelrank==w & dataframe[,which(colnames(dataframe)==cols[type])]==m] )   )))
# # 	xD.m <- unlist(lapply(1:length(xDmn), function(type) colMeans(xDmn[[type]],na.rm=T)))
# 	xD.se <-unlist(lapply(1:length(xDse), function(type)  sapply(1:ncol(xDse[[type]]), function(m) mean(xDse[[type]][,m],na.rm=T))))
# 	return(list(xE.m,xE.se,xD.m,xD.se))
# }
