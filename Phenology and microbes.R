library(UpSetR) #for set plots
library(MCMCglmm) #for models
library(grDevices) #for color


##define useful functions
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

cwos_recs <- read.csv("saved records WOS friday sep 11 2020 - 500scoredWOSrecords.csv",
					stringsAsFactors=F,na.strings = c("NA","")) 
split_recs <- read.csv("saved records WOS friday sep 11 2020 - Split Records.csv", #sep="\t",
					stringsAsFactors=F,na.strings = c("NA","")) 
split_recs$mating[which(is.na(split_recs$mating) & split_recs$WOS.relevance.rank==28) ] <- "unk" # mating system score of unknown was left out for one record
split_recs_scr <- split_recs[!is.na(split_recs$direction.effect.split),] 
	#remove several records that were not possible to score for effect direction, but that were relevant to microbial effects on phenology

################################################################################
####Analysis of all broadly relevant records			
################################################################################


relevant <- cwos_recs$"Relevant..Y.N"=="y"
PMo <- as.numeric(as.factor(cwos_recs$otherPhenMicrLink=="y")) -1 #other link between phenology, plants and microbes
MiP <- as.numeric(as.factor(cwos_recs$MicrobesAffectPhenology=="y")) -1 # microbes influence plant phenology

table(relevant)
table(PMo)
table(MiP)
table(cwos_recs$anyeffect)

wos_rel <- cwos_recs[which(relevant),]
wos_MiP <- cwos_recs[which(MiP==1),]

taxand <- gsub("; ", wos_rel$micrtax,replacement="&")
locand <- gsub("; ", wos_rel$micrloc,replacement="&")
phenand <- gsub("; ", wos_rel$phentrait_of,replacement="&")
#replace abbreviation "prb" from object phenand
phenand[phenand == "germination prb"] <- "germination probability"
phenand[phenand == "germination time&germination prb"] <- "germination time&germination probability"
phenand[phenand == "germination prb&flowering time"] <- "germination probability&flowering time"

pdf("upsetrtax.pdf",width=5,height=5)
upset(fromExpression(table(taxand)),nsets=15, order.by = "freq", mainbar.y.max = 45,
	mainbar.y.label="# with microbe taxonomy intersect", sets.x.label = "# with microbe taxonomy")
dev.off()
pdf("upsetrloc.pdf",width=4.5,height=5)
upset(fromExpression(table(locand)),nsets=10, order.by = "freq", mainbar.y.max=105,
	mainbar.y.label="# with microbes in tissue/s in intersect", sets.x.label = "# with microbes in tissue")
dev.off()
pdf("upsetrtrt.pdf",width=4.5,height=5)
upset(fromExpression(table(phenand)), nsets=15, order.by = "freq", mainbar.y.max = 25.5,
	mainbar.y.label="# measured phenophase/s in intersect", sets.x.label = "# measured phenophase")
dev.off()


################################################################################
####Analysis of phenological timing records
################################################################################

##calculating numbers reported in text

split_recs$anyeffnum <- as.numeric(as.factor(split_recs$anyeffect))-1
#by microbe tax
sum(   tapply(split_recs$anyeffnum[split_recs$micrtax=="MF"],split_recs$WOS.relevance.rank[split_recs$micrtax=="MF"],sum) ==0)#4
length(unique(split_recs$WOS.relevance.rank[split_recs$micrtax=="MF"] ))#20
sum(   tapply(split_recs$anyeffnum[split_recs$micrtax=="bacteria"],split_recs$WOS.relevance.rank[split_recs$micrtax=="bacteria"],sum) ==0)#1
length(unique(split_recs$WOS.relevance.rank[split_recs$micrtax=="bacteria"] )) #25
length(table(split_recs$WOS.relevance.rank[split_recs$micrloc=="multiple"])) #14
length(unique(split_recs$WOS.relevance.rank[split_recs$micrtax=="otherfungi"] ))#18
length(unique(split_recs$WOS.relevance.rank[split_recs$micrtax=="mixed"] ))#16
#by phenophase
sum(   tapply(split_recs$anyeffnum[split_recs$phentrait=="flowering time"],split_recs$WOS.relevance.rank[split_recs$phentrait=="flowering time"],sum) ==0)#43
length(unique(split_recs$WOS.relevance.rank[split_recs$phentrait=="flowering time"]))#6
sum(   tapply(split_recs$anyeffnum[split_recs$phentrait=="germination time"],split_recs$WOS.relevance.rank[split_recs$phentrait=="germination time"],sum) ==0)#2
length(unique(split_recs$WOS.relevance.rank[split_recs$phentrait=="germination time"]))#19
sum(   tapply(split_recs$anyeffnum[split_recs$phentrait=="germination prb"],split_recs$WOS.relevance.rank[split_recs$phentrait=="germination prb"],sum) ==0)#1
length(unique(split_recs$WOS.relevance.rank[split_recs$phentrait=="germination prb"]))#26
sum(   tapply(split_recs$anyeffnum[split_recs$phentrait=="fruiting time"],split_recs$WOS.relevance.rank[split_recs$phentrait=="fruiting time"],sum) ==0)#1
length(unique(split_recs$WOS.relevance.rank[split_recs$phentrait=="fruiting time"]))#7
sum(   tapply(split_recs$anyeffnum[split_recs$phentrait=="senescence time"],split_recs$WOS.relevance.rank[split_recs$phentrait=="senescence time"],sum) ==0)#1
length(unique(split_recs$WOS.relevance.rank[split_recs$phentrait=="senescence time"]))#2
sum(   tapply(split_recs$anyeffnum[split_recs$phentrait%in% c("phyllochron","budburst")],split_recs$WOS.relevance.rank[split_recs$phentrait %in% c("phyllochron","budburst") ],sum) ==0)#1
length(unique(split_recs$WOS.relevance.rank[split_recs$phentrait %in% c("phyllochron","budburst") ])) #9
#
table(split_recs$microbeNPHO)
length(unique(split_recs$WOS.relevance.rank[split_recs$lifeform=="annual" | split_recs$lifeform=="both"]))
length(unique(split_recs$WOS.relevance.rank[split_recs$lifeform=="perennial"]))
#appendix S5, able of plant families
apps5 <- cbind(names(table(split_recs$plant.family)),table(split_recs$plant.family), 
	sapply(sort(unique(split_recs$plant.family)), function(fam) 
		 length(unique(split_recs$WOS.relevance.rank[split_recs$plant.family==fam])) ),
		sapply(sort(unique(split_recs$plant.family)), function(fam) length(which(split_recs$lifeform[split_recs$plant.family==fam]=="annual")) ),
		sapply(sort(unique(split_recs$plant.family)), function(fam) length(which(split_recs$lifeform[split_recs$plant.family==fam]=="perennial")) ),
		sapply(sort(unique(split_recs$plant.family)), function(fam) length(which(split_recs$lifeform[split_recs$plant.family==fam]=="both")) )
	)
	#
colnames(apps5) <- c("Plant Family","Number of Tests","Number of Studies","Tests with annuals","Tests with perennials","Tests with mixed life history species")
write.csv(apps5,"AppendixS5.csv",row.names=F)	

#repeat above alternately for subset data, if exclude records that could not be split: some observed effects, but are not included in model results
split_recs_scr$anyeffnum <- as.numeric(as.factor(split_recs_scr$anyeffect))-1 #



sapply(sort(unique(split_recs_scr$microbeNPHO)), function(z) sapply(sort(unique(split_recs_scr$micrtax)), function(k) sum(split_recs_scr$micrtax==k & split_recs_scr$microbeNPHO==z) ) )
sapply(sort(unique(split_recs_scr$micrloc)), function(z) sapply(sort(unique(split_recs_scr$micrtax)), function(k) sum(split_recs_scr$micrtax==k & split_recs_scr$micrloc==z) ) )
sapply(sort(unique(split_recs_scr$micrloc)), function(z) sapply(sort(unique(split_recs_scr$microbeNPHO)), function(k) sum(split_recs_scr$microbeNPHO==k & split_recs_scr$micrloc==z) ) )


##split early/delay/expanded/narrowed into binary variables
numeric.phen.eff <- rep(NA, length=nrow(split_recs_scr))
numeric.phen.eff[split_recs_scr$direction.effect.split=="earlier"] <- -1
numeric.phen.eff[split_recs_scr$direction.effect.split=="delayed"] <- 1
numeric.phen.eff[split_recs_scr$direction.effect.split=="none"] <- 0
numeric.phen.eff[split_recs_scr$direction.effect.split=="expanded"] <- 1
numeric.phen.eff[split_recs_scr$direction.effect.split=="narrowed"] <- -1
split_recs_scr$numeric.phen.eff <- numeric.phen.eff
split_recs_scr$isEarlyOrNrw <- ifelse(numeric.phen.eff=="-1",1,0)
split_recs_scr$notEarlyOrNrw <- 1-split_recs_scr$isEarlyOrNrw
split_recs_scr$isDelayOrExp <- ifelse(numeric.phen.eff=="1",1,0)
split_recs_scr$notDelayOrExp <- 1-split_recs_scr$isDelayOrExp
split_recs_scr$isSig <- ifelse(numeric.phen.eff=="1" | numeric.phen.eff=="-1",1,0)
split_recs_scr$notSig <- 1- split_recs_scr$isSig 


#rotate levels for evaluating significant differences among all groups
split_recs_scr$micrtaxmfa <- split_recs_scr$micrtax
split_recs_scr$micrtaxmfa[split_recs_scr$micrtaxmfa=="MF"] <- "0MF"
split_recs_scr$micrtaxmixa <- split_recs_scr$micrtax
split_recs_scr$micrtaxmixa[split_recs_scr$micrtaxmixa=="mixed"] <- "0mix"
#
split_recs_scr$strainvcommsa <- split_recs_scr$inoctype
split_recs_scr$strainvcommsa[split_recs_scr$strainvcommsa=="strain"] <- "0strain"
#
split_recs_scr$lifeformpa <- split_recs_scr$lifeform
split_recs_scr$lifeformpa[split_recs_scr$lifeformpa == "perennial"] <- "0perennial"
#
split_recs_scr$mphpa <- split_recs_scr$microbeNPHO
split_recs_scr$mphpa[split_recs_scr$mphpa == "pathogen"] <- "0path"
split_recs_scr$mphba <- split_recs_scr$microbeNPHO
split_recs_scr$mphba[split_recs_scr$mphba == "other beneficial"] <- "0ob"
#
split_recs_scr$micrlocseeda <- split_recs_scr$micrloc
split_recs_scr$micrlocseeda[split_recs_scr$micrlocseeda=="seed"] <- "0seed"
split_recs_scr$micrlocshoota <- split_recs_scr$micrloc
split_recs_scr$micrlocshoota[split_recs_scr$micrlocshoota=="shoot"] <- "0shoot"
#
split_recs_scr$phentraitfba <- split_recs_scr$phentrait
split_recs_scr$phentraitfba[split_recs_scr$phentrait=="floral budset time"] <- "0fbt"
split_recs_scr$phentraitfla <- split_recs_scr$phentrait
split_recs_scr$phentraitfla[split_recs_scr$phentrait=="flowering time"] <- "0flwr"
split_recs_scr$phentraitfra <- split_recs_scr$phentrait
split_recs_scr$phentraitfra[split_recs_scr$phentrait=="fruiting time"] <- "0frt"
split_recs_scr$phentraitga <- split_recs_scr$phentrait
split_recs_scr$phentraitga[split_recs_scr$phentrait=="germination time"] <- "0gt"
split_recs_scr$phentraitfda <- split_recs_scr$phentrait
split_recs_scr$phentraitfda[split_recs_scr$phentrait=="flowering duration"] <- "0fd"

#function to pool studies with few tests into one category; otherwise random effects cannot be estimated, and random effects are important to weight impacts of studies with many tests
wosbin <- function(dataframe, woscolname,nrecs){
	dataframe$dummyWOSrelrank <- dataframe[,woscolname]
	dataframe$dummyWOSrelrank[dataframe$dummyWOSrelran%in%names(table(dataframe$WOS.relevance.rank)[table(dataframe$WOS.relevance.rank)<nrecs])] <- "binned"
	return(dataframe)
}

	##using above function to bin large studies together
	##and removing categories with too few observations to test differences 
	   ###*** NOTE THAT SOMETIMES THERE IS STILL ONLY ONE STUDY PER CATEGORY *** 
flowering <- wosbin(split_recs_scr[split_recs_scr$phentrait=="flowering time",],"WOS.relevance.rank",15)
	#low sample, removing: virus = 9
	flowerTAX    <- wosbin(split_recs_scr[split_recs_scr$phentrait=="flowering time" & split_recs_scr$micrtax !="virus",],"WOS.relevance.rank",15)
	# plants that both self and outcross are removed from mating system analysis, here and for germtime and germprb objects below
	floweringMAT <- wosbin(split_recs_scr[split_recs_scr$phentrait=="flowering time" & split_recs_scr$mating != "unk" & split_recs_scr$mating != "s+o",],"WOS.relevance.rank",15)
	#low sample, removing: phytoH = 6
	floweringMPH <- wosbin(split_recs_scr[split_recs_scr$phentrait=="flowering time" & !split_recs_scr$microbeNPHO%in%c( "unknown","phytohormones"),],"WOS.relevance.rank",15)
	#there are 0 records for seed microbe impacts on in flowering time;
germtime <- wosbin(split_recs_scr[split_recs_scr$phentrait=="germination time",],"WOS.relevance.rank",15)
	#0 records of MF, mycorrizal fungi, on germination time
	#low sample, removing: shoot = 4, 
	germtimeLOC <- wosbin(split_recs_scr[split_recs_scr$phentrait=="germination time" & split_recs_scr$micrloc != "shoot",],"WOS.relevance.rank",15)
	germtimeMAT <- wosbin(split_recs_scr[split_recs_scr$phentrait=="germination time" & split_recs_scr$mating != "unk" & split_recs_scr$mating != "s+o",],"WOS.relevance.rank",15)
	#low sample, removing: phytohormones = 3 pathothens = 5
	germtimeMPH <- wosbin(split_recs_scr[split_recs_scr$phentrait=="germination time" & !split_recs_scr$microbeNPHO %in%c("unknown","pathogen","phytohormones"),],"WOS.relevance.rank",15)
germprb <- wosbin(split_recs_scr[split_recs_scr$phentrait=="germination prb",],"WOS.relevance.rank",15)
	#low sample, removing: MF, size = 2
	germprbTAX <- wosbin(split_recs_scr[split_recs_scr$phentrait=="germination prb" & split_recs_scr$micrtax !="MF",],"WOS.relevance.rank",15)
	germprbMAT <- wosbin(split_recs_scr[split_recs_scr$phentrait=="germination prb" & split_recs_scr$mating != "unk" & split_recs_scr$mating != "s+o",],"WOS.relevance.rank",15)
	#low sample; removing: pathogen = 9; phytohormones is 0
	germprbMPH <- wosbin(split_recs_scr[split_recs_scr$phentrait=="germination prb" & !split_recs_scr$microbeNPHO%in% c("unknown","pathogen"),],"WOS.relevance.rank",15)
	#low sample, removing: shoot = 4
	germprbLOC <- wosbin(split_recs_scr[split_recs_scr$phentrait=="germination prb" & split_recs_scr$micrloc != "shoot",],"WOS.relevance.rank",15)
#low sample traits, removing, maturation time=4, senescense time = 2, peak flower = 9, flower senescense time = 4; floral bud duration = 9
allphen <- wosbin(split_recs_scr[!is.na(split_recs_scr$phentrait) & (split_recs_scr$phentrait%in% c("budburst","floral budset time","flowering time","phyllochron","fruiting time","germination time")),],"WOS.relevance.rank",15)

###more numbers reported in text:
flowering$anyeffnum <- as.numeric(as.factor(flowering$anyeffect))-1
sum(tapply(flowering$anyeffnum,flowering$WOS.relevance.rank,sum) == 0) #6
length(tapply(flowering$anyeffnum,flowering$WOS.relevance.rank,sum))#40
germtime$anyeffnum <- as.numeric(as.factor(germtime$anyeffect))-1
sum(tapply(germtime$anyeffnum,germtime$WOS.relevance.rank,sum) ==0)#2 
length(tapply(germtime$anyeffnum,germtime$WOS.relevance.rank,sum))#19
germprb$anyeffnum <- as.numeric(as.factor(germprb$anyeffect))-1 
sum(tapply(germprb$anyeffnum,germprb$WOS.relevance.rank,sum) ==0)#1
length(tapply(germprb$anyeffnum,germprb$WOS.relevance.rank,sum))#26-1=25
length(unique(c(germprb$WOS.relevance.rank, germtime$WOS.relevance.rank)))#overlap 


##PRIOR
#NOTE R MUST be fixed. J Hadfield states, "The residual variance is not  identifiable in the likelihood for binary data:I have tried to  explain this intuitively in Section 2.6 of the CourseNotes using a  tasteless example of hospital deaths."
priornuup=list(R=list(V= 1, fix=1), G=list(G1=list(V=1, nu=8)))
	#the prior is somewhat strong for the random effects to solve sampling issues when studies have few tests, or tests only in one category 
	#this prior keeps the model from trying too hard to fit the random effect when there isn't enough information in the data 

######BINOMIALS
 traitdiffs <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~phentrait,random = ~dummyWOSrelrank,data=allphen,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	 traitdiffs2 <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~phentraitfba,random = ~dummyWOSrelrank,data=allphen,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	 traitdiffs3 <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~phentraitfla,random = ~dummyWOSrelrank,data=allphen,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	 traitdiffs4 <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~phentraitfra,random = ~dummyWOSrelrank,data=allphen,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	 traitdiffs5 <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~phentraitga,random = ~dummyWOSrelrank,data=allphen,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
 #phyllochron&fruit(.9) > budburst; fruiting & phyllo > floral budburst; germ < flower < phyll&fruit(.9);   floral budset germination  & flwr(.9) < fruiting; flwr fruit phyll and flwr(.9) > germ
 traitDdiffs <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~phentrait,random = ~dummyWOSrelrank,data=allphen,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	 traitDdiffs2 <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~phentraitfba,random = ~dummyWOSrelrank,data=allphen,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=100000,thin=50,burnin=1000)
	 traitDdiffs3 <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~phentraitfla,random = ~dummyWOSrelrank,data=allphen,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	 traitDdiffs4 <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~phentraitfra,random = ~dummyWOSrelrank,data=allphen,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	 traitDdiffs5 <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~phentraitga,random = ~dummyWOSrelrank,data=allphen,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
#germ less than  budburst and less than flowering, and less than floral budset
	# germination probability and flowering duration are left out of this analysis and not compared to each other, though they have enough records,
		#this is because neither trait can be forced into the categories of earlier or delay; and expansion of probability for a life-history transistion is not the same as duration of a stage within a season

flowerbitax <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~micrtax,random = ~dummyWOSrelrank,data=flowerTAX,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	flowerbitax2 <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~micrtaxmfa,random = ~dummyWOSrelrank,data=flowerTAX,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	flowerbitax3 <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~micrtaxmixa,random = ~dummyWOSrelrank,data=flowerTAX,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
flowerbiscom <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~inoctype,random = ~dummyWOSrelrank,data=flowering,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	flowerbiscom2 <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~strainvcommsa,random = ~dummyWOSrelrank,data=flowering,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
flowerbimat <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~mating,random = ~dummyWOSrelrank,data=floweringMAT,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
flowerbilh <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~lifeform,random = ~dummyWOSrelrank,data=flowering,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	flowerbilh2 <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~lifeformpa,random = ~dummyWOSrelrank,data=flowering,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
flowerbimph <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~microbeNPHO,random = ~dummyWOSrelrank,data=floweringMPH,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	flowerbimph2 <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~mphba,random = ~dummyWOSrelrank,data=floweringMPH,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
flowerbiloc <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~micrloc,random = ~dummyWOSrelrank,data=flowering,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	flowerbiloc2 <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~micrlocshoota,random = ~dummyWOSrelrank,data=flowering,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
#
flowerDbitax <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~micrtax,random = ~dummyWOSrelrank,data=flowerTAX,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	flowerDbitax2 <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~micrtaxmfa,random = ~dummyWOSrelrank,data=flowerTAX,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	flowerDbitax3 <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~micrtaxmixa,random = ~dummyWOSrelrank,data=flowerTAX,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
flowerDbiscom <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~inoctype,random = ~dummyWOSrelrank,data=flowering,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	flowerDbiscom2 <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~strainvcommsa,random = ~dummyWOSrelrank,data=flowering,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
flowerDbimat <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~mating,random = ~dummyWOSrelrank,data=floweringMAT,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
flowerDbilh <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~lifeform,random = ~dummyWOSrelrank,data=flowering,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	flowerDbilh2 <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~lifeformpa,random = ~dummyWOSrelrank,data=flowering,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
flowerDbimph <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~microbeNPHO,random = ~dummyWOSrelrank,data=floweringMPH, family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	flowerDbimph2 <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~mphba,random = ~dummyWOSrelrank,data=floweringMPH, family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
flowerDbiloc <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~micrloc,random = ~dummyWOSrelrank,data=flowering,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	flowerDbiloc2 <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~micrlocshoota,random = ~dummyWOSrelrank,data=flowering,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)

germtimebitax <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~micrtax,random = ~dummyWOSrelrank,data=germtime,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	germtimebitax2 <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~micrtaxmixa,random = ~dummyWOSrelrank,data=germtime,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
germtimebiscom <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~inoctype,random = ~dummyWOSrelrank,data=germtime,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	germtimebiscom2 <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~strainvcommsa,random = ~dummyWOSrelrank,data=germtime,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
germtimebimat <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~mating,random = ~dummyWOSrelrank,data=germtimeMAT,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
germtimebilh <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~lifeform,random = ~dummyWOSrelrank,data=germtime,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	germtimebilh2 <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~lifeformpa,random = ~dummyWOSrelrank,data=germtime,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
germtimebimph <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~microbeNPHO,random = ~dummyWOSrelrank,data=germtimeMPH,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
germtimebiloc <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~micrloc,random = ~dummyWOSrelrank,data=germtimeLOC,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	germtimebiloc2 <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~micrlocseeda,random = ~dummyWOSrelrank,data=germtimeLOC,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
#
germtimeDbitax <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~micrtax,random = ~dummyWOSrelrank,data=germtime,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	germtimeDbitax2 <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~micrtaxmixa,random = ~dummyWOSrelrank,data=germtime,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
germtimeDbiscom <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~inoctype,random = ~dummyWOSrelrank,data=germtime,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	germtimeDbiscom2 <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~strainvcommsa,random = ~dummyWOSrelrank,data=germtime,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
germtimeDbimat <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~mating,random = ~dummyWOSrelrank,data=germtimeMAT,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
germtimeDbilh <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~lifeform,random = ~dummyWOSrelrank,data=germtime,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	germtimeDbilh2 <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~lifeformpa,random = ~dummyWOSrelrank,data=germtime,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
germtimeDbimph <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~microbeNPHO,random = ~dummyWOSrelrank,data=germtimeMPH,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
germtimeDbiloc <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~micrloc,random = ~dummyWOSrelrank,data=germtimeLOC,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	germtimeDbiloc2 <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~micrlocseeda,random = ~dummyWOSrelrank,data=germtimeLOC,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)

germprbbitax <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~micrtax,random = ~dummyWOSrelrank,data=germprbTAX,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	germprbbitax2 <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~micrtaxmixa,random = ~dummyWOSrelrank,data=germprbTAX,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
germprbbiscom <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~inoctype,random = ~dummyWOSrelrank,data=germprb,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	germprbbiscom2 <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~strainvcommsa,random = ~dummyWOSrelrank,data=germprb,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
germprbbimat <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~mating,random = ~dummyWOSrelrank,data=germprbMAT,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
germprbbilh <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~lifeform,random = ~dummyWOSrelrank,data=germprb,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	germprbbilh2 <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~lifeformpa,random = ~dummyWOSrelrank,data=germprb,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
germprbbimph <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~microbeNPHO,random = ~dummyWOSrelrank,data=germprbMPH,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
germprbbiloc <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~micrloc,random = ~dummyWOSrelrank,data=germprbLOC,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	germprbbiloc2 <- MCMCglmm(cbind(isEarlyOrNrw, notEarlyOrNrw)~micrlocseeda,random = ~dummyWOSrelrank,data=germprbLOC,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
#
germprbDbitax <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~micrtax,random = ~dummyWOSrelrank,data=germprbTAX,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	germprbDbitax2 <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~micrtaxmixa,random = ~dummyWOSrelrank,data=germprbTAX,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
germprbDbiscom <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~inoctype,random = ~dummyWOSrelrank,data=germprb,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	germprbDbiscom2 <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~strainvcommsa,random = ~dummyWOSrelrank,data=germprb,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
germprbDbimat <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~mating,random = ~dummyWOSrelrank,data=germprbMAT,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
germprbDbilh <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~lifeform,random = ~dummyWOSrelrank,data=germprb,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	germprbDbilh2 <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~lifeformpa,random = ~dummyWOSrelrank,data=germprb,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
germprbDbimph <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~microbeNPHO,random = ~dummyWOSrelrank,data=germprbMPH,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
germprbDbiloc <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~micrloc,random = ~dummyWOSrelrank,data=germprbLOC,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)
	germprbDbiloc2 <- MCMCglmm(cbind(isDelayOrExp, notDelayOrExp)~micrlocseeda,random = ~dummyWOSrelrank,data=germprbLOC,family="multinomial2",prior=priornuup,verbose=F,pr=T,nitt=1000000,thin=50,burnin=1000)

#function to calculate weighted means
weightmns2 <- function(dataframelist, cols,dfwithnames){
	typelist <- lapply(1:length(cols), function(x) sort(unique( dfwithnames[ ,cols[x] ] )))
	wos <- lapply(dataframelist, function(dat) sort(unique(dat$dummyWOSrelrank) ) )
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

cols1 <- c("microbeNPHO","micrtax","micrloc") 
flwrdlist <- list(flowering,flowering,flowering)
germtlist <- list(germtime,germtime,germtime)
germplist <- list(germprb,germprb,germprb)

#calculate weighted meins with function
flwrE.m <-   weightmns2(flwrdlist,cols=cols1,flowering)[[1]][-c(5,10)] 
flwrE.se <-  weightmns2(flwrdlist,cols=cols1,flowering)[[2]][-c(5,10)]
flwrD.m <-   weightmns2(flwrdlist,cols=cols1,flowering)[[3]][-c(5,10)]
flwrD.se <-  weightmns2(flwrdlist,cols=cols1,flowering)[[4]][-c(5,10)]
germTE.m <-  weightmns2(germtlist,cols=cols1,germtime)[[1]][-5] 
germTE.se <- weightmns2(germtlist,cols=cols1,germtime)[[2]][-5] 
germTD.m <-  weightmns2(germtlist,cols=cols1,germtime)[[3]][-5] 
germTD.se <- weightmns2(germtlist,cols=cols1,germtime)[[4]][-5] 
germPE.m <-  weightmns2(germplist,cols=cols1,germprb)[[1]] 
germPE.se <- weightmns2(germplist,cols=cols1,germprb)[[2]] 
germPD.m <-  weightmns2(germplist,cols=cols1,germprb)[[3]] 
germPD.se <- weightmns2(germplist,cols=cols1,germprb)[[4]] 
#subtracting unknown category in MPH (and virus for flowering, which is only 1 sudy with three tests), which is not useful. leaving in germP to fill gap from phytohormones, of which there are no studies in germination probability


#list of significant differences, filled in from model results above
siggroupsFE <- c("a","ab","b", "-", # #mph , at 0.9 "a","ab","b", "-" ; .95 " "," "," ", "-"
				"ab","a","b","b",#tax #at 0.9; at .9 "ab","a","b","b"; .95 "ab","a","b","ab"
				" "," ","-"," ")#, #loc
siggroupsFD <- c("b","ab","a","-", #"b", #mph # at 0.9 "b","ab","a","-" ; .95 " "," "," ","-",
				"b","b","a","ab",#tax # same for .95 and .9
				"ab","b","-","a")#, #loc #at 0.95 or .9 same
siggroupsGtE <- c("b","a","-","-", #mph; n.s. at .9  "b","a","-","-", at .95 " "," ","-","-"
				"a","-","b","a",#tax  #same at .95
				"","","","-")#, #loc	
siggroupsGtD <- c(" "," ","-","-",#, #mph
				"b","-","a","b",#tax # same at .95
				"ab","b","a","-")#,#loc  ; at .9 "ab","b","a","-"; at .95 " "," "," ","-"
siggroupsGpE <- c(" "," ","-","-", #mph ****no phytohormone records***
				"b","-","a","a",#tax; same at .95
				" "," "," ","-")#, #loc
siggroupsGpD <- c(" "," ","-","-", #mph ****no phytohormone records***
				"","-","","",#tax
				"a","a","b","-")#, #loc; at .9: a a b; at. 95 "ab","a","b","-"

#plot main text figure
pdf("means_ses_and prelim fitted diffs from binom slim 90hpdi.pdf",height=6,width=4)
xvals <- c(seq(from=0,to=1,length.out=c(4+4+4) ))
xlab <- c("nutrients","other beneficial","pathogen","phytohormones","bacteria","mycorrhizal fungi","mixed","other fungi","multiple","root","seed","shoot")
par(mfrow=c(3,1))
par(mar=c(1,0,1,0))
par(oma=c(10,4,1,1))
plot(c(0,1)~I(c(0,1)),pch=NA,xlab="",ylab="",bty="n",xaxt="n",ylim=c(-0.3,1.2),yaxt="n")
	points(unlist(flwrE.m)~xvals[-11],pch=16)
	points(unlist(flwrD.m)~I(xvals[-11]+0.01),pch=16,col="gray")
	arrows(x0=xvals[-11] , y0=flwrE.m -flwrE.se, y1 = flwrE.m +flwrE.se,length=0     )
	arrows(x0=xvals[-11]+0.01 , y0=flwrD.m -flwrD.se, y1 = flwrD.m +flwrD.se,length=0  ,col="gray"   )
	abline(v= c(sum(xvals[4:5])/2, sum(xvals[8:9])/2 ),lty=3 )
	text(xvals,y=1.2,siggroupsFE)
	text(xvals,y=-0.2,siggroupsFD,col="gray")
	mtext(side=3,"flowering time")
	axis(side=2,at=c(0,0.25,0.5,0.75,1),labels=c(0," ",0.5, " ",1))
legend(xvals[9],y=1.2,c("earlier","delay"),fill=c("black","gray"),bty="n")
plot(c(0,1)~I(c(0,1)),pch=NA,xlab="",ylab="",bty="n",xaxt="n",ylim=c(-0.3,1.2),yaxt="n")
	points(unlist(germTE.m)~xvals[-6],pch=16)
	points(unlist(germTD.m)~I(xvals[-6]+0.01),pch=16,col="gray")
	arrows(x0=xvals[-6] , y0=germTE.m -germTE.se, y1 = germTE.m +germTE.se,length=0     )
	arrows(x0=xvals[-6]+0.01 , y0=germTD.m -germTD.se, y1 = germTD.m +germTD.se,length=0  ,col="gray"   )
	abline(v= c(sum(xvals[4:5])/2, sum(xvals[8:9])/2),lty=3 )
	text(xvals,y=1.2,siggroupsGtE)
	text(xvals,y=-0.2,siggroupsGtD,col="gray")
	mtext(side=3,"germination time")
	mtext(side=2,"probability of observing effect",line=2)
	axis(side=2,at=c(0,0.25,0.5,0.75,1),labels=c(0," ",0.5, " ",1))
legend(xvals[9],y=1.2,c("earlier","delay"),fill=c("black","gray"),bty="n")
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
legend(xvals[9],y=1.2,c("narrowed","expanded"),fill=c("black","gray"),bty="n")
axis(side=1,at=xvals, labels=xlab,las=2)
mtext("microbe", side = 1, at = mean(xvals[2:3]), line=8.5) 
mtext("microbe", side = 1, at = mean(xvals[6:7]), line=8.5) 
mtext("microbe", side = 1, at = mean(xvals[10:11]), line=8.5)  
mtext("mechanism", side = 1, at = mean(xvals[2:3]), line=9.5) 
mtext("taxonomy", side = 1, at = mean(xvals[6:7]), line=9.5) 
mtext("location", side = 1, at = mean(xvals[10:11]), line=9.5)  
dev.off()

###################################
###Supplemental figures. (replaces variables for easier code, watch out if running piecewise)
##################################


cols1 <- c("inoctype","lifeform","mating") #"inoctype","lifeform",
flwrdlist <- list(flowering,flowering,floweringMAT)
germtlist <- list(germtime,germtime,germtimeMAT)
germplist <- list(germprb,germprb,germprbMAT)

flwrE.m <-   weightmns2(flwrdlist,cols=cols1,flowering)[[1]][-c(8,10)] 
flwrE.se <-  weightmns2(flwrdlist,cols=cols1,flowering)[[2]][-c(8,10)] 
flwrD.m <-   weightmns2(flwrdlist,cols=cols1,flowering)[[3]][-c(8,10)] 
flwrD.se <-  weightmns2(flwrdlist,cols=cols1,flowering)[[4]][-c(8,10)] 
germTE.m <-  weightmns2(germtlist,cols=cols1,germtime)[[1]] [-c(8,10)] 
germTE.se <- weightmns2(germtlist,cols=cols1,germtime)[[2]] [-c(8,10)] 
germTD.m <-  weightmns2(germtlist,cols=cols1,germtime)[[3]] [-c(8,10)] 
germTD.se <- weightmns2(germtlist,cols=cols1,germtime)[[4]] [-c(8,10)] 
germPE.m <-  weightmns2(germplist,cols=cols1,germprb)[[1]] [-c(8,10)] 
germPE.se <- weightmns2(germplist,cols=cols1,germprb)[[2]] [-c(8,10)] 
germPD.m <-  weightmns2(germplist,cols=cols1,germprb)[[3]] [-c(8,10)] 
germPD.se <- weightmns2(germplist,cols=cols1,germprb)[[4]] [-c(8,10)] 

#fill in from model results
siggroupsFE <- c("b","a","ab",#scom  same at .95/.9
 				"a","ab","b",#lh  same at .95/.9
 				"","") 
siggroupsFD <- c("a","b","b",#scom; same at .95/.9; note looks different from model predictions
				"","","",#lh
				" "," ") #mat
siggroupsGtE <- c("","","",#scom
				" "," "," ",#lh 
				"a","b") #mat n.s. at 0.95 a .9 "a","b"
siggroupsGtD <- c("a","b","c",#scom # at 0.95/.9 same
				"a","ab","b",#lh at .95/.9 same
				"","") #mat
siggroupsGpE <- c("a","a","b",#scom # at .95/9 same
				"","","",#lh
				"","") #mat
siggroupsGpD <- c("b","b","a",#scom ; recall is comm, strain, strain mix; same at .95/.9
				" "," "," ",#lh
				"","") # mat 

pdf("means_ses_and prelim fitted diffs from binom supp cats 90hpdi.pdf",height=6,width=3.5)
xvals <- c(seq(from=0,to=1,length.out=c(3+3+2) ))
xlab <- c("community","strain","strain mix","annual","both","perennial","outcrosser","selfer")
par(mfrow=c(3,1))
par(mar=c(1,0,1,0))
par(oma=c(8,4,1,1))
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
legend(xvals[4],y=1.2,c("earlier","delay"),fill=c("black","gray"),bty="n")
plot(c(0,1)~I(c(0,1)),pch=NA,xlab="",ylab="",bty="n",xaxt="n",ylim=c(-0.3,1.2),yaxt="n")
	points(unlist(germTE.m)~xvals,pch=16)
	points(unlist(germTD.m)~I(xvals+0.01),pch=16,col="gray")
	arrows(x0=xvals , y0=germTE.m -germTE.se, y1 = germTE.m +germTE.se,length=0     )
	arrows(x0=xvals+0.01 , y0=germTD.m -germTD.se, y1 = germTD.m +germTD.se,length=0  ,col="gray"   )
	abline(v= c(sum(xvals[3:4])/2, sum(xvals[6:7])/2 ),lty=3 )
	text(xvals,y=1.2,siggroupsGtE)
	text(xvals,y=-0.2,siggroupsGtD,col="gray")
	mtext(side=3,"germination time")
	mtext(side=2,"probability of observing effect",line=2)
	axis(side=2,at=c(0,0.25,0.5,0.75,1),labels=c(0," ",0.5, " ",1))
legend(xvals[4],y=1.2,c("earlier","delay"),fill=c("black","gray"),bty="n")
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
legend(xvals[4],y=1.2,c("narrowed","expanded"),fill=c("black","gray"),bty="n")
axis(side=1,at=xvals, labels=xlab,las=2)
mtext("inoculation", side = 1, at = mean(xvals[2]), line=6.5) 
mtext("type", side = 1, at = mean(xvals[2]), line=7.5) 
mtext("lifeform", side = 1, at = mean(xvals[5]), line=7) 
mtext("mating", side = 1, at = mean(xvals[7:8]), line=6.5)  
mtext("system", side = 1, at = mean(xvals[7:8]), line=7.5)  
dev.off()


phentraitmnses <-  weightmns2(list(allphen),cols="phentrait",allphen)
siggroupstrait <- c("bc","bc","b","a","c","a",#early 0.9   0.95, fruiting not diff from germ
 #phyllochron&fruit(.9) > budburst; fruiting & phyllo > floral budburst; germ < flower < phyll&fruit(.9);   floral budset germination  & flwr(.9) < fruiting;  fruit phyll and flwr(.9) > germ
 				"a","a","a","ab","b","ab") #delay 0.9 and .95 same
 				 #germ less than  budburst and less than flowering, and less than floral budset

 				
pdf("phenophase fitted diffs from binom 90hpdi.pdf",height=4,width=4)
xvals <- c(seq(from=0,to=1,length.out=6 ))
xlab <- sort(unique(allphen$phentrait))
par(mar=c(1,0,1,0))
par(oma=c(7,4,1,1))
plot(c(0,1)~I(c(0,1)),pch=NA,xlab="",ylab="",bty="n",xaxt="n",ylim=c(-0.3,1.2),yaxt="n")
	points(unlist(phentraitmnses[[1]])~xvals,pch=16)
	points(unlist(phentraitmnses[[3]])~I(xvals+0.01),pch=16,col="gray")
	arrows(x0=xvals , y0=phentraitmnses[[1]] -phentraitmnses[[2]], y1 = phentraitmnses[[1]] + phentraitmnses[[2]],length=0     )
	arrows(x0=xvals+0.01 , y0=phentraitmnses[[3]] -phentraitmnses[[4]], y1 = phentraitmnses[[3]] + phentraitmnses[[4]],length=0  ,col="gray"   )
	text(xvals,y=1.2,siggroupstrait[1:6])
	text(xvals,y=-0.2,siggroupstrait[7:12],col="gray")
	mtext(side=2,"probability of observing effect",line=2)
	axis(side=2,at=c(0,0.25,0.5,0.75,1),labels=c(0," ",0.5, " ",1))
legend(xvals[4],y=1.2,c("earlier","delay"),fill=c("black","gray"),bty="n")
axis(side=1,at=xvals, labels=xlab,las=2)
dev.off()


###main figure, re-done for 95% HPDI
cols1 <- c("microbeNPHO","micrtax","micrloc") 
flwrdlist <- list(flowering,flowering,flowering)
germtlist <- list(germtime,germtime,germtime)
germplist <- list(germprb,germprb,germprb)

flwrE.m <-   weightmns2(flwrdlist,cols=cols1,flowering)[[1]][-c(5,10)] 
flwrE.se <-  weightmns2(flwrdlist,cols=cols1,flowering)[[2]][-c(5,10)]
flwrD.m <-   weightmns2(flwrdlist,cols=cols1,flowering)[[3]][-c(5,10)]
flwrD.se <-  weightmns2(flwrdlist,cols=cols1,flowering)[[4]][-c(5,10)]
germTE.m <-  weightmns2(germtlist,cols=cols1,germtime)[[1]][-5] 
germTE.se <- weightmns2(germtlist,cols=cols1,germtime)[[2]][-5] 
germTD.m <-  weightmns2(germtlist,cols=cols1,germtime)[[3]][-5] 
germTD.se <- weightmns2(germtlist,cols=cols1,germtime)[[4]][-5] 
germPE.m <-  weightmns2(germplist,cols=cols1,germprb)[[1]] 
germPE.se <- weightmns2(germplist,cols=cols1,germprb)[[2]] 
germPD.m <-  weightmns2(germplist,cols=cols1,germprb)[[3]] 
germPD.se <- weightmns2(germplist,cols=cols1,germprb)[[4]] 
#subtracting unknown category in MPH (and virus for flowering, which is only 1 sudy with three tests), which is not useful. leaving in germP to fill gap from phytohormones, of which there are no studies in germination probability

siggroupsFE <- c(" "," "," ", "-", # #mph , at 0.9 "a","ab","b", "-" ; .95 " "," "," ", "-"
				"ab","a","b","ab",#tax #at 0.9; at .9 "ab","a","b","b"; .95 "ab","a","b","ab"
				" "," ","-"," ")#, #loc
siggroupsFD <- c(" "," "," ","-", #"b", #mph # at 0.9 "b","ab","a","-" ; .95 " "," "," ","-",
				"b","b","a","ab",#tax # same for .95 and .9
				"ab","b","-","a")#, #loc #at 0.95 or .9 same
siggroupsGtE <- c(" "," ","-","-", #mph; n.s. at .9  "b","a","-","-", at .95 " "," ","-","-"
				"a","-","b","a",#tax  #same at .9/.95
				"","","","-")#, #loc	
siggroupsGtD <- c(" "," ","-","-",#, #mph
				"b","-","a","b",#tax # same at .9/.95
				" "," "," ","-")#,#loc  ; at .9 "ab","b","a","-"; at .95 " "," "," ","-"
siggroupsGpE <- c(" "," ","-","-", #mph ****no phytohormone records***
				"b","-","a","a",#tax; same at .95
				" "," "," ","-")#, #loc
siggroupsGpD <- c(" "," ","-","-", #mph ****no phytohormone records***
				"","-","","",#tax
				"ab","a","b","-")#, #loc; at .9: a a b; at. 95 "ab","a","b","-"
pdf("means_ses_and prelim fitted diffs from binom slim 95hpdi.pdf",height=6,width=4)
xvals <- c(seq(from=0,to=1,length.out=c(4+4+4) ))
xlab <- c("nutrients","other beneficial","pathogen","phytohormones","bacteria","mycorrhizal fungi","mixed","other fungi","multiple","root","seed","shoot")
par(mfrow=c(3,1))
par(mar=c(1,0,1,0))
par(oma=c(10,4,1,1))
plot(c(0,1)~I(c(0,1)),pch=NA,xlab="",ylab="",bty="n",xaxt="n",ylim=c(-0.3,1.2),yaxt="n")
	points(unlist(flwrE.m)~xvals[-11],pch=16)
	points(unlist(flwrD.m)~I(xvals[-11]+0.01),pch=16,col="gray")
	arrows(x0=xvals[-11] , y0=flwrE.m -flwrE.se, y1 = flwrE.m +flwrE.se,length=0     )
	arrows(x0=xvals[-11]+0.01 , y0=flwrD.m -flwrD.se, y1 = flwrD.m +flwrD.se,length=0  ,col="gray"   )
	abline(v= c(sum(xvals[4:5])/2, sum(xvals[8:9])/2 ),lty=3 )
	text(xvals,y=1.2,siggroupsFE)
	text(xvals,y=-0.2,siggroupsFD,col="gray")
	mtext(side=3,"flowering time")
	axis(side=2,at=c(0,0.25,0.5,0.75,1),labels=c(0," ",0.5, " ",1))
legend(xvals[9],y=1.2,c("earlier","delay"),fill=c("black","gray"),bty="n")
plot(c(0,1)~I(c(0,1)),pch=NA,xlab="",ylab="",bty="n",xaxt="n",ylim=c(-0.3,1.2),yaxt="n")
	points(unlist(germTE.m)~xvals[-6],pch=16)
	points(unlist(germTD.m)~I(xvals[-6]+0.01),pch=16,col="gray")
	arrows(x0=xvals[-6] , y0=germTE.m -germTE.se, y1 = germTE.m +germTE.se,length=0     )
	arrows(x0=xvals[-6]+0.01 , y0=germTD.m -germTD.se, y1 = germTD.m +germTD.se,length=0  ,col="gray"   )
	abline(v= c(sum(xvals[4:5])/2, sum(xvals[8:9])/2),lty=3 )
	text(xvals,y=1.2,siggroupsGtE)
	text(xvals,y=-0.2,siggroupsGtD,col="gray")
	mtext(side=3,"germination time")
	mtext(side=2,"probability of observing effect",line=2)
	axis(side=2,at=c(0,0.25,0.5,0.75,1),labels=c(0," ",0.5, " ",1))
legend(xvals[9],y=1.2,c("earlier","delay"),fill=c("black","gray"),bty="n")
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
legend(xvals[9],y=1.2,c("narrowed","expanded"),fill=c("black","gray"),bty="n")
axis(side=1,at=xvals, labels=xlab,las=2)
mtext("microbe", side = 1, at = mean(xvals[2:3]), line=8.5) 
mtext("microbe", side = 1, at = mean(xvals[6:7]), line=8.5) 
mtext("microbe", side = 1, at = mean(xvals[10:11]), line=8.5)  
mtext("mechanism", side = 1, at = mean(xvals[2:3]), line=9.5) 
mtext("taxonomy", side = 1, at = mean(xvals[6:7]), line=9.5) 
mtext("location", side = 1, at = mean(xvals[10:11]), line=9.5)  
dev.off()