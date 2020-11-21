library(UpSetR)

# analyze WOS records for AjB lit review on influence of microbes plant phenology
setwd("~/MicrobesAndPhenology/")

wos_recs <- read.csv("~/saved records WOS friday sep 11 2020 - WOS records.tsv",
					sep="\t",stringsAsFactors=F,na.strings = c("NA",""))

cwos_recs <- read.csv("~/saved records WOS friday sep 11 2020 - Copy of WOS records.tsv",
					sep="\t",stringsAsFactors=F,na.strings = c("NA",""))


relevant <- as.numeric(as.factor(cwos_recs$"Relevant..Y.N"))>=2 & !is.na(cwos_recs$"Relevant..Y.N")
PiM <- as.numeric(as.factor(cwos_recs$"Phenology....microbes...if.only.this..stop.."))==3 & !is.na(cwos_recs$"Phenology....microbes...if.only.this..stop..") # plant phenology influences microbes
PoM <- as.numeric(as.factor(cwos_recs$"phenology.of.micriobes...if.only.this..stop.."))==2 & !is.na(cwos_recs$"phenology.of.micriobes...if.only.this..stop..") #phenology of microbes
MiP <- !is.na(cwos_recs$What.phenological.trait.) # microbes influence plant phenology
MaltSonP <- cwos_recs$"Did.they.measure.selection...on.the.plant.only.."=="y" #phenology of microbes


minyear <- min(cwos_recs$Publication.Year,na.rm=T)
wosbyyear <- sapply(minyear:2020, function(yr) sum(cwos_recs$Publication.Year == yr,na.rm=T) )
wosbyyearR <- sapply(minyear:2020, function(yr) sum(cwos_recs[relevant,]$Publication.Year == yr,na.rm=T) )
wosbyyearPiM <- sapply(minyear:2020, function(yr) sum(cwos_recs[PiM,]$Publication.Year == yr,na.rm=T) )
wosbyyearPoM <- sapply(minyear:2020, function(yr) sum(cwos_recs[PoM,]$Publication.Year == yr,na.rm=T) )
# table(cwos_recs[PoM,]$Publication.Year)
wosbyyearMiP <- sapply(minyear:2020, function(yr) sum(cwos_recs[MiP,]$Publication.Year == yr,na.rm=T) )
wosbyyearSEL <- sapply(minyear:2020, function(yr) sum(cwos_recs[MaltSonP,]$Publication.Year == yr,na.rm=T) )
# table(cwos_recs[MiP,]$Publication.Year)

years <- c(minyear:2020)
pdf("~/WOS_AjB_micr_phen.pdf",height=3.5,width=3.5)
par(mar=c(3,3,1,1))
par(oma=c(0,0,0,0))
plot(wosbyyear~years,type="l",lty=2,xlab="",ylab="")
	lines(wosbyyearR~years,col=rgb(0,0,0))
	lines(wosbyyearPiM~years,col=rgb(0.75,0,0))
	lines(wosbyyearPoM~years,col=rgb(0,0,0.75))
	lines(wosbyyearMiP~years,col=rgb(0.5,0,0.5))
	legend(1990,78,c("all records","relevant","phen -> micr","phen of micr","micr -> phen"),
			 lty=c(2,1,1,1,1),col=rgb( c(0,0,0.75,0,0.5),0,c(0,0,0,0.75,0.5) ) ,bty="n"   ) 
	mtext("WOS records",side=2,line=2)
	mtext("Year",side=1,line=2)
dev.off()

gi
 pdf("~/WOS_AjB_micr_phen_r_sim.pdf",height=3.5,width=4.5)
# png("~/WOS_AjB_micr_phen_r_sim.png",res=600)
par(mar=c(3,3,1,1))
par(oma=c(0,0,0,0))
plot(wosbyyearR~years,type="l",lty=2,xlab="",ylab="",ylim=c(0,20))
# 	lines(wosbyyearPiM+wosbyyearPoM~years,col=rgb(0.75,0,0))
	lines(wosbyyearMiP~years)
	lines(wosbyyearSEL~years,col=rgb(0,0,1))
	legend(1986,20,c("plant microbes & phenology","microbes affects phenology?","microbes alter selection on phenology?"),
			 lty=c(2,1,1 ),col=c(rgb(0,0,0),rgb(0,0,0),rgb(0,0,1)) ,bty="n"   ) 

	mtext("WOS records",side=2,line=2)
	mtext("Year",side=1,line=2)
dev.off()


wos_rel <- cwos_recs[relevant,]

# sim2mult <- wos_rel$microbetax.sim2; sim2mult[ grep(";", wos_rel$microbetax.sim2)] <- "multiple treatment types"
sim2and <- gsub("; ", wos_rel$microbetax.sim2,replacement="&")
simand <- gsub("; ", wos_rel$microbetax.sim,replacement="&")
loc2and <- gsub("; ", wos_rel$Microbelocation.2,replacement="&")
phen.and <- gsub("; ", wos_rel$phen.trait,replacement="&")

pdf("~/upsetrtax.pdf",width=4,height=5)
upset(fromExpression(table(sim2and)),nsets=10, order.by = "freq")
dev.off()
pdf("~/upsetrtaxcomp.pdf",width=4,height=5)
upset(fromExpression(table(simand)),nsets=15, order.by = "freq")
dev.off()
pdf("~/upsetrloc.pdf",width=4,height=5)
upset(fromExpression(table(loc2and)),nsets=10, order.by = "freq")
dev.off()
pdf("~/upsetrtrt.pdf",width=4,height=5)
upset(fromExpression(table(phen.and)), nsets=15, order.by = "freq")
dev.off()



# sim2exp <- as.factor(as.character(unlist(sapply(wos_rel$microbetax.sim2, function(x) unlist(strsplit(x,split="; "))))))
# simexp <- as.factor(as.character(unlist(sapply(wos_rel$microbetax.sim, function(x) unlist(strsplit(x,split="; "))))))
trtexp <- as.factor(as.character(unlist(sapply(wos_rel$phen.trait, function(x) unlist(strsplit(x,split="; "))))))


# barplot(table(wos_rel$Strain.mix.comm.simple))
# barplot(table(wos_rel$wos_rel$microbetax.sim2))


