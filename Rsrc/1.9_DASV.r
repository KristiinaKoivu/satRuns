library(MASS)
library(devtools)
# Run settings (if modifiedSettings is not set to TRUE in batch job script, default settings from Github will be used)
source_url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/Rsrc/settings.r")
if(file.exists("localSettings.r")) {source("localSettings.r")} # use settings in local directory if one exists

# Run functions 
source_url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/Rsrc/functions.r")

###check and create output directories
setwd(generalPath)

yearX <- 3
nSample = 1000 ###number of samples from the error distribution
load(paste0("procData/init",startingYear,"/calST/uniqueData.rdata"))  

####load error models
load(url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/data/inputUncer.rdata"))
load(url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/data/logisticPureF.rdata"))
load(url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/data/step.probit.rdata"))
###load surrMods
load("surErrMods/surMod.rdata")

uniqueData[,BAp:= (ba * pineP/(pineP+spruceP+blp))]
uniqueData[,BAsp:= (ba * spruceP/(pineP+spruceP+blp))]
uniqueData[,BAb:= (ba * blp/(pineP+spruceP+blp))]
uniqueData[,BAp2:= (ba2 * pineP2/(pineP2+spruceP2+blp2))]
uniqueData[,BAsp2:= (ba2 * spruceP2/(pineP2+spruceP2+blp2))]
uniqueData[,BAb2:= (ba2 * blp2/(pineP2+spruceP2+blp2))]

dataSurMod <- uniqueData[,.(h,dbh,BAp,BAsp,BAb,siteType1,
                            siteType2,v2,ba2,h2,dbh2,
                            BAp2,BAsp2,BAb2,segID)] 
setnames(dataSurMod,c("H","D","BAp","BAsp","BAb","st1",
                      "st2","V2","ba2","h2","dbh2",
                      "BAp2","BAsp2","BAb2","segID"))


dataSurMod[,BApPer:=.(BAp/sum(BAp,BAsp,BAb)*100),by=segID]
dataSurMod[,BAspPer:=.(BAsp/sum(BAp,BAsp,BAb)*100),by=segID]
dataSurMod[,BAbPer:=.(BAb/sum(BAp,BAsp,BAb)*100),by=segID]
dataSurMod[,BAtot:=.(sum(BAp,BAsp,BAb)),by=segID]
dataSurMod[,BApPer2:=.(BAp2/sum(BAp2,BAsp2,BAb2)*100),by=segID]
dataSurMod[,BAspPer2:=.(BAsp2/sum(BAp2,BAsp2,BAb2)*100),by=segID]
dataSurMod[,BAbPer2:=.(BAb2/sum(BAp2,BAsp2,BAb2)*100),by=segID]
dataSurMod[,BAtot2:=.(sum(BAp2,BAsp2,BAb2)),by=segID]

nSeg <- nrow(dataSurMod)  ##200
load("stProbMod.rdata")
colnames(stProb) <- paste0("pST",1:5)
dataSurMod <- cbind(dataSurMod[1:nSeg],stProb[1:nSeg,])

# pMvNorm <- matrix(NA,42,nSeg)
# system.time({
#   for(i in 1:nSeg){
#     pMvNorm[,i] <- pSVDA(dataSurMod[i],nSample)
#     if (i %% 100 == 0) { print(i) }
#   }
# })

system.time({
 pMvn <- dataSurMod[1:nSeg, pSVDA(.SD,nSample = nSample,year1=startingYear,
                              year2=year2,tileX=tileX), by = seq_len(nSeg)]
})
save(pMvn,file="pMvn.rdata")
