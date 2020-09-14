# Run settings (if modifiedSettings is not set to TRUE in batch job script, default settings from Github will be used)
source_url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/Rsrc/settings.r")
if(modifiedSettings) {
  source("/scratch/project_2000994/PREBASruns/assessCarbon/Rsrc/mainSettings.r") # in CSC
}

# Run functions 
source_url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/Rsrc/functions.r")


###check and create output directories
setwd(generalPath)

yearX <- 3
nSample = 100 ###number of samples from the error distribution
load(paste0("procData/init",startingYear,"/calST_test/uniqueData.rdata")) 
#load("C:/Users/minunno/GitHub/satRuns/data/inputUncer.rdata")
load("/scratch/project_2000994/PREBASruns/assessCarbon/data/inputUncer.rdata") # in CSC
load("surErrMods/logisticPureF_test.rdata")
load("surErrMods/stProbit_test.rdata")
load("surErrMods/surMod_test.rdata")


uniqueData[,BAp:= (ba * pineP/(pineP+spruceP+blp))]
uniqueData[,BAsp:= (ba * spruceP/(pineP+spruceP+blp))]
uniqueData[,BAb:= (ba * blp/(pineP+spruceP+blp))]

dataSurV <- uniqueData[,.(h,dbh,BAp,BAsp,BAb,siteType1,siteType2,v2,segID)] 
setnames(dataSurV,c("H","D","BAp","BAsp","BAb","st1","st2","V2","segID"))


dataSurV[,BApPer:=.(BAp/sum(BAp,BAsp,BAb)*100),by=segID]
dataSurV[,BAspPer:=.(BAsp/sum(BAp,BAsp,BAb)*100),by=segID]
dataSurV[,BAbPer:=.(BAb/sum(BAp,BAsp,BAb)*100),by=segID]
dataSurV[,BAtot:=.(sum(BAp,BAsp,BAb)),by=segID]



nSeg <- nrow(dataSurV)  ##200
stProbMod <- matrix(NA,nSeg,5)

#Processing time is measured with tictoc
tic("total time taken to run the surrogate model")
for(i in 1:nSeg){
  stProbMod[i,] <- pSTx(dataSurV[i],nSample)
  # if (i %% 100 == 0) { print(i) }
}
toc()
stProbMod <- data.table(stProbMod)

###calculate probit2016
dataSurV[,st:=st1]
probit1 <- predict(step.probit,type='p',dataSurV[1:nSeg,])   ### needs to be changed . We need to calculate with 2016 and 2019 data

###calculate probit2019
dataSurV[,st:=st2]
probit2 <- predict(step.probit,type='p',dataSurV[1:nSeg,])   ### needs to be changed . We need to calculate with 2016 and 2019 data

stProb <- array(NA, dim=c(nSeg,5,3))
stProb[,,1] <- probit1
stProb[,,2] <- probit2
stProb[,,3] <- as.matrix(stProbMod)

stProb <- apply(stProb, c(1,2), mean)

save(stProb,probit1,probit2,stProbMod,file="stProbMod.rdata")

