# Run settings (if modifiedSettings is not set to TRUE in batch job script, default settings from Github will be used)
source_url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/Rsrc/settings.r")
if(exists("modifiedSettings")) {
  if(modifiedSettings) {
    source("/scratch/project_2000994/PREBASruns/assessCarbon/Rsrc/mainSettings.r") # in CSC
  }
}

# Run functions 
source_url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/Rsrc/functions.r")


###check and create output directories
setwd(generalPath)

yearX <- 3
nSample = 1000 ###number of samples from the error distribution
# Load unique data.
# If data is processed in split parts, define to variable split_id which split part to process (in batch job script).
# If splitRun is not needed, the unique data dataset for the whole tile is loaded.
if (splitRun) {
  uniqueData_file <- load(paste0("procData/init",startingYear,"/calST_split/uniqueData", split_id, ".rdata"))
  uniqueData <- get(uniqueData_file)
  rm(uniqueData_file)
} else{
  load(paste0("procData/init",startingYear,"/calST/uniqueData.rdata")) 
}

####load error models
load(url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/data/inputUncer.rdata"))
load(url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/data/logisticPureF.rdata"))
load(url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/data/step.probit.rdata"))
###load surrMods
load("surErrMods/surMod.rdata")


uniqueData[,BAp:= (ba * pineP/(pineP+spruceP+blp))]
uniqueData[,BAsp:= (ba * spruceP/(pineP+spruceP+blp))]
uniqueData[,BAb:= (ba * blp/(pineP+spruceP+blp))]

dataSurMod <- uniqueData[,.(h,dbh,BAp,BAsp,BAb,siteType1,
                            siteType2,v2,ba2,h2,dbh2,segID)] 
setnames(dataSurMod,c("H","D","BAp","BAsp","BAb","st1",
                      "st2","V2","ba2","h2","dbh2","segID"))


dataSurMod[,BApPer:=.(BAp/sum(BAp,BAsp,BAb)*100),by=segID]
dataSurMod[,BAspPer:=.(BAsp/sum(BAp,BAsp,BAb)*100),by=segID]
dataSurMod[,BAbPer:=.(BAb/sum(BAp,BAsp,BAb)*100),by=segID]
dataSurMod[,BAtot:=.(sum(BAp,BAsp,BAb)),by=segID]


nSeg <- nrow(dataSurMod)  ##200
stProbMod <- matrix(NA,nSeg,5)

if(parallelRun){
  # Run surrogate model in parallel. Number of cores used for processing is defined with 
  # mc.cores argument. mc.cores=1 disables parallel processing. For now the result is 
  # a matrix in which all 5 site probability values are in one column. This could be modified 
  # later on to produce a matrix where the 5 values are places in separate columns right away. 
  # That modification would decrease the need of reorganizing the data afterwards (in 
  # lines 73-82).
  stProbMod <- matrix(mclapply(1:nrow(stProbMod), function(i,stProbMod){
    stProbMod[i,] <- pSTx(dataSurMod[i],nSample,startingYear,year2,tileX)
  },stProbMod=stProbMod,mc.cores = coresN))
  
  #Clean and convert the matrix to data table
  stProbMod <- as.data.table(stProbMod)
  setDT(stProbMod)[, paste0("V", 1:5) := tstrsplit(V1, ", ")]
  stProbMod$V1 <- gsub("^.{0,2}", "", stProbMod$V1)
  stProbMod$V5 <- gsub(")", "", stProbMod$V5)
  stProbMod <- stProbMod[, V1:=as.numeric(V1)]
  stProbMod <- stProbMod[, V2:=as.numeric(V2)]
  stProbMod <- stProbMod[, V3:=as.numeric(V3)]
  stProbMod <- stProbMod[, V4:=as.numeric(V4)]
  stProbMod <- stProbMod[, V5:=as.numeric(V5)]
  
} else {
  for(i in 1:nSeg){
    stProbMod[i,] <- pSTx(dataSurMod[i],nSample,startingYear,year2,tileX)
    # if (i %% 100 == 0) { print(i) }
  }
  stProbMod <- data.table(stProbMod)
}

if (splitRun) {
  # If executing a split run, the output is saved at this point.
  # Continue by processing the other parts of split data or by running 1.9_finishSplitRun 
  # if all other parts are processed this far as well.
  save(stProbMod, file = paste0(procDataPath,"init",startingYear,"/calST_split/stProbMod",split_id,".rdata"))
  save(dataSurMod, file = paste0(procDataPath,"init",startingYear,"/calST_split/dataSurMod",split_id,".rdata"))
  print("stPRobMod and dataSurMod of this split_id saved to calST_split folder")
  print("Continue by sending remaining parts of split data to processing OR by combining 
        processed parts and finishing processing of 1.8 with 1.8_finishSplitRun.r")
  

} else {

  ###calculate probit2016
  dataSurMod[,st:=st1]
  step.probit1 <- step.probit[[paste0("y",startingYear)]][[paste0("t",tileX)]]
  probit1 <- predict(step.probit1,type='p',dataSurMod[1:nSeg,])   ### needs to be changed . We need to calculate with 2016 and 2019 data

  ###calculate probit2019
  dataSurMod[,st:=st2]
  step.probit2 <- step.probit[[paste0("y",year2)]][[paste0("t",tileX)]]
  probit2 <- predict(step.probit2,type='p',dataSurMod[1:nSeg,])   ### needs to be changed . We need to calculate with 2016 and 2019 data

  stProb <- array(NA, dim=c(nSeg,5,3))
  stProb[,,1] <- probit1
  stProb[,,2] <- probit2
  stProb[,,3] <- as.matrix(stProbMod)

  stProb <- apply(stProb, c(1,2), mean)

  save(stProb,probit1,probit2,stProbMod,file="stProbMod.rdata")
}


