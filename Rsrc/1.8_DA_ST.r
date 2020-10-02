# Run settings 
library(devtools)
source_url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/Rsrc/settings.r")
if(file.exists("localSettings.r")) {source("localSettings.r")} # use settings file from local directory if one exists

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
  uniqueData_file <- load(paste0("procData/init",startingYear,"/DA",year2,"_split/uniqueData", split_id, ".rdata"))
  uniqueData <- get(uniqueData_file)
  rm(uniqueData_file)
} else{
  load(paste0("procData/init",startingYear,"/DA",year2,"/uniqueData.rdata"))  
}

####load error models
load(url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/data/inputUncer.rdata"))
load(url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/data/logisticPureF.rdata"))
load(url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/data/step.probit.rdata"))
###load surrMods  ###change name
load("surErrMods/surMod.rdata")


uniqueData[,BAp:= (ba * pineP/(pineP+spruceP+blp))]
uniqueData[,BAsp:= (ba * spruceP/(pineP+spruceP+blp))]
uniqueData[,BAb:= (ba * blp/(pineP+spruceP+blp))]

dataSurMod <- uniqueData[,.(segID,h,dbh,BAp,BAsp,BAb,siteType1,
                            siteType2,v2,ba2,h2,dbh2)] 
setnames(dataSurMod,c("segID","H","D","BAp","BAsp","BAb","st1",
                      "st2","V2","ba2","h2","dbh2"))


dataSurMod[,BApPer:=.(BAp/sum(BAp,BAsp,BAb)*100),by=segID]
dataSurMod[,BAspPer:=.(BAsp/sum(BAp,BAsp,BAb)*100),by=segID]
dataSurMod[,BAbPer:=.(BAb/sum(BAp,BAsp,BAb)*100),by=segID]
dataSurMod[,BAtot:=.(sum(BAp,BAsp,BAb)),by=segID]


nSeg <- nrow(dataSurMod)  ##200
stProbMod <- matrix(NA,nSeg,6)
colnames(stProbMod) <- c("segID",paste0("pST",1:5))

if(parallelRun){  ### PARALLEL run
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
  setDT(stProbMod)[, c("segID",paste0("pST", 1:5)) := tstrsplit(V1, ", ")]
  stProbMod$segID <- gsub("^.{0,2}", "", stProbMod$segID)
  stProbMod$pST5 <- gsub(")", "", stProbMod$pST5)
  stProbMod <- stProbMod[, segID:=as.numeric(segID)]
  stProbMod <- stProbMod[, pST1:=as.numeric(pST1)]
  stProbMod <- stProbMod[, pST2:=as.numeric(pST2)]
  stProbMod <- stProbMod[, pST3:=as.numeric(pST3)]
  stProbMod <- stProbMod[, pST4:=as.numeric(pST4)]
  stProbMod <- stProbMod[, pST5:=as.numeric(pST5)]
  stProbMod <- stProbMod[, V1:=NULL]
  
} else {   ### SERIAL run
  for(i in 1:nSeg){
    stProbMod[i,] <- pSTx(dataSurMod[i],nSample,startingYear,year2,tileX)
    # if (i %% 100 == 0) { print(i) }
  }
  stProbMod <- data.table(stProbMod)
}

if (splitRun) {
  save(stProbMod, file = paste0(procDataPath,"init",startingYear,"/calST_split/stProbMod",split_id,".rdata"))
  save(dataSurMod, file = paste0(procDataPath,"init",startingYear,"/calST_split/dataSurMod",split_id,".rdata"))
  
  # If split_id is not highest of the range, processing is stopped here. Continue by sending remaining parts of split data to processing
  
  if(split_id == max(splitRange)){ 
    
    # When the split part with highest id number is processed, stProbMod and dataSurMod of all of the split parts 
    # are combined before continuing processing
    
    # Create empty data tables to where the combined data is added
    stProbMod_bind <- data.table()
    dataSurMod_bind <- data.table()
    
    # Iterate through stProbMod and dataSurMod data tables of the split parts. 
    # Bind split parts to a single data table 
    for (i in splitRange) {
      stProbMod_fileX <- load(paste0("procData/init",startingYear,"/calST_split/stProbMod",i,".rdata"))
      dataSurMod_fileX <- load(paste0("procData/init",startingYear,"/calST_split/dataSurMod",i,".rdata"))
      stProbMod_split <- get(stProbMod_fileX)
      rm(list = stProbMod_fileX)
      rm(stProbMod_fileX)
      dataSurMod_split <- get(dataSurMod_fileX)
      rm(list = dataSurMod_fileX)
      rm(dataSurMod_fileX)
      
      stProbMod_bind <- rbindlist(list(stProbMod_bind, stProbMod_split))
      dataSurMod_bind <- rbindlist(list(dataSurMod_bind, dataSurMod_split))
      
      rm(stProbMod_split, dataSurMod_split)
    }
    
    stProbMod <- stProbMod_bind
    dataSurMod <- dataSurMod_bind
    
    
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
    stProb[,,3] <- as.matrix(stProbMod[,2:6])
    
    stProb <- apply(stProb, c(1,2), mean)
    stProb <- cbind(dataSurMod$segID,stProb)
    colnames(stProb) <- colnames(stProbMod)
    save(stProb,probit1,probit2,stProbMod,file="stProbMod.rdata")
  }
  
  
} else { ### if splitRun = F
  
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
  stProb[,,3] <- as.matrix(stProbMod[,2:6])
  
  stProb <- apply(stProb, c(1,2), mean)
  stProb <- cbind(dataSurMod$segID,stProb)
  colnames(stProb) <- colnames(stProbMod)
  save(stProb,probit1,probit2,stProbMod,file="stProbMod.rdata")
}


