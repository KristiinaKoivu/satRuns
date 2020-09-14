
# Run settings (if modifiedSettings is not set to TRUE in batch job script, default settings from Github will be used)
source_url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/Rsrc/settings.r")
if(modifiedSettings) {
  source("/scratch/project_2000994/PREBASruns/assessCarbon/Rsrc/mainSettings.r") # in CSC
}

  
if (splitRun) {
  source("/scratch/project_2000994/PREBASruns/assessCarbon/Rsrc/mainSettings_split.r")
  # If output is set to be split to smaller parts (splitRun = TRUE), create separate
  # folder for the split data tables.
  mkfldr_split <- paste0("procData/",paste0("init",startingYear,"/calST_split"))
    
  if(!dir.exists(file.path(generalPath, mkfldr_split))) {
    dir.create(file.path(generalPath, mkfldr_split), recursive = TRUE)
  }
} else {
  source("/scratch/project_2000994/PREBASruns/assessCarbon/Rsrc/mainSettings.r") ### in CSC
  # Create folders for outputs.
  setwd(generalPath)
  mkfldr <- paste0("procData/",paste0("init",startingYear,"/calST"))
  if(!dir.exists(file.path(generalPath, mkfldr))) {
    dir.create(file.path(generalPath, mkfldr), recursive = TRUE)
  }
}


###extract CurrClim IDs
rastX <- raster(baRast)
if(testRun){
  extNew <- extent(rastX)
  extNew[2]   <- (extent(rastX)[1] + (extent(rastX)[2] - extent(rastX)[1])*fracTest)
  extNew[4]   <- (extent(rastX)[3] + (extent(rastX)[4] - extent(rastX)[3])*fracTest)
  rastX <- crop(rastX,extNew)
  maxSitesRun <- maxSitesRunTest
}

climID <- raster(climIDpath)
# climIDx <- crop(climID,rastX)
# plot(climIDx)
# plot(rastX,add=T)
# climIDs <- resample(climIDx,rastX,method="ngb")
# writeRaster(climIDs,paste0(rasterPath,"climIDs.tif"),overwrite=T)
# climIDs <- raster(paste0(rastersPath,"climIDs.tif"))
rm(rastX)
gc()


fileNames <- c(baRast,
               blPerRast,
               dbhRast,
               vRast,
               hRast,
               pinePerRast,
               sprucePerRast,
               siteTypeRast,
               siteTypeRast2,
               vRast2)

for(i in 1:length(fileNames)){
  rastX <- raster(fileNames[i])
  if(testRun) rastX <- crop(rastX,extNew)    ####if it is a test run rasters are cropped to a smaller area
  dataX <- data.table(rasterToPoints(rastX))
  if(i==1){
    data.all <- dataX 
  }else{
    data.all <- merge(data.all,dataX)
  }
  print(fileNames[i])
}

###attach weather ID
data.all$climID <- extract(climID,data.all[,.(x,y)])
# dataX <- data.table(rasterToPoints(climIDs))
# data.all <- merge(data.all,dataX)
setnames(data.all,c("x","y","ba","blp","dbh","v","h","pineP","spruceP","siteType1","siteType2","v2","climID"))

##filter data 
data.all <- data.all[!ba %in% baNA]
data.all <- data.all[!blp %in% blPerNA]
data.all <- data.all[!dbh %in% dbhNA]
data.all <- data.all[!v %in% vNA]
data.all <- data.all[!v2 %in% vNA]
data.all <- data.all[!h %in% hNA]
data.all <- data.all[!pineP %in% pinePerNA]
data.all <- data.all[!spruceP %in% sprucePerNA]
data.all <- data.all[!siteType1 %in% siteTypeNA]
data.all <- data.all[!siteType2 %in% siteTypeNA]

####convert data to prebas units
data.all <- data.all[, ba := ba * baConv]
data.all <- data.all[, blp := blp * blPerConv]
data.all <- data.all[, dbh := dbh * dbhConv]
data.all <- data.all[, v := v * vConv]
data.all <- data.all[, v2 := v2 * vConv]
data.all <- data.all[, h := h * hConv]
data.all <- data.all[, pineP := pineP * pinePerConv]
data.all <- data.all[, spruceP := spruceP * sprucePerConv]
data.all <- data.all[, siteType1 := siteType1 * siteTypeConv]
data.all <- data.all[, siteType2 := siteType2 * siteTypeConv]

if(siteTypeX==year2){
  data.all[,siteType:=siteType2]  
}else if(siteTypeX==startingYear){
  data.all[,siteType:=siteType1]  
}else{
  data.all[,siteType:=siteTypeX]  
}
data.all[siteType>5,siteType:=5]
data.all[siteType1>5,siteType1:=5]
data.all[siteType2>5,siteType2:=5]


#####I'm excluding from the runs the areas that have been clearcutted and have ba=0 
# data.all[h==0. & dbh==0 & ba==0,clCut:=1]
data.all[,clCut:=0]
data.all[ba==0,clCut:=1]

###calculate tree density
data.all[clCut==0,N:=ba/(pi*(dbh/200)^2)]

####check where H is below minimum initial height and replace
smallH <- intersect(which(data.all$h < initH), which(data.all$clCut==0))
data.all[smallH, h:=initH]

###check where density is too high and replase stand variables with initial conditions
tooDens <- intersect(which(data.all$N> maxDens), which(data.all$clCut==0))
data.all[tooDens,h:=initH]
data.all[tooDens,ba:=initBA]
data.all[tooDens,dbh:=initDBH]
data.all[tooDens,N:=initN]


data.all[pineP == 0 & spruceP == 0 & blp ==0 & siteType ==1, blp:=1  ]
data.all[pineP == 0 & spruceP == 0 & blp ==0 & siteType <= 3 & siteType > 1, spruceP:=1  ]
data.all[pineP == 0 & spruceP == 0 & blp ==0 & siteType >= 4, pineP:=1  ]

###!!!!!!!!!!!!########careful with this part##########!!!!!!!!#########

####calculate dV
data.all[,dV := v2-v]
data.all[,dVy := (v2-v)/(year2 - startingYear)]

####group pixels by same values
data.all[, segID := .GRP, by = .(ba, blp,dbh, h, pineP, spruceP, siteType1,siteType2, climID,dVy,v2)]
# data.all[clCut==1 ,hist(dVy)]


####Count segID pix
data.all[, npix:=.N, segID]

# uniqueData <- data.table()
####find unique initial conditions
uniqueData <- unique(data.all[clCut==0 & dVy >0,.(ba,blp,dbh,h,pineP,spruceP,siteType1,
                                                  siteType2,N,climID,segID,npix,dVy,v2)])
uniqueData[,uniqueKey:=1:nrow(uniqueData)]
setkey(uniqueData, uniqueKey)
# uniqueData[,N:=ba/(pi*(dbh/200)^2)]
range(uniqueData$N)
uniqueData[,area:=npix*resX^2/10000]

###assign ID to similar pixels
XYsegID <- data.all[,.(x,y,segID)]

###!!!!!!!!!!!!########end careful with this part##########!!!!!!!!#########

# nSamples <- ceiling(dim(uniqueData)[1]/20000)
# sampleID <- 1
# 
# for(sampleID in sampleIDs){
#   set.seed(1)
#   samplesX <- split(uniqueData, sample(1:nSample, nrow(uniqueData), replace=T))
#   sampleX <- ops[[sampleID]]
#   sampleX[,area := N*resX^2/10000]
#   # sampleX[,id:=climID]
# }


nSamples <- ceiling(dim(uniqueData)[1]/maxSitesRun)
set.seed(1)
sampleset <- sample(1:nSamples, nrow(uniqueData),  replace=T)
samples <- split(uniqueData, sampleset) 

# adding sampleID, sampleRow (= row within sample) 
uniqueData[,sampleID:=sampleset]
uniqueData[,sampleRow:=1:length(h),by=sampleID]

segID <- numeric(0)
for(i in 1:nSamples){
  sampleX <- samples[[i]]
  segID <- c(segID,sampleX$segID)
}

save(data.all,file=paste0(procDataPath,"init",startingYear,"/calST/allData.rdata"))         ### All data
save(uniqueData,file=paste0(procDataPath,"init",startingYear,"/calST/uniqueData.rdata"))    ### unique pixel combination to run in PREBAS
save(samples,file=paste0(procDataPath,"init",startingYear,"/calST/samples.rdata"))    ### unique pixel combination to run in PREBAS
save(XYsegID,segID,file=paste0(procDataPath,"init",startingYear,"/calST/XYsegID.rdata"))    ### Coordinates and segID of all pixels

#### If needed (splitRun = TRUE), unique data is split to ten tables here to enable 
#    running 1.9_optST in multiple sections. Running in multiple sections will reduce 
#    total processing time of 1.9_optST.
if (splitRun) {
  
  # Create split_id column which is used in splitting the table. Here the data is split to
  # x amount of parts. NOTICE that the parts are not necessarily equal sized: all but the last on
  # are equal sized but the last one is either equal or few observations smaller depending 
  # on the amount of rows in the original data.
  split_length <- ceiling(nrow(uniqueData)/10)
  uniqueData <- uniqueData[, split_id := NA]
  uniqueData$split_id[1:split_length] <- 1
  uniqueData$split_id[(split_length+1):(split_length*2)] <- 2
  uniqueData$split_id[(2*split_length+1):(split_length*3)] <- 3
  uniqueData$split_id[(3*split_length+1):(split_length*4)] <- 4
  uniqueData$split_id[(4*split_length+1):(split_length*5)] <- 5
  uniqueData$split_id[(5*split_length+1):(split_length*6)] <- 6
  uniqueData$split_id[(6*split_length+1):(split_length*7)] <- 7
  uniqueData$split_id[(7*split_length+1):(split_length*8)] <- 8
  uniqueData$split_id[(8*split_length+1):(split_length*9)] <- 9
  uniqueData$split_id[(9*split_length+1):(nrow(uniqueData))] <- 10
  
  # Split the table to list of defined amount of elements. Splitting is done based on the split_id.
  split_list <- split(uniqueData,uniqueData$split_id)
  
  # Convert the split results to separate data tables
  uniqueData1 <- as.data.table(split_list[[1]])
  uniqueData2 <- as.data.table(split_list[[2]])
  uniqueData3 <- as.data.table(split_list[[3]])
  uniqueData4 <- as.data.table(split_list[[4]])
  uniqueData5 <- as.data.table(split_list[[5]])
  uniqueData6 <- as.data.table(split_list[[6]])
  uniqueData7 <- as.data.table(split_list[[7]])
  uniqueData8 <- as.data.table(split_list[[8]])
  uniqueData9 <- as.data.table(split_list[[9]])
  uniqueData10 <- as.data.table(split_list[[10]])
  
  # Remove split_id column
  uniqueData1 <- uniqueData1[, split_id:=NULL]
  uniqueData2 <- uniqueData2[, split_id:=NULL]
  uniqueData3 <- uniqueData3[, split_id:=NULL]
  uniqueData4 <- uniqueData4[, split_id:=NULL]
  uniqueData5 <- uniqueData5[, split_id:=NULL]
  uniqueData6 <- uniqueData6[, split_id:=NULL]
  uniqueData7 <- uniqueData7[, split_id:=NULL]
  uniqueData8 <- uniqueData8[, split_id:=NULL]
  uniqueData9 <- uniqueData9[, split_id:=NULL]
  uniqueData10 <- uniqueData10[, split_id:=NULL]
  
  # Save splitted tables 
  save(uniqueData1,file=paste0(procDataPath,"init",startingYear,"/calST_split/uniqueData1.rdata"))  
  save(uniqueData2,file=paste0(procDataPath,"init",startingYear,"/calST_split/uniqueData2.rdata"))
  save(uniqueData3,file=paste0(procDataPath,"init",startingYear,"/calST_split/uniqueData3.rdata"))
  save(uniqueData4,file=paste0(procDataPath,"init",startingYear,"/calST_split/uniqueData4.rdata"))
  save(uniqueData5,file=paste0(procDataPath,"init",startingYear,"/calST_split/uniqueData5.rdata"))
  save(uniqueData6,file=paste0(procDataPath,"init",startingYear,"/calST_split/uniqueData6.rdata"))
  save(uniqueData7,file=paste0(procDataPath,"init",startingYear,"/calST_split/uniqueData7.rdata"))
  save(uniqueData8,file=paste0(procDataPath,"init",startingYear,"/calST_split/uniqueData8.rdata"))
  save(uniqueData9,file=paste0(procDataPath,"init",startingYear,"/calST_split/uniqueData9.rdata"))
  save(uniqueData10,file=paste0(procDataPath,"init",startingYear,"/calST_split/uniqueData10.rdata"))
  
}
