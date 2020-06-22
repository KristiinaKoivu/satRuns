library(reshape2)
library(plyr)
library("raster")
library(data.table)
require(sm)
require(rgdal)



### Run settings & functions
source("Rsrc/settings.r")
source("Rsrc/functions.r")

###check prebas version and install if needed
devtools::install_github("ForModLabUHel/Rprebasso", ref=vPREBAS)
require(Rprebasso)


setwd(generalPath)
if(!dir.exists("initPrebas")) {
  dir.create("initPrebas")
}

load(paste0(procDataPath,startingYear,"/samples.rdata"))
nSamples <- length(samples)
sampleIDs <- 1:nSamples

if(testRun){
  sampleID <- 1
  rcpfile="CurrClim"
}


for (rcpfile in weather) { ## ---------------------------------------------
  print(date())
  print(rcpfile)
  if(rcpfile=="CurrClim"){
    load(paste(climatepath, rcpfile,".rdata", sep=""))  
    setnames(dat,"id","climID")
    #####process data considering only current climate###
    # dat <- dat[rday %in% 1:10958] #uncomment to select some years (10958 needs to be modified)
    # maxRday <- max(dat$rday)
    # xday <- c(dat$rday,(dat$rday+maxRday),(dat$rday+maxRday*2))
    # dat = rbind(dat,dat,dat)
    # dat[,rday:=xday]
    
  } else{
    load(paste(climatepath, rcpfile, sep=""))  
  }

  gc()
  
  
  for(sampleID in sampleIDs){
    ## Prepare the same initial state for all harvest scenarios that are simulated in a loop below
    data.sample <- samples[[sampleID]]
    # nSample <- nrow(sampleX)
    # data.sample = sample_data.f(sampleX, nSample)
    totAreaSample <- sum(data.sample$area)
    
    ###check if climID matches
    allclIDs <- unique(dat$climID)
    samClIds <- unique(data.sample$climID)
    if(!all(samClIds %in% allclIDs)){
      opsClim <- samClIds[which(!samClIds %in% allclIDs)]
      dt = data.table(allclIDs, val = allclIDs) # you'll see why val is needed in a sec
      setnames(dt,c("x","val"))
      # setattr(dt, "sorted", "x")  # let data.table know that w is sorted
      setkey(dt, x) # sorts the data
      # binary search and "roll" to the nearest neighbour
      replX <- dt[J(opsClim), roll = "nearest"]
      data.sample$climID <- mapvalues(data.sample$climID,replX[[1]],replX[[2]])
    }
    
    clim = prep.climate.f(dat, data.sample, startingYear, nYears)
    
    # Region = nfiareas[ID==r_no, Region]
    
    initPrebas = create_prebas_input.f(r_no, clim, data.sample, nYears = nYears,
                                       startingYear = startingYear,domSPrun=domSPrun)
    
    
    save(initPrebas,file=paste0(initPrebasPath,rcpfile,"_sample",sampleID,".rdata"))
    rm(initPrebas); gc()
    print(sampleID)
  }
}
  
