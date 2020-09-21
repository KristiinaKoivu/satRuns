library(devtools)
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
mkfldr <- paste0("initPrebas/","init",startingYear,"/st",siteTypeX)
if(!dir.exists(file.path(generalPath, mkfldr))) {
  dir.create(file.path(generalPath, mkfldr), recursive = TRUE)
}

load(paste0(procDataPath,"init",startingYear,"/","st",siteTypeX,"/samples.rdata"))  

if(!parallelRun) {
  nSamples <- length(samples)
  sampleIDs <- 1:nSamples
}

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
  
  if(parallelRun) {
    # Create PREBAS input data from sample data. Process data in parallel with mclapply command.
    # Number of cores used for processing can be defined with mc.cores argument. mc.cores=1 disables 
    # parallel processing. Recommended amount of cores to use with this code is 20.
    initPrebFiles <- mclapply(seq_along(samples), function(x) {
      ## Prepare the same initial state for all harvest scenarios that are simulated in a loop below
      
      data.sample <- samples[[x]]
      sampleID <- names(samples)[x]
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
      
      clim = prep.climate.f(dat, data.sample, startingYear, nYears,startYearWeather)
      
      # Region = nfiareas[ID==r_no, Region]
      
      initPrebas = create_prebas_input.f(r_no, clim, data.sample, nYears = nYears,
                                         startingYear = startingYear,domSPrun=domSPrun)
      
      # if(stXruns){
      save(initPrebas,file=paste0(initPrebasPath,"init",startingYear,"/",
                                  "st",siteTypeX,"/",
                                  rcpfile,"_sample",sampleID,".rdata"))
      # }else{
      # save(initPrebas,file=paste0(initPrebasPath,startingYear,"/",
      # rcpfile,"_sample",sampleID,".rdata"))
      # }
      
      # save(initPrebas,file=paste0(initPrebasPath,startingYear,"/",rcpfile,"_sample",sampleID,".rdata"))
      rm(initPrebas); gc()
      #print(sampleID)
    }, mc.cores = coresN)
    
  } else {
      # Create PREBAS input data from sample data.
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
    
    clim = prep.climate.f(dat, data.sample, startingYear, nYears,startYearWeather)
    
    # Region = nfiareas[ID==r_no, Region]
    
    initPrebas = create_prebas_input.f(r_no, clim, data.sample, nYears = nYears,
                                       startingYear = startingYear,domSPrun=domSPrun)
    
    # if(stXruns){
    save(initPrebas,file=paste0(initPrebasPath,"init",startingYear,"/",
                                "st",siteTypeX,"/",
                                rcpfile,"_sample",sampleID,".rdata"))
    # }else{
    # save(initPrebas,file=paste0(initPrebasPath,startingYear,"/",
    # rcpfile,"_sample",sampleID,".rdata"))
    # }
    
    # save(initPrebas,file=paste0(initPrebasPath,startingYear,"/",rcpfile,"_sample",sampleID,".rdata"))
    rm(initPrebas); gc()
    print(sampleID)
    }
  }
}
  
