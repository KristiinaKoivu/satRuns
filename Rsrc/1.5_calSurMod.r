# Run settings (if modifiedSettings is not set to TRUE in batch job script, default settings from Github will be used)
source_url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/Rsrc/settings.r")
if(modifiedSettings) {
  source("/scratch/project_2000994/PREBASruns/assessCarbon/Rsrc/mainSettings.r") # in CSC
}

# Run functions 
source_url("https://raw.githubusercontent.com/ForModLabUHel/satRuns/master/Rsrc/functions.r")


###check and create output directories
setwd(generalPath)
mkfldr <- paste0("surErrMods/")
if(!dir.exists(file.path(generalPath, mkfldr))) {
  dir.create(file.path(generalPath, mkfldr), recursive = TRUE)
}

yearX <- 3

load(paste0(procDataPath,"init",startingYear,"/calST/samples.rdata"))  


sampleID <- 10
rcpfile="CurrClim"

# resample siteType using a uniform distribution 
samples[[sampleID]]$siteType <- sample(1:5,nrow(samples[[sampleID]]),replace = T)

if(rcpfile=="CurrClim"){
  load(paste(climatepath, rcpfile,".rdata", sep=""))  
  setnames(dat,"id","climID")
} else{
  load(paste(climatepath, rcpfile, sep=""))  
}

## Prepare the same initial state for all harvest scenarios that are simulated in a loop below
data.sample <- samples[[sampleID]]
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

print("model initialized")


###run Model
out <- multiPrebas(initPrebas)$multiOut
Vx <- rowSums(out[,yearX,30,,1])

Vx <- data.table(segID=initPrebas$siteInfo[,1],V3=Vx)
print("runs completed")



####build surrogate model
#library(MASS)
#library(minpack.lm)
### Run settings & functions

# load("C:/Users/minunno/GitHub/satRuns/data/inputUncer.rdata")
load("/scratch/project_2000994/PREBASruns/assessCarbon/data/inputUncer.rdata") # in CSC
# load(paste0(procDataPath,"init",startingYear,"/","st",siteTypeX,"/XYsegID.rdata"))  
# load(paste0("output/init",startingYear,"/","st",siteTypeX,"/CurrClim_sample1.rdata"))  
# load(paste0("procData/init",startingYear,"/","st",siteTypeX,"/uniqueData.rdata"))  
# load(paste0("outDT/init",startingYear,"/","st",siteTypeX,"/V_NoHarv_CurrClimlayerall.rdata"))  
# Vmod2019 <- rowSums(out[,yearX,6,])
# load(paste0("initPrebas/init",startingYear,"/","st",siteTypeX,"/CurrClim_sample1.rdata"))  
dataX <- data.table(cbind(initPrebas$multiInitVar[,3:5,1],initPrebas$multiInitVar[,5,2],
                          initPrebas$multiInitVar[,5,3],initPrebas$siteInfo[,3],Vx$V3))
setnames(dataX,c("H","D","BAp","BAsp","BAb","st","Vmod"))
# if(!all(unique(dataX$st) %in% unique(uniqueData$siteType))) stop("not all siteTypes of the tile are in the sample")

#### Here we use stepwise regression to construct an emulator for volume prediction
# dataX$lnVmod<-log(dataX$Vmod)
# dataX$alpha<-NA
dataX$st <- factor(dataX$st)
dataX[,BAtot:=(BAp+BAsp+BAb)]
dataX[,BAh:=BAtot*H]
dataX[,N:=BAtot/(pi*(D/200)^2)]
b = -1.605 ###coefficient of Reineke
dataX[,SDI:=N *(D/10)^b]
full.model <-lm(Vmod~H+D+SDI+BAh+BAp+BAsp+BAb+st,data=dataX)
step.model <- stepAIC(full.model, direction = "both",
                      trace = FALSE)
# start<-as.vector(full.model$coefficients)
### Anonther option is to use nonlinear regression, which differed in error assumption. 
# full.model0<-lm(lnVmod~H+D+SDI+BAh+BAp+BAsp+BAb+st,data=dataX)
# start<-as.vector(full.model0$coefficients)
# nonlinear<-nlsLM(Vmod~exp(a+b*H+c*D+d*SDI+e*BAh+f*BAp+g*BAsp+h*BAb+
#                           i2*as.numeric(st==2)+
#                           i3*as.numeric(st==3)+
#                           i4*as.numeric(st==4)+
#                           i5*as.numeric(st==5)
#                         ),
#                data=dataX,start = list(a=start[1],
#                                  b=start[2],
#                                  c=start[3],
#                                  d=start[4],
#                                  e=start[5],
#                                  f=start[6],
#                                  g=start[7],
#                                  h=start[8],
#                                  i2=start[9],
#                                  i3=start[10],
#                                  i4=start[11],
#                                  i5=start[12]
#                                  ))
#     
sV <- step.model$fitted.values
# sV2<-predict(nonlinear)
pV <- dataX$Vmod
# plot(sV2,pV,col=2)
plot(sV,pV)
abline(0,1)
# summary(nonlinear)
summary(step.model)
save(step.model,file="surErrMods/surMod.rdata")


    
