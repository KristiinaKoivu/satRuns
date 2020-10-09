###choose PREBAS version
vPREBAS <- "v0.2.x"   #### choose PREBAS verson to run the model  "master"


#####Settings####
testRun = T ####set to TRUE to test the code on a small raster proportion
if(!exists("CSCrun")){
  CSCrun = F ### set to TRUE if you are running on CSC
}
fracTest <- 0.2 ###fraction of test area
maxSitesRun <- 20000
maxSitesRunTest <- 20000
saveVars <- c(1,11:13,17,30,43,44) ####select variables to save
varHD <- FALSE #### if true will vary H and D of pine and spruce using siteType

###library path in CSC project_2000994
if(CSCrun){
  .libPaths(c("/projappl/project_2000994/project_rpackages", .libPaths()))
  libpath <- .libPaths()[1]
}

##load libraries
library(lfda)
library(mvtnorm)
library(reshape2)
library(plyr)
library(raster)
library(data.table)
require(sm)
require(rgdal)
library(raster)
library(rgdal)
library(parallel)
library(MASS)
library(minpack.lm)
library(sf)
library(fasterize)


###check prebas version and install if needed
if(!CSCrun){
  devtools::install_github("ForModLabUHel/Rprebasso", ref=vPREBAS)
}
require(Rprebasso)

####indicate rasterPath and climID path
generalPath <- "C:/Users/minunno/Documents/research/assessCarbon/data/Finland/AC_training_FI_34VEQ/"
rasterPath <- paste0(generalPath,"rasters/")
procDataPath <- paste0(generalPath,"procData/")
outPath <- paste0(generalPath,"output/")
initPrebasPath <- paste0(generalPath,"initPrebas/")
climatepath = "C:/Users/minunno/Documents/research/extarctWeather/inputs/" #### local fm
# climatepath = "/scratch/project_2000994/RCP/" ####on CSC
climIDpath <- "C:/Users/minunno/Documents/research/FinSeg/some stuff/climID10km.tif"
# climIDpath <- "/scratch/project_2000994/PREBASruns/metadata/" ####on CSC


startYearWeather <- 1971 ###1971 for Finnish weather dataBase
startingYear <- 2016  #2019
year2 <- 2019 ###year of the second measurement
yearEnd <- 2019     #2024


####indicate raster files
tileX = "34VEQ"
areaID <- "FI"
baRast <-  paste0(rasterPath,areaID,"_",tileX,"-",startingYear,"_BA_10M_1CHS_8BITS.tif")
blPerRast <- paste0(rasterPath,areaID,"_",tileX,"-",startingYear,"_BLP_10M_1CHS_8BITS.tif")
dbhRast <- paste0(rasterPath,areaID,"_",tileX,"-",startingYear,"_DIA_10M_1CHS_8BITS.tif")
vRast <- paste0(rasterPath,areaID,"_",tileX,"-",startingYear,"_GSV_10M_1CHS_16BITS.tif")
hRast <- paste0(rasterPath,areaID,"_",tileX,"-",startingYear,"_HGT_10M_1CHS_16BITS.tif")
pinePerRast <- paste0(rasterPath,areaID,"_",tileX,"-",startingYear,"_P_pine_10M_1CHS_8BITS.tif")
sprucePerRast <- paste0(rasterPath,areaID,"_",tileX,"-",startingYear,"_P_spruce_10M_1CHS_8BITS.tif")
siteTypeRast <- paste0(rasterPath,areaID,"_",tileX,"-",startingYear,"_SITE_10M_1CHS_8BITS.tif")
siteTypeRast2 <- paste0(rasterPath,areaID,"_",tileX,"-",year2,"_SITE_10M_1CHS_8BITS.tif")
vRast2 <- paste0(rasterPath,areaID,"_",tileX,"-",year2,"_GSV_10M_1CHS_16BITS.tif")
baRast2 <-  paste0(rasterPath,areaID,"_",tileX,"-",year2,"_BA_10M_1CHS_8BITS.tif")
dbhRast2 <- paste0(rasterPath,areaID,"_",tileX,"-",year2,"_DIA_10M_1CHS_8BITS.tif")
hRast2 <- paste0(rasterPath,areaID,"_",tileX,"-",year2,"_HGT_10M_1CHS_16BITS.tif")
pinePerRast2 <- paste0(rasterPath,areaID,"_",tileX,"-",year2,"_P_pine_10M_1CHS_8BITS.tif")
sprucePerRast2 <- paste0(rasterPath,areaID,"_",tileX,"-",year2,"_P_spruce_10M_1CHS_8BITS.tif")
blPerRast2 <- paste0(rasterPath,areaID,"_",tileX,"-",year2,"_BLP_10M_1CHS_8BITS.tif")
mgmtmaskRast <- paste0(rasterPath, areaID, "_", tileX, "_mgmtmask.tif")


# Source of tile-specific settings. Defined in batch job script. When set to TRUE will overwrite the tile-specific 
# settings in this script (lines: 41-49, 53-56, 72-81) with settings from filepath in mySettings variable.
if(exists("tileSettings")){
  if(tileSettings) {
    source(mySettings)
    
    # Indicate raster files
    baRast <-  paste0(rasterPath,areaID,"_",tileX,"-",startingYear,"_BA_10M_1CHS_8BITS.tif")
    blPerRast <- paste0(rasterPath,areaID,"_",tileX,"-",startingYear,"_BLP_10M_1CHS_8BITS.tif")
    dbhRast <- paste0(rasterPath,areaID,"_",tileX,"-",startingYear,"_DIA_10M_1CHS_8BITS.tif")
    vRast <- paste0(rasterPath,areaID,"_",tileX,"-",startingYear,"_GSV_10M_1CHS_16BITS.tif")
    hRast <- paste0(rasterPath,areaID,"_",tileX,"-",startingYear,"_HGT_10M_1CHS_16BITS.tif")
    pinePerRast <- paste0(rasterPath,areaID,"_",tileX,"-",startingYear,"_P_pine_10M_1CHS_8BITS.tif")
    sprucePerRast <- paste0(rasterPath,areaID,"_",tileX,"-",startingYear,"_P_spruce_10M_1CHS_8BITS.tif")
    siteTypeRast <- paste0(rasterPath,areaID,"_",tileX,"-",startingYear,"_SITE_10M_1CHS_8BITS.tif")
    siteTypeRast2 <- paste0(rasterPath,areaID,"_",tileX,"-",year2,"_SITE_10M_1CHS_8BITS.tif")
    vRast2 <- paste0(rasterPath,areaID,"_",tileX,"-",year2,"_GSV_10M_1CHS_16BITS.tif")
    baRast2 <-  paste0(rasterPath,areaID,"_",tileX,"-",year2,"_BA_10M_1CHS_8BITS.tif")
    dbhRast2 <- paste0(rasterPath,areaID,"_",tileX,"-",year2,"_DIA_10M_1CHS_8BITS.tif")
    hRast2 <- paste0(rasterPath,areaID,"_",tileX,"-",year2,"_HGT_10M_1CHS_16BITS.tif")
    pinePerRast2 <- paste0(rasterPath,areaID,"_",tileX,"-",year2,"_P_pine_10M_1CHS_8BITS.tif")
    sprucePerRast2 <- paste0(rasterPath,areaID,"_",tileX,"-",year2,"_P_spruce_10M_1CHS_8BITS.tif")
    blPerRast2 <- paste0(rasterPath,areaID,"_",tileX,"-",year2,"_BLP_10M_1CHS_8BITS.tif")
    mgmtmaskRast <- paste0(rasterPath, tileX, "_mgmtmask.tif")
  }
}


nYears <-  yearEnd - startingYear ## number of simulation years
domSPrun = 0.
mgmtmask = F # switch for masking of management


resX <- 10 ### pixel resolution in meters

### define weather inputs (CurrClim, or climate models)
weather = "CurrClim"

###set harvests
defaultThin = 0.
ClCut = 0.
harvscen = "NoHarv"

####set values for NAs and convert factor for prebas units
baNA <- c(253:255); baConv<- 1
blPerNA <- c(253:255); blPerConv<- 1
dbhNA <- c(253:255); dbhConv <- 1
vNA <- c(65533:65535); vConv <- 1
hNA <- c(65533:65535); hConv <- 0.1
pinePerNA <- c(253:255); pinePerConv <- 1
sprucePerNA <- c(253:255); sprucePerConv <- 1
siteTypeNA <- c(253:255); siteTypeConv <- 1

####settings for sitetype estimation
stXruns <- TRUE
siteTypeX <- startingYear #startingYear #year2 #startingYear #1:5


# Set TRUE to enable running 1.8_optST, 2_InitPreb and 3_runModel in parallel. Set to FALSE, these scripts run as serial.
parallelRun <- FALSE
coresN <- 20L ### Set number of cores to use in parallel run 

# Set whether to split unique data in 1.1_procData_siteType to smaller parts. If
# TRUE, data is split.
splitRun <- FALSE
# Range/number of split parts. NOTICE: Code doesn't adjust number of split parts by merely 
# adjusting this variable. If number of parts needs to be changed, 1.1_procData_siteType 
# needs to be modified.
if(splitRun){
  splitRange <- 1:10
}

####thresholds for variables to reset stand from plantation
maxDens <- 10000
initH <- 1.5
initDBH <- 0.5
initN <- 2200
initBA <- pi*(initDBH/200)^2*initN

#####settings for data extraction
varDT <- c(11:13,30)   ####variables to extract in DT
extrFun <- c("mean","mean","sum","sum")
layerDT <- "all" ###layerID to report in data.tables, if layerDT==all the all layers are considered as sum or mean according to extrFun

#####settings for raster creation
varRast <- varDT  #c(44,30)   ####variables to extract in DT
yearOut <- yearEnd#c(2019,2024)

#####filter model output raster
minX <- c(0,0,0,0)
maxX <- c(70,40,60,1000)


#####end Settings####

