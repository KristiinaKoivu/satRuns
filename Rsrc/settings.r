###choose PREBAS version
vPREBAS <- "v0.2.x"   #### choose PREBAS verson to run the model  "master"


#####Settings####
testRun = T ####set to TRUE to test the code on a small raster proportion
CSCrun = F ### set to TRUE if you are running on CSC
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
library(reshape2)
library(plyr)
library(raster)
library(data.table)
require(sm)
require(rgdal)
library(raster)
library(rgdal)
library(parallel)
devtools::install_github("collectivemedia/tictoc")
library(tictoc)
library(MASS)


###check prebas version and install if needed
devtools::install_github("ForModLabUHel/Rprebasso", ref=vPREBAS)
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

coresN <- 10L ###Set number of cores to use in parallel run 

startYearWeather <- 1971 ###1971 for Finnish weather dataBase
startingYear <- 2016  #2019
year2 <- 2019 ###year of the second measurement
yearEnd <- 2019     #2024
nYears <-  yearEnd - startingYear ## number of simulation years
domSPrun = 0.

resX <- 10 ### pixel resolution in meters

### define weather inputs (CurrClim, or climate models)
weather = "CurrClim"

###set harvests
defaultThin = 0.
ClCut = 0.
harvscen = "NoHarv"


####indicate raster files
baRast <-  paste0(rasterPath,"FI_34VEQ-2016_BA_10M_1CHS_8BITS.tif")
blPerRast <- paste0(rasterPath,"FI_34VEQ-2016_BLP_10M_1CHS_8BITS.tif")
dbhRast <- paste0(rasterPath,"FI_34VEQ-2016_DIA_10M_1CHS_8BITS.tif")
vRast <- paste0(rasterPath,"FI_34VEQ-2016_GSV_10M_1CHS_16BITS.tif")
hRast <- paste0(rasterPath,"FI_34VEQ-2016_HGT_10M_1CHS_16BITS.tif")
pinePerRast <- paste0(rasterPath,"FI_34VEQ-2016_P_pine_10M_1CHS_8BITS.tif")
sprucePerRast <- paste0(rasterPath,"FI_34VEQ-2016_P_spruce_10M_1CHS_8BITS.tif")
siteTypeRast <- paste0(rasterPath,"FI_34VEQ-2016_SITE_10M_1CHS_8BITS.tif")
siteTypeRast2 <- paste0(rasterPath,"FI_34VEQ-2019_SITE_10M_1CHS_8BITS.tif")
vRast2 <- paste0(rasterPath,"FI_34VEQ-2019_GSV_10M_1CHS_16BITS.tif")

####set values for NAs and convert factor for prebas units
baNA <- c(253:255); baConv<- 1
blPerNA <- c(253:255); blPerConv<- 1
dbhNA <- c(253:255); dbhConv <- 1
vNA <- c(65533:65535); vConv <- 1
hNA <- c(65533:65535); hConv <- 0.1
pinePerNA <- c(253:255); pinePerConv <- 1
sprucePerNA <- c(253:255); sprucePerConv <- 1
siteTypeNA <- c(254:255); siteTypeConv <- 1

####settings for sitetype estimation
stXruns <- TRUE
siteTypeX <- startingYear #startingYear #year2 #startingYear #1:5

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

