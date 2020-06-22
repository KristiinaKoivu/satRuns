###choose PREBAS version
vPREBAS <- "v0.2.x"   #### choose PREBAS verson to run the model  "master"


#####Settings####
testRun = T ####set to TRUE to test the code on a small raster proportion
fracTest <- 0.2 ###fraction of test area
maxSitesRun <- 20000
maxSitesRunTest <- 1000
saveVars <- c(1,11:13,17,30,43,44) ####select variables to save


###library path in CSC project_2000994
if(!testRun){
  .libPaths(c("/projappl/project_2000994/project_rpackages", .libPaths()))
  libpath <- .libPaths()[1]
}

####indicate rasterPath and climID path
generalPath <- "C:/Users/minunno/Documents/research/assessCarbon/data/Finland/AC_training_FI_34VEQ/"
climIDpath <- "C:/Users/minunno/Documents/research/FinSeg/some stuff/climID10km.tif"
rasterPath <- paste0(generalPath,"rasters/")
procDataPath <- paste0(generalPath,"procData/")
outPath <- paste0(generalPath,"output/")
initPrebasPath <- paste0(generalPath,"initPrebas/")
# climatepath = "/scratch/project_2000994/RCP/" ####on CSC
climatepath = "C:/Users/minunno/Documents/research/extarctWeather/inputs/" #### local fm


startingYear <- 2016
yearEnd <- 2024
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

####set values for NAs and convert factor for prebas units
baNA <- c(253:255); baConv<- 1
blPerNA <- c(253:255); blPerConv<- 1
dbhNA <- c(253:255); dbhConv <- 1
vNA <- c(65533:65535); vConv <- 1
hNA <- c(65533:65535); hConv <- 0.1
pinePerNA <- c(253:255); pinePerConv <- 1
sprucePerNA <- c(253:255); sprucePerConv <- 1
siteTypeNA <- c(254:255); siteTypeConv <- 1

####thresholds for variables to reset stand from plantation
maxDens <- 10000
initH <- 1.5
initDBH <- 0.5
initN <- 2200
initBA <- pi*(initDBH/200)^2*initN

#####settings for data extraction
varDT <- c(44,30)   ####variables to extract in DT
layerDT <- "tot" ###layerID to report in data.tables, if layerDT==tot the totals of all layers is provided

#####settings for raster creation
varRast <- c(44,30)   ####variables to extract in DT
yearOut <- startingYear + 1:3

#####end Settings####
