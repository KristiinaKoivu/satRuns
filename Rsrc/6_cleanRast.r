source("Rsrc/settings.r")
source("Rsrc/functions.r")
setwd(generalPath)

# raster
# load(paste0(procDataPath,startingYear,"/XYsegID.rdata"))
# crsX <- crs(raster(baRast))
# 
clims <- weather
mans <- harvscen
# 
for(ij in yearOut){
  for(varX in 1:length(varRast)){
    rastName <- paste0("outRast/",startingYear,"/",clims,"_",mans,"_var",varNames[varRast[varX]],
                       "_spec",layerDT,"_yearStart",startingYear,"_yearOut",ij,".tif")
    rastX <- raster(rastName)
    # print(length(which(!is.na(getValues(rastX)))))
    # print(length(which(getValues(rastX)>maxX[varX])))
    rastX[rastX < minX[varX]] <- NA
    rastX[rastX > maxX[varX]] <- NA
    # print(length(which(!is.na(getValues(rastX)))))
    rastNameCleaned <- paste0("outRast/",startingYear,"/",clims,"_",mans,"_var",varNames[varRast[varX]],
                       "_spec",layerDT,"_yearStart",startingYear,"_yearOut",ij,"_cleaned.tif")
    writeRaster(rastX,filename = rastNameCleaned,overwrite=T)
  }
  
  
  # rastName <- paste0("outRast/",startingYear,"/","startV_","startYear",startingYear,
  #                    "_layer",layerDT,".tif")
  # rastX <- raster(rastName)
  # print(length(which(!is.na(getValues(rastX)))))
  # # print(length(which(getValues(rastX)>200)))
  # # rastX[rastX < minX[varX]] <- NA
  # # rastX[rastX > maxX[varX]] <- NA
  # print(length(which(!is.na(getValues(rastX)))))
  # 
}


