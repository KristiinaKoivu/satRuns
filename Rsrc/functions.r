
####function to initialize the model
create_prebas_input.f = function(r_no, clim, data.sample, nYears, startingYear=0,domSPrun=0) { # dat = climscendataset
  #domSPrun=0 initialize model for mixed forests according to data inputs 
  #domSPrun=1 initialize model only for dominant species 
  nSites <- nrow(data.sample)
  
  ###site Info matrix. nrow = nSites, cols: 1 = siteID; 2 = climID; 3=site type;
  ###4 = nLayers; 5 = nSpecies;
  ###6=SWinit;   7 = CWinit; 8 = SOGinit; 9 = Sinit
  
  siteInfo <- matrix(c(NA,NA,NA,160,0,0,20,3,3,413,0.45,0.118),nSites,12,byrow = T)
  #siteInfo <- matrix(c(NA,NA,NA,3,3,160,0,0,20),nSites,9,byrow = T)
  siteInfo[,1] <- data.sample$segID
  siteInfo[,2] <- as.numeric(data.sample[,climID])
  siteInfo[,3] <- data.sample[,siteType]
  
  # litterSize <- matrix(0,3,3)
  # litterSize[1,1:2] <- 30
  # litterSize[1,3] <- 10
  # litterSize[2,] <- 2
  
  ###Initialise model
  # initVardension nSites,variables, nLayers
  # variables: 1 = species; 2 = Age; 3 = H; 4=dbh; 5 = ba; 6 = Hc
  initVar <- array(NA, dim=c(nSites,7,3))
  data.sample[,baP:= (ba * pineP/(pineP+spruceP+blp))]
  data.sample[,baSP:= (ba * spruceP/(pineP+spruceP+blp))]
  data.sample[,baB:= (ba * blp/(pineP+spruceP+blp))]
  data.sample[,dbhP:= dbh]
  data.sample[,dbhSP:= dbh]
  data.sample[,hP:= h]
  data.sample[,hSP:= h]
  
  data.sample[,N:=ba/(pi*(dbh/2)^2/10000)]
  
  areas <- data.sample$area
  
  initVar[,1,] <- as.numeric(rep(1:3,each=nSites))
  initVar[,2,] <- 1.  # round(as.numeric(data.sample[,age]))  ##### set to 1 because we do not know age
  initVar[,3,] <- as.numeric(data.sample[,h])
  # initVar[,3,][which(initVar[,3,]<1.5)] <- 1.5  ####if H < 1.5 set to 1.5
  
  initVar[,4,] <- as.numeric(data.sample[,dbh])
  if(domSPrun==1){
    ##initialize model only for dominant species##
    initVar[,5,] = 0.
    ix = unlist(data.sample[, which.max(c(pineP, spruceP, blp)), by=1:nrow(data.sample)] [, 2])
    for(jx in 1:nSites) initVar[jx,5,ix[jx]] = as.numeric(data.sample[, ba])[jx]
  } else{
    ###initialize model for mixed forest runs
    initVar[,5,1] <- as.numeric(data.sample[,(ba * pineP/(pineP+spruceP+blp))])
    initVar[,5,2] <- as.numeric(data.sample[,(ba * spruceP/(pineP+spruceP+blp))])
    initVar[,5,3] <- as.numeric(data.sample[,(ba * blp/(pineP+spruceP+blp))])
    ####increase spruceP dbh 10% for spruceP sitetype 1:2
    data.sample[pineP>0. & spruceP >0. & siteType<2.5 & baSP > baP,X:=(ba-1.1*baSP-baB)/baP]
    data.sample[pineP>0. & spruceP >0. & siteType<2.5 & baSP > baP,dbhP:=X*dbh]
    data.sample[pineP>0. & spruceP >0. & siteType<2.5 & baSP > baP,dbhSP:=1.1*dbh]
    data.sample[pineP>0. & spruceP >0. & siteType<2.5  & baSP > baP & dbhP<0.5,dbhSP:=((ba-(0.5/dbh)*baP-baB)/baSP)*dbh]
    data.sample[pineP>0. & spruceP >0. & siteType<2.5  & baSP > baP & dbhP<0.5,dbhP:=0.5]
    
    data.sample[pineP>0. & spruceP >0. & siteType<2.5 & baSP <= baP,dbhSP:=dbh * (ba - 0.9*baP - baB)/baSP]
    data.sample[pineP>0. & spruceP >0. & siteType<2.5 & baSP <= baP,dbhP:=pmax(0.9*dbh,0.3)]
    
    ####increase spruceP h 10% for spruceP sitetype 1:2
    data.sample[pineP>0. & spruceP >0. & siteType<2.5 & baSP > baP,X:=(ba-1.1*baSP-baB)/baP]
    data.sample[pineP>0. & spruceP >0. & siteType<2.5 & baSP > baP,hP:=X*h]   
    data.sample[pineP>0. & spruceP >0. & siteType<2.5 & baSP > baP,hSP:=1.1*h]
    data.sample[pineP>0. & spruceP >0. & siteType<2.5 & baSP > baP & hP<1.5,hSP:=((ba-(1.5/h)*baP-baB)/baSP)*h]
    data.sample[pineP>0. & spruceP >0. & siteType<2.5 & baSP > baP & hP<1.5,hP:=1.5]
    
    data.sample[pineP>0. & spruceP >0. & siteType<2.5 & baSP <= baP,hSP:=h * (ba - 0.9*baP - baB)/baSP]
    data.sample[pineP>0. & spruceP >0. & siteType<2.5 & baSP <= baP,hP:=pmax(0.9*h,1.3)]
    
    ####increase spruceP dbh 5% for spruceP sitetype 3
    data.sample[pineP>0. & spruceP >0. & siteType==3 & baSP > baP,X:=(ba-1.05*baSP-baB)/baP]
    data.sample[pineP>0. & spruceP >0. & siteType==3 & baSP > baP,dbhP:=X*dbh]   
    data.sample[pineP>0. & spruceP >0. & siteType==3 & baSP > baP,dbhSP:=1.05*dbh]
    data.sample[pineP>0. & spruceP >0. & siteType==3 & baSP > baP & dbhP<0.5,dbhSP:=((ba-(0.5/dbh)*baP-baB)/baSP)*dbh]
    data.sample[pineP>0. & spruceP >0. & siteType==3 & baSP > baP & dbhP<0.5,dbhP:=0.5]
    
    data.sample[pineP>0. & spruceP >0. & siteType==3 & baSP <= baP,dbhSP:=dbh * (ba - 0.95*baP - baB)/baSP]
    data.sample[pineP>0. & spruceP >0. & siteType==3 & baSP <= baP,dbhP:=pmax(0.95*dbh,0.3)]
    
    ####increase spruceP h 5% for spruceP sitetype 3
    data.sample[pineP>0. & spruceP >0. & siteType==3,X:=(ba-1.05*baSP-baB)/baP]
    data.sample[pineP>0. & spruceP >0. & siteType==3,hP:=X*h]  
    data.sample[pineP>0. & spruceP >0. & siteType==3,hSP:=1.05*h]  
    data.sample[pineP>0. & spruceP >0. & siteType==3 & hP<1.5,hSP:=((ba-(1.5/h)*baP-baB)/baSP)*h]
    data.sample[pineP>0. & spruceP >0. & siteType==3 & hP<1.5,hP:=1.5]
    
    data.sample[pineP>0. & spruceP >0. & siteType==3 & baSP <= baP,hSP:=h * (ba - 0.95*baP - baB)/baSP]
    data.sample[pineP>0. & spruceP >0. & siteType==3 & baSP <= baP,hP:=pmax(0.95*h,1.3)]
    
    ####increase pineP dbh 10% for sitetype >= 4
    data.sample[pineP>0. & spruceP >0. & siteType>3.5 & baP > baSP,X:=(ba-1.1*baP-baB)/baSP]
    data.sample[pineP>0. & spruceP >0. & siteType>3.5 & baP > baSP,dbhSP:=X*dbh]   
    data.sample[pineP>0. & spruceP >0. & siteType>3.5 & baP > baSP,dbhP:=1.1*dbh]   
    data.sample[pineP>0. & spruceP >0. & siteType>3.5 & baP > baSP & dbhSP<0.5,dbhP:=((ba-(0.5/dbh)*baSP-baB)/baP)*dbh]
    data.sample[pineP>0. & spruceP >0. & siteType>3.5 & baP > baSP & dbhSP<0.5,dbhSP:=0.5]
    data.sample[pineP>0. & spruceP >0. & siteType>3.5 & baP > baSP,dbhP:=dbh * (ba - 0.9*baSP - baB)/baP]
    data.sample[pineP>0. & spruceP >0. & siteType>3.5 & baP > baSP,dbhSP:=pmax(0.9*dbh,0.3)]
    ####increase pineP h 10% for sitetype >= 4
    data.sample[pineP>0. & spruceP >0. & siteType>3.5 & baP > baSP,X:=(ba-1.1*baP-baB)/baSP]
    data.sample[pineP>0. & spruceP >0. & siteType>3.5 & baP > baSP,hSP:=X*h]   
    data.sample[pineP>0. & spruceP >0. & siteType>3.5 & baP > baSP,hP:=1.1*h]   
    data.sample[pineP>0. & spruceP >0. & siteType>3.5 & baP > baSP & hSP<1.5,hP:=((ba-(1.5/h)*baSP-baB)/baP)*h]
    data.sample[pineP>0. & spruceP >0. & siteType>3.5 & baP > baSP & hSP<1.5,hSP:=1.5]
    data.sample[pineP>0. & spruceP >0. & siteType>3.5 & baP > baSP,hP:=h * (ba - 0.9*baSP - baB)/baP]
    data.sample[pineP>0. & spruceP >0. & siteType>3.5 & baP > baSP,hSP:=pmax(0.9*h,1.3)]
    initVar[,3,1] <- as.numeric(data.sample[,hP])
    initVar[,3,2] <- as.numeric(data.sample[,hSP])
    initVar[,4,1] <- as.numeric(data.sample[,dbhP])
    initVar[,4,2] <- as.numeric(data.sample[,dbhSP])
    
  }
  
  # initVar[,6,] <- as.numeric(data.sample[,hc])
  
  ###check which BA ==0. and set to 0 the rest of the variables
  NoPine <- which(initVar[,5,1]==0.)
  NoSpruce <- which(initVar[,5,2]==0.)
  NoDecid <- which(initVar[,5,3]==0.)
  
  siteInfo[NoPine,8:9] <- siteInfo[NoPine,8:9] - 1
  siteInfo[NoSpruce,8:9] <- siteInfo[NoSpruce,8:9] - 1
  siteInfo[NoDecid,8:9] <- siteInfo[NoDecid,8:9] - 1

  #siteInfo[NoPine,4] <- siteInfo[NoPine,4] - 1
  #siteInfo[NoSpruce,4] <- siteInfo[NoSpruce,4] - 1
  #siteInfo[NoDecid,4] <- siteInfo[NoDecid,4] - 1
  initVar[NoPine,3:6,1] <- 0.
  initVar[NoSpruce,3:6,2] <- 0.
  initVar[NoDecid,3:6,3] <- 0.
  initVar[NoSpruce,,2] <- initVar[NoSpruce,,3]
  initVar[NoPine,,1:2] <- initVar[NoPine,,2:3]
  
  nLay1 <- which(siteInfo[,8]==1)
  nLay2 <- which(siteInfo[,8]==2)
  initVar[nLay1,3:7,2:3] <- 0
  initVar[nLay2,3:7,3] <- 0
  # initVar[which(initVar[,5,1]==0.),,1] <- initVar[which(initVar[,5,1]==0.),,2]
  # initVar[which(initVar[,5,1]==0.),,2] <- initVar[which(initVar[,5,1]==0.),,3]
  # initVar[which(initVar[,5,1]==0.),1,3] <- 1
  # initVar[which(initVar[,5,1]==0.),3:6,3] <- 0
  
  # if (FALSE) {
  #   dat = dat[climID %in% data.sample[, unique(climID)]]
  #   
  #   if(weather!= "CurrClim.rdata"){
  #     # dat[, pvm:= as.Date('1980-01-01') - 1 + rday ]
  #     # dat[, DOY:= as.numeric(format(pvm, "%j"))]
  #     dat[, Year:= as.numeric(floor(rday/366)+1971)]
  #     dat = dat[Year >= startingYear]
  #     dat[DOY==366, DOY:=365]
  #   }
  #   PARtran = t( dcast(dat[, list(climID, rday, PAR)], rday ~ climID,
  #                      value.var="PAR")[, -1])
  #   TAirtran = t( dcast(dat[, list(climID, rday, TAir)], rday ~ climID,
  #                       value.var="TAir")[, -1])
  #   VPDtran = t( dcast(dat[, list(climID, rday, VPD)], rday ~ climID,
  #                      value.var="VPD")[, -1])
  #   Preciptran = t( dcast(dat[, list(climID, rday, Precip)], rday ~ climID,
  #                         value.var="Precip")[, -1])
  #   CO2tran = t( dcast(dat[, list(climID, rday, CO2)], rday ~ climID,
  #                      value.var="CO2")[, -1])
  # }
  siteInfo[, 2]  = match(as.numeric(siteInfo[, 2]), as.numeric(rownames(clim[[1]])))
  # siteInfo[, 2]  = match(siteInfo[,2], unique(dat$climID))
  
  # defaultThin=as.numeric(1-data.sample[, cons])
  # energyCut <- ClCut <- as.numeric(1-data.sample[, cons])
  ## Set to match climate data years
  initPrebas <- InitMultiSite(nYearsMS = rep(nYears,nSites),siteInfo=siteInfo,
                              # litterSize = litterSize,#pAWEN = parsAWEN,
                              defaultThin = defaultThin,
                              ClCut = ClCut,
                              areas =areas,
                              # energyCut = energyCut, 
                              multiInitVar = as.array(initVar),
                              PAR = clim$PAR[, 1:(nYears*365)],
                              TAir=clim$TAir[, 1:(nYears*365)],
                              VPD=clim$VPD[, 1:(nYears*365)],
                              Precip=clim$Precip[, 1:(nYears*365)],
                              CO2=clim$CO2[, 1:(nYears*365)],
                              yassoRun = 1)
  initPrebas
}




yasso.mean.climate.f = function(dat, data.sample, startingYear, nYears){
  dat = dat[climID %in% data.sample[, unique(climID)]]
  dat[, DOY:=rep(1:365, len=dim(dat)[1])]
  dat[, Year:=rep(1980:2099, each=365)]
  #dat[, Year:= as.numeric(format(pvm, "%Y"))]
  dat = dat[Year >= startingYear & Year <= startingYear+nYears]
  dat[, pvm:= as.Date(paste(Year, '-01-01', sep="")) - 1 + DOY ]
  #dat[, DOY:= as.numeric(format(pvm, "%j"))]
  dat[, Mon:= as.numeric(format(pvm, "%m"))]
  #dat[DOY==366, DOY:=365]
  Tmean = dat[, mean(TAir), by = Year]
  Tsum = dat[, sum(ifelse(TAir>5, TAir-5, 0)), by=.(climID, Year)][, mean(V1), by=Year]
  PAR = dat[, mean(PAR), by = Year]
  VPD = dat[, mean(VPD), by = Year]
  CO2 = dat[, mean(CO2), by = Year]
  Precip = dat[, sum(Precip), by = .(climID, Year)][, mean(V1), by=Year]
  Tampl = dat[, .(mean(TAir)), by = .(climID, Year, Mon)][, (max(V1)-min(V1))/2, by=Year]
  
  out = cbind(Tmean, Precip[, -1], Tampl[, -1], CO2[, -1], PAR[, -1], VPD[, -1], Tsum[, -1])
  colnames(out) = c('Year','Tmean','Precip','Tampl', 'CO2', "PAR", "VPD", "Tsum5")
  out
}


prep.climate.f = function(dat, data.sample, startingYear, nYears){
  dat = dat[climID %in% data.sample[, unique(climID)]]
  if(weather== "CurrClim"){
    dat[, Year:= as.numeric(floor((rday-1)/365)+1971)]
    dat1 = dat[Year >= startingYear]
    if(nYears>length(unique(dat1$Year))){
      nSampleYear <- nYears - length(unique(dat1$Year))
      set.seed(123)
      yearX <- sample(1971:min(startingYear,max(dat$Year)),nSampleYear,replace = F)
      lastYear <- max(dat$Year,startingYear)
      newYears <- lastYear + 1:length(yearX)
      dat2 <- dat[Year %in% yearX,]
      dat2$Year <- newYears[match(dat2$Year,yearX)]
      dat <- rbind(dat1,dat2)
      dat[,DOY:=rep(1:365,nYears),by=climID]
      setorder(dat,climID,Year,DOY)
      dat[,rday:=1:(365*nYears),by=climID]
    }

  }else{
  dat[, pvm:= as.Date('1980-01-01') - 1 + rday ]
  dat[, DOY:= as.numeric(format(pvm, "%j"))]
  dat[, Year:= as.numeric(format(pvm, "%Y"))]
  dat = dat[Year >= startingYear]
  dat[DOY==366, DOY:=365]
  }
  climID = dat[,unique(climID)]
  PARtran = t( dcast(dat[, list(climID, rday, PAR)], rday ~ climID,
                     value.var="PAR")[, -1])
  TAirtran = t( dcast(dat[, list(climID, rday, TAir)], rday ~ climID,
                      value.var="TAir")[, -1])
  VPDtran = t( dcast(dat[, list(climID, rday, VPD)], rday ~ climID,
                     value.var="VPD")[, -1])
  Preciptran = t( dcast(dat[, list(climID, rday, Precip)], rday ~ climID,
                        value.var="Precip")[, -1])
  CO2tran = t( dcast(dat[, list(climID, rday, CO2)], rday ~ climID,
                     value.var="CO2")[, -1])
  list(PAR=PARtran, TAir=TAirtran, VPD=VPDtran, 
       Precip=Preciptran, CO2=CO2tran,climID=climID)
}


##function to compile all data and create data.table 
createDT <- function(climate, management,variable, layer,startingYear){
  
  files <- list.files(path= "output/")#,pattern = paste0("year",startingYear,"_"))
  
  # files2 <- intersect(list.files(path= "output/", pattern = climate), list.files(path= "output/",pattern = management))
  # files=intersect(files,files2)
  
  for (ij in variable) assign(varNames[ij],data.table())
  # segID <- areas <-numeric(0)
  
  for(i in 1:length(files)){
    sampleID <- paste0("sample",i,".")
    
    fileX <- files[grep(sampleID,files,fixed = T)]
    
    load(paste0("output/",fileX))
    
    margin= 1:2#(length(dim(out$annual[,,variable,]))-1)
    if(layer=="tot"){
      for (ij in variable){ 
        if(varNames[ij] %in% dimnames(out)$varX){
          varX <- which(dimnames(out)$varX==varNames[ij])
          assign(varNames[ij],data.table(rbind(eval(parse(text = varNames[ij])),
                                               apply(out[,,varX,],margin,sum))))
        }else{
          print(paste(varNames[ij],"not saved"))
        }
      }
    }else{
      for(ij in variable){
        if(varNames[ij] %in% dimnames(out)$varX){
          varX <- which(dimnames(out)$varX==varNames[ij])
          assign(varNames[ij],data.table(rbind(eval(parse(text = varNames[ij])),
                                               out[,,varX,layer])))
        }else{
          print(paste(varNames[ij],"not saved"))
        }
      } 
    }
    print(i)
  }
  
  for(ij in variable) save(list=varNames[ij],file=paste0("outDT/",varNames[ij],"_",management,"_",climate,
                                                         "StartYear",startingYear,"layer",layer,".rdata"))
}


# this function create raster in tif format from data.tables selecting one year or the average of a time priod if yearOut is a vector of years
createTifFromDT <- function(climate, management, yearOut, varX, layerDT, stYear,XYsegID,crsX=NA){
  simYear <- yearOut - stYear
  fileDT=paste0("outDT/",varNames[varX],"_",management,"_",climate,
                "StartYear",startingYear,"layer",layerDT,".rdata")
  load(fileDT)
  
  outX <- t(get(varNames[varX]))
  if (length(simYear)==1) outX <- outX[simYear,]
  if (length(simYear)>1) outX <- colMeans(outX[simYear,],na.rm = T)
  # outX <- data.table(cbind(segID,areas,outX))
  outX <- data.table(cbind(segID,outX))
  
  setnames(outX,c("segID",varNames[varX]))
  
  setkey(XYsegID,segID)
  setkey(outX,segID)
  outXY <- merge(XYsegID,outX,all = T)
  ###remove coordinates NAs
  outXY <- outXY[!is.na(x)]
  
  ###create raster 
  rastX <- rasterFromXYZ(outXY[,c("x","y",varNames[varX]),with=F])
  crs(rastX) <- crsX
  
  rastName <- paste0("outRast/",climate,"_",management,"_var",varNames[varX],
                     "_spec",layerDT,"_yearStart",stYear,".tif")
  writeRaster(rastX,filename = rastName,overwrite=T)
}

