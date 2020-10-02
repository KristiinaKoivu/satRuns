
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
  initVar[,2,] <- as.numeric(data.sample[,h])*3.3  # round(as.numeric(data.sample[,age]))  ##### set to 1 because we do not know age
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

    if(varHD){ #### if true will vary H and D of pine and spruce using siteType
      ###increase spruceP dbh 10% for spruceP sitetype 1:2
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


prep.climate.f = function(dat, data.sample, startingYear, nYears,startYearWeather){
  dat = dat[climID %in% data.sample[, unique(climID)]]
  if(weather == "CurrClim"){
    dat[, Year:= as.numeric(floor((rday-1)/365)+startYearWeather)]
    dat1 = dat[Year >= startingYear]
    if(nYears>length(unique(dat1$Year))){
      nSampleYear <- nYears - length(unique(dat1$Year))
      set.seed(123)
      yearX <- sample(startYearWeather:min(startingYear,max(dat$Year)),nSampleYear,replace = F)
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
createDT <- function(climate, management,variable, layer,startingYear,funX,siteTypeX){
  
  loadFolder <- paste0(outPath,"init",startingYear,"/",
                       "st",siteTypeX,"/")
  files <- list.files(path= loadFolder)#,pattern = paste0("year",startingYear,"_"))
  startV <- data.table()
  for (ij in variable) assign(varNames[ij],data.table())

  for(i in 1:length(files)){
    sampleID <- paste0("sample",i,".")
    
    fileX <- files[grep(sampleID,files,fixed = T)]
    
    load(paste0(loadFolder,fileX))
    
    if(exists("v0")) startV <- rbind(startV,v0)
    margin= 1:2#(length(dim(out$annual[,,variable,]))-1)
    if(layer=="all"){
      ijX <- 0
      for (ij in variable){ 
        ijX = ijX+1
        if(varNames[ij] %in% dimnames(out)$varX){
          varX <- which(dimnames(out)$varX==varNames[ij])
          if(funX[ijX]=="sum") assign(varNames[ij],data.table(rbind(eval(parse(text = varNames[ij])),
                                               apply(out[,,varX,],margin,sum,na.rm=T))))
          
          if(funX[ijX]=="mean"){
            BAind <- which(saveVars==13)
            if(length(BAind)<1){
              print(paste("BA not saved!! for",varNames[saveVars[varX]],"aritmetic mean was used."))
                    assign(varNames[ij],data.table(rbind(eval(parse(text = varNames[ij])),
                               apply(out[,,varX,],margin,mean,na.rm=T))))
            }else{
              sumBA <- apply(out[,,BAind,],1:2,sum,na.rm=T)
              sumBAs <- array(sumBA,dim = dim(out[,,BAind,]))
              fracBA <- out[,,BAind,]/sumBAs
              wgdVar <- out[,,varX,] * fracBA
              assign(varNames[ij],data.table(rbind(eval(parse(text = varNames[ij])),
                                                   apply(wgdVar,margin,sum,na.rm=T))))
          } }
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
  
  for(ij in variable){
    file2save <- paste0("outDT/init",startingYear,"/","st",siteTypeX,"/",varNames[ij],"_",management,"_",climate,
                        "layer",layer,".rdata")
    save(list=varNames[ij],file=file2save)
  }
  if((dim(startV)[1]) > 0 & siteTypeX==startingYear) save(startV,file=paste0("outDT/","init",startingYear,"/","startV_layerall.rdata"))
}


# this function create raster in tif format from data.tables selecting one year or the average of a time priod if yearOut is a vector of years
createTifFromDT <- function(climate, management, yearOut, varX, layerDT, stYear,XYsegID,crsX=NA){
  simYear <- yearOut - stYear
  fileDT=paste0("outDT/","init",startingYear,"/st",siteTypeX,"/",varNames[varX],"_",management,"_",climate,
                "layer",layerDT,".rdata")
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
  
  rastName <- paste0("outRast/","init",startingYear,"/st",siteTypeX,"/",climate,"_",management,"_var",varNames[varX],
                     "_spec",layerDT,"_yearStart",stYear,"_yearOut",yearOut,".tif")
  writeRaster(rastX,filename = rastName,overwrite=T)
}




####functions used in the siteType calculations
fixBAper <- function(BApers){
  minBA <- min(BApers)
  if(minBA<0) BApers <- BApers - minBA
  BApers <- BApers/sum(BApers)*100
  return(BApers)
}

###function for site type data assimilation
pSTx <- function(segIDx,nSample,year1,year2,tileX){
  mu1 <- errData[[paste0("y",year1)]][[paste0("t",tileX)]]$muFSVda
  sigma1 <- errData[[paste0("y",year1)]][[paste0("t",tileX)]]$sigmaFSVda
  mu2 <- errData[[paste0("y",year2)]][[paste0("t",tileX)]]$muSTda
  sigma2 <- errData[[paste0("y",year2)]][[paste0("t",tileX)]]$sigmaSTda
  set.seed(1234)
  sampleError <- data.table(mvrnorm(nSample*2,mu=mu1,Sigma=sigma1))
  # segIDx <- dataSurV[segID==2]
  sampleX <- data.table()
  sampleX$H <- segIDx$H + sampleError$H
  sampleX$D <- segIDx$D + sampleError$D
  sampleX$BAtot <- segIDx$BAtot + sampleError$G
  sampleX$BApPer <- segIDx$BApPer + sampleError$BAp
  sampleX$BAspPer <- segIDx$BAspPer + sampleError$BAsp
  sampleX$BAbPer <- segIDx$BAbPer + sampleError$BAb
  sampleX <- sampleX[H>1.5]
  sampleX <- sampleX[D>0.5]
  sampleX <- sampleX[BAtot>0.045]
  sampleX <- sampleX[1:min(nSample,nrow(sampleX))]
  if(nrow(sampleX)<nSample){
    sample1 <- sampleX
    set.seed(1234)
    sampleError <- data.table(mvrnorm(nSample*2,mu=mu1,Sigma=sigma1))
    # segIDx <- dataSurV[segID==2]
    sampleX <- data.table()
    sampleX$H <- segIDx$H + sampleError$H
    sampleX$D <- segIDx$D + sampleError$D
    sampleX$BAtot <- segIDx$BAtot + sampleError$G
    sampleX$BApPer <- segIDx$BApPer + sampleError$BAp
    sampleX$BAspPer <- segIDx$BAspPer + sampleError$BAsp
    sampleX$BAbPer <- segIDx$BAbPer + sampleError$BAb
    sampleX <- sampleX[H>1.5]
    sampleX <- sampleX[D>0.5]
    sampleX <- sampleX[BAtot>0.045]
    sampleX <- rbind(sample1,sampleX)
    sampleX <- sampleX[1:min(nSample,nrow(sampleX))]
  }
  
  sampleX[, c("BApPer", "BAspPer", "BAbPer"):=
            as.list(fixBAper(unlist(.(BApPer,BAspPer,BAbPer)))), 
          by = seq_len(nrow(sampleX))]
  
  max.pro.est<-apply(segIDx[, c('BApPer','BAspPer','BAbPer')], 1, which.max)
  segIDx[,max.pro.est:=max.pro.est]
  
  if(year1=="all" & tileX=="all"){
    logistic.model <- logisticPureF$all
  }else{
    logistic.model <- logisticPureF[[paste0("y",year1)]][[paste0("t",tileX)]]
  }
  
  set.seed(1234)
  sampleX$pureF <- runif(min(nSample,nrow(sampleX)),0,1)<predict(logistic.model,type="response",newdata = segIDx)
  if(max.pro.est==1) sampleX[which(pureF),c("BApPer","BAspPer","BAbPer"):=list(100,0,0)]
  if(max.pro.est==2) sampleX[which(pureF),c("BApPer","BAspPer","BAbPer"):=list(0,100,0)]
  if(max.pro.est==3) sampleX[which(pureF),c("BApPer","BAspPer","BAbPer"):=list(0,0,100)]
  
  sampleX[,BAp:=BApPer*BAtot/100]
  sampleX[,BAsp:=BAspPer*BAtot/100]
  sampleX[,BAb:=BAbPer*BAtot/100]
  # sampleX[,st:=segIDx$st]
  sampleX[,V2:=segIDx$V2]
  sampleX[,ba2:=segIDx$ba2]
  sampleX[,h2:=segIDx$h2]
  sampleX[,dbh2:=segIDx$dbh2]
  sampleX[,segID:=segIDx$segID]
  
  # sampleX$lnVmod<-log(sampleX$Vmod)
  # sampleX$st<-factor(sampleX$st,levels = 1:5)     ##!!!!Xianglin
  # sampleX$st <- factor(sampleX$st)
  sampleX[,BAtot:=(BAp+BAsp+BAb)]
  sampleX[,BAh:=BAtot*H]
  sampleX[,N:=BAtot/(pi*(D/200)^2)]
  b = -1.605 ###coefficient of Reineke
  sampleX[,SDI:=N *(D/10)^b]
  
  # full.model<-lm(lnVmod~H+D+lnBAp+lnBAsp+lnBAb+st,data=dataX)
  sampleX$st <- factor(1)
  sampleX[,VsurST1 := pmax(0.,predict(step.modelV,newdata=sampleX))]
  sampleX[,BsurST1 := pmax(0.,predict(step.modelB,newdata=sampleX))]
  sampleX[,DsurST1 := pmax(0.,predict(step.modelD,newdata=sampleX))]
  sampleX[,HsurST1 := pmax(0.,predict(step.modelH,newdata=sampleX))]
  sampleX$st <- factor(2)
  sampleX[,VsurST2 := pmax(0.,predict(step.modelV,newdata=sampleX))]
  sampleX[,BsurST2 := pmax(0.,predict(step.modelB,newdata=sampleX))]
  sampleX[,DsurST2 := pmax(0.,predict(step.modelD,newdata=sampleX))]
  sampleX[,HsurST2 := pmax(0.,predict(step.modelH,newdata=sampleX))]
  sampleX$st <- factor(3)
  sampleX[,VsurST3 := pmax(0.,predict(step.modelV,newdata=sampleX))]
  sampleX[,BsurST3 := pmax(0.,predict(step.modelB,newdata=sampleX))]
  sampleX[,DsurST3 := pmax(0.,predict(step.modelD,newdata=sampleX))]
  sampleX[,HsurST3 := pmax(0.,predict(step.modelH,newdata=sampleX))]
  sampleX$st <- factor(4)
  sampleX[,VsurST4 := pmax(0.,predict(step.modelV,newdata=sampleX))]
  sampleX[,BsurST4 := pmax(0.,predict(step.modelB,newdata=sampleX))]
  sampleX[,DsurST4 := pmax(0.,predict(step.modelD,newdata=sampleX))]
  sampleX[,HsurST4 := pmax(0.,predict(step.modelH,newdata=sampleX))]
  sampleX$st <- factor(5)
  sampleX[,VsurST5 := pmax(0.,predict(step.modelV,newdata=sampleX))]
  sampleX[,BsurST5 := pmax(0.,predict(step.modelB,newdata=sampleX))]
  sampleX[,DsurST5 := pmax(0.,predict(step.modelD,newdata=sampleX))]
  sampleX[,HsurST5 := pmax(0.,predict(step.modelH,newdata=sampleX))]
  
  dx1 <- cbind(sampleX$BsurST1 - segIDx$ba2,sampleX$DsurST1 - segIDx$dbh2,
               sampleX$HsurST1 - segIDx$h2,sampleX$VsurST1 - segIDx$V2)
  dx2 <- cbind(sampleX$BsurST2 - segIDx$ba2,sampleX$DsurST2 - segIDx$dbh2,
               sampleX$HsurST2 - segIDx$h2,sampleX$VsurST2 - segIDx$V2)
  dx3 <- cbind(sampleX$BsurST3 - segIDx$ba2,sampleX$DsurST3 - segIDx$dbh2,
               sampleX$HsurST3 - segIDx$h2,sampleX$VsurST3 - segIDx$V2)
  dx4 <- cbind(sampleX$BsurST4 - segIDx$ba2,sampleX$DsurST4 - segIDx$dbh2,
               sampleX$HsurST4 - segIDx$h2,sampleX$VsurST4 - segIDx$V2)
  dx5 <- cbind(sampleX$BsurST5 - segIDx$ba2,sampleX$DsurST5 - segIDx$dbh2,
               sampleX$HsurST5 - segIDx$h2,sampleX$VsurST5 - segIDx$V2)
  
  pst1 <- mean(dmvnorm(x=dx1, mean=mu2, sigma=sigma2, log=FALSE))
  pst2 <- mean(dmvnorm(x=dx2, mean=mu2, sigma=sigma2, log=FALSE))
  pst3 <- mean(dmvnorm(x=dx3, mean=mu2, sigma=sigma2, log=FALSE))
  pst4 <- mean(dmvnorm(x=dx4, mean=mu2, sigma=sigma2, log=FALSE))
  pst5 <- mean(dmvnorm(x=dx5, mean=mu2, sigma=sigma2, log=FALSE))
  psum <- pst1 +pst2+pst3 +pst4+pst5
  pst1 <- pst1/psum
  pst2 <- pst2/psum
  pst3 <- pst3/psum
  pst4 <- pst4/psum
  pst5 <- pst5/psum
  return(pST=c(segIDx$segID,pst1,pst2,pst3,pst4,pst5)) 
}


###function for structural variables data assimilation 
pSVDA <- function(segIDx,nSample,year1,year2,tileX){
  mu1 <- errData[[paste0("y",year1)]][[paste0("t",tileX)]]$muFSVda
  sigma1 <- errData[[paste0("y",year1)]][[paste0("t",tileX)]]$sigmaFSVda
  mu2 <- errData[[paste0("y",year2)]][[paste0("t",tileX)]]$muFSVda
  sigma2 <- errData[[paste0("y",year2)]][[paste0("t",tileX)]]$sigmaFSVda
  
  pST <- c(segIDx$pST1,segIDx$pST2,segIDx$pST3,segIDx$pST4,segIDx$pST5)
  st <- sample(rep(1:5,round(nSample*pST)),nSample,replace = T)
  set.seed(1234)
  sampleError <- data.table(mvrnorm(nSample*2,mu=mu1,Sigma=sigma1))

  # segIDx <- dataSurV[segID==2]
  sampleX <- data.table()
  sampleX$H <- segIDx$H + sampleError$H
  sampleX$D <- segIDx$D + sampleError$D
  sampleX$BAtot <- segIDx$BAtot + sampleError$G
  sampleX$BApPer <- segIDx$BApPer + sampleError$BAp
  sampleX$BAspPer <- segIDx$BAspPer + sampleError$BAsp
  sampleX$BAbPer <- segIDx$BAbPer + sampleError$BAb
  sampleX <- sampleX[H>1.5]
  sampleX <- sampleX[D>0.5]
  sampleX <- sampleX[BAtot>0.045]
  sampleX <- sampleX[1:min(nSample,nrow(sampleX))]
  if(nrow(sampleX)<nSample){
    sample1 <- sampleX
    set.seed(123)
    sampleError <- data.table(mvrnorm(nSample*2,mu=mu1,Sigma=sigma1))
    # segIDx <- dataSurV[segID==2]
    sampleX <- data.table()
    sampleX$H <- segIDx$H + sampleError$H
    sampleX$D <- segIDx$D + sampleError$D
    sampleX$BAtot <- segIDx$BAtot + sampleError$G
    sampleX$BApPer <- segIDx$BApPer + sampleError$BAp
    sampleX$BAspPer <- segIDx$BAspPer + sampleError$BAsp
    sampleX$BAbPer <- segIDx$BAbPer + sampleError$BAb
    sampleX <- sampleX[H>1.5]
    sampleX <- sampleX[D>0.5]
    sampleX <- sampleX[BAtot>0.045]
    sampleX <- rbind(sample1,sampleX)
    sampleX <- sampleX[1:min(nSample,nrow(sampleX))]
  }

  sampleX[, c("BApPer", "BAspPer", "BAbPer"):=
            as.list(fixBAper(unlist(.(BApPer,BAspPer,BAbPer)))), 
          by = seq_len(nrow(sampleX))]
  
  max.pro.est<-apply(segIDx[, c('BApPer','BAspPer','BAbPer')], 1, which.max)
  segIDx$max.pro.est=max.pro.est
  
  if(year1=="all" & tileX=="all"){
    logistic.model <- logisticPureF$all
  }else{
    logistic.model <- logisticPureF[[paste0("y",year1)]][[paste0("t",tileX)]]
  }
  
  set.seed(1234)
  sampleX$pureF <- runif(min(nSample,nrow(sampleX)),0,1)<predict(logistic.model,type="response",newdata = segIDx)
  if(max.pro.est==1) sampleX[which(pureF),c("BApPer","BAspPer","BAbPer"):=list(100,0,0)]
  if(max.pro.est==2) sampleX[which(pureF),c("BApPer","BAspPer","BAbPer"):=list(0,100,0)]
  if(max.pro.est==3) sampleX[which(pureF),c("BApPer","BAspPer","BAbPer"):=list(0,0,100)]
  
  sampleX[,BAp:=BApPer*BAtot/100]
  sampleX[,BAsp:=BAspPer*BAtot/100]
  sampleX[,BAb:=BAbPer*BAtot/100]
  # sampleX[,st:=segIDx$st]
  sampleX[,V2:=segIDx$V2]
  sampleX[,ba2:=segIDx$ba2]
  sampleX[,h2:=segIDx$h2]
  sampleX[,dbh2:=segIDx$dbh2]
  # sampleX[,BApPer2:=segIDx$BApPer2]
  # sampleX[,BAspPer2:=segIDx$BAspPer2]
  # sampleX[,BAbPer2:=segIDx$BAbPer2]
  sampleX[,segID:=segIDx$segID]
  
  # sampleX$lnVmod<-log(sampleX$Vmod)
  # sampleX$st<-factor(sampleX$st,levels = 1:5)     ##!!!!Xianglin
  # sampleX$st <- factor(sampleX$st)
  sampleX[,BAtot:=(BAp+BAsp+BAb)]
  sampleX[,BAh:=BAtot*H]
  sampleX[,N:=BAtot/(pi*(D/200)^2)]
  b = -1.605 ###coefficient of Reineke
  sampleX[,SDI:=N *(D/10)^b]
  sampleX$st <- st
  sampleX[,rootBAp:=BAp^0.5]
  sampleX[,BAp2:=ba2*BApPer/100]
  sampleX[,BAsp2:=ba2*BAspPer/100]
  sampleX[,BAb2:=ba2*BAbPer/100]
  
  # full.model<-lm(lnVmod~H+D+lnBAp+lnBAsp+lnBAb+st,data=dataX)
  sampleX$st <- factor(sampleX$st)
  sampleX[,Hx := pmax(0.,predict(step.modelH,newdata=sampleX))]
  sampleX[,Dx := pmax(0.,predict(step.modelD,newdata=sampleX))]
  sampleX[,Bx := pmax(0.,predict(step.modelB,newdata=sampleX))]
  # sampleX[,Vx := pmax(0.,predict(step.modelV,newdata=sampleX))]
  sampleX[,Bpx := pmax(0.,predict(step.modelBp,newdata=sampleX))]
  sampleX[,Bspx := pmax(0.,predict(step.modelBsp,newdata=sampleX))]
  sampleX[,Bdx := pmax(0.,predict(step.modelBd,newdata=sampleX))]
  sampleX[,BpPerx := Bpx/(Bpx+Bspx+Bdx)*100]
  sampleX[,BspPerx := Bspx/(Bpx+Bspx+Bdx)*100]
  sampleX[,BdPerx := Bdx/(Bpx+Bspx+Bdx)*100]
  sampleX[,rootBAp:=BAp^0.5]
  
  nax <- unique(which(is.na(sampleX),arr.ind = T)[,1])
  if(length(nax)>0) sampleX <- sampleX[-nax]
  mux <- sampleX[,colMeans(cbind(Hx,Dx,Bx,BpPerx,BspPerx,BdPerx))]
  sigmax <- sampleX[,cov(cbind(Hx,Dx,Bx,BpPerx,BspPerx,BdPerx))]
  
  pMvnormx <- c(mux,as.vector(sigmax))

###sample from second measurement
  mux2 <- c(segIDx$h2 + mu2[3],segIDx$dbh2 + mu2[2],segIDx$BAtot2 + mu2[1],
            segIDx$BApPer2 + mu2[4],segIDx$BAspPer2 + mu2[5],segIDx$BAbPer2 + mu2[6])
  sigmax2 <- sigma2[c(3,2,1,4:6),c(3,2,1,4:6)]
  
  pMvnormx2 <- c(mux2,as.vector(sigmax2))

  
  LL <- solve(sigmax + sigmax2,tol=1e-20)
  sigmaPost <- sigmax %*% LL %*% sigmax2
  muPost <- sigmax2 %*% LL %*% mux + 
    sigmax %*% LL %*% mux2
  
  pMvnormPost <- c(muPost,as.vector(sigmaPost))
  # ss= inv(inv(sigmax)+inv(sigmax2))
  # ff <- sigmaPost %*% (sigmax %^%(-1)) %*% mux + sigmaPost %*% (sigmax2 %^%(-1)) %*% mux2

  # return(list(muPrior=mux,muLik=mux2,muPost=as.vector(muPost),
  #             sigPrior=sigmax,sigLik=sigmax2,sigPost=sigmaPost))
  pars <- c(segIDx$segID,as.vector(pMvnormx),as.vector(pMvnormx2),as.vector(pMvnormPost))
  return(pars)
}


###function for structural variables forcast and uncertainty 
prForUnc <- function(segIDx,nSample,yearUnc,tileX){
  if(yearUnc=="all"){
    muUnc <- errData$all$muFSVda
    sigmaUnc <- errData$all$sigmaFSVda
    step.probitX <- step.probit$all
  }else{
    muUnc <- errData[[paste0("y",yearUnc)]][[paste0("t",tileX)]]$muFSVda
    sigmaUnc <- errData[[paste0("y",yearUnc)]][[paste0("t",tileX)]]$sigmaFSVda
    step.probitX <- step.probit[[paste0("y",yearUnc)]][[paste0("t",tileX)]]
  }

  pST <- predict(step.probitX,type='p',segIDx)   ### needs to be changed . We need to calculate with 2016 and 2019 data
  
  st <- sample(rep(1:5,round(nSample*pST)),nSample,replace = T)
  set.seed(1234)
  sampleError <- data.table(mvrnorm(nSample*2,mu=muUnc,Sigma=sigmaUnc))
  
  # segIDx <- dataSurV[segID==2]
  sampleX <- data.table()
  sampleX$H <- segIDx$H + sampleError$H
  sampleX$D <- segIDx$D + sampleError$D
  sampleX$BAtot <- segIDx$BAtot + sampleError$G
  sampleX$BApPer <- segIDx$BApPer + sampleError$BAp
  sampleX$BAspPer <- segIDx$BAspPer + sampleError$BAsp
  sampleX$BAbPer <- segIDx$BAbPer + sampleError$BAb
  sampleX$BAp <- segIDx$BApPer * sampleX$BAtot
  sampleX$BAsp <- segIDx$BAspPer * sampleX$BAtot
  sampleX$BAb <- segIDx$BAbPer * sampleX$BAtot
  sampleX <- sampleX[H>1.5]
  sampleX <- sampleX[D>0.5]
  sampleX <- sampleX[BAtot>0.045]
  sampleX[,rootBAp:=sampleError$BAp^0.5]
  sampleX <- sampleX[1:min(nSample,nrow(sampleX))]
  if(nrow(sampleX)<nSample){
    sample1 <- sampleX
    set.seed(123)
    sampleError <- data.table(mvrnorm(nSample*2,mu=muUnc,Sigma=sigmaUnc))
    # segIDx <- dataSurV[segID==2]
    sampleX <- data.table()
    sampleX$H <- segIDx$H + sampleError$H
    sampleX$D <- segIDx$D + sampleError$D
    sampleX$BAtot <- segIDx$BAtot + sampleError$G
    sampleX$BApPer <- segIDx$BApPer + sampleError$BAp
    sampleX$BAspPer <- segIDx$BAspPer + sampleError$BAsp
    sampleX$BAbPer <- segIDx$BAbPer + sampleError$BAb
    sampleX$BAp <- segIDx$BApPer * sampleX$BAtot
    sampleX$BAsp <- segIDx$BAspPer * sampleX$BAtot
    sampleX$BAb <- segIDx$BAbPer * sampleX$BAtot
    sampleX <- sampleX[H>1.5]
    sampleX <- sampleX[D>0.5]
    sampleX <- sampleX[BAtot>0.045]
    sampleX <- rbind(sample1,sampleX)
    sampleX <- sampleX[1:min(nSample,nrow(sampleX))]
  }
  
  sampleX[, c("BApPer", "BAspPer", "BAbPer"):=
            as.list(fixBAper(unlist(.(BApPer,BAspPer,BAbPer)))), 
          by = seq_len(nrow(sampleX))]
  
  max.pro.est<-apply(segIDx[, c('BApPer','BAspPer','BAbPer')], 1, which.max)
  segIDx$max.pro.est=max.pro.est
  
  if(yearUnc=="all"){
    logistic.model <- logisticPureF$all
  }else{
    logistic.model <- logisticPureF[[paste0("y",yearUnc)]][[paste0("t",tileX)]]
  }
  
  set.seed(1234)
  sampleX$pureF <- runif(min(nSample,nrow(sampleX)),0,1)<predict(logistic.model,type="response",newdata = segIDx)
  if(max.pro.est==1) sampleX[which(pureF),c("BApPer","BAspPer","BAbPer"):=list(100,0,0)]
  if(max.pro.est==2) sampleX[which(pureF),c("BApPer","BAspPer","BAbPer"):=list(0,100,0)]
  if(max.pro.est==3) sampleX[which(pureF),c("BApPer","BAspPer","BAbPer"):=list(0,0,100)]
  
  sampleX[,BAp:=BApPer*BAtot/100]
  sampleX[,BAsp:=BAspPer*BAtot/100]
  sampleX[,BAb:=BAbPer*BAtot/100]

  sampleX[,segID:=segIDx$segID]
  
  sampleX[,BAtot:=(BAp+BAsp+BAb)]
  sampleX[,BAh:=BAtot*H]
  sampleX[,N:=BAtot/(pi*(D/200)^2)]
  b = -1.605 ###coefficient of Reineke
  sampleX[,SDI:=N *(D/10)^b]
  sampleX$st <- st
  
  sampleX$st <- factor(sampleX$st)
  sampleX[,Hx := pmax(0.,predict(step.modelH,newdata=sampleX))]
  sampleX[,Dx := pmax(0.,predict(step.modelD,newdata=sampleX))]
  sampleX[,Bx := pmax(0.,predict(step.modelB,newdata=sampleX))]
  sampleX[,Bpx := pmax(0.,predict(step.modelBp,newdata=sampleX))]
  sampleX[,Bspx := pmax(0.,predict(step.modelBsp,newdata=sampleX))]
  sampleX[,Bdx := pmax(0.,predict(step.modelBd,newdata=sampleX))]
  sampleX[,BpPerx := Bpx/(Bpx+Bspx+Bdx)*100]
  sampleX[,BspPerx := Bspx/(Bpx+Bspx+Bdx)*100]
  sampleX[,BdPerx := Bdx/(Bpx+Bspx+Bdx)*100]
  
  nax <- unique(which(is.na(sampleX),arr.ind = T)[,1])
  if(length(nax)>0) sampleX <- sampleX[-nax]
  mux <- sampleX[,colMeans(cbind(Hx,Dx,Bx,BpPerx,BspPerx,BdPerx))]
  sigmax <- sampleX[,cov(cbind(Hx,Dx,Bx,BpPerx,BspPerx,BdPerx))]
  
  pMvnormx <- c(mux,as.vector(sigmax))
  
  ###sample from second measurement
  pars <- as.vector(pMvnormx)
  return(pars)
}

