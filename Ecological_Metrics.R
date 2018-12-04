library("spatstat", lib.loc="~/R/win-library/3.3")
library("plyr", lib.loc="~/R/win-library/3.3")
library("doParallel", lib.loc="~/R/win-library/3.3")
library("rLiDAR", lib.loc="~/R/win-library/3.2")
library("aspace", lib.loc="~/R/win-library/3.3")
library("maptools", lib.loc="~/R/win-library/3.3")
library("moments", lib.loc="~/R/win-library/3.3")
library("raster", lib.loc="~/R/win-library/3.3")
library("tripack")
library("RANN")
library("rgeos")
rasterOptions(maxmemory=15e+08,chunksize=15e+07)
options(scipen=999)

enveloper=function(area){
  cl=makeCluster(39)
  registerDoParallel(cl)
  envelopes=foreach(j=2:25, .combine='rbind', .inorder=TRUE) %dopar% {
    library("spatstat", lib.loc="~/R/win-library/3.2")
    library("maptools", lib.loc="~/R/win-library/3.2")
    library("raster", lib.loc="~/R/win-library/3.2")
    sample_area=extent(0,sqrt(area),0,sqrt(area))
    sample_area=as(sample_area, 'SpatialPolygons')
    n=99
    Km = matrix(nrow=n, ncol=513)
    for(i in 1:n){
      Pmc = runifpoint(j,win=sample_area)    # Generate random points
      Kmc = Kest(Pmc, nsim=n, border="border")   # Compute K
      Km[i,] = Kmc$iso             
    }
    list(maxes=apply((Km),2,max), mins=apply((Km),2,min))
  }
  stopCluster(cl)
  return(envelopes)
}
envelopes=enveloper(400)

#taken from r raster v 1.0.0-1
reclasser=function(x){
  x1=x
  nc <- ncol(x1)
  nextclump <- 1
  v <- matrix(0, ncol=ncol(x1)+1, nrow=nrow(x1)+1)
  
  for (r in 1:nrow(x1)) {
    rr <- r + 1
    b <- getValues(x, r)
    b <- which(b != 0)
    for ( cc in b ) {
      vv <- na.omit(c(v[(rr-1),cc], v[rr, cc-1]))
      vm <- max(vv)
      if (vm > 0) {  
        v[rr, cc] <- vm
        vvv <- vv[(vv > 0) & (vv < vm)]
        if (length(vvv) > 0) {
          vvv <- unique(vvv)
          for (i in vvv) {
            v[v==i] <- vm
          }
        }
      } else {
        v[rr, cc] <- nextclump
        nextclump <- nextclump + 1					
      }
    }
  }
  
  v[v==0] <- NA
  rm(x)
  x1 <- setValues(x1, as.vector(t(v[-1,-ncol(v)])))
  u <- na.omit(unique(x1))
  u <- cbind(u, u, 1:length(u))
  x1 <- reclassify(x1, u)
  return(x1)
}

species_predictor=function(crown){
  if (leaves=='off'){
    leavesT=0
  }else{
    leavesT=1
  }
  #crown=basins==basin
  values(crown)[values(crown)==0]=NA
  #crown=tiff*crown
  heights=as.matrix(crown)
  heights <- heights[,colSums(is.na(heights))<nrow(heights)]
  if (is.matrix(heights) & sum(heights, na.rm=TRUE)>0){
    heights <- heights[rowSums(is.na(heights))<ncol(heights),]
  }
  if (is.matrix(heights) & sum(heights, na.rm=TRUE)>0){
    max_coords=which(heights == max(heights, na.rm=TRUE), arr.ind = TRUE)
    max_coords=max_coords[1,]
    east_west=heights[max_coords[1],]
    CV_EW=sd(na.omit(east_west))/mean(na.omit(east_west)) #####
    loc_max_EW=which(na.omit(east_west) == max(na.omit(east_west)))/length(na.omit(east_west)) #####
    loc_max_EW_tot=which(east_west == max(na.omit(east_west)))/length(east_west)  #####
    perc_max_EW=max(na.omit(east_west))/mean(na.omit(east_west))  ####
    x=1:length(na.omit(east_west))
    x2=x^2
    mod=lm(na.omit(east_west)~x+x2)
    coefx_EW=mod$coefficients[1]  ####
    coefx2_EW=mod$coefficients[2]  ####
    
    north_south=heights[,max_coords[2]]
    CV_NS=sd(na.omit(north_south))/mean(na.omit(north_south)) #####
    loc_max_NS=which(na.omit(north_south) == max(na.omit(north_south)))/length(na.omit(north_south)) #####
    loc_max_NS_tot=which(north_south == max(na.omit(north_south)))/length(north_south)  #####
    perc_max_NS=max(na.omit(north_south))/mean(na.omit(north_south))  ####
    x=1:length(na.omit(north_south))
    x2=x^2
    mod=lm(na.omit(north_south)~x+x2)
    coefx_NS=mod$coefficients[1]  ####
    coefx2_NS=mod$coefficients[2]  ####
    
    CV_allNS=sd(rowMeans(heights,na.rm=TRUE))/mean(rowMeans(heights,na.rm=TRUE))  ####
    CV_allEW=sd(colMeans(heights,na.rm=TRUE))/mean(colMeans(heights,na.rm=TRUE))  ####
  } else{
    CV_EW=NA  #####
    loc_max_EW=NA  #####
    loc_max_EW_tot=NA  #####
    perc_max_EW=NA   ####
    coefx_EW=NA  ####
    coefx2_EW=NA  ####
    CV_NS=NA  #####
    loc_max_NS=NA  #####
    loc_max_NS_tot=NA   #####
    perc_max_NS=NA   ####
    coefx_NS=NA   ####
    coefx2_NS=NA   ####
    CV_allNS=NA  ####
    CV_allEW=NA   ####
  }
  
  CVht=sd(na.omit(heights))/mean(heights, na.rm=TRUE)   #######
  max_pt=rasterToPoints(crown==max(heights, na.rm=TRUE))
  max_pt= matrix(max_pt[max_pt[,3]==1,],ncol=3, nrow=1)
  
  surrounding=raster::extract(crown, matrix(c(max_pt[,1],max_pt[,2]), ncol=2, nrow=1), buffer=1.25)
  ppmT=raster::extract(ppm, matrix(c(max_pt[,1],max_pt[,2]), ncol=2, nrow=1), buffer=1.25)
  if (is.na(ppmT)){
    ppmT=mean(getValues(ppm), na.rm=TRUE)
  }else{
    ppmT=mean(na.omit(ppmT[[1]]))}
  surrounding[[1]][is.na(surrounding[[1]])]=0
  diffs=max(surrounding[[1]])-surrounding[[1]]
  mean_difs_perc=mean(diffs[diffs!=0])/max(surrounding[[1]])   ##########
  SD_difs_perc=sd(diffs[diffs!=0])/max(surrounding[[1]])   #######
  
  surrounding=raster::extract(crown, matrix(c(max_pt[,1],max_pt[,2]), ncol=2, nrow=1), buffer=2.1)
  surrounding[[1]][is.na(surrounding[[1]])]=0
  diffs=max(surrounding[[1]])-surrounding[[1]]
  mean_difs2_perc=mean(diffs[diffs!=0])/max(surrounding[[1]])   ##########
  SD_difs2_perc=sd(diffs[diffs!=0])/max(surrounding[[1]])   #######
  
  circularity=(perimeter(tree)^2)/(4*3.14159*area(tree))    ######
  width=extent(tree)[2]-extent(tree)[1]
  height=extent(tree)[4]-extent(tree)[3]
  crown_width=width+height/2
  WidthToHeight=crown_width/max(surrounding[[1]])    ######
  
  pred=1.398347-7.572205*leavesT+0.16013*ppmT-1.936999*mean_difs_perc-1.651575*SD_difs2_perc-1.805987*WidthToHeight-
    2.512163*CV_EW-0.478505*loc_max_EW_tot-1.626035*loc_max_NS_tot+2.278141*perc_max_EW-0.095667*coefx_EW-0.079723*coefx2_EW-
    3.180399*CV_NS+ 2.538378*perc_max_NS-0.113139*coefx_NS-0.104072*coefx2_NS+0.448394*1*15+
    0.007175*coefx_EW*coefx2_EW+0.007507*coefx_NS*coefx2_NS
  pred=1/(1+(1/exp(pred)))
  pred=ifelse(pred>.69,1,0)
  pred=ifelse(is.na(pred),.7,pred)
  pred
}

lclmaxer=function(tiff, plotShape, plotShape_b){
  #First smooth CHM
  #smoothed=focal(tiff, w=matrix(1, nrow=3, ncol=3),fun=function(x){mean(x,na.rm=TRUE)})
  gf=matrix(c( 0.01134374, 0.08381951, 0.01134374, 0.08381951, 0.61934703, 0.08381951,0.01134374, 0.08381951, 0.01134374), nrow=3, ncol=3)
  smoothed=focal(tiff, w=gf, na.rm=TRUE, pad=TRUE)
  
  #######################################################
  #variable radius local maxima
  ##VARIABLE RADIUS DISABLED SINCE IT ONLY FILTERED OUT GOOD TREES 
  ##LOCAL MAXIMA DISABLED SINCE THERE'S NO POINT, WS DOES IT ALREADY  
  out=raster(smoothed)
  rastercells=getValues(smoothed)
  #calculate windows at different widths
  out <- apply(getValuesFocal(smoothed, 1, nrow(smoothed),ngb=3, array=TRUE),1,max)
  # v5 <- apply(getValuesFocal(smoothed, 1, nrow(smoothed), ngb=5, array=TRUE),1,max)
  # v7 <- apply(getValuesFocal(smoothed, 1, nrow(smoothed), ngb=7, array=TRUE),1,max)
  # v9 <- apply(getValuesFocal(smoothed, 1, nrow(smoothed), ngb=9, array=TRUE),1,max)
  # v11 <- apply(getValuesFocal(smoothed, 1, nrow(smoothed), ngb=11, array=TRUE),1,max)
  # v13 <- apply(getValuesFocal(smoothed, 1, nrow(smoothed), ngb=13, array=TRUE),1,max)
  # v15 <- apply(getValuesFocal(smoothed, 1, nrow(smoothed), ngb=15, array=TRUE),1,max)
  #estimate crown width from PEF data
  # rastercells[rastercells<0]=0
  # width=1.46056470+(rastercells*-0.56328814)+(log(rastercells)*2.05763111)+((rastercells^2)*0.01697101)
  # #width=.0108*(rastercells^2)-.1175*rastercells+3.534
  # width=width/.75
  ##use the appropriate window based on pixel hieght
  ##VARIABLE RADIUS DISABLED SINCE IT ONLY FILTERED OUT GOOD TREES (local vs regional width equations?, too much smoothing?)
  # out[]=sapply(1:length(width), simplify=TRUE,function(x){
  #   if (is.na(width[x])){
  #     width[x]
  #   }else if (width[x]<4){
  #     v3[x]
  #   }else if (width[x]<6){
  #     v5[x]
  #   }else if (width[x]<8){
  #     v7[x]
  #   }else if (width[x]<10){
  #     v9[x]
  #   }else if (width[x]<12){
  #     v11[x]
  #   }else if (width[x]<14){
  #     v13[x]
  #   }else if (width[x]>=14){
  #     v15[x]
  #   }
  # })
  #get local maxima where the focal pixel equals the smoothed raster, multiply by tiff for actual height
  trueLM<-(smoothed==out)*tiff
  trueLM[trueLM==0]<-NA
  trueLM[trueLM<1]<-NA
  points=rasterToPoints(trueLM)
  ########################################################
  
  masked_tiff=raster::mask(tiff,plotShape)
  #General Metrics
  all_values=getValues(masked_tiff)
  all_values=all_values[!is.na(all_values)]
  total_mean=mean(all_values)  #####COVARIATE
  total_median=median(all_values)  #####COVARIATE
  total_sd=sd(all_values)  #####COVARIATE
  total_range=diff(range(all_values))  #####COVARIATE
  total_skew=skewness(all_values)  #####COVARIATE
  total_kurt=kurtosis(all_values)  #####COVARIATE
  Perc_Above_P20=length(all_values[all_values>.2*max(all_values)])/length(all_values)  #####COVARIATE
  Perc_Above_P40=length(all_values[all_values>.4*max(all_values)])/length(all_values)  #####COVARIATE
  Perc_Above_P60=length(all_values[all_values>.6*max(all_values)])/length(all_values)  #####COVARIATE
  Perc_Above_P80=length(all_values[all_values>.8*max(all_values)])/length(all_values)  #####COVARIATE
  Perc_Above_P95=length(all_values[all_values>.95*max(all_values)])/length(all_values)  #####COVARIATE
  Perc_Above_5m=length(all_values[all_values>5])/length(all_values)  #####COVARIATE
  Perc_Above_10m=length(all_values[all_values>10])/length(all_values)  #####COVARIATE
  Perc_Above_15m=length(all_values[all_values>15])/length(all_values)  #####COVARIATE
  Perc_Above_20m=length(all_values[all_values>20])/length(all_values)  #####COVARIATE
  Perc_Above_25m=length(all_values[all_values>25])/length(all_values)  #####COVARIATE
  Perc_Above_35m=length(all_values[all_values>30])/length(all_values)  #####COVARIATE
  Percentile_20=quantile(all_values, .2)  #####COVARIATE
  Percentile_40=quantile(all_values, .4)  #####COVARIATE
  Percentile_60=quantile(all_values, .6)  #####COVARIATE
  Percentile_80=quantile(all_values, .8)  #####COVARIATE
  Percentile_95=quantile(all_values, .95)  #####COVARIATE
  
  
  ######################################Watershed ALGORITHM, SORT, FILL####################################### 
  #points_notbuffed=rasterToPoints(raster::mask(trueLM, plotShape))
  Ras_Size=sqrt(length(smoothed))
  cellz=sum(stack((smoothed*-1), ((trueLM*0)+100)*-1), na.rm=TRUE)
  cellz=getValues(cellz)
  cellz[cellz>-3]=NA
  cells=data.frame(cellz=as.vector(cellz))
  id=rownames(cells)
  cells=cbind(id=id,cells)
  cells[(cells$cellz<=-100 & !is.na(cells$cellz)),]$cellz=-100
  cellsS=arrange(cells, cellz)
  cells$basin=NA
  cells$id=as.numeric(as.character(cells$id))
  cellsS$id=as.numeric(as.character(cellsS$id))
  basin_num=1
  
  for (row in 1:nrow(cells)){
    pixel=cellsS[row,]
    if (is.na(pixel$cellz)){
      cells[pixel$id,]$basin=NA
      next
    }
    if (as.numeric(pixel$id)%%Ras_Size<=1 & as.numeric(pixel$id)<Ras_Size & as.numeric(pixel$id) >(Ras_Size^2)-Ras_Size){
      cells[pixel$id,]$basin=NA
    }else{
      #eight neighbot search, queen
      neighbors=c(cells[as.numeric(pixel$id)-(Ras_Size+1),]$basin, cells[as.numeric(pixel$id)-Ras_Size,]$basin, 
                  cells[as.numeric(pixel$id)-(Ras_Size-1),]$basin,cells[as.numeric(pixel$id)-1,]$basin,
                  cells[as.numeric(pixel$id)+1,]$basin, cells[as.numeric(pixel$id)+(Ras_Size-1),]$basin,
                  cells[as.numeric(pixel$id)+Ras_Size,]$basin, cells[as.numeric(pixel$id)+(Ras_Size+1),]$basin)
      #four neighbor search, rook, really doesn't work well
      #      neighbors=c(cells[as.numeric(pixel$id)-Ras_Size,]$basin, 
      #                  cells[as.numeric(pixel$id)-1,]$basin,
      #                  cells[as.numeric(pixel$id)+1,]$basin, 
      #                  cells[as.numeric(pixel$id)+Ras_Size,]$basin)      
      if (sum(neighbors, na.rm=TRUE)==0){
        cells[pixel$id,]$basin=basin_num
        basin_num=basin_num+1
      }else{
        if(length(unique(na.omit(neighbors)))>1 ){
          cells[pixel$id,]$basin=NA
        }else{
          cells[pixel$id,]$basin=as.numeric(names(sort(-table(neighbors)))[1])
        }
      }
    }
  }
  
  #DEVELOP WATERSHEDS INTO RASTER
  basins=raster(smoothed)
  basins[]=cells$basin
  
  #CLEAN WATERSHEDS UP
  #small_basins=subset(freq(basins), freq(basins)[,2]<4)[,1]
  #small_basins=append(small_basins, 0)
  #nas=data.frame(small_basins)
  #nas[,"nas"]=NA
  #basins=reclassify(basins, nas)
  #Clean up basins  
  basinvalues=getValues(basins)
  basinlist=na.omit(unique(basinvalues))
  #all values within a wateshed less than 1/2rds the height of the highest point are removed
  bads=lapply(basinlist,function(x){  
    indicies=which(basinvalues == x)
    heights=(tiff[indicies])
    indicies[heights<=max(na.omit(heights))*.5]})
  bads=unlist(bads)
  basins[bads]=NA
  #reclassify to isolate trees further segmented by the above proceedure
  basins=reclasser(basins)
  #remove lone cells surrounded by NA, may eliminate small trees
  #basins=focal(basins, w=matrix(c(NA,1,NA,1,1,1,NA,1,NA),ncol=3, nrow=3),fun=function(x){
  # if(sum(is.na(x))>7){
  #   NA
  # }
  # else{
  #    x[5]
  # }
  # })
  
  #Canopy Counts and Sizes WITH BUFFER
  basins_NB=raster::mask(basins,plotShape)
  if (length(na.omit(freq(basins_NB))[,2])>0){
    in_trees_B=na.omit(unique(getValues(basins_NB)))
    tree_sizes=data.frame(na.omit(freq(basins)))
    in_tree_sizes=tree_sizes[match(in_trees_B,tree_sizes[,1]),]
    mean_crown_area_B=mean(na.omit(in_tree_sizes)[,2])*.5625 #####COVARIATE
    SD_crown_area_B=sd(na.omit(in_tree_sizes)[,2])*.5625 #####COVARIATE
    max_crown_area_B=max(na.omit(in_tree_sizes)[,2])*.5625 #####COVARIATE
    watershed_count_inplot=nrow(na.omit(in_tree_sizes)) #####COVARIATE
    
    #EACH TREES HEIGHT MEASUREMENTS
    max_hts=NULL
    avg_hts=NULL
    points_notbuffed=NULL
    #range_hts=NULL
    #thirdQ_hts=NULL
    crown_volumesT=NULL
    crown_volumesIN=NULL
    mean_resids=NULL
    pos_resids=NULL
    pos_resids_above=NULL
    species=NULL
    diametersT=NULL
    diameters=NULL
    biomassT=NULL
    biomass=NULL
    perc_in_tally=NULL
    
    for (basin in in_trees_B){
      crown=basins==basin
      perc_in=sum(na.omit(getValues(raster::mask(crown,plotShape))))/sum(na.omit(getValues(crown)))
      perc_in_tally=append(perc_in_tally, perc_in)
      
      crown=crown*tiff
      SPpred=species_predictor(crown)
      species=append(species, SPpred)
      
      
      topH=maxValue(crown)
      max_loc=xyFromCell(crown,which.max(crown))
      
      heights=getValues(crown)
      heights=na.omit(heights[heights>0])
      if (length(heights)==0){
        next
      }  
      
      #my own model, diam-height
      if (basin %in% in_trees_B){
        #HW
        if (SPpred==0){
          DBH=13.7426-0.5355*topH+.0518*topH^2
        }else{#SW
          DBH=9.076091-0.003484*topH+ 0.054784*topH^2
        }
      }
      diametersT=append(diametersT, DBH)
      diameters=append(diameters, DBH*perc_in)
      
      #jenkins generic hardwood/softwood equations
      if (SPpred==1){
        mass=0.5 + ((15000 * (DBH ^ 2.7))/((DBH ^ 2.7) + 364946))
      }else{
        mass=0.5 + ((25000 * (DBH ^ 2.5))/((DBH ^ 2.5) + 246872))}
      biomassT=append(biomassT, mass)
      biomass=append(biomass, mass*perc_in)
      
      total_area=length(heights)
      residuals=NULL
      for (crossec in trunc((min(heights)-1),0):topH){
        profile=length(heights[heights>crossec])/length(heights)
        prof_ht=abs(crossec/topH)
        #Curtis equation from 3100 howland trees v v v v
        residuals=append(residuals, profile-(.9626*(exp(-2.5369*prof_ht^6.9433))))
      }
      
      max_hts=append(max_hts, topH)
      points_notbuffed=rbind(points_notbuffed, max_loc)
      avg_hts=append(avg_hts, mean(heights))
      #range_hts=append(range_hts, range(heights)[2]-range(heights)[1])
      #thirdQ_hts=append(thirdQ_hts, quantile(heights, .75))
      heights=heights-(.66*max(heights))
      heights[heights<0] <- 0
      crown_volumesT=append(crown_volumesT,sum(heights*.75*.75))
      crown_volumesIN=append(crown_volumesIN,sum(heights*.75*.75)*perc_in)
      mean_resids=append(mean_resids,round(mean(residuals),3))
      pos_resids=append(pos_resids,sum(residuals>0)/length(residuals))
      pos_resids_above=append(pos_resids_above,sum(residuals[length(residuals)/2:length(residuals)]>0)/length(residuals[length(residuals)/2:length(residuals)]))
      
    }
    points_notbuffed=cbind(points_notbuffed,extract(tiff,points_notbuffed))
    
    ppp_notbuffed=ppp(points_notbuffed[,1],points_notbuffed[,2], window=as.owin(as(plotShape, 'SpatialPolygons')))
    ppp_buffed=ppp(points_notbuffed[,1],points_notbuffed[,2], window=as.owin(as(plotShape_b, 'SpatialPolygons')))
    in_trees_NB=extract(basins,matrix(points_notbuffed[,1:2], ncol=2))
    in_tree_sizes_NB=tree_sizes[match(in_trees_NB,tree_sizes[,1]),]
    mean_crown_area_NB=mean(na.omit(in_tree_sizes_NB)[,2])*.5625 #####COVARIATE
    SD_crown_area_NB=if (nrow(na.omit(freq(basins_NB)))>1){sd(na.omit(in_tree_sizes_NB)[,2])*.5625 }else{0} #####COVARIATE
    max_crown_area_NB=max(na.omit(in_tree_sizes_NB)[,2])*.5625 #####COVARIATE
    
    #encroaching tree info
    encroaching_tree_count=watershed_count_inplot-ppp_notbuffed$n #####COVARIATE
    tree_sizes=data.frame(na.omit(freq(basins_NB)))
    
    rejects=cbind(attr(ppp_notbuffed, "rejects")$x,attr(ppp_notbuffed, "rejects")$y)
    if (!is.null(rejects)){
      encroaching_trees=extract(basins,cbind(attr(ppp_notbuffed, "rejects")$x,attr(ppp_notbuffed, "rejects")$y))
      encroaching_tree_sizes=tree_sizes[tree_sizes$value %in% encroaching_trees,2]
      encroach_area=sum(na.omit(encroaching_tree_sizes)) #####COVARIATE
      encroach_mean=mean(na.omit(encroaching_tree_sizes)) #####COVARIATE
    }else{
      encroaching_trees=0
      encroach_area=0
      encroach_mean=0
    }
    
    canopy_cover_NB=1-mean(is.na(extract(basins, plotShape)[[1]])) #####COVARIATE
    
    perc_in_tally=sum(na.omit(perc_in_tally))   #####COVARIATE
    biomass_in=sum(na.omit(biomass))   #####COVARIATE
    biomass_tot=sum(na.omit(biomassT))   #####COVARIATE
    crown_volume_tot=sum(na.omit(crown_volumesT))   #####COVARIATE
    crown_volume_in=sum(na.omit(crown_volumesIN))   #####COVARIATE
    basal_area_tot=sum(na.omit(diametersT*diametersT*.00007854))   #####COVARIATE
    basal_area_in=sum(na.omit(diameters*diameters*.00007854))   #####COVARIATE
    diameters=na.omit(diameters)
    QMD=sqrt(sum(diametersT^2)/length(diameters))   #####COVARIATE
    softwoodyness=mean(na.omit(species))   #####COVARIATE
    mean_crown_volumes=mean(crown_volumesT)   #####COVARIATE
    SD_crown_volumes=sd(crown_volumesT)   #####COVARIATE
    mean_top_ht=mean(na.omit(max_hts))   #####COVARIATE
    tot_top_ht=sum(na.omit(max_hts))   #####COVARIATE
    mean_ht=mean(na.omit(avg_hts))   #####COVARIATE
    
    #examine the residuals of this tree's crown profile compared to an average howland tree
    mean_residuals=mean(mean_resids)   #####COVARIATE
    sd_residuals=sd(mean_resids)   #####COVARIATE
    mean_positivity=mean(pos_resids)   #####COVARIATE
    sd_positivity=sd(pos_resids)   #####COVARIATE
    mean_top_positivity=mean(pos_resids_above)   #####COVARIATE
    
    #kurt_crown_volumes=kurtosis(crown_volumes)   #####COVARIATE
    #skew_crown_volumes=skewness(crown_volumes)   #####COVARIATE
  } else{
    perc_in_tally=0   
    biomass_in=0   
    biomass_tot=0   
    crown_volume_tot=0  
    crown_volume_in=0  
    basal_area_tot=0  
    basal_area_in=0  
    diameters=0  
    QMD=0   
    softwoodyness=0   
    mean_crown_volumes=0   
    SD_crown_volumes=0   
    mean_top_ht=0  
    tot_top_ht=0   
    mean_ht=0  
    mean_residuals=0  
    sd_residuals=0   
    mean_positivity=0   
    sd_positivity=0   
    mean_top_positivity=0   
    mean_crown_area_B=0 
    SD_crown_area_B=0 
    max_crown_area_B=0 
    watershed_count_inplot=0 
    mean_crown_area_NB=0 
    SD_crown_area_NB=0 
    max_crown_area_NB=0 
    encroaching_tree_count=0 
    encroach_area=0
    encroach_mean=0 
    canopy_cover_NB=0   
    #more optional but probably meaningless metrics
    #kurt_crown_volumes=0   #####COVARIATE
    #skew_crown_volumes=0   #####COVARIATE
    #skew_crown_area_B=0 #####COVARIATE
    #kurtosis_crown_area_B=0 #####COVARIATE
    #skew_crown_area_NB=0 #####COVARIATE
    #kurtosis_crown_area_NB=0 #####COVARIATE
    #SD_max_ht=0   #####COVARIATE
    #kurt_max_ht=0   #####COVARIATE
    #SD_mean_ht=0   #####COVARIATE
    #mean_range_ht=0  #####COVARIATE
    #SD_range_ht=0   #####COVARIATE
    #kurt_range_ht=0   #####COVARIATE
    #mean_thirdQ_hts=0   #####COVARIATE
    #SD_thirdQ_hts=0   #####COVARIATE
    #med_range_ht=0   #####COVARIATE
  }
  ###############################END OF WATERSHED STUFF##########################################  
  
  #nn dists
  
  nndists_NB=nndist(ppp_notbuffed)
  nndists=nndist(ppp_buffed)
  if (length(nndists_NB)<2){
    mean_dists_NB=0 #####COVARIATE
    SD_dists_NB=0  #####COVARIATE
    #med_dists_NB=0  #####COVARIATE
    skew_dists_NB=0  #####COVARIATE
    kurt_dists_NB=0  #####COVARIATE
  }else{
    mean_dists_NB=mean(nndists_NB) #####COVARIATE
    SD_dists_NB=sd(nndists_NB)  #####COVARIATE
    #med_dists_NB=median(nndists_NB)  #####COVARIATE
    skew_dists_NB=skewness(nndists_NB)  #####COVARIATE
    kurt_dists_NB=kurtosis(nndists_NB)  #####COVARIATE
  }
  if (length(nndists)<2){
    mean_dists_B=0 #####COVARIATE
    SD_dists_B=0  #####COVARIATE
    #med_dists_B=0  #####COVARIATE
    skew_dists_B=0  #####COVARIATE
    kurt_dists_B=0  #####COVARIATE
  }else{
    mean_dists_B=mean(nndists)  #####COVARIATE
    SD_dists_B=sd(nndists)  #####COVARIATE
    #med_dists_B=median(nndists)  #####COVARIATE
    skew_dists_B=skewness(nndists)  #####COVARIATE
    kurt_dists_B=kurtosis(nndists)  #####COVARIATE
  }
  
  #Bon Voroni Growing Space
  if (ppp_notbuffed$n>0){
    all_areas=dirichletWeights(ppp_buffed, exact=FALSE)
    all_areasNB=dirichletWeights(ppp_notbuffed, exact=FALSE)
    mean_growing_space=mean(areas_notbuffed)  #####COVARIATE
    sd_growing_space=sd(areas_notbuffed)   #####COVARIATE
    max_growing_space=max(areas_notbuffed)   #####COVARIATE
    #skew_growing_space=skewness(areas_notbuffed)   #####COVARIATE
    #kurt_growing_space=kurtosis(areas_notbuffed)   #####COVARIATE
    mean_growing_spaceNB=mean(all_areasNB)  #####COVARIATE
    sd_growing_spaceNB=sd(all_areasNB)   #####COVARIATE    
  }else{
    mean_growing_space=400 
    sd_growing_space=0
    #skew_growing_space=0  
    #kurt_growing_space=0 
    max_growing_space=400
    mean_growing_spaceNB=100
    sd_growing_spaceNB=0
  }
  
  ######################################Crown Light Availibility#############################################
  obscuration=NULL
  shadowed=NULL
  neighbor_tm_dists=NULL
  sky_view_areas=NULL
  crown_widths=data.frame(na.omit(freq(basins)))
  crown_widths$width=sqrt((crown_widths$count*.5*.5)/3.1415)*2
  
  #if there are trees in plot..
  if (nrow(points_notbuffed)>1){
    #for each tree in plot..
    for (point in 1:nrow(points_notbuffed)){
      #From 2/3 up the focus tree
      focal_tree=points_notbuffed[point,]
      ht=focal_tree[3]*.666
      #Calculate distances and azis of all other trees in buffered plot
      dists=apply(points[,1:2],1, function(x) sqrt((x[1]-focal_tree[1])^2+(x[2]-focal_tree[2])^2))
      azis=apply(points[,1:2],1, function(x){
        dx = x[1]-focal_tree[1]
        dy = x[2]-focal_tree[2]
        if (dx>0){
          90 - atan(dy/dx)*(180/3.14159)
        }else if (dx<0){
          270 - atan(dy/dx)*(180/3.14159)
        }else if (dx==0){
          360
        }
      }) 
      
      #360+(atan((x[1]-focal_tree[1])/(x[2]-focal_tree[2])))*(180/3.14159))
      
      #subset neighbors which are taller than 2/3 focal tree
      nearby_tree_locs=data.frame(points, dists,azis)
      nearby_tree_locs=subset(nearby_tree_locs, layer>ht)
      #order by distance, remove all trees with distances greater than 10m
      nearby_tree_locs <- nearby_tree_locs[order(nearby_tree_locs$dists),] 
      nearby_tree_locs <- nearby_tree_locs[nearby_tree_locs$dists<10,] 
      #remove the first tree, since that's the focal
      nearby_tree_locs <- nearby_tree_locs[-1,]
      nearby_tree_locs <- na.omit(nearby_tree_locs)
      
      #badly estimated crown width
      obscured_space=NULL
      triangle_space=NULL
      #If there are any neighboring trees at all..
      if (nrow(nearby_tree_locs)>0){
        #for each neighbor..
        for (neighbor in 1:nrow(nearby_tree_locs)){
          #extract which basin represents that tree
          basin_num=extract(basins,matrix(c(nearby_tree_locs$x[neighbor],nearby_tree_locs$y[neighbor]),ncol=2, nrow=1))
          
          #if the point happens to fall outside a basin, assume the tree's a red spruce and do some crappy crown width estimate
          if (is.na(basin_num)){
            bad_dbh=log(1-(((nearby_tree_locs$layer[neighbor]-1.35)/23.13)**1.038))/-.038
            if (is.na(bad_dbh)){
              bad_dbh=50
            }
            lcw=(2.1383*bad_dbh**.47933)/(2.32*bad_dbh**-.02993)
          }else{
            #otherwise get diameter based on measured watershed
            lcw=crown_widths[crown_widths$value==basin_num,3]
          }
          
          #Calculate degrees obscured by neighbor
          obsc=2*atan(lcw/(nearby_tree_locs$dists[neighbor]))*(180/3.14159)
          
          #from 2/3 up neighboring tree, construct a triangle, base the size of the crown, height the height of tree
          #$layer is the same as tree height
          ht_live_crown=nearby_tree_locs$layer[neighbor]*.66
          triangle=((nearby_tree_locs$layer[neighbor]-ht_live_crown)*lcw)/2
          
          #???
          shadow_size=seq(from=round(nearby_tree_locs$azis[neighbor]-.5*round(obsc)),to=round(nearby_tree_locs$azis[neighbor]+.5*round(obsc)))
          shadow_size=sapply(shadow_size, function(x) if(x>360){x-360}else{x})
          shadow_size=unique(sapply(shadow_size, function(x) if(x<0){360+x}else{x}))
          
          #if less than half of the tree's area isn't already shadowing the primary, add that shadow area
          if (sum(is.element(shadow_size, obscured_space))<length(shadow_size)/2){
            triangle_space=append(triangle_space,triangle)
          }
          
          obscured_space=append(obscured_space,shadow_size)
          obscured_space=unique(obscured_space) 
          percent_obscured=length(obscured_space)/360  
          
        }
      }else{triangle_space=0;percent_obscured=0}
      
      shadowed=append(shadowed, sum(triangle_space))
      obscuration=append(obscuration, percent_obscured)
      
      #########Upsidedown pyrimid stuff#############
      tree_center=data.frame(points_notbuffed)
      tree_center=data.frame(t(points_notbuffed[point,1:2]), z=points_notbuffed[point,3]*.66)
      nearby_tree_locs=nearby_tree_locs[order(nearby_tree_locs$azis),]
      #Returns vector slope
      #vectors=data.frame(t(mapply(function(x,y,z){ c(x-tree_center[1], y-tree_center[2], z-tree_center[3], azi)}, nearby_tree_locs$x,nearby_tree_locs$y, nearby_tree_locs$layer)))
      vectors=data.frame(t(mapply(function(x,y,z){ c(x-tree_center[1], y-tree_center[2], z-tree_center[3])}, nearby_tree_locs$x,nearby_tree_locs$y, nearby_tree_locs$layer)))
      neighbor_top_dists=mapply(function(x,y,z){sqrt(x*x + y*y + z*z)}, vectors$x, vectors$y, vectors$z)
      #possible covariate, distance to top..
      if (length(neighbor_top_dists)>0){
        neighbor_tm_dists=append(neighbor_tm_dists, mean(neighbor_top_dists))
        
        vectors$Ddist=mapply(function(z){50/z},vectors$z)
        Dvalues=data.frame(t(mapply(function(x,y,z,dv){c(x*dv, y*dv, z*dv)}, vectors$x, vectors$y, vectors$z, vectors$Ddist)))
        vertexes=data.frame(x=nearby_tree_locs$x+Dvalues[,1],y=nearby_tree_locs$y+Dvalues[,2], azi=nearby_tree_locs$azis )
        for (edge in 1:16){
          azi=(22.5*edge)*(3.14159/180)
          Dy=284*cos(azi)
          Dx=284*sin(azi)
          azi=(22.5*edge)
          coords=c(focal_tree[1]+Dx, focal_tree[2]+Dy,azi)
          vertexes=rbind(vertexes, coords)
        }
        vertexes=vertexes[order(vertexes$azi),]
        vertexes_dists=mapply(function(x,y){sqrt((focal_tree[1]-x)^2 + (focal_tree[2]-y)^2 )}, vertexes$x, vertexes$y)
        vertexes=cbind(vertexes,vertexes_dists)
        vertexes=vertexes[vertexes$vertexes_dists<285,]
        sp_poly <- SpatialPolygons(list(Polygons(list(Polygon(vertexes[,1:2])), ID=1)))
        sky_view_areas=append(sky_view_areas, sp_poly@polygons[[1]]@area)
      }else{sky_view_areas=append(sky_view_areas, 253380)
      }
      
      
    }
  }
  
  if (length(obscuration)>0){
    mean_obscured=mean(obscuration)  #####COVARIATE
    sd_obscured=sd(obscuration)  #####COVARIATE
    mean_shadowed=mean(shadowed)  #####COVARIATE
    sd_shadowed=sd(shadowed)  #####COVARIATE
    kurt_shadowed=kurtosis(shadowed)  #####COVARIATE
    sky_view_areas_m=mean(sky_view_areas)  #####COVARIATE
    sky_view_areas_sd=sd(sky_view_areas)  #####COVARIATE
    #sky_view_areas_kurt=kurtosis(sky_view_areas)  #####COVARIATE
    neighbor_top_dist_m=if (is.null(neighbor_tm_dists)){0}else{mean(neighbor_tm_dists)}  #####COVARIATE
    neighbor_top_dist_sd=if (is.null(neighbor_tm_dists)){0}else{sd(neighbor_tm_dists)}  #####COVARIATE
  }else{
    mean_obscured=0  #####COVARIATE
    sd_obscured=0  #####COVARIATE
    sd_shadowed=0  #####COVARIATE
    kurt_shadowed=0  #####COVARIATE
    mean_shadowed=0  #####COVARIATE
    sky_view_areas_m=253380  #####COVARIATE
    sky_view_areas_sd=0  #####COVARIATE
    #sky_view_areas_kurt=0  #####COVARIATE
    neighbor_top_dist_m=45  #####COVARIATE (45 is the greatest diagonal that can occur with a 20x40m triangle
    neighbor_top_dist_sd=0  #####COVARIATE
  }  
  
  ###################################### END OF Crown Light Availibility#############################################
  
  
  #pointiness, mean and standard deviation
  if (nrow(points_notbuffed)>1){  
    #get cells arund points)
    arounds=data.frame(extract(tiff,points_notbuffed[,1:2], buffer=1.25, cellnumbers=TRUE,df=TRUE, na.rm=TRUE))
    if (nrow(arounds)%%9>0){
      arounds=t(unstack(arounds, layer~ID))
      arounds=ldply(arounds, rbind)
    }else{
      arounds=t(unstack(arounds, layer~ID))
    }
    pointyness=t(apply(arounds, 1, function(x){
      if(is.list(x)){
        c(0, 0, 0, 0)
      }else{
        diffs=max(x)-x
        perc_difs=diffs/x
        c(mean(diffs), sd(diffs), mean(perc_difs), sd(perc_difs))
      }
    }))
    mean_raw_pointy=mean(na.omit(pointyness[,1])) #####COVARIATE
    sd_raw_pointy=mean(na.omit(pointyness[,2])) #####COVARIATE
    mean_pointy=mean(na.omit(pointyness[,3])) #####COVARIATE
    sd_pointy=mean(na.omit(pointyness[,4])) #####COVARIATE
    
  }else if (nrow(points_notbuffed)==1){
    arounds=data.frame(extract(tiff,matrix(points_notbuffed[1,1:2], ncol=2), buffer=1.25, cellnumbers=TRUE, na.rm=TRUE))
    surrounding=na.omit(arounds$value)
    diffs=max(surrounding)-surrounding
    perc_difs=diffs/surrounding
    mean_raw_pointy=mean(na.omit(diffs))
    sd_raw_pointy=sd(na.omit(diffs))
    mean_pointy=mean(na.omit(perc_difs))
    sd_pointy=sd(na.omit(perc_difs))
  }else{
    mean_raw_pointy=0
    sd_raw_pointy=0
    mean_pointy=0
    sd_pointy=0
  }
  ###########################################################################################  
  #Clustered or not..
  kest=Kest(ppp_buffed, correction="border")
  kest$border[is.nan(kest$border)]<-0
  if (ppp_buffed$n==1){csr=1
  }else if(ppp_buffed$n>25){csr=0
  }else {
    env=envelopes[ppp_buffed$n-1,]
    #some below envelope min, none above
    if (sum((kest$border-env$mins)<0) > 0 && sum((env$maxes-kest$border)<0) < 1){
      #evenly spaced
      csr=1
    }
    #some above the envelope max, none below
    if (sum((kest$border-env$mins)<0) < 1 && sum((env$maxes-kest$border)<0) > 0){
      #clustered
      csr=2
    }
    #some above and below the envelope
    if (sum((kest$border-env$mins)<0) > 0 && sum((env$maxes-kest$border)<0) > 0){
      #clustered and spaced at different scales
      csr=3
    }
    #if observed were all above lower envelope, and all less than the upper envelope
    if (sum((kest$border-env$mins)<0) == 0 && sum((env$maxes-kest$border)<0) == 0){
      #randomly distributed
      csr=0
    }
  }
  ####csr IS A COVARIATE
  
  ##########################################################################################
  ####Edge detection####, Lapacian transform
  Laplacian=matrix(-1, 3, 3)
  Laplacian[2, 2] <- 8
  Laplacian <- focal(smoothed,w=Laplacian, pad=TRUE, padValue=NA)
  Laplacian=raster::mask(Laplacian, plotShape)
  
  ####Edge detection####, Sobel transform
  A <- focal(smoothed,w=matrix(c(-1,0,1,-2,0,2,-1,0,1),nrow=3,ncol=3), pad=TRUE, padValue=NA)
  B <- focal(smoothed,w=matrix(c(-1,-2,-1,0,0,0,1,2,1),nrow=3,ncol=3), pad=TRUE, padValue=NA)
  Sobel=abs(A)+abs(B)
  Sobel=raster::mask(Sobel, plotShape)
  
  Sobel[Sobel<maxValue(Sobel)/3] <-NA
  rugosityS=na.omit(getValues(Sobel))
  rugoseS_count=length(rugosityS)  #####COVARIATE
  rugoseS_mean= mean(rugosityS)    #####COVARIATE
  rugoseS_sd= sd(rugosityS)    #####COVARIATE
  rugosityL=na.omit(getValues(Laplacian))
  rugoseL_mean= mean(rugosityL)    #####COVARIATE
  rugoseL_sd= sd(rugosityL)    #####COVARIATE
  
  #distance to edge, Sobel transform
  edge_pts=rasterToPoints(Sobel)
  if (nrow(edge_pts)>0 && nrow(points_notbuffed)>0){
    if (nrow(points_notbuffed)==1){pt2edge= apply(t(as.matrix(points[,1:2])), 1, function(center){min(distances(center, edge_pts[,1:2]))})}
    else{pt2edge= apply(points_notbuffed[,1:2], 1, function(center){min(distances(center, edge_pts[,1:2]))})}
    pt2Sedge_mean=mean(pt2edge)  #####COVARIATE
    pt2Sedge_max=max(pt2edge)  #####COVARIATE
  }else{ 
    pt2edge_mean=0
    pt2edge_max=0
  }
  
  
  #centriod distances
  if (watershed_count_inplot>0){
    #get centroids of clusters (lifted from github
    clump_id <- getValues(basins) 
    xy <- xyFromCell(basins,1:ncell(basins))
    centroids <- data.frame(xy, clump_id, is_clump = basins[] %in% freq(basins, useNA = 'no')[,1])
    centroids <- ddply(centroids[centroids$is_clump == T, ], .(clump_id), summarise, xm = mean(x), ym = mean(y))
    #get distances
    pt2cent= apply(matrix(points_notbuffed[,1:2],ncol=2), 1, function(center){min(distances(center, centroids[,2:3]))})
    mean_off_cent=mean(na.omit(pt2cent))  #####COVARIATE
    sd_off_cent=sd(na.omit(pt2cent))  #####COVARIATE
  }else{
    mean_off_cent=0
    sd_off_cent=0
  }
  
  
  basins=disaggregate(basins, 4)
  edges=boundaries(basins,classes=TRUE, type='inner', directions=4, asNA=TRUE)
  edge_pts=rasterToPoints(edges)
  if (nrow(edge_pts)>0 && nrow(points_notbuffed)>0){
    if (nrow(points_notbuffed)==1){pt2edge= apply(t(as.matrix(points_notbuffed[,1:2])), 1, function(center){min(distances(center, edge_pts[,1:2]))})}
    else{pt2edge= apply(points_notbuffed[,1:2], 1, function(center){min(distances(center, edge_pts[,1:2]))})}
    pt2WSedge_mean=mean(pt2edge)  #####COVARIATE
    pt2WSedge_sd=sd(pt2edge)  #####COVARIATE
    pt2WSedge_max=max(pt2edge)  #####COVARIATE
  }else{ 
    pt2WSedge_mean=0
    pt2WSedge_sd=0
    pt2WSedge_max=0
  }
  #percent occupied pixels based on # of pixels with a 20x20m .75m origional raster
  rugose_WS_perc=length(na.omit(getValues(raster::mask(edges,plotShape))))/11377  #####COVARIATE
  
  ##########################################################################################################  
  ####Get heights####
  LMheights=points_notbuffed[,3]
  #Shawn's competition indicies
  if (length(LMheights)>1){
    #Simpson index with 2m height bins
    x=3
    H_Classes=round(LMheights/x)*x
    frequencies=table(H_Classes)
    SIM_1=frequencies*(frequencies-1)
    SIMH =1 - (sum(SIM_1)/ (length(LMheights)*(length(LMheights)-1)))   #####COVARIATE
    #GINI based on Heights
    n=length(LMheights)
    mu=mean(LMheights)
    N=n * (n - 1)
    ox <- LMheights[order(LMheights)]
    dsum <- drop(crossprod(2 * 1:n - n - 1,  ox))
    GINIH = dsum / (mu * N)   #####COVARIATE
    #TD Diameter Differentation Index
    biggers=lapply(LMheights,function(x){x/LMheights})
    biggers=do.call(rbind,biggers)
    diag(biggers)=NA
    biggers[biggers>1]=NA
    TDI=1-biggers
    TDH=sum(TDI,na.rm=TRUE)/n   #####COVARIATE
    #Clark Evans
    area=400
    mean_nn=mean(nndist(ppp_buffed))
    R_EXP=1/(2*sqrt(ppp_buffed$n/area))
    CEH = mean_nn / R_EXP   #####COVARIATE
  } else{
    SIMH=0   
    GINIH=0   
    TDH=0  
    CEH=0   
  }
  
  if (exists("crown_volumes")){
    if (length(crown_volumesT)>1){
      #Simpson index with 10m volume bins
      x=20
      H_Classes=round(crown_volumesT/x)*x
      frequencies=table(H_Classes)
      SIM_1=frequencies*(frequencies-1)
      SIMV=1 - (sum(SIM_1)/ (length(crown_volumesT)*(length(crown_volumesT)-1)))   #####COVARIATE  
      #GINI based on volume
      n=length(crown_volumesT)
      mu=mean(crown_volumesT)
      N=n * (n - 1)
      ox <- crown_volumesT[order(crown_volumesT)]
      dsum <- drop(crossprod(2 * 1:n - n - 1,  ox))
      GINIV = dsum / (mu * N)   #####COVARIATE
      #TD Diameter Differentation Index
      biggers=lapply(crown_volumesT,function(x){x/crown_volumesT})
      biggers=do.call(rbind,biggers)
      diag(biggers)=NA
      biggers[biggers>1]=NA
      TDI=1-biggers
      TDV=sum(TDI,na.rm=TRUE)/n   #####COVARIATE
    }
  } else{
    SIMV=0   
    GINIV=0  
    TDV=0   
  }
  ################################################################################################  
  #height covariates buffered
  if (length(LMheights)<1){
    height_sd=0  
    height_mean=0  
    height_med=0 
    height_max=0  
    height_min=0  
    #height_skew=0  
    #height_kurt=0  
    height_25=0  
    height_75=0  
    height_90=0  
  }else{
    heightB_sd=sd(LMheights)  #####COVARIATE
    heightB_mean=mean(LMheights)  #####COVARIATE
    heightB_med=median(LMheights)  #####COVARIATE
    heightB_max=max(LMheights)  #####COVARIATE
    heightB_min=min(LMheights)  #####COVARIATE
    #heightB_skew=skewness(LMheights)  #####COVARIATE
    #heightB_kurt=kurtosis(LMheights)  #####COVARIATE
    heightB_25=quantile(LMheights, .25)  #####COVARIATE
    heightB_75=quantile(LMheights, .75)  #####COVARIATE
    heightB_90=quantile(LMheights, .90)  #####COVARIATE
  }
  countB_P30=length(LMheights[LMheights>30])  #####COVARIATE
  countB_P20=length(LMheights[LMheights>20])  #####COVARIATE
  countB_P10=length(LMheights[LMheights>10])  #####COVARIATE
  countB=length(LMheights)  #####COVARIATE
  
  #height covariates not buffered 
  LMheights=points_notbuffed[((points_notbuffed[,1] %in% ppp_notbuffed$x) & (points_notbuffed[,2] %in% ppp_notbuffed$y)),3]
  if (length(LMheights)<1){
    height_sd=0  
    height_mean=0  
    height_max=0  
    height_25=0  
    height_75=0  
  }else{
    height_sd=sd(LMheights)  #####COVARIATE
    height_mean=mean(LMheights)  #####COVARIATE
    height_max=max(LMheights)  #####COVARIATE
    height_25=quantile(LMheights, .25)  #####COVARIATE
    height_75=quantile(LMheights, .75)  #####COVARIATE
  }
  count_P30=length(LMheights[LMheights>30])  #####COVARIATE
  count_P20=length(LMheights[LMheights>20])  #####COVARIATE
  count_P10=length(LMheights[LMheights>10])  #####COVARIATE
  count=length(LMheights)  #####COVARIATE
  
  #LM all heights in 20m cell using the local maxima
  LMheights=points[,3]
  if (length(LMheights)<1){
    heightLM_sd=0  
    heightLM_mean=0  
  }else{
    heightLM_sd=sd(LMheights)  #####COVARIATE
    heightLM_mean=mean(LMheights)  #####COVARIATE
  }
  countLM_P30=length(LMheights[LMheights>30])  #####COVARIATE
  countLM_P20=length(LMheights[LMheights>20])  #####COVARIATE
  countLM_P10=length(LMheights[LMheights>10])  #####COVARIATE
  countLM=length(LMheights)  #####COVARIATE
  
  out_in_ratioC=count/countLM
  out_in_ratioH=height_mean/heightLM_mean
  
  return (c(#general CHM metrics
    total_mean, total_median, total_sd, total_range, total_skew, total_kurt, Perc_Above_P20, Perc_Above_P40,
    Perc_Above_P60, Perc_Above_P80, Perc_Above_P95, Perc_Above_5m, Perc_Above_10m, Perc_Above_15m, Perc_Above_20m,
    Perc_Above_25m, Perc_Above_35m, Percentile_20, Percentile_40, Percentile_60, Percentile_80, Percentile_95,
    #Canopy complexity metrics
    mean_dists_NB, SD_dists_NB, skew_dists_NB, kurt_dists_NB,mean_dists_B, SD_dists_B, skew_dists_B, kurt_dists_B,
    mean_growing_space, sd_growing_space, max_growing_space, mean_growing_spaceNB, sd_growing_spaceNB,
    canopy_cover_NB, csr, rugoseS_count, rugoseS_mean, rugoseS_sd, rugoseL_mean, rugoseL_sd, 
    pt2Sedge_mean, pt2Sedge_max, mean_off_cent, sd_off_cent, SIMH, GINIH, TDH, CEH, SIMV, GINIV, TDV,out_in_ratioC,
    out_in_ratioH,
    #Tree metrics
    mean_crown_area_B, SD_crown_area_B, max_crown_area_B, watershed_count_inplot, mean_crown_area_NB,
    max_crown_area_NB, encroaching_tree_count, encroach_area, encroach_mean, perc_in_tally, 
    biomass_in, biomass_tot, crown_volume_tot, crown_volume_in, basal_area_tot, basal_area_in,
    softwoodyness, mean_crown_volumes, SD_crown_volumes, mean_top_ht, tot_top_ht, mean_ht, 
    mean_residuals,sd_residuals, mean_positivity, sd_positivity, mean_top_positivity, mean_raw_pointy,
    sd_raw_pointy, mean_pointy, sd_pointy, mean_off_cent, sd_off_cent, pt2WSedge_mean, pt2WSedge_sd, 
    pt2WSedge_max, rugose_WS_perc, heightB_mean, heightB_sd, heightB_med, heightB_max, heightB_min,
    heightB_25, heightB_75, heightB_90, countB_P30, countB_P20, countB_P10, countB, height_mean, 
    height_sd, height_max, height_25, height_75, count_P30, count_P20, count_P10, count, heightLM_mean,
    heightLM_sd, countLM_P30, countLM_P20, countLM_P10, countLM,
    #Competition metrics summed
    mean_obscured,sd_obscured,mean_shadowed,sd_shadowed,kurt_shadowed,sky_view_areas_m,sky_view_areas_sd,
    neighbor_top_dist_m,neighbor_top_dist_sd
  ))
}

##########
#High points finder algorithm
#A new method for generating canopy height models from discrete return LiDAR point clouds
#Liu and Dong 2014
pnt_fnder=function(inFile){
  #number of points is twice the expected point density given a 1.5m window
  hmm=nn2(inFile[,1:2],treetype='kd',k=round(nrow(inFile)*2*.017675),searchtype = 'radius', radius=1)
  zs_75_quantile=apply(hmm$nn.idx, 1, function(x){quantile(inFile[x[x>0],3],.75)})
  high_pts=inFile[,][inFile[,3]>=zs_75_quantile,]
  #high_pts[complete.cases(high_pts),]
}

tiner=function(high_pts){
  #make tin
  high_pts=high_pts[!duplicated(high_pts[,1:2]),]
  tin=tryCatch({
    tri.mesh(high_pts[,1], high_pts[,2] )
  }, error=function(e){
    high_pts <- high_pts[sample(nrow(high_pts)),]
    tri.mesh(high_pts[,1], high_pts[,2] )
  }
  )
  #make raster cell centers
  grd=data.frame(expand.grid(x=seq(from=e_b[1]+.375, to=e_b[2], by=0.75),y=seq(from=e_b[3]+.375, to=e_b[4],by=0.75)))
  
  #find triangle corners for each raster cell
  nrst_verticies=apply(grd,1, function(x) tri.find(tin, x[1], x[2]))
  grd=cbind(grd, data.frame(matrix(unlist(nrst_verticies), nrow=length(nrst_verticies), byrow=T)))
  #basic linear inerpolation within each triangle
  cell_values=apply(grd,1,function(x){
    p1=high_pts[as.numeric(x[3]),1:3]
    p2=high_pts[as.numeric(x[4]),1:3]
    p3=high_pts[as.numeric(x[5]),1:3]
    d_p1=sqrt((x[1]-p1[1])^2+(x[2]-p1[2])^2)
    d_p2=sqrt((x[1]-p2[1])^2+(x[2]-p2[2])^2)
    d_p3=sqrt((x[1]-p3[1])^2+(x[2]-p3[2])^2)
    
    total_d=sum(d_p1,d_p2,d_p3, na.rm=TRUE)
    total_p=sum((1-d_p1/total_d), (1-d_p2/total_d), (1-d_p3/total_d), na.rm=TRUE)
    
    #eek
    z_value=sum(p1[3]*((1-d_p1/total_d)/total_p), p2[3]*((1-d_p2/total_d)/total_p), p3[3]*((1-d_p3/total_d)/total_p),na.rm=TRUE)
    z_value
  })
  cell_values=matrix(cell_values, nrow=sqrt(nrow(grd)), ncol=sqrt(nrow(grd)))
  rotate <- function(x) t(apply(x, 2, rev))
  cell_values=rotate(rotate(rotate(cell_values)))
  r <- raster(cell_values, xmn=e_b[1], xmx=e_b[2], ymn=e_b[3], ymx=e_b[4])
  r=raster::mask(r, plotShape_b)
  
  #find rumple coefficent from shawn's paper
  tri_nodes=data.frame(triangles(tin))[,1:3]
  areas=apply(tri_nodes,1, function(x){
    tri_coords=high_pts[x,1:3]
    onetwo=sqrt((tri_coords[1,1]-tri_coords[2,1])^2+(tri_coords[1,2]-tri_coords[2,2])^2+(tri_coords[1,3]-tri_coords[2,3])^2)
    onethree=sqrt((tri_coords[1,1]-tri_coords[3,1])^2+(tri_coords[1,2]-tri_coords[3,2])^2+(tri_coords[1,3]-tri_coords[3,3])^2)
    twothree=sqrt((tri_coords[2,1]-tri_coords[3,1])^2+(tri_coords[2,2]-tri_coords[3,2])^2+(tri_coords[2,3]-tri_coords[3,3])^2)
    sqrt(((onetwo+onethree+twothree)/2)*(((onetwo+onethree+twothree)/2)-onetwo)*(((onetwo+onethree+twothree)/2)-onethree)*(((onetwo+onethree+twothree)/2)-twothree))
  })
  rumple_ratio=sum(areas)/400
  list(r, rumple_ratio)
}

#####################INITIALIZE######################################
#inShape_b=readShapePoly("G:\\7_1_2016_Reorganized_Training\\Plot_Shapefiles\\Baxter_Plots\\Baxter_Ninths_Buffed.shp")
#inShape=readShapePoly("G:\\7_1_2016_Reorganized_Training\\Plot_Shapefiles\\Baxter_Plots\\Baxter_Ninths.shp")

#plots=levels(inShape[1,]$Plot_Subpl)
#plot_fuckedup=levels(inShape_b[1,]$Plot_Subpl)
directory="F:\\Demeritt\\State_LiDAR\\LIDAR"
files=list.files(directory, pattern="*.las")
SiteProd= raster("E:\\Aaron's site productivity\\Tiff_Aaron_Prod.tif")

cl=makeCluster(40)
registerDoParallel(cl)
for (file in 1:length(files)){
  tryCatch({
    inpoints=tryCatch({
      readLAS(paste(directory,files[file], sep="\\"), short=FALSE)
    }, error=function(e){
      NULL
    }
    )
    e=extent(inpoints[,1:2])
    start_pts=expand.grid(x=seq(from=e[1], to=e[2], by=10),y=seq(from=e[3], to=e[4],by=10))
    extents=NULL
    area=100
    for (botpt in 1:nrow(start_pts)){
      plot_ext=extent(start_pts[botpt, 1], start_pts[botpt, 1]+10, start_pts[botpt,2], start_pts[botpt,2]+10)
      plot_ext=plot_ext+10
      extents=append(extents,plot_ext)
    }

    results=NULL
    results=foreach(p=1:length(extents), .combine='rbind', .errorhandling="pass") %dopar% {
      library("raster", lib.loc="~/R/win-library/3.3")
      library("spatstat", lib.loc="~/R/win-library/3.3")
      library("plyr", lib.loc="~/R/win-library/3.3")
      library("doParallel", lib.loc="~/R/win-library/3.3")
      library("aspace", lib.loc="~/R/win-library/3.3")
      library("maptools", lib.loc="~/R/win-library/3.3")
      library("moments", lib.loc="~/R/win-library/3.3")
      library("rLiDAR", lib.loc="~/R/win-library/3.2")
      library("circular", lib.loc="~/R/win-library/3.3")
      library("fields", lib.loc="~/R/win-library/3.3")
      library("RANN")
      library("tripack")
      e_b=extents[p][[1]]
      e=e_b-10
      plotShape_b=as(e_b, 'SpatialPolygons')
      plotShape=as(e, 'SpatialPolygons')
      
      #clip lidar plot
      inFile=inpoints[inpoints[,1]>e_b[1] & inpoints[,1]<e_b[2] & inpoints[,2]>e_b[3] & inpoints[,2]<e_b[4],]
      
      if ((!is.null(inFile)) &(length(plotShape)>0)){
        high_pts=pnt_fnder(inFile)
        ppmr<- raster(plotShape, ncol=(e[2]-e[1]), nrow=(e[4]-e[3]))
        ppm <- rasterize(subset(inFile,inFile[,5]==1)[, 1:2], ppmr, subset(inFile,inFile[,5]==1)[, 5], fun='count')
        ppm[is.na(ppm[])] <- 0  
        ppm <- raster::mask(ppm, plotShape)
        cell_ppm=cellStats(ppm, stat='mean')
        tin <- tiner(high_pts)
        DSM=tin[[1]]
        rumple=tin[[2]]
        grounds=inFile[inFile[,9]==2,]
        tin=tiner(grounds)
        DEM=tin[[1]]
        tiff=DSM-DEM
        crs(DEM)="+proj=utm +zone=19. +datum=NAD83"
        slope=mean(na.omit(getValues(terrain(DEM, opt='slope', unit='degrees'))))
        degrees=na.omit(getValues(terrain(DEM, opt='aspect', unit='degrees')))
        circle=circular(degrees,units="degrees")
        aspect=mean.circular(circle)
        if(aspect<0){
          aspect=360+aspect
        }
        TRI=mean(na.omit(getValues(terrain(DEM, opt='TRI'))))
        
        plot_coord=coordinates(plotShape)
        #WARNING, takes a really long time
        prod=extract(SiteProd,matrix(plot_coord, ncol=2))   
        
        if (length(na.omit(getValues(ppmi)))<95){
          NULL}
        else{
          #cell_ppm=nrow(subset(inFile,inFile[,5]==1))/400  #####COVARIATE
          scan_angle_mean=mean(abs(inFile[,10]))
          scan_angle_sd=sd(abs(inFile[,10]))
          leaves="off"  #####COVARIATE
          my_metrics=lclmaxer(tiff, plotShape,plotShape_b)
          my_metrics[is.na(my_metrics)] <- 0
          c(plot_coord, my_metrics,cell_ppm,scan_angle_mean,scan_angle_sd, leaves, prod, slope, aspect, TRI)
        }
      }
    }
    colnames(results)=c("INV","Aquisiton","Plot",
                        #general CHM metrics
                        "total_mean","total_median","total_sd","total_range","total_skew","total_kurt","Perc_Above_P20","Perc_Above_P40",
                        "Perc_Above_P60","Perc_Above_P80","Perc_Above_P95","Perc_Above_5m","Perc_Above_10m","Perc_Above_15m","Perc_Above_20m",
                        "Perc_Above_25m","Perc_Above_35m","Percentile_20","Percentile_40","Percentile_60","Percentile_80","Percentile_95",
                        #Canopy complexity metrics
                        "mean_dists_NB","SD_dists_NB","skew_dists_NB","kurt_dists_NB","mean_dists_B","SD_dists_B","skew_dists_B","kurt_dists_B",
                        "mean_growing_space","sd_growing_space","max_growing_space","mean_growing_spaceNB","sd_growing_spaceNB",
                        "canopy_cover_NB","csr","rugoseS_count","rugoseS_mean","rugoseS_sd","rugoseL_mean","rugoseL_sd", 
                        "pt2Sedge_mean", "pt2Sedge_max", "mean_off_cent","sd_off_cent","SIMH","GINIH","TDH","CEH","SIMV","GINIV","TDV","out_in_ratioC",
                        "out_in_ratioH",
                        #Tree metrics
                        "mean_crown_area_B","SD_crown_area_B","max_crown_area_B","watershed_count_inplot","mean_crown_area_NB",
                        "max_crown_area_NB","encroaching_tree_count","encroach_area","encroach_mean","perc_in_tally",
                        "biomass_in","biomass_tot","crown_volume_tot","crown_volume_in","basal_area_tot","basal_area_in",
                        "softwoodyness","mean_crown_volumes","SD_crown_volumes","mean_top_ht","tot_top_ht","mean_ht",
                        "mean_residuals","sd_residuals","mean_positivity","sd_positivity","mean_top_positivity","mean_raw_pointy",
                        "sd_raw_pointy","mean_pointy","sd_pointy","mean_off_cent","sd_off_cent","pt2WSedge_mean","pt2WSedge_sd",
                        "pt2WSedge_max","rugose_WS_perc","heightB_mean","heightB_sd","heightB_med","heightB_max","heightB_min",
                        "heightB_25","heightB_75","heightB_90","countB_P30","countB_P20","countB_P10","countB","height_mean",
                        "height_sd","height_max","height_25","height_75","count_P30","count_P20","count_P10","count", "heightLM_mean",
                        "heightLM_sd", "countLM_P30", "countLM_P20", "countLM_P10", "countLM",
                        #Competition metrics summed
                        "mean_obscured","sd_obscured","mean_shadowed","sd_shadowed","kurt_shadowed","sky_view_areas_m","sky_view_areas_sd",
                        "neighbor_top_dist_m","neighbor_top_dist_sd",
                        #Aquisition metrics
                        "cell_ppm","scan_angle_mean","scan_angle_sd", "leaves", "prod")
    
    write.table(results, file=paste(directory,'\\', substr(files[file],1,nchar(files[file])-4),'.csv', sep=""), sep=",", col.names= TRUE, row.names=FALSE)
    print(file/length(files))
  }, error=function(e){cat("ERROR : ", conditionMessage(e), " WITH ", files[file], "              ")})
}
stopCluster(cl) 
