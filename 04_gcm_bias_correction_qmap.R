### Perform quantile-mapping bias correction of daily climate data
### H. Achicanoy - C. Navarro
### CIAT, 2020

# R options
options(warn = -1, scipen = 999)

# Load libraries
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(qmap, ncdf4, raster, tidyverse, compiler, vroom, gtools, fst))

# Paths
OSys <- Sys.info()[1]
root <<- switch(OSys,
                'Linux'   = '/home/jovyan/work/cglabs',
                'Windows' = '//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles')

# Scripts
source(paste0(root,'/scripts/win_parallelization.R'))

# Quantile-mapping bias correction function for available pixels with solar radiation from NASA within a country
BC_Qmap <- function(country   = "Ethiopia",
                    county    = "Arsi",
                    rcp       = "rcp85",
                    gcm       = "ipsl_cm5a_mr",
                    period    = "2021_2045",
                    ncores    = 10)
{
  
  bc_qmap <<- function(df_obs, df_his_gcm, df_fut_gcm){
    if('srad' %in% colnames(df_obs)){
      cat('> Fitting the Qmap function per variable\n')
      prec_fit <- qmap::fitQmap(obs=df_obs$prec, mod=df_his_gcm$prec, method="RQUANT", qstep=0.01, wet.day=TRUE, na.rm=TRUE)
      tmax_fit <- qmap::fitQmap(obs=df_obs$tmax, mod=df_his_gcm$tmax, method="RQUANT", qstep=0.01, wet.day=FALSE, na.rm=TRUE)
      tmin_fit <- qmap::fitQmap(obs=df_obs$tmin, mod=df_his_gcm$tmin, method="RQUANT", qstep=0.01, wet.day=FALSE, na.rm=TRUE)
      srad_fit <- qmap::fitQmap(obs=df_obs$srad, mod=df_his_gcm$srad, method="RQUANT", qstep=0.01, wet.day=FALSE, na.rm=TRUE)
      
      cat('> Doing bias correction per variable historical GCMs\n')
      bc_his_gcm <- df_his_gcm
      bc_his_gcm$prec <- qmap::doQmap(x=df_his_gcm$prec, prec_fit, type="linear")
      bc_his_gcm$tmax <- qmap::doQmap(x=df_his_gcm$tmax, tmax_fit, type="linear")
      bc_his_gcm$tmin <- qmap::doQmap(x=df_his_gcm$tmin, tmin_fit, type="linear")
      bc_his_gcm$srad <- qmap::doQmap(x=df_his_gcm$srad, srad_fit, type="linear")
      
      cat('> Doing bias correction per variable future GCMs\n')
      bc_fut_gcm <- df_fut_gcm
      bc_fut_gcm$prec <- qmap::doQmap(x=df_fut_gcm$prec, prec_fit, type="linear")
      bc_fut_gcm$tmax <- qmap::doQmap(x=df_fut_gcm$tmax, tmax_fit, type="linear")
      bc_fut_gcm$tmin <- qmap::doQmap(x=df_fut_gcm$tmin, tmin_fit, type="linear")
      bc_fut_gcm$srad <- qmap::doQmap(x=df_fut_gcm$srad, srad_fit, type="linear")
      
    } else {
      cat('> Fitting the Qmap function per variable\n')
      prec_fit <- qmap::fitQmap(obs=df_obs$prec, mod=df_his_gcm$prec, method="RQUANT", qstep=0.01, wet.day=TRUE, na.rm=TRUE)
      tmax_fit <- qmap::fitQmap(obs=df_obs$tmax, mod=df_his_gcm$tmax, method="RQUANT", qstep=0.01, wet.day=FALSE, na.rm=TRUE)
      tmin_fit <- qmap::fitQmap(obs=df_obs$tmin, mod=df_his_gcm$tmin, method="RQUANT", qstep=0.01, wet.day=FALSE, na.rm=TRUE)
      
      cat('> Doing bias correction per variable historical GCMs\n')
      bc_his_gcm <- df_his_gcm
      bc_his_gcm$prec <- qmap::doQmap(x=df_his_gcm$prec, prec_fit, type="linear")
      bc_his_gcm$tmax <- qmap::doQmap(x=df_his_gcm$tmax, tmax_fit, type="linear")
      bc_his_gcm$tmin <- qmap::doQmap(x=df_his_gcm$tmin, tmin_fit, type="linear")
      
      cat('> Doing bias correction per variable future GCMs\n')
      bc_fut_gcm <- df_fut_gcm
      bc_fut_gcm$prec <- qmap::doQmap(x=df_fut_gcm$prec, prec_fit, type="linear")
      bc_fut_gcm$tmax <- qmap::doQmap(x=df_fut_gcm$tmax, tmax_fit, type="linear")
      bc_fut_gcm$tmin <- qmap::doQmap(x=df_fut_gcm$tmin, tmin_fit, type="linear")
    }
    bc_data <- list(His = bc_his_gcm,
                    Fut = bc_fut_gcm)
    return(list(bc_data))
  }
  
  cat(paste0(' *** Performing Quantile-mapping bias correction for available pixels with SRAD (NASA) within ',county,' in the period ',period,', using: ',rcp,', GCM: ',gcm,'***\n'))
  cat(paste0('>>> Loading obs data\n'))
  
  obsDir <<- paste0(root,"/data/observational_data/",tolower(country))
  if(!file.exists(paste0(obsDir,'/',tolower(county),'.fst'))){
    his_obs <<- fst::read_fst(paste0(obsDir,'/',tolower(county),'_prec_temp.fst'))
    his_obs <<- his_obs %>%
      tidyr::nest(Climate = c('id','Date','prec','tmax','tmin')) %>%
      dplyr::rename(id = 'id1') %>%
      dplyr::select(id, everything(.))
  } else {
    his_obs <<- fst::read_fst(paste0(obsDir,'/',tolower(county),'.fst'))
    his_obs <<- his_obs %>%
      tidyr::nest(Climate = c('id','Date','prec','tmax','tmin','srad')) %>%
      dplyr::rename(id = 'id1') %>%
      dplyr::select(id, everything(.))
  }
  impute_missings <- function(tbl = clim_data){
    Climate <- 1:nrow(tbl) %>%
      purrr::map(.f = function(i){
        df <- tbl$Climate[[i]]
        if(sum(is.na(df$tmax)) > 0){
          df$tmax[which(is.na(df$tmax))] <- median(df$tmax, na.rm = T)
        }
        if(sum(is.na(df$tmin)) > 0){
          df$tmin[which(is.na(df$tmin))] <- median(df$tmin, na.rm = T)
        }
        if(sum(is.na(df$srad)) > 0){
          df$srad[which(is.na(df$srad))] <- median(df$srad, na.rm = T)
        }
        if(sum(is.na(df$prec)) > 0){
          df$prec[which(is.na(df$prec))] <- median(df$prec, na.rm = T)
        }
        return(df)
      })
    
    return(Climate)
  }
  his_obs$Climate <- impute_missings(tbl = his_obs)
  
  cat(paste0('>>> Loading historical GCM data\n'))
  hisGCMDir <<- paste0(root,"/data/gcm_0_05deg_lat_county/",tolower(country),"/",gcm,"/1971_2000")
  his_gcm <<- fst::read_fst(paste0(hisGCMDir,'/',tolower(county),'.fst'))
  his_gcm <<- his_gcm %>%
    tidyr::nest(Climate = c('id','Date','prec','tmax','tmin','srad')) %>%
    dplyr::rename(id = 'id1') %>%
    dplyr::select(id, everything(.))
  
  cat(paste0('>>> Loading future GCM data\n'))
  futGCMDir <<- paste0(root,"/data/gcm_0_05deg_lat_county/",tolower(country),"/",gcm,"/",period)
  fut_gcm <<- fst::read_fst(paste0(futGCMDir,'/',tolower(county),'.fst'))
  fut_gcm <<- fut_gcm %>%
    tidyr::nest(Climate = c('id','Date','prec','tmax','tmin','srad')) %>%
    dplyr::rename(id = 'id1') %>%
    dplyr::select(id, everything(.))
  
  his_gcm_bc <<- his_gcm
  fut_gcm_bc <<- fut_gcm
  
  cl <- createCluster(ncores, export = list("root","obsDir","his_obs","hisGCMDir","his_gcm","futGCMDir","fut_gcm","bc_qmap","his_gcm_bc","fut_gcm_bc"), lib = list("tidyverse","raster","qmap"))
  
  bc_data <- 1:nrow(his_obs) %>% parallel::parLapply(cl, ., function(i){
    bc_data <<- bc_qmap(df_obs    = his_obs$Climate[[i]],
                       df_his_gcm = his_gcm$Climate[[i]],
                       df_fut_gcm = fut_gcm$Climate[[i]])
    return(bc_data)
  })
  parallel::stopCluster(cl)
  his_gcm_bc$Climate <- bc_data %>% purrr::map(1) %>% purrr::map(1)
  his_gcm_bc <- his_gcm_bc %>% tidyr::unnest(.)
  fut_gcm_bc$Climate <- bc_data %>% purrr::map(1) %>% purrr::map(2)
  fut_gcm_bc <- fut_gcm_bc %>% tidyr::unnest(.)
  
  pDir <- paste0(root,'/data/bc_quantile_0_05deg_lat_county/',tolower(country),'/',gcm,'/1971_2000')
  if(!dir.exists(pDir)){dir.create(pDir, recursive = T)}
  fDir <- paste0(root,'/data/bc_quantile_0_05deg_lat_county/',tolower(country),'/',gcm,'/',period)
  if(!dir.exists(fDir)){dir.create(fDir, recursive = T)}
  
  if('srad' %in% colnames(his_obs$Climate[[1]])){
    outHis <- paste0(pDir,'/',tolower(county),'.fst')
    if(!file.exists(outHis)){
      fst::write_fst(his_gcm_bc,outHis)
    }
    outFut <- paste0(fDir,'/',tolower(county),'.fst')
    if(!file.exists(outFut)){
      fst::write_fst(fut_gcm_bc,outFut)
    }
  } else {
    outHis <- paste0(pDir,'/',tolower(county),'_prec_temp.fst')
    if(!file.exists(outHis)){
      fst::write_fst(his_gcm_bc,outHis)
    }
    outFut <- paste0(fDir,'/',tolower(county),'_prec_temp.fst')
    if(!file.exists(outFut)){
      fst::write_fst(fut_gcm_bc,outFut)
    }
  }
  
}
# Run once
periodList <- c('2021_2045','2041_2065')
rcpList    <- 'rcp85'
# gcmList <- c("ipsl_cm5a_mr","miroc_esm_chem","ncc_noresm1_m","bnu_esm","cccma_canesm2","cmcc_cms","gfdl_esm2g")
gcmList <- c("ipsl_cm5a_mr","miroc_esm_chem","ncc_noresm1_m")
for(p in periodList){
  for(gcm in gcmList){
    BC_Qmap(country='Niger',county='Tillaberi',rcp='rcp85',gcm=gcm,period=p,ncores=5)
  }
}
