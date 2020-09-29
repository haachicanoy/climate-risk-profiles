# Load climate files per county
# By: H. Achicanoy & A. Esquivel
# CIAT, 2020

options(warn = -1, scipen = 999)

suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, raster, vroom, sp))

# Please check the correct county name within: Country_Counts.xlsx
get_gcm_data <- function(country = 'Pakistan',
                         county  = 'Kashmore',
                         iso3    = 'PAK',
                         adm_lvl = 3,
                         gcm_list){
  
  country <<- country
  county  <<- county
  # Paths
  OSys <- Sys.info()[1]
  root <<- switch(OSys,
                  'Linux'   = '/home/jovyan/work/cglabs',
                  'Windows' = '//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles')
  
  # Scripts
  source(paste0(root,'/scripts/win_parallelization.R'))
  
  # Load county shapefile
  shp <<- raster::shapefile(paste0(root,'/data/shps/',country,'/',iso3,'_adm',adm_lvl,'.shp'))
  glue::glue('shp <- shp[shp@data$NAME_{adm_lvl} == county,]') %>%
    as.character %>%
    parse(text = .) %>% 
    eval(expr = ., envir = .GlobalEnv)
  
  # Load id coords
  crd <<- vroom(paste0(root,'/data/id_all_country.csv'), delim = ',')
  crd <<- crd %>%
    dplyr::filter(Country == country)
  pnt <<- crd %>% dplyr::select(x,y) %>% sp::SpatialPoints(coords = .)
  crs(pnt) <<- crs(shp)
  # Filter coordinates that are present in the county
  pnt <<- sp::over(pnt, shp) %>% data.frame %>% dplyr::select(ISO) %>% complete.cases() %>% which()
  crd <<- crd[pnt,]
  crd <<- crd %>%
    dplyr::mutate(county = county)
  
  # GCM climate
  # gcm_list <<- c('ipsl_cm5a_mr','miroc_esm_chem','ncc_noresm1_m')
  gcm_list <<- gcm_list
  periods  <<- c('1971_2000','2021_2045','2041_2065')
  
  for(period in periods){
    
    period <<- period
    
    cl <- createCluster(3, export = list("gcm_list","crd","country","root","period","county"), lib = list("tidyverse","raster"))
    
    gcm_list %>% parallel::parLapply(cl, ., function(gcm){
      gcm <<- gcm
      if(period == '1971_2000'){
        prec_fls <- list.files(paste0(root,'/data/gcm_0_05deg_lat/',tolower(country),'/',gcm,'/',period,'/by-month'), pattern = '^prec', full.names = T)
      } else {
        prec_fls <- list.files(paste0(root,'/data/gcm_0_05deg_lat/',tolower(country),'/',gcm,'/',period,'/rcp85/by-month'), pattern = '^prec', full.names = T)
      }
      prec <- prec_fls %>% purrr::map(., .f=function(x){
        date <- basename(x)
        date <- date %>% gsub('prec_','',.) %>% gsub('.nc','',.) %>% gsub('_','-',.)
        rs   <- raster::stack(x)
        df   <- raster::extract(rs,crd %>% dplyr::select(x,y)) %>% data.frame
        colnames(df) <- paste0(date,'-',c(paste0('0',1:9),10:ncol(df)))
        return(df)
      })
      prec <- do.call(cbind, prec)
      prec$id <- crd$id
      
      if(period == '1971_2000'){
        tmax_fls <- list.files(paste0(root,'/data/gcm_0_05deg_lat/',tolower(country),'/',gcm,'/',period,'/by-month'), pattern = '^tmax', full.names = T)
      } else {
        tmax_fls <- list.files(paste0(root,'/data/gcm_0_05deg_lat/',tolower(country),'/',gcm,'/',period,'/rcp85/by-month'), pattern = '^tmax', full.names = T)
      }
      tmax <- tmax_fls %>% purrr::map(., .f=function(x){
        date <- basename(x)
        date <- date %>% gsub('tmax_','',.) %>% gsub('.nc','',.) %>% gsub('_','-',.)
        rs   <- raster::stack(x)
        df   <- raster::extract(rs,crd %>% dplyr::select(x,y)) %>% data.frame
        colnames(df) <- paste0(date,'-',c(paste0('0',1:9),10:ncol(df)))
        return(df)
      })
      tmax <- do.call(cbind, tmax)
      tmax$id <- crd$id
      
      if(period == '1971_2000'){
        tmin_fls <- list.files(paste0(root,'/data/gcm_0_05deg_lat/',tolower(country),'/',gcm,'/',period,'/by-month'), pattern = '^tmin', full.names = T)
      } else {
        tmin_fls <- list.files(paste0(root,'/data/gcm_0_05deg_lat/',tolower(country),'/',gcm,'/',period,'/rcp85/by-month'), pattern = '^tmin', full.names = T)
      }
      tmin <- tmin_fls %>% purrr::map(., .f=function(x){
        date <- basename(x)
        date <- date %>% gsub('tmin_','',.) %>% gsub('.nc','',.) %>% gsub('_','-',.)
        rs   <- raster::stack(x)
        df   <- raster::extract(rs,crd %>% dplyr::select(x,y)) %>% data.frame
        colnames(df) <- paste0(date,'-',c(paste0('0',1:9),10:ncol(df)))
        return(df)
      })
      tmin <- do.call(cbind, tmin)
      tmin$id <- crd$id
      
      if(period == '1971_2000'){
        srad_fls <- list.files(paste0(root,'/data/gcm_0_05deg_lat/',tolower(country),'/',gcm,'/',period,'/by-month'), pattern = '^rsds', full.names = T)
      } else {
        srad_fls <- list.files(paste0(root,'/data/gcm_0_05deg_lat/',tolower(country),'/',gcm,'/',period,'/rcp85/by-month'), pattern = '^rsds', full.names = T)
      }
      srad <- srad_fls %>% purrr::map(., .f=function(x){
        date <- basename(x)
        date <- date %>% gsub('rsds_','',.) %>% gsub('.nc','',.) %>% gsub('_','-',.)
        rs   <- raster::stack(x)
        df   <- raster::extract(rs,crd %>% dplyr::select(x,y)) %>% data.frame * 0.0864
        colnames(df) <- paste0(date,'-',c(paste0('0',1:9),10:ncol(df)))
        return(df)
      })
      srad <- do.call(cbind, srad)
      srad$id <- crd$id
      
      prec2 <- prec %>% tidyr::pivot_longer(cols = colnames(prec)[-which(colnames(prec)=='id')],
                                            names_to = 'Date',
                                            values_to = 'prec')
      tmax2 <- tmax %>% tidyr::pivot_longer(cols = colnames(tmax)[-which(colnames(tmax)=='id')],
                                            names_to = 'Date',
                                            values_to = 'tmax')
      tmin2 <- tmin %>% tidyr::pivot_longer(cols = colnames(tmin)[-which(colnames(tmin)=='id')],
                                            names_to = 'Date',
                                            values_to = 'tmin')
      srad2 <- srad %>% tidyr::pivot_longer(cols = colnames(srad)[-which(colnames(srad)=='id')],
                                            names_to = 'Date',
                                            values_to = 'srad')
      clim <- dplyr::left_join(x = prec2, y = tmax2, by = c('id','Date'))
      clim <- dplyr::left_join(x = clim, y = tmin2, by = c('id','Date'))
      clim <- dplyr::left_join(x = clim, y = srad2, by = c('id','Date'))
      rm(prec, prec2, tmax, tmax2, tmin, tmin2, srad, srad2)
      
      his_gcm <- tibble::tibble(id      = crd$id,
                                x       = crd$x,
                                y       = crd$y,
                                ISO3    = crd$ISO3,
                                Country = crd$Country,
                                Climate = clim %>% dplyr::group_split(id))
      his_gcm <- his_gcm %>% tidyr::unnest(.)
      out <- paste0(root,'/data/gcm_0_05deg_lat_county/',tolower(country),'/',gcm,'/',period)
      if(!dir.exists(out)){dir.create(out, recursive = T)}
      out <- paste0(out,'/',tolower(county),'.fst')
      if(!file.exists(out)){
        fst::write_fst(his_gcm, out)
      } else {
        cat('File exists\n')
      }
      
    })
    parallel::stopCluster(cl)
    
  }
  return(cat('Process successfully finished\n'))
  
}
# Run once
# gcmList <- c("ipsl_cm5a_mr","miroc_esm_chem","ncc_noresm1_m","bnu_esm","cccma_canesm2","cmcc_cms","gfdl_esm2g")
gcmList <- c("ipsl_cm5a_mr","miroc_esm_chem","ncc_noresm1_m")
get_gcm_data(country='Kenya',county='Turkana',iso3='KEN',adm_lvl=1,gcm_list=gcmList)
