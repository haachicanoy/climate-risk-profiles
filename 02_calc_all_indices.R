# Calculate all agro-climatic indices
# A. Esquivel and H. Achicanoy
# CIAT, 2020

options(warn = -1, scipen = 999)

suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, raster, ncdf4, sf, future, furrr, lubridate, glue, vroom, sp, fst, compiler))

OSys <- Sys.info()[1]
root <<- switch(OSys,
                'Linux'   = '/home/jovyan/work/cglabs',
                'Windows' = '//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles')

calc_indices <- function(country = 'Mozambique',
                         county  = 'Sofala',
                         iso3c   = 'MOZ',
                         adm_lvl = 1,
                         seasons = NULL, # Seasons manually defined
                         n_ssns  = 2,    # 2-seasons automatically defined
                         n_wtts  = 100,  # 100-wettest days
                         gcm     = NULL,
                         period  = '1985_2015',
                         time    = 'past',
                         big_cnt = TRUE,
                         ncores  = 10){
  
  country <<- country
  county  <<- county
  
  # Load county shapefile
  shp <<- raster::shapefile(paste0(root,'/data/shps/',country,'/',iso3c,'_adm',adm_lvl,'.shp'))
  glue::glue('shp <<- shp[shp@data$NAME_{adm_lvl} == county,]') %>%
    as.character %>%
    parse(text = .) %>%
    eval(expr = ., envir = .GlobalEnv)
  
  # Load id coords
  crd <- vroom(paste0(root,'/data/id_all_country.csv'), delim = ',')
  crd <- crd %>%
    dplyr::filter(Country == country)
  pnt <- crd %>% dplyr::select(x,y) %>% sp::SpatialPoints(coords = .)
  crs(pnt) <- crs(shp)
  # Filter coordinates that are present in the county
  pnt <- sp::over(pnt, shp) %>% data.frame %>% dplyr::select(ISO) %>% complete.cases() %>% which()
  crd <- crd[pnt,]
  crd <<- crd
  
  timList <- c('past','future')
  periodList <- c('2021_2045','2041_2065')
  
  source(paste0(root,'/scripts/indices.R'))
  
  # Paths
  obsDir <- paste0(root,'/data/observational_data/',tolower(country))
  futDir <- paste0(root,'/data/bc_quantile_0_05deg_lat_county/',tolower(country),'/',gcm,'/',period)
  outDir <- ifelse(test = time == 'past',
                   yes  = paste0(root,'/results/',country,'/',time),
                   no   = paste0(root,'/results/',country,'/',time,'/',gcm,'/',period))
  
  if(big_cnt){
    out <- paste0(outDir,'/',county,'_',period,'_corrected_idw.fst')
  } else {
    out <- paste0(outDir,'/',county,'_',period,'_corrected.fst')
  }
  
  if(!file.exists(out)){
    # Soil data
    Soil <- fst::read.fst(paste0(root,'/data/soilcp_data.fst')) %>%
      tibble::as_tibble() %>%
      dplyr::select(id, soilcp) %>%
      dplyr::filter(id %in% dplyr::pull(crd, id))
    
    # Read climate data
    if(time == 'past'){
      if(file.exists(paste0(obsDir,'/',tolower(county),'.fst'))){
        clim_data <- fst::read_fst(paste0(obsDir,'/',tolower(county),'.fst'))
        clim_data <- clim_data %>%
          tidyr::nest(Climate = c('id','Date','prec','tmax','tmin','srad')) %>%
          dplyr::rename(id = 'id1') %>%
          dplyr::select(id, everything(.))
      } else {
        if(file.exists(paste0(obsDir,'/',tolower(county),'_prec_temp.fst'))){
          clim_data <- fst::read_fst(paste0(obsDir,'/',tolower(county),'_prec_temp.fst'))
          clim_data <- clim_data %>%
            tidyr::nest(Climate = c('id','Date','prec','tmax','tmin')) %>%
            dplyr::rename(id = 'id1') %>%
            dplyr::select(id, everything(.))
          } else {
          clim_data <- readRDS(paste0(obsDir,'/',tolower(county),'.RDS'))
        }
      }
    } else {
      if(file.exists(paste0(futDir,'/',tolower(county),'.fst'))){
        clim_data <- fst::read_fst(paste0(futDir,'/',tolower(county),'.fst'))
        clim_data <- clim_data %>%
          tidyr::nest(Climate = c('id','Date','prec','tmax','tmin','srad')) %>%
          dplyr::rename(id = 'id1') %>%
          dplyr::select(id, everything(.))
      } else {
        if(file.exists(paste0(futDir,'/',tolower(county),'_prec_temp.fst'))){
          clim_data <- fst::read_fst(paste0(futDir,'/',tolower(county),'_prec_temp.fst'))
          clim_data <- clim_data %>%
            tidyr::nest(Climate = c('id','Date','prec','tmax','tmin')) %>%
            dplyr::rename(id = 'id1') %>%
            dplyr::select(id, everything(.))
        } else {
          clim_data <- readRDS(paste0(futDir,'/',tolower(county),'.RDS'))
        }
      }
    }
    
    # Impute missing data
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
    clim_data$Climate <- impute_missings(tbl = clim_data)
    
    if(big_cnt){
      set.seed(1235)
      sample_n  <- nrow(clim_data)*0.3
      id_sample <- sample(unique(clim_data$id), sample_n) 
      clim_data <- dplyr::filter(clim_data, id %in% id_sample)
    }
    
    run_pixel <- function(id = 362540){
      
      cat(' --- Obtain complete time series per pixel\n')
      tbl <- clim_data$Climate[[which(clim_data$id == id)]]
      tbl <- tbl %>%
        dplyr::mutate(year  = lubridate::year(as.Date(Date)),
                      month = lubridate::month(as.Date(Date)))
      years <- tbl$year %>% unique
      
      cat(' --- Calculate water balance for complete time series\n')
      soilcp <- Soil$soilcp[Soil$id == id]
      watbal_loc <- watbal_wrapper(out_all = tbl, soilcp = soilcp)
      # watbal_loc$IRR <- watbal_loc$Etmax - watbal_loc$prec
      
      tbl <- tbl %>%
        dplyr::mutate(ERATIO = watbal_loc$ERATIO,
                      TAV    = (watbal_loc$tmin + watbal_loc$tmax)/2,
                      GDAY   = ifelse(TAV >= 6 & ERATIO >= 0.35, yes=1, no=0))
      
      cat(' --- Estimate growing seasons from water balance\n')
      
      ### CONDITIONS TO HAVE IN ACCOUNT
      # Length of growing season per year
      # Start: 5-consecutive growing days.
      # End: 12-consecutive non-growing days.
      
      # Run process by year
      lgp_year_pixel <- lapply(1:length(years), function(k){
        
        # Subsetting by year
        watbal_year <- tbl[tbl$year==years[k],]
        
        # Calculate sequences of growing and non-growing days within year
        runsDF <- rle(watbal_year$GDAY)
        runsDF <- data.frame(Lengths=runsDF$lengths, Condition=runsDF$values)
        runsDF$Condition <- runsDF$Condition %>% tidyr::replace_na(replace = 0)
        
        # Identify start and extension of each growing season during year
        if(!sum(runsDF$Lengths[runsDF$Condition==1] < 5) == length(runsDF$Lengths[runsDF$Condition==1])){
          
          LGP <- 0; LGP_seq <- 0
          for(i in 1:nrow(runsDF)){
            if(runsDF$Lengths[i] >= 5 & runsDF$Condition[i] == 1){
              LGP <- LGP + 1
              LGP_seq <- c(LGP_seq, LGP)
              LGP <- 0
            } else {
              if(LGP_seq[length(LGP_seq)]==1){
                if(runsDF$Lengths[i] >= 12 & runsDF$Condition[i] == 0){
                  LGP <- 0
                  LGP_seq <- c(LGP_seq, LGP)
                } else {
                  LGP <- LGP + 1
                  LGP_seq <- c(LGP_seq, LGP)
                  LGP <- 0
                }
              } else {
                LGP <- 0
                LGP_seq <- c(LGP_seq, LGP)
              }
            }
          }
          LGP_seq <- c(LGP_seq, LGP)
          LGP_seq <- LGP_seq[-c(1, length(LGP_seq))]
          runsDF$gSeason <- LGP_seq; rm(i, LGP, LGP_seq)
          LGP_seq <- as.list(split(which(runsDF$gSeason==1), cumsum(c(TRUE, diff(which(runsDF$gSeason==1))!=1))))
          
          # Calculate start date and extension of each growing season by year and pixel
          growingSeason <- lapply(1:length(LGP_seq), function(g){
            
            LGP_ini <- sum(runsDF$Lengths[1:(min(LGP_seq[[g]])-1)]) + 1
            LGP <- sum(runsDF$Lengths[LGP_seq[[g]]])
            results <- data.frame(id=tbl$id %>% unique, year=years[k], gSeason=g, SLGP=LGP_ini, LGP=LGP)
            return(results)
            
          })
          growingSeason <- do.call(rbind, growingSeason)
          if(nrow(growingSeason)>2){
            growingSeason <- growingSeason[rank(-growingSeason$LGP) %in% 1:2,]
            growingSeason$gSeason <- rank(growingSeason$SLGP)
            growingSeason <- growingSeason[order(growingSeason$gSeason),]
          }
          
        } else {
          
          growingSeason <- data.frame(id=tbl$id %>% unique, year=years[k], gSeason = 1:2, SLGP = NA, LGP = NA)
          
        }
        
        print(k)
        return(growingSeason)
        
      })
      lgp_year_pixel <- do.call(rbind, lgp_year_pixel); rownames(lgp_year_pixel) <- 1:nrow(lgp_year_pixel)
      if(length(seasons) == 1){
        lgp_year_pixel <- lgp_year_pixel %>%
          dplyr::filter(gSeason == 1)
      } else {
        if(length(seasons) == 2){
          lgp_year_pixel <- lgp_year_pixel %>%
            dplyr::filter(gSeason %in% 1:2)
        }
      }
      
      cat(' --- Calculate agro-climatic indices for an specific season\n')
      if(!is.null(seasons)){
        indices <- 1:length(seasons) %>%
          purrr::map(.f = function(i){
            season = seasons[[i]]
            # Season across two years
            if(sum(diff(season) < 0) > 0){
              pairs     <- NA; for(j in 1:length(years)-1){pairs[j] <- paste0(years[j:(j+1)], collapse = '-')}
              tbl_list  <- lapply(1:(length(years)-1), function(k){
                df <- tbl %>%
                  dplyr::filter(year %in% years[k:(k+1)])
                df$pairs <- paste0(years[k:(k+1)], collapse = '-')
                df1 <- df %>%
                  dplyr::filter(year == years[k] & month %in% season[1]:12)
                df2 <- df %>%
                  dplyr::filter(year == years[k+1] & month %in% 1:season[length(season)])
                df <- rbind(df1, df2); rm(df1, df2)
                return(df)
              })
              tbl_list <- dplyr::bind_rows(tbl_list)
            } else {
              # Season in one year
              tbl_list  <- lapply(1:length(years), function(k){
                df <- tbl %>%
                  dplyr::filter(year %in% years[k])
                df$pairs <- years[k]
                df <- df %>%
                  dplyr::filter(year == years[k] & month %in% season)
                return(df)
              })
              tbl_list <- dplyr::bind_rows(tbl_list)
            }
            idx <- tbl_list %>%
              dplyr::group_split(pairs) %>%
              purrr::map(.f = function(df){
                idx <- tibble::tibble(CDD  = calc_cddCMP(PREC = df$prec),
                                      P5D  = calc_p5dCMP(PREC = df$prec),
                                      P95  = calc_p95CMP(PREC = df$prec),
                                      NT35 = calc_htsCMP(tmax = df$tmax, t_thresh = 35),
                                      ndws = calc_wsdays(df$ERATIO, season_ini=1, season_end=length(df$ERATIO), e_thresh=0.5))
                return(idx)
              })
            idx <- dplyr::bind_rows(idx)
            idx$year <- tbl_list$pairs %>% unique
            if(sum(diff(season) < 0) > 0){
              idx$year <- substr(x = idx$year, start = 6, stop = 10) %>% as.numeric
            } else {
              idx$year <- idx$year %>% as.numeric
            }
            idx$season <- names(seasons)[i]
            idx$id <- id
            
            return(idx)
            
          })
        indices    <- dplyr::bind_rows(indices)
      } else {
        if(!is.null(n_ssns)){
          # 1. Split by semester
          if(n_ssns == 2){
            tbl_list <- lapply(1:length(years), function(k){
              df <- tbl %>%
                dplyr::filter(year %in% years[k])
              df$pairs <- years[k]
              df1 <- df %>%
                dplyr::filter(year == years[k] & month %in% 1:6)
              df1$season <- 's1'
              df2 <- df %>%
                dplyr::filter(year == years[k] & month %in% 7:12)
              df2$season <- 's2'
              df <- list(s1 = df1, s2 = df2)
              return(df)
            }) %>% purrr::flatten()
          } else { # Split by year
            if(n_ssns == 1){
              tbl_list <- lapply(1:length(years), function(k){
                df <- tbl %>%
                  dplyr::filter(year %in% years[k])
                df$pairs <- years[k]
                return(df)
              })
            }
          }
          # 2. Within each semester identify the n-wettest days per pixel
          wettest_days <- tbl_list %>%
            purrr::map(., .f = function(tbl){
              SummDays <- rsum.lapply(x = tbl$prec, n = n_wtts)
              WetDays  <- SummDays[which.max(cumulative.r.sum(SummDays))]
              WetDays  <- WetDays %>% purrr::map(2) %>% unlist %>% data.frame()
              return(WetDays)
            })
          # 3. Calculate the indices for these previously identified days
          indices <- 1:length(tbl_list) %>%
            purrr::map(.f = function(i){
              df  <- tbl_list[[i]]
              wt  <- wettest_days[[i]] %>% dplyr::pull(.)
              idx <- tibble::tibble(year   = df$year %>% unique(),
                                    season = df$season %>% unique(),
                                    CDD    = calc_cddCMP(PREC = df$prec[wt]),
                                    P5D    = calc_p5dCMP(PREC = df$prec[wt]),
                                    P95    = calc_p95CMP(PREC = df$prec[wt]),
                                    NT35   = calc_htsCMP(tmax = df$tmax[wt], t_thresh = 35),
                                    ndws   = calc_wsdays(df$ERATIO[wt], season_ini=1, season_end=length(df$ERATIO), e_thresh=0.5))
              return(idx)
            })
          indices    <- dplyr::bind_rows(indices)
          indices$id <- id
        }
      }
      all <- dplyr::full_join(x = indices, y = lgp_year_pixel, by = c('id','year')) %>% unique()
      
      return(all)
      
    }
    
    plan(cluster, workers = ncores)
    index_by_pixel <- clim_data %>%
      dplyr::pull(id) %>% 
      furrr::future_map(.x = ., .f = run_pixel) %>% 
      dplyr::bind_rows()
    gc()
    gc(reset = T)
    
    index_by_pixel$Country <- country
    index_by_pixel$county <- county
    index_by_pixel$ISO3 <- iso3c
    index_by_pixel <- dplyr::left_join(x = index_by_pixel, y = clim_data %>% dplyr::select(id, x, y) %>% unique, by = 'id')
    index_by_pixel <- index_by_pixel %>% dplyr::select(id, x, y, ISO3, Country, county, year, season, dplyr::everything(.))
    
    if(!dir.exists(outDir)){dir.create(outDir, recursive = T)}
    fst::write_fst(x = index_by_pixel, path = out)
    cat('>>> File created successfully ...\n')
  } else {
    cat('>>> File exists it is not necessary to create it again\n')
  }
  
}

# countyList <- c('Kaduna',
#                 'Plateau')
# for(i in 1:length(countyList)){
#   calc_indices(country = 'Nigeria',
#                county  = countyList[i],
#                iso3c   = 'NGA',
#                adm_lvl = 1,
#                seasons = list(s1 = 4:10),
#                gcm     = NULL,
#                period  = '1985_2015',
#                time    = 'past',
#                big_cnt = TRUE,
#                ncores  = 20)
# }
# 
# gcmList <- c("ipsl_cm5a_mr","miroc_esm_chem","ncc_noresm1_m")
# periodList <- c('2021_2045','2041_2065')
# for(i in 1:length(countyList)){
#   for(gcm in gcmList){
#     for(period in periodList){
#       calc_indices(country = 'Nigeria',
#                    county  = countyList[i],
#                    iso3c   = 'NGA',
#                    adm_lvl = 1,
#                    seasons = list(s1 = 4:10),
#                    gcm     = gcm,
#                    period  = period,
#                    time    = 'future',
#                    big_cnt = TRUE,
#                    ncores  = 20)
#     }
#   }
# }

countyList <- c('Siaya','Bungoma','Kakamega','Nyandarua')
for(i in 1:length(countyList)){
  calc_indices(country = 'Kenya',
               county  = countyList[i],
               iso3c   = 'KEN',
               adm_lvl = 2,
               seasons = list(s1 = 2:6, s2 = 7:12), # Seasons manually defined
               n_ssns  = NULL,    # 2-seasons automatically defined
               n_wtts  = NULL,  # 100-wettest days
               gcm     = NULL,
               period  = '1985_2015',
               time    = 'past',
               big_cnt = FALSE,
               ncores  = 10)
}

gcmList <- c("ipsl_cm5a_mr","miroc_esm_chem","ncc_noresm1_m")
periodList <- c('2021_2045','2041_2065')
for(i in 1:length(countyList)){
  for(gcm in gcmList){
    for(period in periodList){
      calc_indices(country = 'Kenya',
                   county  = countyList[i],
                   iso3c   = 'KEN',
                   adm_lvl = 2,
                   seasons = list(s1 = 2:6, s2 = 7:12), # Seasons manually defined
                   n_ssns  = NULL,    # 2-seasons automatically defined
                   n_wtts  = NULL,  # 100-wettest days
                   gcm     = gcm,
                   period  = period,
                   time    = 'future',
                   big_cnt = FALSE,
                   ncores  = 10)
    }
  }
}
