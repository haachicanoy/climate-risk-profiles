### Do time series graphs
### H. Achicanoy & A. Esquivel
### Alliance Bioversity-CIAT, 2020

# R options
options(warn = -1, scipen = 999)

# Load libraries
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, fst))

# Paths
OSys <- Sys.info()[1]
root <<- switch(OSys,
                'Linux'   = '/home/jovyan/work/cglabs',
                'Windows' = '//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles')

time_series_plot <- function(country = 'Kenya', county = 'Vihiga'){
  
  country <<- country
  county  <<- county
  cat('>>> Load indices for the historical period\n')
  
  
  past <- fst::fst(paste0(root,"/results/",country,"/past/",county,"_1985_2015_corrected.fst")) %>% data.frame
  
  cat('>>> Load indices for the future period\n')
  futDir  <- paste0(root,'/results/',country,'/future')
  fut_fls <- list.files(futDir, pattern = paste0('^',county,'_[0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9]_corrected.fst'), recursive = T)
  fut_fls <- paste0(futDir,'/',fut_fls)
  future  <- fut_fls %>%
    purrr::map(.f = function(x){df <- fst::fst(x) %>% data.frame; return(df)}) %>%
    do.call(rbind, .)
  
  # Join data
  df <- rbind(past, future)
  # =----------------------------------
  #df <- df %>% filter(id %in% crd$id)
  # =----------------------------------
  if('semester' %in% colnames(df)){
    colnames(df)[which(colnames(df) == 'semester')] <- 'season'
  }
  if('1' %in% df$season){gsub('1', 's1', x = df$season, fixed = T)}
  if('2' %in% df$season){gsub('2', 's2', x = df$season, fixed = T)}
  
  cat('>>> Prepare climatology-based indices: CDD, P5D, P95, NT35, NDWS\n')
  df1_ <- df %>% dplyr::select(CDD,P5D,P95,NT35, ndws, everything() ) %>% #drop_na() %>%
    dplyr::select(year,season,CDD:ndws) %>%
    tidyr::pivot_longer(cols = 'CDD':'ndws', names_to = 'Indices', values_to = 'Value') %>%
    dplyr::group_split(Indices)
  
  cat('>>> Prepare water balance-based indices: SLGP and LGP\n')
  df2_ <- df %>% dplyr::select(gSeason,SLGP, LGP,  everything() ) %>% 
    dplyr::select(year,gSeason:LGP) %>%
    tidyr::pivot_longer(cols = 'SLGP':'LGP', names_to = 'Indices', values_to = 'Value') %>%
    dplyr::group_split(Indices)
  
  # Output folder
  outDir <- paste0(root,'/results/',country,'/graphs/',tolower(county),'/time_series')
  if(!dir.exists(outDir)){dir.create(outDir, recursive = T)}
  
  cat('>>> Creating graphs for climatology-based indicators\n')
  df1_ %>%
    purrr::map(.f = function(tbl){
      idx <- tbl$Indices %>% unique
      df_summ <- tbl %>%
        dplyr::group_by(year,season) %>%
        dplyr::summarise(n      = n(),
                         mean   = mean(Value, na.rm = T),
                         sd     = sd(Value, na.rm = T)) %>%
        dplyr::mutate(sem       = sd/sqrt(n-1),
                      CI_lower  = mean + qt((1-0.95)/2, n - 1) * sem,
                      CI_upper  = mean - qt((1-0.95)/2, n - 1) * sem,
                      Serie     = ifelse(year %in% 1985:2015, 'Past','Fut'),
                      Year      = as.Date(ISOdate(year, 1, 1)),
                      season    = as.factor(season))
      if(length(unique(df_summ$season)) == 1){
        sem.labs <- as_labeller(c('s1' = 'S:1'))
      } else {
        sem.labs <- as_labeller(c('s1' = 'S:1', 's2' = 'S:2', 's3' = 'S:3','s4' = 'S:4'))
      }
      df_summ_ <- df_summ %>%
        dplyr::group_by(season) %>%
        dplyr::group_split(season)
      
      1:length(df_summ_) %>%
        purrr::map(.f = function(i){
          df_summ2 <- df_summ_[[i]]
          plt      <-   df_summ2 %>%
            dplyr::filter(Serie == 'Past') %>%
            ggplot2::ggplot(aes(x = Year, y = mean, colour = season)) +
            ggplot2::geom_line() +
            ggplot2::scale_x_date(date_labels = "%Y", date_breaks = '10 years', limits = as.Date(c('1985-01-01', '2065-01-01'))) +
            #ggplot2::ylim(min(df_summ2$mean)-5, max(df_summ2$mean)+5) +
            ggplot2::geom_line(data = df_summ2 %>% dplyr::filter(Serie == 'Fut'), aes(x = Year, y = mean, colour = season)) +
            ggplot2::geom_ribbon(data = df_summ2 %>% dplyr::filter(Serie == 'Fut'), aes(ymin = CI_lower, ymax = CI_upper, fill = season), color = "grey70", alpha = 0.4) +
            ggplot2::geom_smooth(data = df_summ2, method = "loess", color = "blue", alpha = 0.8, se = FALSE) +
            ggplot2::labs(title    = tbl$Indices %>% unique,
                          subtitle = paste0(country,", ",county)) +
            ggplot2::theme_bw() +
            ggplot2::theme(axis.text       = element_text(size = 17),
                           axis.title      = element_text(size = 20),
                           legend.text     = element_text(size = 17),
                           legend.title    = element_text(size = 20),
                           plot.title      = element_text(size = 20),
                           plot.subtitle   = element_text(size = 17),
                           strip.text.x    = element_text(size = 17),
                           legend.position = 'none') +
            ggplot2::facet_wrap(~season, labeller = labeller(season = sem.labs))
          if(idx == 'CDD'){plt <- plt + ggplot2::ylab('CDD (days)')}
          if(idx == 'P5D'){plt <- plt + ggplot2::ylab('P5D (mm)')}
          if(idx == 'P95'){plt <- plt + ggplot2::ylab('P95 (mm/day)')}
          if(idx == 'NT35'){plt <- plt + ggplot2::ylab('NT35 (days)')}
          if(idx == 'ndws'){plt <- plt + ggplot2::ylab('ndws (days)')}
          ggplot2::ggsave(filename = paste0(outDir,'/ts_',tbl$Indices %>% unique,'_season_',i,'.png'), plot = plt, device = "png", width = 12, height = 6, units = "in")
        })
      
      return(cat('Graphs done\n'))
    })
  
  cat('>>> Creating graphs for water balance-based indicators\n')
  df2_ %>%
    purrr::map(.f = function(tbl){
      idx <- tbl$Indices %>% unique
      tbl <- tbl %>%
        tidyr::drop_na()
      df_summ <- tbl %>%
        dplyr::group_by(year,gSeason) %>%
        dplyr::summarise(n      = dplyr::n(),
                         mean   = mean(Value, na.rm = T),
                         sd     = sd(Value, na.rm = T)) %>%
        dplyr::mutate(sem       = sd/sqrt(n-1),
                      CI_lower  = mean + qt((1-0.95)/2, n - 1) * sem,
                      CI_upper  = mean - qt((1-0.95)/2, n - 1) * sem,
                      Serie     = ifelse(year %in% 1985:2015, 'Past','Fut'),
                      Year      = as.Date(ISOdate(year, 1, 1)),
                      gSeason  = as.factor(gSeason))
      if(length(unique(df_summ$gSeason)) == 1){
        sem.labs <- 'S:1'
        names(sem.labs) <- '1'
      } else {
        sem.labs <- c('S:1','S:2', 'S:3')
        names(sem.labs) <- c('1','2','3')
      }
      df_summ_ <- df_summ %>%
        dplyr::group_by(gSeason) %>%
        dplyr::group_split(gSeason)
      
      1:length(df_summ_) %>%
        purrr::map(.f = function(i){
          df_summ2 <- df_summ_[[i]]
          plt      <- df_summ2 %>%
            dplyr::filter(Serie == 'Past') %>%
            ggplot2::ggplot(aes(x = Year, y = mean, colour = gSeason)) +
            ggplot2::geom_line() +
            ggplot2::scale_x_date(date_labels = "%Y", date_breaks = '10 years', limits = as.Date(c('1985-01-01', '2065-01-01'))) +
            #ggplot2::ylim(min(df_summ2$mean)-5, max(df_summ2$mean)+5) +
            ggplot2::geom_line(data = df_summ2 %>% dplyr::filter(Serie == 'Fut'), aes(x = Year, y = mean, colour = gSeason)) +
            ggplot2::geom_ribbon(data = df_summ2 %>% dplyr::filter(Serie == 'Fut'), aes(ymin = CI_lower, ymax = CI_upper, fill = gSeason), color = "grey70", alpha = 0.4) +
            ggplot2::geom_smooth(data = df_summ2, method = "loess", color = "blue", alpha = 0.8, se = FALSE) +
            ggplot2::labs(title    = tbl$Indices %>% unique,
                          subtitle = paste0(country,", ",county)) +
            ggplot2::theme_bw() +
            ggplot2::theme(axis.text       = element_text(size = 17),
                           axis.title      = element_text(size = 20),
                           legend.text     = element_text(size = 17),
                           legend.title    = element_text(size = 20),
                           plot.title      = element_text(size = 20),
                           plot.subtitle   = element_text(size = 17),
                           strip.text.x    = element_text(size = 17),
                           legend.position = 'none') +
            ggplot2::facet_wrap(~gSeason, labeller = labeller(gSeason = sem.labs))
          if(idx == 'SLGP'){plt <- plt + ggplot2::ylab('SLGP (Day of\nthe year)')}
          if(idx == 'LGP'){plt <- plt + ggplot2::ylab('LGP (days)')}
          ggplot2::ggsave(filename = paste0(outDir,'/ts_',tbl$Indices %>% unique,'_season_',i,'.png'), plot = plt, device = "png", width = 12, height = 6, units = "in")
        })
      
      return(cat('Graphs done\n'))
    })
  
  cat('>>> Graphics done.\n')
}




# =--------------------------------------------------

country <- 'Kenya'
count_i  <- c('Siaya', 'Bungoma', 'Kakamega', 'Nyandarua')

for(i in 1:length(count_i)){
  time_series_plot(country = country, county = count_i[i])}




# =-----------------------------------------
# Tablas...
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

### Do time series graphs
### H. Achicanoy & A. Esquivel
### Alliance Bioversity-CIAT, 2020


ts_table <- function(country = 'Benin', county = 'Borgou'){
  
  country <<- country
  county  <<- county
  cat('>>> Load indices for the historical period\n')
  # past <- fst::fst(paste0(root,"/results/",country,"/past/",county,"_1985_2015_corrected_idw.fst")) %>% data.frame
  past <- fst::fst(paste0(root,"/results/",country,"/past/",county,"_1985_2015_corrected.fst")) %>% data.frame
  
  cat('>>> Load indices for the future period\n')
  futDir  <- paste0(root,'/results/',country,'/future')
  # fut_fls <- list.files(futDir, pattern = paste0('^',county,'_[0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9]_corrected_idw.fst'), recursive = T)
  fut_fls <- list.files(futDir, pattern = paste0('^',county,'_[0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9]_corrected.fst'), recursive = T)
  fut_fls <- paste0(futDir,'/',fut_fls)
  future  <- fut_fls %>%
    purrr::map(.f = function(x){df <- fst::fst(x) %>% data.frame; return(df)}) %>%
    do.call(rbind, .)
  
  # Join data
  df <- rbind(past, future)
  # =----------------------------------
  #df <- df %>% filter(id %in% crd$id)
  # =----------------------------------
  if('semester' %in% colnames(df)){
    colnames(df)[which(colnames(df) == 'semester')] <- 'season'
  }
  if('1' %in% df$season){gsub('1', 's1', x = df$season, fixed = T)}
  if('2' %in% df$season){gsub('2', 's2', x = df$season, fixed = T)}
  
  cat('>>> Prepare climatology-based indices: CDD, P5D, P95, NT35, NDWS\n')
  df1_ <- df %>% dplyr::select(CDD,P5D,P95,NT35, ndws, everything() ) %>%
    dplyr::select(year,season,CDD:ndws) %>%
    tidyr::pivot_longer(cols = 'CDD':'ndws', names_to = 'Indices', values_to = 'Value') %>%
    dplyr::group_split(Indices)
  
  cat('>>> Prepare water balance-based indices: SLGP and LGP\n')
  df2_ <- df %>% dplyr::select(gSeason,SLGP, LGP,  everything() ) %>%
    dplyr::select(year,gSeason,SLGP, LGP) %>%
    mutate(gSeason = case_when(gSeason == 1 ~ 's1', gSeason == 2 ~ 's2', TRUE ~ as.character(gSeason))) %>%
    tidyr::pivot_longer(cols = 'SLGP':'LGP', names_to = 'Indices', values_to = 'Value') %>%
    dplyr::group_split(Indices)
  
  # Output folder
  # outDir <- paste0(root,'/results/',country,'/graphs/',tolower(county),'/time_series')
  # if(!dir.exists(outDir)){dir.create(outDir, recursive = T)}
  
  cat('>>> Creating tables for climatology-based indicators\n')
  st_indx <- df1_ %>%
    purrr::map(.f = function(tbl){
      idx <- tbl$Indices %>% unique
      df_summ <- tbl %>%
        dplyr::group_by(year,season) %>%
        dplyr::summarise(n      = n(),
                         mean   = mean(Value, na.rm = T),
                         sd     = sd(Value, na.rm = T)) %>%
        dplyr::mutate(sem       = sd/sqrt(n-1),
                      CI_lower  = mean + qt((1-0.95)/2, n - 1) * sem,
                      CI_upper  = mean - qt((1-0.95)/2, n - 1) * sem,
                      Serie     = ifelse(year %in% 1985:2015, 'Past','Fut'),
                      Year      = as.Date(ISOdate(year, 1, 1)),
                      season    = as.factor(season))
      if(length(unique(df_summ$season)) == 1){
        sem.labs <- 'S:1'
        names(sem.labs) <- 's1'
      } else {
        sem.labs <- c('S:1','S:2')
        names(sem.labs) <- c('s1','s2')
      }
      
      df_summ_ <- df_summ %>% dplyr::arrange(season) %>% 
        dplyr::select(season, year, Serie, mean, CI_lower, CI_upper) %>% 
        setNames(c('season', 'year', 'time', idx, glue::glue('{idx}_CI_lower'), glue::glue('{idx}_CI_upper')))
      
      return(df_summ_)
    }) %>% 
    reduce(full_join)
  
  max_lvl <- unique(st_indx$season)
  
  cat('>>> Creating graphs for water balance-based indicators\n')
  sd_indx <- df2_ %>%
    purrr::map(.f = function(tbl){
      idx <- tbl$Indices %>% unique
      tbl <- tbl %>%
        tidyr::drop_na()
      df_summ <- tbl %>%
        dplyr::group_by(year,gSeason) %>%
        dplyr::summarise(n      = dplyr::n(),
                         mean   = mean(Value, na.rm = T),
                         sd     = sd(Value, na.rm = T)) %>%
        dplyr::mutate(sem       = sd/sqrt(n-1),
                      CI_lower  = mean + qt((1-0.95)/2, n - 1) * sem,
                      CI_upper  = mean - qt((1-0.95)/2, n - 1) * sem,
                      Serie     = ifelse(year %in% 1985:2015, 'Past','Fut'),
                      Year      = as.Date(ISOdate(year, 1, 1)),
                      gSeason  = as.factor(gSeason))
      if(length(unique(df_summ$gSeason)) == 1){
        sem.labs <- 'S:1'
        names(sem.labs) <- '1'
      } else {
        sem.labs <- c('S:1','S:2')
        names(sem.labs) <- c('1','2')
      }
      # df_summ_ <- df_summ %>% dplyr::group_by(gSeason) %>% dplyr::group_split(gSeason)
      
      df_summ_ <- df_summ %>% dplyr::arrange(gSeason) %>% filter(gSeason %in% max_lvl) %>% 
        dplyr::select(gSeason, year, Serie, mean, CI_lower, CI_upper) %>% 
        setNames(c('season', 'year', 'time', idx, glue::glue('{idx}_CI_lower'), glue::glue('{idx}_CI_upper')))
      
      return(df_summ_)
    }) %>% 
    reduce(full_join)
  
  all_idx <- full_join(st_indx, sd_indx) %>%  
    mutate(Country = country, County = county) %>% 
    dplyr::select(Country, County, season, year, time, everything()) %>%  drop_na(season)
  
  cat('>>> Table\n')
  return(all_idx)}

ts_tbl <- list()
for(i in 1:length(count_i)){
  ts_tbl[[i]] <- ts_table(country = country, county = count_i[i])} 

ts_tbl <- ts_tbl %>% bind_rows()

outDir <- paste0(root,'/results/',country,'/graphs/', country, '.csv')
write_csv(x = ts_tbl, path = outDir)

