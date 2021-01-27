### Do climathology graph - improved version
### A. Esquivel, H. Achicanoy
### CIAT, 2020

# R options
options(warn = -1, scipen = 999)

# Load libraries
suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyverse, lubridate, fst))

# Paths
OSys <- Sys.info()[1]
root <<- switch(OSys,
                'Linux'   = '/home/jovyan/work/cglabs',
                'Windows' = '//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles')

source(paste0(root, "/scripts/indices.R"))

do_climatology <- function(country = 'Kenya',
                           county  = 'Vihiga',
                           seasons = TRUE, # Climatology without any season (manual or automatic)
                           manual  = NULL, # Seasons defined manually e.g. list(s1 = c(11:12,1:4)
                           auto    = list(n = 2)) # Seasons defined by algorithm/automatically
{
  
  input1 <- paste0(root, "/data/observational_data/",tolower(country),"/",tolower(county),".fst")
  input2 <- paste0(root, "/data/observational_data/",tolower(country),"/",tolower(county),"_prec_temp.fst")
  if(!file.exists(input1)){
    site <- fst::read_fst(input2)
    site <- site %>%
      tidyr::nest(Climate = c('id','Date','prec','tmax','tmin')) %>%
      dplyr::rename(id = 'id1') %>%
      dplyr::select(id, everything(.))
  } else {
    site <- fst::read_fst(input1)
    site <- site %>%
      tidyr::nest(Climate = c('id','Date','prec','tmax','tmin','srad')) %>%
      dplyr::rename(id = 'id1') %>%
      dplyr::select(id, everything(.))
  }
  if(nrow(site) > 100){
    set.seed(1235)
    smpl <- sample(x = 1:nrow(site), size = 100, replace = F) %>% sort
  } else {
    smpl <- 1:nrow(site)
  }
  
  all_clmtlgy <- 1:length(smpl) %>%
    purrr::map(.f = function(i){
      clmtlgy <- site[smpl[i],] %>%
        dplyr::pull('Climate') %>%
        .[[1]] %>% 
        dplyr::mutate(Year  = lubridate::year(lubridate::as_date(Date)),
                      Month = lubridate::month(lubridate::as_date(Date))) %>%
        dplyr::group_by(Year, Month) %>%
        dplyr::summarise(Prec = sum(prec, na.rm = T),
                         Tmin = mean(tmin, na.rm = T),
                         Tmax = mean(tmax, na.rm = T)) %>%
        dplyr::group_by(Month) %>%
        dplyr::summarise(Prec = mean(Prec, na.rm = T),
                         Tmin = mean(Tmin, na.rm = T),
                         Tmax = mean(Tmax, na.rm = T))
      return(clmtlgy)
    }) %>%
    dplyr::bind_rows()
  
  avrgs <- all_clmtlgy %>%
    dplyr::group_by(Month) %>%
    dplyr::summarise(Prec = mean(Prec, na.rm = T),
                     Tmin = mean(Tmin, na.rm = T),
                     Tmax = mean(Tmax, na.rm = T))
  
  rlc <- mean(avrgs$Prec) / mean(avrgs$Tmax) * 2
  cols <- c("Tmin" = "blue", "Tmax" = "red")
  
  outDir <- paste0(root,'/results/',country,'/graphs/',tolower(county),'/')
  if(!dir.exists(outDir)){dir.create(outDir, recursive = T)}
  
  plt <- all_clmtlgy %>%
    dplyr::group_by(Month) %>%
    dplyr::summarise(Prec = mean(Prec, na.rm = T),
                     Tmin = mean(Tmin, na.rm = T),
                     Tmax = mean(Tmax, na.rm = T)) %>%
    ggplot2::ggplot(aes(x = Month, y = Prec)) +
    ggplot2::geom_bar(stat="identity", fill = 'lightblue') +
    ggplot2::xlab('Month') +
    ggplot2::ylab('Precipitation (mm)') +
    ggplot2::xlim(0,13)+
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(breaks = 1:12)  +
    ggplot2::geom_line(aes(x = Month, y = Tmin*rlc, colour = 'blue'), size = 1.2) +
    ggplot2::geom_line(aes(x = Month, y = Tmax*rlc, colour = 'red'), size = 1.2) +
    #ggplot2::scale_color_identity(guide = 'legend') +
    ggplot2::theme(axis.text       = element_text(size = 17),
                   axis.title      = element_text(size = 20),
                   legend.text     = element_text(size = 17),
                   legend.title    = element_text(size = 20),
                   plot.title      = element_text(size = 20),
                   plot.subtitle   = element_text(size = 17),
                   strip.text.x    = element_text(size = 17),
                   legend.position = "top") +
    ggplot2::scale_y_continuous(sec.axis = sec_axis(~./rlc, name = expression(atop('Temperature('*~degree*C*')')))) +
    ggplot2::scale_colour_manual(name = "", labels = c("Tmin", "Tmax"), values = c("#00ace6", "#b30000")) 
  
  if(!seasons){
    
    cat(paste0('>>> Create climatology graph for: ',county,'-',country,'without seasons\n'))
    ggplot2::ggsave(filename = paste0(root,'/results/',country,'/graphs/',tolower(county),'/',tolower(county),'_climatology.png'), plot = plt, device = "png", width = 12, height = 6, units = "in")
    
  } else {
    
    cat(paste0('>>> Create climatology graph for: ',county,'-',country,'with seasons plotted\n'))
    
    if(!is.null(auto)){
      
      if(auto$n == 1){
        seasonsInfo <- 1:length(smpl) %>%
          purrr::map(.f = function(i){
            info <- site[smpl[i],] %>%
              dplyr::pull('Climate') %>%
              .[[1]] %>% 
              dplyr::mutate(Year  = lubridate::year(lubridate::as_date(Date)),
                            Month = lubridate::month(lubridate::as_date(Date))) %>%
              dplyr::group_by(Year, add = T) %>%
              dplyr::group_split() %>%
              purrr::map(., .f = function(tbl){
                SummDays <- rsum.lapply(x = tbl$prec, n = 150)
                WetDays  <- SummDays[which.max(cumulative.r.sum(SummDays))]
                WetDays  <- WetDays %>% purrr::map(2) %>% unlist %>% data.frame()
                return(WetDays)
              }) %>%
              dplyr::bind_cols()
            return(info)
          }) %>%
          dplyr::bind_cols()
        tbl_summary <- data.frame(Season = 1,
                                  DIni   = round(rowMeans(seasonsInfo, na.rm = T)[1]),
                                  DEnd   = round(rowMeans(seasonsInfo, na.rm = T)[150])) %>%
          dplyr::mutate(Start = as.Date(DIni,'2000-01-01'),
                        End   = as.Date(DEnd,'2000-01-01'),
                        sMnth = Start %>% ymd() %>% { month(.) + day(.) / days_in_month(.) },
                        eMnth = End %>% ymd() %>% { month(.) + day(.) / days_in_month(.) })
      } else {
        if(auto$n == 2){
          seasonsInfo <- 1:length(smpl) %>%
            purrr::map(.f = function(i){
              info <- site[smpl[i],] %>%
                dplyr::pull('Climate') %>%
                .[[1]] %>% 
                dplyr::mutate(Year     = lubridate::year(lubridate::as_date(Date)),
                              Month    = lubridate::month(lubridate::as_date(Date)),
                              Semester = ifelse(Month %in% 1:6, "1", "2")) %>%
                dplyr::group_by(Semester, add = T) %>%
                dplyr::group_split() %>%
                purrr::map(., .f = function(tbl){
                  info2 <- tbl %>%
                    dplyr::group_by(Year, add = T) %>%
                    dplyr::group_split() %>%
                    purrr::map(., .f = function(tbl2){
                      SummDays <- rsum.lapply(x = tbl2$prec, n = 150)
                      WetDays  <- SummDays[which.max(cumulative.r.sum(SummDays))]
                      WetDays  <- WetDays %>% purrr::map(2) %>% unlist %>% data.frame()
                      return(WetDays)
                    }) %>%
                    dplyr::bind_cols()
                })
              return(info)
            })
          seasonsInfo1 <- seasonsInfo %>% purrr::map(1) %>% dplyr::bind_cols()
          seasonsInfo2 <- seasonsInfo %>% purrr::map(2) %>% dplyr::bind_cols()
          tbl_summary <- data.frame(Season = c('1','2'),
                                    DIni   = c(round(rowMeans(seasonsInfo1, na.rm = T)[1]),
                                               round(rowMeans(seasonsInfo2, na.rm = T)[1])+182),
                                    DEnd   = c(round(rowMeans(seasonsInfo1, na.rm = T)[150]),
                                               round(rowMeans(seasonsInfo2, na.rm = T)[150])+182)) %>%
            dplyr::mutate(Start = as.Date(DIni,'2000-01-01'),
                          End   = as.Date(DEnd,'2000-01-01'),
                          sMnth = Start %>% ymd() %>% { month(.) + day(.) / days_in_month(.) },
                          eMnth = End %>% ymd() %>% { month(.) + day(.) / days_in_month(.) })
          
          alphas <- c(.3,.1)
          seasonL <- c('S:1','S:2')
          if(exists('tbl_summary')){
            for(i in 1:nrow(tbl_summary)){
              plt <- plt +
                ggplot2::annotate("rect", xmin=tbl_summary$sMnth[i], xmax=tbl_summary$eMnth[i], ymin=-Inf, ymax=Inf, alpha=alphas[i], fill='forestgreen') +
                ggplot2::annotate("text", x=(tbl_summary$sMnth[i]+tbl_summary$eMnth[i])/2, y = max(avrgs$Prec)+50, label = seasonL[i], size = 10)
            }
          }
        }
      }
      ggplot2::ggsave(filename = paste0(root,'/results/',country,'/graphs/',tolower(county),'/',tolower(county),'_climatology_auto_seasons.png'), plot = plt, device = "png", width = 12, height = 6, units = "in")
    } else {
      if(!is.null(manual)){
        for(i in 1:length(manual)){
          if(sum(diff(manual[[i]]) < 0) > 0){
            # plt <- plt +
            #   ggplot2::annotate("rect",
            #                     xmin  = manual[[i]][1]-.5,
            #                     xmax  = manual[[i]][which(diff(manual[[i]]) < 0)]+.5,
            #                     ymin  = -Inf,
            #                     ymax  = Inf,
            #                     alpha =.3,
            #                     fill  = "forestgreen")
            # plt <- plt +
            #   ggplot2::annotate("rect",
            #                     xmin  = manual[[i]][which(diff(manual[[i]]) < 0)+1]-.5,
            #                     xmax  = manual[[i]][length(manual[[i]])]+.5,
            #                     ymin  = -Inf,
            #                     ymax  = Inf,
            #                     alpha =.3,
            #                     fill  = "forestgreen")
            
            
            plt <- plt + 
              geom_vline(xintercept = manual[[i]][1]-.5, colour = 'forestgreen', size=1) + 
              geom_vline(xintercept = manual[[i]][length(manual[[i]])]+.5,  colour = 'forestgreen', size=1)
            
            
            
          } else {
            
            
            if(length(manual) == 1){
              
              plt <- plt + 
                geom_vline(xintercept = manual[[i]][1]-.5, colour = 'forestgreen', size=1) + 
                geom_vline(xintercept = manual[[i]][length(manual[[i]])]+.5,  colour = 'forestgreen', size=1)
              
              
            }else{
              colors_vct <- c('forestgreen','#b37b53')
              # plt <- plt +
              #   ggplot2::annotate("rect",
              #                     xmin  = manual[[i]][1]-.5,
              #                     xmax  = manual[[i]][length(manual[[i]])]+.5,
              #                     ymin  = -Inf,
              #                     ymax  = Inf,
              #                     alpha =.3,
              #                     colour  = colors_vct[i], 
              #                     fill =  NULL)
              
              linetype_name <- ifelse( i ==  1, "solid", "dashed")
              
              plt <- plt + 
                geom_vline(xintercept = manual[[i]][1]-.5, colour = colors_vct[i], size=1, linetype = linetype_name) + 
                geom_vline(xintercept = manual[[i]][length(manual[[i]])]+.5,  colour = colors_vct[i], size=1, linetype = linetype_name)
              
              
            }
          }
        }
      }
      ggplot2::ggsave(filename = paste0(root,'/results/',country,'/graphs/',tolower(county),'/',tolower(county),'_climatology_manual_seasons.png'), plot = plt, device = "png", width = 12, height = 6, units = "in")
    }
    
  }
  return(cat('Climatology graph correctly created\n'))
}


## Coding examples
# # No seasons
# do_climatology(country = 'Senegal',
#                county  = 'Kaffrine',
#                seasons = FALSE, # Should be FALSE
#                manual  = NULL,  # Should be NULL
#                auto    = NULL)  # Should be NULL

# # Manual seasons
do_climatology(country = 'Togo',
               county  = 'Maritime',
               seasons = TRUE,           # Should be TRUE
               manual  = list(s1 = 2:6, s2 = 7:12), # Should be something like e.g. list(s1 = c(11:12,1:4)
               auto    = NULL)           # Should be NULL

# # Auto seasons
# do_climatology(country = 'Kenya',
#                county  = 'Vihiga',
#                seasons = TRUE,        # Should be TRUE
#                manual  = NULL,        # Should be NULL
#                auto    = list(n = 2)) # Should be something like e.g. list(n = 2)