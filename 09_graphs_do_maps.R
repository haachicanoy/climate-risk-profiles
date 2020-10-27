# setwd('//dapadfs/workspace_cluster_8/climateriskprofiles/results/all_countrys_maps/Kenya_P/')
# raster::getData('GADM', country='KEN', level=2)
# Ken <- readRDS('//dapadfs/workspace_cluster_8/climateriskprofiles/results/all_countrys_maps/Kenya_P/gadm36_KEN_2_sp.rds')  


# Maps for 
# H. Achicanoy & A. Esquivel
# Alliance Bioversity-CIAT
# July - 2020. 


rm(list = ls())
gc(reset = TRUE)

# =--------------------
# Packages 
options(warn = -1, scipen = 999)

suppressMessages(library(pacman))
suppressMessages(pacman::p_load(tidyr, dplyr, tibble, ggplot2, raster, ncdf4, sf, lubridate, glue, cowsay, fst, ggspatial, vroom, sp, compiler))
# =--------------------
# =----------------------------------
# Identificacion de pixel para ETH
# =----------------------------------
country <- 'Kenya'
count_i  <-    c('Kakamega')

iso3c <- 'KEN'
Big <- 'N'
adm_lvl <- 1



for(i in 1:length(count_i)){
  county  <- count_i[i]
  # -----------------------------------
  C_shp <- county
  co <- tolower(country)
  country1 <- co
  # =--------------------
  
  # Ruta Principal para guardados: 
  root <- '//dapadfs/workspace_cluster_8/climateriskprofiles/'
  
  # Load county shapefile
  if(country == 'India'){
    # India 
    country1 <- raster::shapefile(glue::glue('//dapadfs/workspace_cluster_8/climateriskprofiles/results/India/states/Admin2.shp'))
    shp <- raster::shapefile(glue::glue('//dapadfs/workspace_cluster_8/climateriskprofiles/results/India/states/Admin2.shp'))
    shp <- shp[shp@data$ST_NM %in% C_shp,]
    plot(shp)
    shp@data$ISO <- iso3c
  }else{
    country1 <- readRDS(glue::glue('{root}data/shps/shps_from_R/{country}/gadm36_{iso3c}_{adm_lvl}_sp.rds')) 
    shp <- readRDS(glue::glue('{root}data/shps/shps_from_R/{country}/gadm36_{iso3c}_{adm_lvl}_sp.rds')) 
    shp@data$NAME_1 <- iconv(shp@data$NAME_1,from="UTF-8",to="ASCII//TRANSLIT")
    shp@data$NAME_1 <- case_when(shp@data$NAME_1 == 'Trans Nzoia' ~ 'Trans-Nzoia',
                                 shp@data$NAME_1 == "Murang'a" ~ 'Muranga', 
                                 TRUE ~ shp@data$NAME_1)
    shp <- shp[shp@data$NAME_1 %in% C_shp,]
    plot(shp)
    shp@data$ISO <- iso3c
  }
  
  # Load id coords
  crd <- vroom('//dapadfs/workspace_cluster_8/climateriskprofiles/data/id_all_country.csv')
  crd <- crd %>%
    dplyr::filter(Country == country)
  pnt <- crd %>% dplyr::select(x,y) %>% sp::SpatialPoints(coords = .)
  crs(pnt) <- crs(shp)
  # Filter coordinates that are present in the county
  pnt <- sp::over(pnt, shp) %>% data.frame %>% dplyr::select(ISO) %>% complete.cases() %>% which()
  crd <- crd[pnt,]
  crd <<- crd
  
  
  ##### =-----------------------------------------------------------
  # ### =-----------------------------------------------------------
  # =------------ Funciones...
  
  do_climate_maps <- function(data_split){
    
    ISO3 <- unique(data_split$ISO3); county <- unique(data_split$county)
    semester <-  unique(data_split$semester) ;  Country <- unique(data_split$Country)
    
    # Conditions...
    if(Big =='N'){
      
      median_data <- data_split %>%
        dplyr::group_by(id, county, Country, x, y, ISO3, semester, time) %>%
        dplyr::select(-year) %>%
        dplyr::summarise_all(mean) %>%  dplyr::ungroup() %>%
        dplyr::mutate(NT35 = round(x = NT35, digits = 0), 
                      ndws = round(ndws,0))
      
    }else if(Big == 'B'){
      # =---------------------------------------------------------
      median_data <- data_split 
      
    } else{print('Change big argument... >.<')}
    
    
    limits <- dplyr::select(median_data, CDD, P5D, P95, NT35, ndws) %>% 
      dplyr::summarise_all(.funs = c('min', 'max'))
    
    
    median_data <- median_data %>% dplyr::filter(time == 'future') %>%
      dplyr::rename('CDD_f' = 'CDD'   , 'P5D_f' = 'P5D'  , 'P95_f'= 'P95' , 'NT35_f' = 'NT35', 
                    'ndws_f' = 'ndws' ) %>%
      dplyr::select(-time) %>%
      dplyr::inner_join(dplyr::filter(median_data , time == 'past') %>% dplyr::select(-time), . ) %>%
      dplyr::mutate(CDD_c = CDD_f - CDD, P5D_c = P5D_f - P5D, P95_c = P95_f - P95, NT35_c = NT35_f -NT35, ndws_c = ndws_f - ndws)
    
    
    # Aqui se hace solo la figura base...
    shp_sf <- shp  %>% sf::st_as_sf()
    country <- country1 %>% sf::st_as_sf()
    xlims <- sf::st_bbox(shp_sf)[c(1, 3)]
    ylims <- sf::st_bbox(shp_sf)[c(2, 4)]
    
    b <- ggplot() +
      geom_sf(data = shp_sf, fill = 'red', color = gray(.1)) +
      geom_sf(data = country, fill = NA, color = gray(.5)) +
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.text.x = element_blank(), axis.text.y = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 7, face = "bold"))
    
    
    ggsave(glue::glue('{path}{Country}/graphs/{tolower(county)}/maps/{county}.png') , width = 8, height = 5.5, dpi = 300)
    
    
    #===---------------------------------------------------------
    #=------------------------------------------------------------
    
    index_a <- 'CDD'
    
    a <- ggplot() +
      geom_tile(data = median_data, aes(x = x, y = y, fill = CDD)) +
      geom_sf(data = country, fill = NA, color = gray(.8)) +
      geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
      coord_sf(xlim = xlims, ylim = ylims) +
      labs(fill = glue::glue('{index_a}\n(days)  '), title = 'Historic', x = 'Longitude', y = 'Latitude') +
      scale_fill_viridis_c(limits = c(round(limits$CDD_min,2)-0.1, round(limits$CDD_max,2)+0.1), 
                           guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
      # ggspatial::annotation_scale(location = "br", width_hint = 0.5) +
      # ggspatial::annotation_north_arrow(location = "br", which_north = "true",
      # pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), # 0.2 # 0.3
      # style = north_arrow_fancy_orienteering) +
      scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
      scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
      theme_bw() +
      theme(legend.position = 'bottom')
    
    
    ggsave(glue::glue('{path}{Country}/graphs/{county}/maps/{index_a}_past_S{semester}.png') , width = 8, height = 5.5, dpi = 300)
    
    
    # =- Lo mismo para futuro...
    
    a1 <-  ggplot() +
      geom_tile(data = median_data, aes(x = x, y = y, fill = CDD_f)) +
      geom_sf(data = country, fill = NA, color = gray(.8)) +
      geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
      coord_sf(xlim = xlims, ylim = ylims) +
      labs(fill = glue::glue('{index_a}\n(days)  '), title = 'Future',x = 'Longitude', y = 'Latitude') +
      scale_fill_viridis_c(limits = c(round(limits$CDD_min, 2)-0.1, round(limits$CDD_max,2) + 0.1), 
                           guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
      # ggspatial::annotation_scale(location = "br", width_hint = 0.5) +
      # ggspatial::annotation_north_arrow(location = "br", which_north = "true",
      #                                   pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), # 0.2 # 0.3
      #                                   style = north_arrow_fancy_orienteering) +
      scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
      scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
      theme_bw() +
      theme(legend.position = 'bottom')
    
    ggsave(glue::glue('{path}{Country}/graphs/{county}/maps/{index_a}_future_S{semester}.png') , width = 8, height = 5.5, dpi = 300)
    
    
    a_d <- ggplot() +
      geom_tile(data = median_data, aes(x = x, y = y, fill = CDD_c)) +
      geom_sf(data = country, fill = NA, color = gray(.8)) +
      geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
      coord_sf(xlim = xlims, ylim = ylims) +
      labs(fill = glue::glue('{index_a}\n(days)  '), title = 'Change', x = 'Longitude', y = 'Latitude') +
      scale_fill_gradient2(low = '#000099', mid = 'white', high = '#A50026', 
                           guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
      # ggspatial::annotation_scale(location = "br", width_hint = 0.5) +
      # ggspatial::annotation_north_arrow(location = "br", which_north = "true",
      # pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), # 0.2 # 0.3
      # style = north_arrow_fancy_orienteering) +
      scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
      scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
      theme_bw() +
      theme(legend.position = 'bottom')
    
    ggsave(glue::glue('{path}{Country}/graphs/{county}/maps/Dif_{index_a}_S{semester}.png') , width = 8, height = 5.5, dpi = 300)
    
    
    png(filename = glue::glue('{path}{Country}/graphs/{county}/maps/Dif_{index_a}_{semester}.png'), width=12.5,height=4.5,units="in", res = 300) # width = 1580, height = 720,
    print(gridExtra::grid.arrange(a, a1, a_d, ncol=3,  
                                  top = glue::glue('{Country}, {county}\nS:{semester}',
                                                   bottom =   "Data source: Alliance Bioversity-CIAT")))
    dev.off()
    
    
    #===---------------------------------------------------------
    #=------------------------------------------------------------
    
    # Siguiente indice...
    index_c <- 'P5D'
    
    c <- ggplot() +
      geom_tile(data = median_data, aes(x = x, y = y, fill = P5D)) +
      geom_sf(data = country, fill = NA, color = gray(.8)) +
      geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
      coord_sf(xlim = xlims, ylim = ylims) +
      labs(fill = glue::glue('{index_c}\n(mm)  '), title = 'Historic',x = 'Longitude', y = 'Latitude') +
      scale_fill_gradientn(colours = blues9, limits = c(round(limits$P5D_min, 2)- 0.1, 
                                                        round(limits$P5D_max, 2)+0.1), 
                           guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
      # ggspatial::annotation_scale(location = "br", width_hint = 0.5) +
      # ggspatial::annotation_north_arrow(location = "br", which_north = "true",
      #                                   pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), # 0.2 # 0.3
      #                                   style = north_arrow_fancy_orienteering) +
      scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
      scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
      theme_bw() +
      theme(legend.position = 'bottom')
    
    ggsave(glue::glue('{path}{Country}/graphs/{county}/maps/{index_c}_past_S{semester}.png') , width = 8, height = 5.5, dpi = 300)
    
    # =- Futuro.
    c1 <- ggplot() +
      geom_tile(data = median_data, aes(x = x, y = y, fill = P5D_f)) +
      geom_sf(data = country, fill = NA, color = gray(.8)) +
      geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
      coord_sf(xlim = xlims, ylim = ylims) +
      labs(fill = glue::glue('{index_c}\n(mm)  '), title = 'Future', x = 'Longitude', y = 'Latitude') +
      scale_fill_gradientn(colours = blues9, limits = c(round(limits$P5D_min, 2)- 0.1, 
                                                        round(limits$P5D_max, 2)+0.1), 
                           guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
      # ggspatial::annotation_scale(location = "br", width_hint = 0.5) +
      # ggspatial::annotation_north_arrow(location = "br", which_north = "true",
      #                                   pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), # 0.2 # 0.3
      #                                   style = north_arrow_fancy_orienteering) +
      scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
      scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
      theme_bw() +
      theme(legend.position = 'bottom')
    
    ggsave(glue::glue('{path}{Country}/graphs/{county}/maps/{index_c}_future_S{semester}.png') , width = 8, height = 5.5, dpi = 300)
    
    
    c_d <- ggplot() +
      geom_tile(data = median_data, aes(x = x, y = y, fill = P5D_c)) +
      geom_sf(data = country, fill = NA, color = gray(.8)) +
      geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
      coord_sf(xlim = xlims, ylim = ylims) +
      labs(fill = glue::glue('{index_c}\n(mm)  '), title = 'Change', x = 'Longitude', y = 'Latitude') +
      scale_fill_gradient2(low = '#A50026', mid = 'white', high = '#000099', 
                           guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
      # ggspatial::annotation_scale(location = "br", width_hint = 0.5) +
      # ggspatial::annotation_north_arrow(location = "br", which_north = "true",
      #                                   pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), # 0.2 # 0.3
      #                                   style = north_arrow_fancy_orienteering) +
      scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
      scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
      theme_bw() +
      theme(legend.position = 'bottom')
    
    ggsave(glue::glue('{path}{Country}/graphs/{county}/maps/Dif_{index_c}_S{semester}.png') , width = 8, height = 5.5, dpi = 300)
    
    
    png(filename = glue::glue('{path}{Country}/graphs/{county}/maps/Dif_{index_c}_{semester}.png') , width=12.5,height=4.5,units="in", res = 300)
    print(gridExtra::grid.arrange(c, c1, c_d, ncol=3,  
                                  top = glue::glue('{Country}, {county}\nS:{semester}',
                                                   bottom =   "Data source: Alliance Bioversity-CIAT")))
    dev.off()
    
    
    #===---------------------------------------------------------
    # =------------
    index_d <- 'P95'
    
    d <- ggplot() +
      geom_tile(data = median_data, aes(x = x, y = y, fill = P95)) +
      geom_sf(data = country, fill = NA, color = gray(.8)) +
      geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
      coord_sf(xlim = xlims, ylim = ylims) +
      labs(fill = glue::glue('{index_d}\n(mm)  '), title = 'Historic', x = 'Longitude', y = 'Latitude') +
      scale_fill_gradientn(colours = blues9, limits = c(round(limits$P95_min,2)-0.1, 
                                                        round(limits$P95_max, 2)+0.1), 
                           guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
      # ggspatial::annotation_scale(location = "br", width_hint = 0.5) +
      # ggspatial::annotation_north_arrow(location = "br", which_north = "true",
      #                                   pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), # 0.2 # 0.3
      #                                   style = north_arrow_fancy_orienteering) +
      scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
      scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
      theme_bw() +
      theme(legend.position = 'bottom')
    
    ggsave(glue::glue('{path}{Country}/graphs/{county}/maps/{index_d}_past_S{semester}.png') , width = 8, height = 5.5, dpi = 300)
    
    # =----
    
    d1 <-  ggplot() +
      geom_tile(data = median_data, aes(x = x, y = y, fill = P95_f)) +
      geom_sf(data = country, fill = NA, color = gray(.8)) +
      geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
      coord_sf(xlim = xlims, ylim = ylims) +
      labs(fill = glue::glue('{index_d}\n(mm)  '), title = 'Future', x = 'Longitude', y = 'Latitude') +
      scale_fill_gradientn(colours = blues9, limits =  c(round(limits$P95_min,2)-0.1, 
                                                         round(limits$P95_max, 2)+0.1), 
                           guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
      # ggspatial::annotation_scale(location = "br", width_hint = 0.5) +
      # ggspatial::annotation_north_arrow(location = "br", which_north = "true",
      #                                   pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), # 0.2 # 0.3
      #                                   style = north_arrow_fancy_orienteering) +
      scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
      scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
      theme_bw() +
      theme(legend.position = 'bottom')
    
    
    ggsave(glue::glue('{path}{Country}/graphs/{county}/maps/{index_d}_future_S{semester}.png') , width = 8, height = 5.5, dpi = 300)
    
    
    
    d_d <- ggplot() +
      geom_tile(data = median_data, aes(x = x, y = y, fill = P95_c)) +
      geom_sf(data = country, fill = NA, color = gray(.8)) +
      geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
      coord_sf(xlim = xlims, ylim = ylims) +
      labs(fill = glue::glue('{index_d}\n(mm)  '), title = 'Change', x = 'Longitude', y = 'Latitude') +
      scale_fill_gradient2(low = '#A50026', mid = 'white', high = '#000099', 
                           guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
      scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
      scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
      theme_bw() +
      theme(legend.position = 'bottom')
    
    ggsave(glue::glue('{path}{Country}/graphs/{county}/maps/Dif_{index_d}_S{semester}.png') , width = 8, height = 5.5, dpi = 300)
    
    
    
    png(filename = glue::glue('{path}{Country}/graphs/{county}/maps/Dif_{index_d}_{semester}.png') , width=12.5,height=4.5,units="in", res = 300)
    print(gridExtra::grid.arrange(d, d1, d_d, ncol=3,  
                                  top = glue::glue('{Country}, {county}\nS:{semester}',
                                                   bottom =   "Data source: Alliance Bioversity-CIAT")))
    dev.off()
    
    
    #===---------------------------------------------------------
    #=------------------------------------------------------------
    
    # =------------
    index_e <- 'NT35'
    
    e <- ggplot() +
      geom_tile(data = median_data, aes(x = x, y = y, fill = NT35)) +
      geom_sf(data = country, fill = NA, color = gray(.8)) +
      geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
      coord_sf(xlim = xlims, ylim = ylims) +
      labs(fill = glue::glue('{index_e}\n(days)  '), title = 'Historic', x = 'Longitude', y = 'Latitude') +
      scale_fill_viridis_c(limits = c(round(limits$NT35_min, 2) - 0.1, round(limits$NT35_max, 2)+0.1), 
                           guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
      # ggspatial::annotation_scale(location = "br", width_hint = 0.5) +
      # ggspatial::annotation_north_arrow(location = "br", which_north = "true",
      #                                   pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), # 0.2 # 0.3
      #                                   style = north_arrow_fancy_orienteering) +
      scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
      scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
      theme_bw() +
      theme(legend.position = 'bottom')
    
    ggsave(glue::glue('{path}{Country}/graphs/{county}/maps/{index_e}_past_S{semester}.png') , width = 8, height = 5.5)
    
    # =------------
    
    e1 <- ggplot() +
      geom_tile(data = median_data, aes(x = x, y = y, fill = NT35_f)) +
      geom_sf(data = country, fill = NA, color = gray(.8)) +
      geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
      coord_sf(xlim = xlims, ylim = ylims) +
      labs(fill = glue::glue('{index_e}\n(days)  '), title = 'Future', x = 'Longitude', y = 'Latitude') +
      scale_fill_viridis_c(limits = c(round(limits$NT35_min, 2) - 0.1, round(limits$NT35_max, 2)+0.1), 
                           guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
      scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
      scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
      theme_bw() +
      theme(legend.position = 'bottom')
    
    ggsave(glue::glue('{path}{Country}/graphs/{county}/maps/{index_e}_future_S{semester}.png') , width = 8, height = 5.5)
    
    
    
    e_d <- ggplot() + geom_tile(data = median_data, aes(x = x, y = y, fill = NT35_c)) +
      geom_sf(data = country, fill = NA, color = gray(.8)) +
      geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
      coord_sf(xlim = xlims, ylim = ylims) +
      labs(fill = glue::glue('{index_e}\n(days)  '), title = 'Change', x = 'Longitude', y = 'Latitude') +
      scale_fill_gradient2(low = '#000099', mid = 'white', high = '#A50026', 
                           guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
      scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
      scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
      theme_bw() +
      theme(legend.position = 'bottom')
    
    ggsave(glue::glue('{path}{Country}/graphs/{county}/maps/Dif_{index_e}_S{semester}.png') , width = 8, height = 5.5)
    
    
    png(filename = glue::glue('{path}{Country}/graphs/{county}/maps/Dif_{index_e}_{semester}.png') , width=12.5,height=4.5,units="in", res = 300)
    print(gridExtra::grid.arrange(e, e1, e_d, ncol=3,  
                                  top = glue::glue('{Country}, {county}\nS:{semester}',
                                                   bottom =   "Data source: Alliance Bioversity-CIAT")))
    dev.off()
    
    
    #===---------------------------------------------------------
    #=------------------------------------------------------------
    cowsay::say('here')
    # =------------
    index_f <- 'ndws'
    
    f <- ggplot() +
      geom_tile(data = median_data, aes(x = x, y = y, fill = ndws)) +
      geom_sf(data = country, fill = NA, color = gray(.8)) +
      geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
      coord_sf(xlim = xlims, ylim = ylims) +
      labs(fill = glue::glue('{index_f}\n(days)  '), title = 'Historic', x = 'Longitude', y = 'Latitude') +
      scale_fill_viridis_c(limits = c(round(limits$ndws_min, 2), round(limits$ndws_max, 2)), 
                           guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
      scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
      scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
      theme_bw() + theme(legend.position = 'bottom')
    
    ggsave(glue::glue('{path}{Country}/graphs/{county}/maps/{index_f}_past_S{semester}.png') , width = 8, height = 5.5)
    
    # =------------
    
    f1 <-  ggplot() +
      geom_tile(data = median_data, aes(x = x, y = y, fill = ndws_f)) +
      geom_sf(data = country, fill = NA, color = gray(.8)) +
      geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
      coord_sf(xlim = xlims, ylim = ylims) +
      labs(fill = glue::glue('{index_f}\n(days)  '), title = 'Future', x = 'Longitude', y = 'Latitude') +
      scale_fill_viridis_c(limits = c(round(limits$ndws_min, 2), round(limits$ndws_max, 2)), 
                           guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
      scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
      scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
      theme_bw() +
      theme(legend.position = 'bottom')
    
    ggsave(glue::glue('{path}{Country}/graphs/{county}/maps/{index_e}_future_S{semester}.png') , width = 8, height = 5.5)
    
    
    
    f_d <- ggplot() + geom_tile(data = median_data, aes(x = x, y = y, fill = ndws_c )) +
      geom_sf(data = country, fill = NA, color = gray(.8)) +
      geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
      coord_sf(xlim = xlims, ylim = ylims) +
      labs(fill = glue::glue('{index_f}\n(days)  '), title = 'Change', x = 'Longitude', y = 'Latitude') +
      scale_fill_gradient2(low = '#000099', mid = 'white', high = '#A50026', 
                           guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
      scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
      scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
      theme_bw() +
      theme(legend.position = 'bottom')
    
    ggsave(glue::glue('{path}{Country}/graphs/{county}/maps/Dif_{index_f}_S{semester}.png') , width = 8, height = 5.5)
    
    
    png(filename = glue::glue('{path}{Country}/graphs/{county}/maps/Dif_{index_f}_{semester}.png') , width=12.5,height=4.5,units="in", res = 300)
    print(gridExtra::grid.arrange(f, f1, f_d, ncol=3,  
                                  top = glue::glue('{Country}, {county}\nS:{semester}',
                                                   bottom =   "Data source: Alliance Bioversity-CIAT")))
    dev.off()
    
    
    
    
    return(median_data)}
  
  # =------------
  # Funcion for run basic maps. 
  do_srad_maps <- function(data_split){
    
    ISO3 <- unique(data_split$ISO3); county <- unique(data_split$county) 
    Country <- unique(data_split$Country)
    Big <- unique(data_split$Big)
    
    # Aqui se hace solo la figura base...
    shp_sf <- shp  %>% sf::st_as_sf()
    pais <- country1 %>% sf::st_as_sf()
    xlims <- sf::st_bbox(shp_sf)[c(1, 3)]
    ylims <- sf::st_bbox(shp_sf)[c(2, 4)]
    
    
    if(Big =='N'){
      # gSeason...
      gSeason_i <- data_split %>% dplyr::select( time, x, y,year, gSeason) %>%
        group_by(time, x, y, year) %>% summarise(gSeason = max(gSeason)) %>%
        ungroup() %>% group_by(time, x, y) %>%
        summarise(gSeason = round(mean(gSeason, na.rm = TRUE), 0)) %>% ungroup()
      
      # =----------------------------------------------------------------------
      # SLGP -- LGP 
      two_index <- data_split %>% dplyr::select( time, x, y, gSeason, SLGP, LGP) %>% 
        group_by(time, x, y, gSeason) %>% 
        summarise_all(.f = function(x){round(mean(x, na.rm = TRUE), 0)}) %>% 
        # ungroup(year) %>% 
        arrange(gSeason) %>% 
        tidyr::nest(data = c('x', 'y','SLGP', 'LGP')) %>% 
        drop_na() %>% filter(gSeason < 3) %>% tidyr::unnest() %>% ungroup()
      
    }else if(Big == 'B'){
      # =---------------------------------------------------------
      # gSeason...
      gSeason_i <- data_split %>% dplyr::select( time, x, y, gSeason) %>%
        group_by( time, x, y) %>% summarise(gSeason = max(gSeason, na.rm =  TRUE)) %>%
        ungroup() %>% group_by(time, x, y) %>%
        summarise(gSeason = round(mean(gSeason, na.rm = TRUE), 0)) %>% ungroup()
      
      # =----------------------------------------------------------------------
      # gSeason...
      two_index <- data_split %>% dplyr::select( time, x, y, gSeason, SLGP, LGP) %>% 
        group_by(time, x, y, gSeason) %>% 
        summarise_all(.f = function(x){round(mean(x, na.rm = TRUE), 0)}) %>% 
        # ungroup(year) %>%  
        arrange(gSeason) %>% 
        tidyr::nest(data = c('x', 'y','SLGP', 'LGP')) %>% 
        drop_na() %>% filter(gSeason < 3) %>% tidyr::unnest() %>% ungroup()
    } else{print('Change big argument... >.<')}
    
    
    # =---------------------------------------------------------------------
    # =---------------------------------------------------------------------
    # gSeason
    
    limits_gs <- dplyr::select(gSeason_i, gSeason) %>% summarise_all(.funs = c('min', 'max'))
    
    # # Primero dejaré hechos los de presente... luego repito los de futuro...
    gs <-  ggplot() +
      geom_tile(data = filter(gSeason_i, time == 'past'), aes(x = x, y = y, fill = gSeason)) +
      geom_sf(data = pais, fill = NA, color = gray(.8)) +
      geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
      coord_sf(xlim = xlims, ylim = ylims) +
      labs(fill = glue::glue('gSeason\n(day)  '), title = 'Historic', x = 'Longitude', y = 'Latitude') +
      scale_fill_gradientn(colours = blues9, limits = as.integer(limits_gs), 
                           guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
      scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
      scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
      theme_bw() +
      theme(legend.position = 'bottom')
    
    
    ggsave(glue::glue('{path}{Country}/graphs/{county}/maps/gSeason_past.png') , width = 8, height = 5.5, dpi = 300)
    
    # # =- Futuro.
    gs_f <- ggplot() +
      geom_tile(data = filter(gSeason_i, time == 'future'), aes(x = x, y = y, fill = gSeason)) +
      geom_sf(data = pais, fill = NA, color = gray(.8)) +
      geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
      coord_sf(xlim = xlims, ylim = ylims) +
      labs(fill = glue::glue('gSeason\n(day)  '), title = 'Future', x = 'Longitude', y = 'Latitude') +
      scale_fill_gradientn(colours = blues9, limits = as.integer(limits_gs), 
                           guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
      scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
      scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
      theme_bw() +
      theme(legend.position = 'bottom')
    
    ggsave(glue::glue('{path}{Country}/graphs/{county}/maps/gSeason_future.png') , width = 8, height = 5.5, dpi = 300)
    
    # =--------
    gSeason_dif <- gSeason_i %>% pivot_wider(names_from = time, values_from = gSeason) %>%
      mutate(gSeason = future - past)
    
    
    c <- ggplot() +
      geom_tile(data = gSeason_dif, aes(x = x, y = y, fill = gSeason)) +
      geom_sf(data = pais, fill = NA, color = gray(.8)) +
      geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
      coord_sf(xlim = xlims, ylim = ylims) +
      labs(fill = glue::glue('gSeason\n(days)  '), title = 'Change', x = 'Longitude', y = 'Latitude') +
      scale_fill_gradient2(low = '#A50026', mid = 'white', high = '#000099', 
                           guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
      scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
      scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
      theme_bw() +
      theme(legend.position = 'bottom')
    
    ggsave(glue::glue('{path}{Country}/graphs/{county}/maps/Dif_gSeason.png') , width = 8, height = 5.5, dpi = 300)
    
    
    png(filename = glue::glue('{path}{Country}/graphs/{county}/maps/all_gSeason.png'), width=12.5,height=4.5,units="in", res = 300)
    print(gridExtra::grid.arrange(gs, gs_f, c, ncol=3,  
                                  top = glue::glue('{Country}, {county}',
                                                   bottom =   "Data source: Alliance Bioversity-CIAT")))
    dev.off()
    
    
    # =---------------------------------------------------------------------
    # =---------------------------------------------------------------------
    # Aqui hay que hacer los filtros de acuerdo al grupo ...
    
    # rows_n <- data_split %>% filter(time == 'past', year == '1985') %>% dplyr::select(x, y, id) %>% unique() %>% nrow()
    
    SLGP_dif <- two_index %>% dplyr::select(-LGP)  %>%
      pivot_wider(names_from = time, values_from = SLGP) %>%
      mutate(SLGP = future - past, gSeason = glue::glue('gSeason = {gSeason}'))
    
    limits_two <-   two_index %>% group_by(gSeason) %>%
      dplyr::select(SLGP, LGP) %>% summarise_all(.funs = c('min', 'max')) %>% 
      ungroup()
    
    gS <- unique(two_index$gSeason)
    
    for(i in gS){
      SLGP_p_1 <- ggplot(filter(two_index %>% mutate(gSeason = glue::glue('gSeason = {gSeason}')),
                                time == 'past', gSeason == glue::glue('gSeason = {i}'))) +
        geom_tile(aes(x = x, y = y, fill = SLGP))  +
        geom_sf(data = pais, fill = NA, color = gray(.8)) +
        geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
        coord_sf(xlim = xlims, ylim = ylims) +
        labs(fill = glue::glue('SLGP\n(Day of\nthe year)  '), 
             title = glue::glue('gSeason = {i}; Historic'),
             x = 'Longitude', y = 'Latitude') +
        scale_fill_gradientn(colours = blues9, limits = c(limits_two$SLGP_min[i], limits_two$SLGP_max[i]), 
                             guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
        # ggspatial::annotation_scale(location = "br", width_hint = 0.5) +
        # ggspatial::annotation_north_arrow(location = "br", which_north = "true",
        # pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), # 0.2 # 0.3
        # style = north_arrow_fancy_orienteering) +
        scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
        scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
        theme_bw() +
        theme(legend.position = 'bottom')
      
      ggsave(glue::glue('{path}{Country}/graphs/{county}/maps/SLGP_past_{i}.png') , width = 8, height = 4, dpi = 300)
      
      
      
      SLGP_f_1 <- ggplot(filter(two_index %>% mutate(gSeason = glue::glue('gSeason = {gSeason}')),
                                time == 'future', gSeason == glue::glue('gSeason = {i}'))) +
        geom_tile(aes(x = x, y = y, fill = SLGP))  +
        geom_sf(data = pais, fill = NA, color = gray(.8)) +
        geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
        coord_sf(xlim = xlims, ylim = ylims) +
        labs(fill = glue::glue('SLGP\n(Day of\nthe year)  '), 
             title = glue::glue('gSeason = {i}; Future'), x = 'Longitude', y = 'Latitude') +
        scale_fill_gradientn(colours = blues9, limits = c(limits_two$SLGP_min[i], limits_two$SLGP_max[i]), 
                             guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
        # ggspatial::annotation_scale(location = "br", width_hint = 0.5) +
        # ggspatial::annotation_north_arrow(location = "br", which_north = "true",
        # pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), # 0.2 # 0.3
        # style = north_arrow_fancy_orienteering) +
        scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
        scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
        theme_bw() +
        theme(legend.position = 'bottom')
      
      
      ggsave(glue::glue('{path}{Country}/graphs/{county}/maps/SLGP_future_{i}.png') , width = 8, height = 4, dpi = 300)
      
      
      d <- ggplot() +
        geom_tile(data = filter(SLGP_dif, gSeason == glue::glue('gSeason = {i}')), 
                  aes(x = x, y = y, fill = SLGP))  +
        geom_sf(data = pais, fill = NA, color = gray(.8)) +
        geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
        coord_sf(xlim = xlims, ylim = ylims) +
        labs(fill = glue::glue('SLGP\n(days)  '), 
             title = glue::glue('gSeason = {i}; Change'),x = 'Longitude', y = 'Latitude') +
        scale_fill_gradient2(low = '#A50026', mid = 'white', high = '#000099', 
                             guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
        # ggspatial::annotation_scale(location = "br", width_hint = 0.5) +
        # ggspatial::annotation_north_arrow(location = "br", which_north = "true",
        # pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), # 0.2 # 0.3
        # style = north_arrow_fancy_orienteering) +
        scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
        scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
        theme_bw() +
        theme(legend.position = 'bottom')
      
      ggsave(glue::glue('{path}{Country}/graphs/{county}/maps/Dif_SLGP_{i}.png') , width = 8, height = 4, dpi = 300)
      
      
      png(filename = glue::glue('{path}{Country}/graphs/{county}/maps/all_SLGP_{i}.png') , width=12.5,height=4.5,units="in", res = 300)
      print(gridExtra::grid.arrange(SLGP_p_1, SLGP_f_1, d, ncol=3,  
                                    top = glue::glue('{Country}, {county}',
                                                     bottom =   "Data source: Alliance Bioversity-CIAT")))
      dev.off()
    }
    
    
    #### =---------------------------------------------------------------------------
    
    # =--------------------------------------------------------------------
    # =--------------------------------------------------------------------
    LGP_dif <- two_index %>% dplyr::select(-SLGP)  %>%
      pivot_wider(names_from = time, values_from = LGP) %>%
      mutate(LGP = future - past, gSeason = glue::glue('gSeason = {gSeason}'))
    
    for(i in gS){
      LGP_p <-  ggplot(filter(two_index %>% mutate(gSeason = glue::glue('gSeason = {gSeason}')),
                              time == 'past', gSeason == glue::glue('gSeason = {i}'))) +
        geom_tile(aes(x = x, y = y, fill = LGP))  +
        geom_sf(data = pais, fill = NA, color = gray(.8)) +
        geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
        coord_sf(xlim = xlims, ylim = ylims) +
        labs(fill = glue::glue('LGP\n(days)  '), 
             title = glue::glue('gSeason = {i}; Historic'),x = 'Longitude', y = 'Latitude') +
        scale_fill_gradientn(colours = blues9, limits = c(limits_two$LGP_min[i], limits_two$LGP_max[i]), 
                             guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
        # ggspatial::annotation_scale(location = "br", width_hint = 0.5) +
        # ggspatial::annotation_north_arrow(location = "br", which_north = "true",
        # pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), # 0.2 # 0.3
        # style = north_arrow_fancy_orienteering) +
        scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
        scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
        theme_bw() +
        theme(legend.position = 'bottom')
      
      ggsave(glue::glue('{path}{Country}/graphs/{county}/maps/LGP_past_{i}.png') , width = 8, height = 4, dpi = 300)
      
      
      LGP_f <- ggplot(filter(two_index %>% mutate(gSeason = glue::glue('gSeason = {gSeason}')),
                             time == 'future', gSeason == glue::glue('gSeason = {i}'))) +
        geom_tile(aes(x = x, y = y, fill = LGP))  +
        geom_sf(data = pais, fill = NA, color = gray(.8)) +
        geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
        coord_sf(xlim = xlims, ylim = ylims) +
        labs(fill = glue::glue('LGP\n(days)  '), 
             title = glue::glue('gSeason = {i}; Future'), x = 'Longitude', y = 'Latitude') +
        scale_fill_gradientn(colours = blues9, limits = c(limits_two$LGP_min[i], limits_two$LGP_max[i]), 
                             guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
        # ggspatial::annotation_scale(location = "br", width_hint = 0.5) +
        # ggspatial::annotation_north_arrow(location = "br", which_north = "true",
        # pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), # 0.2 # 0.3
        # style = north_arrow_fancy_orienteering) +
        scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
        scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
        theme_bw() +
        theme(legend.position = 'bottom')
      
      ggsave(glue::glue('{path}{Country}/graphs/{county}/maps/LGP_future_{i}.png') , width = 8, height = 4, dpi = 300)
      
      
      
      e <- ggplot(filter(LGP_dif, gSeason == glue::glue('gSeason = {i}')) ) +
        geom_tile(aes(x = x, y = y, fill = LGP)) +
        geom_sf(data = pais, fill = NA, color = gray(.8)) +
        geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
        coord_sf(xlim = xlims, ylim = ylims) +
        labs(fill = glue::glue('LGP\n(days)  '), 
             title = glue::glue('gSeason = {i}; Change'),x = 'Longitude', y = 'Latitude') +
        scale_fill_gradient2(low = '#A50026', mid = 'white', high = '#000099', 
                             guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
        scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
        scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
        theme_bw() +
        theme(legend.position = 'bottom')
      
      ggsave(glue::glue('{path}{Country}/graphs/{county}/maps/Dif_LGP_{i}.png') , width = 8, height = 4, dpi = 300)
      
      png(filename = glue::glue('{path}{Country}/graphs/{county}/maps/all_LGP_{i}.png') , width=12.5,height=4.5,units="in", res = 300)
      print(gridExtra::grid.arrange(LGP_p, LGP_f, e, ncol=3,  
                                    top = glue::glue('{Country}, {county}',
                                                     bottom =   "Data source: Alliance Bioversity-CIAT")))
      dev.off()
      
    }
    
    
  }
  
  
  ##  =-----------------------------------------
  Clim_graph <- function(historic){
    Country <- unique(historic$Country)
    
    shp_sf <- shp  %>% sf::st_as_sf()
    country <- country1 %>% sf::st_as_sf()
    xlims <- sf::st_bbox(shp_sf)[c(1, 3)]
    ylims <- sf::st_bbox(shp_sf)[c(2, 4)]
    
    # =--------------------------------------------------------
    historic <- historic %>% replace_na(list(z = mean)) %>%
      dplyr::group_by( id, x, y, ISO3, Country) %>% 
      dplyr::summarise_all(~round(. , 1))
    # =--------------------------------------------------------
    
    # prueba <- historic %>% ungroup %>% dplyr::select(prec) %>% 
    #   arrange(prec) %>% slice(1, n()) %>% .$prec
    # 
    # range <- (prueba[2]-prueba[1])/5
    # 
    # prec_p <- tibble(tr =  c('a', 'b', 'c', 'd'), 
    # limits = c(prueba[1]+range, round(prueba[1]+(2*range), 1), round(prueba[1]+(3*range), 1), prueba[2]))
    # labels_l <- c(glue::glue('< {prec_p[2,2]}'), glue::glue('{prec_p[2,2]}-{prec_p[3,2]}'), glue::glue('> {prec_p[3,2]}'))
    
    # historic_p <- historic %>% dplyr::select(id, x, y, ISO3, Country, prec) %>%
    #   mutate(val = case_when( prec < prec_p[2,2] ~ 1, 
    #                                        prec >= prec_p[2,2] & prec < prec_p[3,2]  ~ 2,
    #                                        prec >= prec_p[3,2]  ~ 3,
    #                                        TRUE ~ 4) %>% as.factor())
    
    
    # =--------------------------------------------------------
    
    prec <- ggplot() + geom_tile(data = historic, aes(x = x, y = y, fill = prec )) +
      geom_sf(data = country, fill = NA, color = gray(.8)) +
      geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
      coord_sf(xlim = xlims, ylim = ylims) +
      labs(fill = glue::glue('(mm)  '), title = 'Historical Annual Mean Precipitation (mm/year)',x = 'Longitude', y = 'Latitude') +
      # scale_fill_brewer(palette="Greens",  labels = labels_l, 
      # guide = guide_legend(barwidth = 12, label.theme = element_text(angle = 0))) +
      scale_fill_gradientn(colours = blues9, 
                           guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0))) +
      scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
      scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
      theme_bw() + theme(legend.position = 'bottom')
    
    
    ggsave(glue::glue('{path}{Country}/graphs/{county}/maps/H_prec.png') , width = 8, height = 5.5, dpi = 300)
    
    ##########################################################
    
    # prueba_t <- historic %>% ungroup %>% dplyr::select(tmean) %>% 
    #   arrange(tmean) %>% slice(1, n()) %>% .$tmean
    # 
    # range_t <- (prueba_t[2]-prueba_t[1])/4
    # 
    # temp_p <- tibble(tr =  c('a', 'b', 'c', 'd'), 
    #                  limits = c(prueba_t[1]+range_t, round(prueba_t[1]+(2*range_t), 1), round(prueba_t[1]+(3*range_t), 1), prueba_t[2]))
    # labels_t <- c(glue::glue('< {temp_p[2,2]}'), glue::glue('{temp_p[2,2]}-{temp_p[3,2]}'), glue::glue('> {temp_p[3,2]}'))
    
    # historic_t <- historic %>% dplyr::select(id, x, y, ISO3, Country, tmean) %>%
    #   mutate(val = case_when( tmean < temp_p[2,2] ~ 1, 
    #                           tmean >= temp_p[2,2] & tmean < temp_p[3,2]  ~ 2,
    #                           tmean >= temp_p[3,2]  ~ 3,
    #                           TRUE ~ 4) %>% as.factor())
    
    tmn <- ggplot() + geom_tile(data = historic, aes(x = x, y = y, fill = tmean))+
      geom_sf(data = country, fill = NA, color = gray(.8)) +
      geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
      coord_sf(xlim = xlims, ylim = ylims) +
      labs(fill = expression('('*~degree*C*')'), title = expression('Historical Annual Mean Temperature ('*~degree*C*')'),x = 'Longitude', y = 'Latitude') +
      # scale_fill_brewer(palette="YlOrRd",  labels = labels_t, 
      # guide = guide_legend(barwidth = 12, label.theme = element_text(angle = 0))) +
      scale_fill_gradient(low = "yellow", high = "red",
                          guide = guide_colourbar(barwidth = 12, label.theme = element_text(angle = 0)))+
      scale_y_continuous(breaks = round(ylims, 2), n.breaks = 3) +
      scale_x_continuous(breaks = round(xlims, 2), n.breaks = 3) +
      theme_bw() + theme(legend.position = 'bottom')
    
    ggsave(glue::glue('{path}{Country}/graphs/{county}/maps/H_tmn.png') , width = 8, height = 5.5, dpi = 300)
    
    
    png(filename = glue::glue('{path}{Country}/graphs/{county}/maps/A_Multi_Anual.png'), width=12.5,height=4.5,units="in", res = 300)
    print(gridExtra::grid.arrange(prec, tmn, ncol=2,
                                  top = glue::glue('{Country}, {county}',
                                                   bottom =   "Data source: Alliance Bioversity-CIAT")))
    dev.off()
  }
  
  
  
  
  # =---------------------------------------------------------------
  # Cowsay para radiación. 
  cowsay::say('all maps', by = 'smallcat')
  # =----------------------------------------------------------------
  # Aqui se van a poner los mapas de radiación solar... 
  
  
  # Observed data for each country... 
  tictoc::tic()
  all_climate <-  fst::fst(glue::glue('//dapadfs/workspace_cluster_8/climateriskprofiles/data/observational_data/{co}/{county}.fst')) %>% 
    tibble::as_tibble() %>% 
    dplyr::select(-id1)  %>% 
    tidyr::nest(-id, -x, -y, -ISO3) %>% 
    dplyr::mutate(Country = country) %>% 
    dplyr::rename(climate = 'data') %>% 
    dplyr::select(id, x, y, ISO3, Country, climate)
  tictoc::toc()
  
  
  # historic <- all_climate %>%
  #   dplyr::mutate(summary =  purrr::map(.x = climate, .f = function(z){
  #     z <- z %>% dplyr::mutate(year = lubridate::year(Date), month = lubridate::month(Date),
  #                         tmean = (tmax + tmin)/2 ) %>%  
  #       dplyr::group_by(year, month) %>% 
  #       dplyr::summarise(prec = sum(prec), tmean = mean(tmean)) %>% 
  #       dplyr::ungroup() %>% dplyr::group_by(month) %>% dplyr::select(-year) %>% 
  #       dplyr::summarise_all(.funs = function(x){round( mean(x, na.rm = TRUE), 1)}) %>% 
  #       dplyr::ungroup()})) %>% dplyr::select(-climate) %>% tidyr::unnest() %>% 
  #   dplyr::group_by( id, x, y, ISO3, Country) %>% 
  #   dplyr::summarise_all(~mean(.))
  
  historic <-  all_climate %>%
    dplyr::mutate(summary =  purrr::map(.x = climate, .f = function(z){
      z <- z %>% dplyr::mutate(year = lubridate::year(Date), month = lubridate::month(Date) ,
                               tmean = (tmax + tmin)/2 )  %>%
        dplyr::group_by(year) %>%
        dplyr::summarise(prec = sum(prec), tmean = mean(tmean)) %>%
        dplyr::ungroup() %>% dplyr::select(-year) %>%
        dplyr::summarise_all(.funs = function(x){round( mean(x, na.rm = TRUE), 1)}) %>%
        dplyr::ungroup()})) %>% dplyr::select(-climate) %>% tidyr::unnest()  %>%
    dplyr::group_by( id, x, y, ISO3, Country) %>%
    dplyr::summarise_all(~mean(.))
  
  historic <- historic %>% filter(id %in% crd$id)
  
  path <- '//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles/results/'
  # path <- "D:/OneDrive - CGIAR/Desktop/P_indices_H/Profiles_AA/probando/"
  dir.create(glue::glue('{path}{country}/graphs/{tolower(county)}/maps'),recursive = TRUE)
  # 
  # # graph de clima. 
  Clim_graph(historic)
  
  
  if(Big == 'N'){
    past_c <- fst::fst(glue::glue('{path}{country}/past/{county}_1985_2015_corrected.fst')) %>%
      as_tibble() %>%
      dplyr::mutate(time = 'past')
    
    futDir  <- paste0('//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles/results/',country,'/future')
    fut_fls <- list.files(futDir, pattern = paste0('^',county,'_[0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9]_corrected.fst'), recursive = T)
    fut_fls <- paste0(futDir,'/',fut_fls)
    
    future_c  <- fut_fls %>%
      purrr::map(.f = function(x){df <- fst(x) %>% as_tibble() %>% dplyr::mutate(time = 'future'); return(df)}) %>%
      dplyr::bind_rows()
    
    
    # Climate Index
    data_all <- dplyr::bind_rows(past_c, future_c) %>% dplyr::select(-x, -y) %>%
      dplyr::inner_join(., all_climate) %>% filter(id %in% crd$id) %>% 
      dplyr::mutate(county = county)  %>%
      dplyr::mutate_at(vars(CDD:ndws), ~ifelse(is.na(.), mean(., na.rm = TRUE), .) %>% round()) %>%
      dplyr::mutate(ndws = ifelse(ndws == 0, median(ndws), ndws)) %>%
      tidyr::separate(season, c('ok','semester'), sep = 's') %>%
      dplyr::select(id, ISO3, county, Country, x, y,  time ,semester, year, CDD, P5D, P95, NT35, ndws) %>%
      dplyr::group_split(semester)
    
    data_all %>% purrr::walk(.f = do_climate_maps)
    
    
    # Srad index....
    index_complete <- dplyr::bind_rows(past_c, future_c)  %>% dplyr::select(-x, -y) %>%
      dplyr::inner_join(., all_climate) %>% filter(id %in% crd$id) %>% 
      dplyr::mutate(county = county) %>%
      dplyr::mutate_at(vars(CDD:LGP), ~ifelse(is.na(.), mean(., na.rm = TRUE), .) %>% round()) %>%
      dplyr::mutate(ndws = ifelse(ndws == 0, median(ndws), ndws)) %>%
      dplyr::select(id, ISO3, county, Country, x, y,  time, year, gSeason, SLGP,LGP)  %>%
      dplyr::mutate(Big = 'N')
    
    do_srad_maps(data_split = index_complete)
    
  }else if(Big == 'B'){
    # Si... se tienen los idw.
    past_c <- fst::fst(glue::glue('{path}{country}/past/{county}_1985_2015_idw.fst')) %>%
      as_tibble() %>%
      mutate(time = 'past')
    
    futDir  <- paste0('//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles/results/',country,'/future')
    fut_fls <- list.files(futDir, pattern = paste0('^',county,'_[0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9]_idw.fst'), recursive = T)
    fut_fls <- paste0(futDir,'/',fut_fls)
    
    future_c  <- fut_fls %>%
      purrr::map(.f = function(x){df <- fst(x) %>% as_tibble() %>% mutate(time = 'future'); return(df)}) %>%
      dplyr::bind_rows()
    
    
    # Climate Index
    
    data_all <- dplyr::bind_rows(past_c, future_c) %>% dplyr::select(-x, -y) %>%
      dplyr::inner_join(., all_climate) %>% filter(id %in% crd$id) %>% 
      dplyr::mutate(county = county) %>%
      separate(season, c('ok','semester'), sep = 's') %>%
      dplyr::select(id, ISO3, county, Country, x, y,  time ,semester, CDD, P5D, P95, NT35, ndws) %>%
      unique() %>%
      group_split(semester)
    
    data_all %>% purrr::walk(.f = do_climate_maps)
    
    # Srad index....
    index_complete <- bind_rows(past_c, future_c) %>% dplyr::select(-x, -y) %>%
      inner_join(., all_climate) %>% filter(id %in% crd$id) %>% 
      mutate(county = county) %>%
      dplyr::select(id, ISO3, county, Country, x, y, time, gSeason, SLGP,LGP) %>%
      mutate(Big = Big)
    
    # index_complete %>% purrr::walk(.f = do_srad_maps)
    do_srad_maps(data_split = index_complete)
    
  }
  
}
