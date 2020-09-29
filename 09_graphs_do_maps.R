rm(list = ls())
gc(reset = TRUE)

# =--------------------
# Packages 
library(tidyverse)
library(raster)
library(ncdf4)
library(sf)
library(future)
library(furrr)
library(lubridate)
library(glue)
library(cowsay)
library(fst)
library(ggspatial)


options(warn = -1, scipen = 999)

suppressMessages(library(pacman))
suppressMessages(pacman::p_load(vroom, sp, compiler))
# =--------------------

# =----------------------------------
# Identificacion de pixel para ETH
# =----------------------------------
country <- 'Ethiopia'
county  <- 'Arsi'
country1 <- 'ethiopia'


# =--------------------

# Ruta Principal para guardados: 
root <- '//dapadfs/workspace_cluster_8/climateriskprofiles/data/'

# Load county shapefile
country1 <- raster::shapefile(paste0(root,'/shps/Ethiopia/ETH_adm2.shp'))
shp <- raster::shapefile(paste0(root,'/shps/Ethiopia/ETH_adm2.shp'))
shp <- shp[shp@data$NAME_2 == county,]
plot(shp)

# Load id coords
crd <- vroom(paste0(root,'/id_country.csv'), delim = ',')
crd <- crd %>%
  dplyr::filter(Country == country)
pnt <- crd %>% dplyr::select(x,y) %>% sp::SpatialPoints(coords = .)
crs(pnt) <- crs(shp)
# Filter coordinates that are present in the county
pnt <- sp::over(pnt, shp) %>% data.frame %>% dplyr::select(ID_0:TYPE_2) %>% complete.cases() %>% which()
crd <- crd[pnt,]
crd <<- crd  %>% mutate(county = county)



#=---------------------

# Observed data for each country... 
tictoc::tic()
observacional_data <- read_rds('//dapadfs/workspace_cluster_8/climateriskprofiles/data/observational_data/ethiopia/arsi.RDS')
tictoc::toc()

# =--------------------


path <- '//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles/results/'

past <- fst::fst(glue::glue('{path}{country}/past/{county}_1985_2015.fst')) %>% 
  as_tibble() %>% 
  mutate(time = 'past')
fut_fls <- list.files(path = glue::glue('{path}{country}/future'), pattern = county, recursive = TRUE, full.names = TRUE)

future  <- fut_fls %>%
  purrr::map(.f = function(x){df <- fst(x) %>% as_tibble() %>% mutate(time = 'future'); return(df)}) %>%
  bind_rows()


data_all <- bind_rows(past, future) %>% 
  inner_join(., observacional_data) %>% 
  mutate(county = county) %>% 
  dplyr::select(id, ISO3, county, Country, x, y,  time ,semester, CDD, P5D, P95, NT35)
  

data_all


shp_sf <- shp  %>% sf::st_as_sf()
country <- country1 %>% sf::st_as_sf()
xlims <- sf::st_bbox(shp_sf)[c(1, 3)]
ylims <- sf::st_bbox(shp_sf)[c(2, 4)]




data_all %>% nest(-ISO3, -county, -Country, - semester, -time )





b <- ggplot() +
  geom_sf(data = shp_sf, fill = 'red', color = gray(.1)) +
  geom_sf(data = country, fill = NA, color = gray(.5)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 7, face = "bold")) 


index <- 'CDD'

a <- ggplot() +
  geom_tile(data = median_data, aes(x = x, y = y, fill = CDD)) +
  geom_sf(data = country, fill = NA, color = gray(.8)) +
  geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
  coord_sf(xlim = xlims, ylim = ylims) +
  labs(fill = glue::glue('{index}\n(days)'), x = 'Longitude', y = 'Latitude') + 
  scale_fill_viridis_c() + 
  ggspatial::annotation_scale(location = "br", width_hint = 0.5) +
  ggspatial::annotation_north_arrow(location = "br", which_north = "true", 
                                    pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), # 0.2 # 0.3
                                    style = north_arrow_fancy_orienteering) +
  theme_bw() 


Mapa_cons <-  cowplot::ggdraw() +
  cowplot::draw_plot(a, x = -0.08) +
  cowplot::draw_plot(b, x = 0.68, y = 0.68, width = 0.4, height = 0.3)

ggsave(glue::glue('{output_P}/graphs/{index}_1985_2015.png') , width = 8, height = 5.5)


# ggsave(glue::glue('{output_P}{index}_1985_2015.png') , width = 7.2, height = 5.5)
index <- 'P5D'

c <- ggplot() +
  geom_tile(data = median_data, aes(x = x, y = y, fill = P5D)) +
  geom_sf(data = country, fill = NA, color = gray(.8)) +
  geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
  coord_sf(xlim = xlims, ylim = ylims) +
  labs(fill = glue::glue('{index}\n(mm)'), x = 'Longitude', y = 'Latitude') + 
  scale_fill_gradientn(colours = blues9) + 
  ggspatial::annotation_scale(location = "br", width_hint = 0.5) +
  ggspatial::annotation_north_arrow(location = "br", which_north = "true", 
                                    pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), # 0.2 # 0.3
                                    style = north_arrow_fancy_orienteering) +
  theme_bw() 


Mapa_cons <-  cowplot::ggdraw() +
  cowplot::draw_plot(c, x = -0.08) +
  cowplot::draw_plot(b, x = 0.68, y = 0.68, width = 0.4, height = 0.3)

Mapa_cons
ggsave(glue::glue('{output_P}/graphs/{index}_1985_2015.png') , width = 8, height = 5.5)



# =------------
index <- 'P95'

d <- ggplot() +
  geom_tile(data = median_data, aes(x = x, y = y, fill = P95)) +
  geom_sf(data = country, fill = NA, color = gray(.8)) +
  geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
  coord_sf(xlim = xlims, ylim = ylims) +
  labs(fill = glue::glue('{index}\n(mm)'), x = 'Longitude', y = 'Latitude') + 
  scale_fill_gradientn(colours = blues9) +
  ggspatial::annotation_scale(location = "br", width_hint = 0.5) +
  ggspatial::annotation_north_arrow(location = "br", which_north = "true", 
                                    pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), # 0.2 # 0.3
                                    style = north_arrow_fancy_orienteering) +
  theme_bw() 


Mapa_cons <-  cowplot::ggdraw() +
  cowplot::draw_plot(d, x = -0.08) +
  cowplot::draw_plot(b, x = 0.68, y = 0.68, width = 0.4, height = 0.3)


Mapa_cons
ggsave(glue::glue('{output_P}/graphs/{index}_1985_2015.png') , width = 8, height = 5.5)




# =------------
index <- 'NT35'

e <- ggplot() +
  geom_tile(data = median_data, aes(x = x, y = y, fill = as.factor(NT35))) +
  geom_sf(data = country, fill = NA, color = gray(.8)) +
  geom_sf(data = shp_sf, fill = NA, color = gray(.1)) +
  coord_sf(xlim = xlims, ylim = ylims) +
  labs(fill = glue::glue('{index}\n(days)'), x = 'Longitude', y = 'Latitude') + 
  scale_fill_viridis_d() + 
  ggspatial::annotation_scale(location = "br", width_hint = 0.5) +
  ggspatial::annotation_north_arrow(location = "br", which_north = "true", 
                                    pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), # 0.2 # 0.3
                                    style = north_arrow_fancy_orienteering) +
  theme_bw() 


Mapa_cons <-  cowplot::ggdraw() +
  cowplot::draw_plot(e, x = -0.08) +
  cowplot::draw_plot(b, x = 0.68, y = 0.68, width = 0.4, height = 0.3)

Mapa_cons
ggsave(glue::glue('{output_P}/graphs/{index}_1985_2015.png') , width = 8, height = 5.5)



png(filename = glue::glue('{output_P}/graphs/climate_I_1985_2015.png') , width = 1280, height = 720)
print(gridExtra::grid.arrange(a, c, d, e, layout_matrix = rbind(c(1, 2), c(3, 4))))
dev.off()

