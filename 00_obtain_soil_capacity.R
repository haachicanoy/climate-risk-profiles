# Soil capacity calculation: GIZ climate-hazards
# By: H. Achicanoy
# CIAT, 2020

pacman::p_load(raster, tidyverse, fst, GSIF, vroom)

root_depth <- 60 # cm

coords <- vroom::vroom('//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles/data/id_all_country.csv')
coords <- coords %>% dplyr::select(id,x,y)

soils_root <- '//catalogue/BaseLineData_cluster04/GLOBAL/Biofisico/SoilGrids250m'
orc <- raster::stack(list.files(paste0(soils_root,'/Chemical soil properties/Soil organic carbon content'), pattern = '.tif$', full.names = T) %>% sort())
cec <- raster::stack(list.files(paste0(soils_root,'/Chemical soil properties/Cation exchange capacity (CEC)'), pattern = '.tif$', full.names = T) %>% sort())
phx <- raster::stack(list.files(paste0(soils_root,'/Chemical soil properties/Soil ph in H2O'), pattern = '.tif$', full.names = T) %>% sort())
snd <- raster::stack(list.files(paste0(soils_root,'/Physical soil properties/Sand content'), pattern = '.tif$', full.names = T) %>% sort())
slt <- raster::stack(list.files(paste0(soils_root,'/Physical soil properties/Silt content'), pattern = '.tif$', full.names = T) %>% sort())
cly <- raster::stack(list.files(paste0(soils_root,'/Physical soil properties/Clay content (0-2 micro meter) mass fraction'), pattern = '.tif$', full.names = T) %>% sort())
bld <- raster::stack(list.files(paste0(soils_root,'/Physical soil properties/Bulk density (fine earth)'), pattern = '.tif$', full.names = T) %>% sort())

soil <- raster::stack(orc,cec,phx,snd,slt,cly,bld)
soil_data <- cbind(coords, raster::extract(soil, coords[,c('x','y')]))

soil_data2 <- soil_data %>%
  tidyr::gather(key = 'var', value = 'val', -(1:3)) %>%
  tidyr::separate(col = 'var', sep = '_M_', into = c('var','depth')) %>%
  tidyr::spread(key = 'var', value = 'val') %>%
  dplyr::arrange(id)
soil_data2$depth <- gsub('_250m_ll','',soil_data2$depth)

fst::write_fst(soil_data2, '//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles/data/soil_data.fst')

soil_data2$d <- lapply(1:nrow(soil_data2), function(i){
  y <- GSIF::AWCPTF(SNDPPT = soil_data2$SNDPPT[i],
                    SLTPPT = soil_data2$SLTPPT[i],
                    CLYPPT = soil_data2$CLYPPT[i],
                    ORCDRC = soil_data2$ORCDRC[i],
                    BLD = soil_data2$BLDFIE[i],
                    CEC = soil_data2$CECSOL[i],
                    PHIHOX = soil_data2$PHIHOX[i]/10,
                    h1=-10, h2=-20, h3=-33)
  y <- y$AWCh3
  return(y)
}) %>% base::unlist()

soil_data3 <- soil_data2 %>%
  dplyr::select(id,x,y,depth,d) %>%
  tidyr::spread(key='depth',value='d')

names(soil_data3)[4:ncol(soil_data3)] <- paste0('d.',c(0, 5, 15, 30, 60, 100, 200))
soil_data3$rdepth <- root_depth
soil_data3 <- soil_data3 %>%
  dplyr::select('id','x','y','rdepth',paste0('d.',c(0, 5, 15, 30, 60, 100, 200)))

soil_data4 <- soil_data3 %>%
  dplyr::mutate(aws.25 = (d.0+d.5)/2 * 100) %>%
  dplyr::mutate(aws.100 = (d.5+d.15)/2 * 100) %>%
  dplyr::mutate(aws.225 = (d.15+d.30)/2 * 100) %>%
  dplyr::mutate(aws.450 = (d.30+d.60)/2 * 100) %>%
  dplyr::mutate(aws.800 = (d.60+d.100)/2 * 100) %>%
  dplyr::mutate(aws.1500 = (d.100+d.200)/2 * 100) %>%
  dplyr::select('id','x','y','rdepth',paste0('aws.',c(25,100,225,450,800,1500)))

soilcap_calc_mod <- function(x, minval, maxval) {
  if(!is.na(x[4])){
    rdepth <- max(c(x[4],minval)) #cross check
    rdepth <- min(c(rdepth,maxval)) #cross-check
    wc_df <- data.frame(depth=c(2.5,10,22.5,45,80,150),wc=(x[5:10])*.01)
    if (!rdepth %in% wc_df$depth) {
      wc_df1 <- wc_df[which(wc_df$depth < rdepth),]
      wc_df2 <- wc_df[which(wc_df$depth > rdepth),]
      y1 <- wc_df1$wc[nrow(wc_df1)]; y2 <- wc_df2$wc[1]
      x1 <- wc_df1$depth[nrow(wc_df1)]; x2 <- wc_df2$depth[1]
      ya <- (rdepth-x1) / (x2-x1) * (y2-y1) + y1
      wc_df <- rbind(wc_df1,data.frame(depth=rdepth,wc=ya),wc_df2)
    }
    wc_df <- wc_df[which(wc_df$depth <= rdepth),]
    wc_df$soilthick <- wc_df$depth - c(0,wc_df$depth[1:(nrow(wc_df)-1)])
    wc_df$soilcap <- wc_df$soilthick * wc_df$wc
    soilcp <- sum(wc_df$soilcap) * 10 #in mm
    return(soilcp)
  } else {
    soilcp <- NA
    return(soilcp)
  }
}

soil_data4$soilcp <- apply(soil_data4, 1, FUN=soilcap_calc_mod, minval=45, maxval=100)

fst::write_fst(soil_data4, '//dapadfs.cgiarad.org/workspace_cluster_8/climateriskprofiles/data/soilcp_data.fst')
