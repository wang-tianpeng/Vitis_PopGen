rm(list=ls()); graphics.off()


################ generar objetos de entrada #####################
N <- 10
s <- c(9)
### crear datos de presencias
if(FALSE){
  library(readr)
  library(terra)
  library(tidyverse)
  library(CoordinateCleaner)
  
  
  # berlandieri:
  fil <- "inputs/0037997-251120083545085.csv"
  
  
  lines <- readLines(fil, warn = FALSE)
  length(lines)
  sps <- read_tsv(fil,  guess_max = 10000, show_col_types = FALSE, progress = FALSE)
  sps <- sps[!is.na(sps$species),]
  sps <- sps[!is.na(sps$decimalLatitude),]
  
  sps <- sps[,c("species","countryCode","locality","occurrenceStatus","decimalLatitude","decimalLongitude","coordinateUncertaintyInMeters","coordinatePrecision","elevation","year","identifiedBy")]
  
  #   "champinii",
  #   "monticola",
  #   "arizonica",
  #   "girdiana",
  #   "californica"),sep = " ")
  
  lsp <- paste("Vitis",
               c("riparia","acerifolia","shuttleworthii","aestivalis",
                 "cinerea","labrusca","mustangensis",
                 "vulpina","rotundifolia"),
               sep=" ")
  
  
  limpiar <- data.frame(lon=NA,lat=NA,species=NA)
  ras <- rast("inputs/pres.tif")
  pdf("output/mapas_limpios_2.pdf")
  for(i in 1:length(lsp)){
  
    plot(ras,1,main=lsp[i])
    sp <- sps[sps$species==lsp[i],]
    points(sp[,c("decimalLongitude","decimalLatitude")])
    sp <- as.data.frame(sp)
    sp <- sp[,c("species","decimalLongitude","decimalLatitude")]
    
    r <- rast(resolution = 0.1, crs = "EPSG:4326") # ~10 km
    p <- vect(sp, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")
    
    celda <- cellFromXY(r, crds(p))
    sp$celda <- celda
    
    sp <- sp %>%
      group_by(celda) %>%
      slice_sample(n = 1) %>%
      ungroup() %>% 
      as.data.frame()
    sp <- sp[,2:3]
  
    names(sp) <- c("lon","lat")
    points(sp[,c("lon","lat")],col="red") 
    
    sp$species <- lsp[i]
    sp <- CoordinateCleaner::cc_cap(x = sp,lon = "lon",lat = "lat")
    sp <- CoordinateCleaner::cc_equ(x = sp,lon = "lon",lat = "lat")
    sp <- CoordinateCleaner::cc_cen(x = sp,lon = "lon",lat = "lat")
    sp <- CoordinateCleaner::cc_gbif(x = sp,lon = "lon",lat = "lat")
    sp <- CoordinateCleaner::cc_inst(x = sp,lon = "lon",lat = "lat")
    sp <- CoordinateCleaner::cc_sea(sp,lon = "lon",lat = "lat")
    sp <- CoordinateCleaner::cc_outl(x = sp,lon = "lon",lat = "lat",
                                           species = "species",method = "quantile")
    points(sp[,c("lon","lat")],col="green")
    limpiar <- rbind(limpiar,sp)
    
  }
  limpiar <- limpiar[-1,]
  
  print(summary(factor(limpiar$species)))
  print(summary(factor(sps$species)))
  sps <- limpiar
  
  pres <- rast("inputs/pres.tif")
  sps$temp <- terra::extract(pres,sps[,c("lon","lat")])[,2]
  sps <- sps[!is.na(sps$temp),] 
  sps <- sps[,-ncol(sps)]
  
  save(sps,file = "inputs/especies_input2.R")
  dev.off()
}

if(FALSE){
  source("scripts/functions.R")
  load("inputs/especies_input2.R")
  
  limpiar <- data.frame(lon=NA,lat=NA,species=NA)
  lsp <- unique(sps$species)
  for(i in 1:length(lsp)){
    sp <- sps[which(sps$species==lsp[i]),]
    sp <- vect(sp)
    extens <- ext(sp)
    ecorregions <- terra::vect("~/Desktop/capas_climatica/ecoregions/wwf_terr_ecos.shp") 
    ecorregions <- terra::crop(ecorregions,extens) %>% magrittr::extract(.,sp)
    bioclima <- terra::rast("inputs/pres.tif")
    bioclima <-terra::crop(bioclima,ecorregions) %>% mask(.,ecorregions)
    covarData <- terra::extract(bioclima, sp,ID=FALSE) 
    select_var <- fuzzySim::corSelect(
      data = covarData, var.cols = names(covarData),cor.thresh = 0.8,
      use = "pairwise.complete.obs",select = "VIF")
    select_var <- select_var$selected.vars
    bioclima <- bioclima[[select_var]]
    
  
    
    records <- geom(sp)[,c("x","y")]
    records <- data.frame(species=sp[i],records)
    names(records) <- c("species","longitude","latitude") 
    
    variables <- bioclima
    sp_out <- lsp[i]
    sp_out <- sub(" ","_",sp_out)
    out_dir <- file.path(paste("output/grinell/",lsp[i],sep = ""))
    m <- M_simulationR(data = records, current_variables = variables,starting_proportion = 0.5,sampling_rule = "suitability",
                       max_dispersers = 5, replicates = 3, dispersal_events = 25,
                       output_directory = out_dir,overwrite = T)
    
    plot(bioclima[[1]],col=cols,main=names(bioclima)[1])
    points(records[,2:3],pch=19,cex=0.8)
    points(m$Simulation_occurrences[,2:3],col="green",pch=19,cex=0.5)
    sp <- m$Simulation_occurrences[,c("longitude","latitude")]
    sp$species <- lsp[i]
    names(sp) <- c("lon","lat","species")
    limpiar <- rbind(limpiar,sp)
  }
  limpiar <- limpiar[-1,]
  sps <- limpiar
  
  
  save(sps,file = "inputs/especies_input_limpias2.R")

}




if(FALSE){
  library(terra)
  library(tidyverse)
  bio <- "bio_1"
  pres <- rast("inputs/pres.tif")
  mri_245 <- rast("inputs/mri_245_70.tif")
  mri_585 <- rast("inputs/mri_585_70.tif")
  mpi_245 <- rast("inputs/mpi_245_70.tif")
  mpi_585 <- rast("inputs/mpi_585_70.tif")
  load("inputs/especies_input_limpias.R")
  
  
  sps$pres <- terra::extract(pres,sps[,c("lon","lat")],ID=FALSE) %>% .[,bio]
  sps$mri_245 <- terra::extract(mri_245,sps[,c("lon","lat")],ID=FALSE) %>% .[,bio]
  sps$mri_585 <- terra::extract(mri_585,sps[,c("lon","lat")],ID=FALSE) %>% .[,bio]
  sps$mpi_245 <- terra::extract(mpi_245,sps[,c("lon","lat")],ID=FALSE) %>% .[,bio]
  sps$mpi_585 <- terra::extract(mpi_585,sps[,c("lon","lat")],ID=FALSE) %>% .[,bio]
  
  datos_long <- sps[,c("species","pres","mri_245","mri_585","mpi_245","mpi_585")] %>% 
    pivot_longer(cols=pres:mpi_585,
                 names_to = "escenario",
                 values_to = "clima")
  datos_long$escenario <- factor(datos_long$escenario,
                                 levels = c("pres","mri_245","mri_585","mpi_245","mpi_585"))
  ggplot(datos_long, aes(x = species, y = clima, fill = escenario)) +
    geom_boxplot(position = position_dodge(width = 0.8)) +
    scale_fill_manual(values = c(
      "pres" = "grey",
      "mri_245" = "lightblue",
      "mri_585" = "salmon",
      "mpi_245" = "cyan",
      "mpi_585" = "orange"
    ))+
    labs(x = "Especie", y = bio, fill = "Escenario") +
    theme_classic(base_rect_size = NULL)
  
  

}



############ correr modelos de nicho #########
source("scripts/functions.R")

load("inputs/especies_input_limpias2.R")
presente <- rast("inputs/pres.tif")
# plot(presente$bio_1,col=cols)
# points(sps[,c("lon","lat")],col=factor(sps$species),pch="+")
sps$species <- sub(" ",".",sps$species)
lista_sp <- sps$species %>% unique()
 i <- 1
 for(species in lista_sp[s]){
   pdf(paste("output/salidas_pdf/",species,".pdf",sep = ""))
   print(paste("runing", species, "missing", length(lista_sp)-i))
   
  input <- get_input(species,sps,extens=ext(presente),bioclima=presente)
  dir.create(path = paste("output/",species,sep = ""))
  save(input,file = paste("output/",species,"/input.R",sep = ""))
  if(dir.exists(input$myRespName)){
    unlink(input$myRespName,recursive = T)
  }
  
  # model_sdm: genera el modelo de nicho
  myBiomodEM <- model_sdm(formato_input = input,perc = 0.7,modelo = c("GAM","MAXENT","RF"),N_run = N)

    if(is.character(myBiomodEM)){
    if(myBiomodEM=="not_run"){
      next
    }
    
  }
   
  save(myBiomodEM,file = paste("output/",species,"/model.R",sep = ""))
  
  
  #### proyectar modelo
  myExpl_pres <- input$myExpl
  
  pres <- project_model(myBiomodEM = myBiomodEM,myExpl = myExpl_pres,name = "present",formato_input = input)
  plot(pres,col=cols,main="present")
  
  past_bioclim <- rast("inputs/lgm_tr10.tif")
  past_bioclim <- crop(past_bioclim,ext(input$ecorregions))
  myExpl_past <-past_bioclim[[myBiomodEM@expl.var.names]] %>% mask(.,input$ecorregions)
  lgm_model <- project_model(myBiomodEM = myBiomodEM,myExpl = myExpl_past,name = "lgm",formato_input = input)
  plot(lgm_model,col=cols,main="lgm")
  
  past_bioclim <- rast("inputs/ig_tr10.tif")
  past_bioclim <- crop(past_bioclim,ext(input$ecorregions))
  myExpl_past <-past_bioclim[[myBiomodEM@expl.var.names]] %>% mask(.,input$ecorregions)
  ig_model <- project_model(myBiomodEM = myBiomodEM,myExpl = myExpl_past,name = "ig",formato_input = input)
  plot(ig_model,col=cols,main="ig")
  
  
  
  
  writeRaster(pres,filename =paste("output/",species,"/pres.tif",sep = ""))
  writeRaster(lgm_model,filename =paste("output/",species,"/lgm.tif",sep = ""))
  writeRaster(ig_model,filename =paste("output/",species,"/ig.tif",sep = ""))
  i <- i +1
  dev.off()
}



 
 
 
 
 
 
 stop()
########## analizar datos ########## 
rm(list=ls());graphics.off()
apis_pres <- rast("output/Apis.mellifera/pres.tif") 
pasi_pres<- rast("output/Passiflora.foetida/pres.tif")
xylo_pres <- rast("output/Xylocopa.mexicanorum/pres.tif")
apis_fut <- rast("output/Apis.mellifera/mri_585_70.tif")
pasi_fut <- rast("output/Passiflora.foetida/mri_585_70.tif")
xylo_fut <- rast("output/Xylocopa.mexicanorum/mri_585_70.tif")

ov_pres_apis <- apis_pres+pasi_pres
ov_pres_xylo <- xylo_pres+pasi_pres

ov_fut_apis <- apis_fut+pasi_fut
ov_fut_xylo <- xylo_fut+pasi_fut

ov_pres_apis[ov_pres_apis!=2]<-NA
plot(ov_pres_apis)
pres_apis_pasi <- expanse(ov_pres_apis,"km")[2] %>% as.vector() %>% .[[1]]

ov_fut_apis[ov_fut_apis!=2]<-NA
plot(ov_fut_apis)
fut_apis_pasi <- expanse(ov_fut_apis,"km")[2] %>% as.vector() %>% .[[1]]

ov_pres_xylo[ov_pres_xylo!=2]<-NA
plot(ov_pres_xylo)
pres_xylo_pasi <- expanse(ov_pres_xylo,"km")[2] %>% as.vector() %>% .[[1]]

ov_fut_xylo[ov_fut_xylo!=2]<-NA
plot(ov_fut_xylo)
fut_xylo_pasi <-expanse(ov_fut_xylo,"km")[2] %>% as.vector() %>% .[[1]]


temp <- data.frame(overlap=
             c("pres_apis_pasi","fut_apis_pasi",
               "pres_xylo_pasi","fut_xylo_pasi"),
           area=
             c(pres_apis_pasi,fut_apis_pasi,
               pres_xylo_pasi,fut_xylo_pasi))

barplot(temp$area)
