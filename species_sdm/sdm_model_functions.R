#paquetes
#(apuntes toda america; tss y roc, wmeans, no binario)

library("biomod2")
library("dismo");library("dplyr");library("fuzzySim");
library("magrittr");library(terra); library(geodata)
library("readr");library("sp");library("tools"); library(MetBrewer)
library(grinnell)


# color
cols <- met.brewer("Hiroshige")
cols <- cols[10:1]
cols <- colorRampPalette(cols)(100)



get_input_depracated <- function(species,sps,extens,bioclima){
  occsData <- sps[grep(species,sps$species),c("decimalLongitude","decimalLatitude")]
  
  names(occsData) <- c("lon","lat")
  if(species=="Apis mellifera"){
    occsData <- occsData[sample(1:nrow(occsData),10000,F),]
  }
  
  occsData$species <- species
  occsData <- CoordinateCleaner::cc_cap(x = occsData,lon = "lon",lat = "lat")
  occsData <- CoordinateCleaner::cc_equ(x = occsData,lon = "lon",lat = "lat")
  occsData <- CoordinateCleaner::cc_cen(x = occsData,lon = "lon",lat = "lat")
  occsData <- CoordinateCleaner::cc_gbif(x = occsData,lon = "lon",lat = "lat")
  occsData <- CoordinateCleaner::cc_inst(x = occsData,lon = "lon",lat = "lat")
  occsData <- CoordinateCleaner::cc_sea(occsData,lon = "lon",lat = "lat")
  occsData <- CoordinateCleaner::cc_outl(x = occsData,lon = "lon",lat = "lat",
                                         species = "species",method = "quantile")
  occsData <- occsData[,c("lon","lat")] %>% vect()
  
  
  ecorregions <- terra::vect("~/Desktop/capas_climatica/ecoregions/wwf_terr_ecos.shp") 
  ecorregions <- terra::crop(ecorregions,extens) %>% magrittr::extract(.,occsData)
  bioclima <-terra::crop(bioclima,ecorregions) %>% mask(.,ecorregions)
  covarData <- terra::extract(bioclima, occsData,ID=FALSE) 
  select_var <- fuzzySim::corSelect(
    data = covarData, var.cols = names(covarData),cor.thresh = 0.8,
    use = "pairwise.complete.obs",select = "VIF")
  select_var <- select_var$selected.vars
  bioclima <- bioclima[[select_var]]
  
  plot(bioclima[[1]],col=cols,main=names(bioclima)[1])
  points(occsData)
  
  input <- list(occsData=occsData,myExpl=bioclima,myRespName=species,ecorregions=ecorregions)
  
  
  return(input)
}


get_input <- function(species,sps,extens,bioclima){
  occsData <- sps[grep(species,sps$species),c("lon","lat")]
  
  names(occsData) <- c("lon","lat")
  occsData <- vect(occsData)
  
  ecorregions <- terra::vect("~/Desktop/capas_climatica/ecoregions/wwf_terr_ecos.shp") 
  ecorregions <- terra::crop(ecorregions,extens) %>% magrittr::extract(.,occsData)
  bioclima <-terra::crop(bioclima,ecorregions) %>% mask(.,ecorregions)
  
  covarData <- terra::extract(bioclima, occsData,ID=FALSE) 
  select_var <- fuzzySim::corSelect(
    data = covarData, var.cols = names(covarData),cor.thresh = 0.8,
    use = "pairwise.complete.obs",select = "VIF")
  select_var <- select_var$selected.vars
  bioclima <- bioclima[[select_var]]
  
  plot(bioclima[[1]],col=cols,main=names(bioclima)[1])
  points(occsData)
  
  input <- list(occsData=occsData,myExpl=bioclima,myRespName=species,ecorregions=ecorregions)
  
  
  return(input)
}



model_sdm <- function(formato_input,modelo="MAXENT", perc=0.7, roc_tresh=c(0.5,0.8),N_run=10,override=T){
  myBiomodOption <-BIOMOD_ModelingOptions()
  myBiomodData <- BIOMOD_FormatingData(resp.var = formato_input$occsData,
                                       expl.var = formato_input$myExpl,na.rm = T,
                                       resp.name = formato_input$myRespName,
                                       PA.nb.rep = 1,
                                       PA.nb.absences = 10000,
                                       PA.strategy = 'random')
  myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,
                                      bm.options = myBiomodOption,
                                      models = c(modelo),
                                      modeling.id = 'AllModels',
                                      CV.nb.rep  = N_run,
                                      CV.perc = perc,data.split.perc = perc,
                                      var.import = 3,
                                      metric.eval = c("TSS", "ROC"),
                                      nb.cpu = 4,
                                      do.full.models = FALSE
  )
  myBiomodModelEval <- get_evaluations(myBiomodModelOut)
  write.csv(myBiomodModelEval, file = file.path(formato_input$myRespName, "myBiomodModelEval.csv"),row.names = FALSE)
  
  
  myBiomodEM <- BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
                                        models.chosen = 'all',
                                        em.by = 'all',
                                        metric.select = c('TSS','ROC'),
                                        metric.select.thresh = c(roc_tresh),
                                        var.import = 3,
                                        metric.eval = c('TSS','ROC'),
                                        em.algo = c('EMmean', 'EMwmean'),#,'EMcv', 'EMci', 'EMca'),
                                        EMwmean.decay = 'proportional')
  
  ##trabajar mejor con EMwmean
  if(length(myBiomodEM@em.models_kept)==1){
    if(myBiomodEM@em.models_kept=="none"){
      return("not_run")
      stop()
    }
    
  }
  
  myVarImportEM<-data.frame(get_variables_importance(myBiomodEM))
  write.csv(myVarImportEM, file = file.path(formato_input$myRespName, "myVarImportEM.csv"),row.names = T)
  
  myBiomodEMEval<-get_evaluations(myBiomodEM)
  write.csv(myBiomodEMEval, file = file.path(formato_input$myRespName, "myBiomodEMEval.csv"),
            row.names = FALSE)
  
  return(myBiomodEM)
}




project_model <- function(myBiomodEM,myExpl,name="present",formato_input){
  myBiomodEMProj <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM,
                                               proj.name = name,
                                               new.env = myExpl,
                                               models.chosen = 'all',
                                               metric.binary = 'all');myBiomodEMProj
  
  
  
  CurrentProj_ROC <-terra::rast(file.path(formato_input$myRespName,paste("proj_",name,"/proj_",name,"_",formato_input$myRespName,"_ensemble_ROCbin.tif",sep="")))
  CurrentProj_TSS <-terra::rast(file.path(formato_input$myRespName,paste("proj_",name,"/proj_",name,"_",formato_input$myRespName,"_ensemble_TSSbin.tif",sep="")))
  CurrentProj_prob <-terra::rast(file.path(formato_input$myRespName,paste("proj_",name,"/proj_",name,"_",formato_input$myRespName,"_ensemble.tif",sep="")))
  CurrentProj_ROC <- CurrentProj_ROC[[2]]
  CurrentProj_TSS <- CurrentProj_TSS[[2]]
  CurrentProj_prob <- CurrentProj_prob[[2]]/1000
  CurrentProj <- c(CurrentProj_prob,CurrentProj_ROC,CurrentProj_TSS)
  
  return(CurrentProj)
}
