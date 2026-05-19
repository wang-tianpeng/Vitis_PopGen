library(terra); library(tidyverse);library(ggpubr)

rm(list=ls());graphics.off()


  mask <- rast("inputs/wc2.1_2.5m_bio_1.tif")
  
  ops_tit <- data.frame(red=c("ig","pres","lgm"),
                        ent=c("Inter Glacial (~100 ky)",
                              "Present",
                              "Last Glacial Maximum (~20 ky)"))
  ops_tit <- ops_tit[which(ops_tit$red==tiempo),"ent"]
  ex <- ext(mask)
  ex[1] <- -150
  ex[2] <- -50
  ex[3] <- 0
  ex[4] <- 60
  mask <- crop(mask,ex)
  mask[!is.na(mask)]<-0
  res <- rast("input/wc2.1_30s_bio_1.tif")
  res <- crop(res,ext(mask))
  mask <- resample(x = mask,y =  res, method = "near")
  
  mapas <- list()
  cambio <- 1
  valores <- data.frame(sp=NA,val=NA)
  for(i in 1:length(sps)){
    print(paste("running",sps[i]))
    temp <- rast(paste("output/",sps[i],"/",tiempo,".tif",sep = ""))
    temp[temp==0]<-NA
    temp <- resample(temp, mask, method="bilinear")
    temp[is.na(temp)]<-0
    temp[temp>0]<- 1*cambio
    valores <- rbind(valores,data.frame(sp=sps[i],
                                        val=1*cambio))
    if(plot_all){
      plot(temp+mask,main=sps[i])
    }
    mapas[[i]]<-temp
    names(mapas)[i]<- sps[i]
    cambio <- cambio*10
  }
  
  mapas <- rast(mapas)
  mapas <- sum(mapas)
  mapas <- mapas+mask
  modelo <- as.factor(mapas)
  
  if(plot_all){
    plot(modelo)
  }
  
  nex <- as.data.frame(modelo,xy=T)
  nex <- nex[which(nex$sum!=0),]
  n_ext <- ex
  if(nrow(nex)!=0){
    n_ext[1] <- min(nex$x)-2
    n_ext[2] <- max(nex$x)+2
    n_ext[3] <- min(nex$y)-2
    n_ext[4] <- max(nex$y)+2
    if(plot_all){
      lines(n_ext)
    }
    
    modelo <- crop(modelo,n_ext)
  }else {
    modelo <- crop(modelo,n_ext)
  }
  
  
  valores <- valores[-1,]
  
  valores <- map_dfr(2:length(valores$sp), function(k) {
    combn(valores$sp, k, simplify = FALSE) %>%
      map_dfr(~{
        especies <- .x
        suma_valores <- sum(valores$val[valores$sp %in% especies])
        tibble(
          sp = paste(especies, collapse = " + "),
          val = suma_valores
        )
      })
  }) %>% as.data.frame() %>% rbind(valores,.)
  valores <- rbind(valores,data.frame(sp="Absence",val=0))
  
  normaliza_nombre <- function(nombre) {
    if (nombre == "Absence") return("Absence")
    especies <- unlist(strsplit(nombre, "\\+|\\s+"))
    especies <- especies[especies != ""]
    especies <- sort(especies)
    paste(especies, collapse = "+")
  }
  
  valores$sp_norm <- vapply(valores$sp, normaliza_nombre, character(1))
  
  valores$color <- cols[valores$sp_norm]
  
  ja <- valores
  names(ja) <- c("etiqueta","ID")
  
  levels(modelo) <- ja[,c(2:1)]
  #plot(modelo)
  
  cols <- data.frame(valores,row.names = "sp")
  
  df <- as.data.frame(modelo, xy = TRUE, cells = TRUE, na.rm = FALSE)
  names(df)[which(!names(df) %in% c("x", "y", "cell"))] <- "species"
  cols <- unique(df$species) %>% .[!is.na(.)] %>% .[order(.)] %>% cols[.,]
  ggplot(df) +
    geom_raster(aes(x = x, y = y, fill = species)) +
    coord_equal() +
    #scale_fill_viridis_d(option = "turbo") + 
    scale_fill_manual(values = cols$color,
                      na.value = "white",
                      guide = guide_legend(na.translate = FALSE),) +
    labs(title = ops_tit)+
    #  scale_fill_brewer(palette = "Set3", na.value = "white") +  # o viridis, etc.
    theme_void()+theme(legend.position = "bottom") -> plt
  if(plot_all){
  print(plt)
  }
  
  return(list(df=df,plt=plt,n_ext=n_ext))
}


colores_combinaciones_especies <- function(especies) {
  library(gtools)
  colores <- colors <- c(
    "#66c2a5",  # verde menta
    "#fc8d62",  # naranja claro
    "#8da0cb",  # azul lavanda
    "#e78ac3",  # rosa
    "
    "#ffd92f",  # amarillo
    "#e31a1c",   # rojo fuerte
    "#1f78b4",  # azul fuerte
    "
  )
  
  n <- length(especies)
  todas_combinaciones <- list()
  
  for (i in 1:n) {
    combs <- combinations(n, i, especies)
    etiquetas <- apply(combs, 1, function(x) paste(sort(x), collapse = "+"))
    todas_combinaciones <- c(todas_combinaciones, etiquetas)
  }
  
  n_combs <- length(todas_combinaciones)
  
  colores <- sample(colores,n_combs,replace = F)
  
  names(colores) <- todas_combinaciones
  colores <- c(colores,Absence="#e0e0e0")
  return(colores)
}


pdf("output/test2.pdf")

sps <- list.files(path = "output/")



sps_t <- sps[c(8,5,1)]
cols <- colores_combinaciones_especies(especies = sps_t)
tiempo <- c("pres","ig","lgm")
pres <- mapear(sps = sps_t,tiempo = tiempo[1],cols = cols)
lgm <- mapear(sps = sps_t,tiempo = tiempo[3],cols = cols)
ig <- mapear(sps = sps_t,tiempo = tiempo[2],cols=cols)
ggarrange(pres$plt,lgm$plt,ig$plt,ncol = 1,common.legend = T)


sps_t <- sps[c(10,9,1)]
cols <- colores_combinaciones_especies(especies = sps_t)
tiempo <- c("pres","ig","lgm")
pres <- mapear(sps = sps_t,tiempo = tiempo[1],cols = cols)
lgm <- mapear(sps = sps_t,tiempo = tiempo[3],cols = cols)
ig <- mapear(sps = sps_t,tiempo = tiempo[2],cols=cols)
ggarrange(pres$plt,lgm$plt,ig$plt,ncol = 1,common.legend = T)



sps_t <- sps[c(4,8,7)]
cols <- colores_combinaciones_especies(especies = sps_t)
tiempo <- c("pres","ig","lgm")
pres <- mapear(sps = sps_t,tiempo = tiempo[1],cols = cols)
ig <- mapear(sps = sps_t,tiempo = tiempo[2],cols=cols)
ggarrange(pres$plt,ig$plt,ncol = 1,common.legend = T)

sps_t <- sps[c(2,6)]
cols <- colores_combinaciones_especies(especies = sps_t)
tiempo <- c("pres","ig","lgm")
pres <- mapear(sps = sps_t,tiempo = tiempo[1],cols = cols)
lgm <- mapear(sps = sps_t,tiempo = tiempo[3],cols = cols)
ig <- mapear(sps = sps_t,tiempo = tiempo[2],cols=cols)
ggarrange(pres$plt,lgm$plt,ig$plt,ncol = 1,common.legend = T)


dev.off()









  mask <- rast("input/wc2.1_2.5m_bio_1.tif")
  
  ops_tit <- data.frame(red=c("ig","pres","lgm"),
                        ent=c("Inter Glacial (~100 ky)",
                              "Present",
                              "Last Glacial Maximum (~20 ky)"))
  ex <- ext(mask)
  ex[1] <- -150
  ex[2] <- -50
  ex[3] <- 0
  ex[4] <- 60
  mask <- crop(mask,ex)
  mask[!is.na(mask)]<-0
  res <- rast("input/wc2.1_30s_bio_1.tif")
  res <- crop(res,ext(mask))
  mask <- resample(x = mask,y =  res, method = "near")
  
  mapas <- list()
  cambio <- 1
  valores <- data.frame(sp=NA,val=NA)
  l_sps <- list.files(path = paste("output/",sps,sep = ""),pattern = ".tif",full.names = T)
  for(i in 1:length(l_sps)){
    print(paste("running",l_sps[i]))
    temp <- rast(l_sps[i])
    temp[temp==0]<-NA
    temp <- resample(temp, mask, method="bilinear")
    temp[is.na(temp)]<-0
    temp[temp>0]<- 1*cambio
    valores <- rbind(valores,data.frame(sp=l_sps[i],
                                        val=1*cambio))
    if(plot_all){
      plot(temp+mask,main=l_sps[i])
    }
    mapas[[i]]<-temp
    names(mapas)[i]<- l_sps[i]
    cambio <- cambio*10
  }
  
  mapas <- rast(mapas)
  mapas <- sum(mapas)
  mapas <- mapas+mask
  modelo <- as.factor(mapas)
  
  if(plot_all){
    plot(modelo)
  }
  
  nex <- as.data.frame(modelo,xy=T)
  nex <- nex[which(nex$sum!=0),]
  n_ext <- ex
  if(nrow(nex)!=0){
    n_ext[1] <- min(nex$x)-2
    n_ext[2] <- max(nex$x)+2
    n_ext[3] <- min(nex$y)-2
    n_ext[4] <- max(nex$y)+2
    if(plot_all){
      lines(n_ext)
    }
    
    modelo <- crop(modelo,n_ext)
  }else {
    modelo <- crop(modelo,n_ext)
  }
  
  valores <- data.frame(sp=c("Absence","ig","lgm","pres",
                             "ig+lgm","ig+pres","lgm+pres",
                             "ig+lgm+pres"),
                        val=c(0,1,10,100,11,101,110,111))
  


  valores$cols <- c(
    "#e0e0e0",
    "#66c2a5",  # verde menta
    "#fc8d62",
    "#ffd92f",  # azul lavanda
    "#e78ac3",  # rosa
    "
    "#8da0cb",  # amarillo
    "#e31a1c"  # rojo fuerte
     # morado
  )

  ja <- valores
  names(ja) <- c("etiqueta","ID","cols")
  
  levels(modelo) <- ja[,c(2:1)]
  #plot(modelo)
  

  df <- as.data.frame(modelo, xy = TRUE, cells = TRUE, na.rm = FALSE)
  names(df)[which(!names(df) %in% c("x", "y", "cell"))] <- "species"
  cols <- valores$cols
  ggplot(df) +
    geom_raster(aes(x = x, y = y, fill = species)) +
    coord_equal() +
    #scale_fill_viridis_d(option = "turbo") + 
    scale_fill_manual(values = cols,
                      na.value = "white",
                      guide = guide_legend(na.translate = FALSE),) +
    labs(title = sps)+
    #  scale_fill_brewer(palette = "Set3", na.value = "white") +  # o viridis, etc.
    theme_void()+theme(legend.position = "bottom") -> plt
  if(plot_all){
    print(plt)
  }
  
  return(list(df=df,plt=plt,n_ext=n_ext))
}
plt <- mapear_una(sps = "californica",plot_all = F)
pdf("output/calif.pdf")
print(plt$plt)
dev.off()
