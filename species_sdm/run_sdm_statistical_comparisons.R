#centroids
rm(list=ls()); graphics.off()
library(terra)
library(geodata)
library(tidyverse)
library(maps)

if(FALSE){
  mask <- rast("inputs/wc2.1_2.5m_bio_1.tif")
  ex <- ext(mask)
  ex[1] <- -150
  ex[2] <- -50
  ex[3] <- 0
  ex[4] <- 60
  mask <- crop(mask,ex)
  mask[!is.na(mask)]<-0
  res <- rast(list.files(path = "inputs/wc2.1_30s_bio/",pattern = "tif",full.names = T))
  res <- crop(res,ext(mask))
  names(res) <- sub("wc2.1_30s_","",names(res))
  res <- res[[paste("bio_",1:19,sep = "")]]
  mask <- resample(x = mask,y =  res[[1]], method = "near")
  
  l_sps <- list.files(path = "output/",pattern = "Vitis")
  
}



if(FALSE){
  for(i in 1:length(l_sps)){
    ras <- rast(paste("output/",l_sps[i],"/pres.tif",sep = ""))
    ras_prob <- ras[[1]]
    ras <- ras[[3]]
    ras[ras==0]<-NA
    
    plot(ras)
    geos <- as.data.frame(ras,xy=T)[,c("x","y")]
    pol <- as.polygons(ras, dissolve = TRUE)
    centroid <- centroids(pol)
    points(centroid)
    
    bioclim <- terra::extract(x = res,y = geos,ID=FALSE)
    
    pca <- stats::prcomp(bioclim,retx = T,center = T,scale. = T)
    pca <-pca$x[,1:6] 
    cent_nich <- colMeans(pca)
    pca <- data.frame(pca)
    
    geos$niche_cent <- (pca- cent_nich)^2 %>% rowSums() %>% sqrt()
    
    cent_geo <- geom(centroid)[,c("x","y")]
    cent_geo
    geos$geo_cent <- sqrt((geos$x- cent_geo[1])^2 + (geos$y- cent_geo[2])^2)
    
    ras_niche <- rast(geos[,c("x","y","niche_cent")],type="xyz")
    ras_geo <- rast(geos[,c("x","y","geo_cent")],type="xyz")
    
    ras <- resample(ras, mask, method="bilinear")
    ras_prob <- resample(ras_prob, mask, method="bilinear")
    ras_geo <- resample(ras_geo, mask, method="bilinear")
    ras_niche <- resample(ras_niche, mask, method="bilinear")
    
    ras_salida <- c(ras,ras_prob,ras_geo,ras_niche)
    names(ras_salida)[1:2] <- strsplit(names(ras_salida)[1:2],split = "_",fixed = T) %>% sapply(.,"[",1)
    writeRaster(ras_salida,filename = paste("output/centroids/",l_sps[i],".tif",sep = ""),overwrite=TRUE)
    
    
  }
  
}


if(FALSE){
  mask <- rast("inputs/mask.tif")
  pdf("output/centroids/centroids.pdf")
  lsp <- list.files("output") %>% .[grep("Vitis",.)]
  output <- data.frame(sp=NA,period=NA,lon=NA,lat=NA,val=NA)
  
  for(i in 1:length(lsp)){
    
    ras <- rast(paste("output/",lsp[i],"/pres.tif",sep = ""))
    ras <- ras[[3]]
    ras[ras==0]<-NA
    tp <- values(ras) %>% na.omit() %>% t() %>% as.vector()
    if(length(tp)<2){
      val="FALSE"
    }else{
      val="TRUE"
    }
    
    pol <- as.polygons(ras, dissolve = TRUE)
    if((nrow(values(pol))==0)){
      temp <- data.frame(sp=lsp[i],period="pres",lon=NA,lat=NA,val=val)
    }else{
      centroid <- centroids(pol)
      centroid <- geom(centroid)[,c("x","y")]
      temp <- data.frame(sp=lsp[i],period="pres",lon=centroid["x"],lat=centroid["y"],val=val)
    }
    output <- rbind(output,temp)
    ras_pres <- resample(ras, mask, method="bilinear")
    
    ras <- rast(paste("output/",lsp[i],"/lgm.tif",sep = ""))
    ras <- ras[[3]]
    ras[ras==0]<-NA
    tp <- values(ras) %>% na.omit() %>% t() %>% as.vector()
    if(length(tp)<2){
      val="FALSE"
    }else{
      val="TRUE"
    }
    
    pol <- as.polygons(ras, dissolve = TRUE)
    if((nrow(values(pol))==0)){
      temp <- data.frame(sp=lsp[i],period="lgm",lon=NA,lat=NA,val=val)
    }else{
      centroid <- centroids(pol)
      centroid <- geom(centroid)[,c("x","y")]
      temp <- data.frame(sp=lsp[i],period="lgm",lon=centroid["x"],lat=centroid["y"],val=val)
    }
    output <- rbind(output,temp)
    ras_lgm <- resample(ras, mask, method="bilinear")
    
    
    ras <- rast(paste("output/",lsp[i],"/ig.tif",sep = ""))
    ras <- ras[[3]]
    ras[ras==0]<-NA
    tp <- values(ras) %>% na.omit() %>% t() %>% as.vector()
    if(length(tp)<2){
      val="FALSE"
    }else{
      val="TRUE"
    }
    
    pol <- as.polygons(ras, dissolve = TRUE)
    if((nrow(values(pol))==0)){
      temp <- data.frame(sp=lsp[i],period="ig",lon=NA,lat=NA,val=val)
    }else{
      centroid <- centroids(pol)
      centroid <- geom(centroid)[,c("x","y")]
      temp <- data.frame(sp=lsp[i],period="ig",lon=centroid["x"],lat=centroid["y"],val=val)
    }
    output <- rbind(output,temp)
    ras_ig <- resample(ras, mask, method="bilinear")
    
    ras <- c(ras_pres,ras_lgm,ras_ig)
    names(ras) <- paste(lsp[i],c("pres","lgm","ig"),sep = ":")
    plot(ras)
    
    print(output)
    
    
    
  }
  dev.off()
  
  save(output,file = "output/centroids/cambio.R") 
  
  ggplot(output,mapping = aes(lon,lat,col=sp,shape=period))+
    geom_point()+
    theme_classic()



if(FALSE){
  mask <- rast("inputs/mask.tif")
  temp <- ext(mask)
  temp[1:4] <- c(-124.3875, -55.8125, 19.0791666666667, 59.3458333333333)
  mask <- crop(mask,temp)
  
  mapa <- gadm(country = c("US","Mex"),level = 0,path = tempdir())
  mapa <- crop(mapa,ext(mask))
  
  
  
  l_sps <- list.files(path = "output/",pattern = "Vitis")
  for(i in 1:length(l_sps)){
    ras <- rast(paste("output/",l_sps[i],"/pres.tif",sep = ""))
    ras <- ras[[1]]
    #ras[ras!=1] <- 0
    ras <- resample(ras, mask, method="bilinear")
    ras[is.na(ras)] <- 0
    ras <- ras+mask
    assign(l_sps[i],ras)
  }
  
  
  library(readxl)
  admixture <- read.csv("introgression/VitisProject/ForJonas/samples_SDM_INFO.csv")
  datos <- read_excel(path = "introgression/VitisProject/ForJonas/samples_SDM_INFO.xlsx",sheet = 1)%>% as.data.frame()
  admixture <- data.frame(datos,admixture[,-1])
  admixture <- admixture[-which(is.na(admixture$californica)),]
  admixture <- admixture[!duplicated(admixture[,-c(1:2)]),]
  
  Q <- admixture[,-c(1:2)]
  
  Q <- as.matrix(Q)
  K <- ncol(Q)
  n <- nrow(Q)
  
  eps <- 1e-12  # para evitar log(0)
  maxprop <- apply(Q, 1, max)
  admix_1minusmax <- 1 - maxprop
  
  H_raw <- -rowSums(Q * log(Q + eps))         # en nats
  H_norm <- H_raw / log(K)                    # normalizada a 0-1
  
  simpson_raw <- 1 - rowSums(Q^2)
  simpson_norm <- simpson_raw / (1 - 1/K)     # normaliza a 0-1 (cuando K>1)
  
  effN <- exp(H_raw)  # rango [1, K]
  
  onehots <- diag(1, K, K)
  dist_to_pure <- apply(Q, 1, function(q) {
    dists <- apply(onehots, 1, function(e) sqrt(sum((q - e)^2)))
    min(dists)
  })
  
  df <- data.frame(
    id = admixture$Sample_ID,
    maxprop = maxprop,
    admix_1minusmax = admix_1minusmax,
    H_raw = H_raw,
    H_norm = H_norm,
    simpson_raw = simpson_raw,
    simpson_norm = simpson_norm,
    effN = effN,
    dist_to_pure = dist_to_pure,
    stringsAsFactors = FALSE
  )
  
  if(FALSE){
    ### dos datos
    
    r1 <- raster(Vitis.shuttleworthii)
    r2 <- raster(Vitis.californica)
    
    v1 <- values(r1)
    v2 <- values(r2)
    
    v1[is.na(v1)] <- 0
    v2[is.na(v2)] <- 0
    
    # normalizar
    p1 <- v1 / sum(v1)
    p2 <- v2 / sum(v2)
    
    # Schoener's D
    D <- 1 - 0.5 * sum(abs(p1 - p2))
    
  }
  
  
  ### mas de dos datos
  l_sps
  # rasters[rasters==0]<-NA
  # temp[1:2] <- x
  # temp[3:4] <- y
  
  
  salida_over <- vector(mode = "list",length = nrow(admixture))
  df$D <- NA
  df$I <- NA
  
  temp_admixture <- as.matrix(admixture[,-c(1:2)])
  temp_admixture[which(temp_admixture>0)] <- 1
  
  temp_admixture <- data.frame(ID=admixture$Sample_ID,temp_admixture)
  
  for(i in 1:nrow(admixture)){
    if(!is.na(df$D[i])){
      next
    }
    print(i)
    temp <- admixture[i,-c(1:2)] %>% .[,.!=0] ; print(temp)
    if(length(temp)==1){
      next
    }
    sps <- paste("Vitis",".",names(temp),sep = "")
    sps <- intersect(l_sps,sps)
    if(length(sps)==1 | length(sps)==0){
      next
    }
    rasters <- rast(lapply(sps, get))
    
    vals <- sapply(rasters, function(x) {
      v <- values(x)
      v[is.na(v)] <- 0
      v
    })
    
    row_sum <- rowSums(vals)
    vals <- vals[which(row_sum!=0),]
    
    col_sums <- colSums(vals)
    P <- sweep(vals, 2, col_sums, FUN = "/")  # cada columna suma 1
    
    min_per_cell <- apply(P, 1, min)
    D_min_all <- sum(min_per_cell)   # valor entre 0 y 1
    
    D_pair <- function(p, q) 1 - 0.5 * sum(abs(p - q))
    K <- ncol(P)
    pair_vals <- combn(K, 2, function(idx) {
      D_pair(P[, idx[1]], P[, idx[2]])
    })
    avg_pairwise_D <- mean(pair_vals)
    pairwise_matrix <- matrix(NA, K, K)
    comb <- combn(K, 2)
    for (j in seq_len(ncol(comb))) {
      i1 <- comb[1, j]; i2 <- comb[2, j]
      pairwise_matrix[i1, i2] <- pairwise_matrix[i2, i1] <- pair_vals[j]
    }
    diag(pairwise_matrix) <- 1
    
    geom_per_cell <- apply(P, 1, function(x) prod(x)^(1 / K))
    I_geom_all <- sum(geom_per_cell)  # entre 0 y 1
    
    overs <- list(
      D_min_all = D_min_all,
      avg_pairwise_D = avg_pairwise_D,
      pairwise_matrix = pairwise_matrix,
      I_geom_all = I_geom_all
    )
    
    x <- as.vector(as.matrix(temp_admixture[i, -1]))
    
    tp <- which(apply(temp_admixture[,-1], 1, function(z) all(as.vector(as.matrix(z)) == x)))
    
    for(v in tp){
      salida_over[[v]] <- overs
      names(salida_over)[v] <- admixture$Sample_ID[v]
      
    }
    
    df[tp,"D"] <- avg_pairwise_D
    df[tp,"I"] <- I_geom_all
    
  }
  
  df[,c("id","D","maxprop")]
  
  ja <- df[which(!is.na(df$D)),]
  
  metrics <- c("admix_1minusmax", "H_norm", "simpson_norm", "effN", "dist_to_pure")
  cols <- c("maxprop", "admix_1minusmax", "H_raw", "H_norm", "simpson_raw","effN", "dist_to_pure")
  
  cor_D <- cor(ja[, "I"], ja[, cols], use = "complete.obs")
  
  cor_D
  
  write.csv(df,file = "output/centroids/df.csv")
  
}




############ solo sobrelapes considerando un grado de admixture ########################

if(TRUE){
  corte <- 0.1
  mask <- rast("inputs/mask.tif")
  temp <- ext(mask)
  temp[1:4] <- c(-124.3875, -55.8125, 19.0791666666667, 59.3458333333333)
  mask <- crop(mask,temp)
  
  
  
  l_sps <- list.files(path = "output/",pattern = "Vitis")
  for(i in 1:length(l_sps)){
    ras <- rast(paste("output/",l_sps[i],"/pres.tif",sep = ""))
    ras <- ras[[2]]
    #ras[ras!=1] <- 0
    ras <- resample(ras, mask, method="bilinear")
    ras[is.na(ras)] <- 0
    ras <- ras+mask
    assign(l_sps[i],ras)
  }
  
  
  library(readxl)
  admixture <- read.csv("introgression/VitisProject/ForJonas/samples_SDM_INFO.csv")
  datos <- read_excel(path = "introgression/VitisProject/ForJonas/samples_SDM_INFO.xlsx",sheet = 1)%>% as.data.frame()
  admixture <- data.frame(datos,admixture[,-1])
  admixture <- admixture[-which(is.na(admixture$californica)),]
  admixture <- admixture[!duplicated(admixture[,-c(1:2)]),]
  
  Q <- admixture[,-c(1:2)]
  
  Q <- as.matrix(Q)
  pops <- admixture[,1:2]
  pops$sobrelape <- NA
  ### mas de dos datos
  l_sps
  Q[Q>corte] <- 1
  Q[Q<1]<-0
  for(i in 1:nrow(pops)){
    print(i)
    temp <- Q[i,] %>% .[.>0] ; print(temp)
    if(length(temp)==1){
      pops[i,"sobrelape"] <- 0
      next
    }
    sps <- paste("Vitis",".",names(temp),sep = "")
    sps <- intersect(l_sps,sps)
    if(length(sps)==1 | length(sps)==0){
      next
    }
    rasters <- rast(lapply(sps, get))
    rasters <- sum(rasters)
    rasters[rasters!=2]<-NA
    pops[i,"sobrelape"] <- expanse(rasters,unit="km")[[2]] %>% as.vector()
  }
  
  
  df <- read.csv("output/centroids/df.csv",row.names = "id")
  pops <- data.frame(pops,row.names = "Sample_ID")
  rownames(df)==rownames(pops)
  pops <- data.frame(df,pops)
  write.csv(pops,file = "output/centroids/df_overlaps.csv")
    
    
}

mask <- rast("inputs/mask.tif")
mapear <- function(sps,plot_all=F,n_ext){
  
  
  uno <- get(sps[1])
  dos <- get(sps[2])
  dos[dos==1]<-2
  over <- uno+dos
  modelo <- as.factor(over)
  
  
  valores <- data.frame(sp=sps,val=c(1,2))
  
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
  
 
  
  ja <- valores
  names(ja) <- c("etiqueta","ID")
  
  levels(modelo) <- ja[,c(2:1)]
  #plot(modelo)
  
  df <- as.data.frame(modelo, xy = TRUE, cells = TRUE, na.rm = FALSE)
  names(df)[which(!names(df) %in% c("x", "y", "cell"))] <- "species"
  return(modelo)
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



load("centroids/cambio.R")
l_sps <- list.files(path = "output/",pattern = "Vitis")
for(i in 1:length(l_sps)){
  ras <- rast(paste("output/",l_sps[i],"/pres.tif",sep = ""))
  ras <- ras[[2]]
  #ras[ras!=1] <- 0
  ras <- resample(ras, mask, method="bilinear")
  ras[is.na(ras)] <- 0
  ras <- ras+mask
  assign(l_sps[i],ras)
}

opciones <- combn(l_sps,2)
output <- output[output$period=="pres",]
output <- output[-1,]
output <- data.frame(output,row.names = "sp")

library(geosphere)

salida <- data.frame(sp1=NA, sp2=NA, Area_sp1=NA,Area_sp2=NA,Area_Overlap=NA,dist_cent=NA)
pdf("output/mapas_overlap.pdf")
i <- 1
for(i in 1:ncol(opciones)){
  print(ncol(opciones)-i)
  sps <- opciones[,i]
  mapa <- mapear(sps = sps)
  
  plot(mapa,col=c("lightgrey","#e78ac3","#a6d854","#ffd92f"),main=paste0(sps,collapse = "+"))
  temp <- output[sps,c("lon","lat")]
  km <- distGeo(temp[1,],temp[2,])/1000
  points(temp,pch="+")
  over <- mapa
  over[over!=3] <- NA
  over <- expanse(over,unit="km")
  over <- as.vector(over[[2]])
  
  uno <- mapa
  uno[uno!=1] <- NA
  uno <- expanse(uno,unit="km")
  uno <- as.vector(uno[[2]])
  
  dos <- mapa
  dos[dos!=2] <- NA
  dos <- expanse(dos,unit="km")
  dos <- as.vector(dos[[2]])
  
  temp <- data.frame(sp1=sps[1], sp2=sps[2],
                     Area_sp1=uno+over,
                     Area_sp2=dos+over,
                     Area_Overlap=over,
                     dist_cent=km)
  salida <- rbind(salida,temp)
  
}
dev.off()


salida <- salida[-1,]
write.csv(x = salida,file = "output/centroids/areas_dist.csv")




#################### analisis general ##############


# rm(list=ls())
cols <- MetBrewer::met.brewer("Hiroshige")[10:1]
pal <- colorRampPalette(cols)


pdf("./Vitis_centroids_all.pdf", width = 12, height = 8)
files <- list.files("centroids", pattern = "Vitis.*\\.tif$", full.names = TRUE)

my_xlim <- c(-150, -50)
my_ylim <- c(15, 70)

for(f in files){
  ras <- rast(f)
  
  ras <- crop(ras, crop_ext)
  
  par(mfrow=c(2,2), mar=c(3, 3, 2, 4), oma=c(1, 1, 3, 1))
  
  layer_names <- names(ras)
  
    plot(NULL, xlim=my_xlim, ylim=my_ylim, 
         axes=FALSE, xlab="", ylab="", type="n", xaxs="i", yaxs="i",
         main = layer_names[i])
    if(require(maps)) {
      map("world", fill=TRUE, col="lightgrey", border="white", add=TRUE)
    }
    
    axis(1); axis(2); box()
    
    plot(ras[[i]], col=pal(100), add=TRUE, legend=TRUE,
         plg=list(x="right", bty="n"))
  
  for(i in 3:nlyr(ras)){
    plot(NULL, xlim=my_xlim, ylim=my_ylim, 
         axes=FALSE, xlab="", ylab="", type="n", xaxs="i", yaxs="i",
         main = layer_names[i])
    if(require(maps)) {
      map("world", fill=TRUE, col="lightgrey", border="white", add=TRUE)
    }
    
    axis(1); axis(2); box()
    
    plot(ras[[i]], col=rev(pal(100)), add=TRUE, legend=TRUE,
         plg=list(x="right", bty="n"))

  mtext(basename(f), outer=TRUE, cex=1.5, font=2)
}

dev.off()
par(mfrow=c(1,1))









overlaps <- read.csv("centroids/areas_dist.csv")
flextable::flextable(overlaps[c(1:5,15:20,35:40),-1])









load("centroids/cambio.R")
flextable::flextable(output[c(2:4,14:16,35:40),])


admixture <- read.csv("centroids/df_overlaps.csv")
admixture[sample(1:nrow(admixture),10,F),c("X.1","D","I","sobrelape")]

names(admixture)[14] <- "overlap"

names(admixture)[c(3,5,7,9)]
# "maxprop"     "H_raw"       "simpson_raw" "effN"

names(admixture)[c(3,5,7,9)] <- c( "Max_Prop","Shannon_H","Simpson_Index","Eff_N")

variables <- names(admixture)[c(3,5,7,9)]

pdf("./Comparison_statistics.pdf",width = 21,height = 5)

point_col <- "#528fad" 
line_col <- "#e74c3c"

for(i in variables){
  
  fit <- lm(admixture[,i]~admixture$D) %>% summary() 
  p_val <- fit$coefficients[2,4]
  R <- fit$r.squared %>% round(.,2)
  
  if(p_val < 0.001){
  } else {
  }
  
  a <- ggplot(admixture, mapping = aes_string("D", i)) +
  geom_point(alpha = 0.6, color = point_col, size = 2) +
  geom_smooth(method = "lm", color = line_col, fill = line_col, alpha = 0.2) +
  labs(x = "D (Suitability Overlap)", y = i, subtitle = stats_text) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    plot.subtitle = element_text(face = "italic", size = 10)
  )
  
  fit <- lm(admixture[,i]~admixture$I) %>% summary() 
  p_val <- fit$coefficients[2,4]
  R <- fit$r.squared %>% round(.,2)
  
  if(p_val < 0.001){
  } else {
  }
  
  b <- ggplot(admixture, mapping = aes_string("I", i)) +
  geom_point(alpha = 0.6, color = point_col, size = 2) +
  geom_smooth(method = "lm", color = line_col, fill = line_col, alpha = 0.2) +
  labs(x = "I (Suitability Overlap)", y = "", subtitle = stats_text) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    plot.subtitle = element_text(face = "italic", size = 10)
  )
  
  fit <- lm(admixture[,i]~admixture$overlap) %>% summary() 
  p_val <- fit$coefficients[2,4]
  R <- fit$r.squared %>% round(.,2)
  
  if(p_val < 0.001){
  } else {
  }
  
  c <- ggplot(admixture, mapping = aes_string("overlap", i)) +
  geom_point(alpha = 0.6, color = point_col, size = 2) +
  geom_smooth(method = "lm", color = line_col, fill = line_col, alpha = 0.2) +
  labs(x = "Geographic Overlap", y = "", subtitle = stats_text) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    plot.subtitle = element_text(face = "italic", size = 10)
  )
  
  print(ggpubr::ggarrange(a, b, c, ncol = 3))
}
dev.off()


aba_baba <- read_excel("introgression/VitisProject/ForJonas/samples_SDM_INFO.xlsx",sheet = 4)
aba_baba %>%  filter(`P-value`<0.001) %>% select(`Z-score`) %>% range()
plot(aba_baba$`D statistic`,aba_baba$`Z-score`)


areas_dist <- read.csv("centroids/areas_dist.csv")
aba_baba <- as.data.frame(aba_baba[,c("P2","P3","D statistic")])
aba_baba$P2 <-sub("x_","",aba_baba$P2) %>%  paste("Vitis.",.,sep = "")
aba_baba$P3 <-sub("x_","",aba_baba$P3) %>%  paste("Vitis.",.,sep = "")
aba_baba$dist_cent <- NA
aba_baba$Area_Overlap <- NA
test <- areas_dist[,c("sp1","sp2")]
names(test) <- c("P2","P3")
for(i in 1:nrow(aba_baba)){
  temp <- aba_baba[i,1:2]
  key1 <- apply(temp[, c("P2", "P3")], 1, function(x) paste(sort(x), collapse = "_"))
  key2 <- apply(test[, c("P2", "P3")], 1, function(x) paste(sort(x), collapse = "_"))
  n <- rownames(test[key2 == key1, ]) %>% as.numeric()
  if(length(n)>0){
    print(n)
    aba_baba[i,"dist_cent"] <- areas_dist[n,"dist_cent"]
    aba_baba[i,"Area_Overlap"] <- areas_dist[n,"Area_Overlap"]
    }
  
}


plot(aba_baba$Area_Overlap,aba_baba$`D statistic`)
plot(aba_baba$dist_cent,aba_baba$`D statistic`)
ggplot(aba_baba,mapping = aes(Area_Overlap,`D statistic`))+
  geom_point()+geom_smooth(method = "lm")+theme_classic()->a

ggplot(aba_baba,mapping = aes(dist_cent,`D statistic`))+
  geom_point()+geom_smooth(method = "lm")+theme_classic()->b

ggpubr::ggarrange(a,b,nrow = 2)
