#function
rm(list=ls()); graphics.off()
library(terra); 
library(tidyverse)
library(maps)



Vitis.acerifolia
Vitis.aestivalis
Vitis.arizonica
Vitis.californica
Vitis.champinii
Vitis.cinerea
Vitis.giradiana
Vitis.labrusca
Vitis.monticola
Vitis.mustangensis
Vitis.riparia
Vitis.rotundifolia
Vitis.rupestris
Vitis.shutteworthii
Vitis.vulpina

path <- "./centroids/"
names_species <- c("Vitis.mustangensis","Vitis.acerifolia", "Vitis.riparia")

names_species <- c("Vitis.arizonica")


create_map <- function(path,names_species,extension=FALSE){
  final <- mask <- rast("inputs/mask.tif")
  values <- c()
  cambio <- 1
  for(i in 1:length(names_species)){
    sp <- paste(path,names_species[i],".tif",sep = "")
    ras <- rast(sp)
    ras <- ras[[1]]
    ras <- resample(ras, mask, method="bilinear")
    ras[ras==1]<-(cambio)
    values <- c(values,cambio)
    ras[is.na(ras)] <- 0
    ras <- ras+mask
    final <- c(final,ras)
    cambio <- cambio*10
  }
  final <- final[[2:nlyr(final)]]
  overlap <- sum(final)
  
  overlap <- as.factor(overlap)

  valores <- data.frame(sp=names_species,val=values)
  
  valores <- map_dfr(1:length(valores$sp), function(k) {
    combn(valores$sp, k, simplify = FALSE) %>%
      map_dfr(~{
        especies <- .x
        suma_valores <- sum(valores$val[valores$sp %in% especies])
        tibble(
          sp = paste(especies, collapse = " + "),
          val = suma_valores
        )
      })
  }) %>% as.data.frame()
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
  
  levels(overlap) <- ja[,c(2:1)]
  
  
  if(extension){
    extension <- overlap
    extension[extension==0]<-NA
    ja <- as.data.frame(extension,xy=T)
    x <- range(ja$x)
    y <- range(ja$y)
    extension <- ext(extension)
    extension[c(1,2)] <- x
    extension[c(3,4)] <- y
    overlap <- crop(overlap,extension)
    
  }
  
  return(overlap)
  
}


ras <- create_map(path = path,names_species = names_species)
ras_ext <- create_map(path = path,names_species = names_species,extension = TRUE)

colors <- c(
  "#66c2a5", 
  "#fc8d62", 
  "#8da0cb", 
  "#e78ac3", 
  "#a6d854", 
  "#ffd92f", 
  "#e31a1c", 
  "#1f78b4", 
  "#b15928", 
  "#6a3d9a" 
)


n <- values(ras) %>% unique() %>% as.vector() %>% na.omit() %>% .[.>0] %>% length()
cols <- sample(colors,n,F) %>% c("lightgrey",.)
plot(ras,col=cols)
plot(ras_ext,col=cols)


plot(ras, col=cols, plg=list(x="top"))

plot(ras_ext, col=cols, plg=list(x="bottom", horiz=TRUE, inset=c(0, -0.15)))


# # ---------------------------------------------------------
# # ---------------------------------------------------------



# plot(ras_NA, col=cols, main="Vitis.arizonica in North America Context")

# ---------------------------------------------------------
# ---------------------------------------------------------

plot(ras, col=cols, plg=list(x="top"))

if(require(maps)) {
  map("world", xlim=c(-150, -50), ylim=c(15, 70), fill=TRUE, col="lightgrey", border="white")
  map.axes(xlim=c(-160, -40), ylim=c(10,75))
  map("world", xlim=c(-150, -50), ylim=c(15, 70), fill=TRUE, col="lightgrey", border="white")
  map.axes(xlim=c(-160, -40), ylim=c(10,75))
  plot(ras, col=cols, add=TRUE, legend=TRUE,plg=list(x="top", xjust = 3))
}




# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

all_species_list <- c(
  "Vitis.acerifolia",
  "Vitis.aestivalis",
  "Vitis.arizonica",
  "Vitis.californica",
  "Vitis.champinii",
  "Vitis.cinerea",
  "Vitis.girdiana",
  "Vitis.labrusca",
  "Vitis.monticola",
  "Vitis.mustangensis",
  "Vitis.riparia",
  "Vitis.rotundifolia",
  "Vitis.rupestris",
  "Vitis.shuttleworthii",
  "Vitis.vulpina"
)

output_dir <- "./plots"
if (!dir.exists(output_dir)) dir.create(output_dir)

plot_species_map <- function(raster_data, species_name) {
  
  par(mar = c(3, 3, 3, 1))
  
  plot(NULL, xlim=c(-150, -50), ylim=c(15, 70), 
       axes=FALSE, xlab="", ylab="", type="n", xaxs="i", yaxs="i")
  
  if(require(maps)) {
    map("world", fill=TRUE, col="lightgrey", border="white", add=TRUE)
  }
  
  axis(1); axis(2); box()
  
  my_cols <- c("#00000000", "#1f78b4")
  
  plot(raster_data, col=my_cols, add=TRUE, 
       plg=list(x="top", horiz=TRUE, bty="n", xjust=1))
}

for (sp in all_species_list) {
  message(paste("Processing:", sp, "..."))
  
  tryCatch({
    ras_sp <- create_map(path = path, names_species = c(sp))
    
    file_prefix <- file.path(output_dir, sp)
    
    pdf(paste0(file_prefix, ".pdf"), width = 8, height = 6)
    plot_species_map(ras_sp, sp)
    dev.off()
    
    svg(paste0(file_prefix, ".svg"), width = 8, height = 6)
    plot_species_map(ras_sp, sp)
    dev.off()
    
    png(paste0(file_prefix, ".png"), width = 800, height = 600, res = 100)
    plot_species_map(ras_sp, sp)
    dev.off()
    
  }, error = function(e) {
    message(paste("Error processing", sp, ":", e$message))
  })
}

message("All maps generated in ./plots/ folder.")
