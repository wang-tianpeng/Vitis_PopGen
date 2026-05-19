library(terra)
library(tidyverse)
library(maps)

#################### analisis general ##############


# rm(list=ls())
cols <- MetBrewer::met.brewer("Hiroshige")[10:1]
pal <- colorRampPalette(cols)


pdf("./202601_Vitis_centroids_all.pdf", width = 12, height = 8)
files <- list.files("centroids", pattern = "Vitis.*\\.tif$", full.names = TRUE)

my_xlim <- c(-150, -50)
my_ylim <- c(15, 70)

f  <- "centroids/Vitis.rupestris.tif"
for(f in files){
  ras <- rast(f)
  
  crop_ext <- ext(my_xlim[1], my_xlim[2], my_ylim[1], my_ylim[2])
  ras <- crop(ras, crop_ext)
  
  par(mfrow=c(2,2), mar=c(3, 3, 2, 4), oma=c(1, 1, 3, 1))
  
  layer_names <- names(ras)
  layer_names  <- c("Geographic Distribution", 
                    "SDM Niche Distribution", 
                    "Geographic Centroid", 
                    "SDM Niche Centroid")
   main_name <- varnames(ras)
  for(i in 1:2){
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

  mtext(main_name, outer=TRUE, cex=1.5, font=2)
}

dev.off()
par(mfrow=c(1,1))









pdf("./202602_Vitis_centroids_select_all.pdf", width = 15, height = 4)
files <- list.files("centroids", pattern = "Vitis.*\\.tif$", full.names = TRUE)

my_xlim <- c(-150, -50)
my_ylim <- c(15, 70)

f  <- "centroids/Vitis.rupestris.tif"
for(f in files){
  ras <- rast(f)
  
  crop_ext <- ext(my_xlim[1], my_xlim[2], my_ylim[1], my_ylim[2])
  ras <- crop(ras, crop_ext)
  
  par(mfrow=c(1,3), mar=c(3, 3, 2, 4), oma=c(1, 1, 3, 1))
  
  layer_names <- names(ras)
  layer_names  <- c("Binarized SDM prediction", 
                    "Climatic suitability (SDM prediction)", 
                     "Distance to the geographic centroid map" #, 
                    ) # "SDM Niche Centroid"
   main_name <- varnames(ras)
  for(i in c(2,1)){
    plot(NULL, xlim=my_xlim, ylim=my_ylim, 
         axes=FALSE, xlab="", ylab="", type="n", xaxs="i", yaxs="i",
         main = layer_names[i])
    
    if(require(maps)) {
      map("world", fill=TRUE, col="lightgrey", border="white", add=TRUE)
    }
    
    axis(1); axis(2); box()
    
    plot(ras[[i]], col=pal(100), add=TRUE, legend=TRUE,
         plg=list(x="right", bty="n"))
  
  for(i in 3){
    plot(NULL, xlim=my_xlim, ylim=my_ylim, 
         axes=FALSE, xlab="", ylab="", type="n", xaxs="i", yaxs="i",
         main = layer_names[i])
    
    if(require(maps)) {
      map("world", fill=TRUE, col="lightgrey", border="white", add=TRUE)
    }
    
    axis(1); axis(2); box()
    
    plot(ras[[i]], col=rev(pal(100)), add=TRUE, legend=TRUE,
         plg=list(x="right", bty="n"))

  mtext(main_name, outer=TRUE, cex=1.5, font=2)
}

dev.off()
par(mfrow=c(1,1))








########## END FOR ALL SPECIES PLOTTING


## 1.PLUS plot one species map 

sp_name <- "Vitis.arizonica"
f <- list.files("centroids", pattern = paste0(sp_name, "\\.tif$"), full.names = TRUE)

if(length(f) > 0) {
  ras <- rast(f)
  
  my_xlim <- c(-150, -50)
  my_ylim <- c(15, 70)
  crop_ext <- ext(my_xlim[1], my_xlim[2], my_ylim[1], my_ylim[2])
  
  ras_niche <- crop(ras[[2]], crop_ext)
  
  ras_niche[ras_niche < 0.01] <- NA
  
  col_niche <- colorRampPalette(c("#FEE0D2", "#E76254"))(100) # 1e466e
  
  pdf(paste0("./Single_Niche_", sp_name, ".pdf"), width = 8, height = 6)
  par(mar = c(5, 4, 4, 2))
  
  plot(NULL, xlim=my_xlim, ylim=my_ylim, 
       axes=FALSE, xlab="", ylab="", type="n", xaxs="i", yaxs="i",
       main = paste("Niche Distribution:", sp_name),
       cex.main = 0.8)
  
  if(require(maps)) {
    map("world", fill=TRUE, col="lightgrey", border="white", add=TRUE)
  }
  axis(1, labels = FALSE, tick = FALSE); axis(2, labels = FALSE, tick = FALSE); box()
  
  plot(ras_niche, col=col_niche, add=TRUE, legend=F,
       plg=list(x="bottomleft", bty="n", title="Prob"))
  dev.off()
  
  
  ras_cent <- crop(ras[[4]], crop_ext)
  
  col_cent <- colorRampPalette(c("#E76254","#FEE0D2"))(100)
  
  pdf(paste0("./Single_Niche_Centroid_", sp_name, ".pdf"), width = 8, height = 6)
  par(mar = c(5, 4, 4, 2))
  
  plot(NULL, xlim=my_xlim, ylim=my_ylim, 
       axes=FALSE, xlab="", ylab="", type="n", xaxs="i", yaxs="i",
       main = paste("Niche Centroid Distance:", sp_name),
       cex.main = 0.8)
  
  if(require(maps)) {
    map("world", fill=TRUE, col="lightgrey", border="white", add=TRUE)
  }
  axis(1, labels = FALSE, tick = FALSE); axis(2, labels = FALSE, tick = FALSE); box()
  
  plot(ras_cent, col=col_cent, add=TRUE, legend=F,
       plg=list(x="bottomleft", bty="n", title="Dist"))
  dev.off()
  
} else {
  message("File not found for ", sp_name)
}






library(terra)
library(maps)
library(readr)

# ==============================================================================
# ==============================================================================
sp_name <- "Vitis.arizonica"
sample_file <- "email_FollowUp/Sample_Locations.txt"

my_xlim <- c(-150, -50)
my_ylim <- c(15, 70)
crop_ext <- ext(my_xlim[1], my_xlim[2], my_ylim[1], my_ylim[2])

# ==============================================================================
# ==============================================================================

f <- list.files("centroids", pattern = paste0(sp_name, "\\.tif$"), full.names = TRUE)
if(length(f) == 0) stop(paste("Centroid file not found for", sp_name))

ras <- rast(f)
ras_niche <- crop(ras[[4]], crop_ext)

ras_niche[ras_niche < 0.01] <- NA

if(file.exists(sample_file)){
  ind_loc <- read_tsv(sample_file, show_col_types = FALSE)
  if(!all(c("Longitude", "Latitude") %in% names(ind_loc))){
    stop("Sample file must have 'Longitude' and 'Latitude' columns.")
  }
} else {
  warning("Sample locations file not found!")
  ind_loc <- NULL
}

# ==============================================================================
# ==============================================================================

col_niche <- colorRampPalette(c("#E76254","#FEE0D2"))(100)

output_file <- paste0("./Single_Niche_with_Samples_", sp_name, ".pdf")
pdf(output_file, width = 8, height = 6)

par(mar = c(5, 4, 4, 5)) 

plot(NULL, xlim=my_xlim, ylim=my_ylim, 
     axes=FALSE, xlab="", ylab="", type="n", xaxs="i", yaxs="i",
     main = paste(sp_name, "- Niche & Samples"),
     cex.main = 1.2)

if(require(maps)) {
  map("world", fill=TRUE, col="lightgrey", border="white", add=TRUE)
}

axis(1, labels = FALSE, tick = FALSE); axis(2, labels = FALSE, tick = FALSE); box()
  
plot(ras_niche, col=col_niche, add=TRUE, legend=F,
     plg=list(x="right", inset=c(-0.1, 0), bty="n", title="Prob"))

if(!is.null(ind_loc)) {
  points(ind_loc$Longitude, ind_loc$Latitude, 
         pch = 20,
         cex = 0.8,
         col = "
}

dev.off()

message(paste("Map saved to:", output_file))










### PLOT 2 SPECIES MAPS COMPARISON ###


### PLOT 2 SPECIES MAPS COMPARISON ###

sp1_name <- "Vitis.riparia" # "Vitis.arizonica"
sp2_name <- "Vitis.aestivalis" # "Vitis.girdiana"


sp1_name <- "Vitis.arizonica" # "Vitis.arizonica"
sp2_name <- "Vitis.californica" # "Vitis.girdiana"

sp1_name <- "Vitis.mustangensis" # "Vitis.arizonica"
sp2_name <- "Vitis.acerifolia" # "Vitis.girdiana"
sample_location_doan <- "email_FollowUp/Sample_Locations_DOAN.txt"
ind_loc_doan <- read_tsv(sample_location_doan, show_col_types = FALSE)
ind_loc <- ind_loc_doan

file1 <- list.files("centroids", pattern = paste0(sp1_name, "\\.tif$"), full.names = TRUE)
file2 <- list.files("centroids", pattern = paste0(sp2_name, "\\.tif$"), full.names = TRUE)




ras1 <- rast(file1)[[4]]
ras2 <- rast(file2)[[4]]

# ras1[ras1 < threshold] <- NA
# ras2[ras2 < threshold] <- NA

ras_overlap <- (ras1 + ras2) / 2

my_xlim <- c(-125, -75) # for doaniana
my_ylim <- c(25, 50)

my_xlim <- c(-150, -50)
my_ylim <- c(15, 70)
crop_ext <- ext(my_xlim[1], my_xlim[2], my_ylim[1], my_ylim[2])
ras1 <- crop(ras1, crop_ext)
ras2 <- crop(ras2, crop_ext)
ras_overlap <- crop(ras_overlap, crop_ext)


  
col1_pal <- colorRampPalette(c("#E76254","#FEE0D2"))(100) 

col2_pal <- colorRampPalette(c("#528fad", "#DEEBF7"))(100)
col_overlap_pal <- colorRampPalette(c("#FFE6B7","#FFD06F" ))(10)

pdf(paste0("./Overlap_", sp1_name, "_vs_", sp2_name, ".pdf"), width = 8, height = 6)

par(mar = c(5, 4, 4, 8)) 
plot(NULL, xlim=my_xlim, ylim=my_ylim, 
     axes=FALSE, xlab="", ylab="", type="n", xaxs="i", yaxs="i",
     main = paste("Niche Overlap:", sp1_name, "vs", sp2_name))

map("world", fill=TRUE, col="grey90", border="white", add=TRUE)
axis(1, labels = FALSE, tick = FALSE); 
axis(2, labels = FALSE, tick = FALSE); 
box()

plot(ras1, col="#3FA7D6", add=TRUE,  legend=FALSE) #col2_pal

plot(ras2, col="#DF8277", add=TRUE, legend=FALSE)

# plot(ras_overlap, col=col_overlap_pal, add=TRUE, legend=FALSE)

if(!is.null(ind_loc)) {
  points(ind_loc_doan$Longitude, ind_loc_doan$Latitude, 
         pch = 20,
         cex = 1,
         col = "
}

# plot(ras1, col=col1_pal, legend.only=TRUE, 
#      plg=list(x="topright", inset=c(-0.22, 0), title=sp1_name, title.cex=0.8))

# plot(ras_overlap, col=col_overlap_pal, legend.only=TRUE, 
#      plg=list(x="right", inset=c(-0.22, 0), title="Overlap", title.cex=0.8))

# plot(ras2, col=col2_pal, legend.only=TRUE, 
#      plg=list(x="bottomright", inset=c(-0.22, 0), title=sp2_name, title.cex=0.8))

dev.off()



### PLOT 3: NICHE CENTROID DISTANCE COMPARISON ###

ras1_cent <- rast(file1)[[4]]
ras2_cent <- rast(file2)[[4]]

ras1_cent <- crop(ras1_cent, crop_ext)
ras2_cent <- crop(ras2_cent, crop_ext)

col1_cent_pal <- colorRampPalette(c("#1e466e", "#528fad"  ))(10) # Greens
col2_cent_pal <- colorRampPalette(c("#E76254", "#f7aa58"))(10) # RdPu

pdf(paste0("./Niche_Centroid_", sp1_name, "_vs_", sp2_name, ".pdf"), width = 8, height = 6)

par(mar = c(5, 4, 4, 8)) 
plot(NULL, xlim=my_xlim, ylim=my_ylim, 
     axes=FALSE, xlab="", ylab="", type="n", xaxs="i", yaxs="i",
     main = paste("Niche Centroid Dist:", sp1_name, "vs", sp2_name),
     cex.main = 0.8,
     cex.axis = 0.7)

map("world", fill=TRUE, col="lightgrey", border="white", add=TRUE)
axis(1, labels = FALSE, tick = FALSE); 
axis(2, labels = FALSE, tick = FALSE); 
box()

plot(ras1_cent, col=col1_cent_pal, add=TRUE, alpha=0.9, legend=FALSE)

plot(ras2_cent, col=col2_cent_pal, add=TRUE, alpha=0.9, legend=FALSE)

# plot(ras1_cent, col=col1_cent_pal, legend.only=TRUE, 
#      plg=list(x="topright", inset=c(-0.22, 0), title=sp1_name, title.cex=0.8))

# plot(ras2_cent, col=col2_cent_pal, legend.only=TRUE, 
#      plg=list(x="bottomright", inset=c(-0.22, 0), title=sp2_name, title.cex=0.8))

dev.off()




library(terra)
library(maps)
library(MetBrewer)

# ==============================================================================
# ==============================================================================
sp1_name <- "Vitis.arizonica"
sp2_name <- "Vitis.californica"
my_xlim <- c(-150, -50)
my_ylim <- c(15, 70)
crop_ext <- ext(my_xlim[1], my_xlim[2], my_ylim[1], my_ylim[2])

CEX_MAIN <- 1.0
CEX_LEGEND_TITLE <- 0.9

f1 <- list.files("centroids", pattern = paste0(sp1_name, "\\.tif$"), full.names = TRUE)
f2 <- list.files("centroids", pattern = paste0(sp2_name, "\\.tif$"), full.names = TRUE)

if(length(f1) == 0 || length(f2) == 0) stop("Files not found!")

ras1_all <- rast(f1)
ras2_all <- rast(f2)

# ==============================================================================
# ==============================================================================

ras1_niche <- crop(ras1_all[[2]], crop_ext)
# ras1_niche[ras1_niche < 0.01] <- NA
col_niche_single <- colorRampPalette(c("#FEE0D2", "#E76254"))(100)

ras1_cent <- crop(ras1_all[[4]], crop_ext)
col_cent_single <- colorRampPalette(c("#E76254", "#f89f96"))(100)

ras2_niche <- crop(ras2_all[[2]], crop_ext)
# ras2_niche[ras2_niche < 0.01] <- NA
ras_overlap <- (ras1_niche + ras2_niche) / 2

col_p3_sp1 <- colorRampPalette(c( "#FEE0D2", "#E76254"))(10) # Blueish
col_p3_sp2 <- colorRampPalette(c("#DEEBF7", "#1e466e"))(10) # Reddish
col_p3_over <- colorRampPalette(c("#FFE6B7","#FFD06F"))(10) # Yellowish

ras2_cent <- crop(ras2_all[[4]], crop_ext)
col_p4_sp1 <- colorRampPalette(c("#E76254", "#f7aa58" ))(10) # Blueish
col_p4_sp2 <- colorRampPalette(c("#1e466e", "#528fad"))(10) # Reddish

# ==============================================================================
# ==============================================================================


# ==============================================================================
# ==============================================================================
draw_combined_maps <- function() {
  par(mfrow=c(2,2), oma=c(1,1,2,1), mar=c(2, 2, 3, 5))
#   par(mfrow=c(2,2), oma=c(0.5, 0.5, 1, 0.5), mar=c(1, 1, 2, 8))

  # ------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------
  plot(NULL, xlim=my_xlim, ylim=my_ylim, axes=FALSE, xlab="", ylab="", type="n", xaxs="i", yaxs="i",
       main = paste(sp1_name, "Niche Prob"), cex.main = CEX_MAIN)
  if(require(maps)) map("world", fill=TRUE, col="lightgrey", border="white", add=TRUE)
  box()
  
  plot(ras1_niche, col=col_niche_single, add=TRUE, legend=TRUE,
       plg=list(x="right", inset=c(-0.15, 0), bty="n", title="Prob", title.cex=CEX_LEGEND_TITLE),
       axis.args=list(cex.axis=CEX_AXIS))
  
  # ------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------
  plot(NULL, xlim=my_xlim, ylim=my_ylim, axes=FALSE, xlab="", ylab="", type="n", xaxs="i", yaxs="i",
       main = paste(sp1_name, "Centroid Dist"), cex.main = CEX_MAIN)
  if(require(maps)) map("world", fill=TRUE, col="lightgrey", border="white", add=TRUE)
  box()
  
  plot(ras1_cent, col=col_cent_single, add=TRUE, legend=TRUE,
       plg=list(x="right", inset=c(-0.15, 0), bty="n", title="Dist", title.cex=CEX_LEGEND_TITLE),
       axis.args=list(cex.axis=CEX_AXIS))
  
  # ------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------
  plot(NULL, xlim=my_xlim, ylim=my_ylim, axes=FALSE, xlab="", ylab="", type="n", xaxs="i", yaxs="i",
       main = "Niche Overlap Comparison", cex.main = CEX_MAIN)
  if(require(maps)) map("world", fill=TRUE, col="lightgrey", border="white", add=TRUE)
  box()
  
  plot(ras1_niche, col=col_p3_sp1, add=TRUE, alpha=0.9, legend=FALSE)
  plot(ras2_niche, col=col_p3_sp2, add=TRUE, alpha=0.9, legend=FALSE)
  plot(ras_overlap, col=col_p3_over, add=TRUE, alpha=0.9, legend=FALSE)
  
#   plot(ras1_niche, col=col_p3_sp1, legend.only=TRUE, 
#        plg=list(x="topright", inset=c(-0.25, 0), title=sp1_name, title.cex=CEX_LEGEND_TITLE),
#        axis.args=list(cex.axis=CEX_AXIS))
#   plot(ras_overlap, col=col_p3_over, legend.only=TRUE, 
#        plg=list(x="right", inset=c(-0.25, 0), title="Overlap", title.cex=CEX_LEGEND_TITLE),
#        axis.args=list(cex.axis=CEX_AXIS))
#   plot(ras2_niche, col=col_p3_sp2, legend.only=TRUE, 
#        plg=list(x="bottomright", inset=c(-0.25, 0), title=sp2_name, title.cex=CEX_LEGEND_TITLE),
#        axis.args=list(cex.axis=CEX_AXIS))
  
  # ------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------
  plot(NULL, xlim=my_xlim, ylim=my_ylim, axes=FALSE, xlab="", ylab="", type="n", xaxs="i", yaxs="i",
       main = "Centroid Distance Comparison", cex.main = CEX_MAIN)
  if(require(maps)) map("world", fill=TRUE, col="lightgrey", border="white", add=TRUE)
  box()
  
  plot(ras1_cent, col=col_p4_sp1, add=TRUE, alpha=0.9, legend=FALSE)
  plot(ras2_cent, col=col_p4_sp2, add=TRUE, alpha=0.9, legend=FALSE)
  
#   plot(ras1_cent, col=col_p4_sp1, legend.only=TRUE, 
#        plg=list(x="topright", inset=c(-0.25, 0), title=sp1_name, title.cex=CEX_LEGEND_TITLE),
#        axis.args=list(cex.axis=CEX_AXIS))
#   plot(ras2_cent, col=col_p4_sp2, legend.only=TRUE, 
#        plg=list(x="bottomright", inset=c(-0.25, 0), title=sp2_name, title.cex=CEX_LEGEND_TITLE),
#        axis.args=list(cex.axis=CEX_AXIS))
}

# ==============================================================================
# ==============================================================================

pdf("./Combined_Analysis_Maps.pdf", width = 14, height = 10)
draw_combined_maps()
dev.off()

svg("./Combined_Analysis_Maps.svg", width = 14, height = 10)
draw_combined_maps()
dev.off()

png("./Combined_Analysis_Maps.png", width = 2800, height = 2000, res = 200)
draw_combined_maps()
dev.off()

message("All combined maps (PDF, SVG, PNG) generated successfully.")
