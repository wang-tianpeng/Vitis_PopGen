
### PLOT 2 SPECIES MAPS COMPARISON ###

sp1_name <- "Vitis.riparia" # "Vitis.arizonica"
sp2_name <- "Vitis.aestivalis" # "Vitis.girdiana"


sp1_name <- "Vitis.arizonica" # "Vitis.arizonica"
sp2_name <- "Vitis.californica" # "Vitis.girdiana"

sp1_name <- "Vitis.mustangensis" # "Vitis.arizonica"
sp2_name <- "Vitis.acerifolia" # "Vitis.girdiana"
sp3_name <- "Vitis.doaniana"
sample_location_doan <- "email_FollowUp/Sample_Locations_DOAN.txt"
ind_loc_doan <- read_tsv(sample_location_doan, show_col_types = FALSE)
ind_loc <- ind_loc_doan

file1 <- list.files("centroids", pattern = paste0(sp1_name, "\\.tif$"), full.names = TRUE)
file2 <- list.files("centroids", pattern = paste0(sp2_name, "\\.tif$"), full.names = TRUE)
file3 <- list.files("centroids", pattern = paste0(sp3_name, "\\.tif$"), full.names = TRUE)



ras1 <- rast(file1)[[4]]
ras2 <- rast(file2)[[4]]
ras3 <- rast(file3)[[4]]

# ras1[ras1 < threshold] <- NA
# ras2[ras2 < threshold] <- NA

ras_overlap <- (ras1 + ras2) / 2

my_xlim <- c(-125, -75) # for doaniana
my_ylim <- c(25, 50)

crop_ext <- ext(my_xlim[1], my_xlim[2], my_ylim[1], my_ylim[2])
ras1 <- crop(ras1, crop_ext)
ras2 <- crop(ras2, crop_ext)
ras3 <- crop(ras3, crop_ext)


  
col1_pal <- colorRampPalette(c("#E76254","#FEE0D2"))(100) 

col2_pal <- colorRampPalette(c("#528fad", "#DEEBF7"))(100)
col3_pal <- colorRampPalette(c("#FFD06F","#FFE6B7" ))(10)


pdf(paste0("./Overlap_", sp1_name, "_vs_", sp2_name, sp3_name, ".pdf"), width = 8, height = 6)

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

plot(ras3, col="#FFD06F", add=TRUE, legend=FALSE)


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
