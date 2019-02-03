library(tidyverse)
library(vegan)
library(sf)
library(devtools)
# for 3D plots install this:
# devtools::install_github("AckerDWM/gg3D")
# library(gg3D)
library(rgl)
library(scatterplot3d)



# DATA =========================================================================

# read checklist data
spec_dist <- read_tsv("data/Dist_Legitimate.txt")
# read OGU shapefile
ogu <- st_read("data/mapdata/OGU.shp")

# get number of species per OGU
spec_count <- spec_dist %>%
  group_by(AreaCode) %>%
  summarise(richness = n())
# exclude OGUs with less than 100 species
spec <- spec_dist %>%
  left_join(spec_count, by = "AreaCode") %>%
  filter(richness > 100)


# ORDINATION ===================================================================

# create dissimilarity matrix
xtab <- table(
  spec[["AreaCode"]],
  spec[["Genus"]]) %>%
  unclass()

# use custom distance metric
distfun <- function(comm, method) {
  dist <- betadiver(x = comm, method = "sim")
  return(dist)
}

# run multidimensional scaling
mds_fitted <- metaMDS(
  comm = xtab,
  k = 3,
  distfun = distfun,
  distance = "mixed",
  zerodist = "add",
  noshare = 0.1,
  trymax = 100,
  autotransform = T,
  wascores = T,
  expand = T,
  pc = T)

# values for each species
mds_fitted[["species"]]


# PLOTTING =====================================================================

# extract raw coordinates
mds_ogu <- mds_fitted %>%
  pluck("points") %>%
  as_tibble(rownames = "AreaCode")

# Helper function to rescale mds dimensions to match colour range
rescale_mds <- function(data) {
  # extract MDS valus
  mds <- select(data, MDS1, MDS2, MDS3)
  # rescale MDS values to [0 - 255]
  mds <- ((mds - min(mds)) / max(mds - min(mds))) * 255
  # recombine with labels
  data[["MDS1"]] <- mds[[1]]
  data[["MDS2"]] <- mds[[2]]
  data[["MDS3"]] <- mds[[3]]
  return(data)
}

# create colour scale from MDS values
mds_ogu_rescaled <- rescale_mds(mds_ogu)
rgb_colours <- rgb(
  red = mds_ogu_rescaled[["MDS1"]],
  green = mds_ogu_rescaled[["MDS2"]],
  blue = mds_ogu_rescaled[["MDS3"]],
  maxColorValue = 255)

# make interactive 3D plot
plot3d(
  x = mds_ogu_rescaled[["MDS1"]],
  y = mds_ogu_rescaled[["MDS2"]],
  z = mds_ogu_rescaled[["MDS3"]],
  xlab = "x",
  ylab = "y",
  zlab = "z",
  size = 10,
  col = rgb_colours)
# add text labels
text3d(
  x = mds_ogu_rescaled[["MDS1"]],
  y = mds_ogu_rescaled[["MDS2"]],
  z = mds_ogu_rescaled[["MDS3"]],
  text = mds_ogu_rescaled[["AreaCode"]], adj = 1) 

# static 3d plot
s3d <- scatterplot3d(
  x = mds_ogu_rescaled[["MDS1"]],
  y = mds_ogu_rescaled[["MDS2"]],
  z = mds_ogu_rescaled[["MDS3"]],
  color = rgb_colours,
  xlab = "", ylab = "", zlab = "",
  xlim = c(0,255), ylim = c(0,255), zlim = c(0,255),
  pch = 16, cex.symbols = 1, cex.lab = 0.6, cex.axis = 0.6,
  lty.hplot = 3, type = "h",
  box = FALSE, mar = c(1, 1, 1, 2))
# add border
s3d <- s3d$points3d(
  x = mds_ogu_rescaled[["MDS1"]],
  y = mds_ogu_rescaled[["MDS2"]],
  z = mds_ogu_rescaled[["MDS3"]],
  pch = 21, cex=1, bg = rgb_colours)
# text(9.5, 1, "NMDS axis 2", cex = 0.9, srt = 40)
# text(s3d$xyz.convert(-20,0,365), labels="(a)", cex=2)
