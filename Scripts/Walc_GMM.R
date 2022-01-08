
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtree")


library(Momocs)
library(here)
library(tidyverse)
library(ggtree)



# STEP 1: Import Data

lf <- list.files("WalcOut", full.names = TRUE)

coo <- import_jpg(lf, auto.notcentered = TRUE, 
                  fun.notcentered = NULL, 
                  threshold = 0.5)

database <- read.csv("Walcott.csv")

database$Site <- as.factor(database$Site)
is.factor(database$Site)

data <- Out(coo, fac = database)


mosaic(data, asp = 1, legend = TRUE)




# STEP 2: Normalise the data

shapenorm <- coo_center(data)
shapenorm <- coo_scale(shapenorm)
shapenorm <- coo_close(shapenorm)

stack(shapenorm, title = "Stack: Normalised Walcott Outlines")


shapenorm2 <- shapenorm %>% 
  coo_slidedirection("right") %>% 
  coo_untiltx()

stack(shapenorm2, title = "Stack: Normalised Walcott Outlines")




# STEP 3: Elliptic Fourier Transformation

calibrate_harmonicpower_efourier(shapenorm2, nb.h = 20, plot = FALSE)

calibrate_reconstructions_efourier(shapenorm2, range = 1:20)

calibrate_deviations_efourier(shapenorm2, id = 4)

efashape <- efourier(shapenorm2, nb.h = 17, smooth.it = 0, norm = FALSE)




# STEP 4: Principal Component Analysis

pcashape <- PCA(efashape)

scree(pcashape)

scree_plot(pcashape, nax = 1:6)

PCcontrib(pcashape, nax = 1:3)




#STEP 5: Plot analysis

boxplot(pcashape, ~Site, nax = 1:5) + theme_minimal()

plot_PCA(pcashape,axes = c(1, 2),
         morphospace_position = "range_axes",
         zoom = 1.5,
         chull = FALSE) %>%
  layer_points(cex = 1)

plot_PCA(pcashape, axes = c(1, 2),
         morphospace_position = "full_axes",
         zoom = 1.5,
         chull = FALSE) %>%
  layer_points(cex = 1) 

plot_PCA(pcashape, axes = c(1, 2),
         morphospace_position = "circle",
         zoom = 1.5, 
         chull = FALSE,
         labelpoints = TRUE) %>%
  layer_points(cex = 1)




# STEP 6: HCA & Cladistics

pcascores <- as.data.frame(pcashape$x) %>% 
  rownames_to_column() %>% tibble()

database.expanded <- as.tibble(cbind(database, pcascores))

tree <- hclust(dist(database.expanded[,12:38]))

reference <- 
  database.expanded %>%
  select(rowname) %>%
  as.tibble() %>%
  rename(artefact_id = rowname) %>%
  rownames_to_column() %>%
  rename(id = rowname)

ggtree(tree) %<+% reference +
  geom_tippoint(size = 3) +
  geom_tiplab(aes(label = artefact_id), size = 8, offset = 0.01, align = TRUE) +
  theme(plot.margin = margin(40,40,40,40))



# Cladistics
ggtree(tree, branch.length = "none") %<+% reference +
  geom_tippoint(size = 3) +
  geom_tiplab(aes(label = artefact_id), size = 8, offset = 0.01, align = TRUE) +
  theme(plot.margin = margin(40,40,40,40))




ggsave("Images/HCA.tiff", plot = last_plot(), dpi = 450, height = 350, width = 300, units = "mm")


