
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtree")


library(Momocs)
library(here)
library(tidyverse)
library(ggtree)



# STEP 1: Import Data

lf <- list.files("Outlines", full.names = TRUE)

coo <- import_jpg(lf, auto.notcentered = TRUE, 
                  fun.notcentered = NULL, 
                  threshold = 0.5)

database <- read.csv("HandaxeMaster.csv")

database$Site <- as.factor(database$Site)
is.factor(database$Site)

data <- Out(coo, fac = database)


mosaic(data, ~Site, asp = 1, legend = TRUE)




# STEP 2: Normalise the data

shapenorm <- coo_center(data)
shapenorm <- coo_scale(shapenorm)
shapenorm <- coo_close(shapenorm)

stack(shapenorm, title = "Stack: Normalised Outlines")


shapenorm2 <- shapenorm %>% 
  coo_slidedirection("right") %>% 
  coo_untiltx()   ### insert break where drawn outlines begin (visual only)

stack(shapenorm2, title = "Stack: Normalised Outlines")




# STEP 3: Elliptic Fourier Transformation

calibrate_harmonicpower_efourier(shapenorm2, nb.h = 20, plot = FALSE)

calibrate_reconstructions_efourier(shapenorm2, range = 1:20)

calibrate_deviations_efourier(shapenorm2, id = 4)

efashape <- efourier(shapenorm2, nb.h = 17, smooth.it = 0, norm = FALSE)




# STEP 4: Principal Component Analysis

pcashape <- PCA(efashape)

scree(pcashape)

scree_plot(pcashape, nax = 1:6)

PCcontrib(pcashape, nax = 1:6)




#STEP 5: Plot analysis

boxplot(pcashape, ~Site, nax = 1:5) + theme_minimal()

plot_PCA(pcashape,axes = c(1, 2), ~Site,
         morphospace_position = "circle",
         zoom = 0.8,
         chull = FALSE) %>%
  layer_points(cex = 1)

plot_PCA(pcashape,axes = c(1, 2), ~Site,
         morphospace_position = "range_axes",
         zoom = 0.8,
         chull = FALSE,
         chullfilled = TRUE) %>%
  layer_points(cex = 1)

plot_PCA(pcashape, axes = c(1, 2), ~Site,
         morphospace_position = "full_axes",
         zoom = 0.8,
         chull = FALSE,
         chullfilled = TRUE) %>%
  layer_points(cex = 1)

plot_PCA(pcashape, axes = c(1, 2), ~Site,
         morphospace_position = "circle",
         zoom = 0.8, 
         labelpoints = TRUE,
         chull = FALSE) %>%
  layer_points(cex = 1)




#STEP 6: Discriminant Analysis

dashapefc <- LDA(efashape, ~Site)
dashape95 <- LDA(pcashape, ~Site, retain = 0.95)

dashapefc$CV.correct
dashapefc$CV.ce
dashape95$CV.correct
dashape95$CV.ce

classification_metrics(dashapefc)
classification_metrics(dashape95)

plot_LDA(dashapefc,               
         axes = c(1,2),
         zoom = 0.9,
         chull = FALSE) %>% 
  layer_points(cex = 1) %>% 
  layer_morphospace_LDA(position = "circle")

plot_LDA(dashapefc,
         axes = c(1,2),
         zoom = 0.9,
         chull = FALSE,
         chullfilled = TRUE) %>% 
  layer_points(cex = 1) %>% 
  layer_morphospace_LDA(position = "range_axes")

plot_LDA(dashapefc,
         axes = c(1,2),
         zoom = 1.2,
         chull = FALSE,
         chullfilled = TRUE) %>% 
  layer_points(cex = 1) %>% 
  layer_morphospace_LDA(position = "full_axes")

plot_LDA(dashapefc, 
         axes = c(1,2),
         zoom = 1,
         chull = FALSE) %>% 
  layer_points(cex = 1) %>% 
  layer_morphospace_LDA(position = "circle") %>%
  layer_ellipses()


efashape %>% MANOVA(~Site)
pcashape %>% MANOVA(~Site, retain = 0.95)

pcashape %>% MANOVA_PW(~Site, retain = 0.95)




# STEP 7: Hierarchical Cluster Analysis

pcascores <- as.data.frame(pcashape$x) %>% 
  rownames_to_column() %>% tibble()

database.expanded <- as.tibble(cbind(database, pcascores))

tree <- hclust(dist(database.expanded[,11:34]))

reference <- 
  database.expanded %>%
  select(Site, rowname) %>%
  as.tibble() %>%
  rename(artefact_id = rowname) %>%
  rownames_to_column() %>%
  rename(id = rowname) %>%
  mutate(Site = as.character(Site)) %>%
  mutate(
    Site2 = case_when(
      Site == "Pontnewydd" ~ "Pontnewydd",
      Site == "Lynford" ~ "Lynford",
      Site == "Walcott (Bacton)" ~ "Offshore",
      Site == "Area 240" ~ "Offshore",
      TRUE ~ Site)) %>%
  mutate(Site = as.factor(Site)) %>%
  mutate(Site2 = as.factor(Site2))

ggtree(tree, layout = "circular") %<+% reference +
  geom_tippoint(aes(colour = Site), size = 3) +
  geom_tiplab2(aes(label = artefact_id, colour = Site, align = TRUE), offset = 0.005) +
  theme(plot.margin = margin(40,40,40,40),
        legend.position = "bottom") +
  scale_color_brewer("Site 2", palette = "Set2")

ggtree(tree, layout = "circular") %<+% reference +
  geom_tippoint(aes(colour = Site2), size = 3) +
  geom_tiplab2(aes(label = artefact_id, colour = Site2, align = TRUE), offset = 0.005) +
  theme(plot.margin = margin(40,40,40,40),
        legend.position = "bottom") +
  scale_color_brewer("Site 2", palette = "Set2")




ggsave("Images/HCA.tiff", plot = last_plot(), dpi = 450, height = 350, width = 300, units = "mm")


