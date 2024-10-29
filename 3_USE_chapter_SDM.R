setwd("/home/daniele/Documents/working_files/USE_bookChapter/")

#R libraries
#handle spatial data
library(sf)
library(terra)
#handle data
library(dplyr)
# #plot data
library(ggplot2)
library(viridis)
# collienarity
library(caret) 
library(usdm)
#modelling tools
library(sdm) 
library(USE)
# some metrics  
library(Metrics)
library(ecospat)

# Get density of points in 2 dimensions.
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

outdir <- getwd()

#---- 1. Load climatic data ----
#download country area and bioclimatic variables
ita <- geodata::gadm(country="ITA", 
                     level=0, 
                     path = outdir)
bio <- geodata::worldclim_global(var = 'bio', 
                                     res = 5, 
                                     path = outdir)
bio <- terra::crop(bio, ita, mask = TRUE )
names(bio) <- paste0("bio", seq_along(names(bio)))

##---- 1.1 Find optimal resolution for sampling Italian environmental (here climatic) space ----
PCA_ita <- USE::rastPCA(env.rast = bio, nPC = 2, stand = TRUE)
PCstack <- c(PCA_ita$PCs$PC1, PCA_ita$PCs$PC2)

#Get dataframe from PCstack
PCstack.df <- as.data.frame(PCstack, xy = T, na.rm = T)

#make the PCstack.df spatial
PCstack.sp <- st_as_sf(PCstack.df, coords = c("PC1", "PC2"))

#get optimal resolution of the grid for sampling Italian environmental space
Optres_ita <- USE::optimRes(sdf = PCstack.sp, grid.res = seq(1, 15), cr = 4, showOpt = T)


#---- 2. Load occurrences ----
occ <- readRDS("outputs/gbifClean.RDS") %>%
  dplyr::select( decimalLongitude, decimalLatitude) %>% 
  dplyr::rename(lon = "decimalLongitude", lat = "decimalLatitude") %>%  
  sf::st_as_sf(coords=c("lon", "lat"), crs=4326) %>% 
  mutate(PA = 1)

##---- 2.1 environmental thinning ----
#extract pc scores
occ <- occ %>% 
  bind_cols(terra::extract(PCstack, occ,  xy=TRUE, ID=FALSE)) %>% 
  relocate( geometry, .after = last_col())

occ.env <- occ %>% 
  sf::st_drop_geometry() %>% 
  sf::st_as_sf(coords=c("PC1", "PC2"))

# subset presences using uniformSampling
# 10 presences will be (whenever possible) sampled from each cell of the sampling grid
#same thing will be done to sample testing presences
set.seed(17534)
occ.pres.tr.ts <- USE::uniformSampling(sdf = occ.env, 
                                       grid.res = Optres_ita$Opt_res,
                                       n.tr = 10, 
                                       sub.ts = TRUE, n.ts = 10, 
                                       plot_proc = FALSE)




occ.thin.env <- occ.pres.tr.ts$obs.tr

#make training and testing (presence) datasets spatial 
occ.pres.tr.ts <- lapply(occ.pres.tr.ts, function(x) {
  spdt <- x
  spdt$geometry <- NULL
  spdt.sp <- st_as_sf(spdt, coords = c("x", "y"))
  st_crs(spdt.sp) <- 4326
  return(spdt.sp)
})

lapply(occ.pres.tr.ts, nrow)

#---- 3. Get pseudo-absences within the environmental space with paSampling ----
#get pseudo-absences using paSampling
#here we use all the species presences we obtained from the cleaned gbif dataset (n = 673) to inform the kernel-based filter, so that we can use all information available about locations suitable for the species
#this will allow us safely excluding as many pseudo-absences located in areas of the environmental space where conditions may actually be suitable for Acer (presences) as possible

#however, as we want prevalence to be = 1 in the training set, we set the prev arg as if we were using occ.pres.tr.ts$obs.tr
#(so if we had n= 450 presences)
#basically we 'inflate' prevalence so that we can finally get a number of pseudo-absences comparable with occ.pres.tr.ts$obs.tr

#our prevalence is
myPrev <- nrow(occ)/(nrow(occ.pres.tr.ts$obs.tr)) 
myPrev

set.seed(17893)
occ.abs <- USE::paSampling(env.rast = bio,
                              pres = occ, n.tr = 10,
                              grid.res = Optres_ita$Opt_res, 
                              prev = myPrev,
                              sub.ts = TRUE, 
                              n.ts = 10, 
                              plot_proc = FALSE)


# dispaly thinned occurrences and pseudo-absences in the environmental space
env.pa <- occ.thin.env%>% 
  select(PA) %>% 
  bind_rows(occ.abs$obs.tr %>% 
              select(PA)) 

env.pa <- data.frame(PA=env.pa$PA,
                    st_coordinates(env.pa)) %>% 
  mutate(PA=as.factor(PA), 
         PA=factor(PA, levels=c("1", "0"))) %>% 
  rename(PC1=X, PC2=Y)

PCstack.df$density <- get_density(PCstack.df$PC1, PCstack.df$PC2, n = 100)
p <- ggplot() +
    geom_point(data = PCstack.df, aes(PC1, PC2, color = density), alpha = 0.7) +
    scale_color_viridis(name = "Density of PC-scores", option = "rocket") +
    geom_point(data = env.pa, aes(PC1, PC2, fill = PA), shape=21, size=3) +
    scale_fill_manual(name = "Thinned presence/Pseudo-absence", values = c("1" = "#1E88E5", "0" = "lightgrey")) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      text = element_text(size = 12),
      legend.text = element_text(size = 12)
    )
outname <- paste("figures/", "thinnedPresencePseudoAbs_sdm_", Sys.Date(),".png", sep="")
ggsave(p, filename = outname, width = 18, height = 8, device='png', dpi=320)
  
#make training and testing (absence) datasets spatial 
abs.tr.ts.sp <- lapply(occ.abs, function(x) {
  spdt <- x
  spdt$geometry <- NULL
  spdt.sp <- st_as_sf(spdt, coords = c("x", "y"))
  st_crs(spdt.sp) <- 4326
  return(spdt.sp)
})






# plot in environmental PC space
plot(dt2[, c('PC1', 'PC2')], cex = 0.5, col = gray(0.8))
plot(st_geometry(myGrid.psAbs), add = TRUE) # pseudo-absences
points(pcPres[, c('PC1', 'PC2')], col = 'purple') # presences

#extract climatic data and join training and testing datasets
species_data <- Map( x = occ.pres.tr.ts, y = abs.tr.ts.sp, function(x, y) {
  x <- data.frame(PA = 1, st_coordinates(x))
  y <- data.frame(PA = 0, st_coordinates(y))
  colnames(x)[-1] <- colnames(y)[-1] <- c("long", "lat")
  df <- rbind(x, y)
  df <- data.frame(df, terra::extract(bio, df[c(2, 3)]))
  return(df)
})

#set names of the list including training and testing datasets
names(species_data) <- c("Tr_data", "Ts_data")

#n° of presences and absences
table(species_data$Tr_data$PA)
table(species_data$Ts_data$PA)

#make prevalence in the testing set = 1 (i.e. reduce number of testing presences)
Subs.ts.pres <- species_data$Ts_data
Subs.ts.pres <- Subs.ts.pres[Subs.ts.pres$PA == 1, ]
Subs.ts.pres <- Subs.ts.pres[sample(x = nrow(Subs.ts.pres), size = sum(species_data$Ts_data$PA == 0), replace = T), ]

#substitute testing dataset in FagusEU_data
species_data$Ts_data <- rbind(Subs.ts.pres, species_data$Ts_data[species_data$Ts_data$PA == 0, ])

#check
table(species_data$Ts_data$PA)

#check NAs
sapply(species_data, anyNA) #F F
# saveRDS(species_data, "cleanedOccDf.RDS")

#----4. Modelling ----
##---- 4.1. Step 1 Avoid multicollinearity ----
# Variance inflation factor measures how much the behavior (variance) of an independent variable is influenced, or inflated, by its interaction/correlation with the other independent variables.
# VIF selection let's keep 3 as a threshold (Disclaimer: there is a debate around this)
Vars_to_remove<-usdm::vifstep(species_data$Tr_data[,5:ncol(species_data$Tr_data)], th=3, keep = c("bio1", "bio12"))
Vars_to_remove

#get rid of correlated variables from training and testing datasets
species_data <- lapply(species_data, function(x) {
  x <- x[!colnames(x) %in% Vars_to_remove@excluded]
  return(x)
})

#remove also othe not necessary variables
species_data.f <- lapply(species_data, function(x) {
  x <- x[!colnames(x) %in% c("ID")]
  return(x)
})

#check distribution of values of predictors for presence and absence data
ggarrange(plotlist = lapply(colnames(species_data.f$Tr_data[,-c(2:3)]), function(nm) {
  ggplot(species_data.f$Tr_data[,-c(2:3)], aes_string(x = nm)) +
    geom_density(aes(fill = as.factor(PA)), alpha = .3) +
    scale_fill_viridis_d(name = paste("Pres (1)", "Abs (0)", sep = "\n")) +
    theme_classic()
}),
nrow = 2, ncol = 3)

##---- 4.2. modelling with the sdm package ----
percTesting <- 30
Nreps <- 5
myCores <- 5

# extract enviromental data
d <- sdmData(formula=PA~.  + coords(long+lat), 
             train=species_data.f$Tr_data, 
             test=species_data.f$Ts_data)
class(d)

#model
m1 <- sdm::sdm(PA~.,
          data=d,
          methods=c('glm', 'gbm'),
          replication='sub',
          test.percent=percTesting, 
          n=Nreps, 
          parallelSettings=list("parallel", ncore=myCores))

# Inspect the model
m1

#partial dependency plots
#check id
m1@models$PA$brt

#glm
rcurve(m1, id=1:5)
getVarImp(m1, id=1:5)
#brt
rcurve(m1, id=6:10)
getVarImp(m1, id=6:10)

##---- 4.3. Model performance----
# We can generate the roc curve and compare the results for all models:
roc(m1)

myComp.sdm <- setNames(data.frame(matrix(ncol = 8, nrow = 0)), c("Model", "AUC",  "TSS", "Sensitivity", "Specificity", "CBI", 'RMSE',"R2"  ))

glm_out <- rbind(c("GLM.sdm",
                   mean(unlist(lapply(m1@models$PA$glm, function(x){ mean(x@evaluation$test.dep@statistics$AUC) })), na.rm=T),
                   mean(unlist(lapply(m1@models$PA$glm, function(x){ mean(x@evaluation$test.dep@threshold_based$TSS[2]) })), na.rm=T),
                   mean(unlist(lapply(m1@models$PA$glm, function(x){ mean(x@evaluation$test.dep@threshold_based$sensitivity[2]) })), na.rm=T),
                   mean(unlist(lapply(m1@models$PA$glm, function(x){ mean(x@evaluation$test.dep@threshold_based$specificity[2]) })), na.rm=T),
                   mean(unlist(lapply(m1@models$PA$glm, function(x){ mean(ecospat::ecospat.boyce(fit = x@evaluation$test.dep@predicted, obs=x@evaluation$test.dep@predicted[which(x@evaluation$test.dep@observed==1)],
                                                                                                       nclass=0,
                                                                                                       window.w="default", res=100, PEplot = FALSE)$cor)})), na.rm=T),
                   mean(unlist(lapply(m1@models$PA$glm, function(x){Metrics::rmse(x@evaluation$test.dep@observed,x@evaluation$test.dep@predicted)})), na.rm=T),
                   mean(unlist(lapply(m1@models$PA$glm, function(x){modEvA::RsqGLM(obs=x@evaluation$test.dep@observed, 
                                                                                         pred=x@evaluation$test.dep@predicted, 
                                                                                         plot=FALSE)$Nagelkerke})), na.rm=T)))# (pseudo) R-squared 

brt_out <- rbind(c("BRT.sdm",
                  mean(unlist(lapply(m1@models$PA$brt, function(x){ mean(x@evaluation$test.dep@statistics$AUC) })), na.rm=T),
                  mean(unlist(lapply(m1@models$PA$brt, function(x){ mean(x@evaluation$test.dep@threshold_based$TSS[2]) })), na.rm=T),
                  mean(unlist(lapply(m1@models$PA$brt, function(x){ mean(x@evaluation$test.dep@threshold_based$sensitivity[2]) })), na.rm=T),
                  mean(unlist(lapply(m1@models$PA$brt, function(x){ mean(x@evaluation$test.dep@threshold_based$specificity[2]) })), na.rm=T),
                  mean(unlist(lapply(m1@models$PA$brt, function(x){ mean(ecospat::ecospat.boyce(fit = x@evaluation$test.dep@predicted,
                                                                                                     obs=x@evaluation$test.dep@predicted[which(x@evaluation$test.dep@observed==1)],nclass=0,
                                                                                                     window.w="default", res=100, PEplot = FALSE)$cor) 
                  })), na.rm=T),
                  mean(unlist(lapply(m1@models$PA$brt, function(x){Metrics::rmse(x@evaluation$test.dep@observed,x@evaluation$test.dep@predicted)})), na.rm=T),
                  mean(unlist(lapply(m1@models$PA$brt, function(x){mean(x@object$rsq)})), na.rm=T)))# (pseudo) R-squared 

myComp.sdm[1,] <- glm_out
myComp.sdm[2,] <- brt_out
myComp.sdm[,2:length(myComp.sdm)]=lapply(myComp.sdm[,2:length(myComp.sdm)], as.numeric) #needed to change it up due to chr
myComp.sdm

##---- 4.4. Spatial predictions----
#if you want to get directly the averaged predictions, then
p2m <- predict(m1, 
               newdata = bio,
               overwrite = TRUE, mean=T)


plot(p2m)
# ensemble based on a Weighted averaging that is weighted using TSS statistic
e1 <- ensemble(m1, newdata = bio, 
               overwrite = TRUE, 
               setting = list(method='weighted',stat='tss'))
plot(e1)

##---- 4.5. Convert to binary map----
df <- data.frame(PA = as.data.frame(d)$PA, 
                 sdm::coords(d)
        )
xy <- as.matrix(df[, c("long", "lat")])
p <- terra::extract(e1, xy)
ev <- sdm::evaluates(df$PA, p[,1])

#extract max sens ans spec threshold
th <- ev@threshold_based$threshold[2]
e1pa<-e1
e1pa[] <- ifelse(e1[] >= th, 1, 0)
plot(e1pa, main="Acer pseudoplatanus")

#---- 5. Load chorological distribution map----
acer <- read_sf("Chorological data for the main European woody species/chorological_maps_dataset_20230824/chorological_maps_dataset/Acer pseudoplatanus/shapefiles/Acer_pseudoplatanus_plg.shp")
# no-data countries mask
mask <- geodata::gadm(country=c("ITA", "FRA", "LIE", "CH", "AUT",  "DEU", "SVN", "CZE", "HUN", 
                                "HRV", "MNE","BIH", "MKD", "GRC", "XK"), 
                      level=0, 
                      path = outdir) %>% 
  project("EPSG:4326") 
mask <- st_as_sf(mask) 
acer <- sf::st_intersection(acer, mask %>% filter(GID_0=="ITA"))

myPres <- species_data.f$Tr_data %>% 
  select(PA, long, lat) %>% 
  st_as_sf(coords = c("long", "lat"), crs=4326)

e1pa.df<- terra::as.data.frame(e1pa, xy = TRUE) %>% 
  mutate(ensemble_weighted = as.factor(ensemble_weighted)) %>% 
  as_tibble()
myCols<- c("0" = "lightgrey", "1" ="#1E88E5", "Chorological map"="#004D40")

p <- ggplot()+
  geom_sf(data=mask, col="white", linewidth=0.2)+
  geom_tile(data = e1pa.df, aes(x = x, y = y, fill = ensemble_weighted))+
  geom_sf(data=acer, aes(fill= "Chorological map"), col="#8F9779", alpha=0.5, linewidth=0.8)+ ##708238
  geom_sf(data=myPres,  aes(color= "GBIF Observations"),size=.5)+
  scale_fill_manual(values=myCols)+
  scale_color_manual(values= "#D81B60")+
  labs(x="Longitude",y="Latitude", fill= "Presence-Absence", col="")+
  theme_classic()+
  theme(legend.position = "bottom",  
        text = element_text(size=12),
        legend.text = element_text( size = 12))+
  coord_sf(xlim = c(6.5, 18.8), ylim = c(35, 47.5), expand = FALSE)

outname <- paste("figures/", "AcerPseudoplatanus_sdm_", Sys.Date(),".png", sep="")
ggsave(p, filename = outname, width = 18, height = 8, device='png', dpi=320)
