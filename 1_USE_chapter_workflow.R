# USE book chapter, sensitivity analysis code
# devtools::install_github("danddr/USE")
setwd("/home/daniele/Documents/working_files/USE_bookChapter/")
library(sf)
library(terra)
library(geodata)
library(USE)
library(hypervolume)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
source("SpatialProba.R")
outdir <- getwd()

#---- 1. Prepare data ----
#download country area and bioclimatic variables
ita <- geodata::gadm(country="ITA", 
                     level=0, 
                     path = outdir)
BioData <- geodata::worldclim_global(var = 'bio', 
                                     res = 5, 
                                     path = outdir)
BioData <- terra::crop(BioData, ita, mask = TRUE )
names(BioData) <- paste0("bio", seq_along(names(BioData)))

# set coefficients for virtual species
myCoef <- c(-8, 1.7, -0.1, 0.0009)
names(myCoef) <- c("intercept", "bio1", "quad_bio1", "bio12" )
# simulate virtual species
VS_proba_layer <- SpatialProba(coefs = myCoef,
                               env.rast = BioData, 
                               quadr_term = "bio1", 
                               marginalPlots = TRUE)

# Convert to presence-absence using a random binomial draw
N_VS_layers <- terra::app(VS_proba_layer$rast, fun = function(x) {
  replicate(10, rbinom(n = length(x), 1, x))
})
names(N_VS_layers) <- rep("pa", nlyr(N_VS_layers))

# save figure plots for Figure 1
ggplot2::ggsave("figures/margEff.pdf", VS_proba_layer$margEff, width = 12, height = 8)
pdf("figures/VS_pa.pdf", p, width = 12, height = 8)
plot(N_VS_layers[[1]], plg=list(x="topright"), axes=FALSE)
dev.off()

# Define sensitivity matrix for sample prevalence and buffer size
sens.m <- expand.grid(samplePrev = c(0.25, 0.5, 1),
                      BufferSize = c(15000, 25000, 50000))

# Setting Environmental grid resolution
rpc <- USE::rastPCA(BioData, stand = TRUE)
dt <- na.omit(as.data.frame(rpc$PCs[[c("PC1", "PC2")]], xy = TRUE))
dt <- sf::st_as_sf(dt, coords = c("PC1", "PC2"))
myRes <- USE::optimRes(sdf= dt,
                       grid.res= c(1:10),
                       perc.thr = 20,
                       showOpt = FALSE, 
                       cr=5)

#---- 2. Sensitivity analysis: Hypervolume overlaps and range estimates ----
dir.create("outputs/", recursive=TRUE)
mySensOut <- list()
for(sens in 1:nrow(sens.m)){ #
  # sens=1
  # retrieve sample prevalence
  sample.prev <- unname(sens.m[sens, 1]) 
  # retrieve  buffer size
  buff.dim <- unname(sens.m[sens, 2]) #udm [m]
  message("#-------------------------------------------------------------------------------------------#")
  message("#---------- Processing sample prev ", sample.prev, " and buff size ", buff.dim, " ----------#")
  message("#-------------------------------------------------------------------------------------------#")
    # get the presences from each presence-absence raster
  presList <- lapply(1:nlyr(N_VS_layers), function(x){terra::vect(subset(as.data.frame(N_VS_layers[[x]], xy=TRUE), pa==1), 
                                                      geom = c("x", "y"), 
                                                      crs = "EPSG:4326")})
  # Generate pseudo-absences in the geographical space accordingly to sample prevalence
  prevList <- do.call(c, lapply(presList, nrow))/sample.prev
  # sampling the geographic space: random sampling
  randList <- lapply(prevList, function(x){terra::spatSample(BioData$bio1, method= 'random', size=x, xy=TRUE, na.rm=TRUE )})
  randList <- lapply(randList, setNames, c("x", "y", "pa"))
  randList <- lapply(randList, function(x){x$pa <- 0 ;  return(x)})
  randList <- lapply(randList, function(x){terra::vect(x, geom = c("x", "y"), crs = "EPSG:4326")})
  randList <- Map(rbind, presList, randList)
  # sampling the geographic space: buffer-out random sampling
  buf.surface <- lapply(presList, function(x){
    terra::mask(BioData$bio1, terra::buffer(x, width = buff.dim), inverse=TRUE)
  })
  buffList <- lapply(1:length(N_VS_layers), function(x){terra::spatSample(buf.surface[[x]], method= 'random', size=prevList[x], xy=TRUE, na.rm=TRUE )})
  buffList <- lapply(buffList, setNames, c("x", "y", "pa"))
  buffList <- lapply(buffList, function(x){x$pa <- 0 ;  return(x)})
  buffList <- lapply(buffList, function(x){terra::vect(x, geom = c("x", "y"), crs = "EPSG:4326")})
  buffList <- Map(rbind, presList, buffList)
  # sampling the environmental space: uniform sampling
  uniformList <- lapply(presList, function(x){sf::st_drop_geometry(USE::paSampling(BioData, 
                                                                                     pres=x, 
                                                                                     thres=0.75,
                                                                                     H=NULL,
                                                                                     grid.res=myRes$Opt_res,
                                                                                     prev=sample.prev, 
                                                                                     sub.ts=FALSE,
                                                                                     plot_proc=FALSE)[, c("x", "y")])})
  uniformList <- lapply(uniformList, function(x){x$pa <- 0 ;  return(x)})
  uniformList <- lapply(uniformList, function(x){terra::vect(x, geom = c("x", "y"), crs = "EPSG:4326")})
  uniformList <- Map(rbind, presList, uniformList)
  # put together the the presences and the pseudo-absences for a given combination of prevalence and buffer size
  myPAlist <- list(
    "random"= randList,
    "buffer"= buffList,
    "uniform" = uniformList
  )
  # Estimate Hypervolume overlap and environmental range size 
  message("Starting hypervolumes overlap...")
  myOut_tmp <- list() 
  for(i in 1:length(myPAlist)){
    # i=1
    cat(paste("\nProcessing sampling strategy: ", names(myPAlist)[i]))
    subList <- myPAlist[[i]]
    tmpOut <- list()
    for(j in 1:length(subList)){
      # j=1
      cat(paste("\nProcessing sampling strategy: ", names(myPAlist)[i]), "; VS", j, "of",length(subList) )
      df <- data.frame(pa=subList[[j]]$pa, 
                       terra::extract(x=BioData,
                                      y=subList[[j]], 
                                      ID=FALSE))
      PCAenv <- data.frame(pa=df$pa,princomp(scale(df[,-1]))$scores[,1:2])
      hypAbs <- hypervolume::hypervolume(PCAenv[PCAenv$pa==0,-1],verbose=FALSE)
      hypPres <- hypervolume::hypervolume(PCAenv[PCAenv$pa==1,-1],verbose=FALSE)
      ovrlp <- hypervolume::get_volume(hypervolume::hypervolume_set(hypPres,hypAbs, check.memory=FALSE,verbose=FALSE))[[3]] #3 is intersection
      # compute PC1-2 range   
      rpc <- USE::rastPCA(env.rast = BioData, 
                          nPC = 2, stand = TRUE)
      pcstack <- c(rpc$PCs$PC1, rpc$PCs$PC2)
      tmpVect <- na.omit(terra::extract(pcstack, subList[[j]], df=TRUE))
      rangePC1 <- round(range(tmpVect$PC1)[2] - range(tmpVect$PC1)[1],3)
      rangePC2 <- round(range(tmpVect$PC2)[2] - range(tmpVect$PC2)[1],3)
      tmpOut[[j]] <- data.frame(Sampling = names(myPAlist)[i], 
                              VS = j, 
                              Exp.prev = sample.prev, 
                              buff.size = buff.dim, 
                              overlap = ovrlp, 
                              rangePC1 = rangePC1, 
                              rangePC2 = rangePC2)
    }
    myOut_tmp[[i]] <- do.call(rbind, tmpOut)
  }
  myOut_tmp<-do.call(rbind, myOut_tmp)
  mySensOut[[sens]]<-myOut_tmp
}
mySensOut <- do.call(rbind, mySensOut)
outname <- paste0("outputs/sensitivityOutput_", Sys.Date(), ".RDS")
saveRDS(mySensOut, outname)

#---- 3. Plotting ----
dir.create("figures/", recursive=TRUE)
mySensOut <- readRDS("outputs/sensitivityOutput_2024-10-26.RDS")
mySensOut <- mySensOut %>%
  mutate(buff.size = paste("Buffer radius:",buff.size/1000, "km"), 
         buff.size = factor(buff.size, levels= c("Buffer radius: 15 km", "Buffer radius: 25 km", "Buffer radius: 50 km")), 
         Sampling = recode(Sampling, uniform = "Uniform", random = "Random", buffer = "Buffer-out"), 
         Sampling = factor(Sampling, levels = c("Uniform", "Random", "Buffer-out")), 
         Exp.prev = as.factor(paste("Prevalence:", Exp.prev)))

#Class overlap                              
p <- mySensOut %>%   
  ggplot(aes(y = overlap, x = Sampling, color = Sampling))+
  geom_boxplot()+
  stat_summary(fun = median, geom = "point", size = 2) +
  scale_color_manual(breaks = c("Uniform", "Random", "Buffer-out"),
                     values=c("#0072B2", "#E69F00", "#CC79A7"))+
  labs(x ="",  y = "Overlap", color = "Sampling strategy")+
  scale_fill_brewer(palette = 'Set1')+
  scale_x_discrete(labels = NULL, breaks = NULL) +
  facet_grid(Exp.prev ~ buff.size) +
  theme_light()+
  theme(aspect.ratio = 1.2, 
        legend.background=element_blank(),
        panel.grid = element_blank(),
        legend.position = 'bottom',
        text = element_text(size=14), 
        strip.text = element_text(size=12),
        legend.text = element_text(size=14,angle = 0), 
        legend.title = element_text(size=14),
        legend.key.size = unit(1, 'cm'))
p
outname <- paste("figures/", "classOverlap_boxplot", Sys.Date(),".png", sep="")
ggsave(p, filename = outname, width = 18, height = 8, device='png', dpi=320)

# PCs-range
mySensOut <- mySensOut %>% 
  select(Sampling, Exp.prev, buff.size ,rangePC1, rangePC2) %>% 
  pivot_longer(c("rangePC1", "rangePC2")) %>% 
  rename(range = value, PCs = name) %>% 
  mutate(PCs = ifelse(PCs == "rangePC1", "PC1", "PC2"))

# compute bootstraps
unique_combo <- unique(mySensOut[, c("Sampling", "PCs", "Exp.prev", "buff.size" )])
outRange<-list()
for(i in 1:nrow(unique_combo)){
  # i=1
  df <- mySensOut %>% 
    filter(Sampling==unique_combo[i, "Sampling"]$Sampling,
           PCs==unique_combo[i, "PCs"]$PCs,
           Exp.prev==unique_combo[i, "Exp.prev"]$Exp.prev,
           buff.size==unique_combo[i, "buff.size"]$buff.size) 
  ourRange <- Hmisc::smean.cl.boot(df$range, B=2000) 
  ourRange <- data.frame(Sampling = unique_combo[i, "Sampling"]$Sampling,
                         PCs = unique_combo[i, "PCs"]$PCs,
                         Exp.prev = unique_combo[i, "Exp.prev"]$Exp.prev,
                         buff.size = unique_combo[i, "buff.size"]$buff.size, 
                         Mean = ourRange[1],
                         lci = ourRange[2],
                         uci = ourRange[3]
  )
  outRange[[i]] <- ourRange
  
}
outRange<-do.call(rbind, outRange)

#plot
pc1 <-outRange %>%
  filter(PCs == "PC1") %>% 
  ggplot(aes(y = Mean, x = Sampling, color = Sampling))+
  geom_point()+
  geom_errorbar(aes(ymin=lci, ymax=uci))+
  scale_color_manual(breaks = c("Uniform", "Random", "Buffer-out"),
                     values=c("#0072B2", "#E69F00", "#CC79A7"))+
  labs(x="",  y="Total range (95% CI)", color="Sampling strategy")+
  scale_x_discrete(labels = NULL, breaks = NULL) +
  facet_grid(Exp.prev ~ buff.size) +
  theme_light()+
  theme(aspect.ratio = 1.2, 
        legend.background=element_blank(),
        panel.grid = element_blank(),
        legend.position = 'bottom',
        text = element_text(size=14), 
        strip.text = element_text(size=12),
        legend.text = element_text(size=14,angle = 0), 
        legend.title = element_text(size=14),
        legend.key.size = unit(1, 'cm'))

pc2 <- outRange %>%
  filter(PCs == "PC2") %>% 
  ggplot(aes(y = Mean, x = Sampling, color = Sampling))+
  geom_point()+
  geom_errorbar(aes(ymin=lci, ymax=uci))+
  scale_color_manual(breaks = c("Uniform", "Random", "Buffer-out"),
                     values=c("#0072B2", "#E69F00", "#CC79A7"))+
  labs(x="",  y="Total range (95% CI)", color="Sampling strategy")+
  scale_x_discrete(labels = NULL, breaks = NULL) +
  facet_grid(Exp.prev ~ buff.size) +
  theme_light()+
  theme(aspect.ratio = 1.2, 
        legend.background=element_blank(),
        panel.grid = element_blank(),
        legend.position = 'bottom',
        text = element_text(size=14), 
        strip.text = element_text(size=12),
        legend.text = element_text(size=14,angle = 0), 
        legend.title = element_text(size=14),
        legend.key.size = unit(1, 'cm'))

pp <- ggarrange(pc1, pc2,  
              common.legend = TRUE, 
              legend="bottom", 
              labels = c("A", "B"))
pp
outname <- paste("figures/", "PCs_ranges_", Sys.Date(),".png", sep="")
ggsave(pp, filename = outname, width = 18, height = 8, device='png', dpi=320)