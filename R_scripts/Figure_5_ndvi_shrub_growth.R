# Analyses for Figure 3
# January 2019
# Isla Myers-Smith

#packrat::init()

# Packages ----
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(viridis)

# Library for detrending data - has conflict with Phenex package - must uninstall Phenex before use
library(dplR)

# Libraries for Bayesian shrubring models
library(MCMCglmm)

# Libraries for mixed models
library(nlme)
library(MuMIn)

# Function to extract MCMCglmm model summary outputs ----
clean.MCMC <- function(x) {
  sols <- summary(x)$solutions  # pull out relevant info from model summary
  Gcovs <- summary(x)$Gcovariances
  Rcovs <- summary(x)$Rcovariances
  
  fixed <- data.frame(row.names(sols), sols, row.names = NULL)  # convert to dataframes with the row.names as the first col
  random <- data.frame(row.names(Gcovs), Gcovs, row.names = NULL)
  residual <- data.frame(row.names(Rcovs), Rcovs, row.names = NULL)
  
  names(fixed)[names(fixed) == "row.names.sols."] <- "variable"  # change the columns names to variable, so they all match
  names(random)[names(random) == "row.names.Gcovs."] <- "variable"
  names(residual)[names(residual) == "row.names.Rcovs."] <- "variable"
  
  fixed$effect <- "fixed"  # add ID column for type of effect (fixed, random, residual)
  random$effect <- "random"
  residual$effect <- "residual"
  
  modelTerms <- as.data.frame(bind_rows(fixed, random, residual))  # merge it all together
}

# Function for when models have no random effects
clean.MCMC.2 <- function(x) {
  sols <- summary(x)$solutions  # pull out relevant info from model summary
  Gcovs <- summary(x)$Gcovariances
  Rcovs <- summary(x)$Rcovariances
  
  fixed <- data.frame(row.names(sols), sols, row.names = NULL)  # convert to dataframes with the row.names as the first col
  residual <- data.frame(row.names(Rcovs), Rcovs, row.names = NULL)
  
  names(fixed)[names(fixed) == "row.names.sols."] <- "variable"  # change the columns names to variable, so they all match
  names(residual)[names(residual) == "row.names.Rcovs."] <- "variable"
  
  fixed$effect <- "fixed"  # add ID column for type of effect (fixed, random, residual)
  residual$effect <- "residual"
  
  modelTerms <- as.data.frame(bind_rows(fixed, residual))  # merge it all together
}

getName.MCMC <- function(x) deparse(substitute(x))  # adding the model name

# For models with two random effects
prior1 <- list(R = list(V = 1, nu = 0.002), 
               G = list(G1 = list(V = 1, nu = 1, 
                                  alpha.mu = 0, alpha.v = 10000)))

prior2 <- list(R = list(V = 1, nu = 0.002), 
               G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.v = 10000), 
                        G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.v = 10000)))

# Customised ggplot2 theme function ----
theme_QHI <- function(){
  theme_bw() +
    theme(axis.text = element_text(size = 16), 
          axis.title = element_text(size = 20),
          axis.text.x = element_text(angle = -45, hjust = -0.05),
          axis.line.x = element_line(color="black"), axis.line.y = element_line(color="black"),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),                                          
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),  
          plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
          plot.title = element_text(size=20, vjust=1, hjust=0),
          legend.text = element_text(size=16, face="italic"),          
          legend.title = element_blank(),                              
          legend.position = c(0.9, 0.9), 
          legend.key = element_blank(),
          legend.background = element_rect(color = "black", fill = "transparent", size = 4, linetype="blank"))
}


# MODIS greening trends ---- 

# Update the coordinates for the trends to be the extents of the drone plots:
load("data/greenup.qhi.kanger.RData")

qikiq_greenup <- greenup.all %>% filter(site == "QHI")
kanger_greenup <- greenup.all %>% filter(site == "Kanger")

# Qikiqtaruk peak greenness trends

QHI_greening_trend <- ggplot(qikiq_greenup) + 
  geom_line(aes(x=year, y=ndvi.date.max), alpha = 0.5, size = 0.5, show.legend=T) +
  geom_point(aes(x = year, y = ndvi.date.max), colour = "#41A002", alpha = 0.5, size = 3) +
  geom_smooth(method=lm, aes(x = year, y = ndvi.date.max), colour = "#41A002", fill = "#41A002", alpha = 0.5, show.legend=F, linetype = "dotted") +
  ylab("Annual Max. NDVI\n") +
  xlab("") +
  scale_y_continuous(breaks = c(0.3, 0.4, 0.5, 0.6, 0.7), limits = c(0.3, 0.7)) +
  scale_x_continuous(breaks = c(2000, 2005, 2010, 2015), limits = c(2000, 2017)) +
  theme_QHI()

ggsave("plots/QHI_greening_trend.png", plot = QHI_greening_trend, width = 4, height = 4, units = "in", dpi = 300)

Kanger_greening_trend <- ggplot(kanger_greenup) + 
  geom_line(aes(x=year, y=ndvi.date.max), alpha = 0.5, size = 0.5, show.legend=T) +
  geom_point(aes(x = year, y = ndvi.date.max), colour = "#41A002", alpha = 0.5, size = 3) +
  geom_smooth(method=lm, aes(x = year, y = ndvi.date.max), colour = "#41A002", fill = "#41A002", alpha = 0.5, show.legend=F, linetype = "dotted") +
  ylab("Annual Max. NDVI\n") +
  xlab("") +
  scale_y_continuous(breaks = c(0.3, 0.4, 0.5, 0.6, 0.7), limits = c(0.3, 0.7)) +
  scale_x_continuous(breaks = c(2000, 2005, 2010, 2015), limits = c(2000, 2017)) +
  theme_QHI()

ggsave("plots/Kanger_greening_trend.png", plot = Kanger_greening_trend, width = 4, height = 4, units = "in", dpi = 300)


# Qikiqtaruk shrub growth ----

QHI_shrub_height <- read.csv("data/qhi_point_fraiming_ITEX_1999-2017.csv")
QHI_shrub_height$Height <- as.numeric(as.character(QHI_shrub_height$Height))
QHI_shrub_height <- QHI_shrub_height %>% filter(SPP == "SALPUL" & Height >= 0) %>% group_by(YEAR) %>% summarise(mean_height = mean(Height), height_SD = sd(Height))

QHI_shrub_growth <- ggplot(QHI_shrub_height) + 
  geom_point(aes(x=YEAR, y=mean_height), colour = "#41A002", alpha = 0.5, size = 3) +
  geom_errorbar(aes(x=YEAR, ymin = mean_height-height_SD, ymax=mean_height+height_SD), colour = "#41A002", alpha = 0.5) +
  geom_smooth(method=lm, aes(x=YEAR, y=mean_height), colour = "#41A002", fill = "#41A002", alpha = 0.5, show.legend=F, linetype = "dotted") +
  ylab("Shrub height (cm)\n") +
  xlab("") +
  scale_y_continuous(breaks = c(0, 5, 10, 15), limits = c(0, 15)) +
  scale_x_continuous(breaks = c(2000, 2005, 2010, 2015), limits = c(1998, 2017)) +
  theme_QHI()

ggsave("plots/QHI_shrub_growth.png", plot = QHI_shrub_growth, width = 4, height = 4, units = "in", dpi = 300)


# Qikiqtaruk shrub abundance ----

# The full dataset is archived at:
# Isla H. Myers-Smith, Daskalova, G. N., Bjorkman, A. D. & Thomas, H. J. D. ShrubHub/QikiqtarukHub: QikiqtarukHub_v1.0. (Zenodo, 2018). doi:10.5281/zenodo.2397996
# https://zenodo.org/record/2397996#.XFIU8fx7knM
# https://github.com/ShrubHub/QikiqtarukHub

# 2. Please cite: Myers-Smith, I. H. & et al. 2019. Eighteen years of ecological monitoring reveals multiple lines of evidence for tundra vegetation change. Ecological Monographs. http://onlinelibrary.wiley.com/doi/10.1002/ecm.1351/full

QHI_shrub_cover <- read.csv("data/QHI_cover_1999_2017_ITEX.csv")

#Subset to Salix pulchra
QHI_shrub_cover <- QHI_shrub_cover %>% filter(name == "Salix pulchra") %>% dplyr::select(year, plot, cover)

# calculate mean cover (and using SE rather than SD to match Eric Post's data)
QHI_mean_cover <- QHI_shrub_cover %>% group_by(year) %>%
   summarise(cover_mean = mean(cover),
             cover_SE = sd(cover)/(sqrt(length(cover))))

QHI_shrub_abund <- read.csv("data/qhi_point_fraiming_ITEX_1999-2017.csv")
QHI_shrub_abund <- QHI_shrub_abund %>% filter(SPP == "SALPUL" | SPP == "Salix pulchra") %>% group_by(YEAR, PLOT, SPP) %>% summarise(abund = sum(Abundance)) %>% group_by(YEAR, SPP) %>% summarise(abund_mean = mean(abund), abund_SE = sd(abund)/(sqrt(length(abund))))
QHI_shrub_abund$SPP[QHI_shrub_abund$SPP == "SALPUL"] <- "Salix pulchra"

# Plot
QHI_shrub_abund <- ggplot(QHI_shrub_abund) + 
   geom_point(aes(x = YEAR, y = abund_mean, colour = SPP), alpha = 0.5, size = 3) +
   geom_errorbar(aes(x = YEAR, ymin = abund_mean - abund_SE, ymax= abund_mean + abund_SE, colour = SPP), alpha = 0.5) +
   geom_smooth(method = glm, aes(x = YEAR, y = abund_mean), colour = "#41A002", fill = "#41A002", alpha = 0.5, show.legend=F, linetype = "dotted") +
   ylab("Shrub abundance\n") +
   xlab("") +
   scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100), limits = c(0, 100)) +
   scale_x_continuous(breaks = c(2000, 2005, 2010, 2015), limits = c(1998, 2017)) +
  scale_colour_manual(values = c("#41A002"), labels = c("S. pulchra")) +
  #scale_fill_manual(values = c("#41A002"), labels = c("S. pulchra")) +
  theme_bw() +
  theme(axis.text = element_text(size = 16), 
        axis.title = element_text(size = 20),
        axis.text.x = element_text(angle = -45, hjust = -0.05),
        axis.line.x = element_line(color="black"), axis.line.y = element_line(color="black"),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),                                          
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),  
        plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
        plot.title = element_text(size=20, vjust=1, hjust=0),
        legend.text = element_text(size=14, face="italic"),          
        legend.title = element_blank(),                              
        legend.position = c(0.75, 0.12), 
        legend.key = element_blank(),
        legend.background = element_rect(color = "black", fill = "transparent", size = 4, linetype="blank"))

ggsave("plots/QHI_shrub_abund.png", plot = QHI_shrub_abund, width = 4, height = 4, units = "in", dpi = 300)


# Kanger shrub abundance ----

# Researchers who are interested in using the Kangerlussuaq shrub abundance data for purposes other than reproducing the results of Myers-Smith et al. analyses are encouraged to
# 1. Access the original data from the Arctic Data Center, where the sampling design is described in detail:
# Eric Post and Christian Pedersen. Low Arctic monitoring of plant community composition and dynamics. Arctic Data Center. doi:10.5065/D6542KRH. https://arcticdata.io/catalog/view/doi:10.5065/D6542KRH

# 2. Include the following citation in any resulting publications: Post, E.  2013.  Erosion of community diversity and stability by herbivore removal under warming.  Proc. Roy. Soc. Lond. Ser. B., 208: 20122722

kanger_shrub_abund <- read.csv("data/kanger_abundance.csv")

# Plot
kanger_abund <- ggplot(kanger_shrub_abund, group = Species) + 
   geom_point(aes(x=Year, y=abund_mean, colour = Species), alpha = 0.5, size = 3) +
   geom_errorbar(aes(x=Year, ymin = abund_mean - abund_SE, ymax= abund_mean + abund_SE, colour = Species), alpha = 0.5) +
   geom_smooth(method=loess, aes(x=Year, y=abund_mean, colour = Species, fill = Species), alpha = 0.5, show.legend=F, linetype = "dotted") +
   ylab("Shrub abundance\n") +
   xlab("") +
   scale_y_continuous(breaks = c(0, 20, 40, 60), limits = c(0, 60)) +
   scale_x_continuous(breaks = c(2000, 2005, 2010, 2015), limits = c(1998, 2017)) +
   scale_colour_manual(values = c("#EEB422", "#41A002"), labels = c("B. nana", "S. glauca")) +
   scale_fill_manual(values = c("#EEB422", "#41A002")) +
   theme_bw() +
   theme(axis.text = element_text(size = 16), 
         axis.title = element_text(size = 20),
         axis.text.x = element_text(angle = -45, hjust = -0.05),
         axis.line.x = element_line(color="black"), axis.line.y = element_line(color="black"),
         panel.border = element_blank(),
         panel.grid.major.x = element_blank(),                                          
         panel.grid.minor.x = element_blank(),
         panel.grid.minor.y = element_blank(),
         panel.grid.major.y = element_blank(),  
         plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
         plot.title = element_text(size=20, vjust=1, hjust=0),
         legend.text = element_text(size=14, face="italic"),          
         legend.title = element_blank(),                              
         legend.position = c(0.75, 0.12), 
         legend.key = element_blank(),
         legend.background = element_rect(color = "black", fill = "transparent", size = 4, linetype="blank"))
   
ggsave("plots/kanger_shrub_abund.png", plot = kanger_abund, width = 4, height = 4, units = "in", dpi = 300)


# Qikiqtaruk shrub dendro ----

# The full dataset is archived at:
# Isla H. Myers-Smith, Daskalova, G. N., Bjorkman, A. D. & Thomas, H. J. D. ShrubHub/QikiqtarukHub: QikiqtarukHub_v1.0. (Zenodo, 2018). doi:10.5281/zenodo.2397996
# https://zenodo.org/record/2397996#.XFIU8fx7knM
# https://github.com/ShrubHub/QikiqtarukHub

# Please cite: Myers-Smith, I. H. & et al. 2019. Eighteen years of ecological monitoring reveals multiple lines of evidence for tundra vegetation change. Ecological Monographs. http://onlinelibrary.wiley.com/doi/10.1002/ecm.1351/full

QHI_shrub <- read.csv("data/qhi_dendro_2017.csv")
QHI_dendro_sp <- QHI_shrub %>% filter(Sp == "Salix pulchra") %>% filter(Year > 1998)

# Detach Phenex to allow the dplR package to work
detach("package:phenex", unload=TRUE)

QHI_shrub_detrend <- QHI_dendro_sp %>% group_by(IndivUN) %>% mutate(detrend_rw = detrend.series(rw)$Spline)

QHI_ndvi_detrend <- cbind.data.frame(qikiq_greenup$year, detrend.series(qikiq_greenup$ndvi.date.max)$Spline)
names(QHI_ndvi_detrend) <- c("year", "detrend_ndvi")

QHI_shrub_detrend <- inner_join(QHI_shrub_detrend, QHI_ndvi_detrend, by = c("Year" = "year"))

# Models fit using MCMCglmm
QHI_model <- MCMCglmm(detrend_rw ~ detrend_ndvi, 
                            random = ~Year + detrend_ndvi:IndivUN, data = QHI_shrub_detrend, 
                            family = "gaussian", pr = TRUE, nitt = 100000, burnin = 20000, prior = prior2)

# Models fit using nlme
# QHI_model <- lme(detrend_rw ~ detrend_ndvi, random = list(~detrend_ndvi|IndivUN, ~1|Year), cor=corAR1(), data = QHI_shrub_detrend, method = "REML", control=list(maxIter=1000))

summary(QHI_model)

# Calculates pseudo R2 for nlme models
# r.squaredGLMM(QHI_model)

fit_QHI_model_pred <- as.data.frame(predict.MCMCglmm(QHI_model, QHI_shrub_detrend, interval = "prediction"))
fit_QHI_model <- cbind.data.frame(Year = QHI_shrub_detrend$Year, IndivUN = QHI_shrub_detrend$IndivUN, detrend_rw = QHI_shrub_detrend$detrend_rw, detrend_ndvi = QHI_shrub_detrend$detrend_ndvi, fit = fit_QHI_model_pred$fit, lwr = fit_QHI_model_pred$lwr, upr = fit_QHI_model_pred$upr) %>% arrange(desc(detrend_ndvi)) %>% filter(row_number() %% 10 == 0)

# Mixed model figure

QHI_rw_ndvi <- ggplot() + 
  geom_point(data = QHI_shrub_detrend, aes(x=detrend_ndvi, y=detrend_rw), colour = "#41A002", alpha = 0.5, size = 3) +
  geom_line(data = fit_QHI_model, aes(x=detrend_ndvi, y=fit), colour = "#41A002", size = 1) +
  geom_ribbon(data = fit_QHI_model, aes(x=detrend_ndvi, ymin=lwr, ymax=upr), fill = "#41A002", alpha = 0.5) +
  ylab("Detrended ring width\n") +
  xlab("\nDetrended NDVI") +
  scale_y_continuous(breaks = c(-1, 0, 1, 2, 3, 4), limits = c(-1, 4)) +
  scale_x_continuous(breaks = c(0.8, 0.9, 1.0, 1.1, 1.2), limits = c(0.8, 1.2)) +
  theme_QHI()
  
ggsave("plots/QHI_rw_ndvi.png", plot = QHI_rw_ndvi, width = 4, height = 4, units = "in", dpi = 300)

# Dendro plot

QHI_dendro <- QHI_shrub %>% filter(Year > 1998) %>% group_by(Year, Sp) %>% 
  summarise(mean_rw = mean(rw), rw_SD = sd(rw), rw_SE = (sd(rw)/sqrt(length(unique(IndivUN)))))

QHI_dendro_plot <- ggplot(QHI_dendro, group = Sp) + 
  geom_errorbar(aes(x=Year, ymin = mean_rw-rw_SE, ymax=mean_rw+rw_SE, colour = Sp), alpha = 0.5) +
  geom_path(aes(x=Year, y=mean_rw, colour = Sp), size = 1) +
  geom_point(aes(x=Year, y=mean_rw, colour = Sp), alpha = 0.5, size = 3) +
  ylab("Shrub ring width (mm)\n") +
  xlab("") +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
  scale_x_continuous(breaks = c(2000, 2005, 2010, 2015), limits = c(1998, 2017)) +
  theme_bw() +
  theme(axis.text = element_text(size = 16), 
        axis.title = element_text(size = 20),
        axis.text.x = element_text(angle = -45, hjust = -0.05),
        axis.line.x = element_line(color="black"), axis.line.y = element_line(color="black"),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),                                          
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),  
        plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
        plot.title = element_text(size=20, vjust=1, hjust=0),
        legend.text = element_text(size=14, face="italic"),          
        legend.title = element_blank(),                              
        legend.position = c(0.75, 0.92), 
        legend.key = element_blank(),
        legend.background = element_rect(color = "black", fill = "transparent", size = 4, linetype="blank"))

ggsave("plots/QHI_dendro.png", plot = QHI_dendro_plot, width = 4, height = 4, units = "in", dpi = 300)


# Kanger shrub dendro ----

# Researchers who are interested in using the Kangerlussuaq shrub annual growth ring data for purposes other than reproducing the results of Myers-Smith et al. analyses are encouraged to
# 1. Access the original data from the Arctic Data Center, where the sampling design is described in detail:
# Betula nana: Patrick Sullivan. 2016. Betula nana ring widths. Arctic Data Center. doi:10.18739/A28Q18. https://arcticdata.io/catalog/view/doi:10.18739/A28Q18
# Salix glauca: Patrick Sullivan. 2016. Salix glauca ring widths. Arctic Data Center. doi:10.18739/A24X0Q. https://arcticdata.io/catalog/view/doi:10.18739/A24X0Q
# 2. Include the following citation in any resulting publications: Gamm CM, Sullivan PF, Buchwal A, Dial RJ, Young AB, Watts DA, Cahoon SMP, Welker JM, Post E. 2018. Declining growth of deciduous shrubs in the warming climate of continental western Greenland. Journal of Ecology 106: 640-654.

kanger_dendro1 <- read.csv("data/kanger_rings_betula.csv") %>% mutate(Species = "B_nana")
kanger_dendro2 <- read.csv("data/kanger_rings_salix.csv") %>% mutate(Species = "S_glauca")

# Combine the Betula and Salix data
kanger_dendro <- rbind(kanger_dendro1, kanger_dendro2)
rm(kanger_dendro1, kanger_dendro2)

# Subset to same start as Qikiqtaruk
kanger_dendro_all <- kanger_dendro %>% filter(year > 1998) %>% group_by(year, Species) %>% 
  summarise(mean_rw = mean(ring.width_mm), rw_SD = sd(ring.width_mm), rw_SE = (sd(ring.width_mm)/sqrt(length(unique(disc)))))

kanger_shrub_detrend_salix <- kanger_dendro %>% filter(Species == "S_glauca") %>% group_by(disc) %>% mutate(detrend_rw = detrend.series(ring.width_mm)$Spline)

kanger_shrub_detrend_betula <- kanger_dendro %>% filter(Species == "B_nana") %>% group_by(disc) %>% mutate(detrend_rw = detrend.series(ring.width_mm)$Spline)

kanger_ndvi_detrend <- cbind.data.frame(kanger_greenup$year, detrend.series(kanger_greenup$ndvi.date.max)$Spline)
names(kanger_ndvi_detrend) <- c("year", "detrend_ndvi")

kanger_shrub_salix <- inner_join(kanger_shrub_detrend_salix, kanger_ndvi_detrend, by = c("year" = "year"))
kanger_shrub_betula <- inner_join(kanger_shrub_detrend_betula, kanger_ndvi_detrend, by = c("year" = "year"))

# Models fit using nlme
# kanger_model_salix <- lme(detrend_rw ~ detrend_ndvi, random=list(~detrend_ndvi|disc, ~1|year), cor=corAR1(), data = kanger_shrub_salix, method = "ML", control=list(maxIter=1000))

# kanger_model_betula <- lme(detrend_rw ~ detrend_ndvi, random= list(~detrend_ndvi|disc, ~1|year), cor=corAR1(), data = kanger_shrub_betula, method = "ML", control=list(maxIter=1000))

# Models fit using MCMCglmm
kanger_model_salix <- MCMCglmm(detrend_rw ~ detrend_ndvi, 
                      random = ~year + detrend_ndvi:disc, data = data.frame(kanger_shrub_salix), 
                      family = "gaussian", pr = TRUE, nitt = 100000, burnin = 20000, prior = prior2)

kanger_model_betula <- MCMCglmm(detrend_rw ~ detrend_ndvi, 
                      random = ~year + detrend_ndvi:disc, data = kanger_shrub_betula, 
                      family = "gaussian", pr = TRUE, nitt = 100000, burnin = 20000, prior = prior2)

summary(kanger_model_salix)

# Calculates pseudo R2 for nlme models
# r.squaredGLMM(kanger_model_salix)

summary(kanger_model_betula)

# Calculates pseudo R2 for nlme models
# r.squaredGLMM(kanger_model_betula)

fit_kanger_model_salix_pred <- as.data.frame(predict.MCMCglmm(kanger_model_salix, kanger_shrub_salix, interval = "prediction"))
fit_kanger_model_salix <- cbind.data.frame(year = kanger_shrub_salix$year, disc = kanger_shrub_salix$disc, detrend_rw = kanger_shrub_salix$detrend_rw, detrend_ndvi = kanger_shrub_salix$detrend_ndvi, fit = fit_kanger_model_salix_pred$fit, lwr = fit_kanger_model_salix_pred$lwr, upr = fit_kanger_model_salix_pred$upr) %>% arrange(desc(detrend_ndvi)) %>% filter(row_number() %% 10 == 0)

fit_kanger_model_betula_pred <- as.data.frame(predict.MCMCglmm(kanger_model_betula, kanger_shrub_betula, interval = "prediction"))
fit_kanger_model_betula <- cbind.data.frame(year = kanger_shrub_betula$year, disc = kanger_shrub_betula$disc, detrend_rw = kanger_shrub_betula$detrend_rw, detrend_ndvi = kanger_shrub_betula$detrend_ndvi, fit = fit_kanger_model_betula_pred$fit, lwr = fit_kanger_model_betula_pred$lwr, upr = fit_kanger_model_betula_pred$upr) %>% arrange(desc(detrend_ndvi)) %>% filter(row_number() %% 10 == 0)

# Mixed model figure
kanger_rw_ndvi <- ggplot() +
  geom_point(data = kanger_shrub_betula, aes(x=detrend_ndvi, y=detrend_rw), colour = "#EEB422", alpha = 0.5, size = 3) +
  geom_line(data = fit_kanger_model_betula, aes(x=detrend_ndvi, y=fit), colour = "#EEB422", size = 1) +
  geom_ribbon(data = fit_kanger_model_betula, aes(x=detrend_ndvi, ymin=lwr, ymax=upr), fill = "#EEB422", alpha = 0.5) +
  geom_point(data = kanger_shrub_salix, aes(x=detrend_ndvi, y=detrend_rw), colour = "#41A002", alpha = 0.5, size = 3) +
  geom_line(data = fit_kanger_model_salix, aes(x=detrend_ndvi, y=fit), colour = "#41A002", size = 1) +
  geom_ribbon(data = fit_kanger_model_salix, aes(x=detrend_ndvi, ymin=lwr, ymax=upr), fill = "#41A002", alpha = 0.5) +
  ylab("Detrended ring width\n") +
  xlab("\nDetrended NDVI") +
  scale_y_continuous(breaks = c(-1, 0, 1, 2, 3, 4), limits = c(-1, 4)) +
  scale_x_continuous(breaks = c(0.8, 0.9, 1.0, 1.1, 1.2), limits = c(0.8, 1.2)) +
  theme_QHI()

ggsave("plots/kanger_rw_ndvi.png", plot = kanger_rw_ndvi, width = 4, height = 4, units = "in", dpi = 300)

# Dendro plot
kanger_dendro <- kanger_dendro %>% filter(year > 1998) %>% group_by(year, Species) %>% 
   summarise(mean_rw = mean(ring.width_mm), rw_SD = sd(ring.width_mm), rw_SE = (sd(ring.width_mm)/sqrt(length(unique(disc)))))

kanger_dendro_plot <- ggplot(kanger_dendro, group = Species) + 
   geom_errorbar(aes(x=year, ymin = mean_rw-rw_SE, ymax=mean_rw+rw_SE, colour = Species), alpha = 0.5) +
   geom_path(aes(x=year, y=mean_rw, colour = Species), size = 1) +
   geom_point(aes(x=year, y=mean_rw, colour = Species), alpha = 0.5, size = 3) +
   #geom_smooth(method=lm, aes(x=year, y=mean_rw, colour = Species, fill = Species), alpha = 0.5, show.legend=F, linetype = "dotted") +
   ylab("Shrub ring width (mm)\n") +
   xlab("") +
   scale_y_continuous(breaks = c(0, 0.25, 0.5), limits = c(0, 0.5)) +
   scale_x_continuous(breaks = c(2000, 2005, 2010, 2015), limits = c(1998, 2017)) +
   scale_colour_manual(values = c("#EEB422", "#41A002"), labels = c("B. nana", "S. glauca")) +
   scale_fill_manual(values = c("#EEB422", "#41A002")) +
   theme_bw() +
   theme(axis.text = element_text(size = 16), 
         axis.title = element_text(size = 20),
         axis.text.x = element_text(angle = -45, hjust = -0.05),
         axis.line.x = element_line(color="black"), axis.line.y = element_line(color="black"),
         panel.border = element_blank(),
         panel.grid.major.x = element_blank(),                                          
         panel.grid.minor.x = element_blank(),
         panel.grid.minor.y = element_blank(),
         panel.grid.major.y = element_blank(),  
         plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
         plot.title = element_text(size=20, vjust=1, hjust=0),
         legend.text = element_text(size=14, face="italic"),          
         legend.title = element_blank(),                              
         legend.position = c(0.75, 0.92), 
         legend.key = element_blank(),
         legend.background = element_rect(color = "black", fill = "transparent", size = 4, linetype="blank"))

ggsave("plots/kanger_dendro.png", plot = kanger_dendro_plot, width = 4, height = 4, units = "in", dpi = 300)


# Correlations between rings and MODIS ----

# Combine QHI and Kanger info
greenup_sites <- rbind(kanger_greenup, qikiq_greenup) %>% dplyr::select(Site = site, Year = year, ndvi.date.max)

QHI_mean_cover2 <- QHI_mean_cover %>% mutate(Species = "S_pulchra") %>% rename(Year = year, abund_mean = cover_mean, abund_SE = cover_SE) %>% dplyr::select(Year, Species, abund_mean, abund_SE) %>% mutate(Site = "QHI")
kanger_shrub_abund2 <- kanger_shrub_abund %>% mutate(Site = "Kanger")

cover_sites <- rbind(QHI_mean_cover2, kanger_shrub_abund2)

QHI_shrub_dendro2 <- QHI_dendro %>% mutate(Site = "QHI", Species = "S_pulchra") %>% dplyr::select(Site, Year, mean_rw, rw_SD, Species) %>% as.data.frame()
kanger_dendro2 <- kanger_dendro %>% rename(Year = year) %>% mutate(Site = "Kanger") %>% dplyr::select(Site, Year, mean_rw, rw_SD, Species) %>% as.data.frame()

dendro_sites <- rbind(QHI_shrub_dendro2, kanger_dendro2)

# Create df of rings, cover and MODIS
correlations <- merge(cover_sites, dendro_sites, all = TRUE)
correlations <- merge(correlations, greenup_sites, all = TRUE)

corr <- correlations %>% group_by(Site, Species) %>% mutate(detrend_ndvi = detrend.series(ndvi.date.max)$Spline, detrend_rw = detrend.series(mean_rw)$Spline) %>%
   summarise(
      corr_rw_ndvi = cor(detrend_rw, detrend_ndvi, method = "pearson", use = "pairwise.complete.obs")
   ) %>% as.data.frame()

# Plot for Qikiqtaruk
ggplot(filter(correlations, Site == "QHI")) +
   geom_point(aes(x = mean_rw, y = ndvi.date.max), colour = "#41A002", alpha = 0.5, size = 2) +
   annotate("text", x = 0.25, y = 0.7, label = format(paste("Pearson's R = ", round(corr$corr_rw_ndvi[3],digits = 3)))) +
   theme_QHI()

# Plot for Kanger
ggplot(filter(correlations, Site == "Kanger")) +
   geom_point(aes(x = mean_rw, y = ndvi.date.max, colour = Species), alpha = 0.5, size = 2) +
   annotate("text", x = 0.25, y = 0.7, label = format(paste("Pearson's R = ", round(corr$corr_rw_ndvi[1],digits = 3)))) +
   annotate("text", x = 0.3, y = 0.65, label = format(paste(round(corr$corr_rw_ndvi[2],digits = 3)))) +
   theme_QHI()
