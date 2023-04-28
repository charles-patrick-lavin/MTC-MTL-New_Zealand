# Script for:

#Distinguishing the effects of fisheries and a warming climate on fish populations in Aotearoa, New Zealand

#Charles P. Lavin, Faculty of Biosciences and Aquaculture, Nord University, Bodø, Norway.
#ORCID: 0000-0002-2068-1080

#Daniel Pauly, Sea Around Us, Institute for the Ocean and Fisheries, University of British Columbia, Vancouver, British Columbia, Canada.
#ORCID: 0000-0003-3756-4793

#Donna Dimarchopoulou, Department of Biology, Dalhousie University, Halifax, Nova Scotia, Canada; Department of Biology, Woods Hole Oceanographic #Institution, Woods Hole, Massachusetts, USA.
#ORCID: 0000-0003-3412-3503

#Cui Liang, Key Laboratory of Marine Ecology and Environmental Science, Institute of Oceanology, Chinese Academy of Sciences, Qingdao, China
#ORCID: 0000-0001-6099-4965

#Mark Costello, Faculty of Biosciences and Aquaculture, Nord University, Bodø, Norway.
#ORCID: 0000-0003-2362-0328


# Corresponding author: Charles Lavin, charles.p.lavin@nord.no

# Loading required packages
require(tidyverse)
require(dplyr)
require(plyr)
require(data.table)
require(car)
require(tibble)
require(caTools)
require(rfishbase)
require(worrms)
require(corrplot)
require(ggplot2)
require(ggpubr)
require(ggpmisc)
require(segmented)
require(FactoMineR)
require(factoextra)

# Set working directory to file location
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Loading SSTA data for the New Zealand EEZ
# Extracted from https://psl.noaa.gov/data/gridded/data.kaplan_sst.html
ssta <- read.csv("NZ_EEZ_SSTA_1950_2020.csv")

Fig_3a <- ggplot(data=ssta, aes(x=Year, y=SSTA)) +
  geom_point() + geom_line(cex=0.75) +
  theme_classic() +
  geom_smooth(method="lm", col="red") +
  scale_x_continuous(breaks=seq(1950,2020, by=10)) +
  scale_y_continuous(breaks=seq(-0.75,0.75, by=0.25), limits = c(-0.75,0.75)) +
  labs(y = "SSTA (°C)", x= "Year") +
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=19),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=18)) +
  theme(axis.ticks.length=unit(.25, "cm"))  +
  annotate("text", x = 1950, y = 0.75, size = 7,
           label = "(a)", fontface = 1)

# SSTA linear regression
NZ_SSTA <- lm(SSTA ~ Year, data=ssta)
summary(NZ_SSTA)

# SSTA segmented regression
ssta.seg <- segmented(NZ_SSTA,
                      seg.Z = ~ Year,
                      psi = list(Year = c(2009)))

summary.segmented(ssta.seg)


# Reading in Sea Around Us fisheries data 
# https://www.seaaroundus.org/
EEZ_1_1 <- read.csv("SAU EEZ 554 v50-1.csv") # Mainland New Zealand EEZ, chunk 1
EEZ_1_2 <- read.csv("SAU EEZ 554 v50-2.csv") # Mainland New Zealand EEZ, chunk 2
EEZ_1_3 <- read.csv("SAU EEZ 554 v50-3.csv") # Mainland New Zealand EEZ, chunk 3

EEZ_1 <- rbind(EEZ_1_1, EEZ_1_2, EEZ_1_3)

EEZ_2 <- read.csv("SAU EEZ 555 v50-0.csv") # Kermadec Island EEZ

SAU <- rbind(EEZ_1[1:218957, 2:18], EEZ_2)



nrow(SAU) #227534

final <- SAU

final$year <- as.factor( paste0( as.factor( final$year )))
final$scientific_name <- as.factor( paste0( as.factor( final$scientific_name )))
SAU_final_spp <- as.data.frame(unique(final$scientific_name))


# Extracting species trait info from FishBase
SAU_FB <- rfishbase::estimate(species_list = unique(final$scientific_name))

# selecting spp, trophic level, growth parameter K,
# min depth, max depth, mean temp
SAU_traits <- SAU_FB[,(c(2,5,43))] 
SAU_traits$Species <- as.factor( paste0( as.factor( SAU_traits$Species )))
str(SAU_traits)
SAU_traits <- SAU_traits %>% drop_na(TempPrefMean) %>% distinct()
SAU_traits$Species <- as.factor( paste0( as.factor( SAU_traits$Species )))
str(SAU_traits) #110 species
names(SAU_traits) <- c("scientific_name", "trophic_level", "mean_temp")
# Note: FishBase species info for Centroselachus crepidater
# was added to the database during production of manuscript.
# This species was not included in the original analysis,
# and therefore, has been removed from the species list
# in the present script.

# removing Centroselachus crepidater
SAU_traits <- SAU_traits %>%
  filter(!grepl("Centroselachus crepidater", scientific_name))
list(unique(SAU_traits$scientific_name)) #109 species

# Catch by species / groups
SAU_catch <- final %>%
  dplyr::group_by(year, scientific_name) %>%
  dplyr::summarise(Annual_catch = sum(tonnes))
names(SAU_catch) <- c("Year", "scientific_name", "catch")

# Catch by gear type
SAU_catch_gear <- final %>%
  dplyr::group_by(year, scientific_name, gear_type) %>%
  dplyr::summarise(Annual_catch = sum(tonnes))


#Calculating MTC
SAU_mtc <- inner_join(SAU_catch, SAU_traits, by="scientific_name")
summary(SAU_mtc) #no NAs
str(SAU_mtc)
list(unique(SAU_mtc$scientific_name)) #110 species

#check temp outliers
boxplot(SAU_mtc$mean_temp)
sort( unique( SAU_mtc$mean_temp ) )

#calculating MTC
SAU_mtc_year <- SAU_mtc %>%
  dplyr::group_by(Year) %>%
  dplyr::summarise(MTC = sum(catch*mean_temp) / sum(catch))
str(SAU_mtc_year)
SAU_mtc_year$Year <- as.numeric(as.character(SAU_mtc_year$Year))

# Linear regression, MTC on Year
SAU_MTC_year <- lm(MTC ~ Year, data=SAU_mtc_year)
summary(SAU_MTC_year)

# Pearson's correlation, MTC with SSTA trend
cor.test(SAU_mtc_year$MTC, ssta[1:70,2]) # no significant correlation

# Segmented regression of MTC on Year
sau.seg.mtc <- segmented(SAU_MTC_year,
                         seg.Z = ~ Year,
                         psi = list(Year = c(1983)))

summary.segmented(sau.seg.mtc)

# get breakpoints
sau.seg.mtc$psi #1983
# get the slopes
slope(sau.seg.mtc)

# get the fitted data
seg.fitted.mtc <- fitted(sau.seg.mtc)
seg.model.mtc <- data.frame(Year = SAU_mtc_year$Year, MTC = seg.fitted.mtc)

# MTC
Fig_3b <- ggplot(data=SAU_mtc_year, aes(x=Year, y=MTC)) + geom_point() +
  geom_line(cex=0.75) +
  geom_line(data=seg.model.mtc, aes(x=Year, y=MTC)) +
  theme_classic() +
  scale_x_continuous(breaks=seq(1950,2020, by=10)) +
  scale_y_continuous(breaks=seq(8,18, by=2), limits=c(8,18)) +
  labs(y = "Mean Temp. Catch (°C)", x= "Year") +
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=19),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=18)) +
  theme(plot.margin = margin(0, 1, 0, 0.5,"cm")) +
  theme(axis.ticks.length=unit(.25, "cm")) +
  annotate("text", x = 1957, y = 18, size = 7,
           label = "(b)", fontface = 1)


# Figure 2, catch and # of groups / species

annual_catch <- SAU_catch %>% 
  dplyr::group_by(Year) %>%
  dplyr::summarise(Annual_catch = sum(catch)) %>%
  mutate(fig_catch = Annual_catch*0.000001)

str(annual_catch)
annual_catch$Year <- as.numeric(as.character(annual_catch$Year))

# # of species - includes higher taxonomic entries
annual_group <- SAU_catch %>% 
  dplyr::group_by(Year) %>%
  dplyr::summarise(spp = n_distinct(scientific_name))

annual_group$Year <- as.numeric(as.character(annual_group$Year))

# removing entries not at species level
tax <- SAU_catch %>% filter(grepl(" ", scientific_name)) %>%
  filter(!grepl("Marine", scientific_name)) %>%
  filter(!grepl("Miscellaneous ", scientific_name))

# which groups were removed?
tax_remove <- SAU_catch %>% filter(!grepl(" ", scientific_name))
marine_remove <- SAU_catch %>% filter(grepl("Marine", scientific_name))
misc_remove <- SAU_catch %>% filter(grepl("Miscellaneous", scientific_name))

taxa <- as.data.frame(unique(tax_remove$scientific_name))
names(taxa) <- c("group")
marine <- as.data.frame(unique(marine_remove$scientific_name))
names(marine) <- c("group")
misc <- as.data.frame(unique(misc_remove$scientific_name))
names(misc) <- c("group")

# which species were removed?
full_spp <- as.data.frame(unique(tax$scientific_name))
names(full_spp) <- ("scientific_name")
included_spp <- as.data.frame(unique(SAU_mtc$scientific_name))
names(included_spp) <- ("scientific_name")

spp_removed <- as.data.frame(setdiff(full_spp$scientific_name, included_spp$scientific_name))
names(spp_removed) <- c("group")


final_groups_removed <- rbind(taxa,marine,misc, spp_removed)


# WoRMS classification of removed groups
worms_class <- wm_classification_(name = c("Anguilliformes",
                                           "Bivalvia",
                                           "Bramidae",
                                           "Carangidae",
                                           "Carcharhinidae",
                                           "Carcharhiniformes",
                                           "Centrolophidae",
                                           "Chimaeriformes",
                                           "Echinoidea",
                                           "Elasmobranchii",
                                           "Gempylidae",
                                           "Haliotis",
                                           "Lamnidae",
                                           "Latridae",
                                           "Macrouridae",
                                           "Merlucciidae",
                                           "Moridae",
                                           "Mytilidae",
                                           "Nemadactylus",
                                           "Nototheniidae",
                                           "Pectinidae",
                                           "Perciformes",
                                           "Pinguipedidae",
                                           "Pleuronectiformes",
                                           "Polyprion",
                                           "Rajidae",
                                           "Saccostrea",
                                           "Scombridae",
                                           "Scorpaeniformes",
                                           "Scyliorhinidae",
                                           "Squalidae",
                                           "Squaliformes",
                                           "Trachurus",
                                           "Triakidae",
                                           "Dosinia",
                                           "Alopias",
                                           "Berycidae",
                                           "Bothidae",
                                           "Callorhinchidae",
                                           "Centriscidae",
                                           "Gadidae",
                                           "Istiophoridae",
                                           "Isurus",
                                           "Labridae",
                                           "Mullidae",
                                           "Myctophidae",
                                           "Nephropidae",
                                           "Ommastrephidae",
                                           "Oreosomatidae",
                                           "Pleuronectoidei",
                                           "Scyphozoa",
                                           "Serranidae",
                                           "Teuthida",
                                           "Trachichthyidae",
                                           "Triglidae",
                                           "Xiphiidae",
                                           "Zeidae",
                                           "Clupeidae",
                                           "Monacanthidae",
                                           "Sparidae",
                                           "Coryphaena",
                                           "Sphyrnidae",
                                           "Hydrolagus",
                                           "Decapoda",
                                           "Lamniformes",
                                           "Mugilidae",
                                           "Palinuridae",
                                           "Portunidae",
                                           "Veneridae",
                                           "Beryx",
                                           "Gadiformes",
                                           "Rajiformes",
                                           "Seriola",
                                           "Seriolella",
                                           "Lithodidae",
                                           "Mollusca",
                                           "Scorpaenidae",
                                           "Epinephelus",
                                           "Decapterus",
                                           "Myliobatidae",
                                           "Pterygotrigla",
                                           "Batoidea",
                                           "Dendrobranchiata",
                                           "Gastropoda",
                                           "Solenoceridae",
                                           "Chondrichthyes",
                                           "Clupeiformes",
                                           "Etmopterus",
                                           "Sphyraena",
                                           "Epigonus",
                                           "Benthodesmus",
                                           "Lampris",
                                           "Squalus",
                                           "Kyphosidae",
                                           "Tetraodontidae",
                                           "Sphyrna",
                                           "Trichiurus"))

names(worms_class) <- c("id", "AphialID", "rank", "scientific_name")
worms_class <- worms_class[,3:4]
taxa_rank <- left_join(taxa, worms_class, by=c("group"="scientific_name"))

taxa_rank <- taxa_rank %>% unique()



# Figure 2
annual_spp_tax <- tax %>% 
  dplyr::group_by(Year) %>%
  dplyr::summarise(spp = n_distinct(scientific_name))

annual_catch <- cbind(annual_catch, annual_spp_tax[,2])

annual_catch$Year <- as.numeric(as.character(annual_catch$Year))

Fig_2a = ggplot() +
  geom_area(data=annual_catch, aes(x=Year, y=fig_catch), fill="darkgrey") +
  theme_classic() +
  geom_path(data=annual_group, aes(x=Year, y=spp*0.005767442), cex=1, col="red1") +
  geom_path(data=annual_catch, aes(x=Year, y=spp*0.005767442), cex=1, col="black") +
  scale_x_continuous(breaks = seq(1950, 2020, by = 10)) +
  scale_y_continuous(breaks = seq(0, 1.2, by = 0.2),
                     sec.axis = sec_axis(~./0.005767442)) +
  theme(axis.text = element_text(colour = "black")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.y = element_text(size=19)) +
  theme(axis.text.x = element_text(size=19)) +
  theme(axis.ticks.length=unit(.25, "cm")) +
  annotate("text", x = 1950, y = 1.2, size = 6,
           label = "(a)", fontface = 1) +
  annotate("text", x = 1969.5, y = 1.2, size = 6,
           label = "Full fisheries catch dataset") +
  annotate("text", x = 1960, y = 0.725, size = 5.5,
           label = "Groups & species", col = "red1") +
  annotate("text", x = 1960, y = 0.325, size = 5.5,
           label = "Species only")

# Species and catch included in MTC and MTL analyses
mtc_catch <- SAU_mtc %>%
  dplyr::group_by(Year) %>%
  dplyr::summarise(Annual_catch = sum(catch)) %>%
  mutate(fig_catch = Annual_catch*0.000001)

mtc_spp <- SAU_mtc %>%
  dplyr::group_by(Year) %>%
  dplyr::summarise(spp = n_distinct(scientific_name))

mtc_catch <- cbind(mtc_catch, mtc_spp[,2])

str(mtc_catch)
mtc_catch$Year <- as.numeric(as.character(mtc_catch$Year))

Fig_2b <- ggplot() +
  geom_area(data=mtc_catch, aes(x=Year, y=fig_catch), fill="darkgrey") +
  theme_classic() +
  geom_path(data=mtc_catch, aes(x=Year, y=spp*0.005767442), cex=1, col="black") +
  scale_x_continuous(breaks = seq(1950, 2020, by = 10)) +
  scale_y_continuous(breaks = seq(0, 1.2, by = 0.2),
                     sec.axis = sec_axis(~./0.005767442)) +
  theme(axis.text = element_text(colour = "black")) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.y = element_text(size=19)) +
  theme(axis.text.x = element_text(size=19)) +
  theme(axis.ticks.length=unit(.25, "cm")) +
  annotate("text", x = 1950, y = 1.2, size = 6,
           label = "(b)", fontface = 1) +
  annotate("text", x = 1968.25, y = 1.2, size = 6,
           label = "Fisheries catch analysed") +
  annotate("text", x = 1960, y = 0.3, size = 5.5,
           label = "Species")

Fig_2 <- ggarrange(Fig_2a, Fig_2b,
                   ncol=1, nrow=2,
                   align = c("hv"))

Fig_2 = annotate_figure(Fig_2,
                        left = text_grob("Catch"~(tons~x~10^6), size=26,  rot = 90),
                        right = text_grob("Number of species or groups recorded in catch", size = 26, rot=-90)) +
  theme(plot.margin = margin(0, 0.5, 0, 0.5,"cm"))

#ggsave("Figure_2.tiff", Fig_2,
#       width = 22,
#       height = 28,
#       units = c("cm"),
#       device = c("tiff"))

# Percentage of total SAU catch (from raw data) included
# in MTC analyses (n=109 spp)
(sum(mtc_catch$Annual_catch))/(sum(annual_catch$Annual_catch)) # = 50%


# Now we account for geographical expansion of fisheries offshore and to deeper waters

# Following methods from:
# Liang C, Pauly D (2017) Fisheries impacts
# on China's coastal ecosystems: Unmasking a
# pervasive `fishing down' effect. PLoS ONE 12(3):
# e0173296. doi:10.1371/journal.pone.0173296

# 1 Calculating marine trophic index (MTI) ##################
sau_trophic <- SAU_FB[,(c(2,5))] %>% drop_na()  #selecting spp and Troph

sau_trophic$Species <- as.factor( paste0( as.factor( sau_trophic$Species )))

sau_trophic <- sau_trophic %>% drop_na() %>% distinct()
str(sau_trophic) # 112 spp
names(sau_trophic) <- c("scientific_name", "trophic_level")

# Two species missing temp value, but contain trophic level value from FB
SAU_traits$scientific_name <- as.factor( paste0( as.factor( SAU_traits$scientific_name )))
sau_trophic$scientific_name <- as.factor( paste0( as.factor( sau_trophic$scientific_name )))

setdiff(sau_trophic$scientific_name, SAU_traits$scientific_name)
# "Odax pullus", "Hyporhamphus ihi" (also "Centroselachus crepidater")

# remove these two species from trophic level analyses
nrow(sau_trophic)
sau_trophic <- droplevels(sau_trophic[!sau_trophic$scientific_name == 'Odax pullus',])
sau_trophic <- droplevels(sau_trophic[!sau_trophic$scientific_name == 'Hyporhamphus ihi',])
sau_trophic <- droplevels(sau_trophic[!sau_trophic$scientific_name == 'Centroselachus crepidater',])
nrow(sau_trophic) #109 spp

# filtering SAU data to include species with TL info from FishBase
final_n <- final %>%
  filter(scientific_name %in% sau_trophic$scientific_name)

NZ_catch <- final_n %>%
  dplyr::group_by(year, scientific_name) %>%
  dplyr::summarise(annual_catch = sum(tonnes))

list(unique(NZ_catch$scientific_name))

# Now assign trophic level to each spp

main_MTI <- inner_join(NZ_catch, sau_trophic, by="scientific_name")
summary(main_MTI) #no NAs
str(main_MTI)
list(unique(main_MTI$scientific_name))

# Calculating MTI
NZ_MTI <- main_MTI %>%
  dplyr::group_by(year) %>%
  dplyr::mutate(MTI = sum(annual_catch*trophic_level)/sum(annual_catch)) %>%
  dplyr::mutate(Catch = sum(annual_catch)) %>%
  distinct(year, MTI, Catch)
names(NZ_MTI) <- c("Year", "MTI", "Catch")

NZ_MTI$Year <- as.numeric(as.character(NZ_MTI$Year))
NZ_MTI$MTI <- as.numeric(as.character(NZ_MTI$MTI))
NZ_MTI$Catch <- as.numeric(as.character(NZ_MTI$Catch))

# Linear regression, MTI on Year
SAU_MTI_year <- lm(MTI ~ Year, data=NZ_MTI)
summary(SAU_MTI_year)

fitted.mti.full <- fitted(SAU_MTI_year)
model.mti.full <- data.frame(Year = NZ_MTI$Year, MTI = fitted.mti.full)


# Segmented regression of full MTI and Year
sau.seg.mti <- segmented(SAU_MTI_year,
                         seg.Z = ~ Year,
                         psi = list(Year = c(1995)))

summary.segmented(sau.seg.mti)

# get breakpoints
sau.seg.mti$psi #1995
# get the slopes
slope(sau.seg.mti)

# get the fitted data
seg.fitted.mti <- fitted(sau.seg.mti)
seg.model.mti <- data.frame(Year = NZ_MTI$Year, MTI = seg.fitted.mti)

# MTL
Fig_3c <- ggplot(data=NZ_MTI, aes(x=Year, y=MTI)) + geom_point() +
  geom_line(cex=0.75) +
  geom_line(data=seg.model.mti, aes(x=Year, y=MTI)) +
  theme_classic() +
  scale_x_continuous(breaks=seq(1950,2020, by=10)) +
  scale_y_continuous(breaks=seq(3.6,4.2, by=0.2), limits = c(3.5,4.3)) +
  labs(y = "Mean Trophic Level", x= "Year") +
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=19),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=18)) +
  theme(plot.margin = margin(0, 1, 0, 0.5,"cm")) +
  theme(axis.ticks.length=unit(.25, "cm")) +
  annotate("text", x = 1957, y = 4.3, size = 7,
           label = "(c)", fontface = 1)


# 2 Calculating fishing in balance index (FiB)
FIB_NZ <- as.data.frame(log10(NZ_MTI$Catch*((1/0.1)^NZ_MTI$MTI)) - log10(27176.18*((1/0.1)^3.705097)))

NZ_MTI_FiB <- cbind(NZ_MTI, FIB_NZ) %>%
  mutate(EEZ = "NZ_location")

names(NZ_MTI_FiB) = c("Year", "MTL", "Catch", "FiB", "EEZ")
NZ_MTI_FiB <- NZ_MTI_FiB[,c(1,5,3,2,4)] # re-ordering

# Linear regression of FiB
SAU_FiB_year <- lm(FiB ~ Year, data=NZ_MTI_FiB)
summary(SAU_FiB_year)

fitted.FiB.full <- fitted(SAU_FiB_year)
model.FiB.full <- data.frame(Year = NZ_MTI_FiB$Year, FiB = fitted.FiB.full)

# Segmented regression of  FiB on Year
sau.seg.FiB <- segmented(SAU_FiB_year, #SAU_FiB_year
                         seg.Z = ~ Year,
                         psi = list(Year = c(1990)))

summary.segmented(sau.seg.FiB)

# get breakpoints
sau.seg.FiB$psi #1990
# get the slopes
slope(sau.seg.FiB)

# get the fitted data
seg.fitted.FiB <- fitted(sau.seg.FiB)
seg.model.FiB <- data.frame(Year = NZ_MTI_FiB$Year, FiB = seg.fitted.FiB)

#FiB
Fig_3d <- ggplot(data=NZ_MTI_FiB, aes(x=Year, y=FiB)) +
  geom_point() +
  geom_line(cex=0.75) +
  theme_classic() +
  geom_line(data=seg.model.FiB, aes(x=Year, y=FiB)) +
  scale_x_continuous(breaks=seq(1950,2020, by=10)) +
  labs(y = "FiB Index", x="Year") +
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=19),
        axis.title.x = element_text(size=24),
        axis.title.y = element_text(size=18)) +
  theme(plot.margin = margin(0, 1, 0, 0.5,"cm")) +
  theme(axis.ticks.length=unit(.25, "cm")) +
  annotate("text", x = 1957, y = 2, size = 7,
           label = "(d)", fontface = 1)

#######################################################################
# Utilising RMTI tool from:

# http://www.seaaroundus.org/regional-mti-tools/

# Based on the techniques defined in the following publication:

# Kleisner K, Mansour H, Pauly D (2014) Region-based MTI: resolving geographic
# expansion in the Marine Trophic Index. Mar Ecol Prog Ser 512:185-199

# Load input MTL and FiB data
EEZ <- NZ_MTI_FiB

# Display input file organization
head(EEZ)

# Display dimension of input data
dim(EEZ)

# Specify output file
newpdf<-pdf(file="RMTI_NZ_output.pdf", paper="letter", onefile=T)

# Split input data by EEZ
EEZ1split<-   split( EEZ,EEZ$EEZ)


#################################################
####          Define RMTI function           ####
#################################################
RMTI<-function(EEZ, EEZ1split, newpdf){
  
  par(mfrow=c(2,2), mex=.5, mai=c( 0.4, 0.4, 0.4, 0.1) )
  
  ExpandYears = NULL
  for(j in names(EEZ1split) ) {
    
    EEZ1<-data.frame(EEZ1split[j]  )
    names(EEZ1) <-names(EEZ)
    EEZ1<-EEZ1[ order( EEZ1[,1]),]
    
    CatchSeries<-t(t(EEZ1[,3]))
    CSS<-runmean(CatchSeries, 5, endrule="mean")
    CSSdf<-data.frame(CSS, row.names=NULL)
    EEZ1<- cbind(EEZ1,CSSdf)
    
    MTLSeries<-t(t(EEZ1[,4]))
    MTLS<-runmean(MTLSeries, 5, endrule="mean")
    MTLSdf<-data.frame(MTLS, row.names=NULL)
    EEZ1<- cbind(EEZ1,MTLSdf)
    
    
    pYk=EEZ1$Catch[1]
    
    # loop to calculate potential catch for region 1
    for( i in 1:length(EEZ1$MTL)){
      # Define TL to calculate pYk ranging from [2.5, 4.5]
      Tlevel = seq(max(2.5,min(EEZ1$MTL)-0.5),min(max(EEZ1$MTL)+0.5,4.5),0.1)
      # potential catch per TL
      pYkl = EEZ1$CSS[1]*10^(-EEZ1$MTLS[i] + Tlevel)
      # potential catch per year based on potential catch per TL
      pYk[i] = mean(pYkl)
    }
    y0 = 1
    
    # Save the start year of each region in Y0
    Y0 = EEZ1$Year[y0]
    Country = EEZ1$EEZ[1]
    
    # find year where catch exceeds potential catch
    yn = which(EEZ1$CSS - pYk >= 0*mean(pYk))
    
    ###################################################################
    # start plotting data
    frmbnd = c(min(min(EEZ1$CSS),min(pYk)),max(max(EEZ1$CSS),max(pYk)))
    plot(EEZ1$Year, EEZ1$CSS, ylim=frmbnd, xlab=" ", ylab="Smoothed Catch",  type="l", cex.lab=1.1)
    #lines(EEZ1$Year, pYk, type="l",col=2)
    plot(EEZ1$Year, EEZ1$FIB, xlab=" ", ylab="FiB",  type="l", cex.lab=1.1)
    plot(EEZ1$Year, EEZ1$MTLS,  ylim=c(2,5),xlab=" ", ylab="MTI", type="l", cex.lab=1.1)
    
    
    mtext(as.character(EEZ1$EEZ[1]),outer=T, side=3, line=-4, cex=2)
    mtext("Year", outer=T, side=1, line=-2)
    
    # end of plotting
    ############################################################################
    
    #### split regions ####
    EEZ2<-EEZ1
    N<-length(EEZ1$MTL)
    i=7
    while(1){
      i = i+1
      
      vn = (CSS - pYk >= 0.05*mean(pYk))
      vvn = ifelse(vn=="TRUE", 1, 0)
      zn = as.numeric(runmed(vvn, 5, endrule="median"))
      
      yn = which(diff(zn) > 0)
      
      if (yn[1] >= y0+4 && length(yn)>0){
        if(length(yn)>1){
          yn = yn[1:length(yn)]
        }}
      else{
        yn = NULL
      }
      
      
      if(length(yn) == 0 || i == 10 ){
        EEZ2$cMTL<-MTLS
        
        # append newly detected MTL as column to EEZ1 data
        EEZ1<-cbind(EEZ1,EEZ2$cMTL)
        if(i==8){
          plot(EEZ1$Year, EEZ1[,i], ylim=c(2,5), xlab=" ", ylab="cMTI", type="l", cex.lab=1.1)
        }
        else{
          lines(EEZ1$Year, EEZ1[,i], type="l",col=i-6)
        }
        break
      }
      else{
        # initialize Catch and MTL based on the new year 0
        Catcho = CSS[yn[1]-1]
        MTLo = MTLS[yn[1]-1]
        
        # calculate catch for old region in new time period
        EEZ2$pCatchk<- (10^(MTLo-MTLS))*Catcho
        EEZ2$pCatchk[1:yn[1]-1] = CSS[1:yn[1]-1]
        
        # calculate MTI for old region in new time period
        EEZ2$MTL_stat<-MTLo-log10(CSS/Catcho)
        # set the corresponding catch in the first region to zero
        EEZ2$pCatchk[which(EEZ2$MTL_stat <= 2)] = 0
        EEZ2$MTL_stat[which(EEZ2$MTL_stat <= 2)] = 2
        EEZ2$MTL_stat[which(EEZ2$MTL_stat >= 4.5)] = 4.5
        EEZ2$MTL_stat[1:yn[1]-1] = MTLS[1:yn[1]-1]
        
        # calculate catch and MTI of new region in new time period
        EEZ2$Catchk_new <- CSS - EEZ2$pCatchk
        EEZ2$Catchk_new[which(EEZ2$Catchk_new < 0)] = 0
        EEZ2$MTL_new <- (MTLS*CSS - EEZ2$pCatchk*EEZ2$MTL_stat)/(EEZ2$Catchk_new)
        
        # save corrected MTI in cMTL and attach to EEZ1
        EEZ2$cMTL<-EEZ2$MTL_stat
        EEZ2$cMTL[1:yn[1]-1] = MTLS[1:yn[1]-1]
        EEZ1<-cbind(EEZ1,EEZ2$cMTL)
        
        # Plot new region based MTI
        if(i==8){
          plot(EEZ1$Year, EEZ1[,i], ylim=c(2,5), xlab=" ", ylab="cMTI", type="l", cex.lab=1.1,col=2)
        }
        else{
          lines(EEZ1$Year, EEZ1[,i], type="l",col=i-6)
        }
        
        # calculate MTLS and CSS for new region
        y0 = yn[1]
        Y0 = rbind(Y0, EEZ1$Year[y0])
        Country = rbind(Country, EEZ1$EEZ[1])
        
        EEZ2$MTL=EEZ2$MTL_new
        MTLS<-EEZ2$MTL
        EEZ2$MTL[1:y0-1]=NaN
        MTLS[1:y0-1]=NaN
        EEZ2$Catch=EEZ2$Catchk_new
        CSS<-EEZ2$Catch
        EEZ2$Catch[1:y0-1]=NaN
        CSS[1:y0-1]=0
        
        
      } # end if(length(yn) == 0)
      
    } # end while(1)
    
    
    ExpandYears = rbind(ExpandYears, cbind(Country,Y0));
    
  } # end for loop
  
  return(ExpandYears)
} # end function

# Call RMTI function
ExpandYears = RMTI(EEZ, EEZ1split, newpdf)

dev.off()

# loading regional MTI (adjusted Mean Trophic Level) values
aMTL <- read.csv("aMTL_values.csv")
names(aMTL) <- c("Year", "MTL1", "MTL2", "MTL3")

Fig_3e <- ggplot(data=aMTL, aes(x=Year, y=MTL1)) +
  geom_path(cex=0.75, col="chartreuse2") +
  annotate("text", x = 1957, y = 3.5, size = 4.25,
           label = "Stock") +
  annotate("text", x = 1957.5, y = 3.35, size = 4.25,
           label = "Assemblage 1") +
  geom_path(data= aMTL, aes(x=Year, y=MTL2), cex=0.75, col = "blue1") +
  annotate("text", x = 1955.5, y = 4.15, size = 4.25,
           label = "Stock") +
  annotate("text", x = 1955.5, y = 4.0, size = 4.25,
           label = "Assemblage") +
  annotate("text", x = 1955.5, y = 3.85, size = 4.25,
           label = "2") +
  geom_path(data= aMTL, aes(x=Year, y=MTL3), cex=0.75, col = "red1") +
  annotate("text", x = 1995, y = 4.07, size = 4.25,
           label = "Stock Assemblage 3") +
  theme_classic() +
  scale_x_continuous(breaks=seq(1950,2020, by=10)) +
  scale_y_continuous(breaks=seq(2.6, 4.2, 0.4), limits = c(2.6,4.5)) +
  labs(y = "aMTL", x="Year") +
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=19),
        axis.title.x = element_text(size=24),
        axis.title.y = element_text(size=18)) +
  geom_point(x=1965, y=3.66975, cex=2, shape=1, col="chartreuse2") +
  geom_point(x=1965, y=4.12074, cex=2, shape=1, col ="blue1") +
  geom_point(x=1969, y=3.88778, cex=2, shape=1, col="blue1") +
  geom_point(x=1969, y=4.37216, cex=2, shape=1, col="red1") +
  geom_segment(aes(x=1965, xend=1965, y=3.66975, yend=4.12074), linetype=2) +
  geom_segment(aes(x=1969, xend=1969, y=3.88778, yend=4.37216), linetype=2) +
  theme(plot.margin = margin(0, 1, 0, 0.5,"cm")) +
  theme(axis.ticks.length=unit(.25, "cm")) +
  annotate("text", x = 1957, y = 4.5, size = 7,
           label = "(e)", fontface = 1)

# Figure 3
Fig_3 <- ggarrange(Fig_3a,
                   ggarrange(Fig_3b, Fig_3c, Fig_3d, Fig_3e, ncol = 2, nrow=2, align = c("hv")),
                   nrow = 2, align = c("h"))

#ggsave("Figure_3.tiff", Fig_3,
#       width = 28,
#       height = 30,
#       units = c("cm"),
#       device = c("tiff"))



# Calculating MTC between species' milieu
SAU_spp_milieu <- read.csv("SAU_spp_milieu.csv") #adding spp milieu classification
# (taken from individual species' FishBase page)

SAU_milieu <- inner_join(SAU_mtc, SAU_spp_milieu, by="scientific_name")

SAU_milieu$milieu <- as.factor(SAU_milieu$milieu)

str(SAU_milieu)

MTC_milieu <- SAU_milieu %>%
  dplyr::group_by(Year, milieu) %>%
  dplyr::summarise(MTC = sum(catch*mean_temp)/sum(catch))

Milieu_mean_temp <- SAU_milieu %>%
  dplyr::group_by(milieu) %>%
  dplyr::summarise(sd = sd(mean_temp),
                   mean_temp = mean(mean_temp),
                   n=n(),
                   se=sd/sqrt(n))

MTC_milieu$Year <- as.numeric(as.character(MTC_milieu$Year))
MTC_milieu <- as.data.frame(MTC_milieu)

ggplot(data = MTC_milieu, aes(x=Year, y=MTC)) + geom_path() + 
  facet_wrap(~milieu) +
  stat_poly_line(col="black") +
  stat_poly_eq(use_label(c("eq", "R2"))) +
  scale_x_continuous(breaks=seq(1950,2020, by=10)) +
  theme_classic() +
  theme(axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=30),
        axis.title.y = element_text(size=30))


# Bathydemersal
bathydem.mtc <- lm(MTC ~ Year, data=MTC_milieu %>% filter(milieu == "bathydemersal"))
summary(bathydem.mtc)

# get the fitted data
bathydem.fitted.mtc <- fitted(bathydem.mtc)
bathydem.model.mtc <- data.frame(Year = seq(1950,2019, 1), MTC = bathydem.fitted.mtc)

bathydemersal <- ggplot(data = MTC_milieu %>% filter(milieu == "bathydemersal"), aes(x=Year, y=MTC)) +
  geom_point() +
  geom_line(cex=0.75) +
  geom_line(data=bathydem.model.mtc, aes(x=Year, y=MTC)) +
  #stat_poly_line(col="red") +
  #stat_poly_eq(use_label(c("eq", "R2")), size = 6, label.x = 0.5, label.y = 0.9) +
  scale_x_continuous(breaks=seq(1950,2020, by=10)) +
  scale_y_continuous(breaks=seq(5,25, by=5), limits = c(5,25)) +
  labs(y = "MTC (°C)", x= "Year") +
  theme_classic() +
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=19),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.length=unit(.25, "cm")) +
  annotate("text", x = 1975, y = 25, size = 8, label = "(a) Bathydemersal", col="black") +
  annotate("text", x = 2010, y = 22, size = 6, label = "n = 18", col="black") + 
  annotate("text", x = 1995, y = 19.75, size = 6, parse = T, label = "bar(x)", col="black") +
  annotate("text", x = 2010, y = 19.75, size = 6, label = "= 8.8 °C ± 0.7", col="black")


# Bathypelagic
bathypel.mtc <- lm(MTC ~ Year, data=MTC_milieu %>% filter(milieu == "bathypelagic"))

summary(bathypel.mtc)


bathypel.seg.mtc <- segmented(bathypel.mtc, 
                              seg.Z = ~ Year,
                              psi = list(Year = c(1980)))

summary.segmented(bathypel.seg.mtc)

# get breakpoints
bathypel.seg.mtc$psi # 1980
# get the slopes
slope(bathypel.seg.mtc)

# P-score for significant non-similar slope
pscore.test(bathypel.seg.mtc)

# get the fitted data
bathypel.fitted.mtc <- fitted(bathypel.seg.mtc)
bathypel.model.mtc <- data.frame(Year = seq(1975,2019, 1), MTC = bathypel.fitted.mtc)

bathypelagic <- ggplot(data = MTC_milieu %>% filter(milieu == "bathypelagic"), aes(x=Year, y=MTC)) +
  geom_point() +
  geom_line(cex=0.75) +
  geom_line(data=bathypel.model.mtc, aes(x=Year, y=MTC)) +
  #stat_poly_line(col="red") +
  #stat_poly_eq(use_label(c("eq", "R2")), size = 6, label.x = 0.5, label.y = 0.9) +
  scale_x_continuous(breaks=seq(1950,2020, by=10), limits = c(1950,2020)) +
  scale_y_continuous(breaks=seq(5,25, by=5), limits = c(5,25)) +
  labs(y = "MTC (°C)", x= "Year") +
  theme_classic() +
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=19),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.length=unit(.25, "cm"),
        plot.margin = margin(0, 0.5, 0, 0.5,"cm")) +
  annotate("text", x = 1975, y = 25, size = 8, label = "(b) Bathypelagic", col="black") +
  annotate("text", x = 2010, y = 22, size = 6, label = "n = 13", col="black") + 
  annotate("text", x = 1995, y = 19.75, size = 6, parse = T, label = "bar(x)", col="black") +
  annotate("text", x = 2010, y = 19.75, size = 6, label = "= 7.8 °C ± 1.0", col="black")


# Benthopelagic
benthopel.mtc <- lm(MTC ~ Year, data=MTC_milieu %>% filter(milieu == "benthopelagic"))

summary(benthopel.mtc)

benthopel.seg.mtc <- segmented(benthopel.mtc, 
                               seg.Z = ~ Year,
                               psi = list(Year = c(1977)))

summary.segmented(benthopel.seg.mtc)

# get breakpoints
benthopel.seg.mtc$psi #1977
# get the slopes
slope(benthopel.seg.mtc)

# P-score for significant non-similar slope
pscore.test(benthopel.seg.mtc)

# get the fitted data
benthopel.fitted.mtc <- fitted(benthopel.seg.mtc)
benthopel.model.mtc <- data.frame(Year = unique(MTC_milieu$Year), MTC = benthopel.fitted.mtc)

benthopelagic <- ggplot(data = MTC_milieu %>% filter(milieu == "benthopelagic"), aes(x=Year, y=MTC)) +
  geom_point() +
  geom_line(cex=0.75) +
  geom_line(data=benthopel.model.mtc, aes(x=Year, y=MTC)) +
  #stat_poly_line(col="blue") +
  #stat_poly_eq(use_label(c("eq", "R2")), size = 6, label.x = 0.5, label.y = 0.9) +
  scale_x_continuous(breaks=seq(1950,2020, by=10)) +
  scale_y_continuous(breaks=seq(5,25, by=5), limits = c(5,25)) +
  labs(y = "MTC (°C)", x= "Year") +
  theme_classic() +
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=19),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.length=unit(.25, "cm"),
        plot.margin = margin(0, 0.5, 0, 0.5,"cm")) +
  annotate("text", x = 1974, y = 25, size = 8, label = "(c) Benthopelagic", col="black") +
  annotate("text", x = 2010, y = 22, size = 6, label = "n = 26", col="black") + 
  annotate("text", x = 1992, y = 19.75, size = 6, parse = T, label = "bar(x)", col="black") +
  annotate("text", x = 2008, y = 19.75, size = 6, label = "= 13.6 °C ± 4.7", col="black")


# Demersal
dem.mtc <- lm(MTC ~ Year, data=MTC_milieu %>% filter(milieu == "demersal"))

summary(dem.mtc)

dem.seg.mtc <- segmented(dem.mtc, 
                         seg.Z = ~ Year,
                         psi = list(Year = c(1995)))

summary.segmented(dem.seg.mtc)

# get breakpoints
dem.seg.mtc$psi # 1995
# get the slopes
slope(dem.seg.mtc)

# P-score for significant non-similar slope
pscore.test(dem.seg.mtc)

# get the fitted data
dem.fitted.mtc <- fitted(dem.seg.mtc)
dem.model.mtc <- data.frame(Year = unique(MTC_milieu$Year), MTC = dem.fitted.mtc)

demersal <- ggplot(data = MTC_milieu %>% filter(milieu == "demersal"), aes(x=Year, y=MTC)) +
  geom_point() +
  geom_line(cex=0.75) +
  geom_line(data=dem.model.mtc, aes(x=Year, y=MTC)) +
  #stat_poly_line(col="blue") +
  #stat_poly_eq(use_label(c("eq", "R2")), size = 6, label.x = 0.5, label.y = 0.9) +
  scale_x_continuous(breaks=seq(1950,2020, by=10)) +
  scale_y_continuous(breaks=seq(5,25, by=5), limits = c(5,25)) +
  labs(y = "MTC (°C)", x= "Year") +
  theme_classic() +
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=19),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.length=unit(.25, "cm")) +
  annotate("text", x = 1969, y = 25, size = 8, label = "(d) Demersal", col="black") +
  annotate("text", x = 2010, y = 22, size = 6, label = "n = 16", col="black") + 
  annotate("text", x = 1992, y = 19.75, size = 6, parse = T, label = "bar(x)", col="black") +
  annotate("text", x = 2008, y = 19.75, size = 6, label = "= 14.0 °C ± 2.4", col="black")


# Pelagic-neritic
ner.mtc <- lm(MTC ~ Year, data=MTC_milieu %>% filter(milieu == "pelagic-neritic"))

summary(ner.mtc)

ner.seg.mtc <- segmented(ner.mtc, 
                         seg.Z = ~ Year,
                         psi = list(Year = c(1980)))

summary.segmented(ner.seg.mtc)

# get breakpoints
ner.seg.mtc$psi # 1980
# get the slopes
slope(ner.seg.mtc)

# P-score for significant non-similar slope
pscore.test(ner.seg.mtc)

# get the fitted data
ner.fitted.mtc <- fitted(ner.seg.mtc)
ner.model.mtc <- data.frame(Year = unique(MTC_milieu$Year), MTC = ner.fitted.mtc)


neritic <- ggplot(data = MTC_milieu %>% filter(milieu == "pelagic-neritic"), aes(x=Year, y=MTC)) +
  geom_point() +
  geom_line(cex=0.75) +
  geom_line(data=ner.model.mtc, aes(x=Year, y=MTC)) +
  #stat_poly_line(col="blue") +
  #stat_poly_eq(use_label(c("eq", "R2")), size = 6, label.x = 0.5, label.y = 0.9) +
  scale_x_continuous(breaks=seq(1950,2020, by=10)) +
  scale_y_continuous(breaks=seq(5,25, by=5), limits = c(5,25)) +
  labs(y = "MTC (°C)", x= "Year") +
  theme_classic() +
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=19),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.length=unit(.25, "cm"),
        plot.margin = margin(0, 0.5, 0, 0.5,"cm")) +
  annotate("text", x = 1975, y = 25, size = 8, label = "(e) Pelagic-neritic", col="black") +
  annotate("text", x = 2010, y = 22, size = 6, label = "n = 6", col="black") + 
  annotate("text", x = 1992, y = 19.75, size = 6, parse = T, label = "bar(x)", col="black") +
  annotate("text", x = 2008, y = 19.75, size = 6, label = "= 17.1 °C ± 6.3", col="black")


# Pelagic oceanic
ocean.mtc <- lm(MTC ~ Year, data=MTC_milieu %>% filter(milieu == "pelagic-oceanic"))

summary(ocean.mtc)

ocean.seg.mtc <- segmented(ocean.mtc, 
                           seg.Z = ~ Year,
                           psi = list(Year = c(1965)))

summary.segmented(ocean.seg.mtc)

# get breakpoints
ocean.seg.mtc$psi #1965
# get the slopes
slope(ocean.seg.mtc)

# get the fitted data
ocean.fitted.mtc <- fitted(ocean.seg.mtc)
ocean.model.mtc <- data.frame(Year = unique(MTC_milieu$Year), MTC = ocean.fitted.mtc)

oceanic <- ggplot(data = MTC_milieu %>% filter(milieu == "pelagic-oceanic"), aes(x=Year, y=MTC)) +
  geom_point() +
  geom_line(cex=0.75) +
  geom_line(data=ocean.model.mtc, aes(x=Year, y=MTC)) +
  #stat_poly_line(col="red") +
  #stat_poly_eq(use_label(c("eq", "R2")), size = 6, label.x = 0.5, label.y = 0.1) +
  scale_x_continuous(breaks=seq(1950,2020, by=10)) +
  scale_y_continuous(breaks=seq(5,25, by=5), limits = c(5,25)) +
  labs(y = "MTC (°C)", x= "Year") +
  theme_classic() +
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=19),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.length=unit(.25, "cm"),
        plot.margin = margin(0, 0.5, 0, 0.5,"cm")) +
  annotate("text", x = 1976, y = 25, size = 8, label = "(f) Pelagic-oceanic", col="black") +
  annotate("text", x = 2010, y = 16, size = 6, label = "n = 25", col="black") + 
  annotate("text", x = 1992, y = 13.75, size = 6, parse = T, label = "bar(x)", col="black") +
  annotate("text", x = 2008, y = 13.75, size = 6, label = "= 20.0 °C ± 7.4", col="black")


# Reef-associated
reef.mtc <- lm(MTC ~ Year, data=MTC_milieu %>% filter(milieu == "reef-associated"))

summary(reef.mtc)

# get the fitted data
reef.fitted.mtc <- fitted(reef.mtc)
reef.model.mtc <- data.frame(Year = seq(1950,2019, 1), MTC = reef.fitted.mtc)

reef <- ggplot(data = MTC_milieu %>% filter(milieu == "reef-associated"), aes(x=Year, y=MTC)) +
  geom_point() +
  geom_line(cex=0.75) +
  geom_line(data=reef.model.mtc, aes(x=Year, y=MTC)) +
  #stat_poly_line(col="red") +
  #stat_poly_eq(use_label(c("eq", "R2")), size = 6, label.x = 0.5, label.y = 0.9) +
  scale_x_continuous(breaks=seq(1950,2020, by=10)) +
  scale_y_continuous(breaks=seq(5,25, by=5), limits = c(5,25)) +
  labs(y = "MTC (°C)", x= "Year") +
  theme_classic() +
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=19),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.length=unit(.25, "cm")) +
  annotate("text", x = 1977, y = 25, size = 8, label = "(g) Reef-associated", col="black") +
  annotate("text", x = 2010, y = 22, size = 6, label = "n = 6", col="black") + 
  annotate("text", x = 1992, y = 19.75, size = 6, parse = T, label = "bar(x)", col="black") +
  annotate("text", x = 2008, y = 19.75, size = 6, label = "= 20.4 °C ± 4.5", col="black")



Figure_4 <- ggarrange(bathydemersal, bathypelagic, benthopelagic, demersal, neritic, oceanic, reef,
                      ncol = 2, nrow=4,
                      align = c("hv"))

Fig_4 = annotate_figure(Figure_4,
                        left = text_grob("MTC (°C)", size=26,  rot = 90),
                        bottom = text_grob("Year", size = 26)) +
  theme(plot.margin = margin(0, 0.5, 0, 0.5,"cm"))


#ggsave("Figure_4.tiff", Fig_4,
#       width = 28,
#       height = 36,
#       units = c("cm"),
#       device = c("tiff"))


#Calculating MTC by gear
SAU_mtc_gear <- inner_join(SAU_catch_gear, SAU_traits, by="scientific_name")
summary(SAU_mtc_gear) #no NAs
str(SAU_mtc_gear)
list(unique(SAU_mtc_gear$scientific_name)) #109
list(unique(SAU_mtc_gear$gear_type))

SAU_mtc_gear_class <- SAU_mtc_gear %>%
  mutate(gear_fishing = case_when(gear_type == "small scale gillnets" | gear_type == "small scale seine nets" | gear_type == "recreational fishing gear" | gear_type == "subsistence fishing gear" | gear_type == "small scale longline" | gear_type == "artisanal fishing gear" | gear_type == "small scale hand lines" | gear_type ==  "small scale pots or traps"| gear_type == "small scale purse seine" | gear_type == "small scale lines" | gear_type == "small scale other nets" ~ "Small scale",
                                  gear_type == "longline" ~ "Longline",
                                  gear_type == "unknown class" ~ "Unknown",
                                  gear_type == "bottom trawl" ~ "Bottom trawl",
                                  gear_type == "gillnet" ~ "Gillnet",
                                  gear_type == "hand lines" ~ "Hand lines",
                                  gear_type == "pots or traps" ~ "Pots or traps",
                                  gear_type == "purse seine" ~ "Purse seine",
                                  gear_type == "mixed gear" ~ "Mixed gear",
                                  gear_type == "other" | gear_type == "other nets" ~ "Other",
                                  gear_type == "pelagic trawl" ~ "Pelagic trawl",
                                  gear_type == "pole and line" ~ "Pole and line"))

SAU_mtc_gear_class$gear_fishing <- as.factor(SAU_mtc_gear_class$gear_fishing)
SAU_mtc_gear_class$scientific_name <- as.factor(SAU_mtc_gear_class$scientific_name)
str(SAU_mtc_gear_class)

Gear_spp <- SAU_mtc_gear_class %>% 
  dplyr::group_by(gear_fishing) %>%
  distinct(scientific_name)

Species_num_gear <- Gear_spp %>%
  dplyr::group_by(scientific_name) %>%
  dplyr::summarise(n=n())

Gear_spp_temp <- inner_join(SAU_traits, Gear_spp, by="scientific_name")


Gear_mean_temp <- Gear_spp_temp %>%
  dplyr::group_by(gear_fishing) %>%
  dplyr::summarise(sd = sd(mean_temp),
                   mean_temp = mean(mean_temp),
                   n=n(),
                   se=sd/sqrt(n))


N_spp_by_gear <- SAU_mtc_gear_class %>% 
  dplyr::group_by(gear_fishing) %>%
  dplyr::summarise(spp = n_distinct(scientific_name))


# MTC by gear group
MTC_gear <- SAU_mtc_gear_class %>%
  dplyr::group_by(year, gear_fishing) %>%
  dplyr::summarise(MTC = sum(Annual_catch*mean_temp)/sum(Annual_catch))

MTC_gear$year <- as.numeric(as.character(MTC_gear$year))

ggplot(data = MTC_gear, aes(x=year, y=MTC)) + geom_path() + 
  facet_wrap(~gear_fishing) +
  stat_poly_line(col="black") +
  stat_poly_eq(use_label(c("eq", "R2"))) +
  scale_x_continuous(breaks=seq(1950,2020, by=10)) +
  theme_classic() +
  theme(axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=20),
        axis.title.x = element_text(size=30),
        axis.title.y = element_text(size=30))



# Bottom trawl
bottom.mtc <- lm(MTC ~ year, data=MTC_gear %>% filter(gear_fishing == "Bottom trawl"))

summary(bottom.mtc)

# get the fitted data
bottom.fitted.mtc <- fitted(bottom.mtc)
bottom.model.mtc <- data.frame(year = seq(1951,2019, 1), MTC = bottom.fitted.mtc)

bottom.seg.mtc <- segmented(bottom.mtc, 
                            seg.Z = ~ year,
                            psi = list(year = c(1987)))

summary.segmented(bottom.seg.mtc)

# get breakpoints
bottom.seg.mtc$psi # 1987
# get the slopes
slope(bottom.seg.mtc)

# P-score for significant non-similar slope
pscore.test(bottom.seg.mtc)

# get the fitted data
bottom.fitted.mtc <- fitted(bottom.seg.mtc)
bottom.model.mtc <- data.frame(year = seq(1951,2019, 1), MTC = bottom.fitted.mtc)

bottom <- ggplot(data = MTC_gear %>% filter(gear_fishing == "Bottom trawl"), aes(x=year, y=MTC)) +
  geom_point() +
  geom_line(cex=0.75) +
  geom_line(data=bottom.model.mtc, aes(x=year, y=MTC)) +
  #stat_poly_line(col="red") +
  #stat_poly_eq(use_label(c("eq", "R2")), size = 6, label.x = 0.5, label.y = 0.9) +
  scale_x_continuous(breaks=seq(1950,2020, by=10), limits = c(1950,2020)) +
  scale_y_continuous(breaks=seq(5,25, by=5), limits = c(5,25)) +
  labs(y = "MTC (°C)", x= "Year") +
  theme_classic() +
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=19),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.length=unit(.25, "cm"),
        plot.margin = margin(0, 0.5, 0, 0.5,"cm")) +
  annotate("text", x = 1975, y = 24, size = 8, label = "(a) Bottom trawl", col="black") +
  annotate("text", x = 2010, y = 22, size = 6, label = "n = 83", col="black") + 
  annotate("text", x = 1992, y = 19.75, size = 6, parse = T, label = "bar(x)", col="black") +
  annotate("text", x = 2008, y = 19.75, size = 6, label = "= 12.9 °C ± 5.5", col="black")


# Gillnet
gillnet.mtc <- lm(MTC ~ year, data=MTC_gear %>% filter(gear_fishing == "Gillnet"))

summary(gillnet.mtc)

# get the fitted data
gillnet.fitted.mtc <- fitted(gillnet.mtc)
gillnet.model.mtc <- data.frame(year = seq(1951,2019, 1), MTC = gillnet.fitted.mtc)

gillnet.seg.mtc <- segmented(gillnet.mtc, 
                             seg.Z = ~ year,
                             psi = list(year = c(1999)))

summary.segmented(gillnet.seg.mtc)

# get breakpoints
gillnet.seg.mtc$psi # 1999
# get the slopes
slope(gillnet.seg.mtc)

# P-score for significant non-similar slope
pscore.test(gillnet.seg.mtc)

# get the fitted data
gillnet.fitted.mtc <- fitted(gillnet.seg.mtc)
gillnet.model.mtc <- data.frame(year = seq(1951,2019, 1), MTC = gillnet.fitted.mtc)

gillnet <- ggplot(data = MTC_gear %>% filter(gear_fishing == "Gillnet"), aes(x=year, y=MTC)) +
  geom_point() +
  geom_line(cex=0.75) +
  geom_line(data=gillnet.model.mtc, aes(x=year, y=MTC)) +
  #stat_poly_line(col="red") +
  #stat_poly_eq(use_label(c("eq", "R2")), size = 6, label.x = 0.5, label.y = 0.9) +
  scale_x_continuous(breaks=seq(1950,2020, by=10), limits = c(1950,2020)) +
  scale_y_continuous(breaks=seq(5,25, by=5), limits = c(5,25)) +
  labs(y = "MTC (°C)", x= "Year") +
  theme_classic() +
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=19),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.length=unit(.25, "cm"),
        plot.margin = margin(0, 0.5, 0, 0.5,"cm")) +
  annotate("text", x = 1967, y = 24, size = 8, label = "(b) Gillnet", col="black") +
  annotate("text", x = 2010, y = 22, size = 6, label = "n = 21", col="black") + 
  annotate("text", x = 1992, y = 19.75, size = 6, parse = T, label = "bar(x)", col="black") +
  annotate("text", x = 2008, y = 19.75, size = 6, label = "= 17.3 °C ± 6.9", col="black")


# Hand lines
hand.mtc <- lm(MTC ~ year, data=MTC_gear %>% filter(gear_fishing == "Hand lines"))

summary(hand.mtc)

# get the fitted data
hand.fitted.mtc <- fitted(hand.mtc)
hand.model.mtc <- data.frame(year = seq(1951,2019, 1), MTC = hand.fitted.mtc)

hand <- ggplot(data = MTC_gear %>% filter(gear_fishing == "Hand lines"), aes(x=year, y=MTC)) +
  geom_point() +
  geom_line(cex=0.75) +
  geom_line(data=hand.model.mtc, aes(x=year, y=MTC)) +
  #stat_poly_line(col="red") +
  #stat_poly_eq(use_label(c("eq", "R2")), size = 6, label.x = 0.5, label.y = 0.9) +
  scale_x_continuous(breaks=seq(1950,2020, by=10), limits = c(1950,2020)) +
  scale_y_continuous(breaks=seq(5,25, by=5), limits = c(5,25)) +
  labs(y = "MTC (°C)", x= "Year") +
  theme_classic() +
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=19),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.length=unit(.25, "cm"),
        plot.margin = margin(0, 0.5, 0, 0.5,"cm")) +
  annotate("text", x = 1972, y = 24, size = 8, label = "(c) Hand lines", col="black") +
  annotate("text", x = 2010, y = 22, size = 6, label = "n = 9", col="black") + 
  annotate("text", x = 1992, y = 19.75, size = 6, parse = T, label = "bar(x)", col="black") +
  annotate("text", x = 2008, y = 19.75, size = 6, label = "= 20.7 °C ± 7.9", col="black")


# Longline
long.mtc <- lm(MTC ~ year, data=MTC_gear %>% filter(gear_fishing == "Longline"))

summary(long.mtc)

# get the fitted data
long.fitted.mtc <- fitted(long.mtc)
long.model.mtc <- data.frame(year = seq(1950,2019, 1), MTC = long.fitted.mtc)

long.seg.mtc <- segmented(long.mtc, 
                          seg.Z = ~ year,
                          psi = list(year = c(1969)))

summary.segmented(long.seg.mtc)

# get breakpoints
long.seg.mtc$psi # 1969
# get the slopes
slope(long.seg.mtc)

# P-score for significant non-similar slope
pscore.test(long.seg.mtc)

# get the fitted data
long.fitted.mtc <- fitted(long.seg.mtc)
long.model.mtc <- data.frame(year = seq(1950,2019, 1), MTC = long.fitted.mtc)

long <- ggplot(data = MTC_gear %>% filter(gear_fishing == "Longline"), aes(x=year, y=MTC)) +
  geom_point() +
  geom_line(cex=0.75) +
  geom_line(data=long.model.mtc, aes(x=year, y=MTC)) +
  #stat_poly_line(col="red") +
  #stat_poly_eq(use_label(c("eq", "R2")), size = 6, label.x = 0.5, label.y = 0.9) +
  scale_x_continuous(breaks=seq(1950,2020, by=10), limits = c(1950,2020)) +
  scale_y_continuous(breaks=seq(5,25, by=5), limits = c(5,25)) +
  labs(y = "MTC (°C)", x= "Year") +
  theme_classic() +
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=19),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.length=unit(.25, "cm"),
        plot.margin = margin(0, 0.5, 0, 0.5,"cm")) +
  annotate("text", x = 1972, y = 24, size = 8, label = "(d) Longline", col="black") +
  annotate("text", x = 2010, y = 22, size = 6, label = "n = 22", col="black") + 
  annotate("text", x = 1992, y = 19.75, size = 6, parse = T, label = "bar(x)", col="black") +
  annotate("text", x = 2008, y = 19.75, size = 6, label = "= 19.1 °C ± 7.5", col="black")


# Other
other.mtc <- lm(MTC ~ year, data=MTC_gear %>% filter(gear_fishing == "Other"))

summary(other.mtc)

# get the fitted data
other.fitted.mtc <- fitted(other.mtc)
other.model.mtc <- data.frame(year = seq(1955,2018, 1), MTC = other.fitted.mtc)

other.seg.mtc <- segmented(other.mtc, 
                           seg.Z = ~ year,
                           psi = list(year = c(1964)))

summary.segmented(other.seg.mtc)

# get breakpoints
other.seg.mtc$psi # 1964
# get the slopes
slope(other.seg.mtc)

# P-score for significant non-similar slope
pscore.test(other.seg.mtc)

# get the fitted data
other.fitted.mtc <- fitted(other.seg.mtc)
other.model.mtc <- data.frame(year = seq(1955,2018, 1), MTC = other.fitted.mtc)

other <- ggplot(data = MTC_gear %>% filter(gear_fishing == "Other"), aes(x=year, y=MTC)) +
  geom_point() +
  geom_line(cex=0.75) +
  geom_line(data=other.model.mtc, aes(x=year, y=MTC)) +
  #stat_poly_line(col="red") +
  #stat_poly_eq(use_label(c("eq", "R2")), size = 6, label.x = 0.5, label.y = 0.9) +
  scale_x_continuous(breaks=seq(1950,2020, by=10), limits = c(1950,2020)) +
  scale_y_continuous(breaks=seq(5,25, by=5), limits = c(5,25.97)) +
  labs(y = "MTC (°C)", x= "Year") +
  theme_classic() +
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=19),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.length=unit(.25, "cm"),
        plot.margin = margin(0, 0.5, 0, 0.5,"cm")) +
  annotate("text", x = 1966, y = 9, size = 8, label = "(e) Other", col="black") +
  annotate("text", x = 2010, y = 18, size = 6, label = "n = 8", col="black") + 
  annotate("text", x = 1992, y = 15.75, size = 6, parse = T, label = "bar(x)", col="black") +
  annotate("text", x = 2008, y = 15.75, size = 6, label = "= 24.1 °C ± 3.9", col="black")


# Purse seine
purse.mtc <- lm(MTC ~ year, data=MTC_gear %>% filter(gear_fishing == "Purse seine"))

summary(purse.mtc)

# get the fitted data
purse.fitted.mtc <- fitted(purse.mtc)
purse.model.mtc <- data.frame(year = seq(1951,2019, 1), MTC = purse.fitted.mtc)


purse <- ggplot(data = MTC_gear %>% filter(gear_fishing == "Purse seine"), aes(x=year, y=MTC)) +
  geom_point() +
  geom_line(cex=0.75) +
  geom_line(data=purse.model.mtc, aes(x=year, y=MTC)) +
  #stat_poly_line(col="red") +
  #stat_poly_eq(use_label(c("eq", "R2")), size = 6, label.x = 0.5, label.y = 0.9) +
  scale_x_continuous(breaks=seq(1950,2020, by=10), limits = c(1950,2020)) +
  scale_y_continuous(breaks=seq(5,25, by=5), limits = c(5,25.68)) +
  labs(y = "MTC (°C)", x= "Year") +
  theme_classic() +
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=19),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.length=unit(.25, "cm"),
        plot.margin = margin(0, 0.5, 0, 0.5,"cm")) +
  annotate("text", x = 1973, y = 9, size = 8, label = "(f) Purse seine", col="black") +
  annotate("text", x = 2010, y = 18, size = 6, label = "n = 15", col="black") + 
  annotate("text", x = 1992, y = 15.75, size = 6, parse = T, label = "bar(x)", col="black") +
  annotate("text", x = 2008, y = 15.75, size = 6, label = "= 21.5 °C ± 7.1", col="black")



# Small scale
small.mtc <- lm(MTC ~ year, data=MTC_gear %>% filter(gear_fishing == "Small scale"))

summary(small.mtc)

# get the fitted data
small.fitted.mtc <- fitted(small.mtc)
small.model.mtc <- data.frame(year = seq(1950,2019, 1), MTC = small.fitted.mtc)

small.seg.mtc <- segmented(small.mtc, 
                           seg.Z = ~ year,
                           psi = list(year = c(1989)))

summary.segmented(small.seg.mtc)

# get breakpoints
small.seg.mtc$psi # 1989
# get the slopes
slope(small.seg.mtc)

# P-score for significant non-similar slope
pscore.test(small.seg.mtc)

# get the fitted data
small.fitted.mtc <- fitted(small.seg.mtc)
small.model.mtc <- data.frame(year = seq(1950,2019, 1), MTC = small.fitted.mtc)

small <- ggplot(data = MTC_gear %>% filter(gear_fishing == "Small scale"), aes(x=year, y=MTC)) +
  geom_point() +
  geom_line(cex=0.75) +
  geom_line(data=small.model.mtc, aes(x=year, y=MTC)) +
  #stat_poly_line(col="red") +
  #stat_poly_eq(use_label(c("eq", "R2")), size = 6, label.x = 0.5, label.y = 0.9) +
  scale_x_continuous(breaks=seq(1950,2020, by=10), limits = c(1950,2020)) +
  scale_y_continuous(breaks=seq(5,25, by=5), limits = c(5,25)) +
  labs(y = "MTC (°C)", x= "Year") +
  theme_classic() +
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=19),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.length=unit(.25, "cm"),
        plot.margin = margin(0, 0.5, 0, 0.5,"cm")) +
  annotate("text", x = 1973, y = 24, size = 8, label = "(g) Small scale", col="black") +
  annotate("text", x = 2010, y = 22, size = 6, label = "n = 97", col="black") + 
  annotate("text", x = 1992, y = 19.75, size = 6, parse = T, label = "bar(x)", col="black") +
  annotate("text", x = 2008, y = 19.75, size = 6, label = "= 14.2 °C ± 6.3", col="black")


# Unknown
unknown.mtc <- lm(MTC ~ year, data=MTC_gear %>% filter(gear_fishing == "Unknown"))

summary(unknown.mtc)

# get the fitted data
unknown.fitted.mtc <- fitted(unknown.mtc)
unknown.model.mtc <- data.frame(year = seq(1951,2019, 1), MTC = unknown.fitted.mtc)

unknown.seg.mtc <- segmented(unknown.mtc, 
                             seg.Z = ~ year,
                             psi = list(year = c(1975)))

summary.segmented(unknown.seg.mtc)

# get breakpoints
unknown.seg.mtc$psi # 1975
# get the slopes
slope(unknown.seg.mtc)

# P-score for significant non-similar slope
pscore.test(unknown.seg.mtc)

# get the fitted data
unknown.fitted.mtc <- fitted(unknown.seg.mtc)
unknown.model.mtc <- data.frame(year = seq(1951,2019, 1), MTC = unknown.fitted.mtc)

unknown <- ggplot(data = MTC_gear %>% filter(gear_fishing == "Unknown"), aes(x=year, y=MTC)) +
  geom_point() +
  geom_line(cex=0.75) +
  geom_line(data=unknown.model.mtc, aes(x=year, y=MTC)) +
  #stat_poly_line(col="red") +
  #stat_poly_eq(use_label(c("eq", "R2")), size = 6, label.x = 0.5, label.y = 0.9) +
  scale_x_continuous(breaks=seq(1950,2020, by=10), limits = c(1950,2020)) +
  scale_y_continuous(breaks=seq(5,25, by=5), limits = c(5,25)) +
  labs(y = "MTC (°C)", x= "Year") +
  theme_classic() +
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=19),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.length=unit(.25, "cm"),
        plot.margin = margin(0, 0.5, 0, 0.5,"cm")) +
  annotate("text", x = 1971, y = 24, size = 8, label = "(h) Unknown", col="black") +
  annotate("text", x = 2010, y = 22, size = 6, label = "n = 85", col="black") + 
  annotate("text", x = 1992, y = 19.75, size = 6, parse = T, label = "bar(x)", col="black") +
  annotate("text", x = 2008, y = 19.75, size = 6, label = "= 12.7 °C ± 5.6", col="black")




Figure_5 <- ggarrange(bottom, gillnet, hand, long,
                      other, purse, small, unknown,
                      ncol = 2, nrow=4,
                      align = c("hv"))

Fig_5 = annotate_figure(Figure_5,
                        left = text_grob("MTC (°C)", size=26,  rot = 90),
                        bottom = text_grob("Year", size = 26)) +
  theme(plot.margin = margin(0, 0.5, 0, 0.5,"cm"))


#ggsave("Figure_5.tiff", Fig_5,
#       width = 28,
#       height = 36,
#       units = c("cm"),
#       device = c("tiff"))


# Figure 6

# Catch trends of Pelagic-oceanic
pelagic_oceanic <- SAU_milieu %>% filter(milieu == "pelagic-oceanic") %>%
  filter(scientific_name == "Katsuwonus pelamis" | scientific_name == "Thunnus maccoyii" |
           scientific_name == "Thunnus alalunga" |  scientific_name == "Cetorhinus maximus") %>%
  dplyr::group_by(scientific_name, Year) %>%
  dplyr::summarise(catch= sum(catch))

str(pelagic_oceanic)
pelagic_oceanic <- as.data.frame(pelagic_oceanic)

pelagic_oceanic$Year <- as.numeric(as.character(pelagic_oceanic$Year))
pelagic_oceanic$scientific_name <- as.factor( paste0( as.factor(pelagic_oceanic$scientific_name)))
str(pelagic_oceanic)

ggplot(data=pelagic_oceanic, aes(x=Year, y=catch)) +
  geom_path(aes(color=scientific_name)) +
  theme_classic()

names(pelagic_oceanic) <- c("Species", "Year", "Catch")

pelagic_order <- c("Katsuwonus pelamis", "Thunnus maccoyii", "Thunnus alalunga", "Cetorhinus maximus")


Fig_6c <- ggplot(data=pelagic_oceanic, aes(x=Year, y=Catch*0.001)) +
  geom_path(aes(color=factor(Species, levels = pelagic_order)), cex=1) +
  scale_color_manual(values = c("Katsuwonus pelamis" = "#F8766D",
                                "Thunnus maccoyii" = "green4",
                                "Thunnus alalunga" = "blue3",
                                "Cetorhinus maximus" = "deeppink2")) +
  theme_classic() +
  guides(color=guide_legend(title="Species")) +
  scale_x_continuous(breaks=seq(1960,2020, by=20), limits = c(1950,2020)) +
  scale_y_continuous(breaks=seq(0,20, 5), limits = c(0,20)) +
  theme(axis.ticks.length=unit(.25, "cm")) +
  theme(legend.title=element_text(size=12), 
        legend.text=element_text(size=11, face = "italic")) +
  theme(legend.position = c(0.23, 0.73)) +
  annotate("text", x = 1972, y = 20, size = 5,
           label = "(c) Pelagic-oceanic species", fontface = 1) +
  theme(axis.text.x = element_text(size=19),
        axis.text.y = element_text(size=19),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  theme(plot.margin = margin(0, 0.5, 0, 0,"cm")) +
  labs(y="Catch (tonnes × 1000)", x="Year")


# Purse seine
purse_seine <- SAU_mtc_gear_class %>% filter(gear_fishing == "Purse seine") %>%
  filter(scientific_name == "Katsuwonus pelamis" | scientific_name == "Scomber australasicus" |
           scientific_name == "Thunnus orientalis" |  scientific_name == "Coryphaena hippurus") %>%
  dplyr::group_by(scientific_name, year) %>%
  dplyr::summarise(catch= sum(Annual_catch))


str(purse_seine)
purse_seine <- as.data.frame(purse_seine)

purse_seine$year <- as.numeric(as.character(purse_seine$year))
purse_seine$scientific_name <- as.factor( paste0( as.factor(purse_seine$scientific_name)))
str(purse_seine)

ggplot(data=purse_seine, aes(x=year, y=catch)) +
  geom_path(aes(color=scientific_name)) +
  theme_classic()


names(purse_seine) <- c("Species", "Year", "Catch")

purse_order <- c("Katsuwonus pelamis", "Scomber australasicus", "Thunnus orientalis", "Coryphaena hippurus")


Fig_6d <- ggplot(data=purse_seine, aes(x=Year, y=Catch*0.001)) +
  geom_path(aes(colour=factor(Species, levels = purse_order)), cex=1) +
  scale_color_manual(values = c("Katsuwonus pelamis" = "#F8766D",
                                "Scomber australasicus" = "darkolivegreen3",
                                "Thunnus orientalis" = "darkorange1",
                                "Coryphaena hippurus" = "chocolate4")) +
  theme_classic() +
  guides(color=guide_legend(title="Species")) +
  scale_x_continuous(breaks=seq(1960,2020, by=20), limits = c(1950,2020)) +
  scale_y_continuous(breaks=seq(0,20, 5), limits = c(0,20)) +
  theme(axis.ticks.length=unit(.25, "cm")) +
  theme(legend.title=element_text(size=12), 
        legend.text=element_text(size=11, face = "italic")) +
  theme(legend.position = c(0.79, 0.77)) +
  annotate("text", x = 1970, y = 20, size = 5,
           label = "(d) Purse seine species", fontface = 1) +
  theme(axis.text.x = element_text(size=19),
        axis.text.y = element_text(size=19),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  theme(plot.margin = margin(0, 0, 0, 0,"cm")) +
  labs(y="Catch (tonnes × 1000)", x="Year")

# Benthopelagic
bentho_pelagic <- SAU_milieu %>% filter(milieu == "benthopelagic") %>%
  filter(scientific_name == "Macruronus novaezelandiae" | scientific_name == "Micromesistius australis" |
           scientific_name == "Trachurus declivis" |  scientific_name == "Thyrsites atun") %>%
  dplyr::group_by(scientific_name, Year) %>%
  dplyr::summarise(catch= sum(catch))

str(bentho_pelagic)
bentho_pelagic <- as.data.frame(bentho_pelagic)

bentho_pelagic$Year <- as.numeric(as.character(bentho_pelagic$Year))
bentho_pelagic$scientific_name <- as.factor( paste0( as.factor(bentho_pelagic$scientific_name)))
str(bentho_pelagic)

ggplot(data=bentho_pelagic, aes(x=Year, y=catch)) +
  geom_path(aes(color=scientific_name)) +
  theme_classic()

names(bentho_pelagic) <- c("Species", "Year", "Catch")

bentho_order <- c("Macruronus novaezelandiae", "Micromesistius australis", "Trachurus declivis", "Thyrsites atun")

# 
Fig_6a <- ggplot(data=bentho_pelagic, aes(x=Year, y=Catch*0.001)) +
  geom_path(aes(color=factor(Species, levels = bentho_order)), cex=1) +
  scale_color_manual(values = c("Macruronus novaezelandiae" = "slateblue1",
                                "Micromesistius australis" = "gold2",
                                "Trachurus declivis" = "green3",
                                "Thyrsites atun" = "sienna1")) +
  theme_classic() +
  guides(color=guide_legend(title="Species")) +
  scale_x_continuous(breaks=seq(1960,2020, by=20), limits = c(1950,2020)) +
  scale_y_continuous(breaks=seq(0,360, 60), limits = c(0,360), labels = scales::comma) +
  theme(axis.ticks.length=unit(.25, "cm")) +
  theme(legend.title=element_text(size=12), 
        legend.text=element_text(size=10, face = "italic")) +
  theme(legend.position = c(0.26, 0.74)) +
  annotate("text", x = 1973, y = 360, size = 5,
           label = "(a) Benthopelagic species", fontface = 1) +
  theme(axis.text.x = element_text(size=19),
        axis.text.y = element_text(size=19),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  theme(plot.margin = margin(0, 0.5, 0, 0,"cm")) +
  labs(y="Catch (tonnes × 1000)", x="Year")


# Bottom trawl
bottom_trawl <- SAU_mtc_gear_class %>% filter(gear_fishing == "Bottom trawl") %>%
  filter(scientific_name == "Macruronus novaezelandiae" | scientific_name == "Hoplostethus atlanticus" |
           scientific_name == "Thyrsites atun" |  scientific_name == "Pseudophycis bachus") %>%
  dplyr::group_by(scientific_name, year) %>%
  dplyr::summarise(catch= sum(Annual_catch))


str(bottom_trawl)
bottom_trawl <- as.data.frame(bottom_trawl)

bottom_trawl$year <- as.numeric(as.character(bottom_trawl$year))
bottom_trawl$scientific_name <- as.factor( paste0( as.factor(bottom_trawl$scientific_name)))
str(bottom_trawl)

ggplot(data=bottom_trawl, aes(x=year, y=catch)) +
  geom_path(aes(color=scientific_name)) +
  theme_classic()


names(bottom_trawl) <- c("Species", "Year", "Catch")

bottom_order <- c("Macruronus novaezelandiae", "Hoplostethus atlanticus", "Thyrsites atun", "Pseudophycis bachus")

Fig_6b <- ggplot(data=bottom_trawl, aes(x=Year, y=Catch*0.001)) +
  geom_path(aes(color=factor(Species, levels = bottom_order)), cex=1) +
  scale_color_manual(values = c("Macruronus novaezelandiae" = "slateblue1",
                                "Hoplostethus atlanticus" = "deeppink3",
                                "Thyrsites atun" = "sienna1",
                                "Pseudophycis bachus" = "aquamarine4")) +
  theme_classic() +
  guides(color=guide_legend(title="Species")) +
  scale_x_continuous(breaks=seq(1960,2020, by=20), limits = c(1950,2020)) +
  scale_y_continuous(breaks=seq(0,300, 60), limits = c(0,300), labels = scales::comma) +
  theme(axis.ticks.length=unit(.25, "cm")) +
  theme(legend.title=element_text(size=12), 
        legend.text=element_text(size=10, face = "italic")) +
  theme(legend.position = c(0.2615, 0.74)) +
  annotate("text", x = 1972, y = 300, size = 5,
           label = "(b) Bottom trawl species", fontface = 1) +
  theme(axis.text.x = element_text(size=19),
        axis.text.y = element_text(size=19),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  theme(plot.margin = margin(0, 0.5, 0, 0,"cm")) +
  labs(y="Catch (tonnes × 1000)", x="Year")


Fig_6 <- ggarrange(Fig_6a, Fig_6b,
                   Fig_6c, Fig_6d,
                   ncol=2, nrow=2,
                   align = c("hv"))

Fig_6 = annotate_figure(Fig_6,
                        left = text_grob("Catch"~(tons~x~10^3), size=26,  rot = 90),
                        bottom = text_grob("Year", size = 26)) +
  theme(plot.margin = margin(0, 0.5, 0, 0.5,"cm"))


#ggsave("Figure_6.tiff", Fig_6,
#       width = 28,
#       height = 24,
#       units = c("cm"),
#       device = c("tiff"))


# Calculating slope of each species' annual
# catch across the time series

str(SAU_catch)
SAU_catch$Year <- as.numeric(as.character(SAU_catch$Year))

dt <- data.table(SAU_catch %>% mutate(log_catch = log(catch)))
dt <- na.omit(dt)

slopes <- dt[,list(slope=lm(log_catch~Year)$coef),by=scientific_name]

dat1 <- SAU_catch %>% mutate(log_catch = log(catch))

get.coef <- function(dat1) lm(log_catch ~ Year, dat1)$coefficients

lm_results <- as.data.frame(plyr::ddply(dat1, .(scientific_name), get.coef))

lm_results <- lm_results[,c(1,3)]

SAU_traits_slope <- inner_join(SAU_traits, lm_results, by="scientific_name")

mean_catch <- SAU_catch %>%
  dplyr::group_by(scientific_name) %>%
  dplyr::summarise(sd = sd(catch),
                   mean_catch = mean(catch),
                   n=n(),
                   se=sd/sqrt(n))

SAU_traits_slope_mean <- inner_join(SAU_traits_slope, mean_catch, by="scientific_name")

