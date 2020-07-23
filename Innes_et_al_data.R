library(sp)
library(ncdf4)
library(tidyverse)
library(broom)
library(raster)
library(rgdal)
library(ggmap)
library(car)
library(emmeans)
library(multcomp)

select <- dplyr::select
#########################################################################
#                 extracting climatic data                  ############
#######################################################################

#figuring out days since Jan 1, 1900
monnb <- function(d) { lt <- as.POSIXlt(as.Date(d, origin="1900-01-01")); 
lt$year*12 + lt$mon }
mondf <- function(d1, d2) { monnb(d2) - monnb(d1) }
##https://stackoverflow.com/questions/1995933/number-of-months-between-two-dates

mondf("1900-1-1", "1943-1-1") # =516 days
mondf("1900-1-1", "2017-1-1") # =1404 days


#the climatic data

#datasets were downloaded from https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.02/cruts.1811131722.v4.02/

tmp <- "cru_ts4.02.1901.2017.tmp.dat.nc" #mean temperature

pre <- "cru_ts4.02.1901.2017.pre.dat.nc" #precipitation

pet <- "cru_ts4.02.1901.2017.pet.dat.nc" #potential evapotranspiration 


########################################################################
#dataset containing lat longs of all populations (historical and contemporary)
latlongs <- read_csv("pop_lat_longs.csv") %>%
  dplyr::select(ID, long, lat, Country, Continent, Hemisphere) 

#converting all population latlongs into SpatialPoints class object
points <- latlongs %>% 
  select(-ID, -Country, -Continent, -Hemisphere) %>% 
  SpatialPoints(.) 

## function to extract data 
get_raster_data <- function(netcdf, indices, variable){
  
  raster_data <- brick(netcdf, varname = variable)
  layer <- raster_data[[indices]]
  df <- raster::extract(layer, points, method = "simple", df = TRUE)
  return(df)
  }



###############################################################################
#taking desired time slices from climatic dateset (1943-2017)
df_tmp <- get_raster_data(tmp, 516:1404, variable = "tmp")

df_pre <- get_raster_data(pre, 516:1404, variable = "pre")

df_pet <- get_raster_data(pet, 516:1404, variable = "pet")

#########################################################################
#creating five year averages of mean winter temperature based on hemisphere and time of sampling

#northern hemisphere winter, historical
fiveyr_jan_hmean <- df_tmp %>% 
  select(., ID, matches('X195[0-4]\\.01\\.[0-9]{2}')) %>% 
  mutate(mean_temp = rowMeans(select(., -ID))) 

#southern hemisphere winter, historical
fiveyr_july_hmean <- df_tmp %>% 
  select(., ID, matches('X195[0-4]\\.07\\.[0-9]{2}')) %>% 
  mutate(mean_temp = rowMeans(select(., -ID))) 

#northern hemisphere winter, contemporary
fiveyr_jan_cmean <- df_tmp %>% 
  select(., ID, matches('X201[3-7]\\.01\\.[0-9]{2}')) %>% 
  mutate(mean_temp = rowMeans(select(., -ID))) 

#southern hemisphere winter, contemporary
fiveyr_july_cmean <- df_tmp %>% 
  select(., ID, matches('X201[3-7]\\.07\\.[0-9]{2}')) %>% 
  mutate(mean_temp = rowMeans(select(., -ID)))

################################################################################
#creating five year averages of PET

#historical
fiveyr_hpet <-  df_pet %>% 
  select(., ID, matches('X195[0-4]\\.[0-9]{2}\\.[0-9]{2}')) %>% 
  mutate(pet = rowMeans(select(., -ID)))

#contemporary
fiveyr_cpet <- df_pet %>% 
  select(., ID, matches('X201[3-7]\\.[0-9]{2}\\.[0-9]{2}')) %>% 
  mutate(pet = rowMeans(select(., -ID))) 


################################################################################
#creating five year averages of precipitation

mon_31 <- df_pre %>% 
  select(matches('X[0-9]{4}\\.01|03|05|07|08|10|12\\.[0-9]{2}')) %>%
  map_dfr(function(x) x / 31)

mon_30 <- df_pre %>% 
  select(matches('X[0-9]{4}\\.04|06|09|11\\.[0-9]{2}')) %>%
  map_dfr(function(x) x / 30)

h_mon_28 <- df_pre %>% 
  select(matches('X195[0-1|3-4]\\.02\\.[0-9]{2}')) %>% 
  map_dfr(function(x) x / 28)

c_mon_28 <- df_pre %>% 
  select(matches('X201[3-5|7]\\.02\\.[0-9]{2}')) %>% 
  map_dfr(function(x) x / 28)

mon_29 <- df_pre %>% 
  select(matches('X1952\\.02\\.|2016\\.02\\.')) %>% 
  map_dfr(function(x) x / 29)
  

all_pre <- bind_cols(mon_31, mon_30, h_mon_28, c_mon_28, mon_29) %>% 
  mutate(ID = 1:330)

#historical
fiveyr_hpre <-  all_pre %>% 
  select(., ID, matches('X195[0-4]\\.[0-9]{2}\\.[0-9]{2}')) %>% 
  mutate(precip = rowMeans(select(., -ID)))

#contemporary
fiveyr_cpre <- all_pre %>% 
  select(., ID, matches('X201[3-7]\\.[0-9]{2}\\.[0-9]{2}')) %>% 
  mutate(precip = rowMeans(select(., -ID))) 

################################################################################
#creating aridity indices

h_aridity <- fiveyr_hpre %>% 
  select(., ID, precip) %>% 
  left_join(., select(fiveyr_hpet, ID, pet), by = "ID") %>% 
  mutate(aridity = precip/pet) %>% 
  select(., ID, aridity)
  

c_aridity <- fiveyr_cpre %>% 
  select(., ID, precip) %>% 
  left_join(., select(fiveyr_cpet, ID, pet), by = "ID") %>% 
  mutate(aridity = precip/pet) %>% 
  select(., ID, aridity)


########################################################################
#sanity check to ensure extracted climatic data matches desired locations

#function to plot single layer and overlay points for whole globe
plot_raster_layer <- function(netcdf, index, variable){
  
  raster_data <- brick(netcdf, varname = variable)
  layer <- raster_data[[index]]
  plot(layer)
  plot(points, add = TRUE, pch = 1, cex = .2)
  return(df)
}
quick_plot <- plot_raster_layer(tmp, 1404, variable = "tmp")

#########################################################################
#                 calculating allele frequencies            ############
#######################################################################

########################################################################
#Importing contemporary cyanogenic data collected by SI and MJ

HCN_data <- readr::read_csv("cleanDaday2.0data.csv")


#function to calculate cyanotype frequencies and incorporate climatic data based on hemisphere
calculate_freq <- function(cytype, time, hemisphere){ 
  
  cytype <- enquo(cytype)
  
  HCN_data %>%  
    group_by(ID, city, lat, long, alt, Collector) %>% 
    summarise(total_plants = sum(!is.na(!! cytype)),
              total_pos = sum(!! cytype, na.rm = TRUE),
              freqHCN = total_pos / total_plants) %>% 
    left_join(., select(time, ID, mean_temp), by = "ID") %>% 
    left_join(., latlongs, by = c("ID", "lat", "long")) %>%
    filter(Hemisphere == hemisphere) %>% 
    na.omit() %>% 
    ungroup() -> df
  return(df)
  
  }

#HCN frequency, northern hemisphere
NH_HCNfreq <- calculate_freq(hcn, fiveyr_jan_cmean, "North")

#HCN frequency, southern hemisphere
SH_HCNfreq <- calculate_freq(hcn, fiveyr_july_cmean, "South")

#calculating allele frequencies#####################################################

calculate_allele_freq <- function(allele, time, hemisphere){ 
  
  allele <- enquo(allele)
  
  HCN_data %>%  
    select(ID, city, lat, long, alt, Collector, Ac, Li) %>% 
    group_by(ID, lat, long, alt) %>% 
    summarise(total_plants = sum(!is.na(!! allele)),
              total_pos = sum(!! allele, na.rm = TRUE),
              total_ress = total_plants - total_pos,
              freq = 1 - (total_ress / total_plants)^0.5) %>% 
    left_join(., select(time, ID, mean_temp), by = "ID") %>% 
    left_join(., latlongs, by = c("ID", "lat", "long")) %>%
    filter(Hemisphere == hemisphere) %>% 
    # na.omit() %>% 
    ungroup() -> df
  return(df)
  
}
#enquo() fix: https://dplyr.tidyverse.org/articles/programming.html?fbclid=IwAR0nSWgklzk_nNU6aCx6GVsWymOWLc-VYc56oTQupsVB_Lsq8Vev1EGfPv4

#Ac frequency, northern hemisphere
c_NH_Ac_freq <- calculate_allele_freq(Ac, fiveyr_jan_cmean, "North") %>% 
  rename(freqAc = freq) 

#Ac frequency, southern hemisphere
c_SH_Ac_freq <- calculate_allele_freq(Ac, fiveyr_july_cmean, "South") %>% 
  rename(freqAc = freq) 


#Li frequency, northern hemisphere
c_NH_Li_freq <-  calculate_allele_freq(Li, fiveyr_jan_cmean, "North") %>% 
  rename(freqLi = freq)

#Li frequency, southern hemisphere
c_SH_Li_freq <-  calculate_allele_freq(Li, fiveyr_july_cmean, "South") %>% 
  rename(freqLi = freq)


#putting together freqencies of HCN, Ac and Li ###############################
c_NH_freq_df <- NH_HCNfreq %>% 
  left_join(., select(c_NH_Ac_freq, ID, freqAc), by = 'ID' ) %>% 
  left_join(., select(c_NH_Li_freq, ID, freqLi), by = "ID")

c_SH_freq_df <- SH_HCNfreq %>% 
  left_join(., select(c_SH_Ac_freq, ID, freqAc), by = 'ID' ) %>%
  left_join(., select(c_SH_Li_freq, ID, freqLi), by = "ID")

c_freq_df <- bind_rows(c_SH_freq_df, c_NH_freq_df) %>% 
  mutate(., Locality = city) %>% 
  select(ID, Locality, lat, long, alt, freqHCN, freqAc, freqLi, Collector, Country, 
         Continent, Hemisphere, mean_temp) %>% 
  mutate("Data provenance" = "Present study")


## North American data collected by JS #######################################

NA_data <- read_csv("AllCities_AllPlants.csv")

#standardised distances from urban center used to only include rural populations
rudist <- read_csv("20170128 UAC_Distance.standardize.csv") 

#generating HCN frequencies
NA_HCN <- NA_data %>% 
  select(City, Population, Lat.pop, Long.pop, HCN_Result, Locus.Ac, Locus.Li, Distance) %>% 
  left_join(., select(rudist, City, Distance), by = "City") %>% 
  mutate(Rural = ifelse(Distance.x > Distance.y, 1, 0)) %>% 
  filter(Rural == 1) %>% 
  group_by(City) %>% 
  summarise(total_plants = sum(!is.na(HCN_Result)),
            total_hcn = sum(HCN_Result, na.rm = TRUE),
            freqHCN =  total_hcn/ total_plants) 

#allele frequencies for Ac and Li
each_locus <- function(allele) {
  
  allele <- enquo(allele)
  
  NA_data %>% 
    select(City, Population, Lat.pop, Long.pop, HCN_Result, Locus.Ac, Locus.Li, Distance) %>% 
    left_join(., select(rudist, City, Distance), by = "City") %>% 
    mutate(Rural = ifelse(Distance.x > Distance.y, 1, 0)) %>% 
    filter(Rural == 1) %>% 
    group_by(City) %>% 
    summarise(total_plants = sum(!is.na(!! allele)),
              total_pos = sum(!! allele, na.rm = TRUE),
              total_ress = total_plants - total_pos,
              freq = 1 - (total_ress / total_plants)^0.5) -> df
  return(df)
}

#for Ac
NA_Ac <- each_locus(Locus.Ac)

#for Li
NA_Li <- each_locus(Locus.Li)

#adding Ac and Li allele frequencies to HCN frequencies and incorporating climatic data
NA_freq <- NA_HCN %>% 
  left_join(., select(NA_Ac, freq, City), by = "City") %>% 
  rename(freqAc = freq) %>% 
  left_join(., select(NA_Li, freq, City), by = "City") %>% 
  rename(freqLi = freq) %>% 
  mutate(., ID = 166:181) %>% 
  left_join(., latlongs, by = "ID") %>%
  left_join(., select(fiveyr_jan_cmean, ID, mean_temp), by = "ID") %>% 
  mutate(., Locality = City)

#formatting column orders to match other datasets
NA_freq <- NA_freq %>% 
  mutate(., Collector = "innes") %>% 
  select(ID, Locality, lat, long, freqHCN, freqAc, freqLi, Collector, Country, Continent, 
         Hemisphere, mean_temp) %>% 
  mutate("Data provenance" = "Santangelo et al. 2020")


#North American and New Zealand data from NK and KO ####################

kooyers_NA_data <- read_csv("kooyers_NA_data.csv")

#taking the mean of populations to match "by city" sampling of SI and MJ and adding climatic data
kooyers_NA_freq <- kooyers_NA_data %>% 
  mutate(., Locality = Population) %>% 
  group_by(Locality) %>% 
  summarise(freqAc = mean(`Freq (Ac)`),
            freqLi = mean(`Freq (Li)`),
            freqHCN = mean(`Freq(AcLi)`),
            lat = mean(Latitude), 
            long = mean(Longitude)) %>%
  mutate(., ID = 311:320) %>%
  mutate(., Country = "USA") %>%
  mutate(., Continent = "North America") %>%
  mutate(., Collector = "innes") %>%
  mutate(., Hemisphere = "North") %>%
  left_join(., select(fiveyr_jan_cmean, ID, mean_temp), by = "ID") %>%
  select(ID, Locality, lat, long, freqHCN, freqAc, freqLi, Collector, Country, Continent, 
         Hemisphere, mean_temp) %>% 
  mutate("Data provenance" = "Kooyers et al. 2012")

########################################################################
kooyers_NZ_data <- read_csv("kooyers_NZ_data.csv")

kooyers_NZ_freq <- kooyers_NZ_data %>% 
  mutate(., Locality = Population) %>% 
  group_by(Locality) %>% 
  summarise(freqAc = mean(`Freq (Ac)`),
            freqLi = mean(`Freq (Li)`),
            freqHCN = mean(`Freq(AcLi)`),
            lat = mean(Latitude),
            long = mean(Longitude)) %>%
  mutate(., ID = 321:325) %>%
  mutate(., Country = "New Zealand") %>%
  mutate(., Continent = "Oceania") %>%
  mutate(., Collector = "innes") %>%
  mutate(., Hemisphere = "South") %>%
  left_join(., select(fiveyr_july_cmean, ID, mean_temp), by = "ID") %>%
  select(ID, Locality, lat, long, freqHCN, freqAc, freqLi, Collector, Country, Continent, 
         Hemisphere, mean_temp) %>% 
  mutate("Data provenance" = "Kooyers et al. 2013")

#######################################################################
kooyers_OKTN_data <- read_csv("OK_TN_data.csv")

kooyers_OKTN_freq <- kooyers_OKTN_data %>% 
  mutate(., Locality = Population) %>% 
  group_by(Locality) %>% 
  summarise(freqAc = mean(`Freq. (Ac)`),
            freqLi = mean(`Freq. (Li)`),
            freqHCN = mean(`Freq. (AcLi)`),
            lat = mean(Latitude), 
            long = mean(Longitude)) %>% 
  mutate(., ID = 326:329) %>%
  mutate(., Country = "USA") %>%
  mutate(., Continent = "North America") %>%
  mutate(., Collector = "innes") %>%
  mutate(., Hemisphere = "North") %>%
  left_join(., select(fiveyr_jan_cmean, ID, mean_temp), by = "ID") %>%
  select(ID, Locality, lat, long, freqHCN, freqAc, freqLi, Collector, Country, Continent, 
         Hemisphere, mean_temp) %>% 
  mutate("Data provenance" = "Kooyers et al. 2014")



########################################################################

#historical data extracted from Daday 1954, 1958
daday_data <- read_csv("Daday_data.csv")

#function adding mean winter temperatures to Daday data
add_daday_tmp <- function(hemisphere, time) {

  daday_data %>%
    select(., ID, Locality, lat, long, alt, freqHCN, freqAc, freqLi, Collector) %>%
    left_join(., latlongs, by = c("ID", "lat", "long")) %>%
    left_join(., select(time, ID, mean_temp), by = "ID") %>%
    filter(Hemisphere == hemisphere) -> df
  return(df)
  
}

#northern hemisphere
h_NH_freq_df <- add_daday_tmp('North', fiveyr_jan_hmean)

#southern hemisphere
h_SH_freq_df <- add_daday_tmp("South", fiveyr_july_hmean)

h_freq_df <- bind_rows(h_NH_freq_df, h_SH_freq_df) %>% 
  mutate("Data provenance" = if_else(Continent == "Europe", "Daday 1954", "Daday 1958"))

#Piecing together all data #########################################################

all_freq_df_sup <- bind_rows(c_freq_df, h_freq_df, NA_freq, kooyers_NA_freq, 
                             kooyers_NZ_freq, kooyers_OKTN_freq) %>% 
  filter(!(ID %in% c(2, 14))) %>% 
  filter(!(Country %in% c("India", "Iran", "USSR", "Australia"))) %>% 
  filter(Continent != "Africa")



all_freq_df <- bind_rows(c_freq_df, h_freq_df, NA_freq, kooyers_NA_freq, 
                         kooyers_NZ_freq, kooyers_OKTN_freq) %>% 
  filter(!(ID %in% c(2, 14))) %>% 
  filter(!(Country %in% c("India", "Iran", "USSR", "Australia"))) %>% 
  filter(Continent != "Africa") %>% 
  filter(long >= -98)

#filters are to remove strong outlier populations and historical populations from 
#outside contemporary sampling range

#adding aridity data
h_arid <- all_freq_df %>% 
  filter(Collector == "daday") %>% 
  left_join(., select(h_aridity, ID, aridity), by = "ID") %>% 
  select(ID, Collector, aridity)

c_arid <- all_freq_df %>% 
  filter(Collector == "innes") %>% 
  left_join(., select(c_aridity, ID, aridity), by = "ID") %>% 
  select(ID, Collector, aridity)

arid <- bind_rows(h_arid, c_arid)

all_freq_df <- all_freq_df %>% 
  left_join(., select(arid, ID, aridity), by = "ID")





#####################################################################################
#themes to pretty up figures to follow
ng1 <- theme(aspect.ratio=0.7,panel.background = element_blank(),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             panel.border=element_blank(),
             axis.line.x = element_line(color="black",size=1),
             axis.line.y = element_line(color="black",size=1),
             axis.ticks=element_line(color="black"),
             axis.text=element_text(color="black",size=15),
             axis.title=element_text(color="black",size=1),
             axis.title.y=element_text(vjust=2,size=17),
             axis.title.x=element_text(vjust=0.1,size=17),
             axis.text.x=element_text(size=15),
             axis.text.y=element_text(size=15),
             strip.text.x = element_text(size = 10, colour = "black",face = "bold"),
             strip.background = element_rect(colour="black"),
             legend.position = "top", legend.direction="vertical",
             legend.text=element_text(size=17), legend.key = element_rect(fill = "white"),
             legend.title = element_text(size=17),legend.key.size = unit(1.0, "cm"))

ng2 <- theme(aspect.ratio=0.7,panel.background = element_blank(),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             panel.border = element_blank(),
             axis.line.x = element_line(color="black",size=1),
             axis.line.y = element_line(color="black",size=1),
             axis.ticks=element_line(color="black"),
             axis.text=element_text(color="black",size=15),
             axis.title=element_text(color="black",size=1),
             axis.title.y=element_text(vjust=2,size=17),
             axis.title.x=element_text(vjust=0.1,size=17),
             axis.text.x=element_text(size=15),
             axis.text.y=element_text(size=15),
             strip.text.x = element_text(size = 10, colour = "black",face = "bold"),
             strip.background = element_rect(colour="black"),
             legend.position = "none")


#########################################################################
#                           analyses                        ############
#######################################################################

#manual contrasts for type III ANOVA

options(contrasts = c("contr.sum", "contr.poly"))

##################################################################################

#Assumptions of normality and homogeneity of variance 
par(mfrow = c(1,2))


plot(lm(freqHCN ~ lat + Collector + Continent + Collector:lat + 
          Continent:Collector + lat:Continent, data = all_freq_df))

plot(lm(freqAc ~ lat + Collector + Continent + Collector:lat + 
          Continent:Collector + lat:Continent, data = all_freq_df))

plot(lm(freqLi ~ lat + Collector + Continent + Collector:lat + 
          Continent:Collector + lat:Continent, data = all_freq_df))




plot(lm(freqHCN ~ mean_temp + Collector + Continent + Collector:mean_temp + 
          Continent:Collector + mean_temp:Continent, data = all_freq_df))

plot(lm(freqAc ~ mean_temp + Collector + Continent + Collector:mean_temp + 
          Continent:Collector + mean_temp:Continent, data = all_freq_df))

plot(lm(freqLi ~ mean_temp + Collector + Continent + Collector:mean_temp + 
          Continent:Collector + mean_temp:Continent, data = all_freq_df))

par(mfrow = c(1,1))



##################################################################################

#analysis by latitude of frequencies of HCN, Ac and Li through time

the_lanova <- function(var){
  
  var <- all_freq_df %>%  pull(var)
  lat <- all_freq_df %>% pull("lat")
  sample <- all_freq_df %>% pull("Collector")
  continent <- all_freq_df %>% pull("Continent")
  
  lat <- abs(lat)
  
  mod <- lm(var ~ lat + sample + continent + lat:sample + 
              lat:continent + sample:continent) 
  
  anova_mod <- Anova(mod, type = 3)
  
  print(AIC(mod)) 
  print(anova_mod)
  return(mod)
  
}

#ANOVA results and post hoc analyses by geographic range for HCN
the_lanova("freqHCN")
  
#ANOVA results and post hoc analyses by geographic range for Ac
the_lanova("freqAc")

#ANOVA results and post hoc analyses by geographic range for Li
the_lanova("freqLi")




##################################################################################

#analysis by mean winter temperature of frequencies of HCN, Ac and Li through time

the_tanova <- function(var){
  
  var <- all_freq_df %>% pull(var)
  mean_wt <- all_freq_df %>% pull("mean_temp")
  sample <- all_freq_df %>% pull("Collector")
  continent <- all_freq_df %>% pull("Continent")
  
  mod <- lm(var ~ mean_wt + sample + continent + 
              sample:mean_wt + sample:continent + mean_wt:continent)
  
  anova_mod <- Anova(mod, type = 3)
  
  print(AIC(mod))
  print(anova_mod)
  return(mod)
}

#ANOVA results for HCN
the_tanova("freqHCN")

#ANOVA results Ac
the_tanova("freqAc")

#ANOVA results Li
the_tanova("freqLi")





#########################################################################
#         plotting analyses by lat and temp through time    ############
#######################################################################
library(patchwork)
#plotting regression of cyanogenic frequencies by time of collection with respect to latitude

#for HCN
fig3a <- all_freq_df %>%
  group_by(Collector) %>%
  ggplot(data = ., aes(x = abs(lat), y = freqHCN, shape = Collector, color = Collector, linetype = Collector)) +
  geom_point(size = 2.5) +
  coord_cartesian(ylim = c(0, 1)) +
  geom_smooth(method = 'lm', se =FALSE, color = 'black') +
  scale_color_manual(values = c("darkorange2", "darkorchid4")) +
  labs(x = "", y = "HCN frequency", tag = expression(paste(~bold(A)))) +
  ng2 +
  annotate('text', x = 18.8, y = 0.47 , label = 'paste(beta[Cont] == -0.056)',
           parse = TRUE, size = 5) +
  annotate('text', x = 18.5, y = 0.37 , label = 'paste(beta[Hist] == -0.026)',
           parse = TRUE, size = 5) +
  annotate('text', x = 19.8, y = 0.57 , label = 'italic(P)[Lat~x~time] == 0.861',
           parse = TRUE, size = 5) +
  theme(axis.text.x = element_blank())

#for Ac
fig3b <- all_freq_df %>%
  group_by(Collector) %>%
  ggplot(data = ., aes(x = abs(lat), y = freqAc, shape = Collector, color = Collector, linetype = Collector)) +
  geom_point(size = 2.5) +
  coord_cartesian(ylim = c(0, 1))+
  geom_smooth(method = 'lm', se =FALSE, color = 'black') +
  scale_color_manual(values = c("darkorange2", "darkorchid4")) +
  labs(x = "", y = expression(paste(~italic(Ac)~ "frequency")), 
       tag = expression(paste(~bold(B)))) +
  ng2 + 
  annotate('text', x = 18.8, y = 0.47 , label = 'paste(beta[Cont] == -0.028)',
                 parse = TRUE, size = 5) +
  annotate('text', x = 17.9, y = 0.37 , label = 'paste(beta[Hist] == -0.030)',
           parse = TRUE, size = 5) +
  annotate('text', x = 19.8, y = 0.57 , label = 'italic(P)[Lat~x~time] == 0.278',
           parse = TRUE, size = 5) +
  theme(axis.text.x = element_blank())

#for Li
fig3c <- all_freq_df %>%
  group_by(Collector) %>%
  ggplot(data = ., aes(x = abs(lat), y = freqLi, shape = Collector, color = Collector, linetype = Collector)) +
  geom_point(size = 2.5) +
  coord_cartesian(ylim = c(0, 1))+
  geom_smooth(method = 'lm', se =FALSE, color = 'black') +
  scale_color_manual(values = c("darkorange2", "darkorchid4")) +
  ng2 +
  labs(x = "|Latitude|", y = expression(paste(~italic(Li)~ "frequency")), 
       tag = expression(paste(~bold(C)))) +
  annotate('text', x = 19, y = 0.4 , label = 'paste(beta[Cont] == -0.037)',
           parse = TRUE, size = 5) +
  annotate('text', x = 18.7, y = 0.3 , label = 'paste(beta[Hist] == -0.016)',
           parse = TRUE, size = 5) +
  annotate('text', x = 20, y = 0.5 , label = 'italic(P)[Lat~x~time] == 0.707',
           parse = TRUE, size = 5) 

#patching figures together
fig3 <- fig3a + fig3b + fig3c + plot_layout(nrow = 3, ncol = 1)
#aspect ratio = 600 x 1200


#plotting regression of cyanogenic frequencies by time of collection with respect to mean winter temperature

#native pops

#for HCN
fig4a <- all_freq_df %>%
  filter(Continent == "Europe") %>% 
  group_by(Collector) %>%
  ggplot(data = ., aes(x = mean_temp, y = freqHCN, shape = Collector, color = Collector, linetype = Collector)) +
  geom_point(size = 2.5) +
  coord_cartesian(ylim = c(0, 1)) +
  geom_smooth(method = 'lm', se =FALSE, color = 'black') +
  scale_color_manual(values = c("darkorange2", "darkorchid4")) +
  labs(title = ~bold("Native range"), x = "", y = "HCN frequency", 
       tag = expression(paste(~bold(A)))) +
  ng2 +
  annotate('text', x = -9.8, y = 0.57 , label = 'paste(beta[Cont] == 0.050)',
           parse = TRUE, size = 5) +
  annotate('text', x = -9.7, y = 0.47 , label = 'paste(beta[Hist] == 0.071)',
           parse = TRUE, size = 5) +
  annotate('text', x = -8.5, y = 0.67 , label = 'italic(P)[MWT~x~time] == 0.167',
           parse = TRUE, size = 5) +
  theme(plot.title = element_text(hjust = 0.5, size = 17), axis.text.x = element_blank())

#for Ac
fig4b <-all_freq_df %>%
  filter(Continent == "Europe") %>% 
  group_by(Collector) %>%
  ggplot(data = ., aes(x = mean_temp, y = freqAc, shape = Collector, color = Collector, linetype = Collector)) +
  geom_point(size = 2.5) +
  coord_cartesian(ylim = c(0, 1))+
  geom_smooth(method = 'lm', se =FALSE, color = 'black') +
  scale_color_manual(values = c("darkorange2", "darkorchid4")) +
  labs(x = "", y = expression(paste(~italic(Ac)~ "frequency")), tag = expression(paste(~bold(B)))) +
  ng2 +
  annotate('text', x = -9.8, y = 0.57 , label = 'paste(beta[Cont] == 0.053)',
           parse = TRUE, size = 5) +
  annotate('text', x = -9.9, y = 0.47 , label = 'paste(beta[Hist] == 0.077)',
           parse = TRUE, size = 5) +
  annotate('text', x = -8.85, y = 0.67 , label = 'italic(P)[MWT~x~time] == 0.072',
           parse = TRUE, size = 5) +
  theme(axis.text.x = element_blank())

#for Li
fig4c <- all_freq_df %>%
  filter(Continent == "Europe") %>% 
  group_by(Collector) %>%
  ggplot(data = ., aes(x = mean_temp, y = freqLi, shape = Collector, color = Collector, linetype = Collector)) +
  geom_point(size = 2.5) +
  coord_cartesian(ylim = c(0, 1))+
  geom_smooth(method = 'lm', se =FALSE, color = 'black') +
  scale_color_manual(values = c("darkorange2", "darkorchid4")) +
  labs(x = expression(paste("Mean winter temperature "(degree*C))),
       y = expression(paste(~italic(Li)~ "frequency")), 
       tag = expression(paste(~bold(C)))) +
  ng2 +
  annotate('text', x = -9.6, y = 0.57 , label = 'paste(beta[Cont] == 0.039)',
           parse = TRUE, size = 5) +
  annotate('text', x = -9.7, y = 0.47 , label = 'paste(beta[Hist] == 0.052)',
           parse = TRUE, size = 5) +
  annotate('text', x = -8.85, y = 0.67 , label = 'italic(P)[MWT~x~time] == 0.310',
           parse = TRUE, size = 5)


#introduced pops

#for HCN
fig4d <- all_freq_df %>%
  filter(Continent != "Europe") %>% 
  group_by(Collector) %>%
  ggplot(data = ., aes(x = mean_temp, y = freqHCN, shape = Collector, color = Collector, linetype = Collector)) +
  geom_point(size = 2.5) +
  coord_cartesian(ylim = c(0, 1)) +
  geom_smooth(method = 'lm', se =FALSE, color = 'black') +
  scale_color_manual(values = c("darkorange2", "darkorchid4")) +
  labs(title = ~bold("Introduced range"), x = "", y = "HCN frequency", 
       tag = expression(paste(~bold(D)))) +
  ng2 +
  annotate('text', x = -20.85, y = 0.66 , label = 'paste(beta[Cont] == 0.017)',
           parse = TRUE, size = 5) +
  annotate('text', x = -21, y = 0.56 , label = 'paste(beta[Hist] == 0.026)',
           parse = TRUE, size = 5) +
  annotate('text', x = -18.8, y = 0.76 , label = 'italic(P)[MWT~x~time] == 0.515',
           parse = TRUE, size = 5) +
  theme(plot.title = element_text(hjust = 0.5, size = 17), axis.text.x = element_blank())

#for Ac
fig4e <- all_freq_df %>%
  filter(Continent != "Europe") %>% 
  group_by(Collector) %>%
  ggplot(data = ., aes(x = mean_temp, y = freqAc, shape = Collector, color = Collector, linetype = Collector)) +
  geom_point(size = 2.5) +
  coord_cartesian(ylim = c(0, 1)) +
  geom_smooth(method = 'lm', se =FALSE, color = 'black') +
  scale_color_manual(values = c("darkorange2", "darkorchid4")) +
  labs(x = "", y = expression(paste(~italic(Ac)~ "frequency")), 
       tag = expression(paste(~bold(E)))) +
  ng2 +
  annotate('text', x = -21.3, y = 0.7 , label = 'paste(beta[Cont] == 0.010)',
           parse = TRUE, size = 5) +
  annotate('text', x = -21, y = 0.6 , label = 'paste(beta[Hist] == 0.041)',
           parse = TRUE, size = 5) +
  annotate('text', x = -19, y = 0.8 , label = 'italic(P)[MWT~x~time] == 0.954',
           parse = TRUE, size = 5) +
  theme(axis.text.x = element_blank())

#for Li
fig4f <- all_freq_df %>%
  filter(Continent != "Europe") %>% 
  group_by(Collector) %>%
  ggplot(data = ., aes(x = mean_temp, y = freqLi, shape = Collector, color = Collector, linetype = Collector)) +
  geom_point(size = 2.5) +
  coord_cartesian(ylim = c(0, 1)) +
  geom_smooth(method = 'lm', se =FALSE, color = 'black') +
  scale_color_manual(values = c("darkorange2", "darkorchid4")) +
  labs(x = expression(paste("Mean winter temperature "(degree*C))), 
       y = expression(paste(~italic(Li)~ "frequency")), 
       tag = expression(paste(~bold(F)))) +
  ng2 +
  annotate('text', x = -20.65, y = 0.63 , label = 'paste(beta[Cont] == 0.014)',
           parse = TRUE, size = 5) +
  annotate('text', x = -20.8, y = 0.53 , label = 'paste(beta[Hist] == 0.016)',
           parse = TRUE, size = 5) +
  annotate('text', x = -18.5, y = 0.73 , label = 'italic(P)[MWT~x~time] == 0.665',
           parse = TRUE, size = 5)


#patching figure together
fig4 <- (fig4a / fig4b / fig4c) | (fig4d / fig4e / fig4f)
#aspect ratio = 1400 x 1400

#################################################################################

#plotting changes in mean frequencies of cyanogenesis through time (from lm with latitude)

#for HCN
the_lanova("freqHCN") %>% 
  emmeans(., tukey ~ sample | continent) %>% 
  as.data.frame() %>% 
  ggplot(., aes(x = emmeans.continent, y = emmeans.emmean, shape = 
                  emmeans.sample, colour = emmeans.sample)) +
  geom_point(size = 2, position = position_dodge(width = 0.25)) +
  geom_errorbar(aes(ymin = emmeans.lower.CL, ymax = emmeans.upper.CL),
                position = position_dodge(width = 0.25), width = 0.15) +
  scale_color_manual(values = c("darkorange2", "darkorchid4")) +
  coord_cartesian(ylim = c(0, 1))+
  ylab("Mean HCN") + xlab("Continent") +
  ng1

#for Ac
the_lanova("freqAc") %>% 
  emmeans(., tukey ~ sample | continent) %>% 
  as.data.frame() %>% 
  ggplot(., aes(x = emmeans.continent, y = emmeans.emmean, shape = 
                  emmeans.sample, colour = emmeans.sample)) +
  geom_point(size = 2, position = position_dodge(width = 0.25)) +
  geom_errorbar(aes(ymin = emmeans.lower.CL, ymax = emmeans.upper.CL),
                position = position_dodge(width = 0.25), width = 0.15) +
  scale_color_manual(values = c("darkorange2", "darkorchid4")) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(y = expression(paste(~italic(Li) "Frequency")), 
       x = ("Continent")) +
  ng1

#for Li
the_lanova("freqLi") %>% 
  emmeans(., tukey ~ sample | continent) %>% 
  as.data.frame() %>% 
  ggplot(., aes(x = emmeans.continent, y = emmeans.emmean, shape = 
                  emmeans.sample, colour = emmeans.sample)) +
  geom_point(size = 2, position = position_dodge(width = 0.25)) +
  geom_errorbar(aes(ymin = emmeans.lower.CL, ymax = emmeans.upper.CL),
                position = position_dodge(width = 0.25), width = 0.15) +
  scale_color_manual(values = c("darkorange2", "darkorchid4")) +
  coord_cartesian(ylim = c(0, 1))+
  ylab("Mean") + xlab("Continent") +
  ng1

################################################################################
#plotting changes in mean frequencies of cyanogenesis through time (from lm with mean winter temperature)

#for HCN
the_tanova("freqHCN") %>% 
  emmeans(., tukey ~ sample | continent) %>% 
  as.data.frame() %>% 
  ggplot(., aes(x = emmeans.continent, y = emmeans.emmean, shape = 
                  emmeans.sample, colour = emmeans.sample)) +
  geom_point(size = 2, position = position_dodge(width = 0.25)) +
  geom_errorbar(aes(ymin = emmeans.lower.CL, ymax = emmeans.upper.CL),
                position = position_dodge(width = 0.25), width = 0.15) +
  scale_color_manual(values = c("darkorange2", "darkorchid4")) +
  coord_cartesian(ylim = c(0, 1))+
  ylab("Mean") + xlab("Continent") +
  ng1

#for Ac
the_tanova("freqAc") %>% 
  emmeans(., tukey ~ sample | continent) %>% 
  as.data.frame() %>% 
  ggplot(., aes(x = emmeans.continent, y = emmeans.emmean, shape = 
                  emmeans.sample, colour = emmeans.sample)) +
  geom_point(size = 2, position = position_dodge(width = 0.25)) +
  geom_errorbar(aes(ymin = emmeans.lower.CL, ymax = emmeans.upper.CL),
                position = position_dodge(width = 0.25), width = 0.15) +
  scale_color_manual(values = c("darkorange2", "darkorchid4")) +
  coord_cartesian(ylim = c(0, 1))+
  ylab("Mean") + xlab("Continent") +
  ng1

#for Li
the_tanova("freqLi") %>% 
  emmeans(., tukey ~ sample | continent) %>% 
  as.data.frame() %>% 
  ggplot(., aes(x = emmeans.continent, y = emmeans.emmean, shape = 
                  emmeans.sample, colour = emmeans.sample)) +
  geom_point(size = 2, position = position_dodge(width = 0.25)) +
  geom_errorbar(aes(ymin = emmeans.lower.CL, ymax = emmeans.upper.CL),
                position = position_dodge(width = 0.25), width = 0.15) +
  scale_color_manual(values = c("darkorange2", "darkorchid4")) +
  coord_cartesian(ylim = c(0, 1))+
  ylab("Mean") + xlab("Continent") +
  ng1

################################################################################

#########################################################################
#                       herbivory data                      ############
#######################################################################
#herbivory vs latitude

herb_dmg <- read_csv("herb_dmg.csv")

#calculation of mean herbivory by site and incorporation of HCN frequencies from sites for which we have them
euro_herbivory <- herb_dmg %>%  
  group_by(ID) %>% 
  mutate(mean_herb = log(mean(herb.dmg, na.rm = TRUE)+1)) %>%
  distinct(., ID, .keep_all = TRUE) %>% 
  select(., -herb.dmg, -pop.ind) %>% 
  left_join(., select(fiveyr_jan_cmean, mean_temp, ID), by = "ID") %>% 
  left_join(., select(all_freq_df, freqHCN, ID), by = "ID") %>% 
  na.omit() %>% 
  ungroup()

#ANOVA for the effects of latitude on herbivory with HCN frequency as covariate
herbivory <- function(var){
  
  df <- euro_herbivory
  
  latitude <- df %>% pull(var)
  herbivory <- df %>% pull("mean_herb")
  cyanogenesis <- df %>% pull("freqHCN")
  
  mod <- lm(herbivory ~ cyanogenesis + latitude + latitude:cyanogenesis) 
  
  aov_mod <- anova(mod)
  print(AIC(mod))
  print(mod)
  return(aov_mod)

}

#summary of ANOVA
herbivory("lat")

#creating residuals of mean herbivory to remove effects of HCN frequency on herbivory
euro_herbivory$residualHerb <- resid(lm(mean_herb ~ freqHCN, data = euro_herbivory))

#plotting residuals of mean herbivory against latitude
#fig 5 
euro_herbivory %>% 
  ggplot(data = ., aes(x = lat, y = residualHerb)) +
  geom_point(size = 2.5) +
  geom_smooth(method = 'lm', se =FALSE, color = 'black') +
  labs(x = "Latitude", y = "log(Herbivory) [residuals]") +
  ng1 +
  annotate('text', x = 40.75, y = 1.2 , label = 'paste(~italic(P)[Lat] == 0.005)',
           parse = TRUE, size = 5) +
  annotate('text', x = 42, y = 1 , label = 'paste(~italic(P)[Lat~x~HCN] == 0.006)',
           parse = TRUE, size = 5) +
  annotate('text', x = 40.8, y = 0.8 , label = 'paste(beta[Lat] == 0.014)',
           parse = TRUE, size = 5)



####################################################################################
#supplement
####################################################################################

#Table S1

Table_S1 <- all_freq_df %>% 
  mutate("Mean winter temp." = mean_temp) %>% 
  mutate("Time sampled" = if_else(Collector == "innes", "Contemporary", "Historical")) %>% 
  select(ID, Locality, lat, long, alt, freqHCN, freqAc, freqLi, "Mean winter temp.", 
         "Time sampled", Country, Continent, Hemisphere, "Data provenance") %>% 
  arrange(`Time sampled`, Hemisphere, Continent, Country)


#################################################################################

#Fitting binomial regression
the_bilanova <- function(var){
  
  var <- all_freq_df %>% pull(var)
  lat <- all_freq_df %>% pull("lat")
  sample <- all_freq_df %>% pull("Collector")
  continent <- all_freq_df %>% pull("Continent")
  
  lat <- abs(lat)
  
  mod <- glm(var ~ lat + sample + continent + lat:sample + 
               lat:continent + sample:continent, family = "quasibinomial") 
  
  anova_mod <- Anova(mod, type = 3)
  
  print(AIC(mod)) 
  return(anova_mod)
  
}

#Table S2
the_bilanova("freqHCN")

the_bilanova("freqAc")

the_bilanova("freqLi")


the_bitanova <- function(var){
  
  var <- all_freq_df %>% pull(var)
  mean_wt <- all_freq_df %>% pull("mean_temp")
  sample <- all_freq_df %>% pull("Collector")
  continent <- all_freq_df %>% pull("Continent")
  
  mod <- glm(var ~ mean_wt + sample + continent + mean_wt:sample + 
               mean_wt:continent + sample:continent, family = "quasibinomial") 
  
  anova_mod <- Anova(mod, type = 3)
  
  print(AIC(mod)) 
  return(anova_mod)
  
}

#Table S2
the_bitanova("freqHCN") 

the_bitanova("freqAc") 

the_bitanova("freqLi") 


#################################################################################

# models with three-way interactions

sup_lanova <- function(var){
  
  var <- all_freq_df %>%  pull(var)
  lat <- all_freq_df %>% pull("lat")
  sample <- all_freq_df %>% pull("Collector")
  continent <- all_freq_df %>% pull("Continent")
  
  lat <- abs(lat)
  
  mod <- lm(var ~ lat + sample + continent + lat:sample + 
              lat:continent + sample:continent + lat:sample:continent) 
  
  anova_mod <- Anova(mod, type = 3)
  
  print(AIC(mod)) 
  return(anova_mod)
  
}

#Table S3

#ANOVA results and post hoc analyses by geographic range for HCN
sup_lanova("freqHCN") 

#ANOVA results and post hoc analyses by geographic range for Ac
sup_lanova("freqAc")

#ANOVA results and post hoc analyses by geographic range for Li
sup_lanova("freqLi")


##################################################################################

sup_tanova <- function(var){
  
  var <- all_freq_df %>% pull(var)
  mean_wt <- all_freq_df %>% pull("mean_temp")
  sample <- all_freq_df %>% pull("Collector")
  continent <- all_freq_df %>% pull("Continent")
  
  mod <- lm(var ~ mean_wt + sample + continent + sample:mean_wt + 
              sample:continent + mean_wt:continent + mean_wt:sample:continent)
  
  anova_mod <- Anova(mod, type = 3)
  
  print(AIC(mod))
  return(anova_mod)
}

#Table S4

#ANOVA results for HCN
sup_tanova("freqHCN") 

#ANOVA results for Ac
sup_tanova("freqAc") 

#ANOVA results for Li
sup_tanova("freqLi")


#################################################################################

#Table S5

#post-hoc analyses for HCN
the_lanova("freqHCN") %>% emmeans(., tukey ~ sample|continent)

#post-hoc analyses for Ac
the_lanova("freqAc") %>% emmeans(., tukey ~ sample|continent)

#post-hoc analyses for Li
the_lanova("freqLi") %>% emmeans(., tukey ~ sample|continent)


#Table S6

#post-hoc analyses for HCN
the_tanova("freqHCN") %>% emmeans(., tukey ~ sample|continent)

#post-hoc analyses for Ac
the_tanova("freqAc") %>% emmeans(., tukey ~ sample|continent)

#post-hoc analyses for Li
the_tanova("freqLi") %>% emmeans(., tukey ~ sample|continent)


#################################################################################

#evaluating differences in climatic variables through time
pre_anova <- function(var){
  
  var <- all_freq_df %>% pull(var)
  lat <- all_freq_df %>% pull("lat")
  time <- all_freq_df %>% pull("Collector")
  range <- all_freq_df %>% pull("Continent")
  
  lat <- abs(lat)
  
  mod <- lm(var ~ lat + time + lat:time)
  anova_mod <- Anova(mod, type = 3)
  
  print(AIC(mod))
  print(anova_mod)
  return(mod)
}

#results of ANOVA + post hoc analysis by time for mean winter temperature 
pre_anova("mean_temp") %>% 
  emmeans(., tukey ~ time)

#visualisation of mean winter temp ~ latitude by time (Fig. S1)
all_freq_df %>% 
  group_by(Collector) %>%
  ggplot(data = ., aes(x = abs(lat), y = mean_temp, shape = Collector, color = Collector, linetype = Collector)) +
  geom_point(size = 2.5) +
  geom_smooth(method = 'lm', se =FALSE, color = 'black') +
  scale_color_manual(values = c("darkorange2", "darkorchid4")) +
  ng2 +
  labs(y = expression(paste("Mean winter temperature "(degree*C))),
       x = "|Latitude|")

#################################################################################

#the same test as above with inclusion of ommited populations
pre_supanova <- function(var){
  
  var <- all_freq_df_sup %>% pull(var)
  lat <- all_freq_df_sup %>% pull("lat")
  time <- all_freq_df_sup %>% pull("Collector")
  range <- all_freq_df_sup %>% pull("Continent")
  
  lat <- abs(lat)
  
  mod <- lm(var ~ lat + time + lat:time)
  anova_mod <- Anova(mod, type = 3)
  
  print(AIC(mod))
  print(anova_mod)
  return(mod)
}

#results of ANOVA + post hoc analysis by time for mean winter temperature 
pre_supanova("mean_temp") %>% 
  emmeans(., tukey ~ time)

#visualisation of mean winter temp ~ latitude by time
all_freq_df_sup %>% 
  group_by(Collector) %>%
  ggplot(data = ., aes(x = abs(lat), y = mean_temp, shape = Collector, color = Collector, linetype = Collector)) +
  geom_point(size = 2.5) +
  geom_smooth(method = 'lm', se =FALSE, color = 'black') +
  scale_color_manual(values = c("darkorange2", "darkorchid4")) +
  ng2 +
  labs(y = expression(paste("Mean winter temperature "(degree*C))),
       x = "|Latitude|")


###############################################################################

#Fig. S2

library(plotly)



p <- plot_ly(euro_herbivory, x = ~freqHCN, y = ~lat, z = ~residualHerb, 
             type = "scatter3d", mode = "markers", opacity = 1,
             #line = list(color = "black"),
             marker = list(size = 3.5, color = "black"))

p %>% layout(scene = list(xaxis = list(title = "HCN frequency"),
                          yaxis = list(title = "Latitude"),
                          zaxis = list(title = "log(Herbivory) [residuals]")))


#################################################################################

#results of ANOVA + post hoc analysis by time for drought (aridity)
pre_anova("aridity") %>% 
  emmeans(., tukey ~ time)

#visualisation of aridity ~ latitude split by time
all_freq_df %>% 
  group_by(Collector) %>%
  ggplot(data = ., aes(x = abs(lat), y = aridity, shape = Collector, color = Collector, linetype = Collector)) +
  geom_point(size = 2.5) +
  geom_smooth(method = 'lm', se =FALSE, color = 'black') +
  scale_color_manual(values = c("darkorange2", "darkorchid4")) +
  ng1








