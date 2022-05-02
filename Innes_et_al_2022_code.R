library(sp)
library(ncdf4)
library(tidyverse)
library(broom)
library(raster)
library(rgdal)
library(car)
library(emmeans)
library(multcomp)
library(visreg)
library(patchwork)
library(rsq)
library(interactions)

select <- dplyr::select
#########################################################################
#                 extracting climatic data                  ############
#######################################################################

#figuring out days since Jan 1, 1900 or 2011
monnb <- function(d) { lt <- as.POSIXlt(as.Date(d, origin="1900-01-01")); 
lt$year*12 + lt$mon }
mondf <- function(d1, d2) { monnb(d2) - monnb(d1) }
##https://stackoverflow.com/questions/1995933/number-of-months-between-two-dates

mondf("1900-1-1", "1943-1-1") # 516 months
mondf("1900-1-1", "2017-1-1") # 1404 months

difftime("2018-1-1", "2011-1-1", units = "days") # 2557 days
difftime("2018-9-1", "2011-1-1", units = "days") # 2800 days


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

#taking times slices from 1901-2017 for longer climatic dataset
df_tmp1 <- get_raster_data(tmp, 1:1404, variable = "tmp")

#



#########################################################################

#for temperature by latitude comparison through time
#historical temps, northern hemisphere
jan_hmean <- df_tmp1 %>% 
  select(., ID, matches('X19[0-4][0-9]\\.01\\.[0-9]{2}|X195[0-4]\\.01\\.[0-9]{2}')) %>% 
  mutate(long_temp = rowMeans(select(., -ID))) 

#contemporary temps, northern hemisphere
jan_cmean <- df_tmp1 %>% 
  select(., ID, matches('X195[5-9]\\.01\\.[0-9]{2}|X19[6-9][0-9]\\.01\\.[0-9]{2}|X20[0-9]{2}\\.01\\.[0-9]{2}')) %>% 
  mutate(long_temp = rowMeans(select(., -ID))) 

#historical temps, southern hemisphere
july_hmean <- df_tmp1 %>% 
  select(., ID, matches('X19[0-4][0-9]\\.07\\.[0-9]{2}|X195[0-4]\\.07\\.[0-9]{2}')) %>% 
  mutate(long_temp = rowMeans(select(., -ID))) 

#contemporary temps, southern hemisphere
july_cmean <- df_tmp1 %>% 
  select(., ID, matches('X195[5-9]\\.07\\.[0-9]{2}|X19[6-9][0-9]\\.07\\.[0-9]{2}|X20[0-9]{2}\\.07\\.[0-9]{2}')) %>% 
  mutate(long_temp = rowMeans(select(., -ID))) 

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
  mutate(ID = 1:344)

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
  mutate(aridity2 = precip-pet) %>% 
  select(., ID, aridity, aridity2)
  

c_aridity <- fiveyr_cpre %>% 
  select(., ID, precip) %>% 
  left_join(., select(fiveyr_cpet, ID, pet), by = "ID") %>% 
  mutate(aridity = precip/pet) %>% 
  mutate(aridity2 = precip-pet) %>% 
  select(., ID, aridity, aridity2)


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
    summarise(N = sum(!is.na(!! cytype)),
              total_pos = sum(!! cytype, na.rm = TRUE),
              freqHCN = total_pos / N) %>% 
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
  left_join(., select(c_NH_Li_freq, ID, freqLi), by = "ID") %>% 
  left_join(., select(jan_cmean, ID, long_temp), by = "ID")

c_SH_freq_df <- SH_HCNfreq %>% 
  left_join(., select(c_SH_Ac_freq, ID, freqAc), by = 'ID' ) %>%
  left_join(., select(c_SH_Li_freq, ID, freqLi), by = "ID") %>% 
  left_join(., select(july_cmean, ID, long_temp), by = "ID")

c_freq_df <- bind_rows(c_SH_freq_df, c_NH_freq_df) %>% 
  mutate(., Locality = city) %>% 
  select(ID, Locality, lat, long, alt, freqHCN, freqAc, freqLi, Collector, Country, 
         Continent, Hemisphere, mean_temp, long_temp, N) %>% 
  mutate("Data provenance" = "Present study") %>% 
  mutate("Year of collection" = "2018")


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
  mutate(., Locality = City) %>% 
  rename("N" = "total_plants")

#formatting column orders to match other datasets
NA_freq <- NA_freq %>% 
  mutate(., Collector = "innes") %>%
  left_join(., select(jan_cmean, ID, long_temp), by = "ID") %>% 
  select(ID, Locality, lat, long, freqHCN, freqAc, freqLi, Collector, Country, Continent, 
         Hemisphere, mean_temp, long_temp, N) %>% 
  mutate("Data provenance" = "Santangelo et al. 2020") %>% 
  mutate("Year of collection" = "2016") # except Boston, Montreal, Toronto and New York sampled in 2015


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
            long = mean(Longitude),
            N = sum(N)) %>%
  mutate(., ID = 311:320) %>%
  mutate(., Country = "USA") %>%
  mutate(., Continent = "North America") %>%
  mutate(., Collector = "innes") %>%
  mutate(., Hemisphere = "North") %>%
  left_join(., select(fiveyr_jan_cmean, ID, mean_temp), by = "ID") %>%
  left_join(., select(jan_cmean, ID, long_temp), by = "ID") %>% 
  select(ID, Locality, lat, long, freqHCN, freqAc, freqLi, Collector, Country, Continent, 
         Hemisphere, mean_temp, long_temp, N) %>% 
  mutate("Data provenance" = "Kooyers et al. 2012") %>% 
  mutate("Year of collection" = "2008–2009")

########################################################################
kooyers_NZ_data <- read_csv("kooyers_NZ_data.csv")

kooyers_NZ_freq <- kooyers_NZ_data %>% 
  mutate(., Locality = Population) %>% 
  group_by(Locality) %>% 
  summarise(freqAc = mean(`Freq (Ac)`),
            freqLi = mean(`Freq (Li)`),
            freqHCN = mean(`Freq(AcLi)`),
            lat = mean(Latitude),
            long = mean(Longitude),
            N = sum(N)) %>%
  mutate(., ID = 321:325) %>%
  mutate(., Country = "New Zealand") %>%
  mutate(., Continent = "Oceania") %>%
  mutate(., Collector = "innes") %>%
  mutate(., Hemisphere = "South") %>%
  left_join(., select(fiveyr_july_cmean, ID, mean_temp), by = "ID") %>% 
  left_join(., select(jan_cmean, ID, long_temp), by = "ID") %>% 
  select(ID, Locality, lat, long, freqHCN, freqAc, freqLi, Collector, Country, Continent, 
         Hemisphere, mean_temp, long_temp, N) %>% 
  mutate("Data provenance" = "Kooyers et al. 2013") %>% 
  mutate("Year of collection" = "2010")

#######################################################################
kooyers_OKTN_data <- read_csv("OK_TN_data.csv")

kooyers_OKTN_freq <- kooyers_OKTN_data %>% 
  mutate(., Locality = Population) %>% 
  group_by(Locality) %>% 
  summarise(freqAc = mean(`Freq. (Ac)`),
            freqLi = mean(`Freq. (Li)`),
            freqHCN = mean(`Freq. (AcLi)`),
            lat = mean(Latitude), 
            long = mean(Longitude),
            N = sum(`Total Plants`)) %>% 
  mutate(., ID = 326:329) %>%
  mutate(., Country = "USA") %>%
  mutate(., Continent = "North America") %>%
  mutate(., Collector = "innes") %>%
  mutate(., Hemisphere = "North") %>%
  left_join(., select(fiveyr_jan_cmean, ID, mean_temp), by = "ID") %>%
  left_join(., select(jan_cmean, ID, long_temp), by = "ID") %>% 
  select(ID, Locality, lat, long, freqHCN, freqAc, freqLi, Collector, Country, Continent, 
         Hemisphere, mean_temp, long_temp, N) %>% 
  mutate("Data provenance" = "Kooyers et al. 2014") %>% 
  mutate("Year of collection" = "2012")



########################################################################

#historical data extracted from Daday 1954, 1958
daday_data <- read_csv("Daday_data.csv")

#function adding mean winter temperatures to Daday data
add_daday_tmp <- function(hemisphere, time) {

  daday_data %>%
    select(., ID, Locality, lat, long, alt, freqHCN, freqAc, freqLi, Collector, no_plants) %>%
    rename("N" = "no_plants") %>% 
    left_join(., latlongs, by = c("ID", "lat", "long")) %>%
    left_join(., select(time, ID, mean_temp), by = "ID") %>%
    filter(Hemisphere == hemisphere) -> df
  return(df)
  
}

#northern hemisphere
h_NH_freq_df <- add_daday_tmp('North', fiveyr_jan_hmean) %>% 
  left_join(., select(jan_hmean, ID, long_temp), by = "ID")

#southern hemisphere
h_SH_freq_df <- add_daday_tmp("South", fiveyr_july_hmean) %>% 
  left_join(., select(july_hmean, ID, long_temp), by = "ID")

h_freq_df <- bind_rows(h_NH_freq_df, h_SH_freq_df) %>% 
  select(ID, Locality, lat, long, freqHCN, freqAc, freqLi, Collector, Country, 
         Continent, Hemisphere, mean_temp, long_temp, N) %>% 
  mutate("Data provenance" = if_else(Continent == "Europe", "Daday 1954", "Daday 1958")) %>% 
  mutate("Year of collection" = if_else(Continent == "Europe", "1954", "1958"))

#Piecing together all data #########################################################

all_freq_df_sup <- bind_rows(c_freq_df, h_freq_df, NA_freq, kooyers_NA_freq, 
                             kooyers_NZ_freq, kooyers_OKTN_freq) %>% 
  filter(!(ID %in% c(2, 14))) %>% 
  filter(!(Country %in% c("India", "Iran", "USSR", "Australia"))) %>% 
  filter(Continent != "Africa") %>% 
  mutate("Origin" = if_else(Continent == "Europe", "Native", "Introduced")) %>% 
  mutate(abslat = abs(lat))



all_freq_df <- bind_rows(c_freq_df, h_freq_df, NA_freq, kooyers_NA_freq, 
                         kooyers_NZ_freq, kooyers_OKTN_freq) %>% 
  filter(!(ID %in% c(2, 14))) %>% 
  filter(!(Country %in% c("India", "Iran", "USSR", "Australia", "Great Britain",
                          "Ireland"))) %>% 
 # filter(long >= -98) %>% 
  filter(Continent != "Africa") %>% 
  mutate("Origin" = if_else(Continent == "Europe", "Native", "Introduced")) %>% 
  mutate(abslat = abs(lat))

#filters are to remove strong outlier populations

#adding aridity data
h_arid <- all_freq_df %>% 
  filter(Collector == "daday") %>% 
  left_join(., select(h_aridity, ID, aridity, aridity2), by = "ID") %>% 
  select(ID, Collector, aridity, aridity2)

c_arid <- all_freq_df %>% 
  filter(Collector == "innes") %>% 
  left_join(., select(c_aridity, ID, aridity, aridity2), by = "ID") %>% 
  select(ID, Collector, aridity, aridity2)

arid <- bind_rows(h_arid, c_arid)

all_freq_df <- all_freq_df %>% 
  left_join(., select(arid, ID, aridity, aridity2), by = "ID")




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
#                       analyses & figures                  ############
#######################################################################

# manual contrasts for type III ANOVA

options(contrasts = c("contr.sum", "contr.poly"))

##################################################################################

delta_matrix <- read_csv("delta_matrix.csv")

delta_matrix_native <- delta_matrix %>% filter(Origin == "native")

# for native range
# for change in MWT versus change in HCN 
native_delta_lm <- lm(Delta_HCN ~ Delta_MWT, data = delta_matrix_native)

summary(native_delta_lm)

mean(delta_matrix_native$Delta_HCN)*100

# (Fig. 3A)
native_delta_lm %>%
  ggplot(data = ., aes(x = Delta_MWT, y = Delta_HCN)) +
  geom_point(size = 2.5, alpha = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.75) +
  #geom_smooth(method = 'lm', se =FALSE, colour = "black") +
  #scale_color_manual(values = c("darkgoldenrod2", "darkseagreen4")) +
  labs(x = expression(paste(Delta," Mean winter temperature "(degree*C))), 
       y = expression(paste(Delta," HCN frequency"))) +
  coord_cartesian(ylim = c(-0.6, 0.4), xlim = c(0.75, 3)) +
  ng2 

# for change in MWT versus change in Ac
native_delta_matrix2 <- delta_matrix_native %>% filter(Historical_Ac != 1)

native_delta_lm2 <- lm(Delta_Ac ~ Delta_MWT, data = native_delta_matrix2)

summary(native_delta_lm2) 

mean(native_delta_matrix2$Delta_Ac)*100

# (Fig. 3B)
native_delta_matrix2 %>% 
  ggplot(data = ., aes(x = Delta_MWT, y = Delta_Ac)) +
  geom_point(size = 2.5, alpha = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.75) +
  #geom_smooth(method = 'lm', se =FALSE, colour = "black") +
  #scale_color_manual(values = c("darkgoldenrod2", "darkseagreen4")) +
  labs(x = expression(paste(Delta," Mean winter temperature "(degree*C))), 
       y = expression(paste(Delta, ~italic(Ac)~"frequency"))) +
  coord_cartesian(ylim = c(-0.6, 0.4), xlim = c(0.75, 3)) +
  ng2 

# for change in MWT versus change in Li
native_delta_matrix3 <- delta_matrix_native %>% filter(Historical_Li != 1)

native_delta_lm3 <- lm(Delta_Li ~ Delta_MWT, data = native_delta_matrix3)

summary(native_delta_lm3)

mean(native_delta_matrix3$Delta_Li)*100

# (Fig. 3C)
native_delta_matrix3 %>% 
  ggplot(data = ., aes(x = Delta_MWT, y = Delta_Li)) +
  geom_point(size = 2.5, alpha = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.75) +
  #geom_smooth(method = 'lm', se =FALSE, colour = "black") +
  #scale_color_manual(values = c("darkgoldenrod2", "darkseagreen4")) +
  labs(x = expression(paste(Delta," Mean winter temperature "(degree*C))), 
       y = expression(paste(Delta, ~italic(Li)~"frequency"))) +
  coord_cartesian(ylim = c(-0.6, 0.4), xlim = c(0.75, 3)) +
  ng2 


# for native + introduced
# for change in MWT versus change in HCN
delta_lm <- lm(Delta_HCN ~ Delta_MWT, data = delta_matrix)

summary(delta_lm)

# (Fig. S1A)
delta_matrix %>%
  ggplot(data = ., aes(x = Delta_MWT, y = Delta_HCN, colour = Origin)) +
  geom_point(size = 2.5, alpha = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.75) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.75) +
  #geom_smooth(method = 'lm', se =FALSE, colour = "black") +
  scale_color_manual(values = c("darkgoldenrod2", "darkseagreen4")) +
  labs(x = expression(paste(Delta," Mean winter temperature "(degree*C))), 
       y = expression(paste(Delta," HCN frequency"))) +
  coord_cartesian(xlim = c(-3.7,3), ylim = c(-0.85, 0.55)) +
  ng2 

# for change in MWT versus change in Ac
delta_matrix2 <- delta_matrix %>% filter(Historical_Ac != 1)

delta_lm2 <- lm(Delta_Ac ~ Delta_MWT, data = delta_matrix2)

summary(delta_lm2) 

# (Fig. S1B)
delta_matrix2 %>% 
  ggplot(data = ., aes(x = Delta_MWT, y = Delta_Ac, colour = Origin)) +
  geom_point(size = 2.5, alpha = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.75) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.75) +
  #geom_smooth(method = 'lm', se =FALSE, colour = "black") +
  scale_color_manual(values = c("darkgoldenrod2", "darkseagreen4")) +
  labs(x = expression(paste(Delta," Mean winter temperature "(degree*C))), 
       y = expression(paste(Delta, ~italic(Ac)~"frequency"))) +  
  coord_cartesian(xlim = c(-3.7,3), ylim = c(-0.85, 0.55)) +
  ng2

# for change in MWT versus change in Li
delta_matrix3 <- delta_matrix %>% filter(Historical_Li != 1)

delta_lm3 <- lm(Delta_Li ~ Delta_MWT, data = delta_matrix3)

summary(delta_lm3)

# (Fig. S1C)
delta_matrix3 %>% 
  ggplot(data = ., aes(x = Delta_MWT, y = Delta_Li, colour = Origin)) +
  geom_point(size = 2.5, alpha = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.75) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.75) +
  #geom_smooth(method = 'lm', se =FALSE, colour = "black") +
  scale_color_manual(values = c("darkgoldenrod2", "darkseagreen4")) +
  labs(x = expression(paste(Delta," Mean winter temperature "(degree*C))),
       y = expression(paste(Delta, ~italic(Li)~"frequency"))) +
  coord_cartesian(xlim = c(-3.7,3), ylim = c(-0.85, 0.55)) +
  ng2


################################################################################




#########################################################################
#                       herbivory data                      ############
#######################################################################
#herbivory vs latitude (cyanogenesis) (fig. 4)
##################################################################################
herb_dmg <- read_csv("herb_dmg.csv")



#calculation of mean herbivory by site and incorporation of HCN frequencies from sites for which we have them
euro_herbivory <- herb_dmg %>%  
  group_by(ID) %>% 
  mutate(mean_herb = mean(herb.dmg, na.rm = TRUE)) %>% 
  distinct(., ID, .keep_all = TRUE) %>% 
  select(., -herb.dmg, -pop.ind) %>% 
  left_join(., select(fiveyr_jan_cmean, mean_temp, ID), by = "ID") %>% 
  left_join(., select(all_freq_df, freqHCN, ID), by = "ID") %>% 
  na.omit() %>% 
  ungroup()



#################################################################################

# daily mean temp data for EU from https://cds.climate.copernicus.eu/cdsapp#!/dataset/insitu-gridded-observations-europe?tab=overview

# daily min temp for EU

dtmp_min <- "tn_ens_mean_0.1deg_reg_2011-2020_v23.1e.nc"

df_dmin <- get_raster_data(dtmp_min, 2558:2800, variable = "tn")

#daily max temp for EU

dtmp_max <- "tx_ens_mean_0.1deg_reg_2011-2020_v23.1e.nc"

df_dmax <- get_raster_data(dtmp_max, 2558:2800, variable = "tx")


# average daily temps (min + max)/2

avg_dtmp <- (df_dmax + df_dmin)/2



# setting daily temps to degrees above baseiine of 5C and setting upper threshold of 30C
set_gdd <- function(x){
  y <- case_when(x > 5 & x <= 30 ~ x -5, # Tbase is subtracted from Tavg
                 x > 30 ~ 25, # Tavg above upper threshold is set to threshold and subtracted by Tbase
                 is.na(x) ~ 0, # NAs are converted to 0
                 TRUE ~ 0) # Tavg below 5 do not contribute GDDs (i.e., set to Tbase and then subtracted by Tbase)
}


# Custom function that returns sum of specific range of columns based on column "date"
spec_gdd <- function(x){
  end_date <- as.character(x %>% pull(date_collected))
  sum(select(x, X2018.01.01:!!sym(end_date)))
}

# add in collection dates for each site
collection_dates <- read_csv("HCN_collection_date.csv")


# Create modified dataset
euro_herbivory_gdd <- euro_herbivory %>% 
  left_join(., select(collection_dates, ID, date_collected), by = "ID") %>% 
  left_join(., avg_dtmp, by = "ID") %>% 
  # as.data.frame() %>%   # Need to convert to dataframe first
  mutate(across(X2018.01.01:X2018.08.31, ~as.numeric(.))) %>%  # Ensure columns are numeric
  mutate(across(X2018.01.01:X2018.08.31, ~set_gdd(.))) %>%   # Apply function to all columns from 'X.a' to 'X.c'. Modified columns in place
  # mutate(across(X.a:X.c, ~test_fun(.), .names="{.col}_trans")) %>%   # Apply function to all columns from 'X.a' to 'X.c'. Use this version to add new columns named "<col>_trans"
  mutate(total_gdd = rowSums(across(X2018.01.01:X2018.08.31))) %>%    # Apply function to all columns from 'X.a' to 'X.c'. Modified columns in place
  rowwise() %>%  # Create Row-wise dataframe
  do( X = as_data_frame(.) ) %>%  # Nest entire dataframe as cell value for each row. More details https://dplyr.tidyverse.org/articles/rowwise.html
  ungroup() %>%
  mutate(cum_gdd = map(X, spec_gdd)) %>%  # Map custom funtion that take a row as input and outputs the sum of specific range of columns based on "date"
  unnest(c(X, cum_gdd)) %>% # Unnest nested dataframe column to return to original dataframe format
  select(., -matches("X20[0-9]{2}")) 

# adding closest GDD to missing populations
euro_herbivory_gdd[5, 12] = euro_herbivory_gdd[4, 12] # San Sebastian
euro_herbivory_gdd[47, 12] = euro_herbivory_gdd[46, 12] # Umea
  



write_csv(euro_herbivory_gdd, path = "euro_herbivory.csv")




################################################################################
# herbivory interaction plot (fig. 4)
################################################################################

herb_lat <- lm(log1p(mean_herb) ~ freqHCN + lat + lat:freqHCN + cum_gdd, 
               data = euro_herbivory_gdd)

Anova(herb_lat, type = 3) 

#####

interact_plot(herb_lat, lat, freqHCN, plot.points = TRUE, point.size = 2.5,
              point.alpha = 0.9, data = euro_herbivory_gdd) +
  labs(x = "Latitude", y = "log(Herbivory)") +
  ng2 #+
  # theme(axis.text.x = element_text(size=23),
  #       axis.title.x = element_text(vjust=2,size=25),
  #       axis.title.y = element_text(vjust=2,size=25),
  #       axis.text.y = element_text(size=23))

sim_slopes(herb_lat, lat, freqHCN)


#####
# for finding averages of herbivory across latitudinal gradient

# max lat - mins lat to ind high lat herbivory
max(herb_dmg$lat)-min(herb_dmg$lat) 
# 28.44989

# above distance multiplied by slope plus intercept (y=mx+b)
28.44989*0.05+1.167402
# 2.589897

# transformation to remove log-scaling
exp(2.589897)-1
# 12.3284

# using intercept for low lat herbivory
exp(1.167402)-1
# 2.213633


################################################################################



#################################################################################
#analyses for test of conituned adaptation in the introduced range
#for HCN
#native vs introduced (fig. 5a & 5d)
##################################################################################
nat_intro_HCN <- glm(freqHCN ~ mean_temp + Collector + Origin + mean_temp:Collector + 
                     mean_temp:Origin + Origin:Collector + 
                     Origin:Collector:mean_temp, data = all_freq_df,
                     weights = N, family = "quasibinomial")

Anova(nat_intro_HCN, type = 3, test.statistic = "F")

#qualitatively same between sup and maintext models


sim_slopes(nat_intro_HCN, mean_temp, Collector, Origin, mod2.values = "Native")

sim_slopes(nat_intro_HCN, mean_temp, Collector, Origin, mod2.values = "Introduced")


#w/out sup data: daday: slope = 0.04  S.E = 0.02  P = 0.11
#                innes: slope = 0.11  S.E = 0.02  P = 0.00
#w/ sup data:    daday: slope = 0.04  S.E = 0.01  P = 0.00
#                innes: slope = 0.11  S.E = 0.02  P = 0.00

fig5_HCN <- interact_plot(nat_intro_HCN, mean_temp, Collector, mod2 = Origin, 
                          mod2.values = c("Native", "Introduced"), 
                          plot.points = TRUE, point.size = 2.5, point.shape = TRUE) +
  labs(x = "", y = expression(paste("HCN frequency"))) +
  geom_line(aes(col = Collector)) +
  scale_color_manual(values = c("darkorange2", "darkorchid4")) +
  coord_cartesian(ylim = c(0, 1)) +
  ng2 

# aspect ratio = 12.5" x 5.57" for PDF portait

#################################################################################
#for Ac
#native vs introduced (fig. 5b & 5e)
##################################################################################
nat_intro_Ac <- glm(freqAc ~ mean_temp + Collector + Origin + mean_temp:Collector + 
                    mean_temp:Origin + Origin:Collector + 
                    Origin:Collector:mean_temp, data = all_freq_df, 
                    weights = N, family = "quasibinomial")

Anova(nat_intro_Ac, type = 3, test.statistic = "F") 

#qualitatively same between sup and maintext models


sim_slopes(nat_intro_Ac, mean_temp, Collector, Origin, mod2.values = "Native")

sim_slopes(nat_intro_Ac, mean_temp, Collector, Origin, mod2.values = "Introduced")


#w/out sup data: daday: slope = 0.04  S.E = 0.02  P = 0.04
#                innes: slope = 0.07  S.E = 0.02  P = 0.00
#w/ sup data:    daday: slope = 0.04  S.E = 0.01  P = 0.00
#                innes: slope = 0.06  S.E = 0.02  P = 0.00

fig5_Ac <- interact_plot(nat_intro_Ac, mean_temp, Collector, mod2 = Origin,
                         mod2.values = c("Native", "Introduced"),
                         plot.points = TRUE, point.size = 2.5, point.shape = TRUE) +
  labs(x = "", y = expression(paste(~italic(Ac)~ "frequency"))) +
  geom_line(aes(col = Collector)) +
  scale_color_manual(values = c("darkorange2", "darkorchid4")) +
  coord_cartesian(ylim = c(0, 1)) +
  ng2

# aspect ratio = 12.5" x 5.57" for PDF portait

#################################################################################
#for Li
#native vs introduced (fig. 5c & 5f)
##################################################################################
nat_intro_Li <- glm(freqLi ~ mean_temp + Collector + Origin + mean_temp:Collector + 
                    mean_temp:Origin + Origin:Collector + 
                    Origin:Collector:mean_temp, data = all_freq_df, 
                    weights = N, family = "quasibinomial")

Anova(nat_intro_Li, type = 3, test.statistic = "F") 
#qualitatively same between sup and maintext models


sim_slopes(nat_intro_Li, mean_temp, Collector, Origin, mod2.values = "Native")

sim_slopes(nat_intro_Li, mean_temp, Collector, Origin, mod2.values = "Introduced")

#w/out sup data: daday: slope = 0.03  S.E = 0.02  P = 0.09
#                innes: slope = 0.09  S.E = 0.02  P = 0.00
#w/ sup data:    daday: slope = 0.02  S.E = 0.01  P = 0.06
#                innes: slope = 0.09  S.E = 0.02  P = 0.00



fig5_Li <- interact_plot(nat_intro_Li, mean_temp, Collector, mod2 = Origin,
                         mod2.values = c("Native", "Introduced"),
                         plot.points = TRUE, point.size = 2.5, point.shape = TRUE) +
  labs(x = expression(paste("Mean winter temperature "(degree*C))), 
       y = expression(paste(~italic(Li)~ "frequency"))) +
  geom_line(aes(col = Collector)) +
  scale_color_manual(values = c("darkorange2", "darkorchid4")) +
  coord_cartesian(ylim = c(0, 1)) +
  ng2 

# aspect ratio = 12.5" x 5.57" for PDF portrait

##################################################################################
# patching panels for fig. 5
##################################################################################
fig5_new <- (fig5_HCN / fig5_Ac / fig5_Li)

#aspect ratio = 1400 x 1400

##################################################################################

####################################################################################
#         comparison of just NA and EU           ##################################
##################################################################################
EU_NA_df <- all_freq_df %>% 
  filter(!(Continent %in% c("Asia", "Oceania", "South America")))


EU_NA_df_sup <- all_freq_df_sup %>% 
  filter(!(Continent %in% c("Asia", "Oceania", "South America")))

##################################################################################
#by temp

#for HCN
##################################################################################
EU_vs_NA_HCNtemp <- glm(freqHCN~mean_temp + Collector + Origin + mean_temp:Collector + 
                        mean_temp:Origin + Collector:Origin + Collector:Origin:mean_temp, 
                        weights = N, data = EU_NA_df, family = "quasibinomial")

Anova(EU_vs_NA_HCNtemp, type = 2, test.statistic = "F")

#w/ out mean_temp:Collector:Origin    2.123  1  0.1450946    
#w/ sup mean_temp:Collector:Origin    3.217  1  0.0728821


sim_slopes(EU_vs_NA_HCNtemp, mean_temp, Collector, Origin, mod2.values = "Introduced")

sim_slopes(EU_vs_NA_HCNtemp, mean_temp, Collector, Origin, mod2.values = "Native")



#w/out sup data: daday: slope = 0.07  S.E = 0.04  P = 0.06
#                innes: slope = 0.08  S.E = 0.03  P = 0.00
#w/ sup data:    daday: slope = 0.05  S.E = 0.01  P = 0.00
#                innes: slope = 0.08  S.E = 0.03  P = 0.00

#
interact_plot(EU_vs_NA_HCNtemp, mean_temp, Collector, mod2 = Origin,
                         mod2.values = c("Native", "Introduced"),
                         plot.points = TRUE, point.size = 2.5, point.shape = TRUE) +
  labs(x = expression(paste("Mean winter temperature "(degree*C))), 
       y = expression(paste("HCN frequency"))) +
  geom_line(aes(col = Collector)) +
  scale_color_manual(values = c("darkorange2", "darkorchid4")) +
  coord_cartesian(ylim = c(0, 1)) +
  ng2 




##################################################################################

#for Ac
##################################################################################
EU_vs_NA_Actemp <- glm(freqAc~mean_temp + Collector + Origin + mean_temp:Collector + 
                       mean_temp:Origin + Collector:Origin + Collector:Origin:mean_temp, 
                       weights = N, data = EU_NA_df, family = "quasibinomial")

Anova(EU_vs_NA_Actemp, type = 2, test.statistic = "F")

#qualitatively same between sup and maintext models


sim_slopes(EU_vs_NA_Actemp, mean_temp, Collector, Origin, mod2.values = "Introduced")

sim_slopes(EU_vs_NA_Actemp, mean_temp, Collector, Origin, mod2.values = "Native")


#w/out sup data: daday: slope = 0.06  S.E = 0.04  P = 0.08
#                innes: slope = 0.04  S.E = 0.02  P = 0.08
#w/ sup data:    daday: slope = 0.05  S.E = 0.01  P = 0.00
#                innes: slope = 0.04  S.E = 0.02  P = 0.08

#
interact_plot(EU_vs_NA_Actemp, mean_temp, Collector, mod2 = Origin,
                         mod2.values = c("Native", "Introduced"),
                         plot.points = TRUE, point.size = 2.5, point.shape = TRUE) +
  labs(x = expression(paste("Mean winter temperature "(degree*C))), 
       y = expression(paste(~italic(Ac)~ "frequency"))) +
  geom_line(aes(col = Collector)) +
  scale_color_manual(values = c("darkorange2", "darkorchid4")) +
  coord_cartesian(ylim = c(0, 1)) +
  ng2 



##################################################################################

#for Li
##################################################################################
EU_vs_NA_Litemp <- glm(freqLi~mean_temp + Collector + Origin + mean_temp:Collector + 
                       mean_temp:Origin + Collector:Origin + Collector:Origin:mean_temp, 
                       weights = N, data = EU_NA_df, family = "quasibinomial")

Anova(EU_vs_NA_Litemp, type = 2, test.statistic = "F")

#qualitatively same between sup and maintext models


sim_slopes(EU_vs_NA_Litemp, mean_temp, Collector, Origin, mod2.values = "Introduced")

sim_slopes(EU_vs_NA_Litemp, mean_temp, Collector, Origin, mod2.values = "Native")


#w/out sup data: daday: slope = 0.06  S.E = 0.03  P = 0.05
#                innes: slope = 0.06  S.E = 0.03  P = 0.03
#w/ sup data:    daday: slope = 0.03  S.E = 0.01  P = 0.03
#                innes: slope = 0.06  S.E = 0.03  P = 0.03

#
interact_plot(EU_vs_NA_Litemp, mean_temp, Collector, mod2 = Origin,
              mod2.values = c("Native", "Introduced"),
              plot.points = TRUE, point.size = 2.5, point.shape = TRUE) +
  labs(x = expression(paste("Mean winter temperature "(degree*C))), 
       y = expression(paste(~italic(Li)~ "frequency"))) +
  geom_line(aes(col = Collector)) +
  scale_color_manual(values = c("darkorange2", "darkorchid4")) +
  coord_cartesian(ylim = c(0, 1)) +
  ng2 



##################################################################################
##################################################################################

#################################################################################
#################################################################################
#################################################################################
# regular lms for average changes

# by lat
##################################################################################
# HCN
##################################################################################
lm_lat_HCN <- lm(freqHCN ~ abslat * Collector * Origin, data = all_freq_df, weights = N)

sim_slopes(lm_lat_HCN, abslat, Collector)
##################################################################################
# Ac
##################################################################################
lm_lat_Ac <- lm(freqAc ~ abslat * Collector * Origin, data = all_freq_df, weights = N)

sim_slopes(lm_lat_Ac, abslat, Collector)
##################################################################################
# Li
##################################################################################
lm_lat_Li <- lm(freqLi ~ abslat * Collector * Origin, data = all_freq_df, weights = N)

sim_slopes(lm_lat_Li, abslat, Collector)

##################################################################################
# by temp

upper_quart <- all_freq_df %>% 
  filter(lat >= 48.86)
##################################################################################
# HCN
##################################################################################
lm_temp_HCN <- lm(freqHCN ~ mean_temp * Collector * Origin, data = all_freq_df, weights = N)

sim_slopes(lm_temp_HCN, mean_temp, Collector, Origin, mod2.values = "Native")

sim_slopes(lm_temp_HCN, mean_temp, Collector, Origin, mod2.values = "Introduced")

sim_slopes(lm_temp_HCN, mean_temp, Collector) %>% tidy(.)

lm_temp_HCN %>% 
  emmeans(., tukey ~ Collector)

##################################################################################
# Ac
##################################################################################
lm_temp_Ac <- lm(freqAc ~ mean_temp * Collector * Origin, data = all_freq_df, weights = N)

sim_slopes(lm_temp_Ac, mean_temp, Collector, Origin, mod2.values = "Native")

sim_slopes(lm_temp_Ac, mean_temp, Collector, Origin, mod2.values = "Introduced")

sim_slopes(lm_temp_Ac, mean_temp, Collector) %>% tidy(.)

lm_temp_Ac %>% 
  emmeans(., tukey ~ Collector)

##################################################################################
# Li
##################################################################################
lm_temp_Li <- lm(freqLi ~ mean_temp * Collector * Origin, data = all_freq_df, weights = N)

sim_slopes(lm_temp_Li, mean_temp, Collector, Origin, mod2.values = "Native")

sim_slopes(lm_temp_Li, mean_temp, Collector, Origin, mod2.values = "Introduced")

sim_slopes(lm_temp_Li, mean_temp, Collector) %>% tidy(.)

lm_temp_Li %>%
  emmeans(., tukey ~ Collector)


################################################################################










####################################################################################
#supplementary Table S2
####################################################################################

#Table S2

Table_S1 <- all_freq_df %>% 
  rename("Mean winter temp." = "mean_temp") %>% 
  mutate("Time sampled" = if_else(Collector == "innes", "Contemporary", "Historical")) %>% 
  select(ID, Locality, lat, long, freqHCN, freqAc, freqLi, "Mean winter temp.", N,
         "Time sampled", Country, Continent, Hemisphere, "Data provenance", "Year of collection") %>% 
  arrange(`Time sampled`, Hemisphere, Continent, Country)

write_csv(Table_S1, path = "Table_S1.csv")
#################################################################################


#################################################################################
#evaluating differences in climatic variables through time
pre_anova <- function(var, df){
  
  var <- df %>% pull(var)
  lat <- df %>% pull("lat")
  time <- df %>% pull("Collector")
  range <- df %>% pull("Continent")
  
  lat <- abs(lat)
  
  mod <- lm(var ~ lat + time + lat:time)
  anova_mod <- Anova(mod, type = 3)
  
  print(AIC(mod))
  print(anova_mod)
  return(mod)
}

#################################################################################

#results of ANOVA + post hoc analysis by time for drought (aridity)
pre_anova("aridity", all_freq_df) %>% 
  emmeans(., tukey ~ time)

pre_anova("aridity2", all_freq_df) %>% 
  emmeans(., tukey ~ time)

#visualisation of aridity ~ latitude split by time
all_freq_df %>% 
  group_by(Collector) %>%
  ggplot(data = ., aes(x = abs(lat), y = aridity2, shape = Collector, color = Collector, linetype = Collector)) +
  geom_point(size = 2.5) +
  geom_smooth(method = 'lm', se =FALSE, color = 'black') +
  scale_color_manual(values = c("darkorange2", "darkorchid4")) +
  ng1

AI_HCN <- glm(freqHCN~aridity + Collector + Origin + aridity:Collector + aridity:Origin + 
              Collector:Origin + Collector:aridity:Origin, 
              data = all_freq_df, family = "quasibinomial")

Anova(AI_HCN, type = 3)

AI_Ac <- glm(freqAc~aridity + Collector + Origin + aridity:Collector + aridity:Origin + 
                Collector:Origin + Collector:aridity:Origin, 
              data = all_freq_df, family = "quasibinomial")

Anova(AI_Ac, type = 3)

AI_Li <- glm(freqLi~aridity + Collector + Origin + aridity:Collector + aridity:Origin + 
                Collector:Origin + Collector:aridity:Origin, 
              data = all_freq_df, family = "quasibinomial")

Anova(AI_Li, type = 3)


###################################################################################
the_d2anova <- function(var){
  
  var <- all_freq_df %>%  pull(var)
  AI <- all_freq_df %>% pull("aridity2")
  sample <- all_freq_df %>% pull("Collector")
  continent <- all_freq_df %>% pull("Continent")
  
  mod <- glm(var ~ AI + sample + continent + AI:sample + 
              AI:continent + sample:continent, family = "quasibinomial") 
  
  anova_mod <- Anova(mod, type = 3)
  
  print(AIC(mod)) 
  print(anova_mod)
  return(mod)
  
}

the_d2anova("freqHCN")

the_d2anova("freqAc")

the_d2anova("freqLi")







#################################################################################

