install.packages("ncdf4")

??ncdf4
install.packages("cruts")

install.packages("tidyverse")

install.packages("raster")

library(ncdf4)
library(cruts)
library(tidyverse)

wd <- "cru_ts4.02.1901.2017.tmp.dat.nc"

tmp <- nc_open(wd, write=FALSE, readunlim=TRUE, verbose=FALSE, auto_GMT=TRUE, suppress_dimvals=FALSE )
print(tmp) #name file



?ncvar_get

lon <- ncvar_get(tmp,"lon")
nlon <- dim(lon)
head(lon) #longitudes

lat <- ncvar_get(tmp,"lat",verbose=F)
nlat <- dim(lat)
head(lat) #latitudes
print(c(nlon,nlat))

t <- ncvar_get(tmp,"time")
tunits <- ncatt_get(tmp,"time","units")
nt <- dim(t) #time

dname <- "tmp"
tmp_array <- ncvar_get(tmp,dname)
dlname <- ncatt_get(tmp,dname,"long_name")
dunits <- ncatt_get(tmp,dname,"units")
fillvalue <- ncatt_get(tmp,dname,"_FillValue") #building array
dim(tmp_array)

#tmp_array[, , 1]

#
index_vector <- seq(from = 516, to = 660, by = 1) #each item in list becomes single object
tmp_array_lst <- list() #list for output arrays

##https://www.r-bloggers.com/how-to-write-the-first-for-loop-in-r/
for (i in 1:length(index_vector)){
  
  index = index_vector[i]
  # print(index)
  tmp_array_lst[[i]] <- tmp_array[, , index]
  # print(tmp_array[, , i])
} ##for taking slices of array

#?lapply




#figuring out days since
monnb <- function(d) { lt <- as.POSIXlt(as.Date(d, origin="1900-01-01")); 
lt$year*12 + lt$mon }
mondf <- function(d1, d2) { monnb(d2) - monnb(d1) }
##https://stackoverflow.com/questions/1995933/number-of-months-between-two-dates

mondf("1900-1-1", "1943-1-1") # =516
mondf("1900-1-1", "1955-1-1") # =660


#tidy verse

list.files()
HCN_data <- readr::read_csv("cleanDaday2.0data.csv")
str(HCN_data) # Look at dataset in more detail

Acyan <- HCN_data %>%
  dplyr::filter(hcn == 0) %>%
  select(ID, ind, hcn) #adds coloumns ind + ID + hcn

write_csv(Acyan, path = "test/Acyan_only.csv")
