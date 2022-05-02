# MSc_project

DATA AND FILE OVERVIEW


This README_ERC.txt file was generated on 2022-04-25 by Simon Innes.


GENERAL INFORMATION

1. Title of Dataset: “Data for: Evolution in response to climate in the native and introduced ranges of a globally distributed plant.”

2. Author Information
	Corresponding author
		Name: Simon Innes
		Institution: University of Louisiana at Lafayette
		Email: simon.innesz@gmail.com

	Co-author 
		Name: James Santangelo
		Institution: University of Toronto Mississauga

	Co-author 
		Name: Nicholas Kooyers
		Institution: University of Louisiana at Lafayette

	Co-author 
		Name: Kenneth Olsen
		Institution: Washington University in St-Louis
    
	Senior author 
		Name: Marc Johnson
		Institution: University of Toronto Mississauga

3. Date of data collection: 2008-2018 

4. Geographic location of data collection: see Table S1

5. Funding: JSS was funded by NSERC PGS-D and OGS. Funding for NJK was provided through a National Science Foundation grant (OIA-1920858). This work was additionally funded by an NSERC Discovery grant and a Canada Research Chair (CRC) to MTJJ.


DATA & FILE OVERVIEW

1. File list:
	File 1: “delta_matrix.csv”
	Description: Pairs of historical and contemporary populations within 50 km of each other. 

	File 3: “HCN_collection_date.csv”
	Description: Collection dates for each site in white clover's native range used to calculate growing degree days up until the date of sampling.

	File 4: “AllCities_AllPlants.csv”
	Description: Raw data for all individuals collected from Santangelo et al. 2020.
	
	File 5: “OK_TN_data.csv”
	Description: Population averages for sampling locations from Kooyers et al. 2014. 

	File 6: “herb_dmg.csv”
	Description: Raw data of estimated herbivory for each individual sampled from white clover's native range.

	File 7: “cleanDaday2.0data.csv”
	Description: Raw data for all individuals collected from the present study.

	File 8: “Kooyers_NA_data.csv”
	Description: Population averages for sampling locations from Kooyers et al. 2012.

	File 9: “Kooyers_NZ_data.csv”
	Description: Population averages for sampling locations from Kooyers and Olsen 2013.

	File 11: “Daday_data.csv”
	Description: Population averages for all historical data used in this study.
	
	File 12: “20170128 UAC_Distance.standardize.csv”
	Description: Distances from city centres used to standardize lengths of urban-rural transects from Santangelo et al. 2020.

	File 14: “Table_S2.csv”
	Description: Table of all historical and contemporary sampling locations with associated metadata. 


DATA-SPECIFIC INFORMATION FOR: "delta_matrix.csv"

1. Number of variables: 15

2. Number of cases/rows: 32

3. Variable List: 
	Contemporary_pop: Identification number of contemporary population or pairs of populations described
	Contemporary_HCN: Mean frequency of HCN production for individual or pair of contemporary populations
	Contemporary_Ac: Mean frequency of dominant Ac allele for individual or pair of contemporary populations
	Contemporary_Li: Mean frequency of dominant Li allele for individual or pair of contemporary populations
	Contemporary_MWT: Five year average of mean winter temperature (2013-2017) for each individual or pair of contemporary populations
	Historical_pop: Identification number of historical population described
	Historical_HCN: Mean frequency of HCN production for each historical population
	Historical_Ac: Mean frequency of dominant Ac allele for each historical population
	Historical_Li: Mean frequency of dominant Li allele for each historical population
	Historical_MWT: Five year average of mean winter temperature (1950-1954) for each historical population
	Delta_HCN: Change in the frequency of HCN production between historical and contemporary populations (Contemporary - historical)
	Delta_Ac: Change in the frequency of the dominant Ac allele between historical and contemporary populations (Contemporary - historical)
	Delta_Li: Change in the frequency of the dominant Li allele between historical and contemporary populations (Contemporary - historical)
	Delta_MWT: Change in in MWT between historical and contemporary populations (Contemporary - historical)
	Origin: Categorical descriptor of the origin of the population (native versus introduced; two levels)

4. Missing data codes: 
	none

5. Specialized formats or other abbreviations used: 
	HCN; Hydrogen cyanide
	MWT; Mean winter temperature


DATA-SPECIFIC INFORMATION FOR: "HCN_collection_date.csv"

1. Number of variables: 4

2. Number of cases/rows: 49

3. Variable List: 
	ID: population specific identification number
	country: Country of origin
	city: City of origin
	date_collected: Date the herbivory measurements were made (YYYY-MM-DD)

4. Missing data codes: 
	none

5. Specialized formats or other abbreviations used: 
	none


DATA-SPECIFIC INFORMATION FOR: "AllCities_AllPlants.csv"

1. Number of variables: 13

2. Number of cases/rows: 11,908

3. Variable List: 
	City: City of origin
	Population: Population specific identification number
	Transect: Letter descriptor of each transect in Toronto, ON
	Plant: Plant specific identification number
	Lat.pop: Latitudinal coordinates of population in decimal degrees
	Long.pop: Longitudinal coordinates of population in decimal degrees
	HCN_Result: Presence or absence of HCN production (1 = present; 0 = absent)
	Locus.Li: Presence or absence of dominant Li allele (1 = present; 0 = absent)
	Locus.Ac: Presence or absence of dominant Ac allele (1 = present; 0 = absent)
	Dmg.1: Estimate of herbivory damage on first leaf (not sued in this study)
	Dmg.2: Estimate of herbivory damage on second leaf (not used in this study)
	Dmg.avg: Average of estimated herbivory damage on both leaves (not used in this study)
	Distance: Distance from urban centre for each population in km
	
4. Missing data codes: 
	NA

5. Specialized formats or other abbreviations used: 
	HCN; Hydrogen cyanide


DATA-SPECIFIC INFORMATION FOR: "OK_TN_data.csv"

1. Number of variables: 15

2. Number of cases/rows: 16

3. Variable List: 
	Population: City of origin
	Subpopulation Name: Specific subpopulation name
	Latitude: Latitudinal coordinates of population in decimal degrees
	Longitude: Longitudinal coordinates of population in decimal degrees
	Annual Aridity Index: Aridity index for each population (not used in this study)
	Annual Precipitation (mm): Annual precipitation in mm for each population (not used in this study)
	Minimum Winter Temperature: Minimum winter temperature for each population (not used in this study)
	Total Plants: Total number of plants in each subpopulation
	AcLi: Total number of plants capable of producing HCN
	Acli: Total number of plants with a dominant Ac allele
	acLi: Total number of plants with a dominant Li allele
	acli: Total number of plants homozygous recessive at Ac and Li
	Freq.(AcLi): Frequency of HCN production in each subpopulation
	Freq.(Ac): Frequency of dominant Ac allele in each subpopulation
	Freq.(Li): Frequency of dominant Li allele in each subpopulation

4. Missing data codes: 
	none

5. Specialized formats or other abbreviations used: 
	HCN; Hydrogen cyanide


DATA-SPECIFIC INFORMATION FOR: "herb_dmg.csv"

1. Number of variables: 8

2. Number of cases/rows: 2073

3. Variable List: 
	ID: Population specific identification number
	country: Country of origin
	city: City of origin
	pop.ind: Individual specific identification number
	herb.dmg: Estimated percentage of herbivory to the nearest %
	u_r: Categorical descriptor of urban or rural habitat (two levels)
	lat: Latitudinal coordinates of population in decimal degrees
	long: Longitudinal coordinates of population in decimal degrees

4. Missing data codes: 
	none

5. Specialized formats or other abbreviations used: 
	none


DATA-SPECIFIC INFORMATION FOR: "cleanDaday2.0data.csv"

1. Number of variables: 18

2. Number of cases/rows: 3673

3. Variable List: 
	ID: Population specific identification number
	country: Country of origin
	city: City of origin
	u_r: Categorical descriptor of urban or rural habitat (two levels)
	pop: Population number
	ind: Individual number
	lat: Latitudinal coordinates of population in decimal degrees
	long: Longitudinal coordinates of population in decimal degrees
	alt: Elevation of population in m
	hcn: Presence or absence of HCN production (1 = present; 0 = absent)
	Ac: Presence or absence of dominant Ac allele (1 = present; 0 = absent)
	Li: Presence or absence of dominant Li allele (1 = present; 0 = absent)
	Collector:

4. Missing data codes: 
	NA

5. Specialized formats or other abbreviations used: 
	HCN; Hydrogen cyanide


DATA-SPECIFIC INFORMATION FOR: "kooyers_NA_data.csv"

1. Number of variables: 12

2. Number of cases/rows: 127

3. Variable List: 
	Pop #: Population specific identifier
	Transect: Transect name
	Population: Population name
	Latitude: Latitudinal coordinates of population in decimal degrees
	Longitude: Longitudinal coordinates of population in decimal degrees
	Minimum Winter Temperature: Minimum winter temperature for each population (not used in this study)
	Elevation: Elevation of population in m
	N: Total number of plants in each population
	Freq.(AcLi): Frequency of HCN production in each population
	Freq.(Ac): Frequency of dominant Ac allele in each population
	Freq.(Li): Frequency of dominant Li allele in each population
	N: Number of individuals of individuals used for subsequent analysis (not used in this study)

4. Missing data codes: 
	none

5. Specialized formats or other abbreviations used: 
	HCN; Hydrogen cyanide


DATA-SPECIFIC INFORMATION FOR: "kooyers_NZ_data.csv"

1. Number of variables: 12

2. Number of cases/rows: 15

3. Variable List: 
	Pop #: Population specific identifier
	Transect: Transect name
	Population: Population name
	Latitude: Latitudinal coordinates of population in decimal degrees
	Longitude: Longitudinal coordinates of population in decimal degrees
	Minimum Winter Temperature: Minimum winter temperature for each population (not used in this study)
	Elevation: Elevation of population in m
	N: Total number of plants in each population
	Freq.(AcLi): Frequency of HCN production in each population
	Freq.(Ac): Frequency of dominant Ac allele in each population
	Freq.(Li): Frequency of dominant Li allele in each population
	N: Number of individuals of individuals used for subsequent analysis (not used in this study)

4. Missing data codes: 
	none

5. Specialized formats or other abbreviations used: 
	HCN; Hydrogen cyanide


DATA-SPECIFIC INFORMATION FOR: "Daday_data.csv"

1. Number of variables: 28

2. Number of cases/rows: 143

3. Variable List: 
	Locality: Name of locality for each population
	Country: Country of origin
	lat: Latitudinal coordinates of population in decimal degrees
	long: Longitudinal coordinates of population in decimal degrees
	alt: Elevation of population in m
	AcLi: Frequency of HCN production within populations (0-100)
	Acli: Frequency of dominant Ac allele within populations (0-100)
	acLi: Frequency of dominant Li allele within populations (0-100)
	acli: Frequency of homozygous recessives within populations (0-100)
	no_AcLi: Number of individuals capable of HCN production within populations
	no_Acli: Number of individuals with the dominant Ac allele within populations
	no_acLi: Number of individuals with the dominant Li allele within populations
	no_acli: Number of homozygous recessive individuals within populations
	no_plants: Number of plants within each population
	AcAc: Frequency of homozygous dominant individuals at Ac within populations
	Acac: Frequency of heterozygous individuals at Ac within populations
	acac: Frequency of homozygous recessive individuals at Ac within populations
	LiLi: Frequency of homozygous dominant individuals at Li within populations
	Lili: Frequency of heterozygous individuals at Li within populations
	lili: Frequency of homozygous recessive individuals at li within populations
	freqAc: Frequency of dominant Ac allele within populations
	freqLi: Frequency of dominant Li allele within populations
	freqHCN: Frequency of HCN production within populations
	Continent: Continent of origin
	Collector: Categorical variable describing historical or contemporary sampling (innes = contemporary; daday = historical; two level)
	ID: Population specific identification number
	
4. Missing data codes: 
	NA

5. Specialized formats or other abbreviations used: 
	HCN; Hydrogen cyanide
	

DATA-SPECIFIC INFORMATION FOR: "20170128 UAC_Distance.standardize.csv"

1. Number of variables: 10

2. Number of cases/rows: 19

3. Variable List: 
	City: City of origin
	Center.Lat: Latitudinal coordinates of city centre in decimal degrees
	Center.Long: Longitudinal coordinates of city centre in decimal degrees
	Center.lat.rad: Latitudinal coordinates of city centre in radians
	Center.long.rad: Longitudinal coordinates of city centre in radians
	Border.Lat: Latitudinal coordinates of urban border
	Border.Long: Longitudinal coordinates of urban border
	Border.lat.rad: Latitudinal coordinates of urban border in radians
	Border.long.rad: Longitudinal coordinates of urban border in radians
	Distance: Distance from urban centre to urban border in km

4. Missing data codes: 
	NA

5. Specialized formats or other abbreviations used: 
	HCN; Hydrogen cyanide
	

DATA-SPECIFIC INFORMATION FOR: "Table_S2.csv"

1. Number of variables: 15

2. Number of cases/rows: 220

3. Variable List: 
	ID: Population specific identification number
	Locality: Name of locality for each population
	lat: Latitudinal coordinates of population in decimal degrees
	long: Longitudinal coordinates of population in decimal degrees
	freqHCN: Frequency of HCN production within populations
	freqAc: Frequency of dominant Ac allele within populations
	freqLi: Frequency of dominant Li allele within populations
	Mean winter temp.: Five year average of mean winter temperature (1950-1954 or 2013-2017 for historical and contemporary populations, respectively) for each population
	N: Number of individuals in each population
	Time sampled: Categorical variable describing historical or contemporary sampling (two levels)
	Country: Country of origin
	Continent: Continent of origin
	Hemisphere: Hemisphere of origin
	Data provenance: Origin of the data
	Year of collection: Year the data were collected
	
4. Missing data codes: 
	NA

5. Specialized formats or other abbreviations used: 
	HCN; Hydrogen cyanide


