## Please find the R script and data utilised for the manuscript:
# Distinguishing the effects of fisheries and a warming climate on fish populations in Aotearoa, New Zealand.
## Charles P. Lavin, Daniel Pauly, Donna Dimarchopoulou, Cui Liang, Mark John Costello
## Corresponding author: Charles P. Lavin, charles.p.lavin@nord.no

### Included is the main R script to complete analyses and
### re-create figures (file: 'Lavin_et_al_NZ_MTC_MTL_code.R'),

### As well as data files required to run the R script, including:

#### Sea surface temperature anomaly data (file: 'NZ_EEZ_SSTA_1950_2020.csv')
##### extracted from https://psl.noaa.gov/data/gridded/data.kaplan_sst.html

#### The Sea Around Us fisheries catch data for New Zealand's EEZ (accessed 01 May 2022)
#### (files: 'SAU EEZ 554 v50-0.csv' & 'SAU EEZ 555 v50-0.csv')
##### This data is available at https://www.seaaroundus.org/

#### Species milieu classifications (Excel document), with
#### classifications derived from individual species' FishBase pages
#### (file: 'SAU_spp_milieu.csv')
##### see https://www.fishbase.se/search.php

#### Values calculated for the adjusted Mean Trophic Level (aMTL)
#### (file: 'aMTL_values.csv')
##### aMTL analysis was adapted from the region-based Marine Trophic Index (RMTI)
##### completed by Liang & Pauly (2017)
###### Liang C, Pauly D (2017) Fisheries impacts on China’s coastal ecosystems: Unmasking a pervasive ‘fishing down’ effect. PLoS One 12:e0173296
###### An online tool for completing the RMTI is available at:
###### http://www.seaaroundus.org/regional-mti-tools/
