fishbase2
================
Han lab
8/21/2020

\#\#\#\#\#install packages

``` r
testing_var = 1#0: testing out, not doing as full grid search
nruns = 10
date = "20200821"
time = "1149"
output_name = paste("vert", date, time, sep = "_")
save(output_name, file = "output_name.Rdata")
print(Sys.time())
```

    ## [1] "2020-08-21 11:51:11 EDT"

``` r
id_field = "Species"
rm_list = c("nchar", "Order", "nchar", "haddock_score_sd")
if (testing_var == 1){
  eta = c(0.0001, 0.001, 0.01, 0.1)
  max_depth = c(2,3,4)
  n.minobsinnode = c(2,5)
  nrounds = 100000
} else {
  eta = c(0.0001, 0.001)
  max_depth = c(2)
  n.minobsinnode = c(2)
  nrounds = 100000
}

do_haddock_infected = 1#make 1 to do analysis with over/under infectable species w/ highest haddock score
save(do_haddock_infected, file = "do_haddock_infected.Rdata")

do_score = 1#1 to do regression analysis with haddock score
save(do_score, file = "do_score.Rdata")

do_mammal = 1#whether to do separate analysis for mammals (excluding other verts)
save(do_mammal, file = "do_mammal.Rdata")

do_haddock_median = 0
save(do_haddock_median, file = "do_haddock_median.Rdata")

make_v = 0 #whether to remake V, which takes a little time to get all the data from fishbase

performance_out = NULL#this will store performance metrics
DF_merge_out = NULL#this will have the mean predictions

ACE2_file_name = "ACE2_sequences_fixed_20200820.csv"
```

\#\#cores

\#\#set up function gridSearch.R

\#\#try out grid search using Pima dataset
\#\#<https://www.kaggle.com/kumargh/pimaindiansdiabetescsv?select=pima-indians-diabetes.csv>
\#\#label = “X1”

``` r
DF <- read.csv("../../functions/pima-indians-diabetes.csv")
DF$id_field_vals = seq(1,dim(DF)[1])
load("bootstrapGBM.Rdata")
  nrounds_test = 2000
  min_trees = 0
  k_split = 0.8
  label = "X1"
vars_pred = setdiff(names(DF), label)
OUT <-  bootstrapGBM(DF, label = "X1", vars_pred, k_split, distribution = "bernoulli", eta = 0.5, max_depth = 1, nrounds = nrounds_test, nruns =2, bootstrap = "observed", method = "cv", cv.folds = 5,
                         n.minobsinnode=2,file_label = "test", id_field= "id_field_vals")

temp = OUT[[3]]
```

``` r
source("bootstrapGBM.R")
```

\#\#function to make binary fields into factor

``` r
binary_factor <- function(DF){
  T = DF
  names = names(T)
  binary = NULL
  for (a in 1:length(names)){
    vals = unique(T[,names[a]])
    vals = vals[!is.na(vals)]
    if (length(which(vals==0)) + length(which(vals == -1)) == 2){
      T[,names[a]]=factor(T[,names[a]])
      binary = c(binary, names[a])
      #change to 1 and 0
    }
  }
  T
}
save(binary_factor, file = "binary_factor.Rdata")
```

\#\#function to take the same across rows of categorical variables that
have been 1/0 encoded, where a species may have 1 for more than one
condition of a variable

\#\#function to replace NAs with real values for binary fields

\#\#settings

\#\#look at docs about tables available from fishbase

\#\#read in data and fix species names

\#\#distribution \#\#currently this is ~ FAO areas table (minus “note”
field) e.g. <http://www.fishbase.us/Country/FaoAreaList.php?ID=5537>
\#\#each species may have multiple bounding boxes

Read in the FAO areas (from
<http://www.fao.org/geonetwork/srv/en/main.home?uuid=ac02a460-da52-11dc-9d70-0017f293bd28>
as described by
<http://www.fishbase.us/manual/English/FishbaseThe_FAOAREAS_Table.htm>).
It looks like our data contain both the inland and marine FAOs, so I
read in both and combined them according to a single column of FAO code.

\#\#check out some tables in fishbase \#\#brains: one entry for each
individual fish: BrainWeight, BodyWeight
\#\#<https://www.fishbase.in/manual/fishbasethe_brains_table.htm>

\#\#country: multiple rows per species; for
example:

## countrysub – multiple rows per species

\#\#<https://www.fishbase.de/manual/english/FishBaseThe_Countries_Table.htm>

\#\#get ecology data
\#\#<http://fishbase.us/manual/English/FishbaseThe_ECOLOGY_Table.htm>

\#\#distribution \#\#currently this is ~ FAO areas table (minus “note”
field) e.g. <http://www.fishbase.us/Country/FaoAreaList.php?ID=5537>
\#\#each species may have multiple bounding boxes

\#\#ecosystem – couldn’t find description of this online \#\#multiple
rows per species, one for each ecosystem

\#\#estimate: a table of estimates from some models on trophic levels
\#\#<http://www.fishbase.us/manual/English/FishbaseThe_FOOD_ITEMS_table.htm>

\#\#faoareas, seems to be redundant to countrysub?

\#\#fecundity \#\#sometimes multiple rows per species. could not
\#\#could not locate doc table about fecundity. spawning table seems to
be something different (different fields):
<https://www.fishbase.in/manual/fishbasethe_spawning_table.htm>

\#\#fooditems – including this one
\#\#<http://www.fishbase.org/manual/english/fishbasethe_food_items_table.htm>
\#\#multiple rows per species, for different food types, life stages of
predator, locality, etc.

\#\#genetic – don’t think we want to use this, but including just to see
what it shows

\#\#introductions – species introductions data. for now making one new
feature: the number of records about introductions; it seems that each
row is a different place
\#\#<https://www.fishbase.in/manual/fishbasethe_introduction_table.htm>

\#\#larvae
\#\#<https://www.fishbase.in/manual/fishbasethe_larvae_table.htm> \#\#2
out of the 74 species have multiple records w/ different values.
excluding for now.

\#\#length\_freq; multiple records for some species; excluding for now;
could not find metadata

\#\#length\_length: conversion of length types

\#\#length\_weight: The LENGTH-WEIGHT table presents the a and b values
of over 5,000 length-weight relationships of the form W = a x Lb,
pertaining to about over 2,000 fish species. \#\#multiple records for
some species. \#\#seems like this may only be useful in combination with
length\_length
\#\#<https://www.fishbase.de/manual/FishbaseThe_LENGTH_WEIGHT_Table.htm>

\#\#maturity \#\#multiple records for some species, would need to take
averages if we wanted to use. there are multiple measures of maturity to
choose from.
\#\#<https://www.fishbase.in/manual/fishbasethe_maturity_table.htm>

\#\#morphology
\#\#<https://www.fishbase.in/manual/fishbasethe_morphology_table.htm>
\#\#there are multiple records for some species.

\#\#morphometrics \#\#there are multiple records for some species; to
include we would need to take averages \#\#exclude for now because
couldn’t find documentation

\#\#oxygen
\#\#<https://www.fishbase.in/manual/fishbasethe_oxygen_table.htm>
\#\#there are multiple records for some species (e.g. for different
sexes); to include we would need to take averages \#\#include along with
potentially influencing variables – e.g. salinity, temp, swimming speed,
etc.

\#\#popchar: Table of maximum length (Lmax), weight (Wmax) and age
(tmax)
\#\#<https://www.fishbase.in/manual/fishbasethe_popchar_table.htm>
\#\#there are multiple records for some species; to include we would
need to take averages \#\#

\#\#popgrowth
\#\#<https://www.fishbase.in/manual/fishbasethe_popgrowth_table.htm>
\#\#multiple records for some species, e.g. for different sexes;

\#\#popqb
\#\#<https://www.fishbase.se/manual/english/fishbasethe_popqb_table.htm>
\#\#population-based estimates of food consumption (i.e., estimates that
account for the age structure of populations) \#\#multiple responses for
some species. here there are two measures, popqb and maintenance qb.

\#\#predators
\#\#<https://www.fishbase.se/manual/English/fishbasethe_predators_table.htm>

\#\#ration \#\#�ration� (Rd) pertains to an estimate of daily food
consumption by fish of a specific size
\#\#<https://www.fishbase.in/manual/fishbasethe_ration_table.htm>
\#\#multiple rows for some species

\#\#reproduction
\#\#<https://www.fishbase.in/manual/fishbasethe_reproduction_table.htm>
\#\#only one row per species for these HADDOCK species; adding these
fields

\#\#spawning
\#\#<https://www.fishbase.in/manual/fishbasethe_spawning_table.htm>
\#\#multiple rows per species, for different localities

\#\#speed
\#\#<https://www.fishbase.se/manual/English/PDF/FB_Book_ATorres_Swimming_Speed_RF_JG.pdf>
\#\#<https://www.fishbase.in/manual/fishbasethe_swimming_and_speed_tables.htm>
\#\#multiple records for some species

\#\#stocks
\#\#<https://www.fishbase.in/manual/fishbasethe_stocks_table.htm>
\#\#multiple records for some species, one for each stock

\#\#diet \#\#<https://www.fishbase.in/manual/fishbasethe_diet_table.htm>
\#\#has multiple rows for different stages

\#\#diet\_items – multiple rows per species. seems to be linked with
DietCode to diet table
\#\#<https://www.fishbase.se/manual/English/fishbasethe_food_items_table.htm>

\#\#swimming
\#\#<https://www.fishbase.in/manual/fishbasethe_swimming_and_speed_tables.htm>
\#\#one record per species

\#\#see what coverage is

\#\#remove fields with 0 coverage

\#\#add back haddock fields

\#\#remove fields with near-zero variation

\#\#look for fields in common with other taxa that are not fish and
output to add to datasets from other verts

\#\#add field with AA position 30

\#\#add AA value to rest of fishbase data

\#\#remove fields with near-zero variation again

\#\#combine data Adrian made with rest of fields from fish. remove Order
field

``` r
load("V.Rdata")
alt <- V
 #this is to create development_d
nrows = dim(alt)[1]
development_d = rep(NA, dim(alt)[1])
for (a in 1:nrows){
  tmp = sum(alt$gestation_d[a],alt$incubation_d[a], na.rm =TRUE)
  if (tmp == 0){
    development_d[a] = NA
  } else {
    development_d[a] = tmp
  }
  print(development_d[a])
}
```

    ## [1] 178.12
    ## [1] 65
    ## [1] 70
    ## [1] 26.5
    ## [1] 87
    ## [1] NA
    ## [1] 20
    ## [1] NA
    ## [1] 63.5
    ## [1] 79
    ## [1] NA
    ## [1] 28
    ## [1] 24.5
    ## [1] 618.6127
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] 554
    ## [1] 625.8325
    ## [1] NA
    ## [1] NA
    ## [1] 288
    ## [1] 494.7493
    ## [1] 16.25
    ## [1] NA
    ## [1] 791.495
    ## [1] 800.755
    ## [1] NA
    ## [1] 124.1969
    ## [1] 124.1969
    ## [1] 310.5
    ## [1] 25.5
    ## [1] 358
    ## [1] 12.5
    ## [1] 135.495
    ## [1] 324
    ## [1] 1030
    ## [1] 167
    ## [1] 20
    ## [1] 24.55
    ## [1] 118
    ## [1] 112
    ## [1] 223.47
    ## [1] NA
    ## [1] 23.5
    ## [1] 329.99
    ## [1] 167.99
    ## [1] 68.5
    ## [1] NA
    ## [1] 80
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] 17
    ## [1] NA
    ## [1] 87.5
    ## [1] 220
    ## [1] 12.5
    ## [1] 13.5
    ## [1] 266
    ## [1] 825.0598
    ## [1] 378
    ## [1] 58.91082
    ## [1] 56
    ## [1] 112.8284
    ## [1] NA
    ## [1] 13.9
    ## [1] 356.0956
    ## [1] 78.57
    ## [1] 720.6595
    ## [1] 674.985
    ## [1] NA
    ## [1] 66.06754
    ## [1] 582.05
    ## [1] 27
    ## [1] 29
    ## [1] 129
    ## [1] NA
    ## [1] 185.8298
    ## [1] 49
    ## [1] 20.5
    ## [1] 27.5
    ## [1] 150
    ## [1] NA
    ## [1] NA
    ## [1] 903.4514
    ## [1] NA
    ## [1] 513
    ## [1] 48
    ## [1] 35
    ## [1] 140
    ## [1] NA
    ## [1] 560
    ## [1] 460
    ## [1] 55.75
    ## [1] 76
    ## [1] 639.01
    ## [1] NA
    ## [1] NA
    ## [1] 638.474
    ## [1] 13.5
    ## [1] 1339.68
    ## [1] 126
    ## [1] 135.8952
    ## [1] 330
    ## [1] 330
    ## [1] 344
    ## [1] 18.5
    ## [1] 358.22
    ## [1] NA
    ## [1] 60.195
    ## [1] 70.41506
    ## [1] NA
    ## [1] 18
    ## [1] NA
    ## [1] NA
    ## [1] 32
    ## [1] NA
    ## [1] 45.52264
    ## [1] 113.44
    ## [1] 29.775
    ## [1] 900
    ## [1] 420
    ## [1] 40.6
    ## [1] 38.6
    ## [1] NA
    ## [1] 86.26174
    ## [1] 83.49
    ## [1] NA
    ## [1] NA
    ## [1] 112.395
    ## [1] 29.18
    ## [1] NA
    ## [1] 679.15
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] 207.5911
    ## [1] NA
    ## [1] NA
    ## [1] 27.5
    ## [1] 122.84
    ## [1] 60.436
    ## [1] 180.0906
    ## [1] 399.39
    ## [1] 28
    ## [1] 889.721
    ## [1] 31.74
    ## [1] 447.05
    ## [1] 60.59594
    ## [1] 260.7443
    ## [1] 304.8944
    ## [1] NA
    ## [1] 467.3204
    ## [1] 458
    ## [1] 194
    ## [1] 210.19
    ## [1] 178.96
    ## [1] 13.5
    ## [1] 30.5
    ## [1] 61.5
    ## [1] 51.26188
    ## [1] 50
    ## [1] 41
    ## [1] 29
    ## [1] 62
    ## [1] 24
    ## [1] 20
    ## [1] 210
    ## [1] 464.62
    ## [1] 171.3547
    ## [1] 71.48
    ## [1] NA
    ## [1] 75
    ## [1] 259.42
    ## [1] 126.9143
    ## [1] 145.5246
    ## [1] NA
    ## [1] 91
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] 184.3
    ## [1] 34
    ## [1] 154.82
    ## [1] 42.60956
    ## [1] NA
    ## [1] NA
    ## [1] 165.62
    ## [1] 119.61
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] 123.54
    ## [1] 398.34
    ## [1] 232
    ## [1] 233.33
    ## [1] 320.545
    ## [1] 316
    ## [1] 52
    ## [1] 12.5
    ## [1] 42.00042
    ## [1] 30
    ## [1] NA
    ## [1] 12
    ## [1] 154
    ## [1] 230.2
    ## [1] 14.15
    ## [1] NA
    ## [1] NA
    ## [1] 340.4166
    ## [1] NA
    ## [1] 669.58
    ## [1] 92
    ## [1] 730
    ## [1] 31.9
    ## [1] 25.16
    ## [1] 206
    ## [1] 212.4159
    ## [1] 345
    ## [1] 55.5
    ## [1] 104.3525
    ## [1] NA
    ## [1] 593.84
    ## [1] 13
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] 18
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] 16
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA
    ## [1] NA

``` r
alt$development_d = development_d

alt[which(is.na(alt$longevity_y)) %in% which(!is.na(alt$maximum_longevity_y)), "longevity_y"] <- alt$maximum_longevity_y[which(is.na(alt$longevity_y)) %in% which(!is.na(alt$maximum_longevity_y))] #this is to add max_longevity to places where longevity doesn't exist. I can't test it with the dataset that I have up because all the mammals I'm working with right now either have longevity or don't have both
rm = c("maximum_longevity_y", "incubation_d", "gestation_d")
keep = setdiff(names(alt), rm)
alt = alt[,keep]
V <- alt
save(alt, file = "alt.Rdata")
save(V, file = "V.Rdata")
```

\#add field mass\_specific\_production and remove
“major\_habitat\_type\_breadth” (because it seems redundant to tnc
ecoregion breadth)

\#\#add AA value to rest of vert data

\#\#add WOS hits from R package wosr

``` r
# if (make_v == 1){
  load("V.Rdata")
  dim(V)
```

    ## [1] 278  38

``` r
  W <- read.csv("wos_species_hits.csv")#this version has results only from exact matches
  names(W)[names(W)=="query"]="Species"
  names(W)[names(W)=="rec_cnt"]="WOS_hits"
  
  W$Species = str_replace(W$Species, "TS", "")
  W$Species = str_replace(W$Species, " = ", "")
  W$Species = str_sub(W$Species, 2, -2)
  
  W2 <- read.csv("query_hits_synonyms.csv")#this includes synonyms
  names(W2)[names(W2)=="query"]="Species"
  #add small value to rec_cnt so we can take log
  W2$rec_cnt = W2$rec_cnt+0.0001
  W2$rec_cnt = log(W2$rec_cnt)
  names(W2)[names(W2)=="rec_cnt"]="log_WOS_hits_synonyms"
  
  a = 18
  # hits_synonyms = rep(NA, dim(W2)[1])
  out = NULL
  species = unique(W$Species)
  for (a in 1:length(species)){
    ind_W2 = str_which(W2$Species, species[a])#find the index in W2 that includes the a'th species
    tmp = data.frame(log_WOS_hits_synonyms = sum(W2$log_WOS_hits_synonyms[ind_W2]),
                     Species = species[a]) 
    # print("a")
    # print(a)
    # print(dim(tmp)[1])
    out = rbind(out, tmp)
  }
  V = merge(V, out)
  dim(V)
```

    ## [1] 278  39

``` r
  save(V, file = "V.Rdata")
  table(V$Class)
```

    ## 
    ## Actinopterygii       Amphibia           Aves Elasmobranchii    Holocephali 
    ##             70              4             59              1              1 
    ##       Mammalia       Reptilia 
    ##            126             17

``` r
# }
```

\#\#output: “haddock\_vert\_for\_gbm.csv”

``` r
print(Sys.time())
```

    ## [1] "2020-08-21 11:51:31 EDT"

``` r
load("gridSearch.Rdata")
# save(output_name, file = "output_name.Rdata")
load("V.Rdata")
V$adult_svl_cm[is.nan(V$adult_svl_cm)] <- NA
V$log_adult_body_mass_g[is.nan(V$log_adult_body_mass_g)] <- NA

DF = V
names_df = names(DF)
#remove variables with near zero variation
sp_ind = which(names(DF)=="Species")
nzv = nearZeroVar(DF, freqCut = 95/5, saveMetrics = TRUE)
okay_inds = which(nzv$nzv == FALSE)

near_zero_inds =which(nzv$nzv == TRUE)
print("near zero variation fields")
```

    ## [1] "near zero variation fields"

``` r
print(names_df[near_zero_inds])
```

    ## character(0)

``` r
DF = DF[,c(okay_inds, near_zero_inds)]#include only the columns that have variation; also include adult_svl_cm because it looks like it does have variation

out = DF#now 1/0 encode, to keep every Class
Species = out$Species
sp_ind = which(names(out)=="Species")
dmy <- dummyVars(" ~ .", data = out[,-sp_ind])
out <- data.frame(predict(dmy, newdata = out))
out$Species = Species
DF = out

load("binary_factor.Rdata")
DF = binary_factor(DF)

#find out what haddock_score_median is across all species
haddock_median = median(DF$haddock_score_mean)

haddock_median_and_below = rep(0, dim(DF)[1])
inds = which(DF$haddock_score_mean <= haddock_median)
haddock_median_and_below[inds]= 1
DF$haddock_median_and_below = haddock_median_and_below
DF_gbm = DF
rm = c("haddock_score_sd", "nchar")
keep = setdiff(names(DF_gbm), rm)
DF_gbm = DF_gbm[,keep]
save(DF_gbm, file = "DF_gbm.Rdata")

rm =  c("haddock_score_mean")
keep = setdiff(names(DF),rm)
DF = DF[,keep]
write.csv(DF, file = "haddock_vert_for_gbm.csv", row.names = FALSE)
names(DF)
```

    ##  [1] "ClassActinopterygii"              "ClassAmphibia"                   
    ##  [3] "ClassAves"                        "ClassElasmobranchii"             
    ##  [5] "ClassHolocephali"                 "ClassMammalia"                   
    ##  [7] "ClassReptilia"                    "nchar"                           
    ##  [9] "haddock_score_sd"                 "ForStrat.ground"                 
    ## [11] "ForStrat.understory"              "ForStrat.arboreal"               
    ## [13] "ForStrat.aerial"                  "ForStrat.marine"                 
    ## [15] "Activity.Nocturnal"               "Activity.Crepuscular"            
    ## [17] "Activity.Diurnal"                 "female_maturity_d"               
    ## [19] "male_maturity_d"                  "weaning_d"                       
    ## [21] "log_litterclutch_size_n"          "litters_or_clutches_per_y"       
    ## [23] "log_inter_litterbirth_interval_y" "log_birthhatching_weight_g"      
    ## [25] "log_weaning_weight_g"             "log_adult_body_mass_g"           
    ## [27] "infantMortalityRate_per_year"     "mortalityRateDoublingTime_y"     
    ## [29] "metabolicRate_W"                  "temperature_K"                   
    ## [31] "longevity_y"                      "log_female_body_mass_g"          
    ## [33] "log_male_body_mass_g"             "log_no_sex_body_mass"            
    ## [35] "adult_svl_cm"                     "diet_breadth"                    
    ## [37] "major_habitat_type_breadth"       "tnc_ecoregion_breadth"           
    ## [39] "log_range_size"                   "development_d"                   
    ## [41] "AA_83_Y"                          "AA_30_negative"                  
    ## [43] "log_WOS_hits_synonyms"            "Species"                         
    ## [45] "haddock_median_and_below"

\#find out if score is above or below species that’s been infected that
has highest haddock score

``` r
load("binary_factor.Rdata")
load("DF_gbm.Rdata")
DF = DF_gbm
#read in file where we have added real-world outcomes to supplementary data from Damas et al. 2020

R <- read.csv("Damas_validation_ACE2_species - Data.csv")

R = subset(R, !is.na(transmission.to.conspecifics) )
R = R[,c("transmission.to.conspecifics", "Species")]
R$Species[R$Species == "Peromyscus maniculatus bairdii" ]= "Peromyscus maniculatus"
haddock = rep(NA, dim(R)[1])
species = unique(R$Species)
for (a in 1:length(species)){
  tmp = subset(DF, as.character(Species) == as.character(species[a]))
  haddock[a]= tmp$haddock_score_mean
}  
R$haddock = haddock
R_highest = subset(R, haddock == max(R$haddock))
print("infected species w/ conspecific transmission and highest haddock score")
```

    ## [1] "infected species w/ conspecific transmission and highest haddock score"

``` r
R_highest
```

    ##     transmission.to.conspecifics              Species   haddock
    ## 362                            1 Mesocricetus auratus -129.4977

``` r
haddock_cutoff = R_highest$haddock

haddock_infected_and_below = rep(0, dim(DF)[1])#default is that you have a haddock score that is above species infectable and w/ highest score

inds_equal_and_lower = which(DF$haddock_score_mean <= haddock_cutoff)
haddock_infected_and_below[inds_equal_and_lower] = 1
DF$haddock_infected_and_below = haddock_infected_and_below

rm = which(names(DF) %in% c( "Order",  "nchar", "haddock_score_sd"))

if (length(rm)>0){
  DF = DF[,-rm]
}
DF_gbm = DF
save(DF_gbm, file = "DF_gbm.Rdata")

keep = setdiff(names(DF), "haddock_median_and_below")
keep = setdiff(keep, "haddock_score_mean")

DF = DF[,keep]#remove just for making this dataset for infected

write.csv(DF, file = "haddock_vert_infected_for_gbm.csv", row.names = FALSE)
```

``` r
if (do_haddock_median == 1){
  load("output_name.Rdata")
  label = "haddock_median_and_below"
  output_name = paste(output_name, label)
}
```

\#\#gridsearch w/ median

``` r
if (do_haddock_median == 1){

  load("DF_gbm.Rdata")
  DF = DF_gbm
  rm = which(names(DF) %in% c("haddock_score_mean", "Order", "nchar", "haddock_score_sd", "haddock_infected_and_below"))
  if (length(rm)>0){
    DF = DF[,-rm]
  }
  
  # n.minobsinnode = c(2)
  k_split = 0.8
  distribution = "bernoulli"
  
  label_col_ind = which(names(DF)==label)
  x_col = seq(1:dim(DF)[2])
  x_col = setdiff(x_col, label_col_ind)
  vars = colnames(DF)[x_col]
  vars = setdiff(vars, id_field)#remove species
  
  GRID <- gridSearch(DF = DF, label = label, vars = vars, k_split = k_split, 
                           distribution = distribution, 
                           eta = eta, 
                           max_depth = max_depth, 
                           n.minobsinnode = n.minobsinnode,
                           nrounds = nrounds, 
                           method = "cv", 
                           cv.folds = 5)
  
  hyper_grid = GRID[[1]]
  # print(hyper_grid)
  dev <- GRID[[2]]
  save(GRID, file = paste0("GRID", ".", output_name, ".Rdata"))
  save(hyper_grid, file = paste0("hyper_grid", ".", output_name, ".Rdata"))
  print(Sys.time())
}
```

\#\#make deviance plots

``` r
if (do_haddock_median == 1){
  source("deviance_plots_all.R")
}
```

\#\#find out what happens if we set no lower limit on number of trees

``` r
if (do_haddock_median == 1){

  source("deviance_no_min_trees.R")
  PLTS_no_min
}
```

\#\#make deviance plot just for the “best” parameters requiring at least
10000 trees as optimal number of trees

``` r
if (do_haddock_median == 1){
  source("deviance_min_trees.R")
  PLTS_min
  print("hyper_grid")
  hyper_grid
}
```

\#bootstrapGBM – run with all vertebrates and all fields

``` r
if (do_haddock_median == 1){
  distribution = "bernoulli"
  source("bootstrapGBM.R")
}
```

\#\#look at performance

``` r
if (do_haddock_median == 1){
  source("performance_metrics.R")
}
```

\#\#plot importance

``` r
if (do_haddock_median == 1){
  rm = c("eta", "max_depth", "n.trees", "eval_train", "eval_test")
  keep = setdiff(names(I),rm)
  I = I[,keep]
  
  data_long <- gather(I, key = "var", value = "value", c(2:dim(I)[2]), factor_key=TRUE)
  
  data_long_sum <- data_long %>% group_by(var) %>%
    summarize(mean_importance = mean(value))
  
  print(data_long_sum)
  data_long_sum_nonzero = subset(data_long_sum, mean_importance > 0)
  
  #find the variables that have importance at least 1
  data_long_sum_one = subset(data_long_sum, mean_importance >=1)
  var_one = data_long_sum_one$var
  var_one = c(as.character(var_one),"haddock_median_and_below")#add back label
  var_one = c(var_one, id_field)
  V = read.csv("haddock_vert_for_gbm.csv")
  col_keep = which(names(V) %in% var_one)
  V = V[,col_keep]
  write.csv(V, "haddock_vert_for_gbm_importance_over_one.csv", row.names = FALSE)
  
  data_long_nonzero = subset(data_long, var %in% data_long_sum_nonzero$var)
  
  plot <- ggplot(data = data_long_nonzero, aes(x = reorder(var, -value), y = value))+
    geom_boxplot()+
    theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, color = "black", size = 1), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2), panel.grid.major.y = element_line(color = "grey80"), panel.grid.major.x = element_line(color = "transparent"))+
    xlab("variable")+
    ylab("importance")
  
  
  plot
  ggsave(filename = paste0("importance", output_name, ".jpg"), plot = plot, height = 8, width = 8)
}
```

\#\#partial\_plotR.R – define

``` r
source("partial_plotR.R")
```

\#\#redo analysis with only features that have importance greater than
one

``` r
if (do_haddock_median == 1){
  load("output_name.Rdata")
  output_name = paste0(output_name, label, "importance_over_one")
}
```

\#bootstrapGBM – run with all vertebrates – importance over one

``` r
if (do_haddock_median == 1){
  DF = read.csv("haddock_vert_for_gbm_importance_over_one.csv")
  DF = binary_factor(DF)
  vars = setdiff(names(DF), label)
  vars = setdiff(vars, id_field)
  source("run_bootstrapGBM.R")
}
```

\#\#look at performance – importance over one

``` r
if (do_haddock_median == 1){
  source("performance_metrics.R")
}
```

\#\#plot importance (over one)

``` r
if (do_haddock_median == 1){
  rm = c("eta", "max_depth", "n.trees", "eval_train", "eval_test")
  keep = setdiff(names(I),rm)
  I = I[,keep]
  
  data_long <- gather(I, key = "var", value = "value", c(2:dim(I)[2]), factor_key=TRUE)
  
  data_long_sum <- data_long %>% group_by(var) %>%
    summarize(mean_importance = mean(value))
  
  print(data_long_sum)
  
  data_long_sum_nonzero = subset(data_long_sum, mean_importance > 0)
  
  #find the variables that have importance at least 1
  data_long_sum_one = subset(data_long_sum, mean_importance >=1)
  var_one = data_long_sum_one$var
  
  data_long_nonzero = subset(data_long, var %in% data_long_sum_nonzero$var)
  
  plot <- ggplot(data = data_long_nonzero, aes(x = reorder(var, -value), y = value))+
    geom_boxplot()+
    theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, color = "black", size = 1), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2), panel.grid.major.y = element_line(color = "grey80"), panel.grid.major.x = element_line(color = "transparent"))+
    xlab("variable")+
    ylab("importance")
  
  plot
  ggsave(filename = paste0("importance", output_name, ".jpg"), plot = plot, height = 8, width = 8)
}
```

\#\#redo analysis, this time classifying equal to or below vs. above
haddock score of species with highest haddock score that has also been
found to have conspecific transmission

``` r
if (do_haddock_infected == 1){

  load("output_name.Rdata")
  label = "haddock_infected_and_below"
  output_name = paste(output_name, label)
}
```

\#\#grid search for at/below vs above infected (hamster)

``` r
if (do_haddock_infected == 1){
  
  distribution = "bernoulli"
  DF = read.csv("haddock_vert_infected_for_gbm.csv")
  DF = binary_factor(DF)
  label_col_ind = which(names(DF)==label)
  x_col = seq(1:dim(DF)[2])
  x_col = setdiff(x_col, label_col_ind)
  vars = colnames(DF)[x_col]
  vars = setdiff(vars, id_field)#remove species
  
  GRID <- gridSearch(DF = DF, label = label, vars = vars, k_split = k_split, 
                           distribution = distribution, 
                           eta = eta, 
                           max_depth = max_depth, 
                           n.minobsinnode = n.minobsinnode,
                           nrounds = nrounds, 
                           method = "cv", 
                           cv.folds = 5)
  
  hyper_grid = GRID[[1]]
  dev <- GRID[[2]]
  save(GRID, file = paste0("GRID", ".", output_name, ".Rdata"))
  save(hyper_grid, file = paste0("hyper_grid", ".", output_name, ".Rdata"))
  print(Sys.time())
}
```

    ## [1] "2020-08-21 12:13:48 EDT"

\#\#make deviance plots

``` r
if (do_haddock_infected == 1){
  source("deviance_plots_all.R")
}
```

\#\#find out what happens if we set no lower limit on number of trees

``` r
if (do_haddock_infected == 1){
  source("deviance_no_min_trees.R")
  PLTS_no_min
}
```

    ## [[1]]

![](fishbase2_files/figure-gfm/nolower_inf-1.png)<!-- -->

\#\#make deviance plot just for the “best” parameters requiring at least
10000 trees as optimal number of trees

``` r
if (do_haddock_infected == 1){
  source("deviance_min_trees.R")
  PLTS_min
  print("hyper_grid")
  hyper_grid
}
```

    ## [1] "hyper_grid"
    ## [1] "hyper_grid"

    ##     eta max_depth n.minobsinnode n.trees eval_train eval_test
    ## 4 1e-04         3              5   64108  0.9868766 0.8642473
    ##                                       group
    ## 4 eta:1e-04, max depth:3, min obs in node:5

\#bootstrapGBM – run with all vertebrates and all fields

``` r
if (do_haddock_infected == 1){

  source("run_bootstrapGBM.R")
}
```

    ## [1] "2020-08-21 12:14:17 EDT"
    ## [1] "2020-08-21 12:40:06 EDT"

\#\#look at performance

``` r
if (do_haddock_infected == 1){

  source("performance_metrics.R")
}
```

    ## [1] "observed data, eval train"
    ## [1] 0.9839321
    ## [1] "observed data, eval test"
    ## [1] 0.8557796
    ## [1] "observed data, sd eval test"
    ## [1] 0.04186636
    ## [1] "null data, eval train"
    ## [1] "null data, eval test"
    ## [1] 0.5504032
    ## [1] "true negative"
    ## [1] 147
    ## [1] "true positive"
    ## [1] 107
    ## [1] "false negative"
    ## [1] 13
    ## [1] "false negative species"
    ##  [1] "Ceratotherium simum"     "Chaetura pelagica"      
    ##  [3] "Chlamydotis macqueenii"  "Esox lucius"            
    ##  [5] "Hippocampus comes"       "Homo sapiens"           
    ##  [7] "Leptosomus discolor"     "Merops nubicus"         
    ##  [9] "Neopelma chrysocephalum" "Nipponia nippon"        
    ## [11] "Pipistrellus abramus"    "Tauraco erythrolophus"  
    ## [13] "Xenopus tropicalis"     
    ## [1] "false_positive"
    ## [1] 11
    ## [1] "false positive species"
    ##  [1] "Amphiprion ocellaris"         "Carassius auratus"           
    ##  [3] "Electrophorus electricus"     "Ictalurus punctatus"         
    ##  [5] "Larimichthys crocea"          "Myotis brandtii"             
    ##  [7] "Paguma larvata"               "Rousettus aegyptiacus"       
    ##  [9] "Sinocyclocheilus anshuiensis" "Stegastes partitus"          
    ## [11] "Tachysurus fulvidraco"

    ## Saving 7 x 5 in image

    ## [1] "corrected test eval vert_20200821_1149 haddock_infected_and_below"
    ## [1] 0.8053763

\#\#plot importance w all vars and get vars w importance over one

``` r
if (do_haddock_infected == 1){

  file_in = "haddock_vert_infected_for_gbm.csv"
  file_one = "haddock_vert_infected_for_gbm_importance_over_one.csv"
  source("plot_importance_one.R")
}
```

    ## # A tibble: 41 x 2
    ##    var                 mean_importance
    ##    <fct>                         <dbl>
    ##  1 ClassActinopterygii         3.01   
    ##  2 ClassAmphibia               0      
    ##  3 ClassAves                   0.220  
    ##  4 ClassElasmobranchii         0      
    ##  5 ClassHolocephali            0      
    ##  6 ClassMammalia               0.0342 
    ##  7 ClassReptilia               0.00640
    ##  8 ForStrat.ground             1.11   
    ##  9 ForStrat.understory         0.535  
    ## 10 ForStrat.arboreal           4.59   
    ## # … with 31 more rows

\#\#make correlation plot for variables w/ importance \>1 for above
vs. below infected

``` r
if (do_haddock_infected == 1){

  csv_file = "haddock_vert_infected_for_gbm_importance_over_one.csv"
  
  pdf(file = paste("corr matrix", csv_file, ".pdf"))
  IF<- read.csv(csv_file)
  rm = which(names(IF) == id_field)
  if (length(rm)>0){
    IF = IF[,-rm]
  }
  
  # tmp = subset(IF, ClassActinopterygii == 1 )
  M<-cor(IF, use = "complete")
  # M<-cor(IF, use = "na.or.complete")
  plot <- corrplot(M, method="circle", type = "lower")
  
  dev.off()
}
```

    ## Warning in cor(IF, use = "complete"): the standard deviation is zero

    ## quartz_off_screen 
    ##                 2

\#\#redo analysis with only features that have importance greater than
one – infected

``` r
if (do_haddock_infected == 1){

  load("output_name.Rdata")
  output_name = paste0(output_name, label, "import_over_one_")
}
```

\#bootstrapGBM – run with all vertebrates – importance over one –
infected. use same hyper\_grid

``` r
if (do_haddock_infected == 1){

  distribution = "bernoulli"
  DF = read.csv("haddock_vert_infected_for_gbm_importance_over_one.csv")
  DF = binary_factor(DF)
  rm = "haddock_median_and_below" 
  keep = setdiff(names(DF), rm)
  DF = DF[,keep]
  vars = setdiff(names(DF), label)
  vars = setdiff(vars, id_field)
  
  source("run_bootstrapGBM.R")
}
```

    ## [1] "2020-08-21 12:40:09 EDT"
    ## [1] "2020-08-21 12:58:38 EDT"

\#\#look at performance

``` r
if (do_haddock_infected == 1){

  source("performance_metrics.R")
}
```

    ## [1] "observed data, eval train"
    ## [1] 0.9832021
    ## [1] "observed data, eval test"
    ## [1] 0.8581989
    ## [1] "observed data, sd eval test"
    ## [1] 0.04301192
    ## [1] "null data, eval train"
    ## [1] "null data, eval test"
    ## [1] 0.5508065
    ## [1] "true negative"
    ## [1] 146
    ## [1] "true positive"
    ## [1] 108
    ## [1] "false negative"
    ## [1] 12
    ## [1] "false negative species"
    ##  [1] "Ceratotherium simum"     "Chaetura pelagica"      
    ##  [3] "Chlamydotis macqueenii"  "Esox lucius"            
    ##  [5] "Hippocampus comes"       "Homo sapiens"           
    ##  [7] "Leptosomus discolor"     "Neopelma chrysocephalum"
    ##  [9] "Nipponia nippon"         "Pipistrellus abramus"   
    ## [11] "Tauraco erythrolophus"   "Xenopus tropicalis"     
    ## [1] "false_positive"
    ## [1] 12
    ## [1] "false positive species"
    ##  [1] "Amphiprion ocellaris"         "Carassius auratus"           
    ##  [3] "Electrophorus electricus"     "Ictalurus punctatus"         
    ##  [5] "Larimichthys crocea"          "Myotis brandtii"             
    ##  [7] "Paguma larvata"               "Rousettus aegyptiacus"       
    ##  [9] "Sinocyclocheilus anshuiensis" "Stegastes partitus"          
    ## [11] "Tachysurus fulvidraco"        "Zalophus californianus"

    ## Saving 7 x 5 in image

    ## [1] "corrected test eval vert_20200821_1149haddock_infected_and_belowimport_over_one_"
    ## [1] 0.8073925

\#\#plot importance - infected, vars w importance over one

``` r
if (do_haddock_infected == 1){

  source("plot_importance.R")
}
```

    ## # A tibble: 24 x 2
    ##    var                        mean_importance
    ##    <fct>                                <dbl>
    ##  1 ClassActinopterygii                   3.20
    ##  2 ForStrat.ground                       1.17
    ##  3 ForStrat.arboreal                     4.87
    ##  4 ForStrat.aerial                       1.44
    ##  5 ForStrat.marine                       1.67
    ##  6 female_maturity_d                     1.65
    ##  7 male_maturity_d                       1.98
    ##  8 log_litterclutch_size_n               4.68
    ##  9 litters_or_clutches_per_y             1.68
    ## 10 log_birthhatching_weight_g            4.50
    ## # … with 14 more rows

\#AA\_30\_negative

``` r
load("output_name.Rdata")
label = "AA_30_negative"
output_name = paste0(output_name, label)
distribution = "bernoulli"
```

\#\#now make model to predict “AA\_30\_negative”

``` r
load("DF_gbm.Rdata")
load("gridSearch.Rdata")
DF <- DF_gbm
rm= c(rm_list, "haddock_median_and_below", "haddock_score_mean" ,"haddock_infected_and_below") 
keep = setdiff(names(DF), rm)
DF = DF[,keep]
str_which(names(DF), "haddock")
```

    ## integer(0)

``` r
DF = binary_factor(DF)
write.csv(DF, file = paste0(output_name,".csv"), row.names = FALSE)
vars_gbm = setdiff(names(DF), label)
vars_gbm = setdiff(vars_gbm, id_field)

GRID <- gridSearch(DF = DF, label = label, vars = vars_gbm, k_split = k_split, 
                         distribution = distribution, 
                         eta = eta, 
                         max_depth = max_depth, 
                         n.minobsinnode = n.minobsinnode,
                         nrounds = nrounds, 
                         method = "cv", 
                         cv.folds = 5)

hyper_grid = GRID[[1]]
# print(hyper_grid)
dev <- GRID[[2]]
save(GRID, file = paste0("GRID", ".", output_name, ".Rdata"))
save(hyper_grid, file = paste0("hyper_grid", ".", output_name, ".Rdata"))
print(Sys.time())
```

    ## [1] "2020-08-21 13:30:56 EDT"

\#\#make deviance plots

``` r
source("deviance_plots_all.R")
```

\#\#find out what happens if we set no lower limit on number of trees

``` r
source("deviance_no_min_trees.R")
PLTS_no_min
```

    ## [[1]]

![](fishbase2_files/figure-gfm/dev_30_no_limit-1.png)<!-- -->

\#\#make deviance plot just for the “best” parameters requiring at least
10000 trees as optimal number of trees

``` r
source("deviance_no_min_trees.R")
PLTS_no_min
```

    ## [[1]]

![](fishbase2_files/figure-gfm/dev_best_30-1.png)<!-- -->

``` r
print("hyper_grid")
```

    ## [1] "hyper_grid"

``` r
hyper_grid
```

    ##    eta max_depth n.minobsinnode n.trees eval_train eval_test
    ## 20 0.1         2              5      43  0.9358392 0.9057143
    ##                                      group
    ## 20 eta:0.1, max depth:2, min obs in node:5

\#bootstrapGBM – run with all vertebrates and all fields –

``` r
vars = vars_gbm
source("run_bootstrapGBM.R")
```

    ## [1] "2020-08-21 13:31:56 EDT"
    ## [1] "2020-08-21 13:57:46 EDT"

\#\#look at performance

``` r
source("performance_metrics.R")
```

    ## [1] "observed data, eval train"
    ## [1] 0.9573601
    ## [1] "observed data, eval test"
    ## [1] 0.8590714
    ## [1] "observed data, sd eval test"
    ## [1] 0.04695254
    ## [1] "null data, eval train"
    ## [1] "null data, eval test"
    ## [1] 0.556
    ## [1] "true negative"
    ## [1] 86
    ## [1] "true positive"
    ## [1] 166
    ## [1] "false negative"
    ## [1] 12
    ## [1] "false negative species"
    ##  [1] "Betta splendens"        "Bos mutus"              "Calidris pugnax"       
    ##  [4] "Clupea harengus"        "Cynoglossus semilaevis" "Fundulus heteroclitus" 
    ##  [7] "Gadus morhua"           "Larimichthys crocea"    "Lepisosteus oculatus"  
    ## [10] "Oncorhynchus mykiss"    "Poecilia latipinna"     "Salarias fasciatus"    
    ## [1] "false_positive"
    ## [1] 14
    ## [1] "false positive species"
    ##  [1] "Anolis carolinensis"      "Charadrius vociferus"    
    ##  [3] "Chrysochloris asiatica"   "Danio rerio"             
    ##  [5] "Echinops telfairi"        "Eurypyga helias"         
    ##  [7] "Gallus gallus"            "Mus pahari"              
    ##  [9] "Nothobranchius furzeri"   "Ornithorhynchus anatinus"
    ## [11] "Pipistrellus abramus"     "Rhinolophus pusillus"    
    ## [13] "Vicugna pacos"            "Zonotrichia albicollis"

    ## Saving 7 x 5 in image

    ## [1] "corrected test eval vert_20200821_1149AA_30_negative"
    ## [1] 0.8030714

\#\#plot importance and output csv w importance over one

``` r
load("output_name.Rdata")
label = "AA_30_negative"
output_name = paste0(output_name, label)

file_in = paste0(output_name,".csv")
file_one = paste0(output_name, "importance_over_one", ".csv")
source("plot_importance_one.R")
```

    ## # A tibble: 40 x 2
    ##    var                 mean_importance
    ##    <fct>                         <dbl>
    ##  1 ClassActinopterygii           0.265
    ##  2 ClassAmphibia                 0    
    ##  3 ClassAves                     0.420
    ##  4 ClassElasmobranchii           0    
    ##  5 ClassHolocephali              0    
    ##  6 ClassMammalia                 0.119
    ##  7 ClassReptilia                 0.511
    ##  8 ForStrat.ground               0.999
    ##  9 ForStrat.understory           1.04 
    ## 10 ForStrat.arboreal             0.654
    ## # … with 30 more rows

\#\#redo analysis with only features that have importance greater than
one

``` r
load("output_name.Rdata")
output_name = paste0(output_name, label, "importance_over_one")
```

\#bootstrapGBM – run with all vertebrates – importance over one

``` r
DF = read.csv(paste0(output_name, ".csv"))
DF = binary_factor(DF)
vars = setdiff(names(DF), label)
vars = setdiff(vars, id_field)
print(Sys.time())
```

    ## [1] "2020-08-21 13:57:50 EDT"

``` r
source("run_bootstrapGBM.R")
```

    ## [1] "2020-08-21 13:57:50 EDT"
    ## [1] "2020-08-21 14:14:47 EDT"

\#\#look at performance

``` r
source("performance_metrics.R")
```

    ## [1] "observed data, eval train"
    ## [1] 0.9568706
    ## [1] "observed data, eval test"
    ## [1] 0.8581429
    ## [1] "observed data, sd eval test"
    ## [1] 0.05103154
    ## [1] "null data, eval train"
    ## [1] "null data, eval test"
    ## [1] 0.5644286
    ## [1] "true negative"
    ## [1] 86
    ## [1] "true positive"
    ## [1] 164
    ## [1] "false negative"
    ## [1] 14
    ## [1] "false negative species"
    ##  [1] "Betta splendens"        "Bos mutus"              "Calidris pugnax"       
    ##  [4] "Camelus ferus"          "Clupea harengus"        "Cynoglossus semilaevis"
    ##  [7] "Equus przewalskii"      "Fundulus heteroclitus"  "Gadus morhua"          
    ## [10] "Larimichthys crocea"    "Lepisosteus oculatus"   "Oncorhynchus mykiss"   
    ## [13] "Poecilia latipinna"     "Salarias fasciatus"    
    ## [1] "false_positive"
    ## [1] 14
    ## [1] "false positive species"
    ##  [1] "Anolis carolinensis"      "Charadrius vociferus"    
    ##  [3] "Chrysochloris asiatica"   "Danio rerio"             
    ##  [5] "Echinops telfairi"        "Eurypyga helias"         
    ##  [7] "Gallus gallus"            "Mus musculus"            
    ##  [9] "Mus pahari"               "Ornithorhynchus anatinus"
    ## [11] "Pipistrellus abramus"     "Rhinolophus pusillus"    
    ## [13] "Vicugna pacos"            "Zonotrichia albicollis"

    ## Saving 7 x 5 in image

    ## [1] "corrected test eval vert_20200821_1149AA_30_negativeimportance_over_one"
    ## [1] 0.7937143

\#\#plot importance

``` r
source("plot_importance.R")
```

    ## # A tibble: 19 x 2
    ##    var                        mean_importance
    ##    <fct>                                <dbl>
    ##  1 ForStrat.understory                   1.55
    ##  2 ForStrat.marine                       8.01
    ##  3 Activity.Diurnal                      2.88
    ##  4 female_maturity_d                     3.00
    ##  5 male_maturity_d                       1.54
    ##  6 log_litterclutch_size_n              13.6 
    ##  7 litters_or_clutches_per_y             2.46
    ##  8 log_birthhatching_weight_g            2.86
    ##  9 log_adult_body_mass_g                 8.77
    ## 10 longevity_y                           3.26
    ## 11 log_no_sex_body_mass                  2.65
    ## 12 adult_svl_cm                          3.30
    ## 13 diet_breadth                          5.35
    ## 14 major_habitat_type_breadth            5.86
    ## 15 tnc_ecoregion_breadth                 7.38
    ## 16 log_range_size                       12.9 
    ## 17 development_d                         2.71
    ## 18 AA_83_Y                               4.61
    ## 19 log_WOS_hits_synonyms                 7.29

\#AA\_30\_negative – including haddock score

``` r
load("output_name.Rdata")
label = "AA_30_negative"
output_name = paste0(output_name, label, "haddock_score_mean")
distribution = "bernoulli"
```

\#\#now make model to predict “AA\_30\_negative” – including haddock
score as predictor

``` r
load("DF_gbm.Rdata")
load("gridSearch.Rdata")
DF <- DF_gbm
rm= c(rm_list, "haddock_median_and_below","haddock_infected_and_below") #don't take out haddock score mean
keep = setdiff(names(DF), rm)
DF = DF[,keep]
str_which(names(DF), "haddock")#there should still be one
```

    ## [1] 8

``` r
DF = binary_factor(DF)
write.csv(DF, file = paste0(output_name,".csv"), row.names = FALSE)
vars_gbm = setdiff(names(DF), label)
vars_gbm = setdiff(vars_gbm, id_field)

GRID <- gridSearch(DF = DF, label = label, vars = vars_gbm, k_split = k_split, 
                         distribution = distribution, 
                         eta = eta, 
                         max_depth = max_depth, 
                         n.minobsinnode = n.minobsinnode,
                         nrounds = nrounds, 
                         method = "cv", 
                         cv.folds = 5)

hyper_grid = GRID[[1]]
# print(hyper_grid)
dev <- GRID[[2]]
save(GRID, file = paste0("GRID", ".", output_name, ".Rdata"))
save(hyper_grid, file = paste0("hyper_grid", ".", output_name, ".Rdata"))
print(Sys.time())
```

    ## [1] "2020-08-21 14:37:19 EDT"

\#\#make deviance plots

``` r
source("deviance_plots_all.R")
```

\#\#find out what happens if we set no lower limit on number of trees

``` r
source("deviance_no_min_trees.R")
PLTS_no_min
```

    ## [[1]]

![](fishbase2_files/figure-gfm/dev_30_no_limit_haddock-1.png)<!-- -->

\#\#make deviance plot just for the “best” parameters requiring at least
10000 trees as optimal number of trees, negative at 3

``` r
source("deviance_min_trees.R")
```

    ## [1] "hyper_grid"

``` r
PLTS_min
```

    ## [[1]]

![](fishbase2_files/figure-gfm/dev_best_30_haddock-1.png)<!-- -->

\#bootstrapGBM – run with all vertebrates and all fields – negative at
30, including haddock

``` r
vars = vars_gbm
distribution = "bernoulli"
source("run_bootstrapGBM.R")
```

    ## [1] "2020-08-21 14:38:00 EDT"
    ## [1] "2020-08-21 15:07:56 EDT"

\#\#look at performance

``` r
source("performance_metrics.R")
```

    ## [1] "observed data, eval train"
    ## [1] 0.9898427
    ## [1] "observed data, eval test"
    ## [1] 0.8814286
    ## [1] "observed data, sd eval test"
    ## [1] 0.044083
    ## [1] "null data, eval train"
    ## [1] "null data, eval test"
    ## [1] 0.5613571
    ## [1] "true negative"
    ## [1] 89
    ## [1] "true positive"
    ## [1] 169
    ## [1] "false negative"
    ## [1] 9
    ## [1] "false negative species"
    ## [1] "Betta splendens"        "Camelus ferus"          "Cynoglossus semilaevis"
    ## [4] "Fundulus heteroclitus"  "Larimichthys crocea"    "Lepisosteus oculatus"  
    ## [7] "Oncorhynchus mykiss"    "Poecilia latipinna"     "Salarias fasciatus"    
    ## [1] "false_positive"
    ## [1] 11
    ## [1] "false positive species"
    ##  [1] "Anolis carolinensis"      "Apteryx rowi"            
    ##  [3] "Danio rerio"              "Echinops telfairi"       
    ##  [5] "Eurypyga helias"          "Lipotes vexillifer"      
    ##  [7] "Mus pahari"               "Ornithorhynchus anatinus"
    ##  [9] "Pipistrellus abramus"     "Rhinolophus pusillus"    
    ## [11] "Zonotrichia albicollis"

    ## Saving 7 x 5 in image

    ## [1] "corrected test eval vert_20200821_1149AA_30_negativehaddock_score_mean"
    ## [1] 0.8200714

\#\#plot importance and output csv w importance over one w/ haddock

``` r
file_in = paste0(output_name,".csv")
file_one = paste0(output_name, "importance_over_one", ".csv")
source("plot_importance_one.R")
```

    ## # A tibble: 41 x 2
    ##    var                 mean_importance
    ##    <fct>                         <dbl>
    ##  1 ClassActinopterygii          0.157 
    ##  2 ClassAmphibia                0.0240
    ##  3 ClassAves                    0.674 
    ##  4 ClassElasmobranchii          0     
    ##  5 ClassHolocephali             0     
    ##  6 ClassMammalia                0.243 
    ##  7 ClassReptilia                0.493 
    ##  8 haddock_score_mean          11.4   
    ##  9 ForStrat.ground              0.504 
    ## 10 ForStrat.understory          0.687 
    ## # … with 31 more rows

\#\#bootstrap w/ AA 30, haddock score, import over one

``` r
load("output_name.Rdata")
output_name = paste(output_name, label, "haddock_score_mean", "importance_over_one")

DF = read.csv(file_one)
DF = binary_factor(DF)
rm = which(names(DF) %in% c("Order", "nchar", "haddock_score_sd"))#leave in haddock_score_mean
if (length(rm)>0){
  DF = DF[,-rm]
}
distribution = "bernoulli"

label_col_ind = which(names(DF)==label)
x_col = seq(1:dim(DF)[2])
x_col = setdiff(x_col, label_col_ind)
vars = colnames(DF)[x_col]
vars = setdiff(vars, id_field)#remove species
source("run_bootstrapGBM.R")
```

    ## [1] "2020-08-21 15:08:00 EDT"
    ## [1] "2020-08-21 15:26:36 EDT"

\#\#look at performance – AA 30 one haddock

``` r
source("performance_metrics.R")
```

    ## [1] "observed data, eval train"
    ## [1] 0.989493
    ## [1] "observed data, eval test"
    ## [1] 0.8832857
    ## [1] "observed data, sd eval test"
    ## [1] 0.04569965
    ## [1] "null data, eval train"
    ## [1] "null data, eval test"
    ## [1] 0.5767857
    ## [1] "true negative"
    ## [1] 90
    ## [1] "true positive"
    ## [1] 168
    ## [1] "false negative"
    ## [1] 10
    ## [1] "false negative species"
    ##  [1] "Betta splendens"        "Bos mutus"              "Camelus ferus"         
    ##  [4] "Cynoglossus semilaevis" "Fundulus heteroclitus"  "Larimichthys crocea"   
    ##  [7] "Lepisosteus oculatus"   "Oncorhynchus mykiss"    "Poecilia latipinna"    
    ## [10] "Salarias fasciatus"    
    ## [1] "false_positive"
    ## [1] 10
    ## [1] "false positive species"
    ##  [1] "Anolis carolinensis"      "Danio rerio"             
    ##  [3] "Echinops telfairi"        "Eurypyga helias"         
    ##  [5] "Lipotes vexillifer"       "Mus pahari"              
    ##  [7] "Ornithorhynchus anatinus" "Pipistrellus abramus"    
    ##  [9] "Rhinolophus pusillus"     "Zonotrichia albicollis"

    ## Saving 7 x 5 in image

    ## [1] "corrected test eval vert_20200821_1149 AA_30_negative haddock_score_mean importance_over_one"
    ## [1] 0.8065

\#\#plot importance – AA 30 haddock and get variables w importance over
one

``` r
source("plot_importance.R")
```

    ## # A tibble: 21 x 2
    ##    var                        mean_importance
    ##    <fct>                                <dbl>
    ##  1 haddock_score_mean                   12.1 
    ##  2 ForStrat.marine                       6.82
    ##  3 Activity.Nocturnal                    1.39
    ##  4 Activity.Crepuscular                  1.77
    ##  5 Activity.Diurnal                      2.17
    ##  6 female_maturity_d                     3.46
    ##  7 male_maturity_d                       1.74
    ##  8 log_litterclutch_size_n              15.2 
    ##  9 litters_or_clutches_per_y             3.14
    ## 10 log_birthhatching_weight_g            2.97
    ## # … with 11 more rows

\#\#try doing analysis with mammals

``` r
if (do_mammal == 1){
load("output_name.Rdata")
label = "haddock_infected_and_below"
taxa = "mammals"
output_name = paste(output_name, label, taxa)
load("DF_gbm.Rdata")
DF = DF_gbm

DF= subset(DF, ClassMammalia ==1)
out_names = names(DF)[str_which(names(DF), "Class")]
rm = c(rm_list, out_names, "haddock_median_and_below" )#now remove this outcome variable
keep = setdiff(names(DF), rm)
DF = DF[,keep]

write.csv(DF, file = "mammals_for_gbm.csv", row.names = FALSE)
}
```

\#\#gridsearch w/ infected – mammals

``` r
if (do_mammal == 1){
DF = read.csv("mammals_for_gbm.csv")
DF = binary_factor(DF)
rm = which(names(DF) %in% c("haddock_score_mean", "Order", "nchar", "haddock_score_sd"))
if (length(rm)>0){
  DF = DF[,-rm]
}

# n.minobsinnode = c(2)
k_split = 0.8
distribution = "bernoulli"

label_col_ind = which(names(DF)==label)
x_col = seq(1:dim(DF)[2])
x_col = setdiff(x_col, label_col_ind)
vars = colnames(DF)[x_col]
vars = setdiff(vars, id_field)#remove species

GRID <- gridSearch(DF = DF, label = label, vars = vars, k_split = k_split, 
                         distribution = distribution, 
                         eta = eta, 
                         max_depth = max_depth, 
                         n.minobsinnode = n.minobsinnode,
                         nrounds = nrounds, 
                         method = "cv", 
                         cv.folds = 5)

hyper_grid = GRID[[1]]
# print(hyper_grid)
dev <- GRID[[2]]
save(GRID, file = paste0("GRID", ".", output_name, ".Rdata"))
save(hyper_grid, file = paste0("hyper_grid", ".", output_name, ".Rdata"))
print(Sys.time())
}
```

    ## [1] "2020-08-21 15:36:45 EDT"

\#\#make deviance plots

``` r
if (do_mammal == 1){
source("deviance_plots_all.R")
}
```

\#\#find out what happens if we set no lower limit on number of trees

``` r
if (do_mammal == 1){
source("deviance_no_min_trees.R")
PLTS_no_min
}
```

    ## [[1]]

![](fishbase2_files/figure-gfm/no%20min%20mammals-1.png)<!-- -->

\#\#make deviance plot just for the “best” parameters requiring at least
10000 trees as optimal number of trees –mammals

``` r
if (do_mammal == 1){
source("deviance_min_trees.R")
PLTS_min
}
```

    ## [1] "hyper_grid"

    ## [[1]]

![](fishbase2_files/figure-gfm/dev_best_mamm-1.png)<!-- -->

\#bootstrapGBM – run with mammals

``` r
if (do_mammal == 1){
distribution = "bernoulli"
source("run_bootstrapGBM.R")
}
```

    ## [1] "2020-08-21 15:37:25 EDT"
    ## [1] "2020-08-21 15:46:18 EDT"

\#\#look at performance – mammals

``` r
if (do_mammal == 1){
source("performance_metrics.R")
}  
```

    ## [1] "observed data, eval train"
    ## [1] 0.9952031
    ## [1] "observed data, eval test"
    ## [1] 0.851049
    ## [1] "observed data, sd eval test"
    ## [1] 0.0855862
    ## [1] "null data, eval train"
    ## [1] "null data, eval test"
    ## [1] 0.593007
    ## [1] "true negative"
    ## [1] 66
    ## [1] "true positive"
    ## [1] 52
    ## [1] "false negative"
    ## [1] 6
    ## [1] "false negative species"
    ## [1] "Ceratotherium simum"   "Jaculus jaculus"       "Mesocricetus auratus" 
    ## [4] "Pipistrellus abramus"  "Rhinolophus macrotis"  "Rhinolophus pearsonii"
    ## [1] "false_positive"
    ## [1] 2
    ## [1] "false positive species"
    ## [1] "Oryctolagus cuniculus" "Rhinolophus sinicus"

    ## Saving 7 x 5 in image

    ## [1] "corrected test eval vert_20200821_1149 haddock_infected_and_below mammals"
    ## [1] 0.758042

\#\#plot importance – mammals and get variables w importance over one

``` r
if (do_mammal == 1){
file_in = "mammals_for_gbm.csv"
file_one = "mammals_for_gbm_importance_over_one.csv"
source("plot_importance_one.R")
}
```

    ## # A tibble: 34 x 2
    ##    var                  mean_importance
    ##    <fct>                          <dbl>
    ##  1 ForStrat.ground                0.745
    ##  2 ForStrat.understory            0.136
    ##  3 ForStrat.arboreal              4.02 
    ##  4 ForStrat.aerial                0.230
    ##  5 ForStrat.marine                1.14 
    ##  6 Activity.Nocturnal             0.323
    ##  7 Activity.Crepuscular           1.25 
    ##  8 Activity.Diurnal               1.34 
    ##  9 female_maturity_d              4.03 
    ## 10 male_maturity_d                7.05 
    ## # … with 24 more rows

\#\#bootstrap w/ infected – mammals – import over one

``` r
if (do_mammal == 1){
load("output_name.Rdata")
label = "haddock_infected_and_below"
taxa = "mammals"
output_name = paste(output_name, label, taxa, "importance_over_one")

DF = read.csv("mammals_for_gbm_importance_over_one.csv")
DF = binary_factor(DF)
rm = which(names(DF) %in% c("haddock_score_mean", "Order", "nchar", "haddock_score_sd"))
if (length(rm)>0){
  DF = DF[,-rm]
}
distribution = "bernoulli"

label_col_ind = which(names(DF)==label)
x_col = seq(1:dim(DF)[2])
x_col = setdiff(x_col, label_col_ind)
vars = colnames(DF)[x_col]
vars = setdiff(vars, id_field)#remove species
source("run_bootstrapGBM.R")
}
```

    ## [1] "2020-08-21 15:46:20 EDT"
    ## [1] "2020-08-21 15:54:09 EDT"

\#\#look at performance – mammals –one

``` r
if (do_mammal == 1){
source("performance_metrics.R")
}
```

    ## [1] "observed data, eval train"
    ## [1] 0.9943133
    ## [1] "observed data, eval test"
    ## [1] 0.8531469
    ## [1] "observed data, sd eval test"
    ## [1] 0.08802428
    ## [1] "null data, eval train"
    ## [1] "null data, eval test"
    ## [1] 0.5993007
    ## [1] "true negative"
    ## [1] 66
    ## [1] "true positive"
    ## [1] 52
    ## [1] "false negative"
    ## [1] 6
    ## [1] "false negative species"
    ## [1] "Ceratotherium simum"   "Jaculus jaculus"       "Mesocricetus auratus" 
    ## [4] "Pipistrellus abramus"  "Rhinolophus macrotis"  "Rhinolophus pearsonii"
    ## [1] "false_positive"
    ## [1] 2
    ## [1] "false positive species"
    ## [1] "Oryctolagus cuniculus" "Paguma larvata"

    ## Saving 7 x 5 in image

    ## [1] "corrected test eval vert_20200821_1149 haddock_infected_and_below mammals importance_over_one"
    ## [1] 0.7538462

\#\#plot importance – mammals and get variables w importance over one

``` r
if (do_mammal == 1){
source("plot_importance.R")
}
```

    ## # A tibble: 27 x 2
    ##    var                              mean_importance
    ##    <fct>                                      <dbl>
    ##  1 ForStrat.arboreal                           4.16
    ##  2 ForStrat.marine                             1.13
    ##  3 Activity.Crepuscular                        1.24
    ##  4 Activity.Diurnal                            1.41
    ##  5 female_maturity_d                           4.05
    ##  6 male_maturity_d                             7.08
    ##  7 weaning_d                                   4.26
    ##  8 log_litterclutch_size_n                     5.60
    ##  9 litters_or_clutches_per_y                   1.60
    ## 10 log_inter_litterbirth_interval_y            1.83
    ## # … with 17 more rows

\#\#now do regression analysis

``` r
if (do_score == 1){
  load("output_name.Rdata")
  label = "haddock_score_mean"
  output_name = paste(output_name, label)
  distribution = "gaussian"
}
```

``` r
if (do_score == 1){
  load("DF_gbm.Rdata")
  # load("bootstrapGBM.Rdata")
  load("gridSearch.Rdata")
  load("binary_factor.Rdata")
  DF = DF_gbm
  DF = binary_factor(DF)
  rm = which(names(DF) %in% c("haddock_median_and_below", "haddock_infected_and_below", "Order",  "nchar", "haddock_score_sd"))
  
  if (length(rm)>0){
    DF = DF[,-rm]
  }
  
  write.csv(DF, file = "haddock_score_vert_for_gbm.csv", row.names = FALSE)
  
    
  label_col_ind = which(names(DF)==label)
  x_col = seq(1:dim(DF)[2])
  x_col = setdiff(x_col, label_col_ind)
  vars = colnames(DF)[x_col]
  vars = setdiff(vars, id_field)#remove species
  distribution = "gaussian"
  
  GRID <- gridSearch(DF = DF, label = label, vars = vars, k_split = k_split, 
                           distribution = distribution, 
                           eta = eta, 
                           max_depth = max_depth, 
                           n.minobsinnode = n.minobsinnode,
                           nrounds = nrounds, 
                           method = "cv", 
                           cv.folds = 5)
  
  hyper_grid = GRID[[1]]
  dev <- GRID[[2]]
  save(GRID, file = paste0("GRID", ".", output_name, ".Rdata"))
  save(hyper_grid, file = paste0("hyper_grid", ".", output_name, ".Rdata"))
  print(Sys.time())
  hyper_grid
}
```

    ## [1] "2020-08-21 16:12:51 EDT"

    ##      eta max_depth n.minobsinnode n.trees eval_train eval_test
    ## 1  1e-04         2              2   76889  0.6589906 0.3949209
    ## 2  1e-04         2              5   73534  0.6523962 0.3932300
    ## 3  1e-04         3              2   48305  0.6485872 0.4067840
    ## 4  1e-04         3              5   59021  0.6796365 0.4067741
    ## 5  1e-04         4              2   64024  0.7273870 0.4109896
    ## 6  1e-04         4              5   38152  0.6433715 0.4185736
    ## 7  1e-03         2              2    6134  0.6262726 0.3935995
    ## 8  1e-03         2              5    8549  0.6734341 0.3903337
    ## 9  1e-03         3              2    7348  0.7136463 0.4033448
    ## 10 1e-03         3              5    3835  0.6125445 0.4061408
    ## 11 1e-03         4              2    3529  0.6319883 0.4174513
    ## 12 1e-03         4              5    3853  0.6461973 0.4156531
    ## 13 1e-02         2              2     470  0.5858175 0.4061915
    ## 14 1e-02         2              5    1108  0.7058881 0.3993749
    ## 15 1e-02         3              2     557  0.6657174 0.3972502
    ## 16 1e-02         3              5     377  0.6114259 0.4140017
    ## 17 1e-02         4              2     645  0.7269317 0.3990993
    ## 18 1e-02         4              5     509  0.6906080 0.4216515
    ## 19 1e-01         2              2      89  0.6633685 0.3912643
    ## 20 1e-01         2              5      47  0.5950124 0.3795128
    ## 21 1e-01         3              2      55  0.6622631 0.3898047
    ## 22 1e-01         3              5      56  0.6516021 0.3926344
    ## 23 1e-01         4              2      83  0.7615607 0.4249757
    ## 24 1e-01         4              5      54  0.6933983 0.4248148
    ##                                        group
    ## 1  eta:1e-04, max depth:2, min obs in node:2
    ## 2  eta:1e-04, max depth:2, min obs in node:5
    ## 3  eta:1e-04, max depth:3, min obs in node:2
    ## 4  eta:1e-04, max depth:3, min obs in node:5
    ## 5  eta:1e-04, max depth:4, min obs in node:2
    ## 6  eta:1e-04, max depth:4, min obs in node:5
    ## 7  eta:0.001, max depth:2, min obs in node:2
    ## 8  eta:0.001, max depth:2, min obs in node:5
    ## 9  eta:0.001, max depth:3, min obs in node:2
    ## 10 eta:0.001, max depth:3, min obs in node:5
    ## 11 eta:0.001, max depth:4, min obs in node:2
    ## 12 eta:0.001, max depth:4, min obs in node:5
    ## 13  eta:0.01, max depth:2, min obs in node:2
    ## 14  eta:0.01, max depth:2, min obs in node:5
    ## 15  eta:0.01, max depth:3, min obs in node:2
    ## 16  eta:0.01, max depth:3, min obs in node:5
    ## 17  eta:0.01, max depth:4, min obs in node:2
    ## 18  eta:0.01, max depth:4, min obs in node:5
    ## 19   eta:0.1, max depth:2, min obs in node:2
    ## 20   eta:0.1, max depth:2, min obs in node:5
    ## 21   eta:0.1, max depth:3, min obs in node:2
    ## 22   eta:0.1, max depth:3, min obs in node:5
    ## 23   eta:0.1, max depth:4, min obs in node:2
    ## 24   eta:0.1, max depth:4, min obs in node:5

\#\#make deviance plots

``` r
if (do_score == 1){
  source("deviance_plots_all.R")
}
```

\#\#find out what happens if we set no lower limit on number of trees

``` r
if (do_score == 1){
  source("deviance_no_min_trees.R"  ) 
  PLTS_no_min
}
```

    ## [[1]]

![](fishbase2_files/figure-gfm/dev_no_min_reg-1.png)<!-- -->

\#\#make deviance plot just for the “best” parameters requiring at least
10000 trees as optimal number of trees

``` r
if (do_score == 1){
  source("deviance_min_trees.R")
  PLTS_min
print("hyper_grid")
hyper_grid
}
```

    ## [1] "hyper_grid"
    ## [1] "hyper_grid"

    ##     eta max_depth n.minobsinnode n.trees eval_train eval_test
    ## 6 1e-04         4              5   38152  0.6433715 0.4185736
    ##                                       group
    ## 6 eta:1e-04, max depth:4, min obs in node:5

\#bootstrapGBM – run with all vertebrates and all fields – regression

``` r
if (do_score == 1){
  source("run_bootstrapGBM.R")
}
```

    ## [1] "2020-08-21 16:13:19 EDT"
    ## [1] "2020-08-21 16:36:16 EDT"

\#\#look at performance

``` r
if (do_score == 1){
  source("performance_metrics.R")
}
```

    ## [1] "observed data, eval train"
    ## [1] 0.7187324
    ## [1] "observed data, eval test"
    ## [1] 0.3939597
    ## [1] "observed data, sd eval test"
    ## [1] 0.06171556
    ## [1] "null data, eval train"
    ## [1] "null data, eval test"
    ## [1] -0.006438871
    ## [1] "corrected test eval vert_20200821_1149 haddock_score_mean"
    ## [1] 0.4003986

\#\#plot importance and get vars w import over one

``` r
if (do_score == 1){
  file_in = "haddock_score_vert_for_gbm.csv"
  file_one = "haddock_score_vert_for_gbm_over_one.csv"
  source("plot_importance_one.R")
}
```

    ## # A tibble: 41 x 2
    ##    var                 mean_importance
    ##    <fct>                         <dbl>
    ##  1 ClassActinopterygii         6.06   
    ##  2 ClassAmphibia               0      
    ##  3 ClassAves                   0.160  
    ##  4 ClassElasmobranchii         0      
    ##  5 ClassHolocephali            0      
    ##  6 ClassMammalia               0.141  
    ##  7 ClassReptilia               0.00796
    ##  8 ForStrat.ground             1.73   
    ##  9 ForStrat.understory         0.316  
    ## 10 ForStrat.arboreal           2.22   
    ## # … with 31 more rows

``` r
if (do_score == 1){
  load("output_name.Rdata")
  output_name = paste(output_name, label, "_importance_over_one")
}
```

\#bootstrapGBM – run with all vertebrates – importance over one –
regression use same hyper\_grid

``` r
if (do_score == 1){
  DF = read.csv("haddock_score_vert_for_gbm_over_one.csv")
  DF = binary_factor(DF)
  rm = c("haddock_median_and_below", "haddock_infected_and_below")
  keep = setdiff(names(DF), rm)
  DF = DF[,keep]
  vars = setdiff(names(DF), label)
  vars = setdiff(vars, id_field)
  source("run_bootstrapGBM.R")
}
```

    ## [1] "2020-08-21 16:36:17 EDT"
    ## [1] "2020-08-21 16:52:57 EDT"

\#\#look at performance

``` r
if (do_score == 1){
  source("performance_metrics.R")
}
```

    ## [1] "observed data, eval train"
    ## [1] 0.7238577
    ## [1] "observed data, eval test"
    ## [1] 0.3973572
    ## [1] "observed data, sd eval test"
    ## [1] 0.0593225
    ## [1] "null data, eval train"
    ## [1] "null data, eval test"
    ## [1] -0.006868869
    ## [1] "corrected test eval vert_20200821_1149 haddock_score_mean _importance_over_one"
    ## [1] 0.4042261

\#\#plot importance – regression, import over one

``` r
if (do_score == 1){
  source("plot_importance.R")
}
```

    ## # A tibble: 26 x 2
    ##    var                        mean_importance
    ##    <fct>                                <dbl>
    ##  1 ClassActinopterygii                   6.15
    ##  2 ForStrat.ground                       1.78
    ##  3 ForStrat.arboreal                     2.42
    ##  4 ForStrat.aerial                       1.05
    ##  5 ForStrat.marine                       1.34
    ##  6 Activity.Crepuscular                  1.93
    ##  7 female_maturity_d                     1.42
    ##  8 male_maturity_d                       1.37
    ##  9 log_litterclutch_size_n               6.70
    ## 10 log_birthhatching_weight_g            2.37
    ## # … with 16 more rows

\#\#write outputs – performance\_out and DF\_merge\_out

``` r
load("output_name.Rdata")
write.csv(performance_out, paste0("performance_models", output_name, ".csv"))
write.csv(DF_merge_out, paste0("predictions_models", output_name, ".csv"))
```

\#\#find out how long it takes to run

``` r
print(Sys.time())
```

    ## [1] "2020-08-21 16:52:59 EDT"
