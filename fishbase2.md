fishbase2
================
Han lab
8/18/2020

\#\#\#\#\#install packages

``` r
testing_var = 1#0: testing out, not doing full grid search
nruns = 5
date = "20200818"
time = "1833"
output_name = paste("vert", date, time, sep = "_")
save(output_name, file = "output_name.Rdata")
print(Sys.time())
```

    ## [1] "2020-08-18 18:33:52 EDT"

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
do_haddock_median = 0

do_haddock_infected = 0#make 1 to do analysis with over/under infectable species w/ highest haddock score

do_score = 0#1 to do regression analysis with haddock score

get_fishbase = 0

make_v = 0 #whether to remake V, which takes a little time to get all the data from fishbase

do_mammal = 0#whether to do separate analysis for mammals (excluding other verts)
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

\#add fields and remove “major\_habitat\_type\_breadth” (because it
seems redundant to tnc ecoregion breadth)

\#\#add AA value to rest of vert data

\#\#add WOS hits from R package wosr

``` r
if (make_v == 1){
  load("V.Rdata")
  dim(V)
  
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
  save(V, file = "V.Rdata")
}
```

\#\#output: “haddock\_vert\_for\_gbm.csv”

``` r
print(Sys.time())
```

    ## [1] "2020-08-18 18:34:05 EDT"

``` r
load("gridSearch.Rdata")
# save(output_name, file = "output_name.Rdata")
load("V.Rdata")

# tmp = subset(V, Class == "Actinopterygii")


out = V
Species = out$Species
sp_ind = which(names(out)=="Species")
dmy <- dummyVars(" ~ .", data = out[,-sp_ind])
out <- data.frame(predict(dmy, newdata = out))
out$Species = Species
V = out

names = names(V)
DF = V
Species = DF$Species
##remove variables with near zero variation
sp_ind = which(names(DF)=="Species")
nzv = nearZeroVar(DF, freqCut = 95/5, saveMetrics = TRUE)
okay_inds = which(nzv$nzv == FALSE)

near_zero_inds =which(nzv$nzv == TRUE)
print("near zero variation fields")
```

    ## [1] "near zero variation fields"

``` r
names[near_zero_inds]
```

    ## character(0)

``` r
length(okay_inds)
```

    ## [1] 42

``` r
DF = DF[,okay_inds]#include only the columns that have variation
DF$Species = Species

load("binary_factor.Rdata")
DF = binary_factor(DF)
V = DF

V$adult_svl_cm[is.nan(V$adult_svl_cm)] <- NA
V$log_adult_body_mass_g[is.nan(V$log_adult_body_mass_g)] <- NA

DF = V
#find out what haddock_score_median is across all species
haddock_median = median(V$haddock_score_mean)

haddock_median_and_below = rep(0, dim(DF)[1])
inds = which(DF$haddock_score_mean <= haddock_median)
haddock_median_and_below[inds]= 1
DF$haddock_median_and_below = haddock_median_and_below
DF_gbm = DF
save(DF_gbm, file = "DF_gbm.Rdata")

rm =  c("haddock_score_mean",   "haddock_score_sd")
keep = setdiff(names(DF),rm)
DF = DF[,keep]
write.csv(DF, file = "haddock_vert_for_gbm.csv", row.names = FALSE)
```

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

\#bootstrapGBM – run with all vertebrates and all fields

``` r
if (do_haddock_infected == 1){

  source("run_bootstrapGBM.R")
}
```

\#\#look at performance

``` r
if (do_haddock_infected == 1){

  source("performance_metrics.R")
}
```

\#\#plot importance w all vars and get vars w importance over one

``` r
if (do_haddock_infected == 1){

  file_in = "haddock_vert_infected_for_gbm.csv"
  file_one = "haddock_vert_infected_for_gbm_importance_over_one.csv"
  source("plot_importance_one.R")
}
```

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

\#\#look at performance

``` r
if (do_haddock_infected == 1){

  source("performance_metrics.R")
}
```

\#\#plot importance - infected, vars w importance over one

``` r
if (do_haddock_infected == 1){

  source("plot_importance.R")
}
```

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

\#bootstrapGBM – run with all vertebrates and all fields – regression

``` r
if (do_score == 1){
  source("run_bootstrapGBM.R")
}
```

\#\#look at performance

``` r
if (do_score == 1){
  source("performance_metrics.R")
}
```

\#\#plot importance and get vars w import over one

``` r
if (do_score == 1){
  file_in = "haddock_score_vert_for_gbm.csv"
  file_one = "haddock_score_vert_for_gbm_over_one.csv"
  source("plot_importance_one.R")
}
```

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

\#\#look at performance

``` r
if (do_score == 1){
  source("performance_metrics.R")
}
```

\#\#plot importance – regression, import over one

``` r
if (do_score == 1){
  source("plot_importance.R")
}
```

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

    ## [1] "2020-08-18 18:53:20 EDT"

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
    ## 24 0.1         4              5      36  0.9643926 0.8971429
    ##                                      group
    ## 24 eta:0.1, max depth:4, min obs in node:5

\#bootstrapGBM – run with all vertebrates and all fields –

``` r
vars = vars_gbm
source("run_bootstrapGBM.R")
```

    ## [1] "2020-08-18 18:53:42 EDT"
    ## [1] "2020-08-18 19:10:05 EDT"

\#\#look at performance

``` r
source("performance_metrics.R")
```

    ## [1] "observed data, eval train"
    ## [1] 0.9801585
    ## [1] "observed data, eval test"
    ## [1] 0.8591429
    ## [1] "observed data, sd eval test"
    ## [1] 0.05862889
    ## [1] "null data, eval train"
    ## [1] "null data, eval test"
    ## [1] 0.5675714
    ## [1] "corrected test eval vert_20200818_1833AA_30_negative"
    ## [1] 0.7915714
    ## [1] "true negative"
    ## [1] 85
    ## [1] "true positive"
    ## [1] 168
    ## [1] "false negative"
    ## [1] 9
    ## [1] "false negative species"
    ## [1] "Betta splendens"        "Clupea harengus"        "Cynoglossus semilaevis"
    ## [4] "Fundulus heteroclitus"  "Gadus morhua"           "Larimichthys crocea"   
    ## [7] "Lepisosteus oculatus"   "Oncorhynchus mykiss"    "Salarias fasciatus"    
    ## [1] "false_positive"
    ## [1] 15
    ## [1] "false positive species"
    ##  [1] "Anolis carolinensis"         "Chrysochloris asiatica"     
    ##  [3] "Danio rerio"                 "Delphinapterus leucas"      
    ##  [5] "Echinops telfairi"           "Eurypyga helias"            
    ##  [7] "Gallus gallus"               "Mus pahari"                 
    ##  [9] "Neophocaena asiaeorientalis" "Numida meleagris"           
    ## [11] "Ornithorhynchus anatinus"    "Pipistrellus abramus"       
    ## [13] "Rhinolophus pusillus"        "Vicugna pacos"              
    ## [15] "Zonotrichia albicollis"

    ## Saving 7 x 5 in image

\#\#plot importance and output csv w importance over one

``` r
load("output_name.Rdata")
label = "AA_30_negative"
output_name = paste0(output_name, label)

file_in = paste0(output_name,".csv")
file_one = paste0(output_name, "importance_over_one", ".csv")
source("plot_importance_one.R")
```

    ## # A tibble: 36 x 2
    ##    var                 mean_importance
    ##    <fct>                         <dbl>
    ##  1 ClassActinopterygii          0.426 
    ##  2 ClassAves                    0.903 
    ##  3 ClassMammalia                0.202 
    ##  4 ClassReptilia                0.869 
    ##  5 ForStrat.ground              1.28  
    ##  6 ForStrat.understory          1.64  
    ##  7 ForStrat.arboreal            1.01  
    ##  8 ForStrat.aerial              0.0844
    ##  9 ForStrat.marine              7.68  
    ## 10 Activity.Nocturnal           0.387 
    ## # … with 26 more rows

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

    ## [1] "2020-08-18 19:10:08 EDT"

``` r
source("run_bootstrapGBM.R")
```

    ## [1] "2020-08-18 19:10:08 EDT"
    ## [1] "2020-08-18 19:18:34 EDT"

\#\#look at performance

``` r
source("performance_metrics.R")
```

    ## [1] "observed data, eval train"
    ## [1] 0.9696215
    ## [1] "observed data, eval test"
    ## [1] 0.8677143
    ## [1] "observed data, sd eval test"
    ## [1] 0.05944128
    ## [1] "null data, eval train"
    ## [1] "null data, eval test"
    ## [1] 0.5932857
    ## [1] "corrected test eval vert_20200818_1833AA_30_negativeimportance_over_one"
    ## [1] 0.7744286
    ## [1] "true negative"
    ## [1] 86
    ## [1] "true positive"
    ## [1] 168
    ## [1] "false negative"
    ## [1] 9
    ## [1] "false negative species"
    ## [1] "Betta splendens"        "Calidris pugnax"        "Cynoglossus semilaevis"
    ## [4] "Fundulus heteroclitus"  "Gadus morhua"           "Larimichthys crocea"   
    ## [7] "Lepisosteus oculatus"   "Oncorhynchus mykiss"    "Salarias fasciatus"    
    ## [1] "false_positive"
    ## [1] 14
    ## [1] "false positive species"
    ##  [1] "Anolis carolinensis"      "Charadrius vociferus"    
    ##  [3] "Chrysochloris asiatica"   "Danio rerio"             
    ##  [5] "Delphinapterus leucas"    "Echinops telfairi"       
    ##  [7] "Eurypyga helias"          "Gallus gallus"           
    ##  [9] "Mus pahari"               "Ornithorhynchus anatinus"
    ## [11] "Pipistrellus abramus"     "Rhinolophus pusillus"    
    ## [13] "Vicugna pacos"            "Zonotrichia albicollis"

    ## Saving 7 x 5 in image

\#\#plot importance

``` r
source("plot_importance.R")
```

    ## # A tibble: 20 x 2
    ##    var                        mean_importance
    ##    <fct>                                <dbl>
    ##  1 ForStrat.ground                       1.41
    ##  2 ForStrat.understory                   1.90
    ##  3 ForStrat.arboreal                     1.18
    ##  4 ForStrat.marine                       7.64
    ##  5 Activity.Crepuscular                  1.78
    ##  6 Activity.Diurnal                      2.60
    ##  7 female_maturity_d                     3.46
    ##  8 male_maturity_d                       1.14
    ##  9 development_d                         4.63
    ## 10 log_litterclutch_size_n              12.8 
    ## 11 log_birthhatching_weight_g            1.78
    ## 12 log_adult_body_mass_g                 6.28
    ## 13 longevity_y                           5.84
    ## 14 adult_svl_cm                          2.51
    ## 15 diet_breadth                          4.80
    ## 16 tnc_ecoregion_breadth                 8.99
    ## 17 mass_specific_production              2.49
    ## 18 log_range_size                       15.2 
    ## 19 AA_83_Y                               5.65
    ## 20 log_WOS_hits_synonyms                 7.91

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

    ## [1] 5

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

    ## [1] "2020-08-18 19:35:10 EDT"

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

    ## [1] "2020-08-18 19:35:46 EDT"
    ## [1] "2020-08-18 19:47:00 EDT"

\#\#look at performance

``` r
source("performance_metrics.R")
```

    ## [1] "observed data, eval train"
    ## [1] 0.9855458
    ## [1] "observed data, eval test"
    ## [1] 0.8925714
    ## [1] "observed data, sd eval test"
    ## [1] 0.05020773
    ## [1] "null data, eval train"
    ## [1] "null data, eval test"
    ## [1] 0.6112857
    ## [1] "corrected test eval vert_20200818_1833AA_30_negativehaddock_score_mean"
    ## [1] 0.7812857
    ## [1] "true negative"
    ## [1] 89
    ## [1] "true positive"
    ## [1] 165
    ## [1] "false negative"
    ## [1] 12
    ## [1] "false negative species"
    ##  [1] "Betta splendens"        "Calidris pugnax"        "Camelus ferus"         
    ##  [4] "Clupea harengus"        "Cynoglossus semilaevis" "Fundulus heteroclitus" 
    ##  [7] "Gadus morhua"           "Larimichthys crocea"    "Lepisosteus oculatus"  
    ## [10] "Oncorhynchus mykiss"    "Poecilia latipinna"     "Salarias fasciatus"    
    ## [1] "false_positive"
    ## [1] 11
    ## [1] "false positive species"
    ##  [1] "Anolis carolinensis"      "Apteryx rowi"            
    ##  [3] "Chrysochloris asiatica"   "Echinops telfairi"       
    ##  [5] "Eurypyga helias"          "Mus pahari"              
    ##  [7] "Ornithorhynchus anatinus" "Pipistrellus abramus"    
    ##  [9] "Rhinolophus pusillus"     "Vicugna pacos"           
    ## [11] "Zonotrichia albicollis"

    ## Saving 7 x 5 in image

\#\#plot importance and output csv w importance over one w/ haddock

``` r
file_in = paste0(output_name,".csv")
file_one = paste0(output_name, "importance_over_one", ".csv")
source("plot_importance_one.R")
```

    ## # A tibble: 37 x 2
    ##    var                 mean_importance
    ##    <fct>                         <dbl>
    ##  1 ClassActinopterygii           0.249
    ##  2 ClassAves                     0.971
    ##  3 ClassMammalia                 0.151
    ##  4 ClassReptilia                 0.743
    ##  5 haddock_score_mean           12.5  
    ##  6 ForStrat.ground               1.05 
    ##  7 ForStrat.understory           1.03 
    ##  8 ForStrat.arboreal             0.540
    ##  9 ForStrat.aerial               0.132
    ## 10 ForStrat.marine               7.83 
    ## # … with 27 more rows

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

    ## [1] "2020-08-18 19:47:02 EDT"
    ## [1] "2020-08-18 19:54:52 EDT"

\#\#look at performance – AA 30 one haddock

``` r
source("performance_metrics.R")
```

    ## [1] "observed data, eval train"
    ## [1] 0.9834331
    ## [1] "observed data, eval test"
    ## [1] 0.8962857
    ## [1] "observed data, sd eval test"
    ## [1] 0.0510362
    ## [1] "null data, eval train"
    ## [1] "null data, eval test"
    ## [1] 0.5727143
    ## [1] "corrected test eval vert_20200818_1833 AA_30_negative haddock_score_mean importance_over_one"
    ## [1] 0.8235714
    ## [1] "true negative"
    ## [1] 88
    ## [1] "true positive"
    ## [1] 164
    ## [1] "false negative"
    ## [1] 13
    ## [1] "false negative species"
    ##  [1] "Betta splendens"        "Calidris pugnax"        "Camelus ferus"         
    ##  [4] "Clupea harengus"        "Cynoglossus semilaevis" "Equus przewalskii"     
    ##  [7] "Fundulus heteroclitus"  "Gadus morhua"           "Larimichthys crocea"   
    ## [10] "Lepisosteus oculatus"   "Oncorhynchus mykiss"    "Poecilia latipinna"    
    ## [13] "Salarias fasciatus"    
    ## [1] "false_positive"
    ## [1] 12
    ## [1] "false positive species"
    ##  [1] "Anolis carolinensis"      "Apteryx rowi"            
    ##  [3] "Chrysochloris asiatica"   "Delphinapterus leucas"   
    ##  [5] "Echinops telfairi"        "Eurypyga helias"         
    ##  [7] "Mus pahari"               "Ornithorhynchus anatinus"
    ##  [9] "Pipistrellus abramus"     "Rhinolophus pusillus"    
    ## [11] "Vicugna pacos"            "Zonotrichia albicollis"

    ## Saving 7 x 5 in image

\#\#plot importance – AA 30 haddock and get variables w importance over
one

``` r
source("plot_importance.R")
```

    ## # A tibble: 21 x 2
    ##    var                     mean_importance
    ##    <fct>                             <dbl>
    ##  1 haddock_score_mean                12.6 
    ##  2 ForStrat.ground                    1.32
    ##  3 ForStrat.understory                1.32
    ##  4 ForStrat.marine                    8.06
    ##  5 Activity.Nocturnal                 1.42
    ##  6 Activity.Crepuscular               1.85
    ##  7 Activity.Diurnal                   2.35
    ##  8 female_maturity_d                  4.27
    ##  9 development_d                      3.23
    ## 10 log_litterclutch_size_n           14.0 
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

\#\#make deviance plot just for the “best” parameters requiring at least
10000 trees as optimal number of trees –mammals

``` r
if (do_mammal == 1){
source("deviance_min_trees.R")
PLTS_min
}
```

\#bootstrapGBM – run with mammals

``` r
if (do_mammal == 1){
distribution = "bernoulli"
source("run_bootstrapGBM.R")
}
```

\#\#look at performance – mammals

``` r
if (do_mammal == 1){
source("performance_metrics.R")
}  
```

\#\#plot importance – mammals and get variables w importance over one

``` r
if (do_mammal == 1){
file_in = "mammals_for_gbm.csv"
file_one = "mammals_for_gbm_importance_over_one.csv"
source("plot_importance_one.R")
}
```

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

\#\#look at performance – mammals –one

``` r
if (do_mammal == 1){
source("performance_metrics.R")
}
```

\#\#plot importance – mammals and get variables w importance over one

``` r
if (do_mammal == 1){
source("plot_importance.R")
}
```

\#\#find out how long it takes to run

``` r
print(Sys.time())
```

    ## [1] "2020-08-18 19:54:53 EDT"
