fishbase2
================
Han lab
8/18/2020

\#\#\#\#\#install packages

``` r
testing_var = 1#0: testing out, not doing full grid search
nruns = 10
date = "20200818"
time = "2229"
output_name = paste("vert", date, time, sep = "_")
save(output_name, file = "output_name.Rdata")
print(Sys.time())
```

    ## [1] "2020-08-18 22:29:52 EDT"

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

do_haddock_infected = 1#make 1 to do analysis with over/under infectable species w/ highest haddock score

do_score = 1#1 to do regression analysis with haddock score

get_fishbase = 0

make_v = 0 #whether to remake V, which takes a little time to get all the data from fishbase

do_mammal = 1#whether to do separate analysis for mammals (excluding other verts)
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

    ## [1] "2020-08-18 22:30:07 EDT"

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

    ## [1] "2020-08-18 22:48:16 EDT"

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
    ## 5 1e-04         4              2   64289  0.9874508 0.9228612
    ##                                       group
    ## 5 eta:1e-04, max depth:4, min obs in node:2

\#bootstrapGBM – run with all vertebrates and all fields

``` r
if (do_haddock_infected == 1){

  source("run_bootstrapGBM.R")
}
```

    ## [1] "2020-08-18 22:48:46 EDT"
    ## [1] "2020-08-18 23:10:52 EDT"

\#\#look at performance

``` r
if (do_haddock_infected == 1){

  source("performance_metrics.R")
}
```

    ## [1] "observed data, eval train"
    ## [1] 0.9937008
    ## [1] "observed data, eval test"
    ## [1] 0.8760168
    ## [1] "observed data, sd eval test"
    ## [1] 0.03999423
    ## [1] "null data, eval train"
    ## [1] "null data, eval test"
    ## [1] 0.5638149
    ## [1] "corrected test eval vert_20200818_2229 haddock_infected_and_below"
    ## [1] 0.812202
    ## [1] "true negative"
    ## [1] 151
    ## [1] "true positive"
    ## [1] 109
    ## [1] "false negative"
    ## [1] 10
    ## [1] "false negative species"
    ##  [1] "Ceratotherium simum"     "Chaetura pelagica"      
    ##  [3] "Chlamydotis macqueenii"  "Esox lucius"            
    ##  [5] "Hippocampus comes"       "Merops nubicus"         
    ##  [7] "Neopelma chrysocephalum" "Nipponia nippon"        
    ##  [9] "Pipistrellus abramus"    "Xenopus tropicalis"     
    ## [1] "false_positive"
    ## [1] 7
    ## [1] "false positive species"
    ## [1] "Amphiprion ocellaris"         "Carassius auratus"           
    ## [3] "Electrophorus electricus"     "Ictalurus punctatus"         
    ## [5] "Larimichthys crocea"          "Paguma larvata"              
    ## [7] "Sinocyclocheilus anshuiensis"

    ## Saving 7 x 5 in image

\#\#plot importance w all vars and get vars w importance over one

``` r
if (do_haddock_infected == 1){

  file_in = "haddock_vert_infected_for_gbm.csv"
  file_one = "haddock_vert_infected_for_gbm_importance_over_one.csv"
  source("plot_importance_one.R")
}
```

    ## # A tibble: 37 x 2
    ##    var                 mean_importance
    ##    <fct>                         <dbl>
    ##  1 ClassActinopterygii          2.65  
    ##  2 ClassAves                    0.279 
    ##  3 ClassMammalia                0.0401
    ##  4 ClassReptilia                0.0249
    ##  5 ForStrat.ground              0.966 
    ##  6 ForStrat.understory          0.729 
    ##  7 ForStrat.arboreal            3.22  
    ##  8 ForStrat.aerial              2.75  
    ##  9 ForStrat.marine              1.22  
    ## 10 Activity.Nocturnal           0.277 
    ## # … with 27 more rows

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

    ## [1] "2020-08-18 23:10:54 EDT"
    ## [1] "2020-08-18 23:25:41 EDT"

\#\#look at performance

``` r
if (do_haddock_infected == 1){

  source("performance_metrics.R")
}
```

    ## [1] "observed data, eval train"
    ## [1] 0.9930118
    ## [1] "observed data, eval test"
    ## [1] 0.8778401
    ## [1] "observed data, sd eval test"
    ## [1] 0.03995732
    ## [1] "null data, eval train"
    ## [1] "null data, eval test"
    ## [1] 0.5698457
    ## [1] "corrected test eval vert_20200818_2229haddock_infected_and_belowimport_over_one_"
    ## [1] 0.8079944
    ## [1] "true negative"
    ## [1] 151
    ## [1] "true positive"
    ## [1] 109
    ## [1] "false negative"
    ## [1] 10
    ## [1] "false negative species"
    ##  [1] "Ceratotherium simum"     "Chaetura pelagica"      
    ##  [3] "Chlamydotis macqueenii"  "Esox lucius"            
    ##  [5] "Hippocampus comes"       "Mustela putorius"       
    ##  [7] "Neopelma chrysocephalum" "Nipponia nippon"        
    ##  [9] "Pipistrellus abramus"    "Xenopus tropicalis"     
    ## [1] "false_positive"
    ## [1] 7
    ## [1] "false positive species"
    ## [1] "Amphiprion ocellaris"         "Carassius auratus"           
    ## [3] "Electrophorus electricus"     "Ictalurus punctatus"         
    ## [5] "Larimichthys crocea"          "Paguma larvata"              
    ## [7] "Sinocyclocheilus anshuiensis"

    ## Saving 7 x 5 in image

\#\#plot importance - infected, vars w importance over one

``` r
if (do_haddock_infected == 1){

  source("plot_importance.R")
}
```

    ## # A tibble: 25 x 2
    ##    var                       mean_importance
    ##    <fct>                               <dbl>
    ##  1 ClassActinopterygii                  2.83
    ##  2 ForStrat.arboreal                    3.89
    ##  3 ForStrat.aerial                      2.94
    ##  4 ForStrat.marine                      1.47
    ##  5 female_maturity_d                    1.79
    ##  6 male_maturity_d                      2.45
    ##  7 weaning_d                            1.14
    ##  8 development_d                        4.48
    ##  9 log_litterclutch_size_n              6.54
    ## 10 litters_or_clutches_per_y            1.51
    ## # … with 15 more rows

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

    ## [1] "2020-08-18 23:40:06 EDT"

    ##      eta max_depth n.minobsinnode n.trees eval_train eval_test
    ## 1  1e-04         2              2   69408  0.6210104 0.5266804
    ## 2  1e-04         2              5   52833  0.5807066 0.5153667
    ## 3  1e-04         3              2   48905  0.6291407 0.5213793
    ## 4  1e-04         3              5   60483  0.6622681 0.5291226
    ## 5  1e-04         4              2   38574  0.6255358 0.5133957
    ## 6  1e-04         4              5   50448  0.6697525 0.5257476
    ## 7  1e-03         2              2    6750  0.6167301 0.5267036
    ## 8  1e-03         2              5    6426  0.6085476 0.5212646
    ## 9  1e-03         3              2    4683  0.6219045 0.5228762
    ## 10 1e-03         3              5    5090  0.6353209 0.5235166
    ## 11 1e-03         4              2    5675  0.6900258 0.5287751
    ## 12 1e-03         4              5    5756  0.6909918 0.5311535
    ## 13 1e-02         2              2    1057  0.6806448 0.5433700
    ## 14 1e-02         2              5     455  0.5572571 0.5015936
    ## 15 1e-02         3              2     773  0.6978104 0.5329998
    ## 16 1e-02         3              5     491  0.6268556 0.5221500
    ## 17 1e-02         4              2     476  0.6565509 0.5184193
    ## 18 1e-02         4              5     355  0.6051325 0.5084454
    ## 19 1e-01         2              2      58  0.5771288 0.5088985
    ## 20 1e-01         2              5      99  0.6718390 0.5467305
    ## 21 1e-01         3              2      52  0.6133918 0.5430688
    ## 22 1e-01         3              5      41  0.5787429 0.4819854
    ## 23 1e-01         4              2     110  0.7859531 0.5272194
    ## 24 1e-01         4              5      43  0.6310450 0.4979511
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
    ## 4 1e-04         3              5   60483  0.6622681 0.5291226
    ##                                       group
    ## 4 eta:1e-04, max depth:3, min obs in node:5

\#bootstrapGBM – run with all vertebrates and all fields – regression

``` r
if (do_score == 1){
  source("run_bootstrapGBM.R")
}
```

    ## [1] "2020-08-18 23:40:27 EDT"
    ## [1] "2020-08-18 23:55:11 EDT"

\#\#look at performance

``` r
if (do_score == 1){
  source("performance_metrics.R")
}
```

    ## [1] "observed data, eval train"
    ## [1] 0.6923532
    ## [1] "observed data, eval test"
    ## [1] 0.4198591
    ## [1] "observed data, sd eval test"
    ## [1] 0.08275417
    ## [1] "null data, eval train"
    ## [1] "null data, eval test"
    ## [1] -0.007356416
    ## [1] "corrected test eval vert_20200818_2229 haddock_score_mean"
    ## [1] 0.9272155
    ## [1] "true negative"
    ## [1] 0
    ## [1] "true positive"
    ## [1] 0
    ## [1] "false negative"
    ## [1] 0
    ## [1] "false negative species"
    ## character(0)
    ## [1] "false_positive"
    ## [1] 0
    ## [1] "false positive species"
    ## character(0)

    ## Saving 7 x 5 in image

\#\#plot importance and get vars w import over one

``` r
if (do_score == 1){
  file_in = "haddock_score_vert_for_gbm.csv"
  file_one = "haddock_score_vert_for_gbm_over_one.csv"
  source("plot_importance_one.R")
}
```

    ## # A tibble: 37 x 2
    ##    var                 mean_importance
    ##    <fct>                         <dbl>
    ##  1 ClassActinopterygii         6.99   
    ##  2 ClassAves                   0.389  
    ##  3 ClassMammalia               0.172  
    ##  4 ClassReptilia               0.00550
    ##  5 ForStrat.ground             1.85   
    ##  6 ForStrat.understory         0.426  
    ##  7 ForStrat.arboreal           1.97   
    ##  8 ForStrat.aerial             1.08   
    ##  9 ForStrat.marine             1.56   
    ## 10 Activity.Nocturnal          0.244  
    ## # … with 27 more rows

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

    ## [1] "2020-08-18 23:55:13 EDT"
    ## [1] "2020-08-19 00:07:19 EDT"

\#\#look at performance

``` r
if (do_score == 1){
  source("performance_metrics.R")
}
```

    ## [1] "observed data, eval train"
    ## [1] 0.6994588
    ## [1] "observed data, eval test"
    ## [1] 0.4242069
    ## [1] "observed data, sd eval test"
    ## [1] 0.07988844
    ## [1] "null data, eval train"
    ## [1] "null data, eval test"
    ## [1] -0.006576666
    ## [1] "corrected test eval vert_20200818_2229 haddock_score_mean _importance_over_one"
    ## [1] 0.9307836
    ## [1] "true negative"
    ## [1] 0
    ## [1] "true positive"
    ## [1] 0
    ## [1] "false negative"
    ## [1] 0
    ## [1] "false negative species"
    ## character(0)
    ## [1] "false_positive"
    ## [1] 0
    ## [1] "false positive species"
    ## character(0)

    ## Saving 7 x 5 in image

\#\#plot importance – regression, import over one

``` r
if (do_score == 1){
  source("plot_importance.R")
}
```

    ## # A tibble: 27 x 2
    ##    var                     mean_importance
    ##    <fct>                             <dbl>
    ##  1 ClassActinopterygii                7.09
    ##  2 ForStrat.ground                    1.95
    ##  3 ForStrat.arboreal                  2.17
    ##  4 ForStrat.aerial                    1.12
    ##  5 ForStrat.marine                    1.60
    ##  6 Activity.Crepuscular               2.36
    ##  7 female_maturity_d                  1.39
    ##  8 male_maturity_d                    1.33
    ##  9 development_d                      4.14
    ## 10 log_litterclutch_size_n            7.36
    ## # … with 17 more rows

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

    ## [1] "2020-08-19 00:21:31 EDT"

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
    ## 20 0.1         2              5      46  0.9554577      0.91
    ##                                      group
    ## 20 eta:0.1, max depth:2, min obs in node:5

\#bootstrapGBM – run with all vertebrates and all fields –

``` r
vars = vars_gbm
source("run_bootstrapGBM.R")
```

    ## [1] "2020-08-19 00:22:04 EDT"
    ## [1] "2020-08-19 00:33:15 EDT"

\#\#look at performance

``` r
source("performance_metrics.R")
```

    ## [1] "observed data, eval train"
    ## [1] 0.950647
    ## [1] "observed data, eval test"
    ## [1] 0.8742143
    ## [1] "observed data, sd eval test"
    ## [1] 0.04405013
    ## [1] "null data, eval train"
    ## [1] "null data, eval test"
    ## [1] 0.5725714
    ## [1] "corrected test eval vert_20200818_2229AA_30_negative"
    ## [1] 0.8016429
    ## [1] "true negative"
    ## [1] 84
    ## [1] "true positive"
    ## [1] 164
    ## [1] "false negative"
    ## [1] 13
    ## [1] "false negative species"
    ##  [1] "Betta splendens"        "Bos mutus"              "Calidris pugnax"       
    ##  [4] "Clupea harengus"        "Cynoglossus semilaevis" "Equus przewalskii"     
    ##  [7] "Fundulus heteroclitus"  "Gadus morhua"           "Larimichthys crocea"   
    ## [10] "Lepisosteus oculatus"   "Oncorhynchus mykiss"    "Poecilia latipinna"    
    ## [13] "Salarias fasciatus"    
    ## [1] "false_positive"
    ## [1] 16
    ## [1] "false positive species"
    ##  [1] "Anolis carolinensis"      "Apteryx rowi"            
    ##  [3] "Charadrius vociferus"     "Chrysochloris asiatica"  
    ##  [5] "Danio rerio"              "Delphinapterus leucas"   
    ##  [7] "Echinops telfairi"        "Eurypyga helias"         
    ##  [9] "Gallus gallus"            "Mus musculus"            
    ## [11] "Mus pahari"               "Ornithorhynchus anatinus"
    ## [13] "Pipistrellus abramus"     "Rhinolophus pusillus"    
    ## [15] "Vicugna pacos"            "Zonotrichia albicollis"

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
    ##  1 ClassActinopterygii          0.319 
    ##  2 ClassAves                    0.515 
    ##  3 ClassMammalia                0.157 
    ##  4 ClassReptilia                0.470 
    ##  5 ForStrat.ground              1.04  
    ##  6 ForStrat.understory          1.16  
    ##  7 ForStrat.arboreal            0.792 
    ##  8 ForStrat.aerial              0.0939
    ##  9 ForStrat.marine              6.95  
    ## 10 Activity.Nocturnal           0.661 
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

    ## [1] "2020-08-19 00:33:17 EDT"

``` r
source("run_bootstrapGBM.R")
```

    ## [1] "2020-08-19 00:33:17 EDT"
    ## [1] "2020-08-19 00:40:56 EDT"

\#\#look at performance

``` r
source("performance_metrics.R")
```

    ## [1] "observed data, eval train"
    ## [1] 0.9497271
    ## [1] "observed data, eval test"
    ## [1] 0.8612857
    ## [1] "observed data, sd eval test"
    ## [1] 0.06559409
    ## [1] "null data, eval train"
    ## [1] "null data, eval test"
    ## [1] 0.5735714
    ## [1] "corrected test eval vert_20200818_2229AA_30_negativeimportance_over_one"
    ## [1] 0.7877143
    ## [1] "true negative"
    ## [1] 85
    ## [1] "true positive"
    ## [1] 163
    ## [1] "false negative"
    ## [1] 14
    ## [1] "false negative species"
    ##  [1] "Betta splendens"        "Bos mutus"              "Calidris pugnax"       
    ##  [4] "Clupea harengus"        "Cynoglossus semilaevis" "Equus przewalskii"     
    ##  [7] "Fundulus heteroclitus"  "Gadus morhua"           "Larimichthys crocea"   
    ## [10] "Lepisosteus oculatus"   "Oncorhynchus mykiss"    "Poecilia latipinna"    
    ## [13] "Salarias fasciatus"     "Xiphophorus hellerii"  
    ## [1] "false_positive"
    ## [1] 15
    ## [1] "false positive species"
    ##  [1] "Anolis carolinensis"      "Charadrius vociferus"    
    ##  [3] "Chrysochloris asiatica"   "Danio rerio"             
    ##  [5] "Delphinapterus leucas"    "Echinops telfairi"       
    ##  [7] "Eurypyga helias"          "Gallus gallus"           
    ##  [9] "Mus musculus"             "Mus pahari"              
    ## [11] "Ornithorhynchus anatinus" "Pipistrellus abramus"    
    ## [13] "Rhinolophus pusillus"     "Vicugna pacos"           
    ## [15] "Zonotrichia albicollis"

    ## Saving 7 x 5 in image

\#\#plot importance

``` r
source("plot_importance.R")
```

    ## # A tibble: 19 x 2
    ##    var                        mean_importance
    ##    <fct>                                <dbl>
    ##  1 ForStrat.ground                       1.30
    ##  2 ForStrat.understory                   1.52
    ##  3 ForStrat.marine                       6.99
    ##  4 Activity.Crepuscular                  1.76
    ##  5 Activity.Diurnal                      1.64
    ##  6 female_maturity_d                     3.35
    ##  7 development_d                         6.01
    ##  8 log_litterclutch_size_n              12.8 
    ##  9 litters_or_clutches_per_y             1.66
    ## 10 log_birthhatching_weight_g            2.00
    ## 11 log_adult_body_mass_g                 9.63
    ## 12 longevity_y                           6.81
    ## 13 adult_svl_cm                          4.02
    ## 14 diet_breadth                          4.33
    ## 15 tnc_ecoregion_breadth                 8.72
    ## 16 mass_specific_production              1.50
    ## 17 log_range_size                       14.7 
    ## 18 AA_83_Y                               4.17
    ## 19 log_WOS_hits_synonyms                 7.10

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

    ## [1] "2020-08-19 00:55:42 EDT"

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

    ## [1] "2020-08-19 00:56:12 EDT"
    ## [1] "2020-08-19 01:07:25 EDT"

\#\#look at performance

``` r
source("performance_metrics.R")
```

    ## [1] "observed data, eval train"
    ## [1] 0.9717254
    ## [1] "observed data, eval test"
    ## [1] 0.8851429
    ## [1] "observed data, sd eval test"
    ## [1] 0.03830419
    ## [1] "null data, eval train"
    ## [1] "null data, eval test"
    ## [1] 0.5882143
    ## [1] "corrected test eval vert_20200818_2229AA_30_negativehaddock_score_mean"
    ## [1] 0.7969286
    ## [1] "true negative"
    ## [1] 86
    ## [1] "true positive"
    ## [1] 164
    ## [1] "false negative"
    ## [1] 13
    ## [1] "false negative species"
    ##  [1] "Betta splendens"        "Bos mutus"              "Calidris pugnax"       
    ##  [4] "Camelus ferus"          "Clupea harengus"        "Cynoglossus semilaevis"
    ##  [7] "Fundulus heteroclitus"  "Gadus morhua"           "Larimichthys crocea"   
    ## [10] "Lepisosteus oculatus"   "Oncorhynchus mykiss"    "Poecilia latipinna"    
    ## [13] "Salarias fasciatus"    
    ## [1] "false_positive"
    ## [1] 14
    ## [1] "false positive species"
    ##  [1] "Anolis carolinensis"      "Apteryx rowi"            
    ##  [3] "Chrysochloris asiatica"   "Danio rerio"             
    ##  [5] "Echinops telfairi"        "Eurypyga helias"         
    ##  [7] "Gallus gallus"            "Lipotes vexillifer"      
    ##  [9] "Mus pahari"               "Ornithorhynchus anatinus"
    ## [11] "Pipistrellus abramus"     "Rattus norvegicus"       
    ## [13] "Rhinolophus pusillus"     "Zonotrichia albicollis"

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
    ##  1 ClassActinopterygii           0.146
    ##  2 ClassAves                     0.623
    ##  3 ClassMammalia                 0.137
    ##  4 ClassReptilia                 0.395
    ##  5 haddock_score_mean           10.8  
    ##  6 ForStrat.ground               0.304
    ##  7 ForStrat.understory           0.876
    ##  8 ForStrat.arboreal             0.537
    ##  9 ForStrat.aerial               0.102
    ## 10 ForStrat.marine               8.02 
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

    ## [1] "2020-08-19 01:07:27 EDT"
    ## [1] "2020-08-19 01:15:25 EDT"

\#\#look at performance – AA 30 one haddock

``` r
source("performance_metrics.R")
```

    ## [1] "observed data, eval train"
    ## [1] 0.9707658
    ## [1] "observed data, eval test"
    ## [1] 0.884
    ## [1] "observed data, sd eval test"
    ## [1] 0.04006118
    ## [1] "null data, eval train"
    ## [1] "null data, eval test"
    ## [1] 0.5745714
    ## [1] "corrected test eval vert_20200818_2229 AA_30_negative haddock_score_mean importance_over_one"
    ## [1] 0.8094286
    ## [1] "true negative"
    ## [1] 86
    ## [1] "true positive"
    ## [1] 162
    ## [1] "false negative"
    ## [1] 15
    ## [1] "false negative species"
    ##  [1] "Antrostomus carolinensis" "Betta splendens"         
    ##  [3] "Bos mutus"                "Calidris pugnax"         
    ##  [5] "Camelus ferus"            "Clupea harengus"         
    ##  [7] "Cuculus canorus"          "Cynoglossus semilaevis"  
    ##  [9] "Fundulus heteroclitus"    "Gadus morhua"            
    ## [11] "Larimichthys crocea"      "Lepisosteus oculatus"    
    ## [13] "Oncorhynchus mykiss"      "Poecilia latipinna"      
    ## [15] "Salarias fasciatus"      
    ## [1] "false_positive"
    ## [1] 14
    ## [1] "false positive species"
    ##  [1] "Anolis carolinensis"      "Apteryx rowi"            
    ##  [3] "Chrysochloris asiatica"   "Danio rerio"             
    ##  [5] "Echinops telfairi"        "Etheostoma spectabile"   
    ##  [7] "Eurypyga helias"          "Gallus gallus"           
    ##  [9] "Lipotes vexillifer"       "Mus pahari"              
    ## [11] "Ornithorhynchus anatinus" "Pipistrellus abramus"    
    ## [13] "Rhinolophus pusillus"     "Zonotrichia albicollis"

    ## Saving 7 x 5 in image

\#\#plot importance – AA 30 haddock and get variables w importance over
one

``` r
source("plot_importance.R")
```

    ## # A tibble: 21 x 2
    ##    var                       mean_importance
    ##    <fct>                               <dbl>
    ##  1 haddock_score_mean                  11.2 
    ##  2 ForStrat.marine                      8.48
    ##  3 Activity.Nocturnal                   1.62
    ##  4 Activity.Crepuscular                 2.23
    ##  5 Activity.Diurnal                     2.47
    ##  6 female_maturity_d                    4.00
    ##  7 male_maturity_d                      1.42
    ##  8 development_d                        4.59
    ##  9 log_litterclutch_size_n             16.6 
    ## 10 litters_or_clutches_per_y            2.33
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

    ## [1] "2020-08-19 01:21:55 EDT"

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

    ## [1] "2020-08-19 01:22:34 EDT"
    ## [1] "2020-08-19 01:33:06 EDT"

\#\#look at performance – mammals

``` r
if (do_mammal == 1){
source("performance_metrics.R")
}  
```

    ## [1] "observed data, eval train"
    ## [1] 0.9997628
    ## [1] "observed data, eval test"
    ## [1] 0.8636364
    ## [1] "observed data, sd eval test"
    ## [1] 0.06994949
    ## [1] "null data, eval train"
    ## [1] "null data, eval test"
    ## [1] 0.6283217
    ## [1] "corrected test eval vert_20200818_2229 haddock_infected_and_below mammals"
    ## [1] 0.7353147
    ## [1] "true negative"
    ## [1] 68
    ## [1] "true positive"
    ## [1] 56
    ## [1] "false negative"
    ## [1] 1
    ## [1] "false negative species"
    ## [1] "Ceratotherium simum"
    ## [1] "false_positive"
    ## [1] 0
    ## [1] "false positive species"
    ## character(0)

    ## Saving 7 x 5 in image

\#\#plot importance – mammals and get variables w importance over one

``` r
if (do_mammal == 1){
file_in = "mammals_for_gbm.csv"
file_one = "mammals_for_gbm_importance_over_one.csv"
source("plot_importance_one.R")
}
```

    ## # A tibble: 33 x 2
    ##    var                  mean_importance
    ##    <fct>                          <dbl>
    ##  1 ForStrat.ground                0.650
    ##  2 ForStrat.understory            0.182
    ##  3 ForStrat.arboreal              3.55 
    ##  4 ForStrat.aerial                0.164
    ##  5 ForStrat.marine                0.695
    ##  6 Activity.Nocturnal             0.446
    ##  7 Activity.Crepuscular           0.577
    ##  8 Activity.Diurnal               0.602
    ##  9 female_maturity_d              3.71 
    ## 10 male_maturity_d                6.59 
    ## # … with 23 more rows

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

    ## [1] "2020-08-19 01:33:08 EDT"
    ## [1] "2020-08-19 01:41:52 EDT"

\#\#look at performance – mammals –one

``` r
if (do_mammal == 1){
source("performance_metrics.R")
}
```

    ## [1] "observed data, eval train"
    ## [1] 0.9993281
    ## [1] "observed data, eval test"
    ## [1] 0.8636364
    ## [1] "observed data, sd eval test"
    ## [1] 0.06956002
    ## [1] "null data, eval train"
    ## [1] "null data, eval test"
    ## [1] 0.6311189
    ## [1] "corrected test eval vert_20200818_2229 haddock_infected_and_below mammals importance_over_one"
    ## [1] 0.7325175
    ## [1] "true negative"
    ## [1] 68
    ## [1] "true positive"
    ## [1] 55
    ## [1] "false negative"
    ## [1] 2
    ## [1] "false negative species"
    ## [1] "Ceratotherium simum"  "Rhinolophus macrotis"
    ## [1] "false_positive"
    ## [1] 0
    ## [1] "false positive species"
    ## character(0)

    ## Saving 7 x 5 in image

\#\#plot importance – mammals and get variables w importance over one

``` r
if (do_mammal == 1){
source("plot_importance.R")
}
```

    ## # A tibble: 24 x 2
    ##    var                              mean_importance
    ##    <fct>                                      <dbl>
    ##  1 ForStrat.arboreal                           3.66
    ##  2 female_maturity_d                           3.81
    ##  3 male_maturity_d                             6.74
    ##  4 weaning_d                                   2.91
    ##  5 development_d                               7.17
    ##  6 log_litterclutch_size_n                     3.92
    ##  7 litters_or_clutches_per_y                   2.32
    ##  8 log_inter_litterbirth_interval_y            4.36
    ##  9 log_birthhatching_weight_g                  4.88
    ## 10 log_weaning_weight_g                        4.66
    ## # … with 14 more rows

\#\#find out how long it takes to run

``` r
print(Sys.time())
```

    ## [1] "2020-08-19 01:41:54 EDT"
