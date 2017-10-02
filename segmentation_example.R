######################## SESYNC Research Support: Food and Landscape Diversity ##############################
## Performing Object Based Image Analyses (OBIA) using segments from TerrSet.
## 
## DATE CREATED: 10/02/2017
## DATE MODIFIED: 10/02/2017
## AUTHORS: Benoit Parmentier 
## PROJECT: Food and Landscape Diversity
## ISSUE: 
## TO DO:
##
## COMMIT: initial commit using segments in R
##
## Links to investigate:
##https://stats.idre.ucla.edu/r/dae/logit-regression/
#
#
###################################################
#

###### Library used

library(gtools)                              # loading some useful tools 
library(sp)                                  # Spatial pacakge with class definition by Bivand et al.
library(spdep)                               # Spatial pacakge with methods and spatial stat. by Bivand et al.
library(rgdal)                               # GDAL wrapper for R, spatial utilities
library(gdata)                               # various tools with xls reading, cbindX
library(rasterVis)                           # Raster plotting functions
library(parallel)                            # Parallelization of processes with multiple cores
library(maptools)                            # Tools and functions for sp and other spatial objects e.g. spCbind
library(maps)                                # Tools and data for spatial/geographic objects
library(plyr)                                # Various tools including rbind.fill
library(spgwr)                               # GWR method
library(rgeos)                               # Geometric, topologic library of functions
library(gridExtra)                           # Combining lattice plots
library(colorRamps)                          # Palette/color ramps for symbology
library(ggplot2)                             # Plot package 
library(lubridate)                           # Date utility fuctions
library(dplyr)                               # data manipulation and wrangling
library(ROCR)                                # ROC curve package
library(pROC)                                # prob ROC curve
library(TOC)                                 # TOC and ROC curve package
library(randomForest)                        # random forests
library(lattice)                             # Plot package
library(caret)                               # Modeling with assessment hold outs, CV folds and data splitting
library(gplots)                              # Plot package
library(sf)

#
###### Functions used in this script and sourced from other files

create_dir_fun <- function(outDir,out_suffix=NULL){
  #if out_suffix is not null then append out_suffix string
  if(!is.null(out_suffix)){
    out_name <- paste("output_",out_suffix,sep="")
    outDir <- file.path(outDir,out_name)
  }
  #create if does not exists
  if(!file.exists(outDir)){
    dir.create(outDir)
  }
  return(outDir)
}

#Used to load RData object saved within the functions produced.
load_obj <- function(f){
  env <- new.env()
  nm <- load(f, env)[1]
  env[[nm]]
}

### Other functions ####

#function_modeling <- "roc_weed_risk_functions_08112017d.R" #PARAM 1 #changed this to another file
script_path <- "/nfs/bparmentier-data/Data/projects/FoodandLandscapeDiversity/scripts" #path to script #PARAM 

#source(file.path(script_path,function_modeling)) #source all functions used in this script 1.

############################################################################
####################  Parameters and argument set up ###########

#in_dir <- "/nfs/bparmentier-data/Data/projects/modeling_weed_risk/data"
#/nfs/foodandlandscapediversity-data
out_dir <- "/nfs/bparmentier-data/Data/projects/FoodandLandscapeDiversity/segmentation/"
in_dir <- "/nfs/bparmentier-data/Data/projects/FoodandLandscapeDiversity/segmentation/output_data"

num_cores <- 2 #param 8 #normally I use only 1 core for this dataset but since I want to use the mclappy function the number of cores is changed to 2. If it was 1 then mclappy will be reverted back to the lapply function
create_out_dir_param=TRUE # param 9

out_suffix <-"segmentation_example_10022017" #output suffix for the files and ouptut folder #param 12

infile_data <- "test_0.shp" #segments level zero

model_names <- c("kmeans")

##############################  START SCRIPT  ############################

######### PART 0: Set up the output dir ################

if(is.null(out_dir)){
  out_dir <- in_dir #output will be created in the input dir
  
}
#out_dir <- in_dir #output will be created in the input dir

out_suffix_s <- out_suffix #can modify name of output suffix
if(create_out_dir_param==TRUE){
  out_dir <- create_dir_fun(out_dir,out_suffix)
  setwd(out_dir)
}else{
  setwd(out_dir) #use previoulsy defined directory
}

options(scipen=999)  #remove scientific writing


### PART I READ AND PREPARE DATA #######
#set up the working directory
#Create output directory

segments_sf <- st_read(file.path(in_dir,infile_data))
plot(segments_sf$geometry)

dim(segments_sf)

################################ END OF SCRIPT ###################