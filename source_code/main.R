# ================================================
# Initialize
# ================================================
if( dev.cur() > 1 ) dev.off()
rm( list = ls( envir = globalenv() ), envir = globalenv() )
rm(list = ls(all.names = TRUE))

basedir=dirname(here::here())

# ================================================
# Load
# ================================================
#load parameters
load(paste0(basedir,"/materials/map.RData"))#to keep 
#load functions
source(paste0(basedir,"/source_code/common_functions.R"))
source(paste0(basedir,"/source_code/bootstrap_test.R"))
source(paste0(basedir,"/source_code/T_test_ESS.R"))

# ================================================
# Run
# ================================================
#transform to flock-frame coordinate system
source(paste0(basedir,"/source_code/GPS_to_ECEF.R"))
source(paste0(basedir,"/source_code/ECEF_to_Flock_frame.R"))
source(paste0(basedir,"/source_code/visualization.R"))
##aggregation of the non drone-experienced individuals
source(paste0(basedir,"/source_code/aggregation_analysis.R"))

#statistical analysis
T_test_ESS_analysis()

 #additional bootstrap analysis
full_analysis(birdid="#1", type = "y")
full_analysis(birdid="#1", type = "z")
full_analysis(birdid="average_dist_of_#2-4", type = "ave.dist_to_center_of_NT")
