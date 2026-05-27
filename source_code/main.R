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

# ================================================
# Run
# ================================================
#transform to flock-frame coordinate system
source(paste0(basedir,"/source_code/GPS_to_ECEF.R"))
source(paste0(basedir,"/source_code/ECEF_to_Flock_frame.R"))
source(paste0(basedir,"/source_code/visualization.R"))

#statistical analysis
##position of the drone-experienced individual
full_analysis(birdid="#1", type = "y")
full_analysis(birdid="#1", type = "z")

##aggregation of individuals other than the drone-experienced one
source(paste0(basedir,"/source_code/aggregation_analysis.R"))
full_analysis(birdid="#5(centerpoint)", type = "ave.dist_center_of_three")
