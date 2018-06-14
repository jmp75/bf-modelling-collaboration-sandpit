library(readr)
library(stringr)
library(uchronia)
library(swift)
library(dplyr)

source('./common_functions.r')

########
# Start of the main batch commands
########

root_path <- 'C:/Users/per202/Documents/BF/AUS_Catchments/AUS_Catchments'

# We can and should define a time span for the simulation, and a length of warmup for the calibration:
simul_start <- ISOdate(1990,1,1)
simul_end <- simul_start + lubridate::years(10) - lubridate::days(1)
warmup_days <- 365
model_property <- paste0( model_property_prefix, 'runoff')

p_space <- get_gr4j_parameter_space_df()
# We remove x2 from the feasible parameter space
p_space <- p_space[c(1,3,4),]
p_space <- get_gr4j_parameter_space(p_space)


# # Test:
# cat_id <- "404207"
# blah <- calibrate_catchment(cat_id, root_path, simul_start, simul_end, warmup_days, model_property, obj_id='NSE', model_config=p_space)

## Batch loop.
short_fnames <- list.files(root_path, pattern='*.csv')
cat_ids <- stringr::str_replace_all(short_fnames, '\\.csv', '')
calibration_possible_catchments <- sapply(cat_ids, FUN=catchment_has_runoff)
cat_ids <- cat_ids[calibration_possible_catchments]
# Catchment 403224 has only zeroes for valid observations - this is a problem for NSE (and not meaninful data to calibrate against)
cat_ids <- setdiff(cat_ids, '403224')

univ_func <- function (cat_id) {
  print(paste0("Calibrating ", cat_id))
  calibrate_catchment(cat_id, root_path, simul_start, simul_end, warmup_days, model_property, obj_id='NSE', model_config=p_space)
}

batch_results <- lapply(cat_ids, FUN=univ_func)
batch_results <- dplyr::bind_rows(batch_results)
batch_results$Catchment_ID <- cat_ids
  
out_file = 'c:/tmp/calib_tests.csv'
readr::write_csv(batch_results, out_file)

