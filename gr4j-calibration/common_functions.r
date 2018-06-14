library(readr)
library(stringr)
library(uchronia)
# library(swift)
library(dplyr)


## Define functions to help concise batch optimizations
## Most of this file is function definitions

load_cat_series <- function(cat_id, root_path = 'C:/Users/per202/Documents/BF/AUS_Catchments/AUS_Catchments') {
  csv_file = file.path(root_path,paste0(cat_id, '.csv'))
  stopifnot(file.exists(csv_file))
  d <- readr::read_csv(csv_file)
  d
}

create_lumped_rrmodel <- function(series, simul_start = NA, simul_end = NA, 
    model_id = "GR4J", tstep = "daily", varname_rain = "P", varname_pet = "E") {
  # sSpan <- paste0(simul_start, "/", simul_end)
  ms <- createSubarea(model_id, 1)
  if(is.na(simul_start)) simul_start <- ISOdate(series[1,'year'], series[1,'month'], series[1,'date'])
  s <- as.POSIXct(simul_start, tz='GMT')
  is_negative <- function(val) { val < 0 } 
  rainfall <- uchronia::mkDailySeries(s, series$Rain, isMissingFunc = is_negative)
  pet <- uchronia::mkDailySeries(s, series$Etp, isMissingFunc = is_negative)
  stopifnot(all(!is.na(rainfall)))
  stopifnot(all(!is.na(pet)))
  if(is.na(simul_end)) simul_end <- end(rainfall)
  e <- as.POSIXct(simul_end, tz='GMT')
  setSimulationSpan(ms, s, e)
  setSimulationTimeStep(ms, tstep)
  subAreaName <- getSubareaIds(ms)[1]
  playSubareaInput(ms, input = rainfall, subAreaName, varname_rain)
  playSubareaInput(ms, input = pet, subAreaName, varname_pet)
  return(ms)
}

get_observed_runoff <- function(series, simul_start = NA) {
  if(is.na(simul_start)) simul_start <- ISOdate(series[1,'year'], series[1,'month'], series[1,'date'])
  s <- as.POSIXct(simul_start, tz='GMT')
  is_negative <- function(val) { val < 0 } 
  uchronia::mkDailySeries(s, series$Qobs, isMissingFunc = is_negative)
}

model_property_prefix <- 'subarea.Subarea.'

get_gr4j_parameter_space_df <- function() {
  pSpecGr4j <- joki::getFreeParams('GR4J')
  pSpecGr4j$Value <- c(542.1981111, 0, 7.7403390, 1.2388548)
  pSpecGr4j$Min <- c(1,-30, 1,1)
  pSpecGr4j$Max <- c(3000, 30, 1000, 10)
  pSpecGr4j$Name <- paste0(model_property_prefix, pSpecGr4j$Name)
  return(pSpecGr4j)
}

get_gr4j_parameter_space <- function(p_space_df=get_gr4j_parameter_space_df()) {
  p <- createParameterizer(type='Generic', p_space_df)
  return(p)
}

years_series <- function(tts, subset_start, n_years) {
  stopifnot(subset_start >= start(tts))
  stopifnot(subset_start <= end(tts))
  window(tts, start=subset_start, end=subset_start + lubridate::years(n_years)) 
}

plot_obs_vs_calc <- function(calib_inputs, sim_start = NA, warmup_days = 365, n_years = 3, ylab="runoff (mm)") {
  obs <- calib_inputs$observation
  simulation <- CloneModel_R(calib_inputs$simulation)
  model_config <- calib_inputs$model_config
  applySysConfig(model_config, simulation)
  recordState(simulation, calib_inputs$model_property_id)

  w <- sim_start + lubridate::days(warmup_days)
  e <- w + lubridate::years(n_years)

  setSimulationSpan(simulation, sim_start, e)
  execSimulation(simulation)
  calc <- getRecorded(simulation, calib_inputs$model_property_id)
  joki::plotTwoSeries(obs, calc, ylab=ylab, startTime = w, endTime = e)
}

create_calib_data <- function(simulation, model_property_id, simul_start, warmup_days = 365, simul_end, obs_runoff, obj_id = 'NSE', model_config) {
  simulation <- CloneModel_R(simulation)
  w <- simul_start + lubridate::days(warmup_days)
  setSimulationSpan(simulation, simul_start, simul_end)
  recordState(simulation, model_property_id)
  calib_inputs = list(
    simulation = simulation,
    simul_start = simul_start,
    warmup_to = w,
    simul_end = simul_end,
    obj_id = obj_id,
    observation = obs_runoff,
    model_config = model_config,
    model_property_id = model_property_id
  )
  return(calib_inputs)
}

create_objective <- function(calib_inputs) {
  simulation <- CloneModel_R(calib_inputs$simulation)
  x <- calib_inputs
  setSimulationSpan(simulation, x$simul_start, x$simul_end)
  objective <- createObjective(simulation, x$model_property_id, x$observation, x$obj_id, x$warmup_to, x$simul_end)
  return(objective)
}

calibrate <- function(calib_data, maxIterations = 10000) {
  sceParams <- SCEParameters(nparams = 4)
  urs <- createParameterSampler(0, calib_data$model_config, 'urs')
  term <- getMaxIterationTermination(maxIterations = maxIterations ) 
  obj <- create_objective(calib_data)
  optimizer <- createSceOptimSwift(obj, term, SCEpars=sceParams, urs)
  calibLogger <- setCalibrationLogger(optimizer, '')
  startTime <- lubridate::now()
  calibResults <- executeOptimization(optimizer)
  endTime <- lubridate::now()
  calibWallTime <- endTime-startTime
  print(paste( 'Optimization completed in ', calibWallTime, attr(calibWallTime, 'units')))
  return(list(Results=calibResults, Optimizer=optimizer))
}


calibrate_catchment <- function(cat_id, root_path, simul_start, simul_end, warmup_days, model_property, obj_id, model_config) {
  # Define the base simulation model that will be a building block in the calibration
  d <- load_cat_series(cat_id, root_path=root_path)
  obs_runoff <- get_observed_runoff(d)
  if(all(is.na(obs_runoff))) stop(paste0("Not possible to calibrate without runoff observations for catchment ", cat_id))

  simulation <- create_lumped_rrmodel(d)
  # Apply a default parameter set to the base simulation - notably to set x2 to zero in case it is not calibrated
  applySysConfig(get_gr4j_parameter_space(), simulation)
  # Some of the catchments seem to have no observed runoff data - catch them:
  calib_data <- create_calib_data(simulation, model_property, simul_start, warmup_days, simul_end, obs_runoff, obj_id = obj_id, model_config = model_config)

  ## Viewing interactively:
  # sim_start <- ISOdate(2000,1,1)
  # plot_obs_vs_calc(calib_data, sim_start)

  obj <- create_objective (calib_data)
  # getScore(obj, get_gr4j_parameter_space())
  optim_results <- calibrate(calib_data)
  calibResults <- optim_results$Results
  sortedResults <- sortByScore(calibResults, obj_id)
  results_df <- scoresAsDataFrame(sortedResults)
  # head(results_df)
  return(results_df[1,])
}

# optional, details viz look at calibration logs. 
# optimizer <- optim_results$Optimizer
# d <- getLoggerContent(optimizer)
# d$PointNumber = 1:nrow(d)
# logMh <- mhplot::mkOptimLog(d, fitness = "NSE", messages = "Message", categories = "Category") 
# geomOps <- mhplot::subsetByMessage(logMh)
# print(mhplot::plotParamEvolution(geomOps, pVarIds[1], objLims=c(0,1)))


# We need to weed out some catchment with no data point over the calibration period (or worse none at all...)
# so a couple of functions for it:
runoff_for_catchment <- function(cat_id, do_subset = FALSE) {
  d <- load_cat_series(cat_id, root_path=root_path)
  obs_runoff <- get_observed_runoff(d)
  if(do_subset) obs_runoff <- window(obs_runoff, start=simul_start, end=simul_end)
  obs_runoff
}

# x <- (runoff_for_catchment('225213'))

catchment_has_runoff <- function(cat_id) {
  obs_runoff <- runoff_for_catchment(cat_id, do_subset = TRUE)
  return(!(all(is.na(obs_runoff))))
}

