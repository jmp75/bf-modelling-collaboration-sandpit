library(readr)
library(stringr)

library(uchronia)
library(swift)


## Define functions to help concise batch optimizations

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

get_gr4j_parameter_space <- function() {
  pSpecGr4j <- joki::getFreeParams('GR4J')
  pSpecGr4j$Value <- c(542.1981111, 0, 7.7403390, 1.2388548)
  pSpecGr4j$Min <- c(1,-30, 1,1)
  pSpecGr4j$Max <- c(3000, 30, 1000, 10)
  pSpecGr4j$Name <- paste0(model_property_prefix, pSpecGr4j$Name)
  p <- createParameterizer(type='Generic', pSpecGr4j)
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

create_calib_data <- function(simulation, model_property_id, simul_start, warmup_days = 365, simul_end, obs_runoff, obj_id = 'NSE') {
  simulation <- CloneModel_R(simulation)
  w <- simul_start + lubridate::days(warmup_days)
  setSimulationSpan(simulation, sim_start, simul_end)
  recordState(simulation, model_property_id)
  calib_inputs = list(
    simulation = simulation,
    simul_start = simul_start,
    warmup_to = w,
    simul_end = simul_end,
    obj_id = obj_id,
    observation = obs_runoff,
    model_config = get_gr4j_parameter_space(),
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
  optimizer <- createSceOptimSwift(obj, term, SCEpars=sceParams, urs)
  calibLogger <- setCalibrationLogger(optimizer, '')
  startTime <- lubridate::now()
  calibResults <- executeOptimization(optimizer)
  endTime <- lubridate::now()
  calibWallTime <- endTime-startTime
  print(paste( 'Optimization completed in ', calibWallTime, attr(calibWallTime, 'units')))
  return(list(Results=calibResults, Optimizer=optimizer))
}

## Using functions above. compose a batch calibration workflow. Start with one.

root_path <- 'C:/Users/per202/Documents/BF/AUS_Catchments/AUS_Catchments'

cat_id <- "404207"

d <- load_cat_series(cat_id, root_path=root_path)
simulation <- create_lumped_rrmodel(d)
obs_runoff <- get_observed_runoff(d)
sim_start <- ISOdate(2000,1,1)
plot_obs_vs_calc(calib_inputs, sim_start)

simul_start <- ISOdate(1990,1,1)
simul_end <- simul_start + lubridate::years(10) - lubridate::days(1)
warmup_days <- 365

model_property <- paste0( model_property_prefix, 'runoff')
if(all(is.na(obs_runoff))) stop("Not possible to calibrate without runoff observations")
calib_data <- create_calib_data(simulation, model_property, simul_start, warmup_days, simul_end, obs_runoff, obj_id = 'NSE')
obj <- create_objective (calib_data)
getScore(obj, get_gr4j_parameter_space())

optim_results <- calibrate(calib_data)

calibResults <- optim_results$Results
optimizer <- optim_results$Optimizer

sortedResults <- sortByScore(calibResults, 'NSE')
results_df <- scoresAsDataFrame(sortedResults)
head(results_df)

batch_results <- list()
batch_results[[cat_id]] = results_df[1,]

# optional, details viz look at calibration logs. 
# d <- getLoggerContent(optimizer)
# d$PointNumber = 1:nrow(d)
# logMh <- mhplot::mkOptimLog(d, fitness = "NSE", messages = "Message", categories = "Category") 
# geomOps <- mhplot::subsetByMessage(logMh)
# print(mhplot::plotParamEvolution(geomOps, pVarIds[1], objLims=c(0,1)))

# TODO: batch loop. Note that we need to check observed data as some catchment have no runoff observations.
short_fnames <- list.files(root_path, pattern='*.csv')
cat_ids <- stringr::str_replace_all(short_fnames, '\\.csv', '')
