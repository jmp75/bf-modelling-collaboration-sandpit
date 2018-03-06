import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import spotpy
import sys
import os
import json
import sys
from datetime import datetime


if sys.platform == 'win32':
    root_src = r'C:\src\csiro\stash\bf'
    root_path = r'C:\Users\per202\Documents\BF\AUS_Catchments\AUS_Catchments'
else:
    root_src = '/home/per202/src/csiro/stash/per202/bf'
    root_path = '/home/per202/Data/BF/AUS_Catchments'

sys.path.append(os.path.join(root_src, 'gr4j-sg'))
sys.path.append(os.path.join(root_src, 'reach'))
sys.path.append(os.path.join(root_src, 'watercloud-server'))
sys.path.append(os.path.join(root_src, 'watercloud-server/endpoints'))

import gr4j
from calibration.calib_tasks import *
from calibration.statistics_bf import *
from simulation.runoff import *

#########################
# Time series and dealing with CSV IO
#########################

# Following functions moved to main BF codebase
# def to_daily_time_series(data_array, start_date):
# spotpy and gr4j work, for better or worse, with numpy arrays without date time indices. 
# A function to capture the conversion, which may grow to doing more than as_matrix
# def to_numpy_array(pd_series):
# def concat_pandas_series(colnames, *series):

def load_csv_timeseries(csv_file):
    tseries_df = pd.read_csv(csv_file)
    # TODO
    # stopifnot expected colnames
    start_year = tseries_df['year'][0]
    start_month = tseries_df['month'][0]
    start_day = tseries_df['date'][0]  # date ?? anyway...
    tseries_len = tseries_df.shape[0] 
    ## NOTE: use yyyy-mm-dd for date format. Pandas seems to use the silly US format instead of dd/mm/yyyy even on machines with en_AU locales. 
    start_date = str(start_year) + '-' + str(start_month) + '-' + str(start_day)
    # "Fun" facts: a pd.Series object cannot be at argument to pd.Series, it seems
    # nevertheless there has got to be simplifications possible for the next:
    rainfall_data = to_daily_time_series(to_numpy_array(tseries_df['Rain']), start_date)
    # Sanity check: does time stamps match what I see via excel? 
    # rainfall_data.head()
    # rainfall_data.tail()
    pet_data = to_daily_time_series(to_numpy_array(tseries_df['Etp']), start_date)
    obs_runoff_data = to_daily_time_series(to_numpy_array(tseries_df['Qobs']), start_date)
    obs_runoff_data[obs_runoff_data < 0] = np.nan
    return concat_pandas_series( ['Rain','Etp','Qobs'], rainfall_data, pet_data, obs_runoff_data)

def create_simulation(pd_data, run_start, run_end):
    """
    Create a 'simulation' object (type Gr4jSimulation for now) given daily input climare series and a time span
    
    :pd_data: Pandas data frame with columns 'Rain' and 'Etp'
    :type: Pandas data frame
    
    :run_start: simulation start
    :type: datetime
    
    :run_end: simulation end
    :type: datetime
    
    :return: a Simulation object ready to be run given parameter sets
    :rtype: Gr4jSimulation
    """
    model_inputs = {
        'rainfall' : to_numpy_array(pd_data['Rain'][run_start:run_end]),
        'pet'  : to_numpy_array(pd_data['Etp'][run_start:run_end])
    }
    base_simulation = Gr4jSimulation(model_inputs['rainfall'], model_inputs['pet'], run_start, None)
    # default_runoff = base_simulation.execute_runoff()
    # plt.show(default_runoff)
    return base_simulation

def calibrate_lumped_catchment(base_simulation, objective, param_space, param_fixed, max_iter):
    adapter = LumpedCalibrationDefinition(base_simulation, objective, param_space, param_fixed)
    sampler=spotpy.algorithms.sceua(adapter, alt_objfun = None, save_sim = False, dbname='Grrr', dbformat='ram')
    # I am supposed to be able to have my own logger but the following leads to a fail
    # sampler=spotpy.algorithms.sceua(adapter)
    sampler.sample(max_iter)
    results = sampler.getdata()
    # spotpy.analyser.plot_parametertrace(results)
    # best_pset = spotpy.analyser.get_best_parameterset(results)
    # print(best_pset)
    best = get_best_in_log(results, objective)
    return best

def load_catchment_data(cat_id='226226'):
    csv_file = os.path.join(root_path, cat_id + '.csv')
    pd_data = load_csv_timeseries(csv_file)
    return pd_data

subset_for_calib_stats = ts_subsetter(datetime(1993, 1, 1), datetime(1999, 12, 31))
subset_for_valid_stats = ts_subsetter(datetime(2000, 1, 1), datetime(2009, 12, 31))

def runoff_with_params(base_simulation, parameters):
    base_simulation.apply_parameters(parameters)
    runoff = base_simulation.execute_runoff()
    return runoff

def plot_modelled_observed(simulation, parameters):
    pass

######
# Default arguments
######
# default simulation span
default_run_start = datetime(1990, 1, 1)
default_run_end = datetime(2009, 12, 31)
## Define the feasible parameter space and its complement the fix parameters
default_param_space = [
                spotpy.parameter.Uniform('x1', low=1.0, high=3000.0, optguess=300),
                # spotpy.parameter.Uniform('x2', low=-30.0, high=30.0, optguess=x2),
                spotpy.parameter.Uniform('x3', low=1.0, high=1000.0, optguess=40),
                spotpy.parameter.Uniform('x4', low=1, high=20.0, optguess=1)
            ]
default_param_fixed = {'x2': 0.0}
default_max_calib_iter=1000



def calib_valid_catchment_id(cat_id, run_start = default_run_start, run_end = default_run_end, param_space = default_param_space, param_fixed=default_param_fixed, rep=default_max_calib_iter):
    """
    Perform a calibration and validation for a catchment
    
    :cat_id: catchment id, e.g. '226226'
    :type: string
    
    :run_start: simulation start
    :type: datetime
    
    :run_end: simulation end
    :type: datetime
    
    :param_space: feasible parameter space for the calibration
    :type: a list of spotpy Parameter objects
    
    :param_fixed: fixed gr4j parameters
    :type: dictionary e.g. {'x2': 0.0}
    
    :rep: maximum number of iterations in the calibration
    :type: positive integer
    
    :return: composite of values for model parameters, catchment id, calibration and validation values.
    :rtype: a dictionary
    """
    pd_data = load_catchment_data(cat_id)
    nse_obj = Objective(nse_nan, pd_data['Qobs'], subset_for_calib_stats, is_maximisable=True, name="NSE")
    nse_valid = Objective(nse_nan, pd_data['Qobs'], subset_for_valid_stats, is_maximisable=True, name="NSE")
    # Due to both BF GR4J and probably spotpy not using/supporting 
    # fully fledged pandas time series, we have to hack things a bit to 
    # have a decent subsetting, and put that information into the Objective object.
    nse_obj.simul_start_date = run_start
    nse_valid.simul_start_date = run_start
    base_simulation = create_simulation(pd_data, run_start, run_end)
    best_pset = calibrate_lumped_catchment(base_simulation, nse_obj, param_space, param_fixed, max_iter=rep)
    runoff = runoff_with_params(base_simulation, best_pset)
    valid_obj = nse_valid.objective_statistic(runoff)
    p = dict(best_pset)
    p['NSE_valid'] = valid_obj
    p['cat_id'] = cat_id
    return p

def plot_observed_modelled(simulation, parameters, observation, from_date = default_run_start, to_date = default_run_end):
    r = runoff_with_params(simulation, parameters)
    r = r[from_date:to_date]
    obs = observation[from_date:to_date]
    obs_mod = concat_pandas_series(['mod','obs'], r, obs)
    # matplotlib has no simple default ways to label axes and titles???? 
    # fig, ax = plt.subplots()
    my_plot = plt.plot(obs_mod)
    return(my_plot)


# ###
# # Plotting side by side simulation outputs and observations
# ###

obs_226226 = load_catchment_data('226226')

simulation_226226 = create_simulation(obs_226226, default_run_start, default_run_end)
r = simulation_226226.execute_runoff()
obs_runoff_226226 = obs_226226['Qobs'][default_run_start:default_run_end]
obs_mod = concat_pandas_series(['mod','obs'], r, obs_runoff_226226)
# matplotlib has no simple default ways to label axes and titles???? 
# fig, ax = plt.subplots()
my_plot = plt.plot(obs_mod)
plt.show()

# A data frame as we'd typically get from batch calibrations. 
calib_valid_df = pd.DataFrame.from_records([
    { 'NSE' : 0.672175, 'NSE_valid' : 0.564606, 'cat_id': 226226, 'x1' : 690.653397, 'x2' : 0.0, 'x3' : 54.068796, 'x4' :  1.358228},
    { 'NSE' : 0.772028, 'NSE_valid' : 0.795753, 'cat_id': 403210, 'x1' : 722.697191, 'x2' : 0.0, 'x3' :191.163011, 'x4' :  1.99661},
    { 'NSE' : 0.803792, 'NSE_valid' : 0.647178, 'cat_id': 229661, 'x1' : 711.087971, 'x2' : 0.0, 'x3' :185.898297, 'x4' :  2.98321   }])


parameters = calib_valid_df[['x1','x2','x3','x4']].iloc[[0]].to_dict()

my_plot = plot_observed_modelled(simulation_226226, parameters, obs_runoff_226226, from_date = default_run_start, to_date = default_run_end)
plt.show()

my_plot = plot_observed_modelled(simulation_226226, parameters, obs_runoff_226226, from_date = datetime(1993, 1, 1), to_date = datetime(1999,12, 31))
plt.show()

# ###
# # Batch calibration/verification
# ###
# p = calib_valid_catchment_id('226226')


# cat_ids = ['226226','403210','229661']
# # 225110.csv  225219.csv  226204.csv  229661.csv  403213.csv  403222.csv  403226.csv  403244.csv  404208.csv
# # 225213.csv  226007.csv  226226.csv  403210.csv  403217.csv  403224.csv  403232.csv  404207.csv  405205.csv

# results = []
# for cat_id in cat_ids:
#     try:
#         p = calib_valid_catchment_id(cat_id)
#         results.append(p)
#     except:
#         print("ERROR: calibration failed for " + cat_id)
#         continue #This is not best practice in general...

# calib_valid_df = pd.DataFrame.from_records(results)

