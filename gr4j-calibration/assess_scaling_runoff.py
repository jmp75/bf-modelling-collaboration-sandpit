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
    root_path = '/home/per202/data/bf/AUS_Catchments_Inputs'

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

# def create_simulation(pd_data, run_start, run_end): now in /watercloud-server/simulation/runoff.py

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
    best = get_best_in_log(results, objective.name)
    return best

def load_catchment_data(cat_id='226226'):
    csv_file = os.path.join(root_path, cat_id + '.csv')
    pd_data = load_csv_timeseries(csv_file)
    return pd_data

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

subset_for_calib_stats = ts_subsetter(datetime(1993, 1, 1), datetime(1999, 12, 31))
subset_for_valid_stats = ts_subsetter(datetime(2000, 1, 1), datetime(2009, 12, 31))


def calib_valid_catchment_id(cat_id, run_start = default_run_start, run_end = default_run_end, param_space = default_param_space, param_fixed=default_param_fixed, rep=default_max_calib_iter, objective_fun=nse_log_bias_nan, objective_name="nse_log_bias", scale_obs_factor=1.0):
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
    if scale_obs_factor != 1.0:
        pd_data['Qobs'] = pd_data['Qobs'] * scale_obs_factor
    nse_obj = Objective(objective_fun, pd_data['Qobs'], subset_for_calib_stats, is_maximisable=True, name=objective_name)
    nse_valid = Objective(objective_fun, pd_data['Qobs'], subset_for_valid_stats, is_maximisable=True, name=objective_name)
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
    p[objective_name + '_valid'] = valid_obj
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
# # Batch calibration/verification
# ###

nse_monthly_means = compose_bivariate(mean_monthly, nse_log_bias_nan)
# results = []
# p = calib_valid_catchment_id('603002', rep=200, objective_fun = nse_monthly_means, objective_name="MM_nselb")
# results.append(p)

cat_ids = ['603002', '603003', '603136','604001','604053','605012','606001','606185','606195','606218','607002','607003','607004',
           '607007','607013','607144','607155','607220','608001','608151','608171','609002','609003','609005','609016','609017',
           '609018','609022','610001','610009','610010','610015','610219','611004','611006','611111','612001','612023','612025',
           '612034','612230','613002','614123','614196']


def batch_clibrate(cat_ids, runoff_scaling, output_filename):
    results = []
    counter = 0
    for cat_id in cat_ids:
        counter = counter + 1
        print("INFO: calibration " + str(counter) + " out of " + str(len(cat_ids)) + ": " + cat_id)
        try:
            p = calib_valid_catchment_id(cat_id, rep=3000, objective_fun = nse_monthly_means, objective_name="MM_nselb", scale_obs_factor=runoff_scaling)
            results.append(p)
        except Exception as e:
            print("ERROR: calibration failed for " + cat_id + " " + str(e))
            continue #This is not best practice in general...
    calib_valid_df = pd.DataFrame.from_records(results)
    calib_valid_df.to_csv(output_filename)

batch_clibrate(cat_ids, 1.0, '~/tmp/calib_valid_noscaling.csv')
batch_clibrate(cat_ids, 0.01, '~/tmp/calib_valid_scaling.csv')



csv_file = '~/tmp/calib_valid_noscaling.csv'
parameters_df = pd.read_csv(csv_file)
# https://stackoverflow.com/questions/22231592/pandas-change-data-type-of-series-to-string
parameters_df.cat_id = parameters_df.cat_id.apply(str)
# Otherwise later on:TypeError: invalid type comparison


# def plot_daily_observed_modelled_runoff(cat_id, transform_output_series=no_transform):
cat_id = '603003'



pd_data = load_catchment_data(cat_id)
obs_runoff = subset_for_calib_stats(pd_data['Qobs'])

# nse_log_bias_obj = Objective(nse_log_bias_nan, pd_data['Qobs'], subset_for_calib_stats, is_maximisable=True, name="NSE_log_bias")
# nse_log_bias_valid = Objective(nse_log_bias_nan, pd_data['Qobs'], subset_for_valid_stats, is_maximisable=True, name="NSE_log_bias")
# nse_log_bias_obj.simul_start_date = default_run_start
# nse_log_bias_valid.simul_start_date = default_run_start

base_simulation = create_simulation(pd_data, default_run_start, default_run_end)

# https://stackoverflow.com/questions/19237878/subsetting-a-python-dataframe
# k1 = parameters_df.loc[(parameters_df.cat_id == cat_id), ['Time', 'Product']]
k1 = parameters_df.loc[(parameters_df.cat_id == cat_id)]
model_params = {
        'x1': 350.0,
        'x2':   0.0,
        'x3':  40.0,
        'x4':   0.5
}
ks = list(set(['x1','x2','x3','x4']).intersection(k1.to_dict().keys()))
for k in ks:
    model_params[k] = k1[k][0]


def transform_output_series(s):
    return mean_monthly(subset_for_calib_stats(s))

mod_runoff = runoff_with_params(base_simulation, model_params)

transform_output_series(obs_runoff)

obs_mod_runoff = concat_pandas_series(['mod_runoff','obs_runoff'], transform_output_series(mod_runoff), transform_output_series(obs_runoff))

obs_mod_runoff.cat_id = cat_id

obs_mod_runoff['month'] = np.arange(1, 13)


graph = obs_mod_runoff.plot(title = 'Graph of monthly observed and modelled runoff in mm/day').set(xlabel ='year', ylabel = 'runoff (mm/day))')
return graph

