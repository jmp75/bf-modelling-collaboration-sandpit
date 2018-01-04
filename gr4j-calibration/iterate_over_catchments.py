# Use case

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

#########################
# Time series and dealing with CSV IO
#########################

def to_daily_time_series(data_array, start_date):
    tseries_index = pd.date_range(start_date, periods=len(data_array), freq='D')
    return pd.Series(data_array, index=tseries_index)

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
    rainfall_data = to_daily_time_series(tseries_df['Rain'].as_matrix(), start_date)
    # Sanity check: does time stamps match what I see via excel? 
    # rainfall_data.head()
    # rainfall_data.tail()
    pet_data = to_daily_time_series(tseries_df['Etp'].as_matrix(), start_date)
    obs_runoff_data = to_daily_time_series(tseries_df['Qobs'].as_matrix(), start_date)
    obs_runoff_data[obs_runoff_data < 0] = np.nan
    df = pd.concat([rainfall_data, pet_data, obs_runoff_data], axis=1)
    df.columns = ['Rain','Etp','Qobs'] 
    return df

#########################
# We define in this section some classes that will make their way into BF packages
# For now it is easier to have them inline in this file to match offline use cases
#########################

class runoff_lumped_spot_setup(object):
    """
    An adapter between BF GR4J (or potential other lumped models) and spotpy, to address use cases where:
    - one or more of the 4 model parameters can be optimized 
    - The objective function can be customized
    """
    def __init__(self, base_simulation, obj_calc, param_space, param_fixed=None):
        """
        :base_simulation:    An object of class RunoffYield, or something similar "a la Python"
        :obj_calc:  
        :param_space:  A list of spotpy 'parameter' objects.
        :param_fixed:  A dictionary of values for fixed parameters 
                e.g. {'x2': 0.0} complementing param_space
        """
        self.base_simulation = base_simulation
        self.obj_calc = obj_calc
        self.observed_runoff = obj_calc.observation
        self.param_space = param_space
        self.free_param_names = [ x.name for x in param_space ]
        self.param_fixed = param_fixed
        self.database = list()
    
    def parameters(self):
        return spotpy.parameter.generate(self.param_space)
    
    def simulation(self,vector):
        optim_param_space = dict(zip(self.free_param_names, vector))
        if self.param_fixed is None:
            params = dict(optim_param_space)
        else:
            # params = { **self.param_fixed, **optim_param_space } # python3.5+ way
            params = dict(list(self.param_fixed.items()) + list(optim_param_space.items())) # python 2..
        self.base_simulation.model_params = params
        runoff = self.base_simulation.execute_runoff()
        # runoff = gr4j.gr4j(self.rainfall,self.pet, x1=params['x1'], x2=params['x2'], x3=params['x3'], x4=params['x4']) 
        return runoff
    
    def evaluation(self):
        observations= self.observed_runoff
        return observations
    
    def objectivefunction(self, evaluation, simulation):
        # objectivefunction = -spotpy.objectivefunctions.nashsutcliffe(evaluation=observation,simulation=simulation)
        stat_val = self.obj_calc.bivariate_statistic(observation=evaluation, simulation=simulation)
        if self.obj_calc.is_maximisable:
            return -stat_val
        else:
            return stat_val

    def save(self, objectivefunctions, parameter, simulations):
        self.database.append([])


class Gr4jSimulation:
    """
    A wrapper around a simulation, more convenient and more generic for calibration purposes.
    """
    def __init__(self, rainfall, pet, model_params=None):
        """
            :rainfall:    rainfall series
            :type: numpy array

            :pet:    pet series
            :type: numpy array

            :model_params:  (optional) a dictionary with GR4J model parameters x 1 to 4 . 
        """
        self.rainfall = rainfall
        self.pet = pet
        if model_params is None:
            self.model_params = {
                'x1': 350.0,
                'x2':   0.0,
                'x3':  40.0,
                'x4':   0.5
            }
        else:
            self.model_params = model_params
    def apply_parameters(self, parameters):
        ks = list(set(['x1','x2','x3','x4']).intersection(parameters.keys()))
        for k in ks:
            self.model_params[k] = parameters[k]
    def execute_runoff(self):
        runoff = gr4j.gr4j(self.rainfall,self.pet, **(self.model_params)) 
        return runoff


class Objective:
    """
    A class that wraps bivariate statistics and indicates whether 
    it is a maximizable or minimizable objective in a calibration sense
    """
    def __init__(self, bivariave_stat_function, observation, series_subsetting=None, is_maximisable=False, name='objective'):
        """
            :bivariave_stat_function:    A function that takes two numpy arrays as arguments and returns a double
            :type: bivariate function

            :observation:    A Pandas series of observations to calibrate against
            :type: Pandas time series

            :series_subsetting:    (optional)A function that subsets a pandas time series to a numpy array before stat calculations
                                This function, if present, is applied to both observation and each simulation output.
            :type: univariate function

            :is_maximisable:    Is the bivariave_stat_function such that it is maximisable as an objective
            :type: bool

            :name:    Name of this objective.
            :type: str

        """
        self.biv_func = bivariave_stat_function
        self.is_maximisable = is_maximisable
        self.series_subsetting = series_subsetting
        self.name = name
        self.observation = self.subset_series(observation)
    def subset_series(self, series):
        if self.series_subsetting is None:
            return series.as_matrix() # assumes this is a pandas series
        else:
            return self.series_subsetting(series)
    def objective_statistic(self, simulation):
        return self.bivariate_statistic(self.observation, simulation)
    def bivariate_statistic(self, observation, simulation):
        simul_pd_series = to_daily_time_series(simulation, self.simul_start_date)
        simul_subset = self.subset_series(simul_pd_series)
        return self.biv_func(observation, simul_subset)

def get_best_in_log(log_results, obj_calc):
    """
    A custom processing of spotpy logs to retrieve adequate information for calibration results  
    
    :log_results: Expects an numpy array which should have as first axis an index "like" or "like1". 
    :type: array  
    
    :obj_calc: The objective used in the calibration. 
    :type: Objective  
    
    :return: Best parameter set and objective value
    :rtype: dictionary
    """
    # :maximize: Optional, default=True meaning the highest objectivefunction is taken as best, if False the lowest objectivefunction is taken as best.
    # :type: boolean
    
    maximize=True
    likes=log_results[spotpy.analyser.get_like_fields(log_results)[0]]
    if maximize:
        best=np.nanmax(likes)
    else:
        best=np.nanmin(likes)
    index=np.where(likes==best)
    best_parameter_set = spotpy.analyser.get_parameters(log_results[index])[0]
    parameter_names = spotpy.analyser.get_parameternames(log_results)
    x = dict(zip(parameter_names, best_parameter_set))
    x[obj_calc.name] = likes[index][0]
    return x

def create_simulation(pd_data, run_start, run_end):
    model_inputs = {
        'rainfall' : pd_data['Rain'][run_start:run_end].as_matrix(),
        'pet'  : pd_data['Etp'][run_start:run_end].as_matrix(),
    }
    base_simulation = Gr4jSimulation(model_inputs['rainfall'], model_inputs['pet'], None)
    # default_runoff = base_simulation.execute_runoff()
    # plt.show(default_runoff)
    return base_simulation

def calibrate_lumped_catchment(base_simulation, objective, param_space, param_fixed, max_iter):
    adapter = runoff_lumped_spot_setup(base_simulation, objective, param_space, param_fixed)
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


# Date range for calibration. 
run_start = datetime(1990, 1, 1)
run_end = datetime(2009, 12, 31)

def load_catchment_data(cat_id='226226'):
    csv_file = os.path.join(root_path, cat_id + '.csv')
    pd_data = load_csv_timeseries(csv_file)
    return pd_data

def subset_for_calib_stats(pd_series):
    subset_start = datetime(1993, 1, 1)
    subset_end = datetime(1999, 12, 31)
    return pd_series[subset_start:subset_end].as_matrix()

def subset_for_valid_stats(pd_series):
    subset_start = datetime(2000, 1, 1)
    subset_end = datetime(2009, 12, 31)
    return pd_series[subset_start:subset_end].as_matrix()


## Define the feasible parameter space and its complement the fix parameters
param_space = [
                spotpy.parameter.Uniform('x1', low=1.0, high=3000.0, optguess=300),
                # spotpy.parameter.Uniform('x2', low=-30.0, high=30.0, optguess=x2),
                spotpy.parameter.Uniform('x3', low=1.0, high=1000.0, optguess=40),
                spotpy.parameter.Uniform('x4', low=1, high=20.0, optguess=1)
            ]
param_fixed = {'x2': 0.0}
rep=1000
# rep=10000

def runnoff_with_params(base_simulation, best_pset):
    base_simulation.apply_parameters(best_pset)
    runoff = base_simulation.execute_runoff()
    return runoff

def calib_valid_catchment_id(cat_id):
    pd_data = load_catchment_data(cat_id)
    nse_obj = Objective(spotpy.objectivefunctions.nashsutcliffe, pd_data['Qobs'], subset_for_calib_stats, is_maximisable=True, name="NSE")
    nse_valid = Objective(spotpy.objectivefunctions.nashsutcliffe, pd_data['Qobs'], subset_for_valid_stats, is_maximisable=True, name="NSE")
    # Due to both BF GR4J and probably spotpy not using/supporting 
    # fully fledged pandas time series, we have to hack things a bit to 
    # have a decent subsetting, and put that information into the Objective object.
    nse_obj.simul_start_date = run_start
    nse_valid.simul_start_date = run_start
    base_simulation = create_simulation(pd_data, run_start, run_end)
    best_pset = calibrate_lumped_catchment(base_simulation, nse_obj, param_space, param_fixed, max_iter=rep)
    runoff = runnoff_with_params(base_simulation, best_pset)
    valid_obj = nse_valid.objective_statistic(runoff)
    p = dict(best_pset)
    p['NSE_valid'] = valid_obj
    p['cat_id'] = cat_id
    return p

###
# Batch calibration/verification
###

cat_ids = ['226226','403210']
# 225110.csv  225219.csv  226204.csv  229661.csv  403213.csv  403222.csv  403226.csv  403244.csv  404208.csv
# 225213.csv  226007.csv  226226.csv  403210.csv  403217.csv  403224.csv  403232.csv  404207.csv  405205.csv

results = []
for cat_id in cat_ids:
    try:
        p = calib_valid_catchment_id(cat_id)
        results.append(p)
    except:
        continue #This is not best practice in general...

calib_valid_df = pd.DataFrame.from_records(results)

###########################
# Below is temporary code used for interactive "testing"
# adapter = runoff_lumped_spot_setup(base_simulation, nse_obj, param_space, param_fixed)
# generated_params = adapter.parameters()
# test_values = [x[0] for x in generated_params] 
# simul_outputs = adapter.simulation(test_values)
# observations = adapter.evaluation()
# objective_f = adapter.objectivefunction(pd_data['Qobs'][warmup_to:run_end].as_matrix(), simul_outputs)

# test_values = list(best.values())[0:3] 
# simul_outputs = adapter.simulation(test_values)
# observations = adapter.evaluation()
# calibration_inputs = {'obs_runoff': pd_data['Qobs'][warmup_to:run_end].as_matrix()}
# objective_f = adapter.objectivefunction(calibration_inputs['obs_runoff'], simul_outputs)
