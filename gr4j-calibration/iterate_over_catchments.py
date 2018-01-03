# Use case
# Import CSV data

# import profile

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import spotpy
import sys
import os
import json
import sys


if(sys.platform == 'win32'):
    root_src = r'C:\src\csiro\stash\bf'
    root_path = r'C:\Users\per202\Documents\BF\AUS_Catchments\AUS_Catchments'
else:
    root_src = '/home/per202/src/csiro/stash/per202/bf'
    root_path = '/home/per202/Data/BF/AUS_Catchments'

sys.path.append(os.path.join(root_src, 'gr4j-sg'))
sys.path.append(os.path.join(root_src,'reach'))
sys.path.append(os.path.join(root_src,'watercloud-server'))
sys.path.append(os.path.join(root_src,'watercloud-server/endpoints'))

import gr4j

csv_file = os.path.join(root_path,'226226.csv')

tseries_df = pd.read_csv(csv_file)
# TODO
# stopifnot expected colnames 

start_year = tseries_df['year'][0]
start_month = tseries_df['month'][0]
start_day = tseries_df['date'][0]  # date ?? anyway...

tseries_len = tseries_df.shape[0] 

## NOTE: use yyyy-mm-dd for date format. Pandas seems to use the silly US format for dd/mm/yyyy even on machines with en_AU locales. 
rng = pd.date_range(str(start_year) + '-' + str(start_month) + '-' + str(start_day), periods=tseries_len, freq='D')

# "Fun" facts: a pd.Series object cannot be at argument to pd.Series, it seems
rainfall_data = pd.Series(tseries_df['Rain'].as_matrix(), index=rng)

# Sanity check: does time stamps match what I see via excel? 
rainfall_data.head()
rainfall_data.tail()

pet_data = pd.Series(tseries_df['Etp'].as_matrix(), index=rng)
obs_runoff_data = pd.Series(tseries_df['Qobs'].as_matrix(), index=rng)

obs_runoff_data[obs_runoff_data < 0] = np.nan

# Date range for calibration. 

from datetime import datetime

run_start = datetime(1980, 1, 1)
warmup_to = datetime(1982, 1, 1)
run_end = datetime(1999, 12, 31)


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
        base_simulation.model_params = params
        runoff = base_simulation.execute_runoff()
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

import numpy as np

class Gr4jSimulation:
    def __init__(self, rainfall, pet, model_params=None):
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
    def execute_runoff(self):
        runoff = gr4j.gr4j(self.rainfall,self.pet, **(self.model_params)) 
        return runoff

class Objective:
    """
    A class that wraps bivariate statistics and indicates whether 
    it is a maximizable or minimizable objective in a calibration sense
    """
    def __init__(self, bivariave_stat_function, observation, is_maximisable=False, name='objective'):
        self.biv_func = bivariave_stat_function
        self.is_maximisable = is_maximisable
        self.observation = observation
        self.name = name
    def objective_statistic(self, simulation):
        return self.biv_func(self.observation, simulation)
    def bivariate_statistic(self, observation, simulation):
        return self.biv_func(observation, simulation)


model_inputs = {
    'rainfall' : rainfall_data[run_start:run_end].as_matrix(),
    'pet'  : pet_data[run_start:run_end].as_matrix(),
}
calibration_inputs = {'obs_runoff'  : obs_runoff_data[run_start:run_end].as_matrix() }


base_simulation = Gr4jSimulation(model_inputs['rainfall'], model_inputs['pet'], None)
default_runoff = base_simulation.execute_runoff()
# plt.show(default_runoff)

nse_obj = Objective(spotpy.objectivefunctions.nashsutcliffe, calibration_inputs['obs_runoff'], is_maximisable=True, name="NSE")

param_space = [
                spotpy.parameter.Uniform('x1', low=1.0, high=3000.0, optguess=300),
                # spotpy.parameter.Uniform('x2', low=-30.0, high=30.0, optguess=x2),
                spotpy.parameter.Uniform('x3', low=1.0, high=1000.0, optguess=40),
                spotpy.parameter.Uniform('x4', low=1, high=20.0, optguess=1)
            ]

param_fixed = {'x2': 0.0}

adapter = runoff_lumped_spot_setup(base_simulation, nse_obj, param_space, param_fixed)

generated_params = adapter.parameters()
test_values = [x[0] for x in generated_params] 
simul_outputs = adapter.simulation(test_values)
observations = adapter.evaluation()
objective_f = adapter.objectivefunction(calibration_inputs['obs_runoff'], simul_outputs)

sampler=spotpy.algorithms.sceua(adapter, alt_objfun = None, save_sim = False, dbname='Grrr', dbformat='ram')
# I am supposed to be able to have my own logger but the following leads to a fail
# sampler=spotpy.algorithms.sceua(adapter)

rep=500
rep=2000
# rep=10000
sampler.sample(rep)
results = sampler.getdata()
# spotpy.analyser.plot_parametertrace(results)

# best_pset = spotpy.analyser.get_best_parameterset(results)
# print(best_pset)

# for cat_num in cat_identifiers:

def get_best_in_log(log_results, obj_calc):
    """
    Get the best parameter set of your result array, depending on your first objectivefunction 
    
    :log_results: Expects an numpy array which should have as first axis an index "like" or "like1". 
    :type: array  
    
    :maximize: Optional, default=True meaning the highest objectivefunction is taken as best, if False the lowest objectivefunction is taken as best.
    :type: boolean
    
    :return: Best parameter set and objective value
    :rtype: dictionary
    """
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


best = get_best_in_log(results, nse_obj)
print(best)
test_values = list(best.values())[0:3] 
simul_outputs = adapter.simulation(test_values)
observations = adapter.evaluation()
objective_f = adapter.objectivefunction(calibration_inputs['obs_runoff'], simul_outputs)

