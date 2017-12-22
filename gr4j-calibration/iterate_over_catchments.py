# Use case
# Import CSV data

root_path = r'C:\Users\per202\Documents\BF\AUS_Catchments\AUS_Catchments'

import pandas as pd
import numpy as np
import spotpy
import sys
import os

if(sys.platform == 'win32'):
    sys.path.append(r'C:\src\github_jm\gr4j-sg')
else:
    sys.path.append('/home/per202/src/csiro/stash/per202/bf/gr4j-sg')


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


model_inputs = {
    'rainfall' : rainfall_data[run_start:run_end].as_matrix(),
    'pet'  : pet_data[run_start:run_end].as_matrix(),
    'x1':20.0,
    'x2':0.0,
    'x3':40.0,
    'x4':0.5
}

runoff = gr4j.gr4j(**model_inputs)

import matplotlib.pyplot as plt

# p = plt.plot(runoff)
# plt.show(p)



class gr4j_lumped_spot_setup(object):
    def __init__(self, rainfall, pet, observed_runoff):
        self.params = [spotpy.parameter.Uniform('x1', low=1.0, high=1000.0, optguess=350),
                       spotpy.parameter.Uniform('x2', low=-30.0, high=30.0, optguess=15),
                       spotpy.parameter.Uniform('x3', low=1.0, high=1000.0, optguess=100),
                       spotpy.parameter.Uniform('x4', low=1, high=240.0, optguess=10)
                       ]
        self.observed_runoff = observed_runoff
        self.rainfall = rainfall
        self.pet = pet
    
    def parameters(self):
        return spotpy.parameter.generate(self.params)        
    
    def simulation(self,vector):
        runoff = gr4j.gr4j(self.rainfall,self.pet, x1=vector[0], x2=vector[1], x3=vector[2], x4=vector[3]) 
        return runoff
        
    def evaluation(self):
        observations= self.observed_runoff
        return observations
    
    def objectivefunction(self,simulation,evaluation):
        objectivefunction = -spotpy.objectivefunctions.nashsutcliffe(evaluation=evaluation,simulation=simulation)
        return objectivefunction



calibration_inputs = {'obs_runoff'  : obs_runoff_data[run_start:run_end].as_matrix() }
spotpy_setup=gr4j_lumped_spot_setup(model_inputs['rainfall'], model_inputs['pet'], calibration_inputs['obs_runoff'])

spotpy_setup.objectivefunction(runoff, calibration_inputs['obs_runoff'])

rep=100
sampler=spotpy.algorithms.sceua(spotpy_setup, dbname='GR4J_test', dbformat='ram')

sampler.sample(rep)
results = sampler.getdata()
spotpy.analyser.plot_parametertrace(results)

# The following returns something that looks nonsense to me.
print(spotpy.analyser.get_best_parameterset(results))



# for cat_num in cat_identifiers:
