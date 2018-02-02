
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import json
import sys
from datetime import datetime


# In[2]:


from __future__ import division


# In[3]:


if sys.platform == 'win32':
    root_src = r'C:\src\csiro\stash\bf'
    root_path = r'C:\Users\per202\Documents\BF\AUS_Catchments\AUS_Catchments'
else:
    root_src = '/home/per202/src/csiro/stash/per202/bf'
    root_path = '/home/per202/Data/BF/AUS_Catchments'

root_src = '.'
root_path = '/home/jovyan/work/AUS_Catchments/AUS_Catchments_Inputs'


# In[4]:


sys.path.append(os.path.join(root_src, 'spotpy'))
sys.path.append(os.path.join(root_src, 'gr4j-sg'))
sys.path.append(os.path.join(root_src, 'reach'))
sys.path.append(os.path.join(root_src, 'watercloud-server'))
sys.path.append(os.path.join(root_src, 'watercloud-server/endpoints'))


# In[5]:


import gr4j
import spotpy


# In[6]:


#########################
# Time series and dealing with CSV IO
#########################

def to_daily_time_series(data_array, start_date):
    tseries_index = pd.date_range(start_date, periods=len(data_array), freq='D')
    return pd.Series(data_array, index=tseries_index)

# spotpy and gr4j work, for better or worse, with numpy arrays without date time indices. 
# A function to capture the conversion, which may grow to doing more than as_matrix
def to_numpy_array(pd_series):
    if pd_series is np.array:
        return pd_series
    else:
        return pd_series.as_matrix()  # Assumption...

def concat_pandas_series(colnames, *series):
    df = pd.concat(series, axis=1)
    df.columns = colnames
    return df

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


# In[7]:


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
    def __init__(self, rainfall, pet, start_date, model_params=None):
        """
            :rainfall:    rainfall series
            :type: numpy array
            :pet:    pet series
            :type: numpy array
            :model_params:  (optional) a dictionary with GR4J model parameters x 1 to 4 . 
        """
        self.rainfall = rainfall
        self.pet = pet
        self.start_date = start_date
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
        return to_daily_time_series(runoff, self.start_date)


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
            return to_numpy_array(series) # assumes this is a pandas series
        else:
            return self.series_subsetting(series)
    def objective_statistic(self, simulation):
        return self.bivariate_statistic(self.observation, simulation)
    def bivariate_statistic(self, observation, simulation):
        if simulation is pd.Series:
            simul_pd_series = simulation
        else:
            simul_pd_series = to_daily_time_series(simulation, self.simul_start_date)
        simul_subset = self.subset_series(simul_pd_series)
        return self.biv_func(to_numpy_array(observation), to_numpy_array(simul_subset))

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


# In[8]:


######
# Default arguments
######
# default simulation span
default_run_start = datetime(1990, 1, 1)
default_run_end = datetime(2009, 12, 31)

def load_catchment_data(cat_id):
    csv_file = os.path.join(root_path, cat_id + '.csv')
    pd_data = load_csv_timeseries(csv_file)
    return pd_data

def to_monthly(pd_series):
    monthly_pd_series = pd_series.resample("M").sum()
    return monthly_pd_series

def subset_for_calib_stats(pd_series):
    subset_start = datetime(1993, 1, 1)
    subset_end = datetime(1999, 12, 31)
    x = pd_series[subset_start:subset_end]
    monthly_pd_series = to_monthly(x)
    return monthly_pd_series

def subset_for_valid_stats(pd_series):
    subset_start = datetime(2000, 1, 1)
    subset_end = datetime(2009, 12, 31)
    x = pd_series[subset_start:subset_end]
    monthly_pd_series = to_monthly(x)
    return monthly_pd_series


# In[9]:


def runoff_with_params(base_simulation, parameters):
    base_simulation.apply_parameters(parameters)
    runoff = base_simulation.execute_runoff()
    return runoff

def plot_modelled_observed(simulation, parameters):
    pass

def censor_missing_observations(observed, modelled):
    valid_indices = np.where(~np.isnan(observed))
    return (observed[valid_indices], modelled[valid_indices])

def relative_bias(observed, modelled):
    if observed.sum() <= 0:
        return 'Error'
    if modelled.sum() <= 0:
        return 'Error'
    relative_bias = (modelled.sum() - observed.sum())/observed.sum()
    return relative_bias

def nse_log_bias(observed, modelled):
    """
    NSE - log bias objective function.
    The usefulness of bias constraints in model calibration 
    for regionalisation to ungauged catchments, 18th World IMACS / MODSIM Congress, Cairns, Australia 13-17 July 2009
    https://www.mssanz.org.au/modsim09/I7/viney_I7a.pdf
        .. math::
        nse_log_bias = nse - 5.0 * (abs(log(e_{1 + bias}))) ** 2.5
    :observed: Observed data to be compared with modelled data
    :type: list
    :modelled: Modelled data to be compared with observed data
    :type: list
    :return: Nash-Sutcliffe Effiency with log bias constraint
    :rtype: float
    """
    if len(observed) == len(modelled):
        nse = spotpy.objectivefunctions.nashsutcliffe(observed, modelled)
        bias = relative_bias(observed, modelled)
        return nse - 5.0 * (abs(np.log(1 + bias))) ** 2.5
    else:
        logging.warning("evaluation and simulation lists does not have the same length.")
        return np.nan

# functions dealing with missing values:

def nse_nan(observed, modelled):
    observed, modelled = censor_missing_observations(observed, modelled)
    return spotpy.objectivefunctions.nashsutcliffe(observed, modelled)

def nse_log_bias_nan(observed, modelled):
    observed, modelled = censor_missing_observations(observed, modelled)
    return nse_log_bias(observed, modelled)


# In[10]:


## Define the feasible parameter space and its complement the fix parameters
# Default values and ranges as described in eWater Source
# default_param_space = [
#                 spotpy.parameter.Uniform('x1', low=1.0, high=1500.0, optguess=350.0),
#                 # spotpy.parameter.Uniform('x2', low=-10.0, high=5.0, optguess=0.0),
#                 spotpy.parameter.Uniform('x3', low=1.0, high=500.0, optguess=40.0),
#                 spotpy.parameter.Uniform('x4', low=0.5, high=4.0, optguess=0.5)
#             ]

# Median values and approximate 80% confidence intervals as described in Perrin et al. 2013 
default_param_space = [
                spotpy.parameter.Uniform('x1', low=100.0, high=1200.0, optguess=350.0),
                spotpy.parameter.Uniform('x2', low=-5.0, high=3.0, optguess=0.0),
                spotpy.parameter.Uniform('x3', low=20.0, high=300.0, optguess=90.0),
                spotpy.parameter.Uniform('x4', low=1.1, high=2.9, optguess=1.7)
            ]

default_param_fixed = {}
default_max_calib_iter=1500


# In[11]:


## Catchments of interest

# cat_ids = ['226204']
# old_cat_ids = ['225110','225213','225219','226007','226204','226226','229661','403210','403213','403217','403222','403224','403226','403232','403244','404207','404208','405205']

# new_cat_ids
# cat_ids = ['603002', '603003', '603136','604001','604053','605012','606001','606185','606195','606218','607002','607003','607004',
#            '607007','607013','607144','607155','607220','608001','608151','608171','609002','609003','609005','609016','609017',
#            '609018','609022','610001','610009','610010','610015','610219','611004','611006','611111','612001','612023','612025',
#            '612034','612230','613002','614123','614196']

cat_ids = ['603003']


# In[12]:


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
    nse_log_bias_obj = Objective(nse_log_bias_nan, pd_data['Qobs'], subset_for_calib_stats, is_maximisable=True, name="NSE_log_bias")
    nse_log_bias_valid = Objective(nse_log_bias_nan, pd_data['Qobs'], subset_for_valid_stats, is_maximisable=True, name="NSE_log_bias")
    # Due to both BF GR4J and probably spotpy not using/supporting 
    # fully fledged pandas time series, we have to hack things a bit to 
    # have a decent subsetting, and put that information into the Objective object.
    nse_log_bias_obj.simul_start_date = run_start
    nse_log_bias_valid.simul_start_date = run_start
    base_simulation = create_simulation(pd_data, run_start, run_end)
    best_pset = calibrate_lumped_catchment(base_simulation, nse_log_bias_obj, param_space, param_fixed, max_iter=rep)
    mod_runoff = runoff_with_params(base_simulation, best_pset)
    valid_obj = nse_log_bias_valid.objective_statistic(mod_runoff)
    p = dict(best_pset)
    p['NSE_log_bias_valid'] = valid_obj
    p['cat_id'] = cat_id
    return p


# In[13]:


###
# Batch calibration/verification
###

root_path = '/home/jovyan/work/AUS_Catchments/AUS_Catchments_Inputs'

results = []
for cat_id in cat_ids:
    try:
        print "The calibration and validation for catchment %s" % (cat_id)
        p = calib_valid_catchment_id(cat_id)
        results.append(p)
    except:
        print("ERROR: calibration failed for " + cat_id)
        continue #This is not best practice in general...
        
calib_valid_df = pd.DataFrame.from_records(results)
os.chdir('/home/jovyan/work/AUS_Catchments/AUS_Catchments_Outputs/x2/silo_runoff/monthly/nselogbias')
calib_valid_df.to_csv('calib_valid_silo_monthly_x2_nselogbias_new_cats.csv')


# In[14]:


calib_valid_df


# In[15]:


def no_transform(pd_series):
    return pd_series

def observed_modelled_runoff(cat_id, transform_output_series=no_transform):
    pd_data = load_catchment_data(cat_id)
    obs_runoff = pd_data['Qobs'][default_run_start:default_run_end]
    nse_log_bias_obj = Objective(nse_log_bias_nan, pd_data['Qobs'], subset_for_calib_stats, is_maximisable=True, name="NSE_log_bias")
    nse_log_bias_valid = Objective(nse_log_bias_nan, pd_data['Qobs'], subset_for_valid_stats, is_maximisable=True, name="NSE_log_bias")
    nse_log_bias_obj.simul_start_date = default_run_start
    nse_log_bias_valid.simul_start_date = default_run_start
    base_simulation = create_simulation(pd_data, default_run_start, default_run_end)
    best_pset = calibrate_lumped_catchment(base_simulation, nse_log_bias_obj, default_param_space, default_param_fixed, max_iter = default_max_calib_iter)
    mod_runoff = runoff_with_params(base_simulation, best_pset)
    obs_mod_runoff = concat_pandas_series(['mod_runoff','obs_runoff'], transform_output_series(mod_runoff), transform_output_series(obs_runoff))
    return obs_mod_runoff


# In[16]:


###
# Batch runoff timeseries csv generation
###

# new and old cat_ids areas

areas = {'225110':132000000, '225213':312000000,'225219':570000000,'226007': 206000000,'226204':557000000,'226226':299000000,
         '229661':55000000,'403210':1229000000,'403213':231000000,'403217':181000000,'403222':415000000,'403224':154000000,
         '403226':108000000,'403232':124000000,'403244':122000000,'404207':448000000,'404208':91000000,'405205':106000000,
         '603002':464000000, '603003':235000000, '603136':523000000,'604001':1072000000,'604053':1786000000,'605012':4467000000,
         '606001':474000000,'606185':379000000,'606195':246000000,'606218':398000000,'607002':87000000,'607003':2985000000,
         '607004':661000000,'607007':980000000,'607013':249000000,'607144':470000000,'607155':116000000,'607220':4104000000,
         '608001':156000000,'608151':774000000,'608171':65000000,'609002':641000000,'609003':163000000,'609005':83000000,
         '609016':179000000,'609017':550000000,'609018':544000000,'609022':186000000,'610001':435000000,'610009':208000000,
         '610010':417000000,'610015':163000000,'610219':317000000,'611004':818000000,'611006':631000000,'611111':104000000,
         '612001':1331000000,'612023':59000000,'612025':177000000,'612034':672000000,'612230':169000000,'613002':148000000,
         '614123':59000000,'614196':1418000000}

for cat_id in cat_ids:
    try:
        print "The modelled and observed runoff timeseries for catchment %s" % (cat_id)      
        # Daily timeseries
        ts_daily = observed_modelled_runoff(cat_id)
        daily_obs_runoff_m3_per_s = (ts_daily.obs_runoff * areas[cat_id]) / (1000*24*60*60)
        daily_obs_runoff_m3_per_s.name = 'obs_runoff_m3_per_s'
        daily_mod_runoff_m3_per_s = (ts_daily.mod_runoff * areas[cat_id]) / (1000*24*60*60)
        daily_mod_runoff_m3_per_s.name = 'mod_runoff_m3_per_s'
        ts_daily = pd.concat([ts_daily, daily_obs_runoff_m3_per_s, daily_mod_runoff_m3_per_s], axis = 1)
        os.chdir('/home/jovyan/work/AUS_Catchments/AUS_Catchments_Outputs/x2/silo_runoff/monthly/nselogbias/daily_ts_graphs')
        ts_daily.to_csv('%s_runoff_daily_ts_silo_monthly_x2_nselogbias.csv' % (cat_id))
        # Monthly timeseries
        ts_monthly = to_monthly(ts_daily[['mod_runoff','obs_runoff']])
        os.chdir('/home/jovyan/work/AUS_Catchments/AUS_Catchments_Outputs/x2/silo_runoff/monthly/nselogbias/monthly_ts_graphs')
        ts_monthly.to_csv('%s_runoff_monthly_ts_silo_monthly_x2_nselogbias.csv' % (cat_id))
    except:
        print("ERROR: time series generation failed for " + cat_id)
        continue #This is not best practice in general...


# In[17]:


def plot_daily_observed_modelled_runoff(cat_id, transform_output_series=no_transform):
    pd_data = load_catchment_data(cat_id)
    obs_runoff = pd_data['Qobs'][default_run_start:default_run_end]
    nse_log_bias_obj = Objective(nse_log_bias_nan, pd_data['Qobs'], subset_for_calib_stats, is_maximisable=True, name="NSE_log_bias")
    nse_log_bias_valid = Objective(nse_log_bias_nan, pd_data['Qobs'], subset_for_valid_stats, is_maximisable=True, name="NSE_log_bias")
    nse_log_bias_obj.simul_start_date = default_run_start
    nse_log_bias_valid.simul_start_date = default_run_start
    base_simulation = create_simulation(pd_data, default_run_start, default_run_end)
    best_pset = calibrate_lumped_catchment(base_simulation, nse_log_bias_obj, default_param_space, default_param_fixed, max_iter = default_max_calib_iter)
    mod_runoff = runoff_with_params(base_simulation, best_pset)
    obs_mod_runoff = concat_pandas_series(['mod_runoff','obs_runoff'], transform_output_series(mod_runoff), transform_output_series(obs_runoff))
    graph = obs_mod_runoff.plot(title = 'Graph of monthly observed and modelled runoff in mm/day').set(xlabel ='year', ylabel = 'runoff (mm/day))')
    return graph


# In[18]:


def plot_monthly_observed_modelled_runoff(cat_id, transform_output_series=to_monthly):
    pd_data = load_catchment_data(cat_id)
    obs_runoff = pd_data['Qobs'][default_run_start:default_run_end]
    nse_log_bias_obj = Objective(nse_log_bias_nan, pd_data['Qobs'], subset_for_calib_stats, is_maximisable=True, name="NSE_log_bias")
    nse_log_bias_valid = Objective(nse_log_bias_nan, pd_data['Qobs'], subset_for_valid_stats, is_maximisable=True, name="NSE_log_bias")
    nse_log_bias_obj.simul_start_date = default_run_start
    nse_log_bias_valid.simul_start_date = default_run_start
    base_simulation = create_simulation(pd_data, default_run_start, default_run_end)
    best_pset = calibrate_lumped_catchment(base_simulation, nse_log_bias_obj, default_param_space, default_param_fixed, max_iter = default_max_calib_iter)
    mod_runoff = runoff_with_params(base_simulation, best_pset)
    obs_mod_runoff = concat_pandas_series(['mod_runoff','obs_runoff'], transform_output_series(mod_runoff), transform_output_series(obs_runoff))
    graph = obs_mod_runoff.plot(title = 'Graph of monthly observed and modelled runoff in mm/month').set(xlabel ='year', ylabel = 'runoff (mm/month)')
    return graph


# In[19]:


###
# Batch runoff graphs generation
###

for cat_id in cat_ids:
    try:
        print "The modelled and observed runoff graphs for catchment %s" % (cat_id)
        daily_graph = plot_daily_observed_modelled_runoff(cat_id)
        os.chdir('/home/jovyan/work/AUS_Catchments/AUS_Catchments_Outputs/x2/silo_runoff/monthly/nselogbias/daily_ts_graphs')
        plt.savefig('%s_runoff_daily_graph_silo_monthly_x2_nselogbias.png' % (cat_id))
        plt.show()
        monthly_graph = plot_monthly_observed_modelled_runoff(cat_id)
        os.chdir('/home/jovyan/work/AUS_Catchments/AUS_Catchments_Outputs/x2/silo_runoff/monthly/nselogbias/monthly_ts_graphs')
        plt.savefig('%s_runoff_monthly_graph_silo_monthly_x2_nselogbias.png' % (cat_id))
        plt.show()
    except:
        print("ERROR: time series generation failed for " + cat_id)
        continue #This is not best practice in general...


# In[ ]:




