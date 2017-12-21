# Use case
# Import CSV data

root_path = r'C:\Users\per202\Documents\BF\AUS_Catchments\AUS_Catchments'

import pandas as pd
import numpy as np
import spotpy
import sys

if(sys.platform == 'win32'):
    sys.path.append(r'C:\src\github_jm\gr4j-sg')
else:
    sys.path.append('/home/per202/src/csiro/stash/per202/bf/gr4j-sg')

import gr4j

csv_file = os.path.join(root_path,'225110.csv')

tseries_df = pd.read_csv(csv_file)
# TODO
# stopifnot expected colnames 

start_year = tseries_df['year'][0]
start_month = tseries_df['month'][0]
start_day = tseries_df['date'][0]  # date ?? anyway...

tseries_len = tseries_df.shape[0] 
rng = pd.date_range(str(start_day) + '/' + str(start_month) + '/' + str(start_year), periods=tseries_len, freq='D')

# "Fun" facts: a pd.Series object cannot be at argument to pd.Series, it seems
rainfall_data = pd.Series(tseries_df['Rain'].as_matrix(), index=rng)

# Sanity check: does time stamps match what I see via excel? Nope!
# This will need to be reconciled before doing real runs with it.
rainfall_data.head()
rainfall_data.tail()
# 2012-01-03    3.748100
# 2012-01-04    0.056487
# 2012-01-05    0.068906
# 2012-01-06    1.125100
# 2012-01-07    8.436600
# But the tail end via excel is:
# 2012	6	29	1.1251	0.98921	-9999
# 2012	6	30	8.4366	0.94867	-9999

pet_data = pd.Series(tseries_df['Etp'].as_matrix(), index=rng)
obs_runoff_data = pd.Series(tseries_df['Qobs'].as_matrix(), index=rng)

obs_runoff_data[obs_runoff_data < 0] = np.nan

# Date range for calibration. 

from datetime import datetime

run_start = datetime(1980, 1, 1)
run_end = datetime(1999, 12, 31)

cal_inputs = {}
cal_inputs['rain'] = rainfall_data[run_start:run_end].as_matrix()
cal_inputs['pet'] = pet_data[run_start:run_end].as_matrix()
cal_inputs['obs'] = obs_runoff_data[run_start:run_end].as_matrix()



