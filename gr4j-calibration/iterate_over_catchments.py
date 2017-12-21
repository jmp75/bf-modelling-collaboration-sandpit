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

csv_file = os.path.join(root_path,'225110.csv')

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
run_end = datetime(1999, 12, 31)

cal_inputs = {
    'rainfall' : rainfall_data[run_start:run_end].as_matrix(),
    'pet'  : pet_data[run_start:run_end].as_matrix(),
    # 'obs'  : obs_runoff_data[run_start:run_end].as_matrix(),
    'x1':20.0,
    'x2':0.0,
    'x3':40.0,
    'x4':0.5
}

runoff = gr4j.gr4j(**cal_inputs)

library(lubridate)
ymd('1975-01-01')
