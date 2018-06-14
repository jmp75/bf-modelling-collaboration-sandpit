

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
obs_mod_runoff['cat_id'] = cat_id
obs_mod_runoff['month'] = np.arange(1, 13)

return obs_mod_runoff

graph = obs_mod_runoff.plot(title = 'Graph of monthly observed and modelled runoff in mm/day').set(xlabel ='year', ylabel = 'runoff (mm/day))')
return graph

