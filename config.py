# Meta data with date, collar name, start and end times for each measurement, system height, optionally air pressure and temperature if they are missing from chamber_data-file.
meta_data_path = "Example_files/Metadata.csv"

# Path to chamber_data-files
data_path = 'Example_files/'

# Seconds to remove from beginning of each measurement
remove_from_start = 15

# Minimum length of accetable measurement in seconds
min_length = 80

# Use exponential or linear fits for calculating fluxes, "exponential"/"linear"
fit_type = "linear"

# Quality control

# If PAR standard deviation is greater than this, the measurement is discarded and the flux is not calculated, FEEL FREE TO CHANGE THE LIMIT
par_sd_limit = 100

# If absolute value of flux > 0.1 and if NRMSE > nrmse_limit, the measurement is discarded and the flux is not calculated, FEEL FREE TO CHANGE THE LIMIT
nrmse_limit = 0.1

# Plotting
plotting = True  # If False, no plots are made

# Path to data, data must include date, collar ID, flux and PAR
flux_data_path = 'Results/flux_results.csv'

# Bounds for parameters: min, max
alpha_bounds = -0.1, -0.00000001
GPmax_bounds = -5, -0.0000001
