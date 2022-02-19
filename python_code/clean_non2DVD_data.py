# This program process disdrometer count (non 2DVD) raw data and prepare them for statistical analysis
# The following "rounds" of filtering/processing are implemented:
#   FIRST ROUND (remove zero counts records)
#   SECOND ROUND (remove quiescent minutes)
#   THIRD ROUND (remove outliers of pdf)
#   FOURTH ROUND (remove narrow pdf)
#
# Output 1: clean disdrometer data
# Output 2: summary of the cleaning procedure
#
#
# Author: Massimiliano Ignaccolo
#         massimiliano.ignaccolo [at] tutanota.com. massimiliano.ignaccolo [at] sas.com
# Last modified: Feb 01 2022


# ------------- Necessary Python packages -START
import os
import sys
import pandas as pd
import numpy as np
import argparse

# The disdrorain package is expected to be in the same directory where program is executed
# If this is not the case change value of pack_path accordingly
# pack_path = os.getcwd()
# if pack_path not in sys.path:
# sys.path.append(pack_path)
import disdrorain as dr
# ------------- Necessary Python packages -END

# ARGUMENTS
parser = argparse.ArgumentParser(description='clean impact disdrometer data', epilog="")
# --disdrodata: the disdrometer data to process
parser.add_argument('--disdrodata', action="store", dest='_dd_', default='_NONE_')
# --disdrolimits: the file with the disdrometer class limits (unless standard RD80 are used)
parser.add_argument('--disdrolimits', action="store", dest='_dl_', default='standard')
# --disdroarea: the catchment area of the disdrometer in mm2 (5000 mm2 is default: RD80))
parser.add_argument('--disdroarea', action="store", dest='_da_', default='5000')
# --disdrotimeresolution: the time resolution of the data in seconds (1 min time resolution is default))
parser.add_argument('--disdrotimeresolution', action="store", dest='_dtr_', default='60')
# --disdrooout: the processed didrometer data
parser.add_argument('--disdroout', action="store", dest='_do_', default='_NONE_')
args = parser.parse_args()

print("Executing: ", sys.argv[0])

print('PROCESSING: ', args._dd_)
print()
PI = 3.141592653589793  # approximate value of greek pi
seconds_in_hour = 3600  # for converting rainfall rate in mm/h

if (args._dl_ == "standard"):  # standard
    disdrodata = dr.disdrorain(datapath=args._dd_, instrument_area=float(args._da_), time_interval=float(args._dtr_))
if (args._dl_ != "standard"):  # non 2dvd disdrometer data with specific classlimit
    disdrodata = dr.disdrorain(classpath=args._dl_, datapath=args._dd_, instrument_area=float(args._da_), time_interval=float(args._dtr_))

# ---- preprocess
# define results dataframe
results = pd.DataFrame(columns=['round', 'action_type', 'var', 'value', 'percent_of_total'])
results = results.astype(dtype={'round': 'int64', 'action_type': 'object', 'var':  'object',
                                'value': 'float64', 'percent_of_total': 'float64'})

# ZERO ROUND (Initial Stats)
# define count / rainfall rate dataframe
NR = pd.DataFrame(columns=['N', 'R'])
NR = NR.astype(dtype={'N': 'int64', 'R': 'float64'})
# number of drops trough the catchment area N
NR['N'] = disdrodata.data.sum(axis=1)
# Rainfall rate depends on the 3rd "flux" moment no matter what is the equation for v(D)
R_df_ = disdrodata.flux_moment_calculator([3])
R_df_[np.isnan(R_df_)] = 0
NR['R'] = (PI / 6) * (1 / (disdrodata.instrument_area)) * NR.N * R_df_.M3

new_row = {'round': [0], 'action_type': ["initial stats"], 'var': ["n_records"], 'value': [NR.shape[0]], 'percent_of_total': [100]}
results = pd.concat([results, pd.DataFrame.from_dict(new_row)], ignore_index=True)
new_row = {'round': [0], 'action_type': ["initial stats"], 'var': ["n_drops"], 'value': [NR.N.sum()], 'percent_of_total': [100]}
results = pd.concat([results, pd.DataFrame.from_dict(new_row)], ignore_index=True)
new_row = {'round': [0], 'action_type': ["initial stats"], 'var': ["cumulated rainfall rate"],
           'value': [round(NR.R.sum(), 1)], 'percent_of_total': [100]}
results = pd.concat([results, pd.DataFrame.from_dict(new_row)], ignore_index=True)


# FIRST ROUND (remove zero counts records)
new_row = {'round': [1], 'action_type': ["remove zero count records"], 'var': ["n_records"], 'value': [NR[NR.N > 0].N.count()],
           'percent_of_total': [round(NR[NR.N > 0].N.count()/NR.shape[0]*100, 1)]}
results = pd.concat([results, pd.DataFrame.from_dict(new_row)], ignore_index=True)
new_row = {'round': [1], 'action_type': ["remove zero count records"], 'var': ["n_drops"], 'value': [NR.N.sum()], 'percent_of_total': [100]}
results = pd.concat([results, pd.DataFrame.from_dict(new_row)], ignore_index=True)
new_row = {'round': [1], 'action_type': ["remove zero count records"], 'var': ["cumulated rainfall rate"],
           'value': [round(NR.R.sum(), 1)], 'percent_of_total': [100]}
results = pd.concat([results, pd.DataFrame.from_dict(new_row)], ignore_index=True)

# SECOND ROUND (remove quiescent minutes)
disdrodata2 = disdrodata.remove_counts_below_threshold()
#
NR2 = pd.DataFrame(columns=['N', 'R'])
NR2 = NR2.astype(dtype={'N': 'int64', 'R': 'float64'})
# number of drops trough the catchment area N
NR2['N'] = disdrodata2.data.sum(axis=1)
# Rainfall rate depends on the 3rd "flux" moment no matter what is the equation for v(D)
R_df_ = disdrodata2.flux_moment_calculator([3])
R_df_[np.isnan(R_df_)] = 0
NR2['R'] = (PI / 6) * (1 / (disdrodata2.instrument_area)) * NR2.N * R_df_.M3

new_row = {'round': [2], 'action_type': ["remove quiescent minutes"], 'var': ["n_records"], 'value': [NR2.shape[0]],
           'percent_of_total': [round(NR2.shape[0]/NR.shape[0]*100, 1)]}
results = pd.concat([results, pd.DataFrame.from_dict(new_row)], ignore_index=True)
new_row = {'round': [2], 'action_type': ["remove quiescent minutes"], 'var': ["n_drops"], 'value': [NR2.N.sum()],
           'percent_of_total': [round(NR2.N.sum()/NR.N.sum()*100, 1)]}
results = pd.concat([results, pd.DataFrame.from_dict(new_row)], ignore_index=True)
new_row = {'round': [2], 'action_type': ["remove quiescent minutes"], 'var': ["cumulated rainfall rate"], 'value': [round(NR2.R.sum(), 1)],
           'percent_of_total': [round(NR2.R.sum()/NR.R.sum()*100, 1)]}
results = pd.concat([results, pd.DataFrame.from_dict(new_row)], ignore_index=True)

# THIRD ROUND (remove outliers of pdf)
(disdrodata3, summary) = disdrodata2.outlier_deletion()
#
NR3 = pd.DataFrame(columns=['N', 'R'])
NR3 = NR3.astype(dtype={'N': 'int64', 'R': 'float64'})
# number of drops trough the catchment area N
NR3['N'] = disdrodata3.data.sum(axis=1)
# Rainfall rate depends on the 3rd "flux" moment no matter what is the equation for v(D)
R_df_ = disdrodata3.flux_moment_calculator([3])
R_df_[np.isnan(R_df_)] = 0
NR3['R'] = (PI / 6) * (1 / (disdrodata3.instrument_area)) * NR3.N * R_df_.M3

new_row = {'round': [3], 'action_type': ["remove outliers of pdf"], 'var': ["n_records"], 'value': [NR3.shape[0]],
           'percent_of_total': [round(NR3.shape[0]/NR.shape[0]*100, 1)]}
results = pd.concat([results, pd.DataFrame.from_dict(new_row)], ignore_index=True)
new_row = {'round': [3], 'action_type': ["remove outliers of pdf"], 'var': ["n_drops"], 'value': [NR3.N.sum()],
           'percent_of_total': [round(NR3.N.sum()/NR.N.sum()*100, 1)]}
results = pd.concat([results, pd.DataFrame.from_dict(new_row)], ignore_index=True)
new_row = {'round': [3], 'action_type': ["remove outliers of pdf"], 'var': ["cumulated rainfall rate"], 'value': [round(NR3.R.sum(), 1)],
           'percent_of_total': [round(NR3.R.sum()/NR.R.sum()*100, 1)]}
results = pd.concat([results, pd.DataFrame.from_dict(new_row)], ignore_index=True)

# FOURTH ROUND (remove narrow pdf)
(disdrodata4, disdrodata3_n, summary_nn) = disdrodata3.remove_narrow(_nclmin_=3)

NR4 = pd.DataFrame(columns=['N', 'R'])
NR4 = NR4.astype(dtype={'N': 'int64', 'R': 'float64'})
# number of drops trough the catchment area N
NR4['N'] = disdrodata4.data.sum(axis=1)
# Rainfall rate depends on the 3rd "flux" moment no matter what is the equation for v(D)
R_df_ = disdrodata4.flux_moment_calculator([3])
R_df_[np.isnan(R_df_)] = 0
NR4['R'] = (PI / 6) * (1 / int((disdrodata4.instrument_area))) * NR4.N * R_df_.M3

new_row = {'round': [4], 'action_type': ["remove narrow pdf"], 'var': ["n_records"], 'value': [NR4.shape[0]],
           'percent_of_total': [round(NR4.shape[0]/NR.shape[0]*100, 1)]}
results = pd.concat([results, pd.DataFrame.from_dict(new_row)], ignore_index=True)
new_row = {'round': [4], 'action_type': ["remove narrow pdf"], 'var': ["n_drops"], 'value': [NR4.N.sum()],
           'percent_of_total': [round(NR4.N.sum()/NR.N.sum()*100, 1)]}
results = pd.concat([results, pd.DataFrame.from_dict(new_row)], ignore_index=True)
new_row = {'round': [4], 'action_type': ["remove narrow pdf"], 'var': ["cumulated rainfall rate"], 'value': [round(NR4.R.sum(), 1)],
           'percent_of_total': [round(NR4.R.sum()/NR.R.sum()*100, 1)]}
results = pd.concat([results, pd.DataFrame.from_dict(new_row)], ignore_index=True)


disdrodata4.data.to_csv(args._do_, sep=' ', header=None, index=None)
results.to_csv('summary_'+args._do_, sep=':', index=None)
