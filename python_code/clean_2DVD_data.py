# This program process 2DVD  disdrometer raw data and prepare them for statistical analysis
# The following "rounds" of filtering/processing are implemented:
#   FIRST ROUND (remove drops with off-bound speed)
#   SECOND ROUND (remove quiescent minutes)
#   THIRD ROUND (remove outliers of pdf)
#   FOURTH ROUND (remove narrow pdf)
#
# Output 1: clean disdrometer data
# Output 2: summary of the cleaning procedure
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


# to remove drops ouliers we cast 2DVD dat to PARSIVEL classes
def remove_drops_outlier(r2dvd, rparsi_before, rparsi_after, _summary_):
    _t_ = r2dvd.data.groupby('timestamp').count()
    _t_.reset_index(inplace=True)
    _t_.reset_index(inplace=True)
    _t_.drop(columns=['diameter', 'speed'], inplace=True)
    _t_.set_index('timestamp', inplace=True)
    r2dvd.data.set_index('timestamp', inplace=True)
    _y_ = r2dvd.data.join(_t_)

    for j in range(0, _summary_.shape[0]):
        record = _summary_.iloc[j, 0]
        # print(j,record,type(record))
        rb = rparsi_before.data.iloc[record, 0: 32]
        ra = rparsi_after.data.iloc[record, 0: 32]
        for l in range(0, 32):
            if rb[l] != ra[l]:
                a = rparsi_before.classlimits.loc['left', rb.index[l]]
                b = rparsi_before.classlimits.loc['right', rb.index[l]]
                _y_ = _y_.loc[(_y_['index'] != record) | ((_y_['index'] == record) & ((_y_.diameter < a) | (_y_.diameter >= b))), :]
    return(_y_)


# to remove too narrow pdfs we cast 2DVD dat to PARSIVEL classes
def remove_drops_narrow(r2dvd, rparsi_narrows, _summary_):
    _t_ = r2dvd.data.groupby('timestamp').count()
    _t_.reset_index(inplace=True)
    _t_.reset_index(inplace=True)
    _t_.drop(columns=['diameter', 'speed'], inplace=True)
    _t_.set_index('timestamp', inplace=True)
    recs_toremove = rparsi_narrows.data.index.values.tolist()
    _t_ = _t_.loc[_t_['index'].isin(recs_toremove), :]
    timestamps_toremove = _t_.index.values.tolist()
    _y_ = r2dvd.data.loc[~r2dvd.data['timestamp'].isin(timestamps_toremove), :]
    return(_y_)


# ARGUMENTS
parser = argparse.ArgumentParser(description='clean impact disdrometer data', epilog="")
# --disdrodata: the disdrometer data to process
parser.add_argument('--disdrodata', action="store", dest='_dd_', default='_NONE_')
# --disdrolimits: the file with the PARSIVEL disdrometer class limits
# 2DVD data are cast to PARSIVEL drop count for some rounds
parser.add_argument('--disdroparsivellimits', action="store", dest='_dl_', default='_NONE_')
# --disdroarea: the catchment area of the disdrometer in mm2 (default = 10000 typical area of 2DVD disdro))
parser.add_argument('--disdroarea', action="store", dest='_da_', default='10000')
# --disdrotimeresolution: the time resolution of the data in seconds
parser.add_argument('--disdrotimeresolution', action="store", dest='_dtr_', default='60')
# --disdrooout: the processed didrometer data
parser.add_argument('--disdroout', action="store", dest='_do_', default='_NONE_')
args = parser.parse_args()

#
print('PROCESSING: ', args._dd_)
print()
PI = 3.141592653589793  # approximate value of greek pi
seconds_in_hour = 3600  # for converting rainfall rate in mm/h

disdrodata = dr.disdrorain_2dvd(datapath=args._dd_, instrument_area=float(args._da_), time_interval=float(args._dtr_))

# print(disdrodata.data.head(5))
# ---- preprocess
# define results dataframe
results = pd.DataFrame(columns=['round', 'action_type', 'var', 'value', 'percent_of_total'])
results = results.astype(dtype={'round': 'int64', 'action_type': 'object', 'var':  'object',
                                'value': 'float64', 'percent_of_total': 'float64'})

# ZERO ROUND (Initial Stats)
# define count / rainfall rate dataframe
print("ROUND: 0")
NR = pd.DataFrame(columns=['N', 'R'])
NR = NR.astype(dtype={'N': 'int64', 'R': 'float64'})
# number of drops trough the catchment area N
NR['N'] = disdrodata.data.groupby('timestamp')['diameter'].count().values
# Rainfall rate depends on the 3rd "flux" moment no matter what is the equation for v(D)
R_df_ = disdrodata.flux_moment_calculator([3])
R_df_[np.isnan(R_df_)] = 0
NR['R'] = (PI / 6) * (1 / (disdrodata.instrument_area)) * NR.N * R_df_.M3

new_row = {'round': 0, 'action_type': "initial stats", 'var': "n_records", 'value': NR.shape[0], 'percent_of_total': 100}
results = results.append(new_row, ignore_index=True)
new_row = {'round': 0, 'action_type': "initial stats", 'var': "n_drops", 'value': NR.N.sum(), 'percent_of_total': 100}
results = results.append(new_row, ignore_index=True)
new_row = {'round': 0, 'action_type': "initial stats", 'var': "cumulated rainfall rate",
           'value': round(NR.R.sum(), 1), 'percent_of_total': 100}
results = results.append(new_row, ignore_index=True)

# FIRST ROUND (remove drops with off-bound speed: see krawjeski paper)
print("ROUND: 1")
disdrodata1 = disdrodata.remove_drops_offbound_speed()

# define count / rainfall rate dataframe
NR1 = pd.DataFrame(columns=['N', 'R'])
NR1 = NR1.astype(dtype={'N': 'int64', 'R': 'float64'})
# number of drops trough the catchment area N
NR1['N'] = disdrodata1.data.groupby('timestamp')['diameter'].count().values
# Rainfall rate depends on the 3rd "flux" moment no matter what is the equation for v(D)
R_df_ = disdrodata1.flux_moment_calculator([3])
R_df_[np.isnan(R_df_)] = 0
NR1['R'] = (PI / 6) * (1 / (disdrodata.instrument_area)) * NR1.N * R_df_.M3

new_row = {'round': 1, 'action_type': "remove drops with off-bound speed", 'var': "n_records", 'value': NR1.N.count(),
           'percent_of_total': round(NR1.N.count()/NR.shape[0]*100, 1)}
results = results.append(new_row, ignore_index=True)
new_row = {'round': 1, 'action_type': "remove drops with off-bound speed", 'var': "n_drops", 'value': NR1.N.sum(),
           'percent_of_total': round(NR1.N.sum()/NR.N.sum()*100, 1)}
results = results.append(new_row, ignore_index=True)
new_row = {'round': 1, 'action_type': "remove drops with off-bound speed", 'var': "cumulated rainfall rate",
           'value': round(NR1.R.sum(), 1), 'percent_of_total': round(NR1.R.sum()/NR.R.sum()*100, 1)}
results = results.append(new_row, ignore_index=True)

# SECOND ROUND (remove quiescent minutes)
print("ROUND: 2")
disdrodata2 = disdrodata1.remove_counts_below_threshold()

# define count / rainfall rate dataframe
NR2 = pd.DataFrame(columns=['N', 'R'])
NR2 = NR2.astype(dtype={'N': 'int64', 'R': 'float64'})
# number of drops trough the catchment area N
NR2['N'] = disdrodata2.data.groupby('timestamp')['diameter'].count().values
# Rainfall rate depends on the 3rd "flux" moment no matter what is the equation for v(D)
R_df_ = disdrodata2.flux_moment_calculator([3])
R_df_[np.isnan(R_df_)] = 0
NR2['R'] = (PI / 6) * (1 / (disdrodata2.instrument_area)) * NR2.N * R_df_.M3

new_row = {'round': 2, 'action_type': "remove quiescent minutes", 'var': "n_records", 'value': NR2.shape[0],
           'percent_of_total': round(NR2.shape[0]/NR.shape[0]*100, 1)}
results = results.append(new_row, ignore_index=True)
new_row = {'round': 2, 'action_type': "remove quiescent minutes", 'var': "n_drops", 'value': NR2.N.sum(),
           'percent_of_total': round(NR2.N.sum()/NR.N.sum()*100, 1)}
results = results.append(new_row, ignore_index=True)
new_row = {'round': 2, 'action_type': "remove quiescent minutes", 'var': "cumulated rainfall rate",
           'value': round(NR2.R.sum(), 1), 'percent_of_total': round(NR2.R.sum()/NR.R.sum()*100, 1)}
results = results.append(new_row, ignore_index=True)

# THIRD ROUND (remove outliers of pdf)
print("ROUND: 3")
# cast 2DVD data to parsivel classes
_df_ = disdrodata2.twodvd_to_parsivel()
disdrodata2a = dr.disdrorain(classpath=args._dl_, dataframe=_df_, instrument_area=10000)
(disdrodata2b, summary) = disdrodata2a.outlier_deletion()
_mydata_ = remove_drops_outlier(disdrodata2, disdrodata2a, disdrodata2b, summary)
_mydata_.reset_index(inplace=True)
_mydata_.drop(columns=['index'], inplace=True)
_mydata_.index.rename('record number', inplace=True)
disdrodata3 = dr.disdrorain_2dvd(dataframe=_mydata_, instrument_area=10000)

#
NR3 = pd.DataFrame(columns=['N', 'R'])
NR3 = NR3.astype(dtype={'N': 'int64', 'R': 'float64'})
# number of drops trough the catchment area N
NR3['N'] = disdrodata3.data.groupby('timestamp')['diameter'].count().values
# Rainfall rate depends on the 3rd "flux" moment no matter what is the equation for v(D)
R_df_ = disdrodata3.flux_moment_calculator([3])
R_df_[np.isnan(R_df_)] = 0
NR3['R'] = (PI / 6) * (1 / (disdrodata3.instrument_area)) * NR3.N * R_df_.M3

new_row = {'round': 3, 'action_type': "remove outliers of pdf", 'var': "n_records", 'value': NR3.shape[0],
           'percent_of_total': round(NR3.shape[0]/NR.shape[0]*100, 1)}
results = results.append(new_row, ignore_index=True)
new_row = {'round': 3, 'action_type': "remove outliers of pdf", 'var': "n_drops", 'value': NR3.N.sum(),
           'percent_of_total': round(NR3.N.sum()/NR.N.sum()*100, 1)}
results = results.append(new_row, ignore_index=True)
new_row = {'round': 3, 'action_type': "remove outliers of pdf", 'var': "cumulated rainfall rate",
           'value': round(NR3.R.sum(), 1), 'percent_of_total': round(NR3.R.sum()/NR.R.sum()*100, 1)}
results = results.append(new_row, ignore_index=True)

# FOURTH ROUND (remove narrow pdf)
print("ROUND: 4")
# cast 2DVD data to parsivel classes
_df_ = disdrodata3.twodvd_to_parsivel()
disdrodata3a = dr.disdrorain(classpath=args._dl_, dataframe=_df_, instrument_area=10000)
(disdrodata3b, disdrodata3a_n, summary) = disdrodata3a.remove_narrow(_nclmin_=3)
_mydata_ = remove_drops_narrow(disdrodata3, disdrodata3a_n, summary)
disdrodata4 = dr.disdrorain_2dvd(dataframe=_mydata_, instrument_area=10000)

#
NR4 = pd.DataFrame(columns=['N', 'R'])
NR4 = NR4.astype(dtype={'N': 'int64', 'R': 'float64'})
# number of drops trough the catchment area N
NR4['N'] = disdrodata4.data.groupby('timestamp')['diameter'].count().values
# Rainfall rate depends on the 3rd "flux" moment no matter what is the equation for v(D)
R_df_ = disdrodata4.flux_moment_calculator([3])
R_df_[np.isnan(R_df_)] = 0
NR4['R'] = (PI / 6) * (1 / int((disdrodata4.instrument_area))) * NR4.N * R_df_.M3

new_row = {'round': 4, 'action_type': "remove narrow pdf", 'var': "n_records", 'value': NR4.shape[0],
           'percent_of_total': round(NR4.shape[0]/NR.shape[0]*100, 1)}
results = results.append(new_row, ignore_index=True)
new_row = {'round': 4, 'action_type': "remove narrow pdf", 'var': "n_drops", 'value': NR4.N.sum(),
           'percent_of_total': round(NR4.N.sum()/NR.N.sum()*100, 1)}
results = results.append(new_row, ignore_index=True)
new_row = {'round': 4, 'action_type': "remove narrow pdf", 'var': "cumulated rainfall rate",
           'value': round(NR4.R.sum(), 1), 'percent_of_total': round(NR4.R.sum()/NR.R.sum()*100, 1)}
results = results.append(new_row, ignore_index=True)

# save result to file
disdrodata4.data.to_csv(args._do_, sep=' ', header=None, index=None)
results.to_csv('summary_'+args._do_, sep=':', index=None)
