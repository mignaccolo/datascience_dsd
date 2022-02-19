# This program calculate the accuracy of Local Adaptive Fit Algorithm for a particular site/dataset
#
# Output 1: estimated probability and cumulated probability for error metric associated with LAF
#
# Author: Massimiliano Ignaccolo
#         massimiliano.ignaccolo [at] tutanota.com. massimiliano.ignaccolo [at] sas.com
# Last modified: Feb 01 2022

# ------------- Necessary Python packages -START
import sys
import pandas as pd
import numpy as np
import argparse
import csv
import copy
from grispy import GriSPy

# The disdrorain package is expected to be in the same directory where program is executed
# If this is not the case change value of pack_path accordingly
# pack_path = os.getcwd()
# if pack_path not in sys.path:
# sys.path.append(pack_path)
import disdrorain as dr
#  ------------- Necessary Python packages -END

# ARGUMENTS
parser = argparse.ArgumentParser(description='', epilog="")
# --input: path to file with the adaptive fit results output
parser.add_argument('--input', action="store", dest='_if_', default='_NONE_')
# --output: output file name for the fit accuracy
parser.add_argument('--output', action="store", dest='_of_', default='_NONE_')
# --pdfparameter: the statistical moments for which we want ot calculate the insutu
# LAF fitting accuracy
parser.add_argument('--pdfparameter', action="store", dest='_pp_', default='_NONE_')
args = parser.parse_args()

#
print("Executing: ", sys.argv[0])
print()

# load input file to data frame
_mydf_ = pd.read_csv(args._if_, sep=' ')
# make hard copy
_mydf2_ = copy.deepcopy(_mydf_)
# remove un necessary info
_mydf2_.drop(columns=['mu_r', 'gamma_r', 'predicted_q1', 'predicted_q3', 'predicted_5%',
                      'predicted_95%', 'predicted_r'], inplace=True)
# we first calculate semi IQR and semi 595IQR
_mydf2_.loc[:, 'semiiqr'] = (_mydf2_[f"{args._pp_}_q3"]-_mydf2_[f"{args._pp_}_q1"])/2
_mydf2_.loc[:, 'semi595'] = (_mydf2_[f"{args._pp_}_95%"]-_mydf2_[f"{args._pp_}_5%"])/2
# and then the corresponding coefficient of variations
_mydf2_.loc[:, 'coeffvar'] = np.abs(_mydf2_['semiiqr']/_mydf2_[f"{args._pp_}_median"])
_mydf2_.loc[:, 'coeffvar595'] = np.abs(_mydf2_['semi595']/_mydf2_[f"{args._pp_}_median"])

# We use a bin of 0.05 (5%) to calculate the pdf of the
# coefficient of variation relative to the semi iqr
_mydf2_.loc[:, 'coeffvarbin'] = np.floor(_mydf2_.coeffvar/0.05).astype(int)
_mydf2_.loc[:, 'coeffvar595bin'] = np.floor(_mydf2_.coeffvar595/0.05).astype(int)
_stat_ = _mydf2_.groupby(['coeffvarbin', 'site']).agg({'radius': 'mean', 'mu': 'count'})
_stat_.reset_index(inplace=True)
_stat_.loc[:, 'value'] = _stat_.coeffvarbin*0.05
_stat_.loc[:, 'prob'] = _stat_.mu/_mydf_.shape[0]
_stat_.loc[:, 'variable'] = 'coeffvar'
_stat_.loc[:, 'analysis_type'] = 'pdf'

# We use a bin of 0.01 (1%) to calculate the cdf of the
# coefficient of variation relative to the semi iqr
_mydf2_.loc[:, 'coeffvarbin1'] = np.floor(_mydf2_.coeffvar/0.01).astype(int)
_mydf2_.loc[:, 'coeffvar595bin1'] = np.floor(_mydf2_.coeffvar595/0.01).astype(int)
_stat_cdf = _mydf2_.groupby(['coeffvarbin1', 'site']).agg({'radius': 'mean', 'mu': 'count'})
_stat_cdf.reset_index(inplace=True)
_stat_cdf.loc[:, 'value'] = _stat_cdf.coeffvarbin1*0.01
_stat_cdf.loc[:, 'prob'] = _stat_cdf.mu/_mydf_.shape[0]
_stat_cdf.loc[:, 'variable'] = 'coeffvar'
_stat_cdf.loc[:, 'cumprob'] = _stat_cdf.loc[:, 'prob'].cumsum()
_stat_cdf.loc[:, 'analysis_type'] = 'cdf'

# We use a bin of 0.05 (5%) to calculate the pdf of the
# coefficient of variation relative to the semi 595IQR
_stat595_ = _mydf2_.groupby(['coeffvar595bin', 'site']).agg({'radius': 'mean', 'mu': 'count'})
_stat595_.reset_index(inplace=True)
_stat595_.loc[:, 'value'] = _stat595_.coeffvar595bin*0.05
_stat595_.loc[:, 'prob'] = _stat595_.mu/_mydf_.shape[0]
_stat595_.loc[:, 'variable'] = 'coeffvar595'
_stat595_.loc[:, 'analysis_type'] = 'pdf'

# We use a bin of 0.01 (1%) to calculate the cdf of the
# coefficient of variation relative to the semi 595IQR
_stat595_cdf = _mydf2_.groupby(['coeffvar595bin1', 'site']).agg({'radius': 'mean', 'mu': 'count'})
_stat595_cdf.reset_index(inplace=True)
_stat595_cdf.loc[:, 'value'] = _stat595_cdf.coeffvar595bin1*0.01
_stat595_cdf.loc[:, 'prob'] = _stat595_cdf.mu/_mydf_.shape[0]
_stat595_cdf.loc[:, 'variable'] = 'coeffvar595'
_stat595_cdf.loc[:, 'cumprob'] = _stat595_cdf.loc[:, 'prob'].cumsum()
_stat595_cdf.loc[:, 'analysis_type'] = 'cdf'

# concate all resutls in one data frame
_stat_total = pd.concat([_stat_, _stat595_], ignore_index=True)
_stat_total = pd.concat([_stat_total, _stat_cdf], ignore_index=True)
_stat_total = pd.concat([_stat_total, _stat595_cdf], ignore_index=True)

# keep relevant columns
_stat_total = _stat_total[['analysis_type', 'variable', 'value', 'prob', 'cumprob', 'radius', 'site']]
_stat_total.rename(columns={'radius': 'avg_radius'}, inplace=True)

# save to file
# cumprob values are missing when the analysis type is pdf:
#   we indicate this values with "NaN" in the output file
_stat_total.to_csv(args._of_, sep=' ', na_rep='NaN', index=None)
