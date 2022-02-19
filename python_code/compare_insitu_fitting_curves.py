# This program compare (calculate the discrepancy) the LAF fitted curves between two different sites/datasets
#
# Output 1: estimated discrepancy between the LAF fits of the different datasets (single row format)
# Output 2: estimated discrepancy between the LAF fits of the different datasets (double row format)
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
import copy as cp
import weightedstats as ws

# The disdrorain package is expected to be in the same directory where programy is executed
# If this is not the case change value of pack_path accordingly
# pack_path = os.getcwd()
# if pack_path not in sys.path:
# sys.path.append(pack_path)
import disdrorain as dr
#  ------------- Necessary Python packages -END

# ARGUMENTS
parser = argparse.ArgumentParser(description='', epilog="")
# --fitsiteA; LAF fit dataset for site A
parser.add_argument('--fitsiteA', action="store", dest='_fsA_', default='_NONE_')
# --fitsiteB; LAF fit dataset for site B
parser.add_argument('--fitsiteB', action="store", dest='_fsB_', default='_NONE_')
# --pdfparameter: the statistical moments of interest
parser.add_argument('--pdfparameter', action="store", dest='_pp_', default='_NONE_')
# --output: output file name for the adaptive fit results output
parser.add_argument('--output', action="store", dest='_ou_', default='')
args = parser.parse_args()

#
print("Executing: ", sys.argv[0])
print()

# ##### process site A
# read result of adaptive fitting for site A
_afsA_ = pd.read_csv(args._fsA_, sep=' ')
# drop unnecessary columns
_afsA2_ = _afsA_.loc[:, ['mu_r', 'gamma_r', 'radius', 'mu', 'gamma', f"{args._pp_}_median", 'site']]
# set renormilezed my and gamma as indeces (for join)
_afsA2_.set_index(['mu_r', 'gamma_r'], inplace=True)
# rename columns of interest (for join)
_afsA2_.rename(columns={'radius': 'A_radius', 'mu': 'A_mu', 'gamma': 'A_gamma',
                        f"{args._pp_}_median": f"A_{args._pp_}_median", 'site': 'siteA'}, inplace=True)

# #### process site B
# read result of adaptive fitting for site B
_afsB_ = pd.read_csv(args._fsB_, sep=' ')
# drop unnecessary columns
_afsB2_ = _afsB_.loc[:, ['mu_r', 'gamma_r', 'radius', 'mu', 'gamma', f"{args._pp_}_median", 'site']]
# set renormilezed my and gamma as indeces (for join)
_afsB2_.set_index(['mu_r', 'gamma_r'], inplace=True)
# rename columns of interest (for join)
_afsB2_.rename(columns={'radius': 'B_radius', 'mu': 'B_mu', 'gamma': 'B_gamma',
                        f"{args._pp_}_median": f"B_{args._pp_}_median", 'site': 'siteB'}, inplace=True)

# joined insitu results for both sites
_joined_ = _afsA2_.join(_afsB2_, how='inner')

# calculate semi-absolute distance between the two curves point by point
_joined_.loc[:, 'abs_semi_diff'] = abs(_joined_[f"A_{args._pp_}_median"]-_joined_[f"B_{args._pp_}_median"])/2
# calculate average curve point by point
_joined_.loc[:, 'mean_value'] = abs(_joined_[f"A_{args._pp_}_median"]+_joined_[f"B_{args._pp_}_median"])/2
# calculate weight point by point
_joined_.loc[:, 'weight'] = 1/(_joined_['A_radius']*_joined_['B_radius'])

# ##### L2 RD formalism: calculate point by point contribution to no-weighted/weighted distances
_joined_.loc[:, 'numerator'] = (_joined_['abs_semi_diff'])*(_joined_['abs_semi_diff'])
_joined_.loc[:, 'denominator'] = _joined_['mean_value']*_joined_['mean_value']
_joined_.loc[:, 'numerator_w'] = (_joined_['abs_semi_diff'])*(_joined_['abs_semi_diff'])*_joined_['weight']
_joined_.loc[:, 'denominator_w'] = _joined_['mean_value']*_joined_['mean_value']*_joined_['weight']
# aggregate point by point contribution
num = _joined_[['numerator']].sum()[0]
den = _joined_[['denominator']].sum()[0]
num_w = _joined_[['numerator_w']].sum()[0]
den_w = _joined_[['denominator_w']].sum()[0]

# ##### RD formalism: calculate point by point RD
_joined_.loc[:, 'relative_discrepancy'] = (_joined_.abs_semi_diff)/_joined_.mean_value
# calculate no-weighted/weighted of RD
median_RD = _joined_[['relative_discrepancy']].median()[0]
median_RD_w = ws.weighted_median(_joined_['relative_discrepancy'], _joined_['weight'])

# ##### prepare couple results dataframe
results_couple_only = pd.DataFrame(columns=['site_A', 'site_B', 'support_perc_A', 'support_perc_B',
                                            'L2RD', 'L2RD_w', 'median_RD', 'median_RD_w'])
results_couple_only = results_couple_only.astype(
                      dtype={'site_A': 'object', 'site_B': 'object', 'support_perc_A': 'float64',
                             'support_perc_B': 'float64', 'L2RD': 'float64', 'L2RD_w': 'float64',
                             'median_RD': 'float64', 'median_RD_w': 'float64'})
toappend_dict = {'site_A': _joined_['siteA'].unique(), 'site_B': _joined_['siteB'].unique(),
                 'support_perc_A': [np.round(_joined_.shape[0]/_afsA_.shape[0], 3)],
                 'support_perc_B': [np.round(_joined_.shape[0]/_afsB_.shape[0], 3)],
                 'L2RD': [np.round(np.sqrt(num)/np.sqrt(den), 4)],
                 'L2RD_w': [np.round(np.sqrt(num_w)/np.sqrt(den_w), 4)],
                 'median_RD': [np.round(median_RD, 4)],
                 'median_RD_w': [np.round(median_RD_w, 4)]}
results_couple_only = pd.concat([results_couple_only, pd.DataFrame(toappend_dict)], ignore_index=True)

# ##### prepare both way results dataframe
results_both = pd.DataFrame(columns=['site_A', 'site_B', 'support_perc_A', 'support_perc_B',
                                     'L2RD', 'L2RD_w', 'median_RD', 'median_RD_w'])
results_both = results_both.astype(
               dtype={'site_A': 'object', 'site_B': 'object', 'support_perc_A': 'float64',
                      'support_perc_B': 'float64', 'L2RD': 'float64', 'L2RD_w': 'float64',
                      'median_RD': 'float64', 'median_RD_w': 'float64'})
results_both = cp.deepcopy(results_couple_only)
toappend_dict = {'site_A': _joined_['siteA'].unique(), 'site_B': _joined_['siteB'].unique(),
                 'support_perc_A': [np.round(_joined_.shape[0]/_afsB_.shape[0], 3)],
                 'support_perc_B': [np.round(_joined_.shape[0]/_afsA_.shape[0], 3)],
                 'L2RD': [np.round(np.sqrt(num)/np.sqrt(den), 4)],
                 'L2RD_w': [np.round(np.sqrt(num_w)/np.sqrt(den_w), 4)],
                 'median_RD': [np.round(median_RD, 4)],
                 'median_RD_w': [np.round(median_RD_w, 4)]}
results_both = pd.concat([results_both, pd.DataFrame(toappend_dict)], ignore_index=True)

# output results
results_couple_only.to_csv(args._ou_+'#couples', index=None)
results_both.to_csv(args._ou_+'#both', index=None)
