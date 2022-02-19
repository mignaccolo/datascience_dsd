# This program implement the Local Adaptive Fit Algorithm
#
# Output 1: LAF fit estimate for each grid point in the mu-gamma space
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
# --disdrodatapath: path to directory contain the disdrometer data
parser.add_argument('--disdrodatapath', action="store", dest='_ddp_', default='_NONE_')
# --disdrocatalog: full path to file with data catalog (data catalog contain metdata about the data
parser.add_argument('--disdrocatalog', action="store", dest='_dc_', default='_NONE_')
# --disdroacronym: acronym identifying the data set in the catalog
parser.add_argument('--disdroacronym', action="store", dest='_da_', default='_NONE_')
# --pdfparameter: the statistical moments to fit as a function of mu and gamma
parser.add_argument('--pdfparameter', action="store", dest='_pp_', default='_NONE_')
# --representation: which type of drop size representation to adopt: possible values are "flux", "cloudexpo", "cloudplaw", and "cloud2dvd"
#                   default is "flux"
parser.add_argument('--representation', action="store", dest='_re_', default='flux')
# --renorm_table: path to file containing the renormalization values table
parser.add_argument('--renorm_table', action="store", dest='_rnt_', default='_NONE_')
# --renorm_type: which type of renormalization to adopt:
#                the argument of --renorm_type must be one the values in the
#                column "renorm_type" column of the file with the renormalization values table
parser.add_argument('--renorm_type', action="store", dest='_rnty_', default='')
# --radiusseq
parser.add_argument('--radiusseq', action="store", dest='_rs_', default='_NONE_')
# --occupancy
parser.add_argument('--occupancy', action="store", dest='_oc_', default='20')
# --output: output file name for the adaptive fit results output
parser.add_argument('--output', action="store", dest='_ou_', default='')
args = parser.parse_args()

#
print("Executing: ", sys.argv[0])
print()

catalog = pd.read_csv(args._dc_, sep=',')


# Method calculating the renormalized values of the statistical moments
def _renormalize_(_dfin_, _rntable_, _rnrmtype_, _param_):
    # renormalize the values of statistical moments
    _df_ = copy.deepcopy(_dfin_)
    # obtain renormalization parameters for mean
    mu_min = _rntable_.loc[((_rntable_['renorm_type'] == _rnrmtype_) & (_rntable_['statistical_moment'] == 'mu')), 'xmin'].values[0]
    mu_max = _rntable_.loc[((_rntable_['renorm_type'] == _rnrmtype_) & (_rntable_['statistical_moment'] == 'mu')), 'xmax'].values[0]
    # obtain renormalization parameters for standard deviation
    sigma_min = _rntable_.loc[((_rntable_['renorm_type'] == _rnrmtype_) & (_rntable_['statistical_moment'] == 'sigma')), 'xmin'].values[0]
    sigma_max = _rntable_.loc[((_rntable_['renorm_type'] == _rnrmtype_) & (_rntable_['statistical_moment'] == 'sigma')), 'xmax'].values[0]
    # obtain renormalization parameters for skewness
    gamma_min = _rntable_.loc[((_rntable_['renorm_type'] == _rnrmtype_) & (_rntable_['statistical_moment'] == 'gamma')), 'xmin'].values[0]
    gamma_max = _rntable_.loc[((_rntable_['renorm_type'] == _rnrmtype_) & (_rntable_['statistical_moment'] == 'gamma')), 'xmax'].values[0]
    # obtain renormalization parameters for kurtosis
    kappa_min = _rntable_.loc[((_rntable_['renorm_type'] == _rnrmtype_) & (_rntable_['statistical_moment'] == 'kappa')), 'xmin'].values[0]
    kappa_max = _rntable_.loc[((_rntable_['renorm_type'] == _rnrmtype_) & (_rntable_['statistical_moment'] == 'kappa')), 'xmax'].values[0]
    # obtain renormalization parameters for 5th central moment
    eta_min = _rntable_.loc[((_rntable_['renorm_type'] == _rnrmtype_) & (_rntable_['statistical_moment'] == 'eta')), 'xmin'].values[0]
    eta_max = _rntable_.loc[((_rntable_['renorm_type'] == _rnrmtype_) & (_rntable_['statistical_moment'] == 'eta')), 'xmax'].values[0]

    # organizing renormalization ranges for the paramters that can be fitted
    dict_min = {'sigma': sigma_min, 'kappa': kappa_min, 'eta': eta_min}
    dict_max = {'sigma': sigma_max, 'kappa': kappa_max, 'eta': eta_max}

    # calculate renormalized statistical moments
    _df_.loc[:, 'mu_r'] = (_df_.mu-mu_min)/(mu_max-mu_min)
    _df_.loc[:, 'gamma_r'] = (_df_.gamma-gamma_min)/(gamma_max-gamma_min)
    _df_.loc[:, f"{_param_}_r"] = (_df_[_param_]-dict_min[_param_])/(dict_max[_param_]-dict_min[_param_])

    return _df_


# Method for obtaining values of central moments from re-normalized values of central moments
def _unrenormalize_laf_results_(_dfin_, _rntable_, _rnrmtype_, _param_):
    # renormalize the values of central moments
    _df_ = copy.deepcopy(_dfin_)
    # obtain renormalization parameters for mean
    mu_min = _rntable_.loc[((_rntable_['renorm_type'] == _rnrmtype_) & (_rntable_['statistical_moment'] == 'mu')), 'xmin'].values[0]
    mu_max = _rntable_.loc[((_rntable_['renorm_type'] == _rnrmtype_) & (_rntable_['statistical_moment'] == 'mu')), 'xmax'].values[0]
    # obtain renormalization parameters for standard deviation
    sigma_min = _rntable_.loc[((_rntable_['renorm_type'] == _rnrmtype_) & (_rntable_['statistical_moment'] == 'sigma')), 'xmin'].values[0]
    sigma_max = _rntable_.loc[((_rntable_['renorm_type'] == _rnrmtype_) & (_rntable_['statistical_moment'] == 'sigma')), 'xmax'].values[0]
    # obtain renormalization parameters for skewness
    gamma_min = _rntable_.loc[((_rntable_['renorm_type'] == _rnrmtype_) & (_rntable_['statistical_moment'] == 'gamma')), 'xmin'].values[0]
    gamma_max = _rntable_.loc[((_rntable_['renorm_type'] == _rnrmtype_) & (_rntable_['statistical_moment'] == 'gamma')), 'xmax'].values[0]
    # obtain renormalization parameters for kurtosis
    kappa_min = _rntable_.loc[((_rntable_['renorm_type'] == _rnrmtype_) & (_rntable_['statistical_moment'] == 'kappa')), 'xmin'].values[0]
    kappa_max = _rntable_.loc[((_rntable_['renorm_type'] == _rnrmtype_) & (_rntable_['statistical_moment'] == 'kappa')), 'xmax'].values[0]
    # obtain renormalization parameters for 5th central moment
    eta_min = _rntable_.loc[((_rntable_['renorm_type'] == _rnrmtype_) & (_rntable_['statistical_moment'] == 'eta')), 'xmin'].values[0]
    eta_max = _rntable_.loc[((_rntable_['renorm_type'] == _rnrmtype_) & (_rntable_['statistical_moment'] == 'eta')), 'xmax'].values[0]

    # organizing renormalization ranges for the paramters that can be fitted
    dict_min = {'sigma': sigma_min, 'kappa': kappa_min, 'eta': eta_min}
    dict_max = {'sigma': sigma_max, 'kappa': kappa_max, 'eta': eta_max}

    _df_.loc[:, 'mu'] = (_df_.mu_r*(mu_max-mu_min))+mu_min
    _df_.loc[:, 'gamma'] = (_df_.gamma_r*(gamma_max-gamma_min))+gamma_min

    _df_.loc[:, f"{_param_}_median"] = np.round((fitted.predicted_r*(dict_max[_param_]-dict_min[_param_]))+dict_min[_param_], 6)
    _df_.loc[:, f"{_param_}_5%"] = np.round((fitted['predicted_5%']*(dict_max[_param_]-dict_min[_param_]))+dict_min[_param_], 6)
    _df_.loc[:, f"{_param_}_q1"] = np.round((fitted.predicted_q1*(dict_max[_param_]-dict_min[_param_]))+dict_min[_param_], 6)
    _df_.loc[:, f"{_param_}_q3"] = np.round((fitted.predicted_q3*(dict_max[_param_]-dict_min[_param_]))+dict_min[_param_], 6)
    _df_.loc[:, f"{_param_}_95%"] = np.round((fitted['predicted_95%']*(dict_max[_param_]-dict_min[_param_]))+dict_min[_param_], 6)

    return _df_


# Method implenting the LAF fitting procedure
def _build_model_(_df_, _param_, _listradius_, _occupancy_):
    _data_ = _df_[['mu_r', 'gamma_r']].values
    _params_ = _df_[[_param_+'_r']].values
    _mu_coor = np.mgrid[-0.51:2.01:0.02]
    _gamma_coor = np.mgrid[-0.51:2.01:0.02]
    _centers_ = np.array(np.meshgrid(_mu_coor, _gamma_coor)).T.reshape(-1, 2)

    gsp = GriSPy(_data_)

    results = pd.DataFrame(columns=['mu_r', 'gamma_r', 'predicted_r', 'radius'])
    results = results.astype(dtype={'mu_r': 'float64', 'gamma_r': 'float64', 'predicted_r': 'float64', 'radius': 'float64'})

    for _radius_ in _listradius_:
        print("processing radius:", _radius_)
        # Query for neighbors within upper_radii
        bubble_dist, bubble_ind = gsp.bubble_neighbors(
            _centers_, distance_upper_bound=float(_radius_))

        _squareindex_array_ = []
        _dataindex_array_ = []
        j = 0
        for _item_ in bubble_ind:
            if (len(_item_) >= int(_occupancy_)):
                _squareindex_array_.append(j)
                _dataindex_array_.append(_item_)
            j = j+1

        parammedian = []
        paramq1 = []
        paramq3 = []
        param5 = []
        param95 = []
        for _indexsubset_ in _dataindex_array_:
            _params_sel = np.array([_params_[index] for index in _indexsubset_])
            parammedian.append(np.round(np.median(_params_sel), 6))
            paramq1.append(np.round(np.quantile(_params_sel, 0.25), 6))
            paramq3.append(np.round(np.quantile(_params_sel, 0.75), 6))
            param5.append(np.round(np.quantile(_params_sel, 0.05), 6))
            param95.append(np.round(np.quantile(_params_sel, 0.95), 6))

        _centers_sel = np.array([_centers_[index] for index in _squareindex_array_])
        xy = np.vstack(_centers_sel)
        x = np.round(xy[:, 0], 6)
        y = np.round(xy[:, 1], 6)
        z = np.array(parammedian)
        zq1 = np.array(paramq1)
        zq3 = np.array(paramq3)
        zq5 = np.array(param5)
        zq95 = np.array(param95)

        _temp_ = pd.DataFrame({'mu_r': x, 'gamma_r': y, 'predicted_r': z, 'predicted_5%': zq5,
                               'predicted_q1': zq1, 'predicted_q3': zq3, 'predicted_95%': zq95})
        _temp_.loc[:, 'radius'] = _radius_
        # print(_temp_.head(10))
        # results = results.append(_temp_)
        results = pd.concat([results, _temp_], ignore_index=True)
        _centers_ = _centers_[~np.isin(np.arange(len(_centers_)), _squareindex_array_)]
        # print("3",_radius_,_centers_.shape)

    return results


# load catalog
catalog = pd.read_csv(args. _dc_, sep=',')
# retrieve index of catalog data frame corresponding to the desired acronym
index_acronym = catalog[catalog['ID2'] == args._da_].index.tolist()
# load renormalization parameters`table
renorm_table = pd.read_csv(args._rnt_, sep=',')
# load radii sequence
radius = []
with open(args._rs_, "r") as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for lines in csv_reader:
        radius = lines

# check if representation value is admitted (only: flux, cloudexpo, cloudplaw)
if args._re_ not in ['flux', 'cloudexpo', 'cloudplaw', 'cloud2dvd']:
    print("ERROR: representation value not found!")
    print("possible choices are flux, cloudexpo, cloudplaw, cloud2dvd")
    sys.exit()

# if acronym is part of the catalog we proceed
if index_acronym:
    row = catalog.loc[catalog.ID2 == args._da_, :]  # retrieve row corresponding to acronym
    celllimits_path = args._ddp_+'/'+row['CELLLIMITS'].values[0]  # set file name with cell limits values
    areainstr = row['AREA_INSTRUMENT'].values[0]  # set area of instrument
    instr = row['INSTRUMENT'].values[0]  # set instrument type
    tr = row['TIME_RESOLUTION'].values[0]  # set instrument time resolution
    pathtodata = args._ddp_+'/'+args._da_  # set path to data set

    # if data set is not from 2DVD disdrometer
    if (instr != '2DVD'):
        # if data set is RD80 with standard cell limits division
        if ((instr == 'RD80') & (row['CELLLIMITS'].values[0] == 'standard')):
            disdrodata = dr.disdrorain(datapath=pathtodata, instrument_area=areainstr, time_interval=tr)  # create disdrodata class
        else:
            disdrodata = dr.disdrorain(classpath=celllimits_path, datapath=pathtodata,
                                       instrument_area=areainstr, time_interval=tr)  # create disdrodata class
    else:
        disdrodata = dr.disdrorain_2dvd(datapath=pathtodata, instrument_area=areainstr,
                                        time_interval=tr)  # create disdrodata class for 2dvd data

    # calculate dsd paramters
    if args._re_ == 'flux':
        # calculate central moments in the flux (ground) representation
        dsdpar = disdrodata.flux_phase_space_parameters()
    else:
        # calculate phase space paramter in the cloud representation for on 2DVD data
        if args._re_ == 'cloudexpo':
            # exponential law for drop speed
            dsdpar = disdrodata.cloud_phase_space_parameters(_speed_='expo')
        if args._re_ == 'cloudplaw':
            # expoential law for drop speed
            dsdpar = disdrodata.cloud_phase_space_parameters(_speed_='plaw')
        if args._re_ == 'cloud2dvd':
            # this are 2dvd data and we know the speed of each drop
            dsdpar = disdrodata.cloud_phase_space_parameters()

    # renormalize drop size distribution parameters
    dsdpar_r = _renormalize_(dsdpar, renorm_table, args._rnty_, args._pp_)
    # perform LAF
    fitted = _build_model_(dsdpar_r, args._pp_, radius, args._oc_)

else:
    print("Acronym not found in catalog")
    sys.exit()

# un-renormalize results
results = _unrenormalize_laf_results_(fitted, renorm_table, args._rnty_, args._pp_)
# add site info to results
results.loc[:, 'site'] = args._da_
# save results to file
results.to_csv(args._ou_, sep=' ', index=None)
