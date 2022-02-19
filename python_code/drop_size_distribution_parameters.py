# This program calculates the drop size distribution parameters (statistical moments)
# drop size distribution parameters are
#       N (drop count) and first 6 statistical moments of p(D) for the flux representation
#       N_{V} (number concentration) and first 6 6 statistical moments of f(D) for the cloudexpo or cloudplaw representation
#
# Possible representations of the drop size distribution area available
#   1) "flux" or ground representation: couple (N,p(D))
#   2) "cloudexpo" representation: couple (N_{V},f(D)). If speed of drop is not known than v(D)=9.65-10.3\times e^{-0.6D} is assumed
#   3) "cloudplaw" representation: couple (N_{V},f(D)). If speed of drop is not known than v(D)=3.78\times D^{0.67} is assumed
#   4) "cloud2dvd" representation: couple (N_{V},f(D)). Used only for 2dvd data: drop velocity is the velocity reported by the instrument

# Renormalized values of drop size parameters can be calculated by declaring the path to a renormalization table
#       renormalization apply to the first 5 central moment of p(D) or f(D)
#       no renormalization for N or N_{V}
#       no renormalization for the 6-th statistical moment
#
# Output 1: statistical moments of the pdf of drop diameters
#
# Author: Massimiliano Ignaccolo
#         massimiliano.ignaccolo [at] tutanota.com. massimiliano.ignaccolo [at] sas.com
# Last modified: Feb 01 2022

# ------------- Necessary Python packages -START
import os
import sys
import copy
import pandas as pd
import argparse

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
# --disdrocatalog: full path to file with data catalog (data catalog contain metdata about the dataset)
parser.add_argument('--disdrocatalog', action="store", dest='_dc_', default='_NONE_')
# --disdroacronym: acronym identifying the data set in the catalog
parser.add_argument('--disdroacronym', action="store", dest='_da_', default='_NONE_')
# --representation: which type of drop size representation to adopt: possible values are "flux", "cloudexpo", "cloudplaw", and "cloud2dvd"
#                   default is "flux"
parser.add_argument('--representation', action="store", dest='_re_', default='flux')
# --renorm: default = 'N' -> renormalized central moments are not calculated
#           if different from default then value of argument is assumed to be path to file with
#           renormalization values table. Also renormalized paramters are calculated
parser.add_argument('--renorm', action="store", dest='_rn_', default='N')
# --renorm_type: which type of renormalization to adopt:
#                used only if --renorm is different from 'N'
#                in this case the argument of --renorm_type must be one the values in the
#                column "renorm_type" column of the file with the renormalization values table
parser.add_argument('--renorm_type', action="store", dest='_rnt_', default='_NONE_')
# --output: output file name for the drop size distribution parameters
#           default value is "dsd_paramters"
#           default output is generated in the directory from which code is invoked
parser.add_argument('--output', action="store", dest='_ou_', default='dsd_parameters')
args = parser.parse_args()


#
print("Executing: ", sys.argv[0])
print()


# Method calculating the renormalized values of the statistical moments
def _renormalize_(_dfin_, _rntable_, _rnrmtype_):
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

    # calculate renormalized statistical moments
    _df_.loc[:, 'mu_r'] = (_df_.mu-mu_min)/(mu_max-mu_min)
    _df_.loc[:, 'gamma_r'] = (_df_.gamma-gamma_min)/(gamma_max-gamma_min)
    _df_.loc[:, 'sigma_r'] = (_df_.sigma-sigma_min)/(sigma_max-sigma_min)
    _df_.loc[:, 'kappa_r'] = (_df_.kappa-kappa_min)/(kappa_max-kappa_min)
    _df_.loc[:, 'eta_r'] = (_df_.eta-eta_min)/(eta_max-eta_min)

    return _df_


# load catalog
catalog = pd.read_csv(args._dc_, sep=',')
# retrieve index of catalog data frame corresponding to the desired acronym
index_acronym = catalog[catalog['ID2'] == args._da_].index.tolist()

# check if representation value is admitted (only: flux, cloudexpo, cloudplaw)
if args._re_ not in ['flux', 'cloudexpo', 'cloudplaw', 'cloud2dvd']:
    print("ERROR: representation value not found!")
    print("possible choices are flux, cloudexpo, cloudplaw, cloud2dvd")
    sys.exit()

# if acronym is part of the catalog we calculate the phase space parameters
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
        # calculate phase space paramter in the cloud representation
        if args._re_ == 'cloudexpo':
            # exponential law for drop speed
            dsdpar = disdrodata.cloud_phase_space_parameters(_speed_='expo')
        if args._re_ == 'cloudplaw':
            # expoential law for drop speed
            dsdpar = disdrodata.cloud_phase_space_parameters(_speed_='plaw')
        if args._re_ == 'cloud2dvd':
            # this are 2dvd data and we know the speed of each drop
            dsdpar = disdrodata.cloud_phase_space_parameters()

    # if renormalized value are desired
    if args._rn_ != 'N':
        # load renormalization paramters table
        renorm_table = pd.read_csv(args._rn_, sep=',')
        renorm_type = args._rnt_

        # Add the values of renormalized values of the central moments (no renormalization happens for N / N_{V})
        dsdpar_r = _renormalize_(dsdpar, renorm_table, renorm_type)
else:
    print("ERROR: ID2 not found in catalog")

# output results to file
if args._rn_ != 'N':
    dsdpar_r.to_csv(args._ou_, sep=' ', index=None)
else:
    dsdpar.to_csv(args._ou_, sep=' ', index=None)
