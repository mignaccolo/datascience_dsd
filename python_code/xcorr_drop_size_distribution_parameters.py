# This program calculates the cross correlation among drop size distribution parameters (statistical moments)
# for all possible pdf representations
# drop size distribution parameters are
#       N (drop count) and first 6 statistical moments of p(D) for the flux representation
#       N_{V} (number concentration) and first 6 statistical moments of f(D) for the cloudexpo, cloudplaw, or cloud2dvd representation
#
# Output 1: cross-correlation between statistical moments of the pdf of drop diameters (matrix format)
# Output 2: cross-correlation between statistical moments of the pdf of drop diameters (for plot format)
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

# The disdrorain package is expected to be in the same directory where clean_non2DVD_data.py is executed
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
# --xcorroutmatrix: output file name for cross correlation results in matrix format
#                   default is "xcorrmatrix_dsd_paramters"
parser.add_argument('--xcorroutmatrix', action="store", dest='_dom_', default='xcorrmatrix_dsd_paramters')
# --xcorroutplot: output file name for cross correlation results for plot format
#                 default is "xcorrplot_dsd_paramters"
parser.add_argument('--xcorroutplot', action="store", dest='_dop_', default='xcorrplot_dsd_paramters')
args = parser.parse_args()

#
print("Executing: ", sys.argv[0])
print()


# Method for converting correlation matrix into dataset ideal for plotting
def xcorr_4plot(_xcorr_matrix_, _pdftype_):

    _df_ = pd.DataFrame(columns=['N/Nv_mu', 'N/Nv_sigma', 'N/Nv_gamma', 'N/Nv_kappa', 'N/Nv_eta', 'mu_sigma', 'mu_gamma',
                                 'mu_kappa', 'mu_eta', 'sigma_gamma', 'sigma_kappa', 'sigma_eta', 'gamma_kappa',
                                 'gamma_eta', 'kappa_eta', 'pdftype', 'acronym'])
    _df_ = _df_.astype(dtype={'N/Nv_mu': 'float64', 'N/Nv_sigma': 'float64', 'N/Nv_gamma': 'float64', 'N/Nv_kappa': 'float64',
                              'N/Nv_eta': 'float64', 'mu_sigma': 'float64', 'mu_gamma': 'float64', 'mu_kappa': 'float64',
                              'mu_eta': 'float64', 'sigma_gamma': 'float64', 'sigma_kappa': 'float64', 'sigma_eta': 'float64',
                              'gamma_kappa': 'float64', 'gamma_eta': 'float64', 'kappa_eta': 'float64',
                              'pdftype': 'object', 'acronym': 'object'})

    if _pdftype_ == "flux":
        ndrops = 'N'
    else:
        ndrops = 'Nv'

    _df_['N/Nv_mu'] = _xcorr_matrix_.loc[(_xcorr_matrix_.Variable == f"{ndrops}"), 'mu'].values
    _df_['N/Nv_sigma'] = _xcorr_matrix_.loc[(_xcorr_matrix_.Variable == f"{ndrops}"), 'sigma'].values
    _df_['N/Nv_gamma'] = _xcorr_matrix_.loc[(_xcorr_matrix_.Variable == f"{ndrops}"), 'gamma'].values
    _df_['N/Nv_kappa'] = _xcorr_matrix_.loc[(_xcorr_matrix_.Variable == f"{ndrops}"), 'kappa'].values
    _df_['N/Nv_eta'] = _xcorr_matrix_.loc[(_xcorr_matrix_.Variable == f"{ndrops}"), 'eta'].values
    _df_['mu_sigma'] = _xcorr_matrix_.loc[(_xcorr_matrix_.Variable == 'mu'), 'sigma'].values
    _df_['mu_gamma'] = _xcorr_matrix_.loc[(_xcorr_matrix_.Variable == 'mu'), 'gamma'].values
    _df_['mu_kappa'] = _xcorr_matrix_.loc[(_xcorr_matrix_.Variable == 'mu'), 'kappa'].values
    _df_['mu_eta'] = _xcorr_matrix_.loc[(_xcorr_matrix_.Variable == 'mu'), 'eta'].values
    _df_['sigma_gamma'] = _xcorr_matrix_.loc[(_xcorr_matrix_.Variable == 'sigma'), 'gamma'].values
    _df_['sigma_kappa'] = _xcorr_matrix_.loc[(_xcorr_matrix_.Variable == 'sigma'), 'kappa'].values
    _df_['sigma_eta'] = _xcorr_matrix_.loc[(_xcorr_matrix_.Variable == 'sigma'), 'eta'].values
    _df_['gamma_kappa'] = _xcorr_matrix_.loc[(_xcorr_matrix_.Variable == 'gamma'), 'kappa'].values
    _df_['gamma_eta'] = _xcorr_matrix_.loc[(_xcorr_matrix_.Variable == 'gamma'), 'eta'].values
    _df_['kappa_eta'] = _xcorr_matrix_.loc[(_xcorr_matrix_.Variable == 'kappa'), 'eta'].values
    _df_['pdftype'] = f"{_pdftype_}"
    _df_['acronym'] = args._da_
    _df_['instrument'] = instr

    return _df_


# load catalog
catalog = pd.read_csv(args._dc_, sep=',')
# retrieve index of catalog data frame corresponding to the desired acronym
index_acronym = catalog[catalog['ID2'] == args._da_].index.tolist()

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
    # calculate central moments in the flux (ground) representation
    dsdpar_flux = disdrodata.flux_phase_space_parameters()

    # calculate cross correlation in the flux (ground) representation
    xcorr_dsdpar_flux = dsdpar_flux.corr().reset_index()
    xcorr_dsdpar_flux.rename(columns={'index': 'Variable', 'N': 'N/Nv'}, inplace=True)
    xcorr_dsdpar_flux['pdftype'] = 'flux'
    # prepare for plot dataframe
    _pflux = xcorr_4plot(xcorr_dsdpar_flux, "flux")

    if (instr != '2DVD'):
        # calculate phase space paramter in the cloud representation
        # exponential law for drop speed
        dsdpar_cexpo = disdrodata.cloud_phase_space_parameters(_speed_='expo')
        # expoential law for drop speed
        dsdpar_cplaw = disdrodata.cloud_phase_space_parameters(_speed_='plaw')

        xcorr_dsdpar_cexpo = dsdpar_cexpo.corr().reset_index()
        xcorr_dsdpar_cexpo.rename(columns={'index': 'Variable', 'Nv': 'N/Nv'}, inplace=True)
        xcorr_dsdpar_cexpo['pdftype'] = 'cloudexpo'

        xcorr_dsdpar_cplaw = dsdpar_cplaw.corr().reset_index()
        xcorr_dsdpar_cplaw.rename(columns={'index': 'Variable', 'Nv': 'N/Nv'}, inplace=True)
        xcorr_dsdpar_cplaw['pdftype'] = 'cloudplaw'

        xcorr_dsdpar_summary = xcorr_dsdpar_flux.copy()
        xcorr_dsdpar_summary = pd.concat([xcorr_dsdpar_summary, xcorr_dsdpar_cplaw], ignore_index=True)
        xcorr_dsdpar_summary = pd.concat([xcorr_dsdpar_summary, xcorr_dsdpar_cexpo], ignore_index=True)
        xcorr_dsdpar_summary['site'] = args._da_
        xcorr_dsdpar_summary['instrument'] = instr

        _pcexpo = xcorr_4plot(xcorr_dsdpar_cexpo, "cloudexpo")
        _pcplaw = xcorr_4plot(xcorr_dsdpar_cplaw, "cloudplaw")
        xcorr_dsdpar_summary_xplot = _pflux.copy()
        xcorr_dsdpar_summary_xplot = pd.concat([xcorr_dsdpar_summary_xplot, _pcexpo], ignore_index=True)
        xcorr_dsdpar_summary_xplot = pd.concat([xcorr_dsdpar_summary_xplot, _pcplaw], ignore_index=True)

    else:
        # calculate phase space paramter in the cloud representation
        dsdpar_cloud = disdrodata.cloud_phase_space_parameters()
        xcorr_dsdpar_cloud = dsdpar_cloud.corr().reset_index()
        xcorr_dsdpar_cloud.rename(columns={'index': 'Variable', 'Nv': 'N/Nv'}, inplace=True)
        xcorr_dsdpar_cloud['pdftype'] = 'cloud2dvd'

        xcorr_dsdpar_summary = xcorr_dsdpar_flux.copy()
        xcorr_dsdpar_summary = pd.concat([xcorr_dsdpar_summary, xcorr_dsdpar_cloud], ignore_index=True)
        xcorr_dsdpar_summary['site'] = args._da_
        xcorr_dsdpar_summary['instrument'] = instr

        _pcloud = xcorr_4plot(xcorr_dsdpar_cloud, "cloud2dvd")
        xcorr_dsdpar_summary_xplot = _pflux.copy()
        xcorr_dsdpar_summary_xplot = pd.concat([xcorr_dsdpar_summary_xplot, _pcloud], ignore_index=True)

else:
    print("ERROR: ID2 not found in catalog")

xcorr_dsdpar_summary.to_csv(args._dom_, sep=' ', index=None)
xcorr_dsdpar_summary_xplot.to_csv(args._dop_, sep=' ', index=None)
