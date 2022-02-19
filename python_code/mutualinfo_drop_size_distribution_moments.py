# This program calculates the cross correlation among drop size distribution moments
# for all possible pdf representations
# drop size distribution parameters are
#       N (drop count) and first 6 statistical moments of p(D) for the flux representation
#       N_{V} (number concentration) and first 6 statistical moments of f(D) for the cloudexpo, cloudplaw, or cloud2dvd representation
#
# Author: Massimiliano Ignaccolo
#         massimiliano.ignaccolo [at] tutanota.com. massimiliano.ignaccolo [at] sas.com
# Last modified: Feb 01 2022
#
# Note: RMI calculation might be slow depending on the size of data:
#       program print to standard output its progress
#
# Output 1: rescaled mutual infromation between moments of the pdf of drop diameters (matrix format)
# Output 2: rescaled mutual infromation between moments of the pdf of drop diameters (for plot format)
#
# Author: Massimiliano Ignaccolo
#         massimiliano.ignaccolo [at] tutanota.com. massimiliano.ignaccolo [at] sas.com
# Last modified: Feb 01 2022

# ------------- Necessary Python packages -START
import os
import sys
import copy
import pandas as pd
import numpy as np
import argparse

# The disdrorain package is expected to be in the same directory where program is executed
# The LNC mutual infromation package is expected to be in the same directory where program is executed
# If this is not the case change value of pack_path accordingly
# pack_path = os.getcwd()
# if pack_path not in sys.path:
# sys.path.append(pack_path)
import disdrorain as dr
from lnc import MI
#  ------------- Necessary Python packages -END

# ARGUMENTS
parser = argparse.ArgumentParser(description='', epilog="")
# --disdrodatapath: path to directory contain the disdrometer data
parser.add_argument('--disdrodatapath', action="store", dest='_ddp_', default='_NONE_')
# --disdrocatalog: full path to file with data catalog (data catalog contain metdata about the data
parser.add_argument('--disdrocatalog', action="store", dest='_dc_', default='_NONE_')
# --disdroacronym: acronym identifying the data set in the catalog
parser.add_argument('--disdroacronym', action="store", dest='_da_', default='_NONE_')
# --rmioutmatrix: output file name for rescaled mutual information results in matrix format
#                 default is "rmimatrix_dsd_paramters"
parser.add_argument('--rmioutmatrix', action="store", dest='_dom_', default='rmimatrix_dsd_paramters')
# --rmioutplot: output file name for rescaled mutual information results for plot format
#               default is "rmiplot_dsd_paramters"
parser.add_argument('--rmioutplot', action="store", dest='_dop_', default='rmiplot_dsd_paramters')
args = parser.parse_args()


# Method for calculating the Rescaled Mututal Information
def calculate_RMI(_temp_):
    listv = ['N/Nv', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6']
    mi_matrix = np.zeros(shape=(len(listv), len(listv)))
    for k in range(0, len(listv)):
        for l in range(k, len(listv)):
            if k != l:
                print(listv[k], listv[l])  # info about current couple of paramteter processed
                x = _temp_.loc[:, listv[k]].values.tolist()
                y = _temp_.loc[:, listv[l]].values.tolist()

                xs = sorted(x)
                ys = sorted(y)

                mi_matrix[k][l] = MI.mi_LNC([x, y], k=20, base=np.exp(1), alpha=0.25)/MI.mi_LNC([xs, ys], k=20, base=np.exp(1), alpha=0.25)
                mi_matrix[l][k] = mi_matrix[k][l]
            else:
                mi_matrix[k][l] = 0

    _res_ = pd.DataFrame(mi_matrix)
    _res_.columns = listv
    _res_.index = listv
    _res_.reset_index(inplace=True)

    return _res_


# Method for converting RMI matrix into dataset ideal for plotting
def rmi_4plot(_temp_, _pdftype_='flux'):

    _df_ = pd.DataFrame(columns=['N/Nv_M1', 'N/Nv_M2', 'N/Nv_M3', 'N/Nv_M4', 'N/Nv_M5', 'N/Nv_M6',
                                 'M1_M2', 'M1_M3', 'M1_M4', 'M1_M5', 'M1_M6',
                                 'M2_M3', 'M2_M4', 'M2_M5', 'M2_M6', 'M3_M4', 'M3_M5', 'M3_M6',
                                 'M4_M5', 'M4_M6', 'M5_M6', 'pdftype', 'acronym'])

    _df_ = _df_.astype(dtype={'N/Nv_M1': 'float64', 'N/Nv_M2': 'float64', 'N/Nv_M3': 'float64', 'N/Nv_M4': 'float64',
                              'N/Nv_M5': 'float64', 'N/Nv_M6': 'float64', 'M1_M2': 'float64', 'M1_M3': 'float64',
                              'M1_M4': 'float64', 'M1_M5': 'float64', 'M1_M6': 'float64', 'M2_M3': 'float64', 'M2_M4': 'float64',
                              'M2_M5': 'float64', 'M2_M6': 'float64', 'M3_M4': 'float64', 'M3_M5': 'float64', 'M3_M6': 'float64',
                              'M4_M5': 'float64', 'M4_M6': 'float64', 'M5_M6': 'float64', 'pdftype': 'object', 'acronym': 'object'})

    _df_['N/Nv_M1'] = _temp_.loc[(_temp_.Variable == 'N/Nv'), 'M1'].values
    _df_['N/Nv_M2'] = _temp_.loc[(_temp_.Variable == 'N/Nv'), 'M2'].values
    _df_['N/Nv_M3'] = _temp_.loc[(_temp_.Variable == 'N/Nv'), 'M3'].values
    _df_['N/Nv_M4'] = _temp_.loc[(_temp_.Variable == 'N/Nv'), 'M4'].values
    _df_['N/Nv_M5'] = _temp_.loc[(_temp_.Variable == 'N/Nv'), 'M5'].values
    _df_['N/Nv_M6'] = _temp_.loc[(_temp_.Variable == 'N/Nv'), 'M6'].values
    _df_['M1_M2'] = _temp_.loc[(_temp_.Variable == 'M1'), 'M2'].values
    _df_['M1_M3'] = _temp_.loc[(_temp_.Variable == 'M1'), 'M3'].values
    _df_['M1_M4'] = _temp_.loc[(_temp_.Variable == 'M1'), 'M4'].values
    _df_['M1_M5'] = _temp_.loc[(_temp_.Variable == 'M1'), 'M5'].values
    _df_['M1_M6'] = _temp_.loc[(_temp_.Variable == 'M1'), 'M6'].values
    _df_['M2_M3'] = _temp_.loc[(_temp_.Variable == 'M2'), 'M3'].values
    _df_['M2_M4'] = _temp_.loc[(_temp_.Variable == 'M2'), 'M4'].values
    _df_['M2_M5'] = _temp_.loc[(_temp_.Variable == 'M2'), 'M5'].values
    _df_['M2_M6'] = _temp_.loc[(_temp_.Variable == 'M2'), 'M6'].values
    _df_['M3_M4'] = _temp_.loc[(_temp_.Variable == 'M3'), 'M4'].values
    _df_['M3_M5'] = _temp_.loc[(_temp_.Variable == 'M3'), 'M5'].values
    _df_['M3_M6'] = _temp_.loc[(_temp_.Variable == 'M3'), 'M6'].values
    _df_['M4_M5'] = _temp_.loc[(_temp_.Variable == 'M4'), 'M5'].values
    _df_['M4_M6'] = _temp_.loc[(_temp_.Variable == 'M4'), 'M6'].values
    _df_['M5_M6'] = _temp_.loc[(_temp_.Variable == 'M5'), 'M6'].values
    _df_['pdftype'] = f"{_pdftype_}"
    _df_['acronym'] = args._da_
    _df_['instrument'] = instr

    return _df_


print("Executing: ", sys.argv[0])
print()

# load catalog
catalog = pd.read_csv(args._dc_, sep=',')
# retrieve index of catalog data frame corresponding to the desired acronym
index_acronym = catalog[catalog['ID2'] == args._da_].index.tolist()

if index_acronym:
    row = catalog.loc[catalog.ID2 == args._da_, :]
    celllimits_path = args._ddp_+'/'+row['CELLLIMITS'].values[0]
    areainstr = row['AREA_INSTRUMENT'].values[0]
    instr = row['INSTRUMENT'].values[0]
    tr = row['TIME_RESOLUTION'].values[0]
    pathtodata = args._ddp_+'/'+args._da_

    # if data set is not from 2DVD disdrometer
    if (instr != '2DVD'):
        # if data set is RD80 with standard cell limits division
        if ((instr == 'RD80') & (row['CELLLIMITS'].values[0] == 'standard')):
            disdrodata = dr.disdrorain(datapath=pathtodata, instrument_area=areainstr, time_interval=tr)  # create disdrodata class
        else:
            disdrodata = dr.disdrorain(classpath=celllimits_path, datapath=pathtodata, instrument_area=areainstr,
                                       time_interval=tr)  # create disdrodata class
    else:
        disdrodata = dr.disdrorain_2dvd(datapath=pathtodata, instrument_area=areainstr,
                                        time_interval=tr)  # create disdrodata class for 2dvd data

    # list of pdf's moment to calculate
    list_moments_order = list([1, 2, 3, 4, 5, 6])

    if (instr != '2DVD'):
        # calculate bulk variables: this is used only to find the values of N or Nv
        bvar_vexpo = disdrodata.bulk_variables(_speed_='expo')
        bvar_vplaw = disdrodata.bulk_variables(_speed_='plaw')
        # calculate dsd paramters
        # calculate central moments in the flux (ground) representation
        dsdpar_flux = disdrodata.flux_moment_calculator(list_moments_order)
        dsdpar_flux['N/Nv'] = bvar_vplaw['N']  # it does not matter here if we use plaw or expo dataframe: N is always N
        # calculate phase space paramter in the cloud representation
        # exponential law for drop speed
        dsdpar_cexpo = disdrodata.cloud_moment_calculator(list_moments_order, _speed_='expo')
        dsdpar_cexpo['N/Nv'] = bvar_vexpo['Nv']
        # expoential law for drop speed
        dsdpar_cplaw = disdrodata.cloud_moment_calculator(list_moments_order, _speed_='plaw')
        dsdpar_cplaw['N/Nv'] = bvar_vplaw['Nv']

        print("1/3: processing flux representation")
        _res_flux = calculate_RMI(dsdpar_flux)
        _res_flux.rename(columns={'index': 'Variable'}, inplace=True)
        _res_flux['pdftype'] = 'flux'
        _res_flux_xplot = rmi_4plot(_res_flux)

        print("2/3: processing cloudexpo representation")
        _res_cexpo = calculate_RMI(dsdpar_cexpo)
        _res_cexpo.rename(columns={'index': 'Variable'}, inplace=True)
        _res_cexpo['pdftype'] = 'cloudexpo'
        _res_cexpo_xplot = rmi_4plot(_res_cexpo, _pdftype_="cloudexpo")

        print("3/3: processing cloudplaw representation")
        _res_cplaw = calculate_RMI(dsdpar_cplaw)
        _res_cplaw.rename(columns={'index': 'Variable'}, inplace=True)
        _res_cplaw['pdftype'] = 'cloudplaw'
        _res_cplaw_xplot = rmi_4plot(_res_cplaw, _pdftype_="cloudplaw")

        res_summary = _res_flux.copy()
        res_summary = pd.concat([res_summary, _res_cplaw], ignore_index=True)
        res_summary = pd.concat([res_summary, _res_cexpo], ignore_index=True)
        res_summary['site'] = args._da_
        res_summary['instrument'] = instr

        res_summary_xplot = _res_flux_xplot.copy()
        res_summary_xplot = pd.concat([res_summary_xplot, _res_cplaw_xplot], ignore_index=True)
        res_summary_xplot = pd.concat([res_summary_xplot, _res_cexpo_xplot], ignore_index=True)

    else:
        # calculate bulk variables: this is used only to find the values of N or Nv
        bvar = disdrodata.bulk_variables()
        # calculate dsd paramters
        # calculate central moments in the flux (ground) representation
        list_moments_order = list([1, 2, 3, 4, 5, 6])
        dsdpar_flux = disdrodata.flux_moment_calculator(list_moments_order)
        dsdpar_flux['N/Nv'] = bvar['N']
        # calculate phase space paramter in the cloud representation
        dsdpar_cloud = disdrodata.cloud_moment_calculator(list_moments_order)
        dsdpar_cloud['N/Nv'] = bvar['Nv']

        print("1/2: processing flux representation")
        _res_flux = calculate_RMI(dsdpar_flux)
        _res_flux.rename(columns={'index': 'Variable'}, inplace=True)
        _res_flux['pdftype'] = 'flux'
        _res_flux_xplot = rmi_4plot(_res_flux)

        print("2/2: processing cloud representation")
        _res_cloud = calculate_RMI(dsdpar_cloud)
        _res_cloud.rename(columns={'index': 'Variable'}, inplace=True)
        _res_cloud['pdftype'] = 'cloud2dvd'
        _res_cloud_xplot = rmi_4plot(_res_cloud, _pdftype_="cloud2dvd")

        res_summary = _res_flux.copy()
        res_summary = pd.concat([res_summary, _res_cloud], ignore_index=True)
        res_summary['site'] = args._da_
        res_summary['instrument'] = instr

        res_summary_xplot = _res_flux_xplot.copy()
        res_summary_xplot = pd.concat([res_summary_xplot, _res_cloud_xplot], ignore_index=True)
else:
    print("Acronym not found in catalog")
    sys.exit()

res_summary.to_csv(args._dom_, sep=' ', index=None)
res_summary_xplot.to_csv(args._dop_, sep=' ', index=None)
