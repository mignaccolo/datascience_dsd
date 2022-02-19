# This program perform the PCA analysis of all the drop size paramters (statistical moments)
# for all possible pdf representations
# drop size distribution parameters are
#       N (drop count) and first 6 statistical moments of p(D) for the flux representation
#       N_{V} (number concentration) and first 6 statistical moments of f(D) for the cloudexpo, cloudplaw, or cloud2dvd representation
#
# Output 1: results of PCA analysis of the statistical moments of the pdf of drop diameters
# Output 2: Coefficients of the Principal Components
#
# Author: Massimiliano Ignaccolo
#         massimiliano.ignaccolo [at] tutanota.com. massimiliano.ignaccolo [at] sas.com
# Last modified: Feb 01 2022

# ------------- Necessary Python packages -START
import sys
import pandas as pd
import numpy as np
import argparse
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

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
# --pcavarianceexp: output file name for the variance explained by each component
#                   default is "pca_varexpl_dsd_paramters"
parser.add_argument('--pcavarianceexp', action="store", dest='_pcave_', default='pca_varexpl_dsd_paramters')
# --pcacomponents: output file name for factors identifying each component
#                  default is "pca_compo_dsd_paramters"
parser.add_argument('--pcacomponents', action="store", dest='_pcac_', default='pca_compo_dsd_paramters')
args = parser.parse_args()


# Metohd for calculating the PCA
def PCA_maker(_df0_, ndropspdf='yes', pdftype='flux'):

    listv_pdf = ['mu', 'sigma', 'gamma', 'kappa', 'eta']
    if pdftype == 'flux':
        listv = ['N', 'mu', 'sigma', 'gamma', 'kappa', 'eta']
    else:
        listv = ['Nv', 'mu', 'sigma', 'gamma', 'kappa', 'eta']

    if ndropspdf == 'yes':
        listvar = listv
        _nc_ = 6
        list_index = ['PCA1', 'PCA2', 'PCA3', 'PCA4', 'PCA5', 'PCA6']
        _df_ = _df0_.loc[:, listvar]
    else:
        listvar = listv_pdf
        _nc_ = 5
        list_index = ['PCA1', 'PCA2', 'PCA3', 'PCA4', 'PCA5']
        _df_ = _df0_.loc[:, listvar]

    scaler = StandardScaler()
    scaler.fit(_df_.values)
    _scaled_ = scaler.transform(_df_.values)
    pca = PCA(n_components=_nc_)
    pca.fit(_scaled_)
    data = pca.explained_variance_ratio_.cumsum()
    compo = pca.components_

    list_columns = list()
    for i in range(1, _nc_+1):
        list_columns.append(f"PCA{i}")
    _res_ = pd.DataFrame(data=data.reshape(1, _nc_), columns=list_columns)
    _res_compo = pd.DataFrame(data=compo)
    _res_compo.columns = listvar
    _res_compo.index = list_index
    _res_compo.reset_index(inplace=True)
    if pdftype == 'flux':
        _res_compo.rename(columns={'index': 'component', 'N': 'N/Nv'}, inplace=True)
    else:
        _res_compo.rename(columns={'index': 'component', 'Nv': 'N/Nv'}, inplace=True)

    if ndropspdf == 'yes':
        _res_.loc[:, 'varused'] = 'N_pdf'
        _res_compo.loc[:, 'varused'] = 'N_pdf'
    else:
        _res_.loc[:, 'varused'] = 'pdf'
        _res_compo.loc[:, 'varused'] = 'pdf'

    return (_res_, _res_compo)


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
            disdrodata = dr.disdrorain(classpath=celllimits_path, datapath=pathtodata,
                                       instrument_area=areainstr, time_interval=tr)  # create disdrodata class
    else:
        disdrodata = dr.disdrorain_2dvd(datapath=pathtodata, instrument_area=areainstr,
                                        time_interval=tr)  # create disdrodata class for 2dvd data

    # calculate dsd paramters
    # calculate central moments in the flux (ground) representation
    dsdpar_flux = disdrodata.flux_phase_space_parameters()
    # process flux representation
    (_pca_flux, _pca_compo_flux) = PCA_maker(dsdpar_flux, ndropspdf='yes', pdftype='flux')
    (_pca_flux_noN, _pca_compo_flux_noN) = PCA_maker(dsdpar_flux, ndropspdf='no', pdftype='flux')
    res_pca_flux = pd.concat([_pca_flux, _pca_flux_noN], ignore_index=True)
    res_pca_flux['pdftype'] = 'flux'
    res_pca_compo_flux = pd.concat([_pca_compo_flux, _pca_compo_flux_noN], ignore_index=True)
    res_pca_compo_flux['pdftype'] = 'flux'

    if (instr != '2DVD'):
        # calculate phase space paramter in the cloud representation
        # exponential law for drop speed
        dsdpar_cexpo = disdrodata.cloud_phase_space_parameters(_speed_='expo')
        # expoential law for drop speed
        dsdpar_cplaw = disdrodata.cloud_phase_space_parameters(_speed_='plaw')

        # process cloudexpo representation
        (_pca_cplaw, _pca_compo_cplaw) = PCA_maker(dsdpar_cplaw, ndropspdf='yes', pdftype='cloudplaw')
        (_pca_cplaw_noN, _pca_compo_cplaw_noN) = PCA_maker(dsdpar_cplaw, ndropspdf='no', pdftype='cloudplaw')
        res_pca_cplaw = pd.concat([_pca_cplaw, _pca_cplaw_noN], ignore_index=True)
        res_pca_cplaw['pdftype'] = 'cloudplaw'
        res_pca_compo_cplaw = pd.concat([_pca_compo_cplaw, _pca_compo_cplaw_noN], ignore_index=True)
        res_pca_compo_cplaw['pdftype'] = 'cloudplaw'

        # process cloudplaw representation
        (_pca_cexpo, _pca_compo_cexpo) = PCA_maker(dsdpar_cexpo, ndropspdf='yes', pdftype='cloudexpo')
        (_pca_cexpo_noN, _pca_compo_cexpo_noN) = PCA_maker(dsdpar_cexpo, ndropspdf='no', pdftype='cloudexpo')
        res_pca_cexpo = pd.concat([_pca_cexpo, _pca_cexpo_noN], ignore_index=True)
        res_pca_cexpo['pdftype'] = 'cloudexpo'
        res_pca_compo_cexpo = pd.concat([_pca_compo_cexpo, _pca_compo_cexpo_noN], ignore_index=True)
        res_pca_compo_cexpo['pdftype'] = 'cloudexpo'

        res_pca = res_pca_flux.copy()
        res_pca = pd.concat([res_pca, res_pca_cplaw], ignore_index=True)
        res_pca = pd.concat([res_pca, res_pca_cexpo], ignore_index=True)
        res_pca['site'] = args._da_
        res_pca['instrument'] = instr

        res_pca_compo = res_pca_compo_flux.copy()
        res_pca_compo = pd.concat([res_pca_compo, res_pca_compo_cplaw], ignore_index=True)
        res_pca_compo = pd.concat([res_pca_compo, res_pca_compo_cexpo], ignore_index=True)
        res_pca_compo['site'] = args._da_
        res_pca_compo['instrument'] = instr
    else:
        dsdpar_cloud = disdrodata.cloud_phase_space_parameters()

        # process cloudplaw representation
        (_pca_cloud, _pca_compo_cloud) = PCA_maker(dsdpar_cloud, ndropspdf='yes', pdftype='cloud2dvd')
        (_pca_cloud_noN, _pca_compo_cloud_noN) = PCA_maker(dsdpar_cloud, ndropspdf='no', pdftype='cloud2dvd')
        res_pca_cloud = pd.concat([_pca_cloud, _pca_cloud_noN], ignore_index=True)
        res_pca_cloud['pdftype'] = 'cloud2dvd'
        res_pca_compo_cloud = pd.concat([_pca_compo_cloud, _pca_compo_cloud_noN], ignore_index=True)
        res_pca_compo_cloud['pdftype'] = 'cloud2dvd'

        res_pca = res_pca_flux.copy()
        res_pca = pd.concat([res_pca, res_pca_cloud], ignore_index=True)
        res_pca['site'] = args._da_
        res_pca['instrument'] = instr

        res_pca_compo = res_pca_compo_flux.copy()
        res_pca_compo = pd.concat([res_pca_compo, res_pca_compo_cloud], ignore_index=True)
        res_pca_compo['site'] = args._da_
        res_pca_compo['instrument'] = instr
else:
    print("Acronym not found in catalog")
    sys.exit()

res_pca.to_csv(args._pcave_, sep=' ', index=None)
res_pca_compo.to_csv(args._pcac_, sep=' ', index=None)
