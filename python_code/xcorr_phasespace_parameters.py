import sys
import pandas as pd 
import numpy as np
import argparse
# path to where package is located (change accordingly)
pack_path = '/home/massimiliano/work_vivo/pioggia/disdrorain'

if pack_path not in sys.path:
    sys.path.append(pack_path)

import disdrorain as dr


# Arguments
parser = argparse.ArgumentParser(description='', epilog="")
parser.add_argument('--disdrodatapath', action="store", dest='_ddp_', default='_NONE_')
parser.add_argument('--disdrocatalog', action="store", dest='_dc_', default='_NONE_')
parser.add_argument('--disdroacronym', action="store", dest='_da_', default='_NONE_')
parser.add_argument('--xcorroutmatrix', action="store", dest='_dom_', default='_NONE_')
parser.add_argument('--xcorroutplot', action="store", dest='_dop_', default='_NONE_')
args = parser.parse_args()

#
print("Executing: ", sys.argv[0])
print()
#PI = 3.141592653589793  # approximate value of greek pi
#seconds_in_hour = 3600  # for converting rainfall rate in mm/h

catalog=pd.read_csv(args._dc_,sep=',')

index_acronym=catalog[catalog['ACRONYM']==args._da_].index.tolist()
if index_acronym:
	row = catalog.loc[catalog.ACRONYM==args._da_,:]
	celllimits_path = args._ddp_+'/'+row['CELLLIMITS'].values[0]
	areainstr = row['AREA_INSTRUMENT'].values[0]
	instr = row['INSTRUMENT'].values[0]
	tr = row['TIME_RESOLUTION'].values[0]
	pathtodata=args._ddp_+'/'+args._da_

	if (row['CELLLIMITS'].values[0]!='NONE'):
		if row['CELLLIMITS'].values[0]!='standard':
			disdrodata=dr.disdrorain(classpath=celllimits_path,datapath=pathtodata,instrument_area=areainstr,time_interval=tr)
		if row['CELLLIMITS'].values[0]=='standard':
			disdrodata=dr.disdrorain(datapath=pathtodata,instrument_area=areainstr,time_interval=tr)
		psp_flux = disdrodata.flux_phase_space_parameters()
		psp_cloud_plaw = disdrodata.cloud_phase_space_parameters(_speed_='plaw')
		psp_cloud_expo = disdrodata.cloud_phase_space_parameters(_speed_='expo')

		xcorr_psp_flux=psp_flux.corr().reset_index()
		xcorr_psp_flux.rename(columns={'index' : 'Variable', 'N' :'N/Nv'},inplace=True)
		xcorr_psp_flux['pdftype']='flux'

		xcorr_psp_cloud_plaw=psp_cloud_plaw.corr().reset_index()
		xcorr_psp_cloud_plaw.rename(columns={'index' : 'Variable', 'Nv' :'N/Nv'},inplace=True)
		xcorr_psp_cloud_plaw['pdftype']='cloud_plaw'

		xcorr_psp_cloud_expo=psp_cloud_expo.corr().reset_index()
		xcorr_psp_cloud_expo.rename(columns={'index' : 'Variable', 'Nv' :'N/Nv'},inplace=True)
		xcorr_psp_cloud_expo['pdftype']='cloud_expo'

		xcorr_psp_summary=xcorr_psp_flux.copy()
		xcorr_psp_summary=xcorr_psp_summary.append(xcorr_psp_cloud_plaw)
		xcorr_psp_summary=xcorr_psp_summary.append(xcorr_psp_cloud_expo)
		xcorr_psp_summary['site'] = args._da_
		xcorr_psp_summary['instrument'] =instr


		row_xcorr_psp_summary_xplot = pd.DataFrame(columns=['N/Nv_mu', 'N/Nv_sigma', 'N/Nv_gamma', 'N/Nv_kappa', 'N/Nv_eta', 'mu_sigma','mu_gamma','mu_kappa','mu_eta',
			'sigma_gamma','sigma_kappa','sigma_eta','gamma_kappa', 'gamma_eta','kappa_eta','pdftype','acronym'])
		row_xcorr_psp_summary_xplot = row_xcorr_psp_summary_xplot.astype(dtype = {'N/Nv_mu': 'float64', 'N/Nv_sigma': 'float64', 'N/Nv_gamma': 'float64', 
			'N/Nv_kappa': 'float64', 'N/Nv_eta': 'float64', 'mu_sigma': 'float64','mu_gamma': 'float64', 'mu_kappa': 'float64', 'mu_eta': 'float64',
			'sigma_gamma': 'float64', 'sigma_kappa': 'float64', 'sigma_eta': 'float64', 'gamma_kappa': 'float64', 'gamma_eta': 'float64',
			'kappa_eta': 'float64', 'pdftype' : 'object', 'acronym' : 'object'})
		xcorr_psp_summary_xplot = row_xcorr_psp_summary_xplot.copy() 



		row_xcorr_psp_summary_xplot['N/Nv_mu'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='flux') & (xcorr_psp_summary.Variable=='N'),'mu'].values
		row_xcorr_psp_summary_xplot['N/Nv_sigma'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='flux') & (xcorr_psp_summary.Variable=='N'),'sigma'].values
		row_xcorr_psp_summary_xplot['N/Nv_gamma'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='flux') & (xcorr_psp_summary.Variable=='N'),'gamma'].values
		row_xcorr_psp_summary_xplot['N/Nv_kappa'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='flux') & (xcorr_psp_summary.Variable=='N'),'kappa'].values
		row_xcorr_psp_summary_xplot['N/Nv_eta'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='flux') & (xcorr_psp_summary.Variable=='N'),'eta'].values
		row_xcorr_psp_summary_xplot['mu_sigma'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='flux') & (xcorr_psp_summary.Variable=='mu'),'sigma'].values
		row_xcorr_psp_summary_xplot['mu_gamma'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='flux') & (xcorr_psp_summary.Variable=='mu'),'gamma'].values
		row_xcorr_psp_summary_xplot['mu_kappa'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='flux') & (xcorr_psp_summary.Variable=='mu'),'kappa'].values
		row_xcorr_psp_summary_xplot['mu_eta'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='flux') & (xcorr_psp_summary.Variable=='mu'),'eta'].values
		row_xcorr_psp_summary_xplot['sigma_gamma'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='flux') & (xcorr_psp_summary.Variable=='sigma'),'gamma'].values
		row_xcorr_psp_summary_xplot['sigma_kappa'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='flux') & (xcorr_psp_summary.Variable=='sigma'),'kappa'].values
		row_xcorr_psp_summary_xplot['sigma_eta'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='flux') & (xcorr_psp_summary.Variable=='sigma'),'eta'].values
		row_xcorr_psp_summary_xplot['gamma_kappa'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='flux') & (xcorr_psp_summary.Variable=='gamma'),'kappa'].values
		row_xcorr_psp_summary_xplot['gamma_eta'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='flux') & (xcorr_psp_summary.Variable=='gamma'),'eta'].values
		row_xcorr_psp_summary_xplot['kappa_eta'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='flux') & (xcorr_psp_summary.Variable=='kappa'),'eta'].values
		row_xcorr_psp_summary_xplot['pdftype'] = 'flux'
		row_xcorr_psp_summary_xplot['acronym'] = args._da_
		xcorr_psp_summary_xplot = xcorr_psp_summary_xplot.append(row_xcorr_psp_summary_xplot)

		row_xcorr_psp_summary_xplot['N/Nv_mu'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='cloud_plaw') & (xcorr_psp_summary.Variable=='Nv'),'mu'].values
		row_xcorr_psp_summary_xplot['N/Nv_sigma'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='cloud_plaw') & (xcorr_psp_summary.Variable=='Nv'),'sigma'].values
		row_xcorr_psp_summary_xplot['N/Nv_gamma'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='cloud_plaw') & (xcorr_psp_summary.Variable=='Nv'),'gamma'].values
		row_xcorr_psp_summary_xplot['N/Nv_kappa'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='cloud_plaw') & (xcorr_psp_summary.Variable=='Nv'),'kappa'].values
		row_xcorr_psp_summary_xplot['N/Nv_eta'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='cloud_plaw') & (xcorr_psp_summary.Variable=='Nv'),'eta'].values
		row_xcorr_psp_summary_xplot['mu_sigma'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='cloud_plaw') & (xcorr_psp_summary.Variable=='mu'),'sigma'].values
		row_xcorr_psp_summary_xplot['mu_gamma'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='cloud_plaw') & (xcorr_psp_summary.Variable=='mu'),'gamma'].values
		row_xcorr_psp_summary_xplot['mu_kappa'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='cloud_plaw') & (xcorr_psp_summary.Variable=='mu'),'kappa'].values
		row_xcorr_psp_summary_xplot['mu_eta'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='cloud_plaw') & (xcorr_psp_summary.Variable=='mu'),'eta'].values
		row_xcorr_psp_summary_xplot['sigma_gamma'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='cloud_plaw') & (xcorr_psp_summary.Variable=='sigma'),'gamma'].values
		row_xcorr_psp_summary_xplot['sigma_kappa'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='cloud_plaw') & (xcorr_psp_summary.Variable=='sigma'),'kappa'].values
		row_xcorr_psp_summary_xplot['sigma_eta'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='cloud_plaw') & (xcorr_psp_summary.Variable=='sigma'),'eta'].values
		row_xcorr_psp_summary_xplot['gamma_kappa'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='cloud_plaw') & (xcorr_psp_summary.Variable=='gamma'),'kappa'].values
		row_xcorr_psp_summary_xplot['gamma_eta'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='cloud_plaw') & (xcorr_psp_summary.Variable=='gamma'),'eta'].values
		row_xcorr_psp_summary_xplot['kappa_eta'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='cloud_plaw') & (xcorr_psp_summary.Variable=='kappa'),'eta'].values
		row_xcorr_psp_summary_xplot['pdftype'] = 'cloud_plaw'
		row_xcorr_psp_summary_xplot['acronym'] = args._da_
		xcorr_psp_summary_xplot = xcorr_psp_summary_xplot.append(row_xcorr_psp_summary_xplot)

		row_xcorr_psp_summary_xplot['N/Nv_mu'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='cloud_expo') & (xcorr_psp_summary.Variable=='Nv'),'mu'].values
		row_xcorr_psp_summary_xplot['N/Nv_sigma'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='cloud_expo') & (xcorr_psp_summary.Variable=='Nv'),'sigma'].values
		row_xcorr_psp_summary_xplot['N/Nv_gamma'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='cloud_expo') & (xcorr_psp_summary.Variable=='Nv'),'gamma'].values
		row_xcorr_psp_summary_xplot['N/Nv_kappa'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='cloud_expo') & (xcorr_psp_summary.Variable=='Nv'),'kappa'].values
		row_xcorr_psp_summary_xplot['N/Nv_eta'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='cloud_expo') & (xcorr_psp_summary.Variable=='Nv'),'eta'].values
		row_xcorr_psp_summary_xplot['mu_sigma'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='cloud_expo') & (xcorr_psp_summary.Variable=='mu'),'sigma'].values
		row_xcorr_psp_summary_xplot['mu_gamma'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='cloud_expo') & (xcorr_psp_summary.Variable=='mu'),'gamma'].values
		row_xcorr_psp_summary_xplot['mu_kappa'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='cloud_expo') & (xcorr_psp_summary.Variable=='mu'),'kappa'].values
		row_xcorr_psp_summary_xplot['mu_eta'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='cloud_expo') & (xcorr_psp_summary.Variable=='mu'),'eta'].values
		row_xcorr_psp_summary_xplot['sigma_gamma'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='cloud_expo') & (xcorr_psp_summary.Variable=='sigma'),'gamma'].values
		row_xcorr_psp_summary_xplot['sigma_kappa'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='cloud_expo') & (xcorr_psp_summary.Variable=='sigma'),'kappa'].values
		row_xcorr_psp_summary_xplot['sigma_eta'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='cloud_expo') & (xcorr_psp_summary.Variable=='sigma'),'eta'].values
		row_xcorr_psp_summary_xplot['gamma_kappa'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='cloud_expo') & (xcorr_psp_summary.Variable=='gamma'),'kappa'].values
		row_xcorr_psp_summary_xplot['gamma_eta'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='cloud_expo') & (xcorr_psp_summary.Variable=='gamma'),'eta'].values
		row_xcorr_psp_summary_xplot['kappa_eta'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='cloud_expo') & (xcorr_psp_summary.Variable=='kappa'),'eta'].values
		row_xcorr_psp_summary_xplot['pdftype'] = 'cloud_expo'
		row_xcorr_psp_summary_xplot['acronym'] = args._da_
		xcorr_psp_summary_xplot = xcorr_psp_summary_xplot.append(row_xcorr_psp_summary_xplot)
		xcorr_psp_summary_xplot['instrument'] = instr
		xcorr_psp_summary_xplot['isparsi'] = args._da_.find('parsi')

	else:
		disdrodata=dr.disdrorain_2dvd(datapath=pathtodata,instrument_area=areainstr,time_interval=tr)

		psp_flux = disdrodata.flux_phase_space_parameters()
		psp_cloud = disdrodata.cloud_phase_space_parameters()

		xcorr_psp_flux=psp_flux.corr().reset_index()
		xcorr_psp_flux.rename(columns={'index' : 'Variable', 'N' :'N/Nv'},inplace=True)
		xcorr_psp_flux['pdftype']='flux'

		xcorr_psp_cloud =psp_cloud.corr().reset_index()
		xcorr_psp_cloud.rename(columns={'index' : 'Variable', 'Nv' :'N/Nv'},inplace=True)
		xcorr_psp_cloud['pdftype']='cloud'

		xcorr_psp_summary=xcorr_psp_flux.copy()
		xcorr_psp_summary=xcorr_psp_summary.append(xcorr_psp_cloud)
		xcorr_psp_summary['site'] = args._da_
		xcorr_psp_summary['instrument'] =instr


		row_xcorr_psp_summary_xplot = pd.DataFrame(columns=['N/Nv_mu', 'N/Nv_sigma', 'N/Nv_gamma', 'N/Nv_kappa', 'N/Nv_eta', 'mu_sigma','mu_gamma','mu_kappa','mu_eta',
			'sigma_gamma','sigma_kappa','sigma_eta','gamma_kappa', 'gamma_eta','kappa_eta','pdftype','acronym'])
		row_xcorr_psp_summary_xplot = row_xcorr_psp_summary_xplot.astype(dtype = {'N/Nv_mu': 'float64', 'N/Nv_sigma': 'float64', 'N/Nv_gamma': 'float64', 
			'N/Nv_kappa': 'float64', 'N/Nv_eta': 'float64', 'mu_sigma': 'float64','mu_gamma': 'float64', 'mu_kappa': 'float64', 'mu_eta': 'float64',
			'sigma_gamma': 'float64', 'sigma_kappa': 'float64', 'sigma_eta': 'float64', 'gamma_kappa': 'float64', 'gamma_eta': 'float64',
			'kappa_eta': 'float64', 'pdftype' : 'object', 'acronym' : 'object'})
		xcorr_psp_summary_xplot = row_xcorr_psp_summary_xplot.copy() 



		row_xcorr_psp_summary_xplot['N/Nv_mu'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='flux') & (xcorr_psp_summary.Variable=='N'),'mu'].values
		row_xcorr_psp_summary_xplot['N/Nv_sigma'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='flux') & (xcorr_psp_summary.Variable=='N'),'sigma'].values
		row_xcorr_psp_summary_xplot['N/Nv_gamma'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='flux') & (xcorr_psp_summary.Variable=='N'),'gamma'].values
		row_xcorr_psp_summary_xplot['N/Nv_kappa'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='flux') & (xcorr_psp_summary.Variable=='N'),'kappa'].values
		row_xcorr_psp_summary_xplot['N/Nv_eta'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='flux') & (xcorr_psp_summary.Variable=='N'),'eta'].values
		row_xcorr_psp_summary_xplot['mu_sigma'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='flux') & (xcorr_psp_summary.Variable=='mu'),'sigma'].values
		row_xcorr_psp_summary_xplot['mu_gamma'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='flux') & (xcorr_psp_summary.Variable=='mu'),'gamma'].values
		row_xcorr_psp_summary_xplot['mu_kappa'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='flux') & (xcorr_psp_summary.Variable=='mu'),'kappa'].values
		row_xcorr_psp_summary_xplot['mu_eta'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='flux') & (xcorr_psp_summary.Variable=='mu'),'eta'].values
		row_xcorr_psp_summary_xplot['sigma_gamma'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='flux') & (xcorr_psp_summary.Variable=='sigma'),'gamma'].values
		row_xcorr_psp_summary_xplot['sigma_kappa'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='flux') & (xcorr_psp_summary.Variable=='sigma'),'kappa'].values
		row_xcorr_psp_summary_xplot['sigma_eta'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='flux') & (xcorr_psp_summary.Variable=='sigma'),'eta'].values
		row_xcorr_psp_summary_xplot['gamma_kappa'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='flux') & (xcorr_psp_summary.Variable=='gamma'),'kappa'].values
		row_xcorr_psp_summary_xplot['gamma_eta'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='flux') & (xcorr_psp_summary.Variable=='gamma'),'eta'].values
		row_xcorr_psp_summary_xplot['kappa_eta'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='flux') & (xcorr_psp_summary.Variable=='kappa'),'eta'].values
		row_xcorr_psp_summary_xplot['pdftype'] = 'flux'
		row_xcorr_psp_summary_xplot['acronym'] = args._da_
		xcorr_psp_summary_xplot = xcorr_psp_summary_xplot.append(row_xcorr_psp_summary_xplot)

		row_xcorr_psp_summary_xplot['N/Nv_mu'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='cloud') & (xcorr_psp_summary.Variable=='Nv'),'mu'].values
		row_xcorr_psp_summary_xplot['N/Nv_sigma'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='cloud') & (xcorr_psp_summary.Variable=='Nv'),'sigma'].values
		row_xcorr_psp_summary_xplot['N/Nv_gamma'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='cloud') & (xcorr_psp_summary.Variable=='Nv'),'gamma'].values
		row_xcorr_psp_summary_xplot['N/Nv_kappa'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='cloud') & (xcorr_psp_summary.Variable=='Nv'),'kappa'].values
		row_xcorr_psp_summary_xplot['N/Nv_eta'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='cloud') & (xcorr_psp_summary.Variable=='Nv'),'eta'].values
		row_xcorr_psp_summary_xplot['mu_sigma'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='cloud') & (xcorr_psp_summary.Variable=='mu'),'sigma'].values
		row_xcorr_psp_summary_xplot['mu_gamma'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='cloud') & (xcorr_psp_summary.Variable=='mu'),'gamma'].values
		row_xcorr_psp_summary_xplot['mu_kappa'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='cloud') & (xcorr_psp_summary.Variable=='mu'),'kappa'].values
		row_xcorr_psp_summary_xplot['mu_eta'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='cloud') & (xcorr_psp_summary.Variable=='mu'),'eta'].values
		row_xcorr_psp_summary_xplot['sigma_gamma'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='cloud') & (xcorr_psp_summary.Variable=='sigma'),'gamma'].values
		row_xcorr_psp_summary_xplot['sigma_kappa'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='cloud') & (xcorr_psp_summary.Variable=='sigma'),'kappa'].values
		row_xcorr_psp_summary_xplot['sigma_eta'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='cloud') & (xcorr_psp_summary.Variable=='sigma'),'eta'].values
		row_xcorr_psp_summary_xplot['gamma_kappa'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='cloud') & (xcorr_psp_summary.Variable=='gamma'),'kappa'].values
		row_xcorr_psp_summary_xplot['gamma_eta'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='cloud') & (xcorr_psp_summary.Variable=='gamma'),'eta'].values
		row_xcorr_psp_summary_xplot['kappa_eta'] = xcorr_psp_summary.loc[(xcorr_psp_summary.pdftype=='cloud') & (xcorr_psp_summary.Variable=='kappa'),'eta'].values
		row_xcorr_psp_summary_xplot['pdftype'] = 'cloud'
		row_xcorr_psp_summary_xplot['acronym'] = args._da_
		xcorr_psp_summary_xplot = xcorr_psp_summary_xplot.append(row_xcorr_psp_summary_xplot)
		xcorr_psp_summary_xplot['instrument'] = instr
		xcorr_psp_summary_xplot['isparsi'] = args._da_.find('parsi')
else:
	print("Acronym not found in catalog")
	sys.exit()

xcorr_psp_summary.to_csv(args._dom_,sep=' ',index=None)
xcorr_psp_summary_xplot.to_csv(args._dop_,sep=' ',index=None)
