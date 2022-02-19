# This program calculates the drop size distribution moments: 1st,2nd, 3rd, 4th, 5th, 6th (change values of list_moments_order if necessary)
# 3 possible representations of the drop size distribution area available
# 	1) "flux" or ground representation: couple (N,p(D))
# 	2) "cloudexpo" representation: couple (N_{V},f(D)). If speed of drop is not known than v(D)=9.65-10.3\times e^{-0.6D} is assumed
# 	3) "cloudplaw" representation: couple (N_{V},f(D)). If speed of drop is not known than v(D)=3.78\times D^{0.67} is assumed

# ------------- Necessary Python packages -START
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
# --disdrodatapath	:	path to directory contain the disdrometer data
# --disdrocatalog	:	full path to file with data catalog (data catalog contain metdata about the dataset)
# --disdroacronym	:	acronym identifying the data set in the catalog
# --representation	:	which type of drop size representation to adopt: possible values are "flux", "cloudexpo", "cloudplaw
#					:	default is "flux"
# --output			:	output file name for the drop size distribution parameters
#						default value is "dsd_moments"
#						default output is generated in the directory from which code is invoked
parser = argparse.ArgumentParser(description='', epilog="")
parser.add_argument('--disdrodatapath', action="store", dest='_ddp_', default='_NONE_')
parser.add_argument('--disdrocatalog', action="store", dest='_dc_', default='_NONE_')
parser.add_argument('--disdroacronym', action="store", dest='_da_', default='_NONE_')
parser.add_argument('--representation', action="store", dest='_re_', default='flux')
parser.add_argument('--output', action="store", dest='_ou_', default='dsd_moments')
args = parser.parse_args()


#
print("Executing: ", sys.argv[0])
print()

# load catalog
catalog = pd.read_csv(args._dc_, sep=',')
# retrieve index of catalog data frame corresponding to the desired acronym
index_acronym = catalog[catalog['ID2'] == args._da_].index.tolist()

# check if to calculate also renormalized values

# check if representation value is admitted (only: flux, cloudexpo, cloudplaw)
if args._re_ not in ['flux', 'cloudexpo', 'cloudplaw']:
	print("ERROR: representation value not found!")
	print("possible choices are flux, cloudexpo, cloudplaw")
	sys.exit()

# if acronym is part of the catalog we calculate the phase space parameters 
if index_acronym:
	row = catalog.loc[catalog.ID2 == args._da_,:]  # retrieve row corresponding to acronym
	celllimits_path = args._ddp_+'/'+row['CELLLIMITS'].values[0]  # set file name with cell limits values 
	areainstr = row['AREA_INSTRUMENT'].values[0]  # set area of instrument 
	instr = row['INSTRUMENT'].values[0]  # set instrument type
	tr = row['TIME_RESOLUTION'].values[0]  # set instrument time resolution
	pathtodata = args._ddp_+'/'+args._da_  # set path to data set

	# if data set is not from 2DVD disdrometer
	if (instr != '2DVD'):
		# if data set is RD80 with standard cell limits division
		if ((instr == 'RD80') & (row['CELLLIMITS'].values[0] == 'standard')):
			disdrodata = dr.disdrorain(datapath=pathtodata,instrument_area=areainstr,time_interval=tr)  # create disdrodata class 
		else:
			disdrodata = dr.disdrorain(classpath=celllimits_path,datapath=pathtodata,instrument_area=areainstr,time_interval=tr)  # create disdrodata class  
	else:
		disdrodata = dr.disdrorain_2dvd(datapath=pathtodata,instrument_area=areainstr,time_interval=tr)  # create disdrodata class for 2dvd data

	# calculate dsd paramters
	list_moments_order = list([1, 2, 3, 4, 5, 6])
	if args._re_ == 'flux':
		# calculate central moments in the flux (ground) representation
		dsdpar = disdrodata.flux_moment_calculator(list_moments_order)
	else:
		# calculate phase space paramter in the cloud representation
		if args._re_ == 'cloudexpo':
			# exponential law for drop speed
			dsdpar = disdrodata.cloud_moment_calculator(list_moments_order,_speed_='expo')
		if args._re_ == 'cloudplaw':
			# expoential law for drop speed
			dsdpar = disdrodata.cloud_moment_calculator(list_moments_order,_speed_='plaw')
else:
	print("ERROR: ID2 not found in catalog")

# output results to file
dsdpar.to_csv(args._ou_,sep=' ',index=None)
