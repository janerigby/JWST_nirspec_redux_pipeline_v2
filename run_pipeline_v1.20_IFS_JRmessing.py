# from Brian Welch, reducing jwst nirspec ifs spectra using v1.20 of the jwst pipeline
# this is all copied over from my notebook nirspec_pipeline_rough
# I put it in this script so I can easily run many iterations at once
import numpy as np
import glob
import os

###############################################
# Set up CRDS path and server environment variables
home = "/Users/bwelch/Documents/"
os.environ["CRDS_PATH"] = home + "crds_cache/jwst_ops"
os.environ["CRDS_SERVER_URL"] = "https://jwst-crds.stsci.edu"
#os.environ["CRDS_CONTEXT"] = "jwst_1536.pmap" # probably unnecessary, but included just to be sure
os.environ["CRDS_CONTEXT"] = "jwst_1466.pmap"


import zipfile
import urllib.request

import json
import asdf

from astropy.io import fits
from astropy.utils.data import download_file
import astropy.units as u
from astropy import wcs
from astropy.wcs import WCS
from astropy.visualization import ImageNormalize, ManualInterval, LogStretch, LinearStretch, AsinhStretch

import matplotlib.pyplot as plt
import matplotlib as mpl

# The calwebb_spec and spec3 pipelines
from jwst.pipeline import Spec2Pipeline
from jwst.pipeline import Spec3Pipeline

# the level1 pipeline:
from jwst.pipeline import Detector1Pipeline

# data models
from jwst import datamodels

# association file utilities
from jwst.associations import asn_from_list as afl # Tools for creating association files
from jwst.associations.lib.rules_level2_base import DMSLevel2bBase # Definition of a Lvl2 association file
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base # Definition of a Lvl3 association file

#from stcal import dqflags

from scipy.interpolate import interp1d

targetlist = ['SGAS1110','SGAS2111','SGAS1050','SGAS1527','SGAS1429','COSMIC.EYE',
				'SGAS1723','SGAS1226','SPT0418-47','SPT2147-50', 'SGAS1402m28']

def run_pipeline(target, rundet1=True, runspec2=True, runspec3=True, runspec3_PRISM=False, clean_data=True):
	'''
	runspec3_PRISM == use non-linear wavelength solution for prism observations
		only runs the prism data through spec3 pipeline (things break trying to do all data here)
	clean_data = whether or not to run the "cleaning" step before stage 3 (turn off to save time when only running spec3)
	'''
	print(f'Starting NIRSpec IFU pipeline for {target}')
	# Modify the paths to the relevant directories on your machine
	# ------------------------------------------------------------
	# 1) point to where the jwst pipeline config files are located
	home = "/Users/bwelch/Documents/" # for B. Welch
	#home = "/Users/tahutch1/programs/jwst-drp/" # for T. Hutchison

	# 2) point to where you keep your uncal data
	# for B. Welch
	# LEGGOS targets:
	if target == 'SGAS1110':
		input_path = '/Users/bwelch/Documents/data/leggos/sgas1110/nirspec/full_uncal/'
	if target == 'SGAS2111':
		input_path = '/Users/bwelch/Documents/data/leggos/sgas2111/nirspec/full_uncal/' 
	if target == 'SGAS1050':
		input_path = '/Users/bwelch/Documents/data/leggos/sgas1050/nirspec/full_uncal/'
	if target == 'SGAS1527':
		input_path = '/Users/bwelch/Documents/data/leggos/sgas1527/nirspec/full_uncal' 
	if target == 'SGAS1429':
		input_path = '/Users/bwelch/Documents/data/leggos/sgas1429/nirspec/full_uncal/'
	if target == 'COSMIC.EYE':
		input_path = '/Users/bwelch/Documents/data/leggos/cosmic_eye/nirspec/full_uncal/'
					# TEMPLATES targets:
	if target == 'SGAS1723':
		input_path = '/Users/bwelch/Documents/data/templates/sdss1723/full_uncal_data/'
	if target == 'SGAS1226':
		input_path = '/Users/bwelch/Documents/data/templates/sdss1226/nirspec/MAST_2023-01-03T0916/JWST/' 
	if target == 'SPT0418-47':
		input_path = '/Users/bwelch/Documents/data/templates/spt0418/MAST_2022-10-10T1412/JWST/'
	if target == 'SPT2147-50':
		input_path = '/Users/bwelch/Documents/data/templates/spt2147/MAST_2022-11-02T1127/JWST/'
					# other stuff
	if target == 'SPARKLER':
		input_path = '/Users/bwelch/Documents/data/sparkler/nirspec/MAST_2024-05-30T1445/JWST/' 
	if target == 'SUNBURST':
		input_path = '/Users/bwelch/Documents/data/sgas/sunburst/nirspec/MAST_2023-04-11T0820/JWST/'
	if target == 'J1757132':
		input_path = '/Users/bwelch/Documents/data/NIRSpec_PSF/full_uncal'
	if target == 'SGAS1402m28':
		input_path = '/Users/jrrigby1/SCIENCE/JWST_Data/LEGGOS/SGAS1402m28/Raw_data/'

	# 3) point to where you want your processed outputs to live
	# for B. Welch
	# LEGGOS targets: 
	if target == 'SGAS1110':
		output_path = '/Users/bwelch/Documents/data/leggos/sgas1110/nirspec/pmap1466/'
	if target == 'SGAS2111':
		output_path = '/Users/bwelch/Documents/data/leggos/sgas2111/nirspec/pmap1466/'
	if target == 'SGAS1050':
		output_path = '/Users/bwelch/Documents/data/leggos/sgas1050/nirspec/pmap1466/'
	if target == 'SGAS1527':
		output_path = '/Users/bwelch/Documents/data/leggos/sgas1527/nirspec/pmap1466/'
	if target == 'SGAS1429':
		output_path = '/Users/bwelch/Documents/data/leggos/sgas1429/nirspec/pmap1466/'
	if target == 'COSMIC.EYE':
		output_path = '/Users/bwelch/Documents/data/leggos/cosmic_eye/nirspec/pmap1466/'
					# TEMPLATES targets: 
	if target == 'SGAS1723':
		output_path = '/Users/bwelch/Documents/data/templates/sdss1723/nirspec/pmap1466/'
	if target == 'SGAS1226':
		output_path = '/Users/bwelch/Documents/data/templates/sdss1226/nirspec/pmap1466/'
	if target == 'SPT0418-47':
		output_path = '/Users/bwelch/Documents/data/templates/spt0418/nirspec/pmap1466/'
	if target == 'SPT2147-50':
		output_path = '/Users/bwelch/Documents/data/templates/spt2147/nirspec/pmap1466/'
					# Other stuff
	if target == 'SPARKLER':
		output_path = '/Users/bwelch/Documents/data/sparkler/nirspec/pmap1466/'
	if target == 'SUNBURST':
		output_path = '/Users/bwelch/Documents/data/sgas/sunburst/nirspec/pmap1466/'
	if target == 'J1757132':
		output_path = '/Users/bwelch/Documents/data/NIRSpec_PSF/pmap1466/'
	if target == 'SGAS1402m28':
		output_path = '/Users/jrrigby1/SCIENCE/JWST_Data/LEGGOS/SGAS1402m28/Redux_v1.20.2/'


	if os.path.exists(output_path) == False: # if folder doesn't exist
		print('Creating folder ' + output_path)
		os.system('mkdir ' + output_path) # creates the folder


	# test by TAH to see if the files are accessible

	files = glob.glob(input_path + '**/*_uncal.fits', recursive=True) #list the uncalibrated (level 1b) files.
	files = sorted(files)

	ifu,sky = 0,0

	for exposure in files:
		test = fits.open(exposure)
		head = test[0].header
		if head['BKGDTARG'] == False: ifu += 1
		elif head['BKGDTARG'] == True: sky += 1

	print(f'Number of IFU: {ifu}, Number of sky: {sky}')


	# checking that the file system is in place for these data
	# if not, creating the folders
	folders_L2a = ['L2a/','L2a/sci/','L2a/bkg/']

	for folder in folders_L2a:
		if os.path.exists(output_path + folder) == False: # if folder doesn't exist
			print('Creating folder ' + output_path + folder)
			os.system('mkdir ' + output_path + folder) # creates the folder


	# Run the pipeline, splitting outputs into "sci" and "bkg" output folders
	# Leakcals will be in the same folders as their associated observations
	if rundet1 == True:
		print(f'Running Detector1 Pipeline for {target}')
		for i,exposure in enumerate(files): 
			print(f'Running exposure {i} / {len(files)}', end='\r')
			det1 = Detector1Pipeline()
			# set output directory based on sci vs bkg exposure
			head = fits.open(exposure)[0].header
			if head['BKGDTARG'] == False:
				det1.output_dir = output_path + folders_L2a[1]
			elif head['BKGDTARG'] == True:
				det1.output_dir = output_path + folders_L2a[2]
			else: print('not target')
			if head['IS_IMPRT'] == True:
				print('Skipping leak cal file: ',exposure)
				continue
			# now we set other parameters:
		det1.save_results = True
		det1.jump.maximum_cores = 'half'
		det1.jump.expand_large_events = True
		det1.ramp_fit.maximum_cores = 'half'
		# NSClean step moved to Detector1 pipeline
		# throws an error if there's nothing on NRS2, so
		if 'M' in head['GRATING']:
			if 'nrs2' in exposure:
				det1.clean_flicker_noise.skip = True
		else:
			det1.clean_flicker_noise.skip = False
			det1.clean_flicker_noise.fit_method = 'fft'
			det1.clean_flicker_noise.mask_science_regions = True

			det1.run(exposure)

		print('Done Detector1 for {target}')

	# checking that the file system is in place for these data
	# if not, creating the folders
	folders_L2b = ['L2b/','L2b/sci/','L2b/bkg/']

	for folder in folders_L2b:
		if os.path.exists(output_path + folder) == False: # if folder doesn't exist
			print('Creating folder ' + output_path + folder)
			os.system('mkdir ' + output_path + folder) # creates the folder

	# define some path names for easy access
	use_leak = False
	process_background = False

	level1_sci_dir = os.path.join(output_path, folders_L2a[1]) 
	level2outputdir = os.path.join(output_path, folders_L2b[0])

	ratefiles = sorted(glob.glob(os.path.join(level1_sci_dir,'*rate.fits')))

	# And now we run the pipeline with the asn files we just made! 
	if runspec2 == True:
		print(f'Starting Spec2 Pipeline for {target}')
		# this first part is kinda not needed, but its still here in case we want to reduce the background files at some point
		if (use_leak == True) or (process_background == True):
			if target == 'SGAS1723':
				asnfiles_sci = glob.glob(level2outputdir+'jw'+program_id+'-o'+obs_num_sci[0]+'*asn.json') 
				asnfiles_sci += glob.glob(level2outputdir+'jw'+program_id+'-o'+obs_num_sci[1]+'*asn.json')
			else:
				asnfiles_sci = glob.glob(level2outputdir+'jw'+'*asn.json')
			
			for asn in asnfiles_sci:
				#print(asn)
				spec2 = Spec2Pipeline()
				spec2.output_dir = output_path + folders_L2b[1]
				spec2.save_results = True
				spec2.pixel_replace.skip = False
				spec2.pixel_replace.algorithm = 'mingrad'
				# NSClean in Spec2 is deprecated, moved to Detector1 pipeline 
				#spec2.nsclean.skip = False
				#spec2.nsclean.save_mask = True
				spec2.run(asn)
			
			if process_background == True:
				asnfiles_bkg = glob.glob(level2outputdir+'jw'+program_id+'-o'+obs_num_bkg+'*asn.json')
				for asn in asnfiles_bkg:
					spec2 = Spec2Pipeline()
					spec2.output_dir = output_path + folders_L2b[2]
					spec2.save_results = True
					spec2.run(asn)

		else:
			for i,file in enumerate(ratefiles):
				print(f'Running file {i} / {len(ratefiles)}',end='\r')
				head = fits.open(file)[0].header
				if 'M' in head['GRATING']: 
					if 'nrs2' in file: continue
					#print(file)
					spec2 = Spec2Pipeline()
					spec2.output_dir = output_path + folders_L2b[1]
					spec2.save_results = True
					# skip cube building and 1d extraction, since we don't really use those anyway
					spec2.cube_build.skip = True 
					spec2.extract_1d.skip = True
					spec2.pixel_replace.skip = False
					spec2.pixel_replace.algorithm = 'mingrad'
					# NSClean in Spec2 is deprecated, moved to Detector1 pipeline 
					#spec2.nsclean.skip = False
					#spec2.nsclean.save_mask = True
					spec2.run(file)


		print(f'Done Spec2 Pipeline for {target}')


	# checking that the file system is in place for these data
	# if not, creating the folders
	folders_L3 = ['L3/']#, 

	for folder in folders_L3:
		if os.path.exists(output_path + folder) == False: # if folder doesn't exist
			print('Creating folder ' + output_path + folder)
			os.system('mkdir ' + output_path + folder) # creates the folder


	#orig_calfiles = glob.glob(level2outputdir + '/sci/*cal.fits') # cal files output from L2 pipeline

	# set different thresholds for each target - done by-eye w/ L2 ratefiles
	if target == 'SGAS1723': 
		calmax_blue = 2000
		calmax_red = 350
		calmax = 2000
	if target == 'SGAS1110': calmax = 500 # for g140h grating
	if target == 'SGAS2111': calmax = 50
	if target == 'SGAS1050': calmax = 1000 # edited up to 1000 from 500 on July 14 2026
	if target == 'SGAS1429': calmax = 400
	if target == 'SGAS1527': calmax = 300 # edited up to 300 from 150 on July 14 2026
	if target == 'COSMIC.EYE': calmax = 200
	if target == 'SGAS1226': calmax = 50
	if target == 'SPT0418-47': calmax = 40
	if target == 'SPT2147-50': calmax = 30
	if target == 'SPARKLER': calmax = 15

	calfiles = glob.glob(level2outputdir + '/sci/*cal.fits') # cal files output from L2 pipeline

	#for file in calfiles:
	#	outfile = file[:-5] + '2.fits'
	#	if target == 'SUNBURST':
	#		if ('_02101_' in file) or ('_10101_' in file):
	#			cut_cal(file, outfile, max_threshold=3000, min_threshold=-10)
	#		elif ('_04101_' in file) or ('_12101_' in file):
	#			cut_cal(file, outfile, max_threshold=1800, min_threshold=-10)
	#		elif ('_06101_' in file) or ('_14101_' in file):
	#			cut_cal(file, outfile, max_threshold=2300, min_threshold=-10)
	#		else:
	#			cut_cal(file, outfile, max_threshold=2000, min_threshold=-10)
	#	else:
	#		cut_cal(file, outfile, max_threshold=calmax, min_threshold=-10)
	#print(outfile)

	if process_background == True:
		bkg_calfiles = glob.glob(level2outputdir + '/bkg/*cal.fits') # cal files output from L2 pipeline
		for file in bkg_calfiles:
			outfile = file[:-5] + '2.fits'
			cut_cal(file, outfile, max_threshold=calmax, min_threshold=-10)
			#print(outfile)

	input_spec2 = output_path + folders_L2b[1]
	clean_dir = output_path + folders_L2b[1] # just chuck them all in the same directory. Fuck me up fam. 
	#clean_data = True # turn to True to flag additional pixels - False will not set additional DQ flags
	if clean_data:
		for file in glob.glob(input_spec2+'*cal2.fits'):
			run_clean_step(file, clean_dir, expand_flat=False)

	if target == 'COSMIC.EYE':
		target = 'COSMICEYE' # for better saving purposes

	# Make association files 
	calfiles = glob.glob(output_path + folders_L2b[1] + '*cal2.fits') # changed in this version to save fewer intermediate steps

	#bkgfiles = glob.glob(output_path + folders_L2b[2] + '*x1d.fits')
	if process_background == True:
		bkgfiles_cal = glob.glob(output_path + folders_L2b[2] + '*cal2.fits') # for background testing only

	version = 'nobg' 
	#if process_background == True:
	#    version = 'standard'


	if version == 'standard':
		if len(folders_L3) == 1:
			fold = folders_L3[0]
			asnfile = os.path.join(output_path, fold, 'L3asn.json')
		else:
			fold = folders_L3[1]
			asnfile = os.path.join(output_path, fold, 'L3asn.json')
			lev3asnname = 'Level3_' + target + '_ALLOBS'
			writel3asn(calfiles, bkgfiles, asnfile, lev3asnname) 
		
	if version == 'bgonly':
		fold = folders_L3[0]
		asnfile = os.path.join(output_path, fold, 'L3asn.json')
		lev3asnname = 'Level3_' + target + '_BGONLY'
		writel3asn(bkgfiles_cal, None, asnfile, lev3asnname)

	if version == 'nobg':
		if len(folders_L3) == 1:
			fold = folders_L3[0]
			asnfile = os.path.join(output_path, fold, 'L3asn.json')
		else:
			fold = folders_L3[3]
			asnfile = os.path.join(output_path, fold, 'L3asn.json')
			lev3asnname = 'Level3_' + target + '_NOBG_IFUALIGN_XY0p1'
			writel3asn(calfiles, None, asnfile, lev3asnname)
			# I guess this is deprecated since I put it all in the loop below to be able to do multiple resolutions
			# oh well
					
	# make multiple cubes at different resolutions:
	if runspec3 == True:
		reslist = [0.1, 0.05]
		print(f'Starting Spec3 Pipeline for {target}')
		for res in reslist:
			print(f'Starting res: {res}')
			if target == 'SUNBURST':
				calpath = os.path.join(output_path, folders_L2b[1]) # for convenience
				res_str = str(res).split('.')[-1]
				calfiles_p1 = glob.glob(calpath + '*_02101_*cal2.fits') + glob.glob(calpath + '*_10101_*cal2.fits')
				calfiles_p2 = glob.glob(calpath + '*_04101_*cal2.fits') + glob.glob(calpath + '*_12101_*cal2.fits')
				calfiles_p3 = glob.glob(calpath + '*_06101_*cal2.fits') + glob.glob(calpath + '*_14101_*cal2.fits')
				pos1asnfile = os.path.join(output_path, folders_L3[0], 'P1asn.json')
				lev3asnnameP1 = 'Level3_'+target+'_P1_NOBG_SKYALIGN_XY0p'+res_str
				writel3asn(calfiles_p1, None, pos1asnfile, lev3asnnameP1)
				pos23asnfile = os.path.join(output_path, folders_L3[0], 'P2P3asn.json')
				lev3asnnameP23 = 'Level3_'+target+'_P2P3_NOBG_SKYALIGN_XY0p'+res_str
				writel3asn(calfiles_p2+calfiles_p3, None, pos23asnfile, lev3asnnameP23)
				asnlist = [pos1asnfile,pos23asnfile]
				offsetfileP1 = os.path.join(output_path, fold, 'offsetsP1.asdf')
				make_offset_file(target, calfiles_p1, offsetfileP1)
				offsetfileP23 = os.path.join(output_path, fold, 'offsetsP2P3.asdf')
				make_offset_file(target, calfiles_p2+calfiles_p3, offsetfileP23)
			else:
				fold = folders_L3[0]
				asnfile = os.path.join(output_path, fold, 'L3asn.json')
				res_str = str(res).split('.')[-1]
				lev3asnname = 'Level3_' + target + '_NOBG_NIRCAMALIGN_XY0p'+res_str 
				writel3asn(calfiles, None, asnfile, lev3asnname)
				offsetfile = os.path.join(output_path, fold, 'offsets.asdf')
				make_offset_file(target, calfiles, offsetfile)
				spec3 = Spec3Pipeline()
				spec3.output_dir = output_path + fold # use folder defined in cell above - output to same dir as asn file
				# Outlier detection
				spec3.outlier_detection.skip = False
				spec3.outlier_detection.ifu_second_check = True
				# cube building parameters
				#spec3.cube_build.coord_system = 'ifualign'
				spec3.cube_build.coord_system = 'skyalign'
				spec3.cube_build.scalexy = res
				###############################################
			if target == 'COSMICEYE':
				spec3.cube_build.nspax_x = 59 
				spec3.cube_build.nspax_y = 65 
				spec3.cube_build.ra_center  = 323.8028721 
				spec3.cube_build.dec_center = -1.0287456  
			if target == 'SGAS1110':
				spec3.cube_build.nspax_x = 59
				spec3.cube_build.nspax_y = 79
				spec3.cube_build.ra_center  = 167.5832686
				spec3.cube_build.dec_center = 64.9974530
			if target == 'SGAS1429':
				spec3.cube_build.nspax_x = 55
				spec3.cube_build.nspax_y = 51
				spec3.cube_build.ra_center  = 217.4787892
				spec3.cube_build.dec_center = 12.0439266
			if target == 'SGAS1527':
				spec3.cube_build.nspax_x = 53
				spec3.cube_build.nspax_y = 53
				spec3.cube_build.ra_center  = 231.9388828
				spec3.cube_build.dec_center = 6.8719792
			if target == 'SGAS1050':
				spec3.cube_build.nspax_x = 109
				spec3.cube_build.nspax_y = 79
				spec3.cube_build.ra_center  = 162.6633975
				spec3.cube_build.dec_center = 0.2910689
			if target == 'SGAS2111':
				spec3.cube_build.nspax_x = 79
				spec3.cube_build.nspax_y = 75
				spec3.cube_build.ra_center  = 317.8276240
				spec3.cube_build.dec_center = -1.2412612
			if target != 'SUNBURST':
				spec3.cube_build.offset_file = offsetfile
							###############################################
			spec3.extract_1d.skip = True
			spec3.save_results = True # DON'T FORGET THIS OR YOU'LL WASTE SEVERAL HOURS FOR NAUGHT!
			if target == 'SUNBURST':
				for i,asn in enumerate(asnlist):
					##### COMMENTED ONLY FOR ONE RUN - ADD BACK IN FOR NEXT PRODUCTION VERSION
					#if i == 0:
					#	# P1
					#	spec3.cube_build.nspax_x = 67
					#	spec3.cube_build.nspax_y = 61
					#	spec3.cube_build.ra_center  = 237.5198779
					#	spec3.cube_build.dec_center = -78.1830012
					#	spec3.cube_build.offset_file = offsetfileP1
					#if i == 1:
					#	# P2+P3
					#	spec3.cube_build.nspax_x = 65
					#	spec3.cube_build.nspax_y = 75
					#	spec3.cube_build.ra_center  = 237.5008791
					#	spec3.cube_build.dec_center = -78.1864364
					#	spec3.cube_build.offset_file = offsetfileP23
					spec3.run(asn)
					# run bg subtraction:
					do_bgsub_step(target, output_path, folders_L3, lev3asnnameP1)
					do_bgsub_step(target, output_path, folders_L3, lev3asnnameP23)
			else:
				spec3.run(asnfile)
				# run bg subtraction:
				do_bgsub_step(target, output_path, folders_L3, lev3asnname)
				print(f'Done res: {res}')

		if target == 'SGAS1110': # handle wopr:
			calfiles2 = glob.glob(output_path + folders_L2b[1] + 'jw03843005*cal2.fits') # wopr only = obs 005
			for res in reslist:
				print(f'Starting wopr-only res: {res}')
				fold = folders_L3[0]
				asnfile = os.path.join(output_path, fold, 'L3asn.json')
				res_str = str(res).split('.')[-1]
				lev3asnname = 'Level3_' + target + '_NOBG_NIRCAMALIGN_WOPRONLY_XY0p'+res_str
				writel3asn(calfiles2, None, asnfile, lev3asnname)
				offsetfile2 = os.path.join(output_path, fold, 'offsets_wopr.asdf')
				make_offset_file(target, calfiles2, offsetfile2)
				spec3 = Spec3Pipeline()
				spec3.output_dir = output_path + fold # use folder defined in cell above - output to same dir as asn file
				# Outlier detection
				spec3.outlier_detection.skip = False
				spec3.outlier_detection.ifu_second_check = True
				# cube building parameters
				#spec3.cube_build.coord_system = 'ifualign'
				spec3.cube_build.coord_system = 'skyalign'
				spec3.cube_build.scalexy = res
				###############################################
				spec3.cube_build.nspax_x = 59
				spec3.cube_build.nspax_y = 79
				spec3.cube_build.ra_center  = 167.5832686
				spec3.cube_build.dec_center = 64.9974530
				spec3.cube_build.offset_file = offsetfile2
				###############################################
				spec3.save_results = True # DON'T FORGET THIS OR YOU'LL WASTE SEVERAL HOURS FOR NAUGHT!
				spec3.run(asnfile)
				# run bg subtraction:
				do_bgsub_step(target, output_path, folders_L3, lev3asnname, special=True)
				print(f'Done wopr-only res: {res}')

		print(f'Done Spec3 Pipeline for {target}')
					# and now, run spec3 cube building with non-linear prism wavelength solutions:
	if runspec3_PRISM:
		# get calfiles for prism only:
		calfiles_prism = []
		for file in calfiles:
			with fits.open(file) as hdu:
				if hdu[0].header['GRATING'] == 'PRISM':
					calfiles_prism.append(file)
					# print(calfiles_prism)
					asnfile = os.path.join(output_path, fold, 'L3PRISMasn.json')
					reslist = [0.1]#, 0.03]
		for res in reslist:
			print(f'Starting non-linear prism, res {res}')
			res_str = str(res).split('.')[-1]
			lev3asnname = 'Level3_' + target + '_NOBG_MULTICUBE_XY0p'+res_str
			writel3asn(calfiles_prism, None, asnfile, lev3asnname)
			spec3 = Spec3Pipeline()
			spec3.output_dir = output_path + fold # use folder defined in cell above - output to same dir as asn file
			# Outlier detection
			spec3.outlier_detection.skip = False
			spec3.outlier_detection.ifu_second_check = True
			# cube building parameters
			spec3.cube_build.coord_system = 'ifualign'
			spec3.cube_build.scalexy = res
			spec3.cube_build.output_type = 'multi' ###################################### ONLY USE THIS FOR PRISM STUFF ###############################################
			spec3.save_results = True # DON'T FORGET THIS OR YOU'LL WASTE SEVERAL HOURS FOR NAUGHT!
			spec3.run(asnfile)
			print(f'Done res {res}')





###############################################
## HELPER FUNCTIONS FOR RUNNING THE PIPELINE ##
###############################################


# BELOW COPIED FROM DAVID LAW'S MIRI MRS NOTEBOOK: 
# https://github.com/STScI-MIRI/MRS-ExampleNB/blob/main/Flight_Notebook1/MRS_FlightNB1.ipynb
# 
# Define a useful function to write out a Lvl3 association file from an input list
# Note that any background exposures have to be of type x1d.
def writel3asn(scifiles, bgfiles, asnfile, prodname):
	# Define the basic association of science files
	asn = afl.asn_from_list(scifiles, rule=DMS_Level3_Base, product_name=prodname)

	# Add background files to the association
	if bgfiles:
		nbg=len(bgfiles)
		for ii in range(0,nbg):
			asn['products'][0]['members'].append({'expname': bgfiles[ii], 'exptype': 'background'})
		
	# Write the association to a json file
	_, serialized = asn.dump()
	with open(asnfile, 'w') as outfile:
		outfile.write(serialized)

def make_offset_file(target, calfiles, offsetfile):
	if target == 'COSMICEYE':
		prism_dra  = 0.04427999992913101
		prism_ddec = 0.13320000000023313
		g235m_dra  = -0.016199999959098932
		g235m_ddec = -0.16164000000005174
	if target == 'SGAS1429':
		prism_dra  = 0.0748
		prism_ddec = 0.1256
		g235m_dra  = 0.0748
		g235m_ddec = 0.1256 
	if target == 'SGAS1110':
		prism_dra  = 0.4 - 0.2
		prism_ddec = 0.01368
		g235m_dra  = 0.4 - 0.2
		g235m_ddec = 0.01368
		g140h_dra  = 0.11 - 0.05
		g140h_ddec = -0.14
	if target == 'SGAS1527':
		prism_dra  = 0.154
		prism_ddec = 0.065
		g235m_dra  = 0.104 
		g235m_ddec = 0.103
	if target == 'SGAS1050':
		prism_dra  = 0.09
		prism_ddec = 0.1184
		g235m_dra  = 0.084
		g235m_ddec = 0.1116
	if target == 'SGAS2111':
		prism_dra  = -0.039
		prism_ddec = -0.134
		g235m_dra  = -0.039
		g235m_ddec = -0.134
	if target == 'SUNBURST':
		g235m_dra  = -1.25#6.37
		g235m_ddec = 0.48
		g140h_dra  = -1.25#6.37
		g140h_ddec = 0.48
		ra_offsets  = []
		dec_offsets = []
		filelist    = []
		#print(calfiles)
	for file in calfiles:
		head = fits.open(file)[0].header
		grating = head['GRATING']
		if grating == 'PRISM':
			ra_offsets.append(prism_dra)
			dec_offsets.append(prism_ddec)
		if grating == 'G235M':
			ra_offsets.append(g235m_dra)
			dec_offsets.append(g235m_ddec)
		if grating == 'G140H':
			ra_offsets.append(g140h_dra)
			dec_offsets.append(g140h_ddec)
		if grating == 'G235H':
			ra_offsets.append(g235m_dra) # use g235m also for g235h - they're close enough by eye
			dec_offsets.append(g235m_ddec)
			filelist.append(file.split('/')[-1])
			tree = {
				"units": str(u.arcsec),
				"filename":filelist,
				"raoffset": ra_offsets,
				"decoffset": dec_offsets
			}
	with asdf.AsdfFile(tree) as af:
		af.write_to(offsetfile)


def cut_cal(infile, outfile, max_threshold=20, min_threshold=-1):
	hdu=fits.open(infile)
	sci=hdu['SCI'].data
	dq=hdu['DQ'].data

	dnubit=dqflags.interpret_bit_flags('DO_NOT_USE', mnemonic_map=datamodels.dqflags.pixel)
	indx=np.where((dq & dnubit) != 0)
	sci[indx]=np.nan

	indx=np.where((sci > max_threshold) | (sci < min_threshold))
	sci[indx]=np.nan
	dq[indx] = np.bitwise_or(dq[indx], dnubit)

	hdu['SCI'].data=sci
	hdu.writeto(outfile, overwrite=True)



def run_clean_step(file, clean_dir, expand_flat=False):
	DO_NOT_USE = datamodels.dqflags.pixel["DO_NOT_USE"]
	NON_SCIENCE= datamodels.dqflags.pixel["NON_SCIENCE"]
	NO_SAT_CHECK = datamodels.dqflags.pixel["NO_SAT_CHECK"] # 2^21 = 2097152
	UNRELIABLE_FLAT = datamodels.dqflags.pixel["UNRELIABLE_FLAT"] # 2^25 = 33554432
	OTHER_BAD_PIXEL = datamodels.dqflags.pixel['OTHER_BAD_PIXEL'] # 2^30 = 1073741824
	TELEGRAPH_PIXEL = datamodels.dqflags.pixel['TELEGRAPH'] # 2^15 = 32768
	MSA_FAILED_OPEN = datamodels.dqflags.pixel['MSA_FAILED_OPEN'] # 2^29 = 536870912
	FLUX_ESTIMATED = datamodels.dqflags.pixel['FLUX_ESTIMATED'] # 2^28 = 268435456
	UNRELIABLE_DARK = datamodels.dqflags.pixel['UNRELIABLE_DARK'] # 2^23 = 8388608

	filename = os.path.basename(file)
	#filename = filename.replace('cal2','cal3') # don't update the names for this run, to save some disk space

	new_file = clean_dir + filename  # write to a different directory to keep things straight
	test_file = filename.replace('.fits','test.fits') # this is just for testing and looking at where the new bad pixels have 
	# been flagged
	new_test_file = clean_dir + test_file
	print('test file', new_test_file)

	print('new file', new_file)
	input_cal = datamodels.IFUImageModel(file)
	input_test = input_cal.copy()

	non_science = np.bitwise_and(input_cal.dq, NON_SCIENCE).astype(bool)
	bad = np.bitwise_and(input_cal.dq, DO_NOT_USE).astype(bool)
	num_science = np.where(~non_science & ~bad)
	num_science = len(num_science[0])
	print('number of good science pixels at the start', num_science)
	loc_non_science = np.where(np.bitwise_and(input_cal.dq, NON_SCIENCE).astype(bool))

	data1 = input_cal.data
	dq1 = input_cal.dq
	dq1_donotuse = np.where(np.bitwise_and(dq1, DO_NOT_USE).astype(bool))
	dq_test =  dq1.copy() * 0 # initialize to 0 
	dq_test[dq1_donotuse] = 1 # set up bad pixels
	dq_test[loc_non_science] = 1 

	uflat = np.bitwise_and(dq1, UNRELIABLE_FLAT).astype(bool)
	bad =  np.bitwise_and(dq1, DO_NOT_USE).astype(bool)
	non_science =  np.bitwise_and(dq1, NON_SCIENCE).astype(bool)
	test1 = np.where( uflat & ~bad & ~non_science)
	num = len(test1[0])

	# flag new unreliable flat 
	dq_test[test1] = 3  
	##

	if num > 0:
		print('Number of pixels with UNREIABLE_FLAT but no DO_NOT_USE', num)
		print('% valid pixels',(num/num_science)* 100)
		input_cal.dq[test1] = np.bitwise_or(input_cal.dq[test1], datamodels.dqflags.pixel['DO_NOT_USE'])
		input_cal.data[test1] = np.nan


	# NO_SAT_CHECK
	data2 = input_cal.data
	dq2 = input_cal.dq

	no_sat = np.bitwise_and(dq2, NO_SAT_CHECK).astype(bool)
	bad =  np.bitwise_and(dq2, DO_NOT_USE).astype(bool)
	non_science =  np.bitwise_and(dq2, NON_SCIENCE).astype(bool) # check if it is a NON-Science pixel
	test2 = np.where( no_sat & ~bad & ~non_science)
	num = len(test2[0])

	# set test DQ
	dq_test[test2] = 5
	##

	if num > 0:
		print('Number of pixels with NO_SAT_CHECK but no DO_NOT_USE ', num)
		print('% valid pixels',(num/num_science)* 100)
		input_cal.dq[test2] = np.bitwise_or(input_cal.dq[test2], datamodels.dqflags.pixel['DO_NOT_USE']) 
		input_cal.data[test2] = np.nan


	# OTHER_BAD_PIXEL 
	data3 = input_cal.data
	dq3 = input_cal.dq

	other_bad = np.bitwise_and(dq3, OTHER_BAD_PIXEL).astype(bool)
	bad =  np.bitwise_and(dq3, DO_NOT_USE).astype(bool)
	non_science =  np.bitwise_and(dq3, NON_SCIENCE).astype(bool)
	test3 = np.where( other_bad & ~bad & ~non_science)
	num = len(test3[0])

	# set test DQ
	dq_test[test3] = 7
	##


	if num > 0:
		print('Number of pixels with OTHER_BAD_PIXEL but not DO_NOT_USE', num)
		print('% valid pixels',(num/num_science)* 100)
		input_cal.dq[test3] = np.bitwise_or(input_cal.dq[test3], datamodels.dqflags.pixel['DO_NOT_USE'])
		input_cal.data[test3] = np.nan

	# TELEGRAPH Pixel
	data4 = input_cal.data
	dq4 = input_cal.dq

	# set test DQ
	tp = np.bitwise_and(dq4, TELEGRAPH_PIXEL).astype(bool)
	bad =  np.bitwise_and(dq4, DO_NOT_USE).astype(bool)
	non_science =  np.bitwise_and(dq4, NON_SCIENCE).astype(bool)
	test4 = np.where( tp & ~bad & ~non_science)
	num = len(test4[0])

	#
	dq_test[test4] = 9
	##


	if num > 0:  
		print('Number of pixels with Random Telegraph  but no DO_NOT_USE', num)
		print('% valid pixels',(num/num_science)* 100)
		input_cal.dq[test4] = np.bitwise_or(input_cal.dq[test4], datamodels.dqflags.pixel['DO_NOT_USE'])
		input_cal.data[test4] = np.nan


	# MSA_FAILED_OPEN Pixel
	data5 = input_cal.data
	dq5 = input_cal.dq

	# set test DQ
	msa = np.bitwise_and(dq5, MSA_FAILED_OPEN).astype(bool)
	bad =  np.bitwise_and(dq5, DO_NOT_USE).astype(bool)
	non_science =  np.bitwise_and(dq5, NON_SCIENCE).astype(bool)
	test5 = np.where( msa & ~bad & ~non_science)
	num = len(test5[0])

	#
	dq_test[test5] = 11
	##


	if num > 0:  
		print('Number of pixels with Failed Open Shutter but no DO_NOT_USE', num)
		print('% valid pixels',(num/num_science)* 100)
		input_cal.dq[test5] = np.bitwise_or(input_cal.dq[test5], datamodels.dqflags.pixel['DO_NOT_USE'])
		input_cal.data[test5] = np.nan


	# FLUX_ESTIMATED Pixel
	data6 = input_cal.data
	dq6 = input_cal.dq

	# set test DQ
	fe = np.bitwise_and(dq6, FLUX_ESTIMATED).astype(bool)
	bad =  np.bitwise_and(dq6, DO_NOT_USE).astype(bool)
	non_science =  np.bitwise_and(dq6, NON_SCIENCE).astype(bool)
	test6 = np.where( fe & ~bad & ~non_science)
	num = len(test6[0])

	#
	dq_test[test6] = 13
	##


	if num > 0:  
		print('Number of pixels with Flux Estimated  but no DO_NOT_USE', num)
		print('% valid pixels',(num/num_science)* 100)
		input_cal.dq[test6] = np.bitwise_or(input_cal.dq[test6], datamodels.dqflags.pixel['DO_NOT_USE'])
		input_cal.data[test6] = np.nan


	# UNRELIABLE_DARK Pixel
	data7 = input_cal.data
	dq7 = input_cal.dq

	# set test DQ
	ud = np.bitwise_and(dq7, UNRELIABLE_DARK).astype(bool)
	bad =  np.bitwise_and(dq7, DO_NOT_USE).astype(bool)
	non_science =  np.bitwise_and(dq7, NON_SCIENCE).astype(bool)
	test7 = np.where( ud & ~bad & ~non_science)
	num = len(test7[0])

	#
	dq_test[test7] = 14
	##


	if num > 0:  
		print('Number of pixels with Unreliable Dark but no DO_NOT_USE', num)
		print('% valid pixels',(num/num_science)* 100)
		input_cal.dq[test7] = np.bitwise_or(input_cal.dq[test7], datamodels.dqflags.pixel['DO_NOT_USE'])
		input_cal.data[test7] = np.nan


	# This one may not be need. You should check
	# now expand the unrelaible_flat flag to be 1 pixel more along the edge of the slices
	#expand_flat = True
	if expand_flat:
		dq_expand = input_cal.dq
		uflat = np.where(np.bitwise_and(dq_expand, UNRELIABLE_FLAT).astype(bool))
		nflat = len(uflat[0])
		# loop over the pixels with UNRELIABLE FLAT set. Check y+1, y-1 (not sure which edge of the slice we
		# on) and then these pixels to UNRELIABLE_FLAT
		print('Looping over ', nflat)

		for i in range (nflat):
			ix = uflat[1][i]
			iy = uflat[0][i]
			input_cal.dq[iy+1, ix] = np.bitwise_or(input_cal.dq[iy+1,ix], datamodels.dqflags.pixel['DO_NOT_USE'])
			input_cal.dq[iy+1, ix] = np.bitwise_or(input_cal.dq[iy+1,ix], datamodels.dqflags.pixel['UNRELIABLE_FLAT'])
					
			input_cal.dq[iy-1, ix] = np.bitwise_or(input_cal.dq[iy-1,ix], datamodels.dqflags.pixel['DO_NOT_USE'])
			input_cal.dq[iy-1, ix] = np.bitwise_or(input_cal.dq[iy-1,ix], datamodels.dqflags.pixel['UNRELIABLE_FLAT'])
					
			input_cal.dq[iy-2, ix] = np.bitwise_or(input_cal.dq[iy-2,ix], datamodels.dqflags.pixel['DO_NOT_USE'])
			input_cal.dq[iy-2, ix] = np.bitwise_or(input_cal.dq[iy-2,ix], datamodels.dqflags.pixel['UNRELIABLE_FLAT'])
					
			input_cal.dq[iy+2, ix] = np.bitwise_or(input_cal.dq[iy+2,ix], datamodels.dqflags.pixel['DO_NOT_USE'])
			input_cal.dq[iy+2, ix] = np.bitwise_or(input_cal.dq[iy+2,ix], datamodels.dqflags.pixel['UNRELIABLE_FLAT'])
					

	input_cal.save(new_file)
	input_test.dq = dq_test
	#input_test.save(new_test_file) # you can view the DQ plane in DS9 and easily see additional pixels that are flagged as bad
	print(test_file)




def exp_bg_sub(cubefile, bgfile, outfile=None):
	# first, load the expected background data:
	exp = np.loadtxt(bgfile)
	expwl = np.array([line[0] for line in exp])
	exptot = np.array([line[1] for line in exp])
	# and make an interp1d object sdo we can get the background at any given point:
	bginterp = interp1d(expwl, exptot)
	# now load the data cube:
	hdu = fits.open(cubefile)
	data = hdu[1].data
	head1 = hdu[1].header
	cubewl = np.arange(head1['CRVAL3'],
				 head1['CRVAL3']+(head1['CDELT3']*len(data)),
				 head1['CDELT3'])
	if len(cubewl) > len(data):
		cubewl = cubewl[:-1]
		# next, evaluate background on data wavelength points:
		bgcalculated = bginterp(cubewl)
		bgcube = np.tile(bgcalculated, (data.shape[2],data.shape[1],1)).T # cubeify!
		# and subtract!
		#print(len(data), len(cubewl))
		bgsub_cube = data - bgcube
		# save or return result:
	if outfile:
		hdu[1].data = bgsub_cube
		hdu[1].header['comment'] = 'Expected Background Subtracted'
		hdu.writeto(outfile, overwrite=True)
		hdu.close()
	else:
		hdu.close()
		return bgsub_cube


def do_bgsub_step(target, output_path, folders_L3, lev3asnname, special=False):
	## Background Subtraction:
	# starting with LEGGOS targets:
	if target == 'SGAS1110':
		if special == True:
			expbg_file = '/Users/bwelch/Documents/data/leggos/sgas1110/nirspec/expected_bg/background.txt'
			g235hcube = os.path.join(output_path, folders_L3[0], lev3asnname+'_g235h-f170lp_s3d.fits')
			g235hout = g235hcube.replace('NOBG','BGSUB')
			exp_bg_sub(g235hcube, expbg_file, g235hout)
		else:
			expbg_file = '/Users/bwelch/Documents/data/leggos/sgas1110/nirspec/expected_bg/background.txt'
			g235hcube = os.path.join(output_path, folders_L3[0], lev3asnname+'_g235h-f170lp_s3d.fits')
			g235mcube = os.path.join(output_path, folders_L3[0], lev3asnname+'_g235m-f170lp_s3d.fits')
			g140cube = os.path.join(output_path, folders_L3[0], lev3asnname+'_g140h-f100lp_s3d.fits')
			g235hout = g235hcube.replace('NOBG','BGSUB')
			g235mout = g235mcube.replace('NOBG','BGSUB')
			g140out = g140cube.replace('NOBG','BGSUB')
			exp_bg_sub(g140cube, expbg_file, g140out)
			exp_bg_sub(g235hcube, expbg_file, g235hout)
			exp_bg_sub(g235mcube, expbg_file, g235mout)
	if target == 'SGAS2111':
		expbg_file = '/Users/bwelch/Documents/data/leggos/sgas2111/nirspec/expected_bg/background.txt'
		g235cube = os.path.join(output_path, folders_L3[0], lev3asnname+'_g235m-f170lp_s3d.fits')
		g235out = g235cube.replace('NOBG', 'BGSUB')
		exp_bg_sub(g235cube, expbg_file, g235out)
	if target == 'SGAS1050':
		expbg_file = '/Users/bwelch/Documents/data/leggos/sgas1050/nirspec/expected_bg/background.txt'
		g235cube = os.path.join(output_path, folders_L3[0], lev3asnname+'_g235m-f170lp_s3d.fits')
		g235out = g235cube.replace('NOBG', 'BGSUB')
		exp_bg_sub(g235cube, expbg_file, g235out)
	if target == 'SGAS1429':
		expbg_file = '/Users/bwelch/Documents/data/leggos/sgas1429/nirspec/expected_bg/background.txt'
		g235cube = os.path.join(output_path, folders_L3[0], lev3asnname+'_g235m-f170lp_s3d.fits')
		g235out = g235cube.replace('NOBG', 'BGSUB')
		exp_bg_sub(g235cube, expbg_file, g235out)
	if target == 'SGAS1527':
		expbg_file = '/Users/bwelch/Documents/data/leggos/sgas1527/nirspec/expected_bg/background.txt'
		g235cube = os.path.join(output_path, folders_L3[0], lev3asnname+'_g235m-f170lp_s3d.fits')
		g235out = g235cube.replace('NOBG', 'BGSUB')
		exp_bg_sub(g235cube, expbg_file, g235out)
	if (target == 'COSMIC.EYE') or (target == 'COSMICEYE'):
		expbg_file = '/Users/bwelch/Documents/data/leggos/cosmic_eye/nirspec/expected_bg/background.txt'
		g235cube = os.path.join(output_path, folders_L3[0], lev3asnname+'_g235m-f170lp_s3d.fits')
		g235out = g235cube.replace('NOBG', 'BGSUB')
		exp_bg_sub(g235cube, expbg_file, g235out)
# and now for the TEMPLATES targets:
	if target == 'SGAS1723':
		expbg_file = '/Users/bwelch/Documents/data/templates/sdss1723/nirspec/expected_bg/background.txt'
		g395cube = os.path.join(output_path, folders_L3[0], lev3asnname+'_g395h-f290lp_s3d.fits')
		g140cube = os.path.join(output_path, folders_L3[0], lev3asnname+'_g140h-f100lp_s3d.fits')
		g395out = g395cube.replace('NOBG','BGSUB')
		g140out = g140cube.replace('NOBG','BGSUB')
		exp_bg_sub(g140cube, expbg_file, g140out)
		exp_bg_sub(g395cube, expbg_file, g395out)
	if target == 'SGAS1226':
		expbg_file = '/Users/bwelch/Documents/data/templates/sdss1226/nirspec/expected_bg/background.txt'
		g235cube = os.path.join(output_path, folders_L3[0], lev3asnname+'_g235h-f170lp_s3d.fits')
		g235out = g235cube.replace('NOBG', 'BGSUB')
		exp_bg_sub(g235cube, expbg_file, g235out)
	if target == 'SPT0418-47':
		expbg_file = '/Users/bwelch/Documents/data/templates/spt0418/nirspec/expected_bg/background.txt'
		g395cube = os.path.join(output_path, folders_L3[0], lev3asnname+'_g395m-f290lp_s3d.fits')
		g395out = g395cube.replace('NOBG', 'BGSUB')
		exp_bg_sub(g395cube, expbg_file, g395out)
	if target == 'SPT2147-50':
		expbg_file = '/Users/bwelch/Documents/data/templates/spt2147/nirspec/expected_bg/background.txt'
		g395cube = os.path.join(output_path, folders_L3[0], lev3asnname+'_g395m-f290lp_s3d.fits')
		g395out = g395cube.replace('NOBG', 'BGSUB')
		exp_bg_sub(g395cube, expbg_file, g395out)
# and other stuff: 
	if target == 'SPARKLER':
		expbg_file = '/Users/bwelch/Documents/data/sparkler/nirspec/expected_bg/background.txt'
		g235cube = os.path.join(output_path, folders_L3[0], lev3asnname+'_g140m-f070lp_s3d.fits')
		g235out = g235cube.replace('NOBG', 'BGSUB')
		exp_bg_sub(g235cube, expbg_file, g235out)
	if target == 'SUNBURST':
		expbg_file = '/Users/bwelch/Documents/data/sgas/sunburst/nirspec/expected_bg/background.txt'
		g140cube = os.path.join(output_path,folders_L3[0],lev3asnname+'_g140h-f100lp_s3d.fits')
		g140out = g140cube.replace('NOBG','BGSUB')
		g235cube = os.path.join(output_path,folders_L3[0],lev3asnname+'_g235h-f170lp_s3d.fits')
		g235out = g235cube.replace('NOBG','BGSUB')
		exp_bg_sub(g140cube, expbg_file, g140out)
		exp_bg_sub(g235cube, expbg_file, g235out)



