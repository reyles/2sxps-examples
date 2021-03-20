#Given a 2SXPS source name, performs photometry on the UVOT data at that sources position.
from __future__ import division
import numpy as np
import time
import warnings
from numpy import *
import astropy.units as u
from astropy import coordinates
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
import os
import pandas as pd
from builtins import any

#This takes a filter, does the photometry and returns a dataframe of the sources within the 2SXPS error region
def photometry(src,ob,filter):
    print('Current filter: %s' % filter)
    #Makes lists to return data
    ob_list, mjd_list, mjderr_list, filter_list, ra_list, dec_list, mag_list, magerr_list, ul_list = [],[],[],[],[],[],[],[],[]
    #Calculates the time
    image = fits.open('data/%s/%s/uvot/image/sw%su%s_sk.img.gz' % (src,ob,ob,filter))
    header = image[1].header
    start_time = header['DATE-OBS']
    end_time = header['DATE-END']
    start_mjd = Time(start_time, format='isot', scale='utc').mjd
    end_mjd = Time(end_time, format='isot', scale='utc').mjd
    obs_mjd = (start_mjd + end_mjd) / 2
    obs_mjderr = obs_mjd - start_mjd
    #Does the initial photometry pass
    os.system('uvotsource image=data/%s/%s/uvot/image/sw%su%s_sk.img.gz+1 srcreg=src.reg bkgreg=bkg.reg sigma=5 outfile=phot.fits clobber=yes >/dev/null 2>&1'
              % (src,ob,ob,filter))
    #Determines whether there was a 5-sigma detection
    phot = fits.open('phot.fits')
    phot_data = phot[1].data
    if phot_data['AB_MAG'] < phot_data['AB_MAG_LIM']:
        ob_list.append(ob)
        mjd_list.append(obs_mjd)
        mjderr_list.append(obs_mjderr)
        filter_list.append(filter)
        ra_list.append(phot_data['RA'][0])
        dec_list.append(phot_data['DEC'][0])
        mag_list.append(phot_data['AB_MAG'][0])
        magerr_list.append(phot_data['AB_MAG_ERR'][0])
        ul_list.append(0)
    #If there wasn't a detection, run it again to find a 3-sigma limit
    else:
        os.system('uvotsource image=data/%s/%s/uvot/image/sw%su%s_sk.img.gz+1 srcreg=src.reg bkgreg=bkg.reg sigma=3 outfile=phot.fits clobber=yes >/dev/null 2>&1'
                  % (src,ob,ob,filter))
        phot = fits.open('phot.fits')
        phot_data = phot[1].data
        if phot_data['AB_MAG'] < phot_data['AB_MAG_LIM']:
            #3-sigma background limit higher than detection magnitude: take 3-sigma detection as upper limit
            ob_list.append(ob)
            mjd_list.append(obs_mjd)
            mjderr_list.append(obs_mjderr)
            filter_list.append(filter)
            ra_list.append(phot_data['RA'][0])
            dec_list.append(phot_data['DEC'][0])
            mag_list.append(phot_data['AB_MAG'][0])
            magerr_list.append(phot_data['AB_MAG_ERR'][0])
            ul_list.append(2)
        else:
            #3-sigma background limit lower than detection magnitude: take the background limit as an upper limit
            ob_list.append(ob)
            mjd_list.append(obs_mjd)
            mjderr_list.append(obs_mjderr)
            filter_list.append(filter)
            ra_list.append(phot_data['RA'][0])
            dec_list.append(phot_data['DEC'][0])
            mag_list.append(phot_data['AB_MAG_LIM'][0])
            magerr_list.append(0)
            ul_list.append(1)
            
    #Zips the lists to make a dataframe
    photometry_df = pd.DataFrame(list(zip(ob_list,mjd_list, mjderr_list, filter_list, ra_list, dec_list, mag_list, magerr_list, ul_list)),
                                 columns = ['Ob','MJD','MJD_err','Filter','RA','Dec','Mag','Mag_err','Upper_limit'])
    #print(photometry_df)
    return photometry_df

#This takes a filter, detects all the sources in it and returns the dataframe
def detect(src,ob,filter,src_pos,src_poserr):
    print('Current filter: %s' % filter)
    #Does the photometry
    os.system('uvotdetect infile=data/%s/%s/uvot/image/sw%su%s_sk.img.gz+1 outfile=sources.fits expfile=data/%s/%s/uvot/image/sw%su%s_ex.img.gz threshold=5 clobber=YES calibrate=YES >/dev/null 2>&1'
              % (src,ob,ob,filter,src,ob,ob,filter))
    #Calculates the time
    image = fits.open('data/%s/%s/uvot/image/sw%su%s_sk.img.gz' % (src,ob,ob,filter))
    header = image[1].header
    start_time = header['DATE-OBS']
    end_time = header['DATE-END']
    start_mjd = Time(start_time, format='isot', scale='utc').mjd
    end_mjd = Time(end_time, format='isot', scale='utc').mjd
    obs_mjd = (start_mjd + end_mjd) / 2
    obs_mjderr = obs_mjd - start_mjd
    #Figures out which sources are close enough to the source
    uvotdetect = fits.open('sources.fits')
    uvot_data = uvotdetect[1].data
    #Adds them into lists
    mjd_list, mjderr_list, filter_list, ra_list, raerr_list, dec_list, decerr_list, majaxis_list, mag_list, magerr_list = [],[],[],[],[],[],[],[],[],[]
    for i in range(uvot_data.shape[0]):
            mjd_list.append(obs_mjd)
            mjderr_list.append(obs_mjderr)
            filter_list.append(filter)
            ra_list.append(uvot_data['RA'][i])
            raerr_list.append(uvot_data['RA_ERR'][i])
            dec_list.append(uvot_data['DEC'][i])
            decerr_list.append(uvot_data['DEC_ERR'][i])
            majaxis_list.append(uvot_data['PROF_MAJOR'][i])
            mag_list.append(uvot_data['MAG'][i])
            magerr_list.append(uvot_data['MAG_ERR'][i])
    #Zips the lists to make a dataframe
    detect_df = pd.DataFrame(list(zip( mjd_list, mjderr_list, filter_list, ra_list, raerr_list, dec_list, decerr_list, majaxis_list, mag_list, magerr_list)),
                                 columns = ['MJD','MJD_err','Filter','RA','RA_err','Dec','Dec_err','Major_axis','Mag','Mag_err'])
    #print(detect_df)
    return detect_df

#Finds a random point far away from detected sources to centre the background region around
def find_bkg(detect_df,dist):
    #Converts detections dataframe into a SkyCoord catalogue for cross matching
    detect_cat = SkyCoord(ra=list(detect_df['RA'])*u.degree, dec=list(detect_df['Dec'])*u.degree)
    min_sep = 0
    i = 0
    #Generates random ra/dec positions till suitably far away from other sources. Normally done first time so MCMC etc are overkill
    while min_sep < dist * u.arcsec:
        rand_ra = np.random.uniform(np.min(list(detect_df['RA'])),np.max(list(detect_df['RA'])))
        rand_dec = np.random.uniform(np.min(list(detect_df['Dec'])),np.max(list(detect_df['Dec'])))
        rand_pos = SkyCoord(ra=rand_ra*u.degree, dec=rand_dec*u.degree)
        d2d = rand_pos.separation(detect_cat)
        min_sep = np.min(d2d)
        i = i + 1

    return rand_ra,rand_dec

#The below sets it up and does the loop to do the photometry
#Enter source name
src = 'J105534.3-012614'

#Look up source's coordinates and error region
fullcat = pd.read_csv('../cats_2sxps/2SXPS_Sources.csv', low_memory=False)
full_2sxps = fullcat.loc[(fullcat['IAUName'] == str('2SXPS ' + src))]
print('Data read!')

#Convert to SkyCoord position
src_ra = float(full_2sxps['RA'])
src_dec = float(full_2sxps['Decl'])
src_poserr = float(full_2sxps['Err90'])
src_pos = SkyCoord(ra=src_ra*u.degree, dec=src_dec*u.degree)


#Figure out the files to read
obs = os.listdir('data/%s' % src)
n_obs = len(obs)
n_complete = 0

print('Number of obs: %s' % n_obs)

#Creates a blank dataframe to store the detected source
detect_master = pd.DataFrame(columns = ['MJD','MJD_err','Filter','RA','RA_err','Dec','Dec_err','Major_axis','Mag','Mag_err'])

#Go through the obs and detect sources
for ob in obs:
    try:
        files = os.listdir('data/%s/%s/uvot/image' % (src,ob))
        #Grab the detections for each filter
        if any('wh' in f for f in files) == True:
            detect_master = pd.concat([detect_master,detect(src,ob,'wh',src_pos,src_poserr)])
        if any('uu' in f for f in files) == True:
            detect_master = pd.concat([detect_master,detect(src,ob,'uu',src_pos,src_poserr)])
        if any('bb' in f for f in files) == True:
            detect_master = pd.concat([detect_master,detect(src,ob,'bb',src_pos,src_poserr)])
        if any('vv' in f for f in files) == True:
            detect_master = pd.concat([detect_master,detect(src,ob,'vv',src_pos,src_poserr)])
        if any('w1' in f for f in files) == True:
            detect_master = pd.concat([detect_master,detect(src,ob,'w1',src_pos,src_poserr)])
        if any('w2' in f for f in files) == True:
            detect_master = pd.concat([detect_master,detect(src,ob,'w2',src_pos,src_poserr)])
        if any('m2' in f for f in files) == True:
            detect_master = pd.concat([detect_master,detect(src,ob,'m2',src_pos,src_poserr)])
    except:
        print('No image files')
    n_complete = n_complete + 1
    detect_master.to_csv('photometry/%s_uvot_detections.csv' % src, index=False)
    print('Obs completed: %s/%s' % (n_complete,n_obs))

#Figure out the mean parameters for sources close to the expected source position
detect_cat = SkyCoord(ra=list(detect_master['RA'])*u.degree, dec=list(detect_master['Dec'])*u.degree)
d2d = src_pos.separation(detect_cat)
catalogmsk = d2d < src_poserr * u.arcsec

#Calculates mean position and size
ra_calc, dec_calc, axis_calc = [], [], []
for i in range(len(catalogmsk)):
    if catalogmsk[i] == True:
         ra_calc.append(detect_master['RA'].iloc[i])
         dec_calc.append(detect_master['Dec'].iloc[i])
         axis_calc.append(detect_master['Major_axis'].iloc[i])

mean_ra = np.mean(ra_calc)
mean_dec = np.mean(dec_calc)
mean_axis = np.mean(axis_calc)

if np.isnan(mean_ra) == True:
    mean_ra = src_ra
if np.isnan(mean_dec) == True:
    mean_dec = src_dec
if mean_axis < 5 or np.isnan(mean_axis) == True:
    mean_axis = 5
    
print(mean_ra,mean_dec,mean_axis)

#Turns the source position into a region for the photometry
f = open('src.reg', 'w')
f.write('fk5;circle(%s,%s,2.5")' % (mean_ra,mean_dec))
f.close()

#Find a background region
dist = mean_axis * 2.4
bkg_ra, bkg_dec = find_bkg(detect_master,dist)

f = open('bkg.reg', 'w')
f.write('fk5;circle(%s,%s,%s")' % (mean_ra,mean_dec,dist))
f.close()
       
#Cleanup
os.system('rm sources.fits')

photometry_master = pd.DataFrame(columns = ['Ob','MJD','MJD_err','Filter','RA','Dec','Mag','Mag_err','Upper_limit'])

n_complete = 0

#Go through the obs again and perform final photometry
for ob in obs:
    try:
        files = os.listdir('data/%s/%s/uvot/image' % (src,ob))
        #Grab the photometry for each filter
        if any('wh' in f for f in files) == True:
            photometry_master = pd.concat([photometry_master,photometry(src,ob,'wh')])
        if any('uu' in f for f in files) == True:
            photometry_master = pd.concat([photometry_master,photometry(src,ob,'uu')])
        if any('bb' in f for f in files) == True:
            photometry_master = pd.concat([photometry_master,photometry(src,ob,'bb')])
        if any('vv' in f for f in files) == True:
            photometry_master = pd.concat([photometry_master,photometry(src,ob,'vv')])
        if any('w1' in f for f in files) == True:
            photometry_master = pd.concat([photometry_master,photometry(src,ob,'w1')])
        if any('w2' in f for f in files) == True:
            photometry_master = pd.concat([photometry_master,photometry(src,ob,'w2')])
        if any('m2' in f for f in files) == True:
            photometry_master = pd.concat([photometry_master,photometry(src,ob,'m2')])
    except:
        print('No image files')
    n_complete = n_complete + 1
    #print(photometry_master)
    photometry_master = photometry_master.sort_values(['MJD'])
    photometry_master.to_csv('photometry/%s_uvot_photometry.csv' % src, index=False)
    print('Obs completed: %s/%s' % (n_complete,n_obs))
    
os.system('rm *.reg')
os.system('rm phot.fits')