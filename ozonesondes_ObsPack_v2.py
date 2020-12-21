#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 12:07:56 2020

@author: rgryan

CODE TO CONVERT WOUDC OZONESONDE FILES TO OBSPACK .NC FORMAT,
FOR USE IN THE GEOS-CHEM OBSPACK PACKAGE

WOUDC OZONESONDE FILES AVAILABLE AT https://woudc.org/data/
"""

import pandas as pd
import numpy as np
from netCDF4 import Dataset
import glob
import os

## Path to ozonesondes
path = '/Volumes/GoogleDrive/My Drive/Documents/Postdoc/rockets/'
osp = path+'ozonesondes/woudc_2019_ozonesondes/'
osp = path+'ozonesondes/woudc_2019_ozonesondes/test/'

## Using glob, this reads all .csv sonde files in the folder
files = glob.glob(osp+'*.csv')

## function to calculate datetime from sonde start datetime and 'duration'(seconds) field
def calc_DT(seconds):
        return pd.to_datetime(date+' '+startTime) + pd.Timedelta(str(seconds)+' seconds')
def calc_DT2(dtarray):
    try:
        y,mo,d,h,mi,s = dtarray[0],dtarray[1],dtarray[2],dtarray[3],dtarray[4],dtarray[5]
        return pd.to_datetime(str(d)+'/'+str(mo)+'/'+str(y)+' '+str(h)+':'+str(mi)+':'+str(s), dayfirst=True)
    except ValueError:
        pd.to_datetime(np.NaN)
#%%
## Perform the readin and conversion on each file
for osf in files[:]:
    ## Use 'Readlines' to find the header of the file and location lat/lon
    f = open(osf, "r")
    c,blanks,footer,l =0, 0, 0, [] # Counters for finding information in file
    for line in f:
        l.append(line)
        c=c+1
        if line == '\n': 
            c, blanks = c-1, blanks+1
            #print('    newLine')
        elif line[:9] == '#LOCATION': ## Find index of location information
            ilocation = c+1+blanks 
            #print('    Location')
        elif line[:6] == '#PLATF':    ## Find index of station information
            istation = c+1+blanks
        elif line[:10] == '#TIMESTAMP': ## Find timestamp information
            idate = c+1+blanks
        elif line[:10]  == '#PROFILE\n':## Find start of ozone sonde profile
            header=c
        elif line == '#PROFILE_UNCERTAINTY\n': ## Find index of appended uncertainty info, if present
            footer = c+blanks-1
        elif line == '#PRELAUNCH\n':
            print('yes')
            footer = c+blanks-1
        else:
            continue
    
    ## Get metadata from the early lines of the sonde file
    lat, lon = l[ilocation].split(',')[0], l[ilocation].split(',')[1]
    country, station = l[istation].split(',')[3], l[istation].split(',')[2]
    date, startTime = l[idate].split(',')[1], l[idate].split(',')[2]
    
    print('<><><><><> '+date+', '+station+', '+country+' <><><><><>')

    ## Then use Pandas to read in the sonde profile
    if footer>0:
        print('yes again')
        df = pd.read_csv(osf, header=header, skipfooter=len(l)-footer)
        #df = df[:footer]
    else:
        df = pd.read_csv(osf, header=header)
    df['lat'], df['lon'] = lat, lon ## Appending to the df makes the lat-lon the right dimension for .nc input
    try:
        df = df.astype('float', errors='ignore')
        df = df[df['Duration']>-0.0001] ## Remove "pre-sonde" information
        df['CT_sampling_strategy'] = 4  ## 4=instantaneous sampling in geoschem
        df['Date_Time'] = list(map(calc_DT, df['Duration']))
        df.set_index(pd.DatetimeIndex(df['Date_Time'], dayfirst=True), inplace=True)
        df['year'],df['month'],df['day'] = df.index.year,df.index.month,df.index.day   ## Get breakdown of time components as required for
        df['hour'],df['minute'],df['second'] = df.index.hour,df.index.minute,df.index.second
        
        ## Determine if ozonesonde goes over midnight or not
        if df.index.date[0] == df.index.date[-1]:
            print('--->> same start and end date, carry on :) ')
            dfl, dates = [df], [date]
        else:
            print('--->> goes over midnight, splitting by time! ')
            dfl = [df.loc[date+' 00:00':date+' 23:59'],
                   df.loc[str(df.index.date[-1])+' 00:00':str(df.index.date[-1])+' 23:59']]
            dates = [date, str(df.index.date[-1])]
        
        for j in range(len(dfl)):
        
            ## Put the datetime info into arrays suitable for appending to the .nc
            dtarray = [ [ dfl[j]['year'].values[i], dfl[j]['month'].values[i],  
             dfl[j]['day'].values[i], dfl[j]['hour'].values[i],
             dfl[j]['minute'].values[i], 
             dfl[j]['second'].values[i]] for i in range(len(dfl[j]))]
        
            ## netcdf file name
            netcdfname = 'obspack_O3sonde_'+dates[j]
            obs_id = netcdfname
        
            ## A 200-character array is required for obspack ID, achieve this with the created obs_id
            ##     and a masked array
            chars = []
            [chars.append(obs_id[i]) for i in range(len(obs_id)) ]
            chars_ = np.pad(chars, (0,200-len(obs_id)), 'constant', constant_values=('0'))
            chars_masked = np.ma.masked_where(chars_=='0', chars_)
        
            ## Check if a .nc file has already been made for this location and time
            if os.path.isfile(osp+netcdfname+'.nc'):
                print('-----------> '+obs_id+'.nc exists already, data will be appended')
            
                ncfile = Dataset(osp+netcdfname+'.nc', 'a')
                Time = ncfile.variables['time_components'][:]
                #Time = Time[~np.isnan(Time)]
                Lat,Lon,Alt  = ncfile.variables['latitude'][:],ncfile.variables['longitude'][:],ncfile.variables['altitude'][:]
                ncfile.close(); print('-----------> old file closed')
                os.remove(osp+netcdfname+'.nc') ; print('-----------> old file removed')
            
                # Make a transfer dataframe 'tdf'
                tdf = pd.DataFrame()
                tdf['Date_Time'] = [ calc_DT2(Time[i]) for i in range(len(Time)) ]  
                tdf['lat'], tdf['lon'], tdf['GPHeight'] = Lat,Lon,Alt
                tdf['CT_sampling_strategy'] = 4
                tdf.set_index(pd.DatetimeIndex(tdf['Date_Time'], dayfirst=True), inplace=True)
                tdf['year'],tdf['month'],tdf['day'] = tdf.index.year,tdf.index.month,tdf.index.day   ## Get breakdown of time components as required for
                tdf['hour'],tdf['minute'],tdf['second'] = tdf.index.hour,tdf.index.minute,tdf.index.second                                            
                tdf = pd.concat([dfl[j],tdf])
                tdf = tdf[tdf.index.notnull()]
                tdf.sort_values(by='Date_Time',inplace=True)
                
                ncfile = Dataset(osp+obs_id+'.nc',mode='w',format='NETCDF4_CLASSIC') 
                print('-----------> new file made')
                obs_dim = ncfile.createDimension('obs', None)
                cal_dim = ncfile.createDimension('calendar_components', 6)
                str_dim = ncfile.createDimension('string_of_200chars', 200)
                
                Lat = ncfile.createVariable('latitude', np.float32, ('obs',))
                Lat.units = 'degrees_north'
                Lat[:] = tdf['lat'].values

                Lon = ncfile.createVariable('longitude', np.float32, ('obs',))
                Lon.units = 'degrees_east'
                Lon[:] = tdf['lon'].values

                Alt = ncfile.createVariable('altitude', np.float32, ('obs',))
                Alt.units = 'meters'
                Alt.long_name = 'sample altitude in meters above sea level'
                Alt[:] = tdf['GPHeight'].values

                Time=ncfile.createVariable('time_components', np.int32, ('obs', 'calendar_components'))
                Time.fill_value = -9
                Time.long_name = "Calendar time components as integers.  Times and dates are UTC." 
                Time.order = "year, month, day, hour, minute, second" 
                Time.comment = "Calendar time components as integers.  Times and dates are UTC."
                Time[:] = [ [tdf['year'].values[i], tdf['month'].values[i],  
                           tdf['day'].values[i], tdf['hour'].values[i],
                           tdf['minute'].values[i], tdf['second'].values[i]] for i in range(len(tdf))] 

                ID = ncfile.createVariable('obspack_id', 'S1', ('obs', 'string_of_200chars'))
                ID.long_name = "Unique ObsPack observation id" ;
                ID.comment = "Unique observation id string that includes obs_id, dataset_id and obspack_num." 
                ID[:] = np.tile(chars_masked, (len(df),1))

                samp = ncfile.createVariable('CT_sampling_strategy', np.int32, ('obs',))
                samp.fill_value = -9
                samp.long_name = "model sampling strategy" ;
                samp.values = "How to sample model. 1=4-hour avg; 2=1-hour avg; 3=90-min avg; 4=instantaneous" ;
                samp[:] = tdf['CT_sampling_strategy'].values

                ncfile.close() ; print('-----------> new file  finished and closed')
            
            else: ## Create the netcdf file with required dimensions and variables

                ncfile = Dataset(osp+obs_id+'.nc',mode='w',format='NETCDF4_CLASSIC') 
                #print(ncfile)
                obs_dim = ncfile.createDimension('obs', None)
                cal_dim = ncfile.createDimension('calendar_components', 6)
                str_dim = ncfile.createDimension('string_of_200chars', 200)
                
                Lat = ncfile.createVariable('latitude', np.float32, ('obs',))
                Lat.units = 'degrees_north'
                Lat[:] = dfl[j]['lat'].values

                Lon = ncfile.createVariable('longitude', np.float32, ('obs',))
                Lon.units = 'degrees_east'
                Lon[:] = dfl[j]['lon'].values

                Alt = ncfile.createVariable('altitude', np.float32, ('obs',))
                Alt.units = 'meters'
                Alt.long_name = 'sample altitude in meters above sea level'
                Alt[:] = dfl[j]['GPHeight'].values

                Time=ncfile.createVariable('time_components', np.int32, ('obs', 'calendar_components'))
                Time.fill_value = -9
                Time.long_name = "Calendar time components as integers.  Times and dates are UTC." 
                Time.order = "year, month, day, hour, minute, second" 
                Time.comment = "Calendar time components as integers.  Times and dates are UTC."
                Time[:] = dtarray 

                ID = ncfile.createVariable('obspack_id', 'S1', ('obs', 'string_of_200chars'))
                ID.long_name = "Unique ObsPack observation id" ;
                ID.comment = "Unique observation id string that includes obs_id, dataset_id and obspack_num." 
                ID[:] = np.tile(chars_masked, (len(df),1))

                samp = ncfile.createVariable('CT_sampling_strategy', np.int32, ('obs',))
                samp.fill_value = -9
                samp.long_name = "model sampling strategy" ;
                samp.values = "How to sample model. 1=4-hour avg; 2=1-hour avg; 3=90-min avg; 4=instantaneous" ;
                samp[:] = df['CT_sampling_strategy'].values

                # first print the Dataset object to see what we've got
                #print(ncfile)
                # close the Dataset.
                ncfile.close(); print('-----------> ncfile made for '+date)
        
    except KeyError or TypeError:
        print(' -------X--->  Data incomplete/cannot be processed at '+station+', '+country)


#%%
