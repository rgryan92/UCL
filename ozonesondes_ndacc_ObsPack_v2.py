#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 14:39:48 2020

@author: rgryan

READS IN OZONE SONDE FILES DOWNLOADED FROM THE NDACC WEBSITE http://www.ndaccdemo.org/instruments/sonde
AND CONVERTS THEM TO INPUT NETCDF FILES FOR THE GEOS-CHEM OBSPACK DIAGNOSTIC
"""

import pandas as pd
import numpy as np
from netCDF4 import Dataset
import glob
import os

## PATH TO O3 SONDE FILES
path = '/Volumes/GoogleDrive/My Drive/Documents/Postdoc/rockets/'
osp = path+'ozonesondes/ndacc_2019_ozonesondes/'
outpath = path+'ozonesondes/obsPack_input/' # PATH FOR SAVED OUTPUT

# LIST OF NDACC STATIONS FOR WHICH YOU'RE CONVERTING
stations = ["dumont d'urville", "South Pole", "Belgrano", "Boulder", "Hilo",
             "Santa Cruz", "Natal Brazil" , "Neumayer", "Ny-Aalesund", 
             "Samoa", "Ittoqqortoormiit", "Wallops Flight Facility"]

## function to calculate datetime from sonde start datetime and 'duration'(seconds) field
def calc_DT(seconds):
        return startDT + pd.Timedelta(str(seconds)+' seconds')
def calc_DT2(dtarray):
    try:
        y,mo,d,h,mi,s = dtarray[0],dtarray[1],dtarray[2],dtarray[3],dtarray[4],dtarray[5]
        return pd.to_datetime(str(d)+'/'+str(mo)+'/'+str(y)+' '+str(h)+':'+str(mi)+':'+str(s), dayfirst=True)
    except ValueError:
        pd.to_datetime(np.NaN)
        
## USE "GLOB" TO READ IN ALL FILES IN YOUR OZONE SONDE DIRECTORY
##    AN OBSPACK OUTPUT FILE IS PRODUCED FOR EACH DAY (UTC TIME) 
##    ON WHICH THERE'S AN OZONE SONDE. SO, IF THE SONDE TIME GOES
##    OVER MIDNIGHT UTC, IT WILL BE SPLIT BETWEEN TWO OBSPACK 
##    FILES.

files = glob.glob(osp+'*.*')
successful, unsuccessful = [],[]
for osf in files[:]:
    #print(osf)
    
    ## Use 'Readlines' to find the header of the file and location lat/lon
    f = open(osf, "r")
    c,flag,footer,l =0, 0, 0, [] # Counters for finding information in file
    for line in f:
        l.append(line)  
        c = c+1
        try:
            if list(filter(None, line.split(' ')))[0] == 's':
                header=c-1
                flag=1
            elif list(filter(None, line.split(' ')))[0] == 'hPa':
                header=c-1
                flag=1
        except IndexError:
            continue
        
        if line[:-1] in stations:
            #print( '    Station is '+line[:-1] )
            station, iGeo = line[:-1], c
    
    if flag>0:
        successful.append(1)
        columns = list(filter(None, l[header-1].split(' ')))
        if columns[-1] == '\n':
            columns = columns[:-1]
        else:
            columns = columns
        #print(columns)
        
        top = list(filter(None, l[0].split(' ')))
        if top[-1] == '\n':
            top = top[:-1]
        else:
            top = top
        startDT = pd.to_datetime( top[-3:-2][0] + ' '+ top[-2:-1][0][:8],dayfirst=True) 
        endDT   = pd.to_datetime(top[-2:-1][0][8:]+' '+top[-1:][0][:8] ,dayfirst=True)
        gf = pd.read_csv(osf, delim_whitespace=True, header=header)
        gf.columns = columns
        
        df = pd.DataFrame()
        if 'ElapTime' in columns:
            df['Duration'] = gf['ElapTime']
        else:
            df['Duration'] = gf['Time']
        
        if 'Lat' and 'Lon' in columns:
            df['lat'], df['lon'] = gf['Lat'], gf['Lon']
        elif 'GPSLon' and 'GPSLat' in columns:
            df['lat'], df['lon'] = gf['GPSLat'], gf['GPSLon']
        else:
            geo = list(filter(None, l[iGeo].split(' ')))
            df['lat'], df['lon'] = float(geo[4]), float(geo[3])
            
        if 'Alt' in columns:
            df['alt'] = gf['Alt']
        elif 'GPSAlt' in columns:
            df['alt'] = gf['GPSAlt']
        elif 'GHgt' in columns:
            df['alt'] = gf['GHgt']
        else:
            df['alt'] = gf['GPSHgt']
            
        if df['alt'][len(df)-1] < 100: # then height is most likely in km...
            df['alt'] = 1000*df['alt']
            
        df['Date_Time'] = list(map(calc_DT, df['Duration']))
        df.set_index(pd.DatetimeIndex(df['Date_Time'], dayfirst=True), inplace=True)
        df['year'],df['month'],df['day'] = df.index.year,df.index.month,df.index.day   ## Get breakdown of time components as required for
        df['hour'],df['minute'],df['second'] = df.index.hour,df.index.minute,df.index.second
        df['CT_sampling_strategy'] = 4
        
        ##########
        ## Determine if ozonesonde goes over midnight or not
        
        if top[-3:-2][0] == top[-2:-1][0][8:]:
            print('    -√--> Same start & end day')
            dfl, dates = [df], [top[-3:-2][0]]
        else:
            print('    -X--> goes over midnight')
            dfl = [df.loc[top[-3:-2][0]+' 00:00':top[-3:-2][0]+' 23:59'],
                   df.loc[top[-2:-1][0][8:]+' 00:00':top[-2:-1][0][8:]+' 23:59']]
            dates = [top[-3:-2][0], top[-2:-1][0][8:]]
        
        for j in range(len(dfl)):
        
            ## Put the datetime info into arrays suitable for appending to the .nc
            dtarray = [ [ dfl[j]['year'].values[i], dfl[j]['month'].values[i],  
             dfl[j]['day'].values[i], dfl[j]['hour'].values[i],
             dfl[j]['minute'].values[i], 
             dfl[j]['second'].values[i]] for i in range(len(dfl[j]))]
            
            if len(str(dtarray[0][1])) < 2:
                m_ = '0'+str(dtarray[0][1])
            else:
                m_ = str(dtarray[0][1])
            if len(str(dtarray[0][2])) < 2:
                d_ = '0'+str(dtarray[0][2])
            else:
                d_ = str(dtarray[0][2]) 
            savedate = str(dtarray[0][0])+'-'+m_+'-'+d_
            ## netcdf file name
            netcdfname = 'obspack_O3sonde_'+savedate
            obs_id = netcdfname
        
            ## A 200-character array is required for obspack ID, achieve this with the created obs_id
            ##     and a masked array
            chars = []
            [chars.append(obs_id[i]) for i in range(len(obs_id)) ]
            chars_ = np.pad(chars, (0,200-len(obs_id)), 'constant', constant_values=('0'))
            chars_masked = np.ma.masked_where(chars_=='0', chars_)
        
            ## Check if a .nc file has already been made for this location and time
            if os.path.isfile(outpath+netcdfname+'.nc'):
                print('-----------> '+obs_id+'.nc exists already, data will be appended')
            
                ncfile = Dataset(outpath+netcdfname+'.nc', 'a')
                Time = ncfile.variables['time_components'][:]
                #Time = Time[~np.isnan(Time)]
                Lat,Lon,Alt  = ncfile.variables['latitude'][:],ncfile.variables['longitude'][:],ncfile.variables['altitude'][:]
                ncfile.close(); print('-----------> old file closed')
                os.remove(outpath+netcdfname+'.nc') ; print('-----------> old file removed')
            
                # Make a transfer dataframe 'tdf'
                tdf = pd.DataFrame()
                tdf['Date_Time'] = [ calc_DT2(Time[i]) for i in range(len(Time)) ]  
                tdf['lat'], tdf['lon'], tdf['alt'] = Lat,Lon,Alt
                tdf['CT_sampling_strategy'] = 4
                tdf.set_index(pd.DatetimeIndex(tdf['Date_Time'], dayfirst=True), inplace=True)
                tdf['year'],tdf['month'],tdf['day'] = tdf.index.year,tdf.index.month,tdf.index.day   ## Get breakdown of time components as required for
                tdf['hour'],tdf['minute'],tdf['second'] = tdf.index.hour,tdf.index.minute,tdf.index.second                                            
                tdf = pd.concat([dfl[j],tdf])
                tdf = tdf[tdf.index.notnull()]
                tdf.sort_values(by='Date_Time',inplace=True)
                
                ncfile = Dataset(outpath+obs_id+'.nc',mode='w',format='NETCDF4_CLASSIC') 
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
                Alt[:] = tdf['alt'].values

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

                ncfile = Dataset(outpath+obs_id+'.nc',mode='w',format='NETCDF4_CLASSIC') 
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
                Alt[:] = dfl[j]['alt'].values

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
                ncfile.close(); print('-----------> ncfile made for '+savedate)
        ##########
        
        print('----> Yay! Successful for '+osf[len(osp):])
    else:
        unsuccessful.append(1)
        print('--X-> Boo. Unsuccessful for '+osf[len(osp):])

