'''
>>> Code to utilise OTD-LIS lightning climatology for the purposes of comparison to 
         TROPOMI NO2 product
>>> Rob Ryan, UCL, April 2021

====================================================
GET SEASONAL CLIMATOLOGY FROM MONTHLY OTD-LIS Lightning
====================================================
climpath = path to otdlis lightning climatology, set by default to the path in UCL 
                  Myriad group directory
season = 'jja', 'djf, 'son' or 'mam' to indicate season, for the purposes of comparing to 
               TROPOMI seasonal averages. 
xres, yres = latitudinal and longitudinal resolution of the input array, allowing lightning
                   climatology to be regridded if necessary. 
Returns: nlat x nlon seasonally averaged lightning climatology array

In addition to standard imported packages Requires GDAL regridding code "regrid_gdal" 
   written by Rob Ryan
====================================================
'''

import numpy as np
from osgeo import gdal
from osgeo import gdal_array
from osgeo import osr
import netCDF4
import regrid_gdal as rg 

def get_seasonal_otdlis( season, xres, yres, climpath='/shared/ucl/depts/uptrop/Data/Satellites/' \
                                       +'OTDLIS/Climatology/LISOTD_HRMC_V2.3.2015.nc.nc4'
                                     ):
    
    if season == 'jja':
        months = [6,7,8]
    elif season == 'son':
        months = [9, 10, 11]
    elif season == 'djf':
        months = [1,11,12]
    elif season == 'mam':
        months = [3,4,5]
    else:
        print('   ** Error: season not properly defined ** ')

    s, lst = season, []
    
    d = netCDF4.Dataset(climpath)
    lat, lon = d.variables['Latitude'][:], d.variables['Longitude'][:]
    
    for m in months:
        mclim = rg.gdal_regrid(   np.array(d.variables['HRMC_COM_FR'][:,:, m-1]).transpose(), 
                                                lon, lat, 180/len(lat), 360/len(lon), xres, yres, plot=False  )

        lst.append( mclim[0] )
        
    return np.nanmean( lst, axis=0 )

'''
====================================================
Mask_lightning: masks an input array with a given threshold of climatological lightning
====================================================
arr = input array to be masked (numpy array)
method = '<' or '>', i.e. less than or greater than the threshold
threshold = threshold to mask by. 
                    Can be a float x, i.e. mask values < or > x
                    Can be 'mean', i.e. mask < or > mean
                    Can be 'percentile' i.e. mask < or > nth percentile
percentile = n if masking by nth percentile
climpath = path to otdlis lightning climatology, set by default to the path in UCL 
                  Myriad group directory

season = 'jja', 'djf, 'son' or 'mam' to indicate season, for the purposes of comparing to 
               TROPOMI seasonal averages. Passed to "get_seasonal_otdlis" function.
xres, yres = latitudinal and longitudinal resolution of the input array, allowing lightning
                   climatology to be regridded if necessary. Passed to "get_seasonal_otdlis" function.
Returns: original array where values not satisfying masking criteria are set to NaN
====================================================
'''
def mask_lightning( arr, method, threshold, season, xres, yres, percentile=90,
                                climpath='/shared/ucl/depts/uptrop/Data/Satellites/' \
                                       +'OTDLIS/Climatology/LISOTD_HRMC_V2.3.2015.nc.nc4'
                                 ):
    data = get_seasonal_otdlis(season, xres, yres )
    
    if threshold == 'mean':
        Th = np.nanmean( data )
    elif threshold == 'percentile':
        Th = np.percentile( data, percentile )
    else:
        Th = threshold
        
    if method == '<':
        M = np.ma.masked_where(data < Th, data )
    else:
        M = np.ma.masked_where(data > Th, data )
        
    arr[M.mask] = np.nan
    
    return arr