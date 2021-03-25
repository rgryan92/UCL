import numpy as np
import xarray as xr
import regrid_gdal as rg

def apply_landmask( arr2mask, xres, yres, keep='ocean' ):
    path = '/shared/ucl/depts/uptrop/uclgcfs/gcgrid/gcdata/ExtData/HEMCO/'+\
           'OLSON_MAP/v2019-02/Olson_2001_Land_Type_Masks.025x025.generic.nc'
    landmask = xr.open_dataset(path)

    mLat, mLon, mask = landmask['lat'], landmask['lon'], landmask['LANDTYPE00'][0,:,:]

    mask_ =  rg.gdal_regrid(np.array(mask).transpose(), mLon, mLat, 0.25,0.25,
                        2., 2.5, plot=False )

    if keep == 'ocean':
        return np.ma.masked_where(mask_[0] < 0.9, arr2mask)
    else:
        return np.ma.masked_where(mask_[0] > 0.9, arr2mask)
