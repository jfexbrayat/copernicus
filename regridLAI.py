"""
01/05/2018 - JFE
This script includes a function to regrid Copernicus LAI data
to the desired spatial resolution.
"""

from netCDF4 import Dataset
import numpy as np
import glob, os

def regridLAI(path2orig,path2dest,latres,lonres,variables = ['LAI','LAI_ERR']):
    """
    This function regrids a Copernicus LAI file onto a regular grid with same 
    extent and resolution provided in the arguments.


    - path2orig : path to original file holding the data
    - path2dest : path to new file with regridded info
    - latres    : latitude resolution of target grid
    - lonres    : longitude resolution of target grid
    - variables : which variables to regrid
    """

    # open file
    nc = Dataset(path2orig)

    # get lat/lon from file
    origlat = nc.variables['lat'][:]
    origlon = nc.variables['lon'][:]

    # get resolution
    origlatres = np.abs(origlat[1]-origlat[0])
    origlonres = np.abs(origlon[1]-origlon[0])

    # adjust the lat/lon to centre of cell
    origlat = nc.variables['lat'][:]-origlatres/2.
    origlon = nc.variables['lon'][:]+origlonres/2.

    # define scanning window size
    sizelat = np.abs(np.round(latres/origlatres).astype('i'))
    sizelon = np.abs(np.round(lonres/origlonres).astype('i'))

    # define destination lat / lon arrays
    destlat = np.arange(80-latres/2.,-60,-latres)
    destlon = np.arange(-180+lonres/2.,180,lonres)

    #calculate grid cell area
    areas = np.zeros([origlat.size,origlon.size])
    for la,latval in enumerate(origlat):  
        areas[la]= (6371e3)**2 * ( np.deg2rad(0+origlonres/2.)-np.deg2rad(0-origlonres/2.) ) * (np.sin(np.deg2rad(latval+origlatres/2.))-np.sin(np.deg2rad(latval-origlatres/2.)))

    #create the destination file
    if len(glob.glob(path2dest)) > 0:
        os.remove(glob.glob(path2dest)[0])

    ncdest = Dataset(path2dest,'w')

    ncdest.createDimension('lon',size=destlon.size)
    ncdest.createDimension('lat',size=destlat.size)

    ncdest.createVariable('lat','d',dimensions=('lat'))
    ncdest.variables['lat'][:] = destlat
    ncdest.variables['lat'].units='degrees_north'

    ncdest.createVariable('lon','d',dimensions=('lon'))
    ncdest.variables['lon'][:] = destlon
    ncdest.variables['lon'].units='degrees_east'

    for varname in variables:
        ncdest.createVariable(varname,'d',dimensions=('lat','lon'))
        ncdest.variables[varname].missing_value = -9999.
        ncdest.variables[varname].long_name = nc.variables[varname].long_name
      

        # iterate grid
        target = np.zeros([destlat.size,destlon.size]) - 9999.
        for la, latval in enumerate(destlat):  
            for lo, lonval in enumerate(destlon):
                print la, lo
                slcarea = areas[(la*sizelat):((la+1)*sizelat),(lo*sizelon):((lo+1)*sizelon)]
                slcdata = nc.variables['LAI'][0,(la*sizelat):((la+1)*sizelat),(lo*sizelon):((lo+1)*sizelon)]
                if 'mask' in dir(slcdata):
                    if slcdata.mask.sum() <= 0.5*slcdata.size:
                        target[la,lo] = (slcdata*slcarea).sum()/(~slcdata.mask*slcarea).sum()
                else:
                    target[la,lo] = (slcdata*slcarea).sum()/slcarea.sum()
        ncdest.variables[varname][:] = target

    nc.close();ncdest.sync();ncdest.close()
    print "no problem till now"    
    return 0
    
if __name__ == "__main__":

    path2file = '/disk/scratch/local.2/jexbraya/dummy/copernicus/LAI_201801130000_GLOBE_PROBAV_V1.5.1/c_gls_LAI_201801130000_GLOBE_PROBAV_V1.5.1.nc'
    regridLAI(path2file,'test.nc',1,1)

