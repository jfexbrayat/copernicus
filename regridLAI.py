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
    - mask      : to only consider pixels according to a mask of e.g. land cover
    """

    if len(glob.glob(path2dest)) > 0:
        print 'Destination file "%s" already exists, please remove before proceeding' % (glob.glob(path2dest)[0])

    else:
        print "Regridding file " + path2orig.split('/')[-1]
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

        ncdest = Dataset(path2dest,'w')

        #create dimensions and corresponding variables
        ncdest.createDimension('lat',size=destlat.size)
        ncdest.createVariable('lat','d',dimensions=('lat'))
        ncdest.variables['lat'][:] = destlat
        ncdest.variables['lat'].units='degrees_north'
        ncdest.variables['lat'].latitude='latitude'

        ncdest.createDimension('lon',size=destlon.size)
        ncdest.createVariable('lon','d',dimensions=('lon'))
        ncdest.variables['lon'][:] = destlon
        ncdest.variables['lon'].units='degrees_east'
        ncdest.variables['lon'].latitude='longitude'

        #query for time as it is not in the 300m files
        if 'time' in nc.dimensions:
            ncdest.createDimension('time',size=None)
            ncdest.createVariable('time','d',dimensions=('time'))
            ncdest.variables['time'][:] = nc.variables['time'][:]
            ncdest.variables['time'].units = nc.variables['time'].units
            ncdest.variables['time'].long_name = nc.variables['time'].long_name

        ncdest.createVariable('fraction','d',dimensions=nc.variables[variables[0]].dimensions)
        ncdest.variables['fraction'].units = '%'
        ncdest.variables['fraction'].long_name = 'fraction of regridded cell which had data at original resolution'

        for va,varname in enumerate(variables):
            print "Regridding variable %s ... " % (varname)
            ncdest.createVariable(varname,'d',dimensions=nc.variables[varname].dimensions)
            ncdest.variables[varname].missing_value = -9999.
            ncdest.variables[varname].long_name = nc.variables[varname].long_name

            # iterate grid
            target = np.zeros([destlat.size,destlon.size]) - 9999.
            if va == 0:
                fraction = np.zeros(target.shape)-9999.
            counter = 0
            for la, latval in enumerate(destlat):
                for lo, lonval in enumerate(destlon):
                    counter+=1
                    print '\rRegridding pixel %i / %i' % (counter, len(destlat)*len(destlon)),
                    slcarea = areas[(la*sizelat):((la+1)*sizelat),(lo*sizelon):((lo+1)*sizelon)]
                    if 'time' in nc.dimensions:
                        slcdata = nc.variables[varname][0,(la*sizelat):((la+1)*sizelat),(lo*sizelon):((lo+1)*sizelon)]
                    else:
                        slcdata = nc.variables[varname][(la*sizelat):((la+1)*sizelat),(lo*sizelon):((lo+1)*sizelon)]
                    if 'mask' in dir(slcdata):
                        if slcdata.mask.sum() != slcdata.size:
                            target[la,lo] = (slcdata*slcarea).sum()/(~slcdata.mask*slcarea).sum()
                    # first pass, extract the fraction of pixel with valid data
                        if va == 0:
                          #  fraction[la,lo] = (~slcdata.mask).sum()/(float(slcdata.size))
                            fraction[la,lo] = (~slcdata.mask*slcarea).sum()/(slcarea.sum())
                    else:
                        target[la,lo] = (slcdata*slcarea).sum()/slcarea.sum()
                        if va == 0:
                            fraction[la,lo]=1.
            ncdest.variables[varname][:] = np.expand_dims(target,0)

            if va == 0:
                ncdest.variables['fraction'][:] = np.expand_dims(fraction*100.,0)
        ncdest.sync();ncdest.close()
        nc.close()
    return 0

if __name__ == "__main__":

    import sys
    from sklearn.externals.joblib import Parallel, delayed

    res = sys.argv[1]

    path2files = glob.glob('/disk/scratch/local.2/copernicus/LAI/*/*nc');path2files.sort()
   # path2files = path2files[:1]

  #  print len(path2files)

    path2dest = []
    for fname in path2files:
        dummy = fname.split('/')
        destfile = dummy[-1].split('.')[0]+'_%sx%s.nc' % (res,res)
        path2dest.append('/'.join(fname.split('/')[:5])+'/LAI_%sx%s/%s' % (res,res,destfile))

    Parallel(n_jobs = 10)(delayed(regridLAI)(src,target,float(res),float(res)) for src,target in zip(path2files,path2dest))
