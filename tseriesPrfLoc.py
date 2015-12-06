#!/usr/bin/env python
#Copyright 2014, by the California Institute of Technology. 
#ALL RIGHTS RESERVED. 
#United States Government Sponsorship acknowledged. 
#Any commercial use must be negotiated with the 
#Office of Technology Transfer at the California Institute of Technology.
# A J Mannucci, Jan 2014.
# 4/21/14: Changes to the multimodel handler (ModelHandler class).
# 5/9/14: Parser bug, can't handle negative values for alt input.
# 1/25/15: Further development. Output the variable. 

from __future__ import print_function
import gitm
import gitmLib # This is JPL
import glob
import numpy as np
import sys, os
import ionotime as IT
from netCDF4 import Dataset

# Test command line
# tseriesPrfLoc.py -i "GitmRunDir/3DUSR_t110425*.bin" -o runODensQuietAmiPrfLoc.nc -v 'O Density' --lat 60.0 --lon 235.0.
# tseriesPrfLoc.py -i "../runs/apr2011storm_s008D110429_0100_110501_0000.nc" -o runODensTiegcmPrfLoc.nc -v 'O Density' --lat 60.0 --lon 235.0
# run -G -d -b tseriesPrfLoc.py -i "runDir2/3DUSR_t110429_*bin" -o run52OdensPrfLoc.nc -v 'O Density' --lat 60.0 --lon 235.0
if (__name__ == '__main__'):
    import argparse
    p = argparse.ArgumentParser(description='Calculate time series of a vertical profiles for a single location.')
    p.add_argument('-i','--infile',action='store',dest='infile',default=None,type=str,metavar='<Input wildcard>',help="Wildcard spec for retrieving times series of model files. All files with this wildcard are retrieved. The wild card is typically of the form: header*.bin for GITM and header*.nc for TIEGCM. Sort order of file names should also be time order. ")
    p.add_argument('-o','--outfile',action='store',dest='outfile',default=None,type=str,metavar='<output>',help="Name for output netCDF file containing time series.")
    p.add_argument('-v','--var',action='store',dest='var',default=None,type=str,metavar='<variable>',help="Variable name to extract. Use pullgitm.py to see variable names in the binary output. Enclose in quotes.")
    p.add_argument('--lat',action='store',dest='lat',default=[35.0],type=float,metavar='<lat>',help="Latitude. Degrees.")
    p.add_argument('--lon',action='store',dest='lon',default=[235.0],type=float,metavar='<lon>',help="Longitude. Degrees.")
#    p.add_argument('--mintime',action='store',dest='mintime',default=-1.0e36,type=float,metavar='<Earliest time>',help="Minimum time to return a result. J2000.")
#    p.add_argument('--maxtime',action='store',dest='maxtime',default=1.0e36,type=float,metavar='<Latest time>',help="Maximum time to return a result. J2000.")
    
    args = p.parse_args()
    outfile = args.outfile
    infile = args.infile
    if (infile == None):
        sys.stderr.write("ERROR: tseriesPrfLoc.py: too few arguments. Need an input file name. Use -h for help.\n")
        sys.exit(1)

    # Open the model.
    mObj = gitmLib.ModelHandler(infile)
    # Get all the times within
    times = np.array(mObj.modelTimes)
    ntimes = len(times)
    # Get all the altitudes within. Variable grids are [lon,lat,alt]
    altval = gitmLib.collapseDim(mObj.getVariable('Altitude',times[0]),2)
    nalt = len(altval)
    arr = np.zeros((ntimes,nalt))
    for (it,t) in enumerate(times): # Times are J2000
        sys.stdout.write("tseriesPrfLoc.py: Processing time "+IT.from_J2000(t)+"\n")
        for (ia,alt) in enumerate(altval): # Altitudes. 
            # Pull out value.
            [locval, ialt] = gitmLib.locValue(mObj, t, args.var, alt, \
                                            args.lat*np.pi/180.0, \
                                            args.lon*np.pi/180.0)
            # Fill array
            arr[it,ia] = locval
    # Write output to a file.
    try:
        
        root = Dataset(outfile, 'w', clobber=True)
    except:
        sys.stderr.write("FATAL ERROR: tseriesPrfLoc: error opening file "+outfile+". Terminating\n")
        sys.exit(1)
    root.description = "Produced by tseriesPrfLoc. Vertical profiles at a location versus time"
    root.latitude = args.lat
    root.longitude = args.lon
    root.physicalvar = args.var
    root.createDimension('time',ntimes)
    root.createDimension('altitude',nalt)
    tvar = root.createVariable('time','f8',('time',))
    avar = root.createVariable('altitude','f8',('altitude',))
    var = root.createVariable('variable','f8',('time','altitude'))
    tvar[:] = times.copy()
    avar[:] = altval.copy()
    var[:] = arr.copy()
    root.close()

#    fout.write("#Variable {:s} at {:6.2f}  Longitude and {:6.2f} Latitude and Altitude {:6.2f}\n".format(args.var, args.lon, args.lat, alt))
#    fout.write("#TimeJ2000       TimeStr              LocVal\n")
