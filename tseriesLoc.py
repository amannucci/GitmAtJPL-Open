#!/usr/bin/env python 
#Copyright 2014, by the California Institute of Technology. 
#ALL RIGHTS RESERVED. 
#United States Government Sponsorship acknowledged. 
#Any commercial use must be negotiated with the 
#Office of Technology Transfer at the California Institute of Technology.
# A J Mannucci, Jan 2014.
# 4/21/14: Changes to the multimodel handler (ModelHandler class).
# 5/9/14: Parser bug, can't handle negative values for alt input.

import gitm
import gitmLib # This is JPL
import glob
import numpy as np
import sys, os
import ionotime as IT

# Test command line
# tseriesLoc.py -i "GitmRunDir/3DUSR_t110425*.bin" -o runQuietAmiLocAltInt.txt -v 'Electron Density' --lat 60.0 --lon 235.0 -a -1.
# tseriesLoc.py -i "GitmRunDir/3DUSR_t110425*.bin" -o runQuietAmiLocAltInt.txt -v 'Electron Density' --lat 60.0 --lon 235.0 -a -1.
# tseriesLoc.py -i "../runs/apr2011storm_s008D110429_0100_110501_0000.nc" -o runTiegcmLocAltInt.txt -v 'Electron Density' --lat 60.0 --lon 235.0 -a -1.
# tseriesLoc.py -i "GitmRunDir/3DUSR_t110425*.bin" -o runODensQuietAmiLocAltInt.txt -v 'O Density' --lat 60.0 --lon 235.0 -a -1.
# tseriesLoc.py -i "../runs/apr2011storm_s008D110429_0100_110501_0000.nc" -o runODensTiegcmLocAltInt.txt -v 'O Density' --lat 60.0 --lon 235.0 -a -1.
if (__name__ == '__main__'):
    import argparse
    p = argparse.ArgumentParser(description='Calculate time series of a single variable at a single location. Altitude is given, or altitude integrated. ')
    p.add_argument('-i','--infile',action='store',dest='infile',default=None,type=str,metavar='<Input wildcard>',help="Wildcard spec for retrieving times series of model files. All files with this wildcard are retrieved. The wild card is typically of the form: header*.bin for GITM and header*.nc for TIEGCM. Sort order of file names should also be time order. ")
    p.add_argument('-o','--outfile',action='store',dest='outfile',default=None,type=str,metavar='<output>',help="Name for output text file containing time series.")
    p.add_argument('-v','--var',action='store',dest='var',default=None,type=str,metavar='<variable>',help="Variable name to average. Use pullgitm.py to see variable names in the binary output. Enclose in quotes.")
    p.add_argument('--lat',action='store',dest='lat',default=[35.0],type=float,metavar='<lat>',help="Latitude. Degrees.")
    p.add_argument('--lon',action='store',dest='lon',default=[235.0],type=float,metavar='<lon>',help="Longitude. Degrees.")
    p.add_argument('-a','--alt',action='store',dest='alt',default=350.0,type=float,metavar='<altitude>',help="Altitude of desired quantity. Use negative value for altitude average. km.")
    p.add_argument('--mintime',action='store',dest='mintime',default=-1.0e36,type=float,metavar='<Earliest time>',help="Minimum time to return a result. J2000.")
    p.add_argument('--maxtime',action='store',dest='maxtime',default=1.0e36,type=float,metavar='<Latest time>',help="Maximum time to return a result. J2000.")
    
    args = p.parse_args()
    alt = args.alt
    if (args.infile == None):
        sys.stderr.write("ERROR: tseriesLoc.py: too few arguments. Need an input file name. Use -h for help.\n")
        sys.exit(1)

    fout = open(args.outfile, 'w')
    fout.write("#Variable {:s} at {:6.2f}  Longitude and {:6.2f} Latitude and Altitude {:6.2f}\n".format(args.var, args.lon, args.lat, alt))
    fout.write("#TimeJ2000       TimeStr              LocVal\n")

    # Open the model.
    mObj = gitmLib.ModelHandler(args.infile)
    # Get all the times within
    times = np.array(mObj.modelTimes)
    for t in times: # Times are J2000
        sys.stdout.write("tseriesLoc.py: Processing time "+IT.from_J2000(t)+"\n")
        # Pull out value.
        [locval, ia] = gitmLib.locValue(mObj, t, args.var, alt*1000.0, \
                                        args.lat*np.pi/180.0, \
                                        args.lon*np.pi/180.0)
        # Write to output file.
        dstr = IT.from_J2000(t, "YYYYMMDD_HHMMSS")
        #
        fout.write("{0:15.3f}  {1:s}    {2:10.5g}\n".format(t, dstr, float(locval)))
    fout.close()
