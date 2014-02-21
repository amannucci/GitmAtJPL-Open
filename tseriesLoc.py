#!/usr/bin/env python 
#Copyright 2014, by the California Institute of Technology. 
#ALL RIGHTS RESERVED. 
#United States Government Sponsorship acknowledged. 
#Any commercial use must be negotiated with the 
#Office of Technology Transfer at the California Institute of Technology.
# A J Mannucci, Jan 2014.

import gitm
import gitmLib # This is JPL
import glob
import numpy as np
import sys, os
import ionotime as IT

# Test command line
# tseriesLoc.py -i runDir/3DUSR_t110425 -o runQuietAmiLocAltInt.txt -v 'Electron Density' --lat 60.0 --lon 235.0 -a -1.
if (__name__ == '__main__'):
    import argparse
    p = argparse.ArgumentParser(description='Calculate time series of a single variable at a single location. Altitude is given, or altitude integrated. ')
    p.add_argument('-i','--infile',action='store',dest='infile',default=None,type=str,metavar='<Input header>',help="Header for retrieving times series of GITM files. All files with this header are retrieved. The wild card is: header*.bin. Sort order of names should also be time order. ")
    p.add_argument('-o','--outfile',action='store',dest='outfile',default=None,type=str,metavar='<output>',help="Name for output text file containing time series.")
    p.add_argument('-v','--var',action='store',dest='var',default=None,type=str,metavar='<variable>',help="Variable name to average. Use pullgitm.py to see variable names in the binary output. Enclose in quotes.")
    p.add_argument('--lat',action='store',dest='lat',default=[35.0],type=float,metavar='<lat>',help="Latitude. Degrees.")
    p.add_argument('--lon',action='store',dest='lon',default=[235.0],type=float,metavar='<lon>',help="Longitude. Degrees.")
    p.add_argument('-a','--alt',action='store',dest='alt',default=350.0,type=float,metavar='<altitude>',help="Altitude of desired quantity. Use negative value for altitude average. km.")
    
    args = p.parse_args()
    if (args.infile == None):
        sys.stderr.write("ERROR: tseriesLoc.py: too few arguments. Need an input file name. Use -h for help.\n")
        sys.exit(1)
    
    flist = glob.glob(args.infile+'*.bin')
    if (len(flist) < 1):
        sys.stderr.write("ERROR: tseriesLoc.py: File not found: {:s}*.bin\n".format(args.infile))
        sys.exit(1)

    flist.sort() # Sort order is time order. 
    
    fout = open(args.outfile, 'w')
    fout.write("#Variable {:s} at {:6.2f}  Longitude and {:6.2f} Latitude and Altitude {:6.2f}\n".format(args.var, args.lon, args.lat, args.alt))
    fout.write("#TimeJ2000       TimeStr                         LocVal\n")
    for fs in flist:
        sys.stdout.write("tseriesLoc.py: Processing file "+fs+"\n")
        binf = gitm.GitmBin(fs)
        if (not (args.var in binf.keys())):
            sys.stderr.write("ERROR: tseriesLoc.py: Requested variable >>>"+var+"<<< not in file "+fs+". Aborting.\n")
            sys.exit(1)
        dt = binf['time'] # python datetime object
        secJ2000 = IT.to_J2000((dt.year, dt.month,\
                                dt.day, dt.hour, \
                                dt.minute, \
                                dt.second+\
                                1.0e-6*dt.microsecond))
        # Pull out value.
        [locval, ia] = gitmLib.locValue(binf, args.var, args.alt*1000.0, \
                                        args.lat*np.pi/180.0, \
                                        args.lon*np.pi/180.0)
        # Write to output file.
        dstr = dt.strftime('%Y-%m-%d/%H:%M:%S.%f')
        #
        fout.write("{:15.3f}  {:s}    {:10.3e}\n".format(secJ2000, dstr, locval))
    fout.close()
