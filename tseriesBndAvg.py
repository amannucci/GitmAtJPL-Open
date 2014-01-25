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

if (__name__ == '__main__'):
    import argparse
    p = argparse.ArgumentParser(description='Calculate time series of boundary variables from GITM output. A boundary is defined in terms of a single longitude and a latitude range, or a single latitude and a longitude and/or solar local time range. A mode input determines which boundary is calculated. Default mode is latitude boundary (Lat).')
    p.add_argument('-i','--infile',action='store',dest='infile',default=None,type=str,metavar='<Input header>',help="Header for retrieving times series of GITM files. All files with this header are retrieved. The wild card is: header*.bin. Sort order of names should also be time order. ")
    p.add_argument('-o','--outfile',action='store',dest='outfile',default=None,type=str,metavar='<output>',help="Name for output text file containing average values.")
    p.add_argument('-v','--var',action='store',dest='var',default=None,type=str,metavar='<variable>',help="Variable name to average. Use pullgitm.py to see variable names in the binary output. Enclose in quotes.")
    p.add_argument('-l','--latband',nargs=2, action='store',dest='lat',default=[60.0,90.0],type=float,metavar='<lat>',help="Latitude band to average over, lower, upper. Degrees. Only first value is used for mode Lat.")
    p.add_argument('--lonband',nargs=2, action='store',dest='lon',default=[235.0,290.0],type=float,metavar='<lon>',help="Longitude band to average over, lower, upper. Degrees. Only first value us used for mode Lon.")
    p.add_argument('--localt',nargs=2, action='store',dest='localt',default=None,type=float,metavar='<local-t>',help="Solar local time band to average over, lower, upper. Hours. Ignored for mode Lon.")
    p.add_argument('-a','--alt',action='store',dest='alt',default=350.0,type=float,metavar='<altitude>',help="Altitude of desired quantity. Use negative value for altitude average. km.")
    p.add_argument('-m','--mode',action='store',dest='mode',default='Lat',type=str,metavar='<mode>',help="Mode to use for the boundary. Use Lat or Lon. Lat means there is one latitude. Lon means there is one longitude. For Lon mode, the local times are ignored.")
    
    args = p.parse_args()
    if (args.infile == None):
        sys.stderr.write("ERROR: tseriesBndAvg.py: too few arguments. Need an input file name. Use -h for help.\n")
        sys.exit(1)
    
    flist = glob.glob(args.infile+'*.bin')
    if (len(flist) < 1):
        sys.stderr.write("ERROR: tseriesBndAvg.py: File not found: {:s}*.bin\n".format(args.infile))
        sys.exit(1)

    flist.sort() # Sort order is time order. 
    
    fout = open(args.outfile, 'w')
    fout.write("#Variable {:s} over {:6.2f}-{:6.2f} Lat at Altitude {:6.2f}\n".format(args.var, args.lat[0], args.lat[1], args.alt))
    fout.write("#Lon range: {:6.2f}-{:6.2f}      LT range: {:6.2f}-{:6.2f}\n")
    fout.write("#TimeJ2000       TimeStr                         LonAvg                 LTavg\n")
    for fs in flist:
        sys.stdout.write("tseriesRegAvg.py: Processing file "+fs+"\n")
        binf = gitm.GitmBin(fs)
        if (not (args.var in binf.keys())):
            sys.stderr.write("ERROR: tseriesRegAvg.py: Requested variable >>>"+var+"<<< not in file "+fs+". Aborting.\n")
            sys.exit(1)
        dt = binf['time'] # python datetime object
        secJ2000 = IT.to_J2000((dt.year, dt.month,\ # Seconds past J2000
                                dt.day, dt.hour, \
                                dt.minute, \
                                dt.second+\
                                1.0e-6*dt.microsecond))
        # Average over requested latitude and longitude or local time band.
        [avglon, avglt, ia] = gitmLib.RegAverage(binf, args.var, args.alt*1000.0, args.lat[0]*np.pi/180.0, args.lat[1]*np.pi/180.0,\
                                                 args.lon[0]*np.pi/180.0, args.lon[1]*np.pi/180.0, \
                                                 args.localt[0], args.localt[1])
        # Write to output file.
        dstr = dt.strftime('%Y-%m-%d/%H:%M:%S.%f')
        #
        fout.write("{:15.3f}  {:s}  {:12.5e}  {:12.5e}\n".format(secJ2000, dstr, avglon, avglt))
    fout.close()
