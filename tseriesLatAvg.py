#!/usr/bin/env python 

import gitm
import gitmLib # This is JPL
import glob
import numpy as np
import sys, os
import ionotime as IT

if (__name__ == '__main__'):
    import argparse
    p = argparse.ArgumentParser(description='Calculate latitude-averaged time series from GITM output.')
    p.add_argument('-i','--infile',action='store',dest='infile',default=None,type=str,metavar='<Input header>',help="Header for retrieving times series of GITM files. All files with this header are retrieved. The wild card is: header*.bin. Sort order of names should also be time order. ")
    p.add_argument('-o','--outfile',action='store',dest='outfile',default=None,type=str,metavar='<output>',help="Name for output text file.")
    p.add_argument('-v','--var',action='store',dest='var',default=None,type=str,metavar='<variable to plot>',help="Variable name to average. Use pullgitm.py to see variable names in the binary output. Enclose in quotes.")
    p.add_argument('-l','--latband',nargs=2, action='store',dest='lat',default=[60.0,90.0],type=float,metavar='<lat>',help="Latitude band to average over, lower, upper. Degrees. ")
    p.add_argument('-a','--alt',action='store',dest='alt',default=350.0,type=float,metavar='<altitude>',help="Altitude of desired quantity. km.")
    
    args = p.parse_args()
    if (args.infile == None):
        sys.stdout.write("ERROR: tseriesLatAvg.py: too few arguments. Use -h for help.\n")
        sys.exit(1)
    
    flist = glob.glob(args.infile+'*.bin')
    if (len(flist) < 1):
        sys.stderr.write("ERROR: tseriesLatAvg.py: File not found: {:s}*.bin\n".format(args.infile))
        sys.exit(1)

    flist.sort() # Sort order is time order. 
    
    fout = open(args.outfile, 'w')
    fout.write("#Variable {:s} over {:6.2f}-{:6.2f} Lat at Altitude {:6.2f}\n".format(args.var, args.lat[0], args.lat[1], args.alt))
    fout.write("#TimeJ2000       TimeStr                         Value\n")
    for fs in flist:
        sys.stdout.write("tseriesLatAvg.py: Processing file "+fs+"\n")
        binf = gitm.GitmBin(fs)
        if (not (args.var in binf.keys())):
            sys.stderr.write("ERROR: tseriesLatAvg.py: Requested variable >>>"+var+"<<< not in file "+fs+". Aborting.\n")
            sys.exit(1)
        # Average over requested latitude band.
        [avg,ia] = gitmLib.latAverage(binf, args.var, args.alt*1000.0, args.lat[0]*np.pi/180.0, args.lat[1]*np.pi/180.0)
        # Write to output file.
        dt = binf['time']
        secJ2000 = IT.to_J2000((dt.year, dt.month,\
                                dt.day, dt.hour, \
                                dt.minute, \
                                dt.second+\
                                1.0e-6*dt.microsecond))
        dstr = dt.strftime('%Y-%m-%d/%H:%M:%S.%f')
        #
        fout.write("{:15.3f}  {:s}  {:12.5e}\n".format(secJ2000, dstr, avg))
    fout.close()
