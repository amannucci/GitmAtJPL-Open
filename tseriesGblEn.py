#!/usr/bin/env python
from __future__ import print_function
import gitmLib
import ionotime as IT
import sys
# Global average energy, according to Xing's email of 12/10/2015. 

def main():
    import argparse
    p = argparse.ArgumentParser(description='Calculates global energy for GITM. Writes to stdout.')
    p.add_argument('infiles', action='store',type=str, metavar='<Input wildcard>',help="Wildcard spec for retrieving times series of model files. All files with this wildcard are retrieved. The wild card is typically of the form: header*.bin for GITM and header*.nc for TIEGCM. Sort order of file names should also be time order. Enclose in quotes. Not to be expanded by shell. Use -G option in ipython.")
    p.add_argument('-v','--var',action='store',dest='var',nargs="+",default='Electron Density',type=str,metavar='<variable>',help="Variable name to extract. Use pullgitm.py to see variable names in the binary output. Enclose in quotes.")
    
    args = p.parse_args()
    var = args.var # Could be a list. 
    filespec = args.infiles

    nvar = len(var)
    
    mObj = gitmLib.ModelHandler(filespec)
#    print("#Variable: ",var)
    #Formatting stuff
    s = ""
    for i in var:
        i = i.replace(" ","_")
        s = s+i+"  "
    print("#Date                 Time2000        "+s)
    for t in mObj.modelTimes:
        date = IT.from_J2000(t,'YYYYMMDD_HHMMSS')
        a = []
        for v in var:
            a.append(gitmLib.gblAverage(mObj, t, v))
        sys.stdout.write("{:s}  {:15.1f}  ".format(date, t))
        for i in a:
            sys.stdout.write("{:15.4e}".format(i))
        sys.stdout.write("\n")
               
if __name__ == "__main__":
    main()
