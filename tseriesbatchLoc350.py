#!/usr/bin/env python

import sys
import mylib

# New runs. <> is the abbrlist entry for the quantity of interest. 
Locinputs = {0:['run09_20110424_weimer','run09QuietWei110424_Loc1_<>_Alt350.txt','350.0'],\
          1:['run10_20110428_weimer','run10StormWei110428_Loc1_<>_Alt350.txt','350.0'],\
          2:['run12_20110424_amie','run12QuietAmi110424_Loc1_<>_Alt350.txt','350.0'],\
          3:['run13_20110428_amie','run13StormAmi110428_Loc1_<>_Alt350.txt','350.0'] }

varlist = ["'Electron Density'","'O Density'", "'N2 Density'", "'Neutral Velocity (north)'", "'Neutral Temperature'", "'Ion Velocity (up)'"]
abbrlist = ['eDens','ODens','N2dens','NW','Temp','IonVert']

# Real events
for i in Locinputs.keys():
    cstr = Locinputs[i]
    for j,var in enumerate(varlist):
       # Sub abbreviation into cstr[1] (Output file)
       outfile = cstr[1].replace('<>',abbrlist[j])
       cmd='tseriesLoc.py -i /tec/xingmeng/RunsGITM/RunsEvent/'+cstr[0]+'/IO2/3DUSR -o '+outfile+" -v "+var+" -lat 22.8 -lon 234.0 -a "+cstr[2]
       sys.stdout.write("Executing: \n"+cmd+"\n")
       #mylib.ex(cmd)
