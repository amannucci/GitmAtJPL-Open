#!/usr/bin/env python
# For agu Fall 2015 invited talk.
#
from __future__ import print_function
import gitmLib
import numpy as np
import mylib
import sys
import ionotime as IT

fname = 'agu2015/3DUSR_t120508_*.bin'
alt = 300.0
lat = 35.0
lon = 270.0

# Open the model.
mObj = gitmLib.ModelHandler(fname)

times = mObj.modelTimes

for (i,t) in enumerate(times):
    s = IT.from_J2000(t)
    datestr = s[0:4]+'-'+s[4:6]+'-'+s[6:8]
    l = IT.from_J2000(t,'LIST')
    ut = float(l[3])+(float(l[4])/60.0)
    slt = mylib.solar_local_time(t, lon)
    v = 'O+ Density'
    [oPlusDens, ia] = gitmLib.locValue(mObj, t, v, alt*1000.0, \
                                    lat*np.pi/180.0, \
                                    lon*np.pi/180.0)
    v = 'Electron Density'
    [elecDens, ia] = gitmLib.locValue(mObj, t, v, alt*1000.0, \
                                    lat*np.pi/180.0, \
                                    lon*np.pi/180.0)
    v = 'O Density'
    [ODens, ia] = gitmLib.locValue(mObj, t, v, alt*1000.0, \
                                    lat*np.pi/180.0, \
                                    lon*np.pi/180.0)

    # Calculate production based on the standard number that Xing provided:
    prodStd = ODens*3.0e-7

    # Compare O+ to Ne:

    # Get neutral temperature:
    v = 'Neutral Temperature'
    [neutTemp, ia] = gitmLib.locValue(mObj, t, v, alt*1000.0, \
                                    lat*np.pi/180.0, \
                                    lon*np.pi/180.0)

    # Calc O2 and N2 losses
    v = 'O2 Density'
    [o2Dens, ia] = gitmLib.locValue(mObj, t, v, alt*1000.0, \
                                    lat*np.pi/180.0, \
                                    lon*np.pi/180.0)
    v = 'N2 Density'
    [n2Dens, ia] = gitmLib.locValue(mObj, t, v, alt*1000.0, \
                                    lat*np.pi/180.0, \
                                    lon*np.pi/180.0)
    o2Loss = elecDens * o2Dens * \
             ( 2.82e-11 - 7.74e-12*neutTemp/300.0 \
             + 1.073e-12*(neutTemp/300.0)**2 \
             - 5.17e-14*(neutTemp/300.0)**3 \
             + 9.65e-16*(neutTemp/300.0)**4 ) * 1.e-6

    n2Loss = elecDens * n2Dens * ( 1.533e-12 - 5.92e-13*(neutTemp/300) \
              + 8.6e-14*(neutTemp/300)**2 )*1.0e-6

    # Flux estimation. Here we are interested in plasma flux. Consult
    # Mathematica notebook gitminterpretation.nb for the formula for the
    # divergence of the flux. Also see README in agu2015 area. 
    # The first
    # question is how to get neighboring voxels. See how locValue works.
    # It works with the more general model object, and can retrieve all
    # lats and lons.
    # Presumably we are dealing with a single time. Let's get it.

    latval = gitmLib.collapseDim(mObj.getVariable('Latitude',t), 1)
    lonval = gitmLib.collapseDim(mObj.getVariable('Longitude',t), 0)
    altval = gitmLib.collapseDim(mObj.getVariable('Altitude',t),2)
    # Lat/lon seem to be on a regular grid, but not altitude. However,
    # the ordering is trivial. So NN are just found by adding/subtracting
    # one from the index.
    indlat = mylib.findNearest(lat*np.pi/180.0, latval); lat0 = latval[indlat]
    indlon = mylib.findNearest(lon*np.pi/180.0, lonval); lon0 = lonval[indlon]
    indalt = mylib.findNearest(alt*1000.0, altval); alt0 = altval[indalt]
    # Now get the NN. This is a rectangular grid. Notation is key.
    latPl = latval[indlat+1]; latMi = latval[indlat-1]
    lonPl = lonval[indlon+1]; lonMi = lonval[indlon-1]
    altPl = altval[indalt+1]; altMi = altval[indalt-1]
    # Average centered values with NN values. I want values at each
    # face, and then difference across the voxel for approx. derivative.
    # The flux is Ne*bV. The quantities are E,N,V directions. That is:
    # del-lon, del-lat, del-r.
    # First fetch the arrays.
    EDensArr = mObj.getVariable('Electron Density', t)
    ionVnorth = mObj.getVariable('Ion Velocity (north)', t)
    ionVeast = mObj.getVariable('Ion Velocity (east)', t)
    ionVvert = mObj.getVariable('Ion Velocity (up)', t)
    # Compute flux
    fluxN = EDensArr*ionVnorth
    fluxE = EDensArr*ionVeast
    fluxV = EDensArr*ionVvert
    # North component of flux, compute face values.
    [fluxN_Eplfc, fluxN_Emifc, \
    fluxN_Nplfc, fluxN_Nmifc, \
    fluxN_Aplfc, fluxN_Amifc] = gitmLib.finDiff(fluxN, indlon, indlat, indalt)
    # East component of flux
    [fluxE_Eplfc, fluxE_Emifc, \
    fluxE_Nplfc, fluxE_Nmifc, \
    fluxE_Aplfc, fluxE_Amifc] = gitmLib.finDiff(fluxE, indlon, indlat, indalt)
    # Vertical component of flux
    [fluxV_Eplfc, fluxV_Emifc, \
    fluxV_Nplfc, fluxV_Nmifc, \
    fluxV_Aplfc, fluxV_Amifc] = gitmLib.finDiff(fluxV, indlon, indlat, indalt)
    # Compute voxel centered values
    fluxN0 = fluxN[indlon, indlat, indalt]
    fluxE0 = fluxE[indlon, indlat, indalt]
    fluxV0 = fluxV[indlon, indlat, indalt]

    # Now compute an approx to the divergence of the flux. Spherical coord
    # See Ma notebook. 
    # First I need to compute the voxel size in all three directions. Use
    # same approach as for the functions.
    latPlfc = (latPl - lat0)/2.0 + lat0
    latMifc = (lat0 - latMi)/2.0 + latMi
    lonPlfc = (lonPl - lon0)/2.0 + lon0
    lonMifc = (lon0 - lonMi)/2.0 + lonMi
    altPlfc = (altPl - alt0)/2.0 + alt0
    altMifc = (alt0 - altMi)/2.0 + altMi
    distLat = alt0*(latPlfc - latMifc) # Does not account for rad change across
    # voxel. Uses voxel central radius. 
    distLon = alt0*(lonPlfc - lonMifc)
    distAlt = altPlfc - altMifc
    # Compute needed derivatives.
    fluxEderE = (fluxE_Eplfc - fluxE_Emifc)/distLon
    fluxNderN = (fluxN_Nplfc - fluxN_Nmifc)/distLat
    fluxVderV = (fluxV_Aplfc - fluxV_Amifc)/distAlt
    # Now compose the entire expression. See Ma notebook.
    # Csc = 1/np.sin and Sec = 1/np.cos
    #
    div = (1.0/np.sin(lat0))*(1.0/alt0)*(fluxEderE + np.sin(lat0)*fluxV0 + \
          + np.cos(lat0)*fluxN0) + fluxVderV + \
          (1/alt0)*(fluxV0 + fluxNderN)
    # Let's also compute a cartesian approximation to the divergence
    divC = fluxEderE + fluxNderN + fluxVderV

    # Let's write out the terms to the stdout.
    print("#Date       Time  SLT   Edens    O+      Prod     o2Loss   n2Loss    Div    DivC", file=sys.stdout)
    print("{:s} {:5.2f} {:5.2f} {:8.2E} {:8.2E} {:8.2E} {:8.2E} {:8.2E} {:8.2E} {:8.2e}".format(datestr, ut, slt, elecDens, \
            oPlusDens, prodStd, o2Loss, n2Loss, div, divC))

