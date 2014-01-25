#
# Copyright 2014, by the California Institute of Technology.
# ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged.
# Any commercial use must be negotiated with the
# Office of Technology Transfer at the California Institute of Technology.
#
# AJ Mannucci, JPL.

import copy

import gitm
import mylib
import numpy as np
import sys

def binDiff(bin1, bin2):
    '''Differences two GitmBin objects to produce difference
    values. Does not difference the following keys:
    dLon, Altitude, Longitude, LT, dLat, time, Latitude which are
    ASSUMED to be correctly the same for each bin object. 
    Inputs:
    bin1, bin2: the two objects
    Returns
    bin2-bin1
    '''
    exclude_list = ['dLon','Altitude','Longitude','LT','dLat',\
                    'time','Latitude']
    # Perform a shallow copy and then form differences of all
    # the keys. Perhaps a deep copy is better?
    bindiff = copy.copy(bin1)
    for ikey in bindiff.keys():
        if (not (ikey in exclude_list)):
            bindiff[ikey] = bin2[ikey] - bin1[ikey]
    return bindiff

def binReplace(bin, var, newval):
    '''binreplace will scale values in a binf files, for specified key.
    Does a deep copy (correct?)
    Inputs:
    bin: the gitm binary object.
    scale: multiply by this
    var: the text key variable to multiply.
    Returns:
    modified bin.
    '''
    bin[var] = copy.deepcopy(newval)
    return bin

def latAverage(binf, var, alt, latrange):
    '''Compute latitude-averaged quantity at a specific altitude
    or altitude integrated.
    Inputs:
    binf: the binary gitm object (GitmBin)
    var: the variable name
    alt: the altitude or interest. Function will compute the index for
      a nearby altitude. Use a negative value to obtain height-integrated.
      Meters.
    latrange: Lower lat boundary[0], upper lat boundary[1]. Radians. 
    Returns:
    [Scalar average, alt index used] alt index used if height integrated
    is -1.
    '''
    # Get the array. [lon,lat,alt]
    data = np.array(binf[var])
    nalt = len(data[0,0,:])
    # If a specific height, get the altitude index.
    # We assume altitude is independent of lon/lat.
    latval = binf['Latitude'][0,:,0]
    lonval = binf['Longitude'][:,0,0]
    indlat, = np.where((latval >= latrange[0]) & (latval < latrange[1]))
    # For some reason, certain longitudes are repeated in circle
    # sense. Restrict to 0->2pi
    indlon, = np.where((lonval >=0) & (lonval < 2.*np.pi))
    if (alt > 0):
        indalt = mylib.findNearest(alt, np.array(binf['Altitude'][0,0,:]))
        avg = 0.0
        for i in indlon:
            avg = avg + np.mean(data[i, indlat, indalt])
        return [avg/len(indlon), indalt]
    else: # Height integrated.
        avg = 0.0
        for i in indlon:
            for j in range(nalt):
                avg = avg + np.mean(data[i, indlat, j])
        return [avg/(len(indlon)*nalt),-1]

def latRegAverage(binf, var, alt, ut, latrange, lonrange, ltrange):
    '''Compute regionally-averaged quantity at a specific altitude
    or altitude integrated. Region is defined by longitude and local time (generally, these will
    be distinct regions).
    Inputs:
    binf: the binary gitm object (GitmBin)
    var: the variable name
    alt: the altitude or interest. Function will compute the index for
      a nearby altitude. Use a negative value to obtain height-integrated.
      Meters.
    ut: Univeral time of the observations. Secs past J2000.
    latrange: Lower lat boundary, upper lat boundary, 2 element list. Radians.
    lonrange: longitude range as for lat. Radians.
    ltrange: solar local time range. Hours.
    Returns:
    [Lon avg, LT avg, alt index used] alt index used. If height integrated
    is -1.
    '''
    # Get the array. [lon,lat,alt]
    data = np.array(binf[var])
    nalt = len(data[0,0,:])
    # If a specific height, get the altitude index.
    # We assume altitude is independent of lon/lat.
    latval = binf['Latitude'][0,:,0]
    lonval = binf['Longitude'][:,0,0]
    indlat, = np.where((latval >= latrange[0]) & (latval < latrange[1]))
    # Restrict in longitude. This requires a function call. 
    indlon, = mylib.lonindices(lonrange[0], lonrange[1], lonval)
    # Now determine the longitudes that fit the LT criterion.
    lonlt1 = mylib.lt2lon(ltrange[0], ut)
    lonlt2 = mylib.lt2lon(ltrange[1], ut)
    indlonlt, = mylib.lonindices(lonlt1, lonlt2, lonval)
    if (alt > 0):
        indalt = mylib.findNearest(alt, np.array(binf['Altitude'][0,0,:]))
        avg = 0.0
        for i in indlon:
            avg = avg + np.mean(data[i, indlat, indalt])
        avglt = 0.0
        for i in indlonlt:
            avglt = avglt + np.mean(data[i, indlat, indalt])
        return [avg/len(indlon), avglt/len(indlonlt), indalt]
    else: # Height integrated.
        avg = 0.0
        for i in indlon:
            for j in range(nalt):
                avg = avg + np.mean(data[i, indlat, j])
        avglt = 0.0
        for i in indlonlt:
            for j in range(nalt):
                avglt = avglt + np.mean(data[i, indlat, j])
        return [avg/(len(indlon)*nalt), avglt/(len(indlonlt)*nalt), -1]

def boundaryAvg(binf, var, alt, ut, latrange, lonrange, ltrange,mode='Lat'):
    '''Calculates a variable along a latitude (parallel) or longitude
    (meridian). For longitudes, can sub local time instead. 
    Returns the average along the boundary. Usually used for transport
    quantities, e.g. winds. 
    Inputs:
    binf: the binary gitm object (GitmBin)
    var: the variable name
    alt: the altitude or interest. Function will compute the index for
      a nearby altitude. Use a negative value to obtain height-integrated.
      Meters.
    ut: Univeral time of the observations. Secs past J2000.
    latrange: Lower, upper lat boundary. Radians
    lonrange: longitude range. Radians.
    ltrange: solar local time range in hours.
    mode: 'Lon', 'Lat', means that this coordinate is fixed. Only the first
      element of range is used for the fixed coord. Default is Lat. 
    Returns:
    [lon avg, lat avg, LT avg, alt index used] alt index used.
      If height integrated is -1.
      LT and Lon averages are only computed for mode "Lat", Lat avg
      only computed for mode "Lon".
    '''
    # Get the array. [lon,lat,alt]
    data = np.array(binf[var])
    altitudes = np.array(binf['Altitude'][0,0,:])
    nalt = len(altitudes)
    # If a specific height, get the altitude index.
    # We assume altitude is independent of lon/lat.
    latval = np.array(binf['Latitude'][0,:,0])
    lonval = np.array(binf['Longitude'][:,0,0])
    if (mode == 'Lat'):
        avglat = 0.0 # This is not computed. 
        # Get latitude index. 
        indlat = mylib.findNearest(latrange[0], latval)
        # Restrict in longitude. This requires a function call. 
        indlon = mylib.lonindices(lonrange[0]*180.0/np.pi, \
                                  lonrange[1]*180.0/np.pi, \
                                  lonval*180.0/np.pi)
        # Determine the longitudes that fit the LT criterion.
        lonlt1 = mylib.lt2lon(ltrange[0], ut)*np.pi/180.0
        if (lonlt1 < 0.0): lonlt1 = lonlt1 + 360.0
        lonlt2 = mylib.lt2lon(ltrange[1], ut)*np.pi/180.0
        if (lonlt2 < 0.0): lonlt2 = lonlt2 + 360.0
        indlonlt = mylib.lonindices(lonlt1*180.0/np.pi,
                                    lonlt2*180.0/np.pi, \
                                    lonval*180.0/np.pi)
        if (alt > 0):
            indalt = mylib.findNearest(alt, altitudes)
            avglon = np.mean(data[indlon, indlat, indalt])
            avglt = np.mean(data[indlonlt, indlat, indalt])
        else: # Height integrated
            avglon = 0.0
            for i in range(nalt):
                avglon = avglon + np.mean(data[indlon, indlat, i])
            avglon = avglon/float(nalt)
            avglt = 0.0
            for i in range(nalt):
                avglt = avglt + np.mean(data[indlon, indlat, i])
            avglt = avglt/float(nalt)
            indalt = -1
    elif (mode == 'Lon'):
        # For this case, there is no LT average calculated, or lon.
        avglt = 0.0
        avglon = 0.0
        # Average over a latitude range, for a single longitude.
        # Lon index
        indlon = mylib.findNearest(lonrange[0], \
                      np.array(binf['Longitude'][:,0,0]))
        # Find indices for latitude range.
        indlat = np.where((latval >= latrange[0]) & (latval < latrange[1]))
        if (alt > 0):
            indalt = mylib.findNearest(alt, \
                     np.array(binf['Altitude'][0,0,:]))
            avglat = np.mean(data[indlon, indlat, indalt])
        else: # Height integrated
            avglat = 0.0
            for i in range(nalt):
                avglat = avglat + np.mean(data[indlon, indlat, i])
            avglat = avglat/float(nalt)
            indalt = -1
    else:
        sys.stderr.write("ERROR: gitmLib:boundaryAvg: Incorrect mode passed in ("+mode+"). Allowed values are Lat or Lon. \n")
        sys.exit(1)
    return [avglon, avglat, avglt, indalt]
