#Copyright 2014, by the California Institute of Technology. 
#ALL RIGHTS RESERVED. 
#United States Government Sponsorship acknowledged. 
#Any commercial use must be negotiated with the 
#Office of Technology Transfer at the California Institute of Technology.
# A J Mannucci.
from __future__ import print_function
import datetime as DT
import os
import numpy as np
import ionotime as IT
import sys
from datetime import datetime
import tomolib
import shutil

def dt2j2000(dt):
    """
    This function takes in a python datetime object and computes
    seconds past J2000 (Jan 1, 2000 12UT) for that object.
    To display this message, type help(mylib.dt2j2000).
    """
    return IT.to_J2000((dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second))

def doy2dt(year, doy):
    """
    This function takes in a year and doy and returns
    corresponding python datetime object.
    """
    ref = DT.datetime(year, 1, 1)
    refdays = ref.toordinal() - 1
    nowdays = refdays + doy
    return DT.datetime.fromordinal(nowdays)

def sec2cal(tsec): # OBSOLETE: use ionotime instead (from_J2000)
    """
    This function takes in seconds past J2000 and
    returns the ut hr for the day. It only works past
    Jan 1, 2000. It does not account for leap seconds and
    so is approximate. See ionotime library for more
    precise data. Based on routine sec2uthrfrac.m (matlab)
    """
    tadj = tsec - 43200.0 # Adj to midnight Jan 2, 2000
    secsday = np.mod(tadj, 86400.0) # get remainder to days
    hrs = (secsday - np.mod(secsday, 3600.0))/3600.0 # convt to hours
    mins = np.floor(np.mod(secsday, 3600.0)/60.0) # minutes in hour
    uthrfrac = hrs + (mins/60.0) # hour fraction
    return uthrfrac

def floatOrMiss(val, missingValue=-999.):
    try: val = float(val)
    except: val = missingValue
    return val

# This version will check for PRN number and return that as a float.
def floatOrMissPRN(val, missingValue=-999.):
    try: val = float(val)
    except: 
        if (val[0:3] == 'GPS'):
            val = float(val[3:5])
        else:
            val = missingValue
    return val

def dataLines(file, commentChar='#'):
    """Remove comment lines (beginning with #) from a text data file."""
    for line in open(file, 'r'):
        li = line.strip()
        if len(li) > 0:
            if li[0] != commentChar: yield line
        else:
            yield line

def dataLinesX(file, commentChar='#', filter=[False,0]):
    """Remove comment lines (beginning with #) from a text data file.
    Special version: only return lines where the column filter[0] is
    equal to the value filter[1]. Columns start at 1."""
    for line in open(file, 'r'):
        li = line.strip()
        if len(li) > 0:
            if ((li[0] != commentChar) and filter[0]):
                if (float(li.split()[filter[0]-1]) == filter[1]):
                    yield line
        else:
            yield line

# Luke's access to the command line
# Examples:
# print ex("cat sample.file | cl 1 2 | head -20")
# (prints the output as a single string)
# 
# print eex("cat sample.file | cl 1 2 | head -20")
# (prints the command, then the output, as a single string)
# 
#print exa("   ")
#(prints the output as an ARRAY if you want to manipulate it)
# 
#standard usage:
# 
#for line in exa("commands"):
#    (contents)

def  ex    (command      ): return "".join(os.popen(command).readlines())
def exa    (command      ): return os.popen(command).readlines()
def eex    (command      ): return command+"\n"+ex(command)

# Function that uses linear fitting to extend an array to
# an arbitrary point.
def extend(x, y, xval):
    '''Function that uses linear fitting to extend an array to
    an arbitrary point. Inputs:
    x: x values. Fit a line to these.
    y: y values (same size as x).
    xval: produce a y value at this x value. This has to be outside the
    range of x.
    '''
    if ((xval >= x[0]) and (xval <= x[-1])):
        sys.stderr.write('mylib.extend: ERROR: xval not outside range of x. Returning -9999.0\n')
        return -9999.0
    pl = np.polyfit(x,y,1)
    return np.polyval(pl, xval)

# Function that sets value of first point in array
# to first non-bad value in the array.
def firstgood(arr, dirflag, bad_value):
    '''firstgood: sets first values in array to first non-bad value
    in the array. Returns a new array.
    arr: the array that will the data.
    dirflag: 0 if bad_value is a lower bound, 1 if upper.
    bad_value: the bad value.
    '''
    if (dirflag == 0):
        if (arr[0] < (bad_value + 0.01)):
            barr = arr
            ind, = np.where(barr > bad_value)
            barr[0:ind[0]] = barr[ind[0]]
            return barr
        else:
            return arr
    if (dirflag == 1):
        if (arr[0] > (bad_value - 0.01)):
            barr = arr
            ind, = np.where(barr < bad_value)
            barr[0:ind[0]] = barr[ind[0]]
            return barr
        else:
            return arr

# Function that removes bad values from an array.
def clean_arrays(badValue, dirflag, refarray, *args):
    '''clean_arrays will remove bad values from the passed in arrays.
    Inputs:
    badValue: value to compare. "Bad" means larger or smaller than this value.
    dirflag:  0 if bad_value is a lower bound, 1 if upper. I.e. if 0, then
    good values are all larger than bad_value. If 1, then good values
    are smaller than bad_value.
    refarray: The array containing the bad values.
    *args: additional arrays of the same size as refarray that will also
    have the bad values removed according to criteria applied to refarry.
    Returns:
    list of arrays that are cleaned up according to criterion applied to
    refarray.
    '''
    narrays = len(args)
    # Get the indices that define good data. 
    if (dirflag == 0):
        ind, = np.where(refarray > badValue)
    if (dirflag == 1):
        ind, = np.where(refarray < badValue)
    refarray = refarray[ind]
    returnlist = [refarray]
    for i in range(len(args)):
        returnlist.append(args[i][ind])
    return returnlist

# Functions that sets x and y tick sizes with axes objects
def set_xlabel_size(ax, fsize):
    ''' Uses axes objects to set font size of labels.
    Replaces PLT.xticks(fontsize = fsize)
    ax: axis object
    fsize: numerical font size in pts.
    '''
#    print('Inside set xlabel size, ax = ', ax)
#    print('Inside set xlabel the labels are: ', ax.get_xticklabels())
    tl = ax.get_xticklabels()
#    print('The first one is: ', tl[0])
    tt = [str(i.get_text()) for i in ax.get_xticklabels()]
#    print('set xlabels', tt)
    ax.set_xticklabels(tt, fontsize = fsize)

def set_ylabel_size(ax, fsize):
    ''' Uses axes objects to set font size of labels.
    Replaces PLT.xlabels(fontsize = fsize)
    ax: axis object
    fsize: numerical font size in pts.
    '''
    tt = [str(i.get_text()) for i in ax.get_yticklabels()]
    #    print('set ylabels', tt)
    ax.set_yticklabels(tt, fontsize = fsize)

# list2grid is a miraculous function that reads in three 1-d arrays
# and creates three 2d arrays suitable for image plotting with contourf
# or contour.
# The technique is first done in plotoccplane2.py
#
def list2grid(x, y, z):
    '''Function that reads in a 2d data set that is linearly
    defined, and creates three 2D data sets that can be used
    in contourf, or Basemap.contourf.
    x: x values, 1D
    y: y values, 1D
    z: z values, 1D

    Returns the tuple:
    (X, Y, Z) where these are the desired 2D arrays.
    '''
    # Create unique sorted versions of the two dimensions
    hkeys = np.sort(np.unique(x))
    vkeys = np.sort(np.unique(y))

    # Create the mappings between the hdim, vdim values and
    # their place in the keys arrays
    hkeydict = {}
    for i in np.arange(len(hkeys)):
        hkeydict[hkeys[i]] = i
    vkeydict = {}
    for i in np.arange(len(vkeys)):
        vkeydict[vkeys[i]] = i

    # Use meshgrid to create the appropriate data structures
    # for contourf
    X,Y = np.meshgrid(hkeys, vkeys)

    # Create the z-array for plotting in contourf.
    # Note "transposing" of v and h
    Z = np.zeros((len(vkeys), len(hkeys)), dtype = 'f8')
    # The x,y arrays are in "visual" format.
    # The number of rows is the number of elements of vkeys (rows are vertical).
    # The number columns is the number of elements of hkeys (columns are horiz). 
    for i in np.arange(len(z)):
        i1 = hkeydict[x[i]]
        i2 = vkeydict[y[i]]
        Z[i2, i1] = z[i]
    return (X, Y, Z)

# Function that detrends an evenly space array. Originally from retr.py
def detrend_quadratic(y):
    '''A function that acts like the mlab functions detrend_none and
    detrend_linear, that returns residuals to a quadratic fit to the
    input data. Assumes that x values are evenly spaced. 
    '''
    # Create x values.
    x = np.arange(0.0,len(y))
    p = np.polyfit(x,y,2)
    resid = y - np.polyval(p, x)
    return resid
# Cubic detrend.
def detrend_cubic(y):
    '''A function that acts like the mlab functions detrend_none and
    detrend_linear, that returns residuals to a cubic fit to the
    input data. Assumes that x values are evenly spaced. 
    '''
    # Create x values.
    x = np.arange(0.0,len(y))
    p = np.polyfit(x,y,3)
    resid = y - np.polyval(p, x)
    return resid
# Order supplied detrend.
def detrend_poly(y,order):
    '''A function that acts like the mlab functions detrend_none and
    detrend_linear, that returns residuals to a nth-order fit to the
    input data. Assumes that x values are evenly spaced.
    Inputs:
    y: the data series at regular intervals
    order: order of polynomial to detrend. For linear, this = 1.
    '''
    # Create x values.
    x = np.arange(0.0,len(y))
    p = np.polyfit(x,y,3)
    resid = y - np.polyval(p, x)
    return resid

#Data limiter. For interpolation.
def data_limiter(time_inner, data_series, time_outer):
    '''
    data_limiter will condition the time_inner and data_series so that
    it can always be interpolated using time_outer. time_outer
    will set the upper and lower range for time_inner and the corresponding
    data series. time_outer will
    be adjusted to be monotonic increasing or decreasing according to
    the behaviour of time_inner, then clipping is performed.
    No extrapolation is done, so that depending on cadences, time_inner
    will cover a smaller time range than time_outer.
    Note that "monotonicity" is determined only by the first and last
    members of the time_inner series.
    Inputs:
    time_inner: the time series that may extend past time_outer.
    data_series: data associated with time_inner.
    time_outer: time series that defines the range for the return time
    series and data series.
    Outputs:
    [clipped_time: to fit within time_outer,
    clipped_data: associated with clipped_time (from data_series),
    ind_lo, ind_hi: the indices that are used to clip the input time series]
    '''
    # Determine increasing or decreasing.
    if (len(time_inner) == 0):
        print('ERROR in data_limiter call (mylib). Empty time_inner.', file=sys.stderr)
        sys.exit(1)
    if (time_inner[-1] > time_inner[0]): # Increasing.
        inc_flag = True
        if (time_outer[-1] < time_outer[0]): # Flip it.
            time_outer = np.flipud(time_outer)
    else: # decreasing
        inc_flag = False
        if (time_outer[-1] > time_outer[0]): # Flip it
            time_outer = np.flipud(time_outer)
    if (inc_flag): # time_inner is increasing
        # Check limit at starting indices
        ind, = np.where(time_inner < time_outer[0])
        if (len(ind) == 0):
            # Null set.
            ind_lo = 0
        else:
            ind_lo = ind[-1]+1
        # Check limit at ending indices
        ind, = np.where(time_inner > time_outer[-1])
        if (len(ind) == 0):
            # Null set
            ind_hi = len(time_inner)-1
        else:
            ind_hi = ind[0]-1
    else: # time_inner is decreasing
        # Check limit at starting indices
        ind, = np.where(time_inner > time_outer[0])
        if (len(ind) == 0):
            ind_lo = 0
        else:
            ind_lo = ind[-1]+1
        # Ending indices
        ind, = np.where(time_inner < time_outer[-1])
        if (len(ind) == 0):
            ind_hi = len(time_inner)-1
        else:
            ind_hi = ind[0]-1
    return [time_inner[ind_lo:ind_hi+1], data_series[ind_lo:ind_hi+1], ind_lo, ind_hi]

def flag_repeats(data, prec):
    '''This function returns the indices of repeated data at the beginning
    of a time series. It returns the index that signals where the data
    have stopped repeating.
    Input: some data series (1d numpy array).
    prec: precision to use to detect repeats. Could be zero.
    Returns: index at which non-repeated data starts
    '''
    ind, = np.where(np.abs(np.diff(data)) > prec)
    return ind[0]

def solar_local_time(ut, lon):
    '''This function takes in UTC and longitude and converts it to
    solar local time hours.
    ut: Univeral time in seconds past J2000
    lon: longitude in degrees, 0-360, or =-180.0 -> 180.0
    Returns:
    Solar local time in decimal hours.
    '''
    if (lon < 0.0): lon = lon + 360.0
    (yr,mm,dd,hh,mm,ss) = IT.from_J2000(ut, "LIST")
    uthrs = hh + (mm/60.0) + (ss/3600.0)
    x = (lon/15.0) + (uthrs)
    if (x > 24.0): x = x - 24.0
    return x

def lt2lon(slt, ut):
    '''This function takes in local_time (solar, hours) and UT (hours)
    and converts to longitude, from 0.0 to 360.0
    Inputs: solar_local_time (hours)
            ut in seconds past J2000
    Returns:
      Longitude (degrees) 0-360
    '''
    (yy,mo,dd,hh,mn,sec) = IT.from_J2000(ut,'LIST')
    uthr = hh+(mn/60.0)+(sec/3600.0)
    lon = 15.0*(slt - uthr)
    if (lon < 0.0): lon = lon + 360.0
    return lon
    
def lon_average(lon1, lon2):
    '''The lon_average function will take two equal length longitude
    time series and form an average longitude, taking into account
    that there could be a boundary crossing. Longitude convention is
    -180 to 180.0. Degress assumed. 
    Inputs: lon1 and lon2, two longitude numpy vectors, -180 to 180
    Returns: average lon.
    '''
    # If the difference between the two is >= 180.0, we need to
    # reset.
    avg = np.zeros(len(lon1))
    ind, = np.where(np.abs(lon1 - lon2) >= 180.0)
    for i in range(len(ind)):
        j = ind[i]
        avg[j] = (lon2[j] + 360.0 + lon1[j])/2.0
        if (avg[j] > 360.0):
            avg[j] = avg[j] - 360.0
    ind, = np.where(np.abs(lon1 - lon2) < 180.0)
    avg[ind] = (lon1[ind] + lon2[ind])/2.0
    return avg

def lon2ew(lon):
    '''Converts 0-360 lons to -180->180 lons. Input: lon value in degrees.'''
    if (lon > 180.0): lon = lon - 360.0
    return lon

def lonslices(lonarr,latarr, tol=90.0):
    '''Creates multiple lon arrays from a single lon array. Returns
    the multiple lon arrays in a list. Creates a new list when the 0/360
    boundary is crossed.
    Will correspondingly break up the latarr to keep pace. 
    Inputs: lonarr in degrees, latarr in degrees (matches lonarr)
    tol: Sep of longitude that defines a new slice. Default 90 degrees.
    Output: dictionary
    {'lon':[lonout1, lonout2, ...]; 'lat':[latout1, latout2,...]}
    output lon arrays that don't have
    a large break due to the crossing.
    '''
    lon = lonToZero360(lonarr)
    londiff = np.abs(lon[1:] - lon[0:-1])
    ind, = np.where(londiff > tol)
    if (len(ind) == 0):
        return {'lon':[lon], 'lat':[latarr]}
    outlon = [lon[0:ind[0]+1]]
    outlat = [latarr[0:ind[0]+1]]
    i=0
    for i in range(1,len(ind)):
        outlon.append(lon[ind[i-1]+1:ind[i]+1])
        outlat.append(latarr[ind[i-1]+1:ind[i]+1])
    outlon.append(lon[ind[i]+1:])
    outlat.append(latarr[ind[i]+1:])
    return {'lon':outlon,'lat':outlat}

def lonindices(lonmin, lonmax, lonarr):
    '''Returns indices into lonarr for which lons are <= lonmax
    and >= lonmin. Deals with the whole boundary crossing thing.
    Inputs:
    lonmin, lonmax: degrees.
    lonarr: degrees. Must be sorted in ascending order, and unique.
    Returns: indices as an array.
    '''
    # Convert everything to 0-360.
    lonminZ = lonToZero360(lonmin)
    lonmaxZ = lonToZero360(lonmax)
    lonarrZ = np.array([map(lonToZero360, lonarr)]).flatten()
    # Handle crossed boundary cases.
    arrcross = False
    if (lonarrZ[-1] < lonarrZ[0]): arrcross = True
    # Handle all the cases. 
    if ((not arrcross) and (lonmaxZ > lonminZ)):
        ind, = np.where((lonarrZ >= lonminZ) & (lonarrZ <= lonmaxZ))
    elif ((not arrcross) and (lonmaxZ < lonminZ)):
        ind1, = np.where((lonarrZ >= lonminZ) & (lonarrZ <= 360.0))
        ind2, = np.where((lonarrZ >= 0.0) & (lonarrZ <= lonmaxZ))
        ind = np.append(ind1.flatten(),ind2.flatten())
    elif ((arrcross) and (lonmaxZ > lonarrZ[-1]) and (lonminZ < lonarrZ[-1])):
        ind1, = np.where((lonarrZ >= lonminZ) & (lonarrZ <= lonarrZ[-1]))
        ind2, = np.where((lonarrZ >= lonarrZ[0]) & (lonarrZ <= lonmaxZ))
        ind = np.append(ind1.flatten(),ind2.flatten())
    elif ((arrcross) and (lonmaxZ < lonminZ)):
        ind1, = np.where((lonarrZ >= lonminZ) & (lonarrZ <= 360.0))
        ind2, = np.where((lonarrZ >= 0.0) & (lonarrZ <= lonmaxZ))
        ind = np.append(ind1.flatten(),ind2.flatten())          
    elif ((arrcross) and (lonmaxZ > lonminZ)):
        ind, = np.where((lonarrZ >= lonminZ) & (lonarrZ <= lonmaxZ))
    return ind

def latLonList(latrange, lonrange, latarr, lonarr):
    '''Creates a list of lat/lon pairs that exist within the specified
    ranges based on the passed-in lat/lon arrays. Inclusive.
    Returns indices into the arrays as a list. 
    '''

    # Make sure I fulfill conditions of lonindices.
    if (not arrOrderCheck(lonarr)):
        print("FATAL ERROR: mylib: latlonlist: input array is not meeting requirements", file=sys.stderr)
        sys.exit(1)
    lonind = lonindices(lonrange[0], lonrange[1], lonarr)
    latind = np.where((latarr <= latrange[1]) & \
                      (latarr >= latrange[0]))
    return [latind, lonind]
        
def dump(fname, headers=None, arrays=None):
    '''
    mylib.dump: debugging function that dumps data from arrays
    into a file. Appends to that file. 
    Inputs:
    fname: Name of file. Opens file for appending.
    headers: array of text strings for column headers.
    arrays: list of numpy arrays to write. All are of same size.
    Prints a time stamp also.
    '''
    ncol = len(headers)
    fout = open(fname,'a')
    now = str(datetime.now())
    fout.write('#'+now+'\n')
    fout.write('#')
    map(lambda w: fout.write(w+' '), headers) # Write headers
    fout.write('\n')
    for i in range(len(arrays[0])):
        st = ''
        for j in range(len(arrays)):
            st = st+'{:17.10E} '.format(arrays[j][i])
        fout.write(st+'\n')
    fout.close()
            
def rd_txt(fname):
    '''Wrapper for the standard approach to reading numeric columns.
    Returns the fields variable. Each variable is accessed via:
    var = fields[:,col_no]. col_no starts at zero. 
    Input: fname: file name string.
    FUTURE: add a keyword that returns text fields. 
    '''
    fields = np.array([map(floatOrMiss, line.split()) for line in dataLines(fname)])
    return fields

def interpretflags(numflags, txt):
    '''Simple utility that interprets the string 1,0,1,1 etc. as a set of
    flags. Useful for setting options on the command line.
    String can be any length. Not all options need be specified (can
    be shorter) which means ignored flags are False.
    Inputs:
    numflags: number of returned flags.
    txt: the string 1,0,0 etc. to set the flags.
    Returns:
    Numpy array [True, False, True, True] etc. that is numflags long. 
    '''
    x = np.array(map(int, txt.split(',')), dtype='bool')
    # Pad out to numflags
    if (numflags < np.size(x)): numflags = np.size(x)
    return np.append(x,np.zeros(((numflags-np.size(x)),),dtype=bool))

def lonToZero360(lon):
    '''Takes a lon value or array and converts it to 0->360 range if needed.
    Input: lon (degrees)
    Output: lon (degrees) as an array (possibly of length 1).'''
    if (np.isscalar(lon)):
        if (lon < 0.0): lon = lon+360.0
        return lon
    # Input is an array
    ind = np.where(lon < 0.0)
    lon[ind] = lon[ind] + 360.0
    return lon    

def geo2magdipole(geolat, geolon):
    '''Takes in geographic lat lon in degrees and returns
    a list with lat/lon in degrees magnetic dipole, according to
    the 1965 formula listed by Russell at
    http://www-ssc.igpp.ucla.edu/personnel/russell/papers/gct1.html/
    NOTE: One of the array elements is missing a minus sign. Could not find
    original publication.
    Validate page: http://www.spenvis.oma.be/help/background/coortran/coortran.html.
    Also: http://wdc.kugi.kyoto-u.ac.jp/igrf/gggm/
    This does NOT match the CGM calculation of ModelWeb.
    '''
    # Create the transformation matrix.
    a = np.zeros((3,3))
    a[0,:] = np.array([0.33907, -0.91964, -0.19826])
    a[1,:] = np.array([0.93826, 0.34594, 0.0])
    a[2,:] = np.array([0.06859,   -0.18602,   0.98015])
    #a = np.transpose(a)
    # Form the unit vector in the geo direction.
    ugeo = tomolib.xyz(1.0, geolat*np.pi/180.0, geolon*np.pi/180.0)
    # Transform to geomag
    umag = np.dot(a, ugeo)
    # Transform back to lat/lon in degrees
    return map(lambda x: x*180.0/np.pi, tomolib.latlonPt(umag))

def sublos(azarr, elarr, stalat, stalon, rearth=6370.0e3, h=350.0e3):
    '''sublos routine calculates subionospheric point. Derived from
    IDL sublos routine.
    Inputs:
    azarr, elarr: Input arrays of az/el in radians.
    stalat, stalon: station latitude and longitude in radians.
    Keywords:
    rearth: radius of earth in km
    h: height of ionosphere in km
    Output:
    [sublatarr, sublonarr] radians
    '''
    sublat = np.zeros(len(azarr)); sindlon = np.zeros(len(azarr))
    sublon = np.zeros(len(azarr)); cosdlon = np.zeros(len(azarr))
    dlon = np.zeros(len(azarr))
    b = np.arccos(np.cos(elarr)/(1.0 + h/rearth)) - elarr
    sinlat = np.sin(stalat)*np.cos(b) + np.cos(stalat)*np.sin(b)*np.cos(azarr)
    sublat = np.arcsin(sinlat)

    ind, = np.where((sublat >= np.pi/2.0) | (sublat <= -np.pi/2.0))
    if (len(ind) > 0):
        sublon[ind] = stalon
    ind, = np.where((sublat <= np.pi/2.0) & (sublat >= -np.pi/2.0))
    if (len(ind) > 0):
        sindlon[ind] = np.sin(b[ind])*np.sin(azarr[ind])/np.cos(sublat[ind])
        cosdlon[ind] = (np.cos(b[ind]) * np.cos(stalat) - \
                        np.sin(b[ind]) * np.sin(stalat) * \
                        np.cos(azarr[ind])) / np.cos(sublat[ind])
        dlon[ind] = np.arctan2(sindlon[ind],cosdlon[ind])
        sublon[ind] = stalon + dlon[ind]
    return [sublat, sublon]

def form_series(npts, lolim, hilim):
    '''form_series creates a data series of npts length between lolim and hilim,
    inclusive.
    Inputs:
    npts: number of points in series.
    lolim: lowest value in series.
    hilim: highest value in series. (> lolim is assumed).
    Returns:
    The series.
    '''
    delta = (hilim - lolim)/float(npts-1)
    return lolim + delta*np.arange(npts)

def ltindices(lower, upper, series):
    '''This function takes in a lower and upper local time and
    returns indices into series that find data in that range. Must
    account for 24 hour crossing.
    Inputs:
    lower, upper: lower and upper limits of LT/MLT in hours.
    series: the data series of LT/MLT.
    Returns: the index into the series.'''
    if (upper > lower): # No crossing of 24 hours assumed.
        ind, = np.where((series >= lower) & (series < upper))
        return ind
    elif (upper < lower): # Crossing of 24 hours occurs.
        ind, = np.where(((series >= lower) & (series <= 24.0)) | \
                        ((series >= 0.0) & (series < upper)))
        return ind
        
def ensure_dir(file):
    '''ensure_dir checks for a path name for the file. If the path
    does not exist, it will be created leading to the file.
    No return value.
    Inputs:
    file: the file to check the path for.
    '''
    d = os.path.dirname(file)
    if not os.path.isdir(d):
        os.makedirs(d)
    return
        
def hyphendate(yyyymmdd):
    '''Creates a hyphenated form of date YYYYMMDD -> YYYY-MM-DD.
    Returns hyphenated form.
    '''
    return yyyymmdd[0:4]+'-'+yyyymmdd[4:6]+'-'+yyyymmdd[6:8]

def unzipfile(fs):
    '''Un gzips a file if it is zipped and returns the unzipped
    file name.
    Input: fs, the file name:
    Returns: unzipped file name.
    All names fully qualified (not just base name)
    '''
    if (fs[-3:] == '.gz'):
        ex('gunzip -f '+fs)
        return fs[0:-3]
    elif (os.path.isfile(fs+'.gz')):
        ex('gunzip -f '+fs+'.gz')
        return fs
    else:
        return fs

def idzipfile(fs):
    '''This will return the file name of the existing file. What
    is passed in may or may not have the .gz ending. What is returned
    is what actually exists.
    '''
    if (fs[-3:] == '.gz'):
        if (os.path.isfile(fs)):
            return fs
        elif (os.path.isfile(fs[0:-3])):
            return fs[0:-3]
        else:
            return ''
    else:
        if (os.path.isfile(fs)):
            return fs
        elif (os.path.isfile(fs+'.gz')):
            return fs+'.gz'
        else:
            return ''
            

def chooseDateFormat(yyyymmdd, rootdir, filepath):
    '''chooseDateFormat will return the format for a date, either
    yyyymmdd or yyyy-mm-dd, depending on which points to the file
    filepath. If the file is not found, returns a blank string.
    Inputs:
    yyyymmdd- string date in said format.
    rootdir- File is in rootdir/date/filepath
    filepath- path to the file name, using either hyphenated date or not.
    Returns:
    A string with correct date format to file, or blank if file cannot be found.
    '''
    # Try non-hyphenated date.
    pathtofile = rootdir+'/'+yyyymmdd+'/'+filepath
    if (os.path.isfile(pathtofile)):
        return yyyymmdd
    else:
        hdate = hyphendate(yyyymmdd)
        pathtofile = rootdir+'/'+hdate+'/'+filepath
        if (os.path.isfile(pathtofile)):
            return hdate
        else:
            return ''

def dateToIonexName(prefix, yymmdd):
    '''Takes a date string (YYYYMMDD or YYMMDD) and creates an IONEX file name from it, 
       based on prefix. 
       Input: 
           prefix (e.g. 'jpli', or 'jplg', etc.
           yymmdd date str either YYMMDD or YYYYMMDD.
       Returns: Ionex name.
    '''
    t = IT.to_J2000(yymmdd)
    s = IT.from_J2000(t,'YYDOY')
    return prefix+s[2:5]+'0.'+s[0:2]+'i'

def IonexNameToDate(prefix, name):
    '''Takes an IONEX file name and returns the YYYYMMDD format, and
    seconds past J2000. Uses prefix to find doy (scans past prefix).
    Assumes prefix is unique. Also works for .nc file names. 
    '''
    if (name.find(prefix) < 0):
        print("ERROR: IonexNameToDate: could not find file prefix in name. Returning blanks.", file=sys.stderr)
        return []
    i = name.find(prefix)+len(prefix)
    doy = name[i:i+3]
    yr = name[i+3+2:i+3+4]
    t = IT.to_J2000(yr+doy,'YYDOY')
    d = IT.from_J2000(t,'YYYYMMDD')
    return [hyphendate(d), t]

def findNearest(val, arrval):
    '''Finds index of passed in array that most closely matches
    input value.
    Inputs:
    val: value to match
    arrval: the array
    Returns:
    Index into arrval. arrval[index] is closest matching value.
    '''
    return int((np.abs(arrval - val)).argmin())

def ionex2ncFileName(fn):
    '''Converts an ionex file name to a netCDF file name. Just appends
    .nc. Assumes file is not compressed.
    Inputs:
    fn: the ionex name.
    Returns:
    the netCDF file name.
    '''
    return fn+".nc"

class DumpFile(object):
    '''Class for handling dump files. Mainly, to read in time sequential
    order without storing all of the data.'''
    def __init__(self, filename):
        '''Init method for DumpFile class. Pass in the file name to
        read. Will create the generator object that is subsequently used.'''
        self.filename = filename
        self.gen = dataLines(self.filename)
        self.tstart = 1.0e27
        self.tend = -1.0e27

    def rdtimelimits(self, start, stop, cols):
        '''Input time limits in J2000 (start, stop) and cols which
        is an array of 0-starting column numbers. A list of data
        within the specified time window, for those columns, will
        be returned as a list of numpy 1-d arrays.
        Note: sequential calling. Once a certain time is specified,
        the function will not rewind to earlier times. MUST BE CALLED
        IN TIME SEQUENTIAL ORDER.'''
        lineslist = []
        try:
            line = next(self.gen)
        except StopIteration:
            return []
        tcurrent = float(line.split()[0])
        # Reach the begin time
        while (tcurrent < start):
            try:
                line = next(self.gen)
            except StopIteration:
                break
            tcurrent = float(line.split()[0])
        # We are now ready to start reading
        choke = False
        while ((tcurrent >= start) and (tcurrent < stop)):
            lineslist.append(line)
            try:
                line = next(self.gen)
            except StopIteration:
                choke = True
                break
            tcurrent = float(line.split()[0])
        if (not choke): lineslist.append(line)
        fields = np.array([map(floatOrMiss, line.split()) for line in lineslist])
        retval = []
        for i in range(len(cols)):
            retval.append(fields[:,cols[i]])

        return retval

    def rdnlines(self, n, cols):
        '''Input number of lines to be read, and cols which
        is an array of 0-starting column numbers. A list of data,
        for those columns, will
        be returned as a list of numpy 1-d arrays.
        Note: sequential calling. Picks up where it was called last.'''
        lineslist = []
        try:
            line = next(self.gen)
        except StopIteration:
            return []
        # We are now ready to start reading
        nline = 0
        choke = False
        while (nline < n-1):
            lineslist.append(line)
            try:
                line = next(self.gen)
            except StopIteration:
                choke = True
                break
            nline = nline + 1
        if (not choke): lineslist.append(line)
        fields = np.array([map(floatOrMiss, line.split()) for line in lineslist])
        retval = []
        for i in range(len(cols)):
            retval.append(fields[:,cols[i]])
                        
        return retval

def matlabLDiv(a, b):
    '''Implements the matlab \ operator on arrays (maybe works
    on matrices also).
    Solves the "left divide" equation: a * x = b for x
    Inputs:
    a,b arrays
    Returns: a \ b
    Results depend on dimensions of matrices. See documentation for
    matlab's MLDIVIDE operator \ .
    Theory is on the Numpy for Matlab user's page.
    http://wiki.scipy.org/NumPy_for_Matlab_Users
    See also this site, but beware of caveats
    http://stackoverflow.com/questions/1001634/array-division-translating-from-matlab-to-python
    '''
    import scipy.linalg as LA
    # Check if a is square.
    try:
        (r,c) = np.shape(a)
    except ValueError: # In case a is of dim=1, cannot unpack two vals.
        return LA.solve(a,b)
    else:
        if (r == c): # Square
            return LA.solve(a,b)
        else:
            return LA.lstsq(a,b)
    
def matlabRDiv(a, b):
    '''Carries out matlab / operator on arrays (maybe works
    on matrices also).
    Solves the "right divide" equation: x*a = b for x.
    Inputs: a,b arrays
    Returns: a / b
    Results depend on dimensions of matrices. See documentation for
    matlab's MRDIVIDE operator / .
    Theory is on the Numpy for Matlab user's page.
    http://wiki.scipy.org/NumPy_for_Matlab_Users
    '''
    return np.transpose(matlabLDiv(b.T, a.T))

def getTmpDir():
    '''Provides a directory for temporary files.
    No inputs. Returns directory name as a string without trailing slash.
    '''
    td = os.environ['TMPDIR']
    if (td[-1] == '/'):
        return os.environ['TMPDIR'][0:-1]
    else:
        return os.environ['TMPDIR']

def createShadowFile(fs):
    '''Copies file to a temp area for further file manipulation and
    read. Typically might us TMPDIR.
    Input: file name to be copied.
    Returns: location of copied file.
    '''
    td = getTmpDir()
    target = td+'/'+os.path.basename(fs)
    try:
        shutil.copyfile(fs, target) # Copies file
        return target
    except:
        return ''
    
def meanUT(utarr):
    '''Computes the mean of a UT array.
    Takes account of the fact that UT can bridge the 24 hour discontinuity,
    assuming that the UT range in the array is < 12 hours. 
    Input: array of UT hours.
    Returns: mean hour of array, meaningfully interpreted.
    NOTE: utarr is modified. 
    '''
    inarr = utarr.copy()
    tmin = np.min(inarr)
    tmax = np.max(inarr)
    if (np.abs(tmax - tmin) > 12.0):
        ind, = np.where(inarr < 12.0)
        inarr[ind] = inarr[ind] + 24.0
    m = np.mean(inarr)
    if (m > 24.0): m = m - 24.0
    return m

def arrOrderCheck(arr):
    '''Checks that an array is in ascending order and unique. If not, returns
    False.
    '''
    [tmparr, tmpind] = np.unique(arr, return_index=True)
    if (np.array_equal(tmpind, np.arange(len(arr)))):
        return True
    else:
        return False
    
