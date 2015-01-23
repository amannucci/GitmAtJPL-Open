#
# Copyright 2014, by the California Institute of Technology.
# ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged.
# Any commercial use must be negotiated with the
# Office of Technology Transfer at the California Institute of Technology.
#
# AJ Mannucci, JPL.

from __future__ import print_function
import copy
import glob
import os

import gitm
import mylib
import ionotime as IT
import numpy as np
import sys
#import scipy.io.netcdf as CDF
from netCDF4 import Dataset

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
        return [avg/len(indlon),-1]

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
        return [avg/len(indlon), avglt/len(indlonlt), -1]

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
            avglt = 0.0
            for i in range(nalt):
                avglt = avglt + np.mean(data[indlon, indlat, i])
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
            indalt = -1
    else:
        sys.stderr.write("ERROR: gitmLib:boundaryAvg: Incorrect mode passed in ("+mode+"). Allowed values are Lat or Lon. \n")
        sys.exit(1)
    return [avglon, avglat, avglt, indalt]

def locValue(mObject, time, var, alt, lat, lon):
    '''Returns time series of values at lat/lon. Altitude or integrated. 
    Inputs:
    mObject: an instance of a ModelHandler object (gitmLib)
    time: J2000 time that variable is desired.
    var: the variable name
    alt: the altitude or interest. Function will compute the index for
      a nearby altitude. Use a negative value to obtain height-integrated.
      Meters.
    lat: latitude of interest. Will find nearest. Radians.
    lon: longitude of interest. Will find nearest. Radians.
    Returns:
    [value, alt index used] alt index used if height integrated
    is -1.
    '''
    # Get the array. [lon,lat,alt]
    data = mObject.getVariable(var, time)
    # If a specific height, get the altitude index.
    # We assume altitude is independent of lon/lat.
    latval = collapseDim(mObject.getVariable('Latitude',time), 1)
    lonval = collapseDim(mObject.getVariable('Longitude',time), 0)
    indlat = mylib.findNearest(lat, latval)
    indlon = mylib.findNearest(lon, lonval)
    if (alt > 0):
        altval = collapseDim(mObject.getVariable('Altitude',time),2)
        indalt = mylib.findNearest(alt, altval)
        return [data[indlon, indlat, indalt],indalt]
    else: # Height integrated.
        return [np.sum(data[indlon, indlat, :]),-1.0]

class TgcmBin(object):
    '''A Tiegcm class'''
    def __init__(self, filename):
        try:
            self.root = Dataset(filename,'r')
        except:
            sys.stderr.write("FATAL ERROR: gitmLib: TgcmBin. File open error: {:s}\n",format(filename))
            sys.exit(1)
            
    def __new__(self, filename):
        '''This new method is over-written so I can return a value.
        That follows the gitmBin convention.
        '''
        try:
            self.root = Dataset(filename, 'r')
        except RuntimeError:
            sys.stderr.write("FATAL ERROR: gitmLib: TgcmBin. File open error: {:s}\n".format(filename))
            sys.exit(1)
        tgcm_newdict = {}
        for i in self.root.variables.keys():
#            if (len(np.shape(self.root.variables[i])) == 0):
#                continue
            #print(np.shape(cls.f.variables[i]))
            tgcm_newdict[i] = self.root.variables[i][:].copy()
        self.root.close()
        return tgcm_newdict
        
def modelType(filename):
    '''Returns 'gitm' or 'tgcm' string depending on model type,
    as deduced from the file name.
    '''
    if (filename[filename.rfind('.'):] == '.bin'):
        return 'gitm'
    elif (filename[filename.rfind('.'):] == '.nc'):
        return 'tgcm'
    else:
        sys.stdout.write("FATAL ERROR: gitmLib.modelType, unknown file ending, cannot determine model type.\n")
        sys.exit(1)

def ModelBin(filename, *args, **kwargs):
    '''Generic interface to GitmBin and TgcmBin objects for reading
    their respective file types.
    Inputs:
    filename. Uses this to determine the model type.
    Optional args: only passed to GitmBin object.
    '''
    ftype = modelType(filename)
    if (ftype == 'gitm'):
        return gitm.GitmBin(filename, *args, **kwargs)
    elif (ftype == 'tgcm'):
        return TgcmBin(filename)
    else:
        sys.stderr.write("FATAL ERROR: gitmLib.ModelBin: Unsupported file type.\n")
        sys.exit(1)

class ModelHandler(object):
    '''This is the generic model handler that extracts time series for the
    GITM and TIEGCM models, and performs variable name translation.
    GITM variable names (gitm.py) are viewed as primary. If the model is
    TIEGCM, translation must occur, unless that variable name is
    unique to TIEGCM.
    '''
    # The idea behind the translation layer is as follows:
    # Calling ModelHandler returns an instance of the class (standard).
    # Methods of the class do the work. This includes a method that
    # returns the standard GitmBin class instance in case that wants to be
    # used.
    # filespec can include wild-carding. filename is a single name.
    def __init__(self, filespec, *args, **kwargs):
        self.modelFileNames = expandFileSpec(filespec)
        self.modelType = modelType(self.modelFileNames[0])
        self.timeDict = self.fillModelTimeDictionary(filespec)
        self.varNotTimed = fillNotTimedList(self.modelType)
        self.modelTimes = np.sort(self.timeDict.keys())
        self.modelNameDict = modelTransTable(self.modelType)
        self.modelUnitConv = modelUnitConverter(self.modelType)
        self.gitmArgs = args
        self.gitmKwargs = kwargs
        self.lastModelFileName = "ThisIsNeverAModelFileName.XXX"
        self.modelObject = 0 # This will become a GitmBin object or its tiegcm equivalent
        return

    def fillModelTimeDictionary(self, filespec):
        '''Returns the time/filename dictionary, based on
        a filename input, which includes wildcards.
        Expands filename into a list of files, then queries each
        file for the times within. Associates
        certain times with specific file names.
        '''
        if (self.modelType == 'tgcm'):
            return getTgcmTime(filespec)
        elif (self.modelType == 'gitm'):
            return getGitmTime(filespec)
        else:
            sys.stderr.write("FATAL ERROR: ModelHandler.fillModelTimeDictionary: unrecognized model type: {:s}\n".format(self.modelType))
            sys.exit(1)

    def getVariable(self, varname, time):
        '''Retrieves a specific model variable from a model file.
        Returns it as an appropriate numpy array, e.g. in 3 dimensions.
        Inputs: varname - uses GITM conventions (these may evolve, but have
        to be meaningful to a GITM file. See utility pullgitm.py for
        retrieving variable names from a file.
        time - is in seconds past J2000. If it does not match the times in
        the output file, nearest time will be used (no interpolation).
        Returns: the variable array as a numpy array.
        '''
        if (self.modelType == 'tgcm'):
            filename = time2FileName(time, self.timeDict)
            cacheHit = checkModelCache(filename, self.lastModelFileName)
            if (not cacheHit): self.modelObject = 0 # Save memory
            t = getTgcmObject(filename, time, \
                              self.lastModelFileName, \
                              self.modelObject)
            self.modelObject = t
            # Since a tgcm object can contain multiple times, we have to select the
            # appropriate array for a specific time.
            var = getTgcmVar(varname, self.modelNameDict, time, t, self.varNotTimed, self.modelUnitConv)
            # Update cache stuff.
            self.lastModelFileName = updateCache(filename)
            return var
        elif (self.modelType == 'gitm'):
            filename = time2FileName(time, self.timeDict)
            cacheHit = checkModelCache(filename, self.lastModelFileName)
            if (not cacheHit): self.modelObject = 0 # Recover memory
            t = getGitmObject(varname, time, self.timeDict, \
                              self.gitmArgs, self.gitmKwargs, \
                              self.lastModelFileName, \
                              self.modelObject)
            self.modelObject = t
            filename = time2FileName(time, self.timeDict)
            if (not (varname in t.keys())):
                sys.stderr.write("ERROR: gitmLib: ModelHandler: getVariable: Requested variable >>>"+varname+"<<< not in file "+filename+". Aborting.\n")
                sys.exit(1)
            # Update cache stuff.
            self.lastModelFileName = updateCache(filename)
            return np.array(t[varname])

######End of ModelHandler

def getTgcmTime(filespec):
    '''Pull out times from a series of TGCM output file specified by
    filespec (wildcard) and fill a dictionary associating times with
    file names. Times are in seconds past J2000 (floats).
    Input: wild card spec for the files.
    Output: Dictionary with keys as the times, and values as file
    file names hosting that time
    '''
    # Expand the filespec into a list.
    flist = glob.glob(filespec)
    if (len(flist) < 1):
        sys.stderr.write("FATAL ERROR: gitmLib.py:getTgcmTime: Files not found: {:s}\n".format(filespec))
        sys.exit(1)
    flist.sort()
    d = {}
    for fs in flist:
        # netCDF4 seems to hang if file not found, so...
        if (os.path.exists(fs)):
            root = Dataset(fs, mode='r')
        else:
            sys.stdout.write("FATAL ERROR: gitmLib.py: getTgcmTime: file not found: {:s}".format(fs))
            sys.exit(1)
        # Recover all times in current file.    
        mtime = root.variables['mtime'][:]
        myear = root.variables['year'][:]
        # Do the date computations.
        for (i,t) in enumerate(mtime):
            year = str(myear[i])
            doy = '{:03d}'.format(t[0])
            date2000 = IT.to_J2000(year+doy,'YYYYDOY')
            start2000 = date2000+60.0*t[2]+3600.0*t[1]
            d[start2000] = fs
    root.close()
    return d    

def getGitmTime(filespec):
    '''Pull out times from a series of GITM output files specified by
    filespec (wildcard) and fill a dictionary associating times with
    file names. Times are in seconds past J2000 (floats).
    Input: wild card spec for the files.
    Output: Dictionary with keys as the times, and values as
    file names hosting that time
    '''
    # Expand the filespec into a list.
    flist = glob.glob(filespec)
    if (len(flist) < 1):
        sys.stderr.write("FATAL ERROR: gitmLib.py:getGitmTime: Files not found: {:s}\n".format(filespec))
        sys.exit(1)
    flist.sort()
    d = {}
    for fs in flist:
        # Recover all times in current file. One time per file, as per
        # GITM convention.
        binf = gitm.GitmBin(fs)
        dt = binf['time'] # python date time object
        secJ2000 = IT.to_J2000((dt.year, dt.month,\
                        dt.day, dt.hour, \
                        dt.minute, \
                        dt.second+\
                        1.0e-6*dt.microsecond))
        d[secJ2000] = fs
    return d
 
def getTgcmObject(filename, time, lastModelFileName,\
                  modelObject):
    '''Retrieves a TIEGCM variable at a specific time. For names of variables,
    use pullgitm.py on a GITM file. GITM variables names are
    used if there are common variables shared by GITM and TIEGCM. If a variable
    only exists in a TIEGCM file, then the TIEGCM name is used. See translation table
    in modelTransTable.
    Inputs:
    varname: text name of the variable.
    time: time variable wanted, seconds past J2000. Will pull in closest
      time.
    tdict: a dictionary of time:filename needed to retrieve the data.
    nameDict: the name translation dictionary.
    lastModelFileName: for checking cache. If unchanged, don't
       open a new file.
    modelObject: this is what is returned also, but if there is a cache hit, it just
       returns what is passed in. This is the model dictionary containing the variables
       as arrays.
    Associated with ModelHandler class.
    Returns: the variable as a numpy array.
    '''
    # Check for cache hit or miss.
    cacheHit = checkModelCache(filename, lastModelFileName)
    # Retrieve the data.
    if (not cacheHit):
        return TgcmBin(filename)
    else:
        return modelObject

def getTgcmVar(varname, modelNameDict, time, tgcmObject, varNotTimed, modelUnitConv):
    '''Returns a TIEGCM data array for a specific time, based on the passed in object
    tgcmObject (which is really a dictionary). May need to recover additional variables
    to translate from units such as "mass mixing ratio". 
    Inputs: varname, variable name, matching GITM convention if possible.
    time: J2000
    tgcmObject: the dictionary containing the data for a specific time.
    varNotTimed: whether the variable has a time dimension or not.
    modelUnitConv: dictionary for unit conversions
    Returns: numpy array of that variable for a specific time.
    '''
    try:
        tgcmname = modelNameDict[varname]
    except KeyError:
        tgcmname = varname
    if (varname in varNotTimed):
        return swapTgcmArrays(tgcmObject[tgcmname])
    else:
        # Figure out time index. FUTURE: do this once at ModelHandler initialization.
        timename = modelNameDict['time']
        modeltimes = tgcmObject[timename]
        modelyears = tgcmObject['year']
        timeJ2000 = np.zeros(len(modeltimes))
        for (i,t) in enumerate(modeltimes):        
            year = str(modelyears[i])
            doy = '{:03d}'.format(t[0])
            date2000 = IT.to_J2000(year+doy,'YYYYDOY')
            timeJ2000[i] = date2000+60.0*t[2]+3600.0*t[1]
        indtime = mylib.findNearest(time, timeJ2000)
    # Get unit conversion factor. In the dictionary modelUnitConv, > 0 means
    # it's a simple multiplicative factor.
    # < 0 means if may require conversions such as from mass mixing ratio.
    # This is performed within the subroutine.
    unitfac = tgcmUnitConv(modelUnitConv, varname, tgcmObject, indtime)
    # unitfac might be a scalar or full array.
    return swapTgcmArrays(unitfac*tgcmObject[tgcmname][indtime,:])
           
def getGitmObject(varname, time, tdict, args_in, kwargs_in, \
                  lastModelFileName, modelObject):
    '''Retrieves a GITM variable at a specific time. For names of variables,
    use pullgitm.py on the file of interest.
    Inputs:
    varname: text name of the variable.
    time: time variable wanted, seconds past J2000. Will pull in closest
      time.
    tdict: a dictionary of time:filename needed to retrieve the data.
    args_in, kwargs_in: the arguments and keywords that would be used in the call
      to GitmBin, initializing the class (per U Mich's gitm.py library).
    Associated with ModelHandler class.
    Returns: the variable as a numpy array.
    '''
    # Get the file name containing the time of interest.
    filename = time2FileName(time, tdict)
    # Check for cache hit or miss.
    cacheHit = checkModelCache(filename, lastModelFileName)
    # Retrieve the data.
    if (not cacheHit):
        if (kwargs_in == {}):
            modelObject = gitm.GitmBin(filename, args_in)
        else:
            modelObject = gitm.GitmBin(filename, args_in, kwargs_in)

    return modelObject
    
def modelTransTable(modeltype):
    '''Returns a dictionary that translates from GITM variable names
    to names appropriate to modeltype. Key: GITM name, value: TIEGCM.
    Input: modeltype, e.g. 'gitm' or 'tgcm'.
    Returns: the dictionary gitmname:modeltype_name.
    For GITM, use pullgitm.py (and gitm.py library), and for TIEGCM
    use ncdump on the TIEGCM netCDF files or see TIEGCM manual.
    '''
    if (modeltype == 'tgcm'):
        # Notes: Chemical heating (gitm)/total heating (tgcm)?
        return {'N2 Density':'N2D', 'O2 Density':'O2', \
                'Conduction':'', 'Altitude':'ZG', \
               'O Cooling':'', 'Potential':'POTEN', 'Ion Velocity (north)':'VI_EXB', \
               'Latitude':'lat', 'Longitude':'lon', 'Altitude':'ZG', \
                'EUV Heating':'', 'CO2 Cooling':'CO2_COOL', \
                'Ion Velocity (up)':'WI_EXB', 'Electron Density':'NE', \
                'LT':'', 'Neutral Velocity (East)':'UN', 'NO Density':'NO', \
                'Neutral Velocity (north)':'VN', 'O Density':'O1', \
                'Chemical Heating':'','Neutral Velocity (up)':'', \
                'dLon':'', 'Ion Velocity (east)':'UI_EXB', \
                'Auroral Heating':'', 'Longitude':'lon', 'NO Cooling':'NO_COOL',\
                'Neutral Temperature':'TN', 'Joule Heating':'QJOULE', \
                'Photoelectron Heating':'', 'Neutral Density':'DEN', \
                'dLat':'', 'time':'mtime', 'Region 2 Current':'', \
                'Region 1 Current':''}
    else:
        return {} # No dictionary if GITM is the model.

def modelUnitConverter(modelType):
    '''Provides a multiplicative factor to bring TIEGCM variables to their
    GITM equivalents in units. = 1 if there in no GITM equivalent.
    TIEGCM can also be in "mass mixing ratio" units, meaning total density
    (g/cm3) is needed to calculate density in GITM units of 1/m^3. These special
    cases are encoded as negative conversion factors, with the value telling
    the code how to convert (see unitConv function).
    -(small number) means multiply total density by mass mixing ratio, etc.
    The small number is the relevant mass per particle.
    -1000.0 is for energy conversion (erg/g/s versus ??? for GITM, ask Xing)
    A value of None means a GITM variable without TIEGCM equivalent. 
    TIEGCM value * this unit converter value = GITM value.
    Remarks:
    CO2 Cooling: erg/g/s according to TIEGCM docs, don't know what GITM unit is.
      Similarly for NO Cooling, Joule Heating, etc.
    '''
    d2r = np.pi/180.0
    if (modelType == 'tgcm'):
        return {'N2 Density':-14.00674, 'Conduction':None, 'Altitude':0.01, \
               'O Cooling':None, 'Potential':1.0, 'Ion Velocity (north)':0.01, \
               'Latitude':d2r, 'Longitude':d2r, 'Altitude':0.01, \
                'EUV Heating':None, 'CO2 Cooling':1.0, \
                'Ion Velocity (up)':0.01, 'Electron Density':1.0e6, \
                'LT':None, 'Neutral Velocity (East)':0.01, 'NO Density':-1.0, \
                'Neutral Velocity (north)':0.01, 'O Density':-15.9994, \
                'O2 Density':-31.988, \
                'Chemical Heating':None,'Neutral Velocity (up)':0.01, \
                'dLon':None, 'Ion Velocity (east)':0.0, \
                'Auroral Heating':None, 'Longitude':d2r, 'NO Cooling':1.0,\
                'Neutral Temperature':1.0, 'Joule Heating':-1000.0,
                'Photoelectron Heating':None, 'Neutral Density':1000.0, \
                'dLat':None, 'time':1.0, 'Region 2 Current':1.0, \
                'Region 1 Current':1.0}
    else:
        return {} # No dictionary if GITM is the model.

def tgcmUnitConv(converterDict, varname, tgcmObject, indtime):
    '''Will apply the unit conversion to create GITM like variables.
    Requires special cases, because sometimes unit conversion involves more
    than a simple multiplication, e.g. TIEGCM units may be mass mixing ratios,
    whereas GITM may use absolute mass densities. Returns conversion factor
    appropriate to the time of the data.
    Inputs:
    converterDict: the dictionary containing conversion factors and codes.
    varname: standard variable name
    tgcmObject: The full tgcm object (dictionary) in case fields need to
      be accessed to effect the conversion.
    indtime: time index into tgcmObject.
    Returns:
    Conversion scalar.
    '''
    if (converterDict[varname] > 0.0):
        # Simple scalar
        return converterDict[varname]
    elif ((converterDict[varname] < 0.0) and \
          (converterDict[varname] > -1000.0)):
        # Need to multiply tgcm quantity by total density.
        # Gets total mass density. Factor of 1.0e6 converts from
        # g/cm3 to 1/m3. To get number density, need molecular/atomic
        # weight.
        amu = 1.6726e-24 # grams, atomic mass unit
        molmass = -converterDict[varname] # Depends on species.
        den = tgcmObject['DEN'][indtime,:]
        return (den/(molmass*amu))*1.0e6
    elif (converterDict[varname] == -1000.0):
        # This is for energy conversion. Just use 1.0 for now until
        # this is better understood.
        return 1.0
        
def fillNotTimedList(modelType):
    '''Fills a list that returns variable names that are not time dependent.
    Names are GITM names, except when
    there is no GITM equivalent for a TGCM name. If a name is not listed here, it
    is assumed to have time dependence.
    Returns an empty list for GITM. Should not be needed for models that have one
    time per file.
    '''
    if (modelType == 'tgcm'):
        return ['Latitude','Longitude','time',\
                'lev','ilev','mlon','mlat','mlev','imlev',\
                'timestep','mag','p0','p0_model','grav','LBC']
    else:
        return [] # Not used for GITM model.

def time2FileName(time, tdict):
    '''Converts a time to a file name containing that time, according to
    ModelHandler conventions. Inputs are the time (J2000) and the ModelHandler
    time/filename dictionary. Returns the file name.
    '''
    times = np.array(tdict.keys())
    files = np.array(tdict.values())
    return files[mylib.findNearest(time, times)]

def checkModelCache(filename, lastModelFileName):
    '''Logic to check if there is a cache hit for the model objects.
    Inputs: time J2000, filename (containing that time), previous model time,
    and previous model file name (text).
    Note that some model objects have multiple times per object. GITM has one
    time per object/file.
    Return True if cache is hit, False otherwise.
    '''
    if (filename != lastModelFileName):
        return False
    else:
        return True

def updateCache(filename):
    '''
    '''
    return filename

def expandFileSpec(filespec):
    '''Creates a file list from a potentially wild-carded filespec.
    Returns the list.
    '''
    return glob.glob(filespec)

def swapTgcmArrays(arr):
    '''Swaps arrays from tgcm ordering to GITM ordering. (Absence of time
    dimension is assumed). If 2D, lat/lon (tiegcm) -> lon/lat. If 3D,
    alt/lat/lon (tiegcm) should be lon/lat/alt. If 1D, it's possible that
    it's not 1D in gitm.
    Input: arr to be processed. (tiegcm style)
    Returns: swapped arr.'''
    # Swap axes to meet GITM conventions. Leave 1-d and scalars alone.
    # 2-d arrays, if lat/lon should become lon/lat.
    # 3-d arrays, if alt/lat/lon should become lon/lat/alt.
    # For now, assume 2d is lat/lon and 3d is alt/lat/lon
    sh = np.shape(arr)
    if (len(sh) <= 1):
        return arr
    if (len(sh) == 2):
        return np.swapaxes(arr, 0, 1)
    elif (len(sh) == 3):
        return np.swapaxes(arr, 0, 2)
    else:
        sys.stderr.write("FATAL ERROR: gitmLib: swapTgcmArrays: cannot interpret array size. Terminating.\n")
        sys.exit(1)

def collapseDim(arr,dim):
    '''Handles the fact that some models have lower dimensional arrays
    than others. Specifically, lat/lon in GITM is a full 3D array, whereas
    in tgcm, it's one d.
    Inputs:
    arr: the 1d or 3d array.
    dim: the dimension we want.
    Returns: 1d array.
    '''
    if (len(np.shape(arr)) == 3):
        if (dim == 0):
            return arr[:,0,0]
        elif (dim == 1):
            return arr[0,:,0]
        elif (dim == 2):
            return arr[0,0,:]
    elif (len(np.shape(arr)) == 2):
        if (dim == 0):
            return arr[:,0]
        elif (dim == 1):
            return arr[0,:]
    elif (len(np.shape(arr)) == 1):
        return arr
        
    
if (__name__ == "__main__"):
    # For testing new model handling paradigm.
    gitmname = 'GitmRunDir/3DUSR_t110425_000000.bin'
    tgcmname = 'tgcm.nc'
#    d = TgcmBin(tgcmname)
    gitmO = ModelHandler(gitmname)
    j2000 = IT.to_J2000((2011,4,25,1,0,0.0))
    ed = gitmO.getVariable('Electron Density', j2000)
    print(ed)
    tgcmO = ModelHandler(tgcmname)
    j2000 = IT.to_J2000((2011,4,15,6,0,0.0))
    ed2 = tgcmO.getVariable('Electron Density', j2000)
    print(ed2)
# More testing needed:
# Caching.
# Non-time based variables, such as 'Latitude'
# Once this is worked out, we have the basic machinery, but there
# remain differences, such as the order of the arrays, etc.
    # Now get neutral density
    nd = tgcmO.getVariable('Neutral Density', j2000)
    # Now try another time.
    nd2 = tgcmO.getVariable('Neutral Density', j2000+86400.0)
    # Now get altitude coords.
    zg = tgcmO.getVariable('Altitude', j2000)
    ilev = tgcmO.getVariable('ilev', j2000)
    refPres = tgcmO.getVariable('p0', j2000)
    print(refPres)
    refPres2 = tgcmO.getVariable('p0', j2000+86400.0)
    print(refPres2)
