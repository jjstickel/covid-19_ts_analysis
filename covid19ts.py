"""
Module to collect reading and processing of Johns Hopkins CSSE COVID-19
timeseries data 

https://github.com/CSSEGISandData/COVID-19)

"""

## TODO:
# - check and correct for double counting of cities? e.g. New York City

import numpy as np
import pandas as pd
#import regularsmooth as ds # my local copy
import scikits.datasmooth as ds # pip installed

popfile = ("API_SP.POP.TOTL_DS2_en_csv_v2_887275/"
           "API_SP.POP.TOTL_DS2_en_csv_v2_887275_2018.csv")

##### scaling factor for cases, i.e., "x per mult" ######
mult=1e6
critlow = 10*1e-6 # for time zero, using confirmed
#critlow = 0.1*1e-6 # for time zero, using deaths


def covid19_global(countries, websource=True, JHCSSEpath=None, file_pop=popfile, mult=mult,
                   critlow=critlow):
    """
    Read in Johns Hopkins CSSE COVID-19 timeseries data for the countries
    specified and perform these operations:
    - normalize cases to be per capita
    - smooth the cases data
    - set time-zero for each location and shift elapsed time in days
    - determine rates (i.e., the derivative) for cases

    Optional arguments:
    websource:      whether to get time-series data from the web source or local files
    JHCSSEpath:     path to the local files if websource=False
    file_pop:       path to the population file (default is provided)
    mult:           scaling factor for the data, default is 1e6
    critlow:        cutoff criteria for time zero, default is 10e-6
    """
    
    # read in Johns Hopkins' data tables
    pathname = setup_path(websource, JHCSSEpath)
    # read in global data
    data_confirmed = pd.read_csv(pathname + "time_series_covid19_confirmed_global.csv")
    data_deaths = pd.read_csv(pathname + "time_series_covid19_deaths_global.csv")
    data_recovered = pd.read_csv(pathname + "time_series_covid19_recovered_global.csv")
    # presume these files all have the same "fixed" structure -- could put in a
    # check at some point
    header = data_confirmed.columns.values
    dates = pd.to_datetime(header[4:])

    # read in population data (source: https://data.worldbank.org/indicator/sp.pop.totl)
    data_pop = pd.read_csv(file_pop, header=4)

    corona = dict()
    corona["locs"] = countries # I'm not sure if I like this approach... JJS 4/15/20
    for country in countries:
        ctryd = dict()
        corona[country] = ctryd
        ctryd['name'] = country
        ctryd['population'] = read_pop(data_pop, country)
        ctryd['cnf'] = read_cases(data_confirmed, country)
        ctryd['dth'] = read_cases(data_deaths, country)
        ctryd['rec'] = read_cases(data_recovered, country)

    corona["mult"] = mult
    corona["critlow"] = critlow
    corona["dates"] = dates

    perform_operations(corona)
    
    return corona


def covid19_US(locations, websource=True, JHCSSEpath=None, mult=mult, critlow=critlow):
    """
    Read in Johns Hopkins CSSE COVID-19 timeseries data for the US locations
    specified and perform these operations:
    - normalize cases to be per capita
    - smooth the cases data
    - set time-zero for each location and shift elapsed time in days
    - determine rates (i.e., the derivative) for cases
    
    `locations` argument should be a list of location names in the form of
    `state` OR `county, state`.
    
    Optional arguments:
    websource:      whether to get time-series data from the web source or local files
    JHCSSEpath:     path to the local files if websource=False
    mult:           scaling factor for the data, default is 1e6
    critlow:        cutoff criteria for time zero, default is 10e-6

    """
    
    # read in Johns Hopkins' data tables
    pathname = setup_path(websource, JHCSSEpath)
    # read in US data
    data_confirmed = pd.read_csv(pathname + "time_series_covid19_confirmed_US.csv")
    data_deaths = pd.read_csv(pathname + "time_series_covid19_deaths_US.csv")
    # presume these files all have a known "fixed" structure -- could put in a
    # check at some point
    header = data_confirmed.columns.values
    dates = pd.to_datetime(header[11:]) # note:  data_deaths has an extra population column

    corona = dict()
    corona["locs"] = locations # I'm not sure if I like this approach... JJS 4/15/20
    for loc in locations:
        locd = dict()
        corona[loc] = locd
        locd['name'] = loc
        # handle data handling differently by whether location is a state or
        # county--if a state, then need to read in all the state entries and
        # sum
        if "," in loc: # it's a county
            idx = loc.find(",")
            county = loc[:idx]
            state = loc[idx+2:]
            loc_bool = np.logical_and(data_deaths["Admin2"] == county,
                                      data_deaths["Province_State"] == state)
            if loc_bool.sum() != 1:
                raise Exception("%s is not in the data files, or occurs more than once"
                                % loc)
            #locd["population"] = data_deaths[loc_bool]["Population"].values[0]
            #locd["cnf"] = data_confirmed[loc_bool].iloc[0,11:].values
            #locd["dth"] = data_deaths[loc_bool].iloc[0,12:].values
        else: # it's a state
            loc_bool = data_deaths["Province_State"] == loc
            if loc_bool.sum() < 1:
                raise Exception("%s is not in the data files" % loc)
        # these lines should work even for county location with single line...
        locd["population"] = data_deaths[loc_bool]["Population"].sum()
        locd["cnf"] = data_confirmed[loc_bool].iloc[:, 11:].values.sum(axis=0)
        locd["dth"] = data_deaths[loc_bool].iloc[:, 12:].values.sum(axis=0)

    corona["mult"] = mult
    corona["critlow"] = critlow
    corona["dates"] = dates

    perform_operations(corona)
    
    return corona


def setup_path(websource, JHCSSEpath):
    if websource:
        pathname = "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/"
    elif JHCSSEpath is None:
        raise ValueError("If websource=False, you need to provide a path for the JH CSSE data")
    else:
        pathname = JHCSSEpath
    return pathname

def read_cases(data, country):
    """
    Select country data values. The country data must be in a single row for the
    current implementation.
    """
    # get the row for the country
    c_bool = data["Country/Region"] == country
    c_data = data[c_bool]
    if c_data.shape[0] is 1:
        return np.squeeze(c_data.iloc[:,4:].values)
    else:
        raise Exception("Not implemented:  there is more than one row (or no rows) of %s data." % country)
    
def read_pop(data, country):
    """ Get the population for a country """
    # There may be some name mismatches that need correction in the population
    # file -- please create an issue or pull request when you find them
    row = data[ data["Country Name"]==country ]
    if row.size == 0:
        raise Exception("%s is not in the population data file" % country)
    return row.loc[:, "2018"].values[0]

# normalization and smoothing
def per_capita(locd, mult=mult, criteria=critlow):
    """
    calculate per capita cases (confirmed, deaths, and recovered)

    """
    pop = locd["population"]
    locd["cnf_pc"] = locd["cnf"]/pop*mult
    locd["dth_pc"] = locd["dth"]/pop*mult
    if "rec" in locd.keys():
        locd["rec_pc"] = locd["rec"]/pop*mult
    return

def time_zero(dates, locd, mult=mult, criteria=critlow):
    """
    shift elapsed time to specified number of _smoothed per-capita cases_
    """
    criteval = locd["cnf_pc_h"] > criteria*mult # use confirmed metric for time zero
    #criteval = locd["dth_pc"] > criteria*mult # use death metric for time zero
    idx0 = np.nonzero(criteval)[0][0]
    locd["idx0"] = idx0
    locd["days"] = (dates - dates[idx0]).astype('timedelta64[D]').values
    return

def perform_operations(corona):
    """
    perform the smoothing, per capita, time shift, and derivative operations
    """
    
    locs = corona["locs"]
    dates = corona["dates"]
    ndays = dates.size

    days = (dates - dates[0]).astype('timedelta64[D]').values
    lmbd=5e-5 # smoothing parameter
    for loc in locs:
        locd = corona[loc]
        # compute per-capita
        per_capita(locd)
        # smoothing
        if False:
            # smoothing, unconstrained
            locd['cnf_pc_h'] = ds.smooth_data(days, locd['cnf_pc'], d=2, lmbd=lmbd)
            locd['dth_pc_h'] = ds.smooth_data(days, locd['dth_pc'], d=2, lmbd=lmbd)
        else:
            # smoothing, with constraints
            # constrain cases to be non-negative
            Aones = -np.eye(ndays)
            bzero = np.zeros((ndays,1))
            # constrain cases to be always increasing
            D = -ds.derivative_matrix(days,1)
            bdzero = np.zeros((ndays-1,1))
            Aiq = np.vstack((Aones,D))
            biq = np.vstack((bzero, bdzero))
            locd['cnf_pc_h'] = ds.smooth_data_constr(days, locd['cnf_pc'], 2, lmbd,
                                                      (Aiq,biq))
            locd['dth_pc_h'] = ds.smooth_data_constr(days, locd['dth_pc'], 2, lmbd,
                                                      (Aiq,biq))
        # set time-zero and create elapsed time
        time_zero(dates, locd)

    # determine rates
    for loc in locs:
        locd = corona[loc]
        # take the derivative
        locd['cnf_rate'] = deriv1_fd(locd['cnf_pc_h'], locd['days'], central=True)
        locd['dth_rate'] = deriv1_fd(locd['dth_pc_h'], locd['days'], central=True)


def is_even(val):
    return np.mod(val, 2) == 0

def deriv1_fd(y,x=None,h=None,xmid=False,central=False,e=0.0):
    """
    Computes the approximate derivative of the y vs x data by finite
    difference.

    Parameters
    ___________________
    y : 1D array
        y data values of array length N
    x : 1D array
        monotonically increasing x data values of array length N; need
        not be provided if equally spaced and h is provided.
    h : float
        Uniform grid spacing of x.  If x is not equally spaced, x
        should be provided and h should be None.
    xmid : boolean
        If true, then midpoint values of x are returned along with the
        the derivative values
    central : boolean
        Performs a central difference, correct for unequally spaced x
        values. The endpoints are calculate by forward and backward
        difference so that yp is the same length as x.
    e : float
        Backward/forward adjustment to the central difference
        calculation: -1<e<0 shift to backward differnce, 0<e<1 shift
        to forward difference.
        
    Returns
    _______
    yp : 1D array
        The derivative values; array length N-1 (length N if
        central=True)
    xm : 1D array (optional)
        The midpoint values of x; array length N-1

    see also:  numpy.gradient
    """
    # Jonathan Stickel 2010, 2011, 2012, 2013
    
    # TODO:
    # - in 'deriv1_fd', calculate second-order accurate forward and backward
    # difference for endpoints (currently only first-order accurate)

    if xmid and central:
        raise Exception('Both xmid and central cannot both be True')

    # calculate direct differences
    if h is not None:
        # this calculation could potentially be avoided if doing central difference
        yp = np.diff(y)/h
    else:
        yp = np.diff(y)/np.diff(x.astype(float))
    
    if central:
        yd21 = np.diff(y[::2]) # every-other difference of y[i+1] - y[i-1]
                               # starting at i = 1
        yd22 = np.diff(y[1::2]) # every-other difference of y[i+1] - y[i-1]
                                # starting at i = 2
        if is_even(y.size):
            # interleave the two 1d arrays
            yd2 = np.reshape(np.column_stack((yd21,yd22)),-1)
        else:
            yd22 = np.hstack((yd22,np.nan))
            yd2 = np.reshape(np.column_stack((yd21,yd22)),-1)
            yd2 = yd2[:-1]
        if h is not None:
            ypc = yd2/(2*h)
        else:
            xd21 = np.diff(x[::2])
            xd22 = np.diff(x[1::2])
            if is_even(x.size):
                xd2 = np.reshape(np.column_stack((xd21,xd22)),-1)
            else:
                xd22 = np.hstack((xd22,np.nan))
                xd2 = np.reshape(np.column_stack((xd21,xd22)),-1)
                xd2 = xd2[:-1]
            ypc = yd2/xd2.astype(float) # central difference
        # paramater e can be used to shift towards backward/forward difference
        ypc = (1-e)*yp[:-1] + (1+e)*yp[1:] - ypc
        return np.hstack( ( yp[0], ypc, yp[-1] ) )
    else:
        if xmid:
            if x is None:
                raise Exception('x must be provided if xmid=True')
            if h is not None:
                xm = x[:-1] + h/2.0
            else:
                xm = ( x[:-1] + x[1:] )/2.0
            return yp,xm
        else:
            return yp

