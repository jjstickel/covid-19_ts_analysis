"""
Module to collect reading and processing of Johns Hopkins CSSE COVID-19
timeseries data 

https://github.com/CSSEGISandData/COVID-19)

"""


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
    websource:      whether to get time-series data from the web source or local file
    JHCSSEpath:     path to the local file if websource=False
    file_pop:       path to the population file (default is provided)
    mult:           scaling factor for the data, default is 1e6
    critlow:        cutoff criteria for time zero, default is 10e-6
    """
    
    # read in Johns Hopkins' data tables
    if websource:
        pathname = "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/"
    elif JHCSSEpath is None:
        raise ValueError("If websource=False, you need to provide a path for the JH CSSE data")
    else:
        pathname = JHCSSEpath
    # read in global data
    data_confirmed = pd.read_csv(pathname + "time_series_covid19_confirmed_global.csv")
    data_deaths = pd.read_csv(pathname + "time_series_covid19_deaths_global.csv")
    data_recovered = pd.read_csv(pathname + "time_series_covid19_recovered_global.csv")
    # presume these files all have the same "fixed" structure -- could put in a
    # check at some point
    header = data_confirmed.columns.values
    dates = pd.to_datetime(header[4:])
    ndays = dates.size

    # read in population data (source: https://data.worldbank.org/indicator/sp.pop.totl)
    data_pop = pd.read_csv(file_pop, header=4)

    corona = dict()
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
    
    days = (dates - dates[0]).astype('timedelta64[D]').values
    lmbd=5e-5 # smoothing parameter
    for country in countries:
        ctryd = corona[country]
        # compute per-capita
        per_capita(ctryd)
        # smoothing
        if False:
            # smoothing, unconstrained
            ctryd['cnf_pc_h'] = ds.smooth_data(days, ctryd['cnf_pc'], d=2, lmbd=lmbd)
            ctryd['dth_pc_h'] = ds.smooth_data(days, ctryd['dth_pc'], d=2, lmbd=lmbd)
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
            ctryd['cnf_pc_h'] = ds.smooth_data_constr(days, ctryd['cnf_pc'], 2, lmbd,
                                                      (Aiq,biq))
            ctryd['dth_pc_h'] = ds.smooth_data_constr(days, ctryd['dth_pc'], 2, lmbd,
                                                      (Aiq,biq))
        # set time-zero and create elapsed time
        time_zero(dates, ctryd)

    # determine rates
    for country in countries:
        ctryd = corona[country]
        # take the derivative
        ctryd['cnf_rate'] = deriv1_fd(ctryd['cnf_pc_h'], ctryd['days'], central=True)
        ctryd['dth_rate'] = deriv1_fd(ctryd['dth_pc_h'], ctryd['days'], central=True)

    return corona

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
def per_capita(ctryd, mult=mult, criteria=critlow):
    """
    calculate per capita cases (confirmed, deaths, and recovered)

    """
    pop = ctryd["population"]
    ctryd["cnf_pc"] = ctryd["cnf"]/pop*mult
    ctryd["dth_pc"] = ctryd["dth"]/pop*mult
    ctryd["rec_pc"] = ctryd["rec"]/pop*mult
    return

def time_zero(dates, ctryd, mult=mult, criteria=critlow):
    """
    shift elapsed time to specified number of _smoothed per-capita cases_
    """
    criteval = ctryd["cnf_pc_h"] > criteria*mult # use confirmed metric for time zero
    #criteval = ctryd["dth_pc"] > criteria*mult # use death metric for time zero
    idx0 = np.nonzero(criteval)[0][0]
    ctryd["idx0"] = idx0
    ctryd["days"] = (dates - dates[idx0]).astype('timedelta64[D]').values
    return




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

