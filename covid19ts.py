"""
Module with functions for reading and processing COVID-19 timeseries data 

https://github.com/CSSEGISandData/COVID-19
https://covidactnow.org
https://data.worldbank.org/indicator/sp.pop.totl
https://github.com/kjhealy/fips-codes/blob/master/county_fips_master.csv 
"""

# old usage data sets
# - https://www.census.gov/data/tables/time-series/demo/popest/2010s-state-total.html (US populations)
# - http://www.healthdata.org/ (hospital capacity)
# - https://www.census.gov/data/tables/time-series/demo/popest/2010s-state-total.html


## TODO:
# - (see list in `analyze_covid_time_series.py`)

import numpy as np
import pandas as pd
import json
import urllib.request
import os.path
import datetime
#import regularsmooth as ds # my local copy
import scikits.datasmooth as ds # pip installed

popfile = ("API_SP.POP.TOTL_DS2_en_csv_v2_2763937/"
           "API_SP.POP.TOTL_DS2_en_csv_v2_2763937.csv")

# my covid act now api key after registering
apikey = "d4e3714137fb41eb93d0d17fab7e9387"

##### scaling factor for cases, i.e., "x per mult" ######
#multval=1e6
multval=1e4


def covid19_global(countries, websource=True, JHCSSEpath=None, file_pop=popfile,
                   mult=multval, lmbd=5e-5, dbf=None, nsub=1):
    """
    Read in Johns Hopkins CSSE COVID-19 timeseries data for the countries
    specified and perform these operations:
    - normalize data to be per capita
    - smooth the data
    - set time-zero for each location and shift elapsed time in days
    - determine rates (i.e., the derivative) for data

    Optional arguments:
    websource:      whether to get time-series data from the web source or local files
    JHCSSEpath:     path to the local files if websource=False
    file_pop:       path to the population file (default is provided)
    mult:           scaling factor for the data, default is 1e6
    lmbd:           smoothing parameter
    dbf:            days before latest record to include in the analysis
    nsub:           subsampling period
    """
    
    # read in Johns Hopkins' data tables
    data = {}
    webbase = "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/"
    for dataname in ["confirmed", "deaths", "recovered"]:
        filename = "time_series_covid19_" + dataname + "_global.csv"
        webpath = webbase + filename
        if websource:
            print("Using web data file %s" % webpath)
            data[dataname] = pd.read_csv(webpath)
        else:
            filepath = JHCSSEpath + filename
            try:
                mdate = datetime.date.fromtimestamp(os.path.getmtime(filepath))
                today = datetime.date.today()
                if (today-mdate).days > 1:
                    raise ValueError("File is older than 1 day") # will this ever print?
                else:
                    print("Using existing file %s" % filepath)
            except:
                print("Downloading %s" % filepath)
                urllib.request.urlretrieve(webpath, filepath)
            data[dataname] = pd.read_csv(filepath)

    # convert some pesky NaN values (for missing labels) to strings
    for key, dataset in data.items():
        dataset["Province/State"] = dataset["Province/State"].values.astype("str")

    # presume these files all have the same "fixed" structure -- could put in a
    # check at some point
    header = data["confirmed"].columns.values
    dates = pd.to_datetime(header[4:])
    dates = np.flip(np.flip(dates)[:dbf:nsub])

    # read in population data (source: https://data.worldbank.org/indicator/sp.pop.totl)
    data_pop = pd.read_csv(file_pop, header=4)

    corona = dict()
    corona["locs"] = countries # I'm not sure if I like this approach... JJS 4/15/20
    for country in countries:
        ctryd = dict()
        corona[country] = ctryd
        ctryd['name'] = country
        ctryd['population'] = read_pop(data_pop, country)
        arr = read_cases(data["confirmed"], country)
        ctryd['cnf'] = np.flip(np.flip(arr)[:dbf:nsub])
        arr = read_cases(data["deaths"], country)
        ctryd['dth'] = np.flip(np.flip(arr)[:dbf:nsub])
        arr = read_cases(data["recovered"], country)
        ctryd['rec'] = np.flip(np.flip(arr)[:dbf:nsub])

    corona["mult"] = mult
    corona["dates"] = dates

    perform_operations(corona, lmbd, mult)
    
    return corona


def covid19_US(locations, websource=True, JHCSSEpath=None, mult=multval, lmbd=5e-5,
               dbf=None, nsub=1):
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
    lmbd:           smoothing parameter
    dbf:            days before latest record to include in the analysis
    nsub:           subsampling period
    """
    
    # read in Johns Hopkins' data tables
    data = {}
    webbase = "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/"
    for dataname in ["confirmed", "deaths"]:
        filename = "time_series_covid19_" + dataname + "_US.csv"
        webpath = webbase + filename
        if websource:
            print("Using web data file %s" % webpath)
            data[dataname] = pd.read_csv(webpath)
        else:
            filepath = JHCSSEpath + filename
            try:
                mdate = datetime.date.fromtimestamp(os.path.getmtime(filepath))
                today = datetime.date.today()
                if (today-mdate).days > 1:
                    raise ValueError("File is older than 1 day") # will this ever print?
                else:
                    print("Using existing file %s" % filepath)
            except:
                print("Downloading %s" % filepath)
                urllib.request.urlretrieve(webpath, filepath)
            data[dataname] = pd.read_csv(filepath)
    
    # presume these files all have a known "fixed" structure -- could put in a
    # check at some point
    header = data["confirmed"].columns.values
    dates = pd.to_datetime(header[11:])
    dates = np.flip(np.flip(dates)[:dbf:nsub])

    corona = dict()
    corona["locs"] = locations # I'm not sure if I like this approach... JJS 4/15/20
    for loc in locations:
        locd = dict()
        corona[loc] = locd
        locd['name'] = loc
        # perform data handling differently by whether location is a state or
        # county--if a state, then need to read in all the state entries and
        # sum
        if "," in loc: # it's a county
            idx = loc.find(",")
            county = loc[:idx]
            state = loc[idx+2:]
            loc_bool = np.logical_and(data["deaths"]["Admin2"] == county,
                                      data["deaths"]["Province_State"] == state)
            if loc_bool.sum() != 1:
                raise Exception("%s is not in the data files, or occurs more than once"
                                % loc)
        else: # it's a state
            loc_bool = data["deaths"]["Province_State"] == loc
            if loc_bool.sum() < 1:
                raise Exception("%s is not in the data files" % loc)
        # these lines should work even for county location with single line...
        locd["population"] = data["deaths"][loc_bool]["Population"].sum()
        arr = data["confirmed"][loc_bool].iloc[:, 11:].values.sum(axis=0)
        locd["cnf"] = np.flip(np.flip(arr)[:dbf:nsub])
        arr = data["deaths"][loc_bool].iloc[:, 12:].values.sum(axis=0)
        locd["dth"] = np.flip(np.flip(arr)[:dbf:nsub])

    corona["mult"] = mult
    corona["dates"] = dates

    perform_operations(corona, lmbd, mult)
    
    return corona


def covid19_can(locations, lastday, websource=True, sourcepath=None, mult=multval, lmbd=5e-5,
                dbf=None, nsub=1):
    """Read in Covid Act Now COVID-19 timeseries data for the US locations
    specified and perform these operations:
    - normalize data to be per capita
    - smooth the data
    - set time-zero for each location and shift elapsed time in days
    - determine rates (i.e., the derivative) 
    
    `locations` argument should be a list of location names in the form of
    `US`, `state`, OR `county, state`. [county locations still TBD 8/26/21]
 
    `lastday` is arbitrary reference day, but intended to be the last day of collected data
    
    Optional arguments:
    websource:      whether to get time-series data from the web or local files
    path:           path to the local files if websource=False; filenames are
                    expected to be in the form `[CD].timeseries.json`, where
                    `CD` is one of `US`, the two letter code for the state, or
                    FIPS county code.
    mult:           scaling factor for the data, default is 1e6
    lmbd:           smoothing parameter
    nsub:           subsampling period

    """
    # to read state codes
    codetable = pd.read_json("data/state_codes.json")
    codetable.set_index("State", inplace=True)

    # read fips codes
    # encoding problem, determined by opening in a text-editor -- can also be fixed by
    # opening in LibreOffice and resaving
    fipstable = pd.read_csv("data/county_fips_master.csv", encoding="iso8859_15")
    fipstable.set_index("long_name", inplace=True)

    # # hospital bed capacities -- should be in location-specific files now, or
    # # at least icu-beds seem to be, 8/26/21
    # hosp_data = pd.read_csv("Summary_stats_all_locs.csv")
    # hosp_data.set_index("location_name", inplace=True)
   
    # read state data
    nlocs = len(locations)
    ### FIXME, extract state part of location, and use FIPS codes for counties,
    ### from `county_fips_maser.csv`
    # codes = codetable.loc[locations]["Code"].values # for an array, no loop
    codes = []
    for loc in locations:
        if "," in loc: # it's a county
            idx = loc.find(",")
            county = loc[:idx]
            state = loc[idx+2:]
            #stateAbr = codetable.loc[state]["Code"] # will expect state code instead
            # get fips code as a string with 5 digits, prepending 0 if needed
            fipscode = "%05i" % fipstable.loc[county + " " + state]["fips"]
            codes.append(fipscode)
        elif loc == "US":
            codes.append("US")
        else:
            codes.append( codetable.loc[loc]["Code"] )

    data_locs = dict()
    for i in range(nlocs):
        if codes[i]=="US":
            webbase = "https://api.covidactnow.org/v2/country/"
        elif len(codes[i])==2: # a state
            webbase = "https://api.covidactnow.org/v2/state/"
        else:
            webbase = "https://api.covidactnow.org/v2/county/"
        webpath = webbase + codes[i] + ".timeseries.json?apiKey=" + apikey
        #print(webpath)
        if websource:
            print("Using web data file %s" % webpath[:-40]) # strip the apikey
            data_locs[locations[i]] = json.load(urllib.request.urlopen(webpath))
        else:
            filepath = sourcepath + codes[i] + ".timeseries.json"
            try:
                mdate = datetime.date.fromtimestamp(os.path.getmtime(filepath))
                today = datetime.date.today()
                if (today-mdate).days > 1:
                    raise ValueError("File is older than 1 day") # will this ever print?
                else:
                    print("Using existing file %s" % filepath)
            except:
                print("Downloading %s" % filepath)
                urllib.request.urlretrieve(webpath, filepath)
            data_locs[locations[i]] = json.load(open(filepath))  

    # set of relevant keys, from the set in "actualsTimeseries"
    # ['cases', 'deaths', 'positiveTests', 'negativeTests',
    #    'contactTracers', 'hospitalBeds', 'icuBeds', 'newCases',
    #    'newDeaths', 'vaccinesAdministeredDemographics',
    #    'vaccinationsInitiatedDemographics', 'date', 'vaccinesDistributed',
    #    'vaccinationsInitiated', 'vaccinationsCompleted',
    #    'vaccinesAdministered']
    # note that `hospitalBeds` and `icuBeds` expand further into
    # ['capacity', 'currentUsageTotal', 'currentUsageCovid','typicalUsageRate']
    keys = ['cases', 'deaths', 'positiveTests', 'negativeTests',
            'hospCapacity', 'hospTotal', 'hospCovid', 'hospTypical',
            'icuCapacity', 'icuTotal', 'icuCovid', 'icuTypical',
            'vaccinesDistributed', 'vaccinationsInitiated', 'vaccinationsCompleted',
            'vaccinesAdministered']

    # initiate combined data dictionary
    corona = dict()
    corona["locs"] = locations
    corona["lastday"] = lastday
    corona["mult"] = mult

    # process data
    for loc in locations:
        locd = dict()
        locd["name"] = loc
        data = data_locs[loc]

        # get the datafram of the timeseries data
        data_df = pd.DataFrame(data["actualsTimeseries"])
        # remove the time-series objects from data (not entirely necessary, but
        # will save some ram and maybe some computational cost)
        for key in ["actualsTimeseries", "metricsTimeseries", "riskLevelsTimeseries"]:
            data.pop(key)
        # expand hospital and icu columns
        hosp_df = pd.DataFrame(data_df["hospitalBeds"].values.tolist())
        hosp_df.rename(columns={"capacity": "hospCapacity",
                               "currentUsageTotal": "hospTotal",
                               "currentUsageCovid": "hospCovid",
                                "typicalUsageRate": "hospTypical"}, inplace=True)
        icu_df = pd.DataFrame(data_df["icuBeds"].values.tolist())
        icu_df.rename(columns={"capacity": "icuCapacity",
                               "currentUsageTotal": "icuTotal",
                               "currentUsageCovid": "icuCovid",
                               "typicalUsageRate": "icuTypical"}, inplace=True)
        # remove hospital and icu bed columns that are actually elements of dictionaries
        data_df = data_df.drop(columns=["hospitalBeds", "icuBeds"])
        # add back hospital and icu bed columns
        data_df = pd.concat([data_df, hosp_df, icu_df], axis=1)
        # convert None values to NaN
        data_df.fillna(value=np.nan, inplace=True)

        # get population
        pop = data["population"]
        locd["population"] = pop
        
        # hospital capacity seems to be intended but missing from the data
        # files -- however, icu beds are there, so will use that, JJS 8/26/21
        # capacity = hosp_data["all_bed_capacity"][loc] 
        
        # get and process dates -- US and each sub-location may have different dates
        dates = pd.to_datetime(data_df["date"].values)
        dates = np.flip(np.flip(dates)[:dbf:nsub])
        ndays = dates.size
        days = (dates - lastday).astype('timedelta64[D]')
        locd["dates"] = dates
        locd["days"] = days
        
        # interpolate nan values, subsample, and transfer meaningful data
        for key in keys:
            # arr = data_df[key].values # without interpolation
            arr = data_df[key].interpolate().values
            locd[key] = np.flip(np.flip(arr)[:dbf:nsub])

        # compute per capita for metrics *I* want to look at (more may be added)
        locd["cnf_pc"] = locd["cases"]/pop*mult
        locd["dth_pc"] = locd["deaths"]/pop*mult
        # hospital bed capacity doesn't seem to be populated in the data, so
        # only getting per-capita hospital usage, JJS 8/27/21
        locd["hosp_total_pc"] = locd["hospTotal"]/pop*mult
        locd["hosp_covid_pc"] = locd["hospCovid"]/pop*mult
        locd["icu_capacity"] = locd["icuCapacity"]
        locd["icu_total_pc"] = locd["icuTotal"]/pop*mult
        locd["icu_total_frac"] = locd["icuTotal"]/locd["icu_capacity"]
        locd["icu_covid_pc"] = locd["icuCovid"]/pop*mult
        locd["icu_covid_frac"] = locd["icuCovid"]/locd["icu_capacity"]
        locd["pos_pc"] = locd["positiveTests"]/pop*mult
        locd["neg_pc"] = locd["negativeTests"]/pop*mult
        locd["test_frac"] = (locd["positiveTests"] + locd["negativeTests"])/pop
        locd["vacc_init_frac"] = locd["vaccinationsInitiated"]/pop
        locd["vacc_full_frac"] = locd["vaccinationsCompleted"]/pop
        
        # CFR - calculate from smooth values instead?
        locd["cfr"] = locd["deaths"]/locd["cases"]
        
        # smoothing, with constraints
        # TODO: rework `perform_operations` to use it with this data also? JJS
        # 8/27/21
        # constrain cases to be non-negative
        Aones = -np.eye(ndays)
        bzero = np.zeros((ndays,1))
        # constrain cases to be always increasing
        D = -ds.derivative_matrix(days,1)
        bdzero = np.zeros((ndays-1,1))
        Aiq = np.vstack((Aones,D))
        biq = np.vstack((bzero, bdzero))
        # filter for nan 
        cnf_pc = locd["cnf_pc"].copy()
        keep = ~np.isnan(cnf_pc)
        cnf_pc = cnf_pc[keep]
        dayscnf = days[keep]
        dth_pc = locd["dth_pc"].copy()
        keep = ~np.isnan(dth_pc)
        dth_pc = dth_pc[keep]
        daysdth = days[keep]
        
        locd['cnf_pc_h'] = ds.smooth_data_constr(dayscnf, cnf_pc, 2, lmbd,
                                                 (Aiq,biq), xhat=days)
        locd['dth_pc_h'] = ds.smooth_data_constr(daysdth, dth_pc, 2, lmbd,
                                                 (Aiq,biq), xhat=days)

        # determine rates, i.e., take the derivative
        locd['cnf_rate'] = deriv1_fd(locd['cnf_pc_h'], days, central=True)
        locd['dth_rate'] = deriv1_fd(locd['dth_pc_h'], days, central=True)

        corona[loc] = locd
    return corona


# # this was needed for covid-tracking-project and is no longer used, JJS 8/27/21
# def int_to_date(dateint):
#     """ convert integer of the form yyyymmdd to date object """
#     datestr = str(dateint)
#     return np.datetime64(datestr[:4] + "-" + datestr[4:6] + "-" + datestr[6:8])


def read_cases(data, country):
    """
    Select country data values. The country data must be in a single row for the
    current implementation. Not used for US locations.
    """
    # get the row for the country
    c_bool = np.logical_and(data["Province/State"] == "nan",
                            data["Country/Region"] == country)
    c_data = data[c_bool]
    if c_data.shape[0] == 1:
        return np.squeeze(c_data.iloc[:,4:].values)
    else:
        raise Exception("Not implemented:  there is more than one row (or no rows) of %s data." % country)
    
def read_pop(data, country):
    """ Get the population for a country. Not used for US locations. """
    # There may be some name mismatches that need correction in the population
    # file -- please create an issue or pull request when you find them
    row = data[ data["Country Name"]==country ]
    if row.size == 0:
        raise Exception("%s is not in the population data file" % country)
    return row.loc[:, "2020"].values[0]

# normalization
def per_capita(locd, mult):
    """
    calculate per capita cases (confirmed, deaths, and recovered)

    """
    pop = locd["population"]
    locd["cnf_pc"] = locd["cnf"]/pop*mult
    locd["dth_pc"] = locd["dth"]/pop*mult
    if "rec" in locd.keys():
        locd["rec_pc"] = locd["rec"]/pop*mult
    return


def perform_operations(corona, lmbd, mult):
    """
    Perform the smoothing, per capita, time shift, and derivative operations.
    The argument `lmbd` is the smoothing parameter

    """
    
    locs = corona["locs"]
    dates = corona["dates"]
    ndays = dates.size

    days = (dates - dates[-1]).astype('timedelta64[D]').values
    corona["days"] = days
    
    for loc in locs:
        locd = corona[loc]
        # compute per-capita
        per_capita(locd, mult)
        # smoothing
        if False:
            # smoothing, unconstrained
            locd['cnf_pc_h'] = ds.smooth_data(days, locd['cnf_pc'], d=2, lmbd=lmbd)
            locd['dth_pc_h'] = ds.smooth_data(days, locd['dth_pc'], d=2, lmbd=lmbd)
        else:
            ## could pull this section out into a function
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
        #time_zero(dates, locd)

    # determine rates
    for loc in locs:
        locd = corona[loc]
        # take the derivative
        locd['cnf_rate'] = deriv1_fd(locd['cnf_pc_h'], days, central=True)
        locd['dth_rate'] = deriv1_fd(locd['dth_pc_h'], days, central=True)


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

