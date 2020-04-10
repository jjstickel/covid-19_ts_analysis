"""
Read in Johns Hopkins' COVID-19 timeseries data (source:
https://github.com/CSSEGISandData/COVID-19) and perform some basic analysis

"""

# Jonathan Stickel, 2020

# on 3/23/20, J-H switched to a new set of csv files

# TODO:

# - make a module for all the data processing that can be used for both
#   scripting and notebook
# - Process countries with multiple entries. Will need to make sure that the
#   sum provides the correct result
# - estimate current number sick (confirmed - deaths - recovered)
# - (long term) switch from dict to class
# - add ability to read US data file (by state, city)

import numpy as np
import pandas as pd
from matplotlib.pyplot import *
import matplotlib.dates as mdates
from datetime import datetime
from decimal import Decimal

from findiffjs import deriv1_fd
#import regularsmooth as ds # my local copy
import scikits.datasmooth as ds # pip installed


## User input -- put up to 5 countries of interest in this list. Must be the
## same name used in the JH global files, and (at the moment), it must be a
## single entry in the file (e.g., China has multiple entries and will cause an
## Exception)
countries = ["US", "Italy", "Spain", "Iran"]
# flag for using web or local address for the Johns Hopkins CSSE data
websource = False

# read in Johns Hopkins' data tables
if websource:
    pathname = "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/"
else:
    pathname = "../JH_COVID-19/csse_covid_19_data/csse_covid_19_time_series/"
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
file_pop = "API_SP.POP.TOTL_DS2_en_csv_v2_887275/API_SP.POP.TOTL_DS2_en_csv_v2_887275_2018.csv"
data_pop = pd.read_csv(file_pop, header=4)

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
    row = data_pop[ data_pop["Country Name"]==country ]
    if row.size == 0:
        raise Exception("%s is not in the population data file" % country)
    return row.loc[:, "2018"].values[0]

nctry = len(countries)
corona = dict()
for country in countries:
    ctryd = dict()
    corona[country] = ctryd
    ctryd['name'] = country
    ctryd['population'] = read_pop(data_pop, country)
    ctryd['cnf'] = read_cases(data_confirmed, country)
    ctryd['dth'] = read_cases(data_deaths, country)
    ctryd['rec'] = read_cases(data_recovered, country)

##### scaling factor for cases, i.e., "x per mult" ######
mult=1e6
critlow = 10*1e-6 # for time zero, using confirmed
#critlow = 0.1*1e-6 # for time zero, using deaths

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
        ctryd['cnf_pc_h'] = ds.smooth_data_constr(days, ctryd['cnf_pc'], 2, lmbd, (Aiq,biq))
        ctryd['dth_pc_h'] = ds.smooth_data_constr(days, ctryd['dth_pc'], 2, lmbd, (Aiq,biq))
    # set time-zero and create elapsed time
    time_zero(dates, ctryd)
    
# compute exponential fit 
def expfit(t, y):
    """
    fit data to exponential function
    """
    k, lna = np.polyfit(t, np.log(y), 1)
    return np.exp(lna), k

crithigh = 200*1e-6 # upper bound for exponential fit, for fitting confirmed
#crithigh = 2*1e-6 # upper bound for exponential fit, for fitting deaths
for country in countries:
    ctryd = corona[country]
    criteval = ctryd['cnf_pc'] < crithigh*mult
    idx0 = ctryd['idx0']
    idx1 = np.nonzero(criteval)[0][-1]
    a, k = expfit(ctryd['days'][idx0:idx1], ctryd['cnf_pc'][idx0:idx1])
    dt2in = np.log(2)/k
    print("initial doubling time for %s was %g days" % (country, dt2in))
    ctryd['a'] = a
    ctryd['k'] = k
    ctryd['dt2in'] = dt2in
    ctryd['cnf_expfit'] = a*np.exp(k*ctryd['days'])
    
# determine rates
for country in countries:
    ctryd = corona[country]
    # take the derivative
    ctryd['cnf_rate'] = deriv1_fd(ctryd['cnf_pc_h'], ctryd['days'], central=True)
    ctryd['dth_rate'] = deriv1_fd(ctryd['dth_pc_h'], ctryd['days'], central=True)

### plotting ###
ion()
clr = ['C%g' % i for i in range(10)]
sbl = ["o", "s", "v", "d", "h"]
savefigs = False

# get inverse human readable form for t=0 criterium
clinv_tup = Decimal("%.1g"%(1/critlow)).as_tuple()
clinv_dig = clinv_tup.digits[0]
clinv_exp = clinv_tup.exponent

N = 0

# total confirmed, not scaled
N+=1; figure(N)
#figure(1, figsize=(4,3))
clf()
for i in range(nctry):
    ctryd = corona[countries[i]]
    population = ctryd['population']
    plot(dates, ctryd['cnf'], sbl[i]+clr[i], mfc='none', mew=1.5, label=ctryd["name"])
    #plot(dates, ctryd['cnf_expfit']*population/mult, '-'+clr[i], lw=1.5)
    plot(dates, ctryd['cnf_pc_h']*population/mult, '--'+clr[i])
ax = gca()
ax.get_xaxis().set_major_formatter(mdates.DateFormatter("%m/%d"))
setp(ax.get_xticklabels(), rotation=30, ha="right")
cnf_max = max([corona[country]['cnf'].max() for country in countries])
axis(xmin = datetime(year=2020, month=2, day = 15), ymin=0-cnf_max*0.1, ymax = cnf_max*1.1)
xlabel("date [m/d]")
ylabel("confirmed")
legend(loc='best')
if savefigs:
    savefig("confirmed.pdf", bbox_inches="tight")
    savefig("confirmed.png", bbox_inches="tight")

# confirmed per capita
N+=1; figure(N)
clf()
for i in range(nctry):
    ctryd = corona[countries[i]]
    days = ctryd["days"]
    plot(days, ctryd["cnf_pc"], sbl[i]+clr[i], mfc='none', mew=1.5, label=ctryd["name"])
    plot(days, ctryd["cnf_expfit"], '-'+clr[i], lw=1.5)
    plot(days, ctryd["cnf_pc_h"], '--'+clr[i])
scaled_max = max([corona[country]['cnf_pc'].max() for country in countries])
axis(xmin=-5, ymin=0-scaled_max*0.1, ymax=scaled_max*1.1)
xlabel('days since %i confirmed per $10^%i$' % (clinv_dig, clinv_exp))
#xlabel('days since %i deaths per $10^%i$' % (clinv_dig, clinv_exp))
#ylabel("confirmed per %1.0e" % mult)
ylabel("confirmed per $10^%i$" % np.log10(mult))
legend(loc='best')
if savefigs:
    savefig("confirmed_scaled.pdf", bbox_inches="tight")
    savefig("confirmed_scaled.png", bbox_inches="tight")

# log confirmed per capita
N+=1; figure(N)
clf()
for i in range(nctry):
    ctryd = corona[countries[i]]
    semilogy(ctryd["days"], ctryd["cnf_pc"], sbl[i]+clr[i], mfc='none', mew=1.5,
             label=ctryd["name"])
    semilogy(ctryd["days"], ctryd["cnf_expfit"], '-'+clr[i], lw=1.5)
    semilogy(ctryd["days"], ctryd["cnf_pc_h"], '--'+clr[i])
axis(xmin=-5, ymin=1e-3)
xlabel('days since %i confirmed per $10^%i$' % (clinv_dig, clinv_exp))
#xlabel('days since %i deaths per $10^%i$' % (clinv_dig, clinv_exp))
#ylabel("confirmed per %1.0e" % mult)
ylabel("confirmed per $10^%i$" % np.log10(mult))
legend(loc='best')

# rate confirmed per capita
N+=1; figure(N)
clf()
for i in range(nctry):
    ctryd = corona[countries[i]]
    plot(ctryd["days"], ctryd["cnf_rate"], "-"+sbl[i]+clr[i], mfc='none', mew=1.5,
         label=ctryd["name"])
axis(xmin = -5)
xlabel('days since %i confirmed per $10^%i$' % (clinv_dig, clinv_exp))
#xlabel('days since %i deaths per $10^%i$' % (clinv_dig, clinv_exp))
ylabel("rate [confirmed per $10^%i$ / day]" % np.log10(mult))
legend(loc="best")
if savefigs:
    savefig("confirmed_rate.png", bbox_inches="tight")

# total deaths, not scaled
N+=1; figure(N)
#figure(1, figsize=(4,3))
clf()
for i in range(nctry):
    ctryd = corona[countries[i]]
    population = ctryd['population']
    plot(dates, ctryd['dth'], sbl[i]+clr[i], mfc='none', mew=1.5, label=ctryd["name"])
    plot(dates, ctryd['dth_pc_h']*population/mult, '--'+clr[i])
ax = gca()
ax.get_xaxis().set_major_formatter(mdates.DateFormatter("%m/%d"))
setp(ax.get_xticklabels(), rotation=30, ha="right")
dth_max = max([corona[country]['dth'].max() for country in countries])
axis(xmin = datetime(year=2020, month=2, day = 15), ymin=0-dth_max*0.1, ymax = dth_max*1.1)
xlabel("date [m/d]")
ylabel("deaths")
legend(loc='best')
if savefigs:
    savefig("temp.pdf", bbox_inches="tight")
 
# deaths per capita
N+=1; figure(N)
clf()
for i in range(nctry):
    ctryd = corona[countries[i]]
    days = ctryd["days"]
    plot(days, ctryd["dth_pc"], sbl[i]+clr[i], mfc='none', mew=1.5, label=ctryd["name"])
    plot(days, ctryd["dth_pc_h"], '--'+clr[i])
scaled_max = max([corona[country]['dth_pc'].max() for country in countries])
axis(xmin=-5, ymin=0-scaled_max*0.1, ymax=scaled_max*1.1)
xlabel('days since %i confirmed per $10^%i$' % (clinv_dig, clinv_exp))
#xlabel('days since %i deaths per $10^%i$' % (clinv_dig, clinv_exp))
#ylabel("confirmed per %1.0e" % mult)
ylabel("deaths per $10^%i$" % np.log10(mult))
legend(loc='best')
if savefigs:
    savefig("temp.pdf", bbox_inches="tight")
    
# rate deaths per capita
N+=1; figure(N)
clf()
for i in range(nctry):
    ctryd = corona[countries[i]]
    plot(ctryd["days"], ctryd["dth_rate"], "-"+sbl[i]+clr[i], mfc='none', mew=1.5,
         label=ctryd["name"])
axis(xmin = -5)
xlabel('days since %i confirmed per $10^%i$' % (clinv_dig, clinv_exp))
#xlabel('days since %i deaths per $10^%i$' % (clinv_dig, clinv_exp))
ylabel("rate [deaths per $10^%i$ / day]" % np.log10(mult))
legend(loc="best")
if savefigs:
    savefig("temp.png", bbox_inches="tight")

# case fatality ratio
N+=1; figure(N)
clf()
for i in range(nctry):
    ctryd = corona[countries[i]]
    plot(ctryd["days"], ctryd["dth_pc_h"]/ctryd["cnf_pc_h"]*100, "-"+sbl[i]+clr[i], mfc='none',
         mew=1.5, label=ctryd["name"])
axis(xmin = -5)
xlabel('days since %i confirmed per $10^%i$' % (clinv_dig, clinv_exp))
#xlabel('days since %i deaths per $10^%i$' % (clinv_dig, clinv_exp))
ylabel("case-fatality ratio")
legend(loc="best")
if savefigs:
    savefig("temp.png", bbox_inches="tight")

# confirmed per-capita cases and deaths on one graph
N+=1; figure(N)
clf()
ax1 = subplot(111)
ax2 = ax1.twinx()
for i in range(nctry):
    ctryd = corona[countries[i]]
    days = ctryd["days"]
    ax1.plot(days, ctryd["cnf_pc"], "-"+sbl[i]+clr[i], mfc='none', mew=1.5,
             label=ctryd["name"])
    ax2.plot(days, ctryd["dth_pc"], ":"+sbl[i]+clr[i], mfc='none', mew=1.5, ms=3)
    #plot(days, ctryd["dth_pc_h"], '--'+clr[i])
cnf_max = max([corona[country]['cnf_pc'].max() for country in countries])
dth_max = max([corona[country]['dth_pc'].max() for country in countries])
ax1.axis(xmin=-5, ymin=0-cnf_max*0.1, ymax=cnf_max*1.1)
ax2.axis(xmin=-5, ymin=0-dth_max*0.1, ymax=dth_max*1.1)
ax1.set_xlabel('days since %i confirmed per $10^%i$' % (clinv_dig, clinv_exp))
#ax1.set_xlabel('days since %i deaths per $10^%i$' % (clinv_dig, clinv_exp))
ax1.set_ylabel("confirmed per $10^%i$" % np.log10(mult))
ax2.set_ylabel("deaths per $10^%i$" % np.log10(mult))
handles, labels = ax1.get_legend_handles_labels()
legend(handles, labels, loc='best')
if savefigs:
    savefig("temp.pdf", bbox_inches="tight")

