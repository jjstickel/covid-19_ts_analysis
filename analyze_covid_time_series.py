"""
Read in Johns Hopkins CSSE COVID-19 timeseries data (source:
https://github.com/CSSEGISandData/COVID-19) and perform some basic analysis

"""

# Jonathan Stickel, 2020

# on 3/23/20, J-H switched to a new set of csv files

# TODO:

#   scripting and notebook
# - Process countries with multiple entries. Will need to make sure that the
#   sum provides the correct result
# - estimate current number sick (confirmed - deaths - recovered)
# - (long term) switch from dict to class
# - add ability to read US data file (by state, city)

import numpy as np
from decimal import Decimal
from datetime import datetime
from matplotlib.pyplot import *
import matplotlib.dates as mdates
from covid19ts import covid19_global


## User input -- put up to 5 countries of interest in this list. Must be the
## same name used in the JH global files, and (at the moment), it must be a
## single entry in the file (e.g., China has multiple entries and will cause an
## Exception)
countries = ["US", "Italy", "Spain", "Iran"]

JHCSSEpath = "../JH_COVID-19/csse_covid_19_data/csse_covid_19_time_series/"

corona = covid19_global(countries, websource=False, JHCSSEpath=JHCSSEpath)
mult = corona["mult"]
critlow = corona["critlow"]
nctry = len(countries)
dates = corona["dates"]

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
    
### plotting ###
ion()
clr = ['C%g' % i for i in range(10)]
sbl = ["o", "s", "v", "d", "x"]
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
ylabel("case-fatality ratio [%]")
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

