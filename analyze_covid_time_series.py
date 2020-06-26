"""
Read in Johns Hopkins CSSE COVID-19 timeseries data (source:
https://github.com/CSSEGISandData/COVID-19) and perform some basic analysis

"""

# Jonathan Stickel, 2020

# on 3/23/20, J-H switched to a new set of csv files

# TODO:
# - Process countries with multiple entries. Will need to make sure that the
#   sum provides the correct result
# - (long term) switch from dict to class

import numpy as np

from covid19ts import covid19_global, covid19_US
import covid_plots as cvp

## User input -- put up to 5 countries of interest in this list. Must be the
## same name used in the JH global files, and (at the moment), it must be a
## single entry in the file (e.g., China has multiple entries and will cause an
## Exception)
countries = ["US", "Italy", "Sweden", "Russia", "Brazil"]
#countries = ["US", "Sweden", "Denmark", "Norway"]

# now also process US locations
#US_locs = ["Colorado", "Florida", "Texas", "Arizona"]
#US_locs = ["Colorado", "South Dakota", "Minnesota", "Wisconsin"]
US_locs = ["Colorado", "Florida", "New York", "Arizona"]
#US_locs = ["Colorado", "Washington", "California", "New York"]
#US_locs = ["Colorado", "New York", "New York, New York"]

dbf = 90
saveplots = False

JHCSSEpath = "../JH_COVID-19/csse_covid_19_data/csse_covid_19_time_series/"

lmbd = 5e-5
corona = covid19_global(countries, websource=False, JHCSSEpath=JHCSSEpath, lmbd=lmbd)

mult = corona["mult"]
critlow = corona["critlow"]
nctry = len(countries)
dates = corona["dates"]

coronaUS = covid19_US(US_locs, websource=False, JHCSSEpath=JHCSSEpath, lmbd=lmbd)
if not np.alltrue(corona["dates"] == coronaUS["dates"]):
    raise ValueError("the dates from the global and US files do not match")
nUSloc = len(US_locs)

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

# estimate "active" cases, then recovered and compare to "recovered" data
# estimate of recovery time from:
# - https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(20)30243-7/fulltext
# - https://towardsdatascience.com/visual-notes-from-singapores-first-100-fully-recovered-covid-19-patients-aad7f2e1d0a0
rectime = 14 # days, 1 day = 1 data point
for country in countries:
    ctryd = corona[country]
    # calculate active cased based on JH data
    ctryd["acv_pc"] = ctryd["cnf_pc"] - ctryd["dth_pc"] - ctryd["rec_pc"]
    # estimate active cases based on an "average" recovery time
    ctryd["acvest_pc"] = ctryd["cnf_pc"].copy()
    ctryd["acvest_pc"][rectime:] = ctryd["acvest_pc"][rectime:] - ctryd["cnf_pc"][:-rectime]
    # recovered estimate -- not sure of the value of this...
    ctryd["recest_pc"] = ctryd["cnf_pc"] - ctryd["dth_pc"] - ctryd["acv_pc"]
for loc in US_locs:
    locd = coronaUS[loc]
    locd["acvest_pc"] = locd["cnf_pc"].copy()
    locd["acvest_pc"][rectime:] = locd["acvest_pc"][rectime:] - locd["cnf_pc"][:-rectime]


# plotting
cvp.critlow_readable(corona) # provide convenient readable terms for time labeling
N = 1
cvp.total_global_plot(corona, N)
N+=1
cvp.per_capita_global_plot(corona, N, savefigs=saveplots, days_before=dbf)
N+=1
cvp.rate_global_plot(corona, N, savefigs=saveplots, days_before=dbf)
N+=1
cvp.active_CFR_global_plot(corona, N, savefigs=saveplots, days_before=dbf)
N+=1
cvp.exp_fit_confirmed_plot(corona, N)
N+=1
cvp.confirmed_deaths_simul_global_plot(corona, N)
N+=1
cvp.per_capita_US_plot(coronaUS, corona, N, savefigs=saveplots, days_before=dbf)
N+=1
cvp.rate_US_plot(coronaUS, corona, N, savefigs=saveplots, days_before=dbf)
N+=1
cvp.active_CFR_US_plot(coronaUS, corona, N, savefigs=saveplots, days_before=dbf)
