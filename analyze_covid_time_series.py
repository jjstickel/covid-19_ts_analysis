"""Read in Johns Hopkins CSSE COVID-19 timeseries data (source:
https://github.com/CSSEGISandData/COVID-19) and Covid Tracking project
(https://covidtracking.com/) and perform some basic analysis

"""

# Jonathan Stickel, 2020, 2021

# on 3/23/20, J-H switched to a new set of csv files
# on 7/15/20, switched to using Covid Tracking Project for US data
# after March 2021, Covid Tracking Project stopped updating their data, so no longer using it

# US populations from
# https://www.census.gov/data/tables/time-series/demo/popest/2010s-state-total.html

# TODO:
# - Process countries with multiple entries. Will need to make sure that the
#   sum provides the correct result
# - (long term) switch from dict to class
# - check and correct for double counting of cities? e.g. New York City
# - switch to use other databases, especially for state data with
#   hospitalizations and vaccinations

import numpy as np

from covid19ts import covid19_global, covid19_US, covid19_ctp
import covid_plots as cvp

from matplotlib.pyplot import *
ion()

## User input -- put up to 7 countries of interest in this list. Must be the
## same name used in the JH global files, and (at the moment), it must be a
## single entry in the file (e.g., China has multiple entries and will cause an
## Exception)
#countries = ["US", "Italy", "Spain", "Sweden", "Brazil"]
#countries = ["US", "Italy", "Spain", "Germany", "Sweden", "Brazil", "India"]
countries = ["US", "Sweden", "Denmark", "Norway"]

# US states, up to 6; US will also be added automatically
#US_locs = ["Colorado", "California", "Arizona", "Florida", "Wisconsin", "South Dakota"]
US_locs = ["Colorado", "New York", "Arizona", "Florida", "California", "South Dakota"]
#US_locs = ["Colorado", "New York", "Florida", "Wisconsin", "North Dakota", "South Dakota"]
#US_locs = ["Colorado", "Washington", "California", "New York"]
#US_locs = ["Colorado", "New York", "New York, New York"]

dbf = None
#dbf = 150
nsub = 7 # subsample every `nsub` points
if (nsub > 14):
    raise Warning("Subsampling of %g is too large for estimating active cases" % nsub)

saveplots = False

JHCSSEpath = "../JH_COVID-19/csse_covid_19_data/csse_covid_19_time_series/"

lmbd = 5e-5
mult = 1e4
#mult = 1e6
corona = covid19_global(countries, websource=False, JHCSSEpath=JHCSSEpath, lmbd=lmbd,
                        mult=mult, dbf=dbf, nsub=nsub)

#mult = corona["mult"]
#critlow = corona["critlow"]
nctry = len(countries)
nUSloc = len(US_locs)

dates = corona["dates"]
lastday = dates[-1]  

## not using this data set for the time being -- will be interesting to check
## for differences with the new US data set
coronaUS = covid19_US(US_locs, websource=False, JHCSSEpath=JHCSSEpath, lmbd=lmbd, mult=mult,
                      dbf=dbf, nsub=nsub)
if not np.alltrue(corona["dates"] == coronaUS["dates"]):
    raise ValueError("the dates from the global and US files do not match")

# create function to analyze US data from the COVID Tracking Project; collect
# that data and analysis in it's own dict
# OBSOLETE -- COVID Tracking Project stopped collating data March 2021
#coronaUS_ctp = covid19_ctp(US_locs, lastday, websource=False,
#                           sourcepath="../covidtracking/", lmbd=lmbd, mult=mult, dbf=dbf)

# estimate "active" cases; since data for recovered cases is so unreliable,
# just this estimate is used
# estimate of recovery time from:
# - https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(20)30243-7/fulltext
# - https://towardsdatascience.com/visual-notes-from-singapores-first-100-fully-recovered-covid-19-patients-aad7f2e1d0a0
rectime = 14 # days, 1 day = 1 original data point
recsamp = np.int(14/nsub)
for country in countries:
    ctryd = corona[country]
    ctryd["acvest_pc"] = ctryd["cnf_pc"].copy()
    ctryd["acvest_pc"][:recsamp] = np.nan
    ctryd["acvest_pc"][recsamp:] = ctryd["acvest_pc"][recsamp:] - ctryd["cnf_pc"][:-recsamp]
for loc in coronaUS["locs"]:
    locd = coronaUS[loc]
    locd["acvest_pc"] = locd["cnf_pc"].copy()
    locd["acvest_pc"][:recsamp] = np.nan
    locd["acvest_pc"][recsamp:] = locd["acvest_pc"][recsamp:] - locd["cnf_pc"][:-recsamp]
    

# plotting
N = 1
#cvp.total_global_plot(corona, N)
#N+=1
cvp.per_capita_global_plot(corona, lastday, N, savefigs=saveplots, days_before=dbf)
N+=1
cvp.rate_global_plot(corona, lastday, N, savefigs=saveplots, days_before=dbf)
N+=1
cvp.active_CFR_global_plot(corona, lastday, N, savefigs=saveplots, days_before=dbf)
#N+=1
#cvp.exp_fit_confirmed_plot(corona, N)
#N+=1
#cvp.confirmed_deaths_simul_global_plot(corona, N)
N+=1
cvp.per_capita_US_plot(coronaUS, lastday, N, savefigs=saveplots, days_before=dbf)
N+=1
cvp.rate_US_plot(coronaUS, lastday, N, savefigs=saveplots, days_before=dbf)
# N+=1
# cvp.active_hosp_US_plot(coronaUS_ctp, lastday, N, savefigs=saveplots, days_before=dbf,
#                         capacity=False)
#N+=1
#cvp.tests_CFR_US_plot(coronaUS, lastday, N, savefigs=saveplots, days_before=dbf)
# N+=1
# cvp.hosp_cap_deaths_US_plot(coronaUS_ctp, lastday, N, savefigs=saveplots, days_before=dbf)


# # custom analysis
# co = coronaUS_ctp['Colorado']
# poslast = co["positive"][-1] - co["positive"][-2]
# teslast = co["totalTestResults"][-1] - co["totalTestResults"][-2]
# percpos = poslast/teslast
# print("current fraction of positive tests in CO = %g" % percpos)
# posidaily = np.diff(co["positive"])
# testdaily = np.diff(co["totalTestResults"])
# posfracdaily = posidaily/testdaily                    
# newconfirmed = np.diff(co["positive"])/co["population"]

# N+=1
# figure(N)
# clf()
# days = co["days"]
# plot(days[1:], newconfirmed*1e5, "--", lw=1, label="new confirmed")
# plot(days[1:], posfracdaily*100, lw=1, label="new positive test percent")
# plot(days, co["hspcur_pc"]*10, lw=2, label="hopitalizations")
# legend(loc='best')
# axis(xmin=-dbf)#,ymax=30)
# xlabel("days before %s" % lastday.date())
# ylabel("percent or per 100,000")
# title("Colorado")
# if saveplots: savefig("plots/testing_hosp_CO.pdf", bbox_inches='tight')
