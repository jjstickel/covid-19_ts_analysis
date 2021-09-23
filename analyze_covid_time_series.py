"""
Read in COVID-19 timeseries data and perform some basic analysis.  See
header of `covid19ts.py` or `analyze_covid_time_series_notebook.md` for list of
data sources

"""

# Jonathan Stickel, 2020, 2021

# on 3/23/20, J-H switched to a new set of csv files
# on 7/15/20, switched to using Covid Tracking Project for US data
# after March 2021, Covid Tracking Project stopped updating their data, so no longer using it

# TODO:
# - Process countries with multiple entries. Will need to make sure that the
#   sum provides the correct result
# - switch to top-bottom plotting for most plots?
# - (long term) switch from dict to class

import numpy as np

from covid19ts import covid19_global, covid19_US, covid19_can
import covid_plots as cvp

from matplotlib.pyplot import *
ion()

## User input -- put up to 7 countries of interest in this list. Must be the
## same name used in the JH global files, and (at the moment), it must be a
## single entry in the file (e.g., China has multiple entries and will cause an
## Exception)
#countries = ["US", "Italy", "Spain", "Sweden", "Brazil"]
countries = ["US", "Italy", "Spain", "Germany", "Sweden", "Brazil", "India"]
#countries = ["US", "Sweden", "Denmark", "Norway"]

# US locations, up to 7 (`US`, States, and counties in the form `[name] County, [ST]`
# where ST is the state code)
US_locs = ["US", "Colorado", "Idaho", "New York", "Florida", "Arizona", "Alabama"]
#US_locs = ["Colorado", "New York", "Florida", "Wisconsin", "North Dakota", "South Dakota"]
#US_locs = ["US", "Colorado", "Jefferson County, CO", "Douglas County, CO", "Denver County, CO", "Boulder County, CO"]
#US_locs = ["US", "Colorado", "Jefferson County, CO", "Larimer County, CO",
#           "Alaska", "Anchorage Municipality, AK"]

#dbf = None
dbf = 550
nsub = 7 # subsample every `nsub` points
if (nsub > 14):
    raise Warning("Subsampling period of %g is too large (>14) for estimating active cases" % nsub)

saveplots = False

#JHCSSEpath = "../JH_COVID-19/csse_covid_19_data/csse_covid_19_time_series/" # github 
JHCSSEpath = "../JH_COVID-19/" # direct download

lmbd = 5e-5
mult = 1e4
#mult = 1e6
corona = covid19_global(countries, websource=False, JHCSSEpath=JHCSSEpath, lmbd=lmbd,
                        mult=mult, dbf=dbf, nsub=nsub)

#mult = corona["mult"]
#critlow = corona["critlow"]
#nctry = len(countries)
nUSloc = len(US_locs)

dates = corona["dates"]
lastday = dates[-1]

# ## not using this data set for the time being -- will be interesting to check
# ## for differences with the other US data set; cannot include `US` as one of the locations
# coronaUS = covid19_US(US_locs[1:], websource=False, JHCSSEpath=JHCSSEpath,
#                       lmbd=lmbd, mult=mult, dbf=dbf, nsub=nsub)
# if not np.alltrue(corona["dates"] == coronaUS["dates"]):
#     raise ValueError("the dates from the global and US files do not match")

# analyze US data from the Covid Act Now database
coronaUS_can = covid19_can(US_locs, lastday, websource=False, sourcepath="../covidactnow/",
                           lmbd=lmbd, mult=mult, dbf=dbf, nsub=nsub)

# estimate "active" cases; since data for recovered cases is so unreliable,
# just this estimate is used
# estimate of recovery time from:
# - https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(20)30243-7/fulltext
# - https://towardsdatascience.com/visual-notes-from-singapores-first-100-fully-recovered-covid-19-patients-aad7f2e1d0a0
rectime = 14 # days, 1 day = 1 original data point
recsamp = int(14/nsub)
for country in countries:
    ctryd = corona[country]
    ctryd["acvest_pc"] = ctryd["cnf_pc_h"].copy()
    ctryd["acvest_pc"][:recsamp] = np.nan
    ctryd["acvest_pc"][recsamp:] = ctryd["acvest_pc"][recsamp:] - ctryd["cnf_pc"][:-recsamp]
    # some anomlies in the case data can lead to negative active cases for this
    # crude estimation
    ctryd["acvest_pc"][ctryd["acvest_pc"] < 0] = 0
for loc in coronaUS_can["locs"]:
    locd = coronaUS_can[loc]
    locd["acvest_pc"] = locd["cnf_pc_h"].copy()
    locd["acvest_pc"][:recsamp] = np.nan
    locd["acvest_pc"][recsamp:] = locd["acvest_pc"][recsamp:] - locd["cnf_pc"][:-recsamp]
    locd["acvest_pc"][locd["acvest_pc"] < 0] = 0

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
cvp.per_capita_US_plot(coronaUS_can, lastday, N, savefigs=saveplots, days_before=dbf)
N+=1
cvp.rate_US_plot(coronaUS_can, lastday, N, savefigs=saveplots, days_before=dbf)
N+=1
cvp.active_hosp_US_plot(coronaUS_can, lastday, N, savefigs=saveplots, days_before=dbf)
N+=1
cvp.icu_US_plot(coronaUS_can, lastday, N, savefigs=saveplots, days_before=dbf)
N+=1
cvp.tests_vacc_US_plot(coronaUS_can, lastday, N, savefigs=saveplots, days_before=dbf)


#N+=1
# this plot isn't very interesting anymore -- make a vaccination plot instead of tests
#cvp.tests_CFR_US_plot(coronaUS_can, lastday, N, savefigs=saveplots, days_before=dbf)
# N+=1
# cvp.hosp_cap_deaths_US_plot(coronaUS_can, lastday, N, savefigs=saveplots, days_before=dbf)


# # custom analysis
# co = coronaUS_can['Colorado']
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
