---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.11.4
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

# COVID-19 time-series analysis

Jonathan Stickel, 2020-2021

This purpose of this notebook (and repository) is to make available a set of time-series
analyses of COVID-19 data. Data sources are:

- https://github.com/CSSEGISandData/COVID-19 (COVID-19 data)
- https://covidactnow.org (COVID-19 US data)
- https://data.worldbank.org/indicator/sp.pop.totl (world populations)
- https://www.cdc.gov/ (total US deaths)

Scroll down to see the data plots and analysis. 

To run the entire notebook (via Binder link below or on your local computer), click `Cell` from the menu above and then `Run All`. By default, the COVID-19 data is grabbed from web sources directly. See arguments to functions for specifying files in your local path.

Dependencies are:

- `numpy`
- `pandas`
- `matplotlib`
- `cvxopt`
- `datetime`
- `scikit.datasmooth` (can be pip installed)

Also, `covid19ts` and `covid_plots` are part of this repository.

Open this link in a new tab (right-click, "open link in new tab") for a web-based "live" notebook (it may take awhile to load, be patient):
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/jjstickel/covid-19_ts_analysis.git/master?filepath=analyze_covid_time_series_notebook.md)


```python
# import modules
import numpy as np
# locally defined modules
from covid19ts import covid19_global, covid19_can
import covid_plots as cvp
```

# User input

Put up to 7 countries of interest in the `countries` list. Must be the same name used in the J-H global files, and (at the moment), it must be a single entry in the file (e.g., China has multiple entries and will cause an Exception). 

```python
countries = ["US", "Italy", "Spain", "Germany", "Sweden", "Brazil", "Mexico"]
```

US locations, up to 7 (`US`, States, and counties in the form `[name] County, [ST]` where ST is the state code).

```python
US_locs = ["US", "Colorado", "Idaho", "New York", "Florida", "Arizona", "Alabama"]
# example of counties -- analyzed separately below
#US_locs = ["US", "Colorado", "Jefferson County, CO", "Douglas County, CO", "Denver County, CO", "Boulder County, CO"]
```

Read in COVID-19 timeseries data for the locations specified and perform these operations:
- normalize cases to be per capita
- smooth the cases data
- set time-zero for each location and shift elapsed time in days
- determine rates (i.e., the derivative) for cases

```python
# days before today (`dbf`) to analyze, and subsampling by `nsub`; more data takes a little more processing time;
# use `dbf = None` and `nsub = 1` to use all data
dbf = 550 
nsub = 7 # subsample every `nsub` points
if (nsub > 14):
    raise Warning("Subsampling period of %g is too large (>14) for estimating active cases" % nsub)
# global data
lmbd = 5e-5 # smoothing parameter, larger means more smooth
corona = covid19_global(countries, lmbd=lmbd, dbf=dbf, nsub=nsub)
# extract common variables for ease-of-use
mult = corona["mult"]
nctry = len(countries)
dates = corona["dates"]
lastday = dates[-1]
# US data
coronaUS_can = covid19_can(US_locs, lastday, lmbd=lmbd, dbf=dbf, nsub=nsub)
```

*Estimate* "active" cases by presuming all confirmed cases have recovered or died in an aeverage number of days. Estimate of recovery time from:
- https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(20)30243-7/fulltext
- https://towardsdatascience.com/visual-notes-from-singapores-first-100-fully-recovered-covid-19-patients-aad7f2e1d0a0

Unfortunately, recovered data is really poor, and so directly calculating active cases is not useful.

```python
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
```

# Plotting

```python
# plot setup
cvp.rcParams.update({'font.size': 14})
cvp.fw = 8
cvp.fh = 6
cvp.ms = 5
```

# Per capita cases (confirmed and deaths)


## Global

```python
cvp.per_capita_global_plot(corona, lastday, days_before=dbf)
```

The US has a lot more confirmed per-capita cases than many countries. This could be attributed to more testing. Deaths continue to rise with recent uptick presumably due to Delta variant. See rate plots below.


## US

```python
cvp.per_capita_US_plot(coronaUS_can, lastday, days_before=dbf)
```

US local per capita data. 


# Per capita growth rates (confirmed and deaths)


## Global

```python
cvp.rate_global_plot(corona, lastday, days_before=dbf)
```

Growth rate is the derivative of the cases (i.e., instantaneous slope for each day). Rates have gone up and down over time. A flat rate means linear growth. Exponential growth only happened very early in the pandemic and has been cyclical since, despite all the media buzz.


## US

```python
cvp.rate_US_plot(coronaUS_can, lastday, days_before=dbf)
```

US local rate data. 


# Active cases, CFR, hospitalizations, and vaccinations. 


## Global

```python
cvp.active_CFR_global_plot(corona, lastday, days_before=dbf)
```

Have we peaked? A curve of active cases help us answer this. While an initial peak ocurred long ago, therecontinue to be more waves.

The "case fatality ratio", or *CFR*, is an indication of how deadly a disease is. It is only an indication because it is limited by how many actual cases are measured and *confirmed*. Here, we see that the US is doing pretty good compared to other countries. Generally, more testing increases the denominator of the ratio and makes the CFR *lower*. (Note: the CFR is commonly called the case fatality *rate*. The use of the word rate here is technically incorrect---rate refers to something changing over *time*. [More info here](https://ourworldindata.org/coronavirus?fbclid=IwAR3zOvtt7gqkhitoHJ_lXDr3eDeE_JPtfukpOkY94PSaBm_hmrMvWCXWFpg#what-do-we-know-about-the-risk-of-dying-from-covid-19))


## US

```python
cvp.active_hosp_US_plot(coronaUS_can, lastday, days_before=dbf)
```

The Covid Act Now data has hospitalizations by US state. I find it informative to plot hospitalizations next to active cases. Waves of hospitalizations follow closely active cases.

```python
cvp.icu_US_plot(coronaUS_can, lastday, days_before=dbf)
```

ICU bed usage per capita (left) and percent of total beds (right). Other than Alabama, I haven't noticed many states truly running out of ICU beds (at least not observed by these data).

```python
cvp.tests_vacc_US_plot(coronaUS_can, lastday, days_before=dbf)
```

COVID-19 testing (left) and completed vaccinations (right). There is a slight inverse correlation between vaccinations and ICU bed usage (and hospitalizations as shown in the figure above). Overall, there is not a large spread in vaccinations between states compared to the spread in the COVID-19 incidence data (hospitalizations, deaths).

```python
cvp.deaths_persp_US_plot(coronaUS_can, lastday, days_before=dbf)
```

How bad is COVID-19 really? Deaths are plotted with the total yearly US deaths in 2018. Total COVID-19 deaths so far (per capita) are about 20% of yearly US deaths, and yearly deaths are less than 1% of the population in any given year, so COVID-19 mortality is so far about 0.2% of the population. While not trivial, it is nothing close to disasters of the pre-modern era. The Black Plague mortality is estimated to be about 50%, and famines during the middle ages resulted in 10-25% mortality, sometimes for several years in a row (Wikipedia).


# Colorado counties

```python
US_locs = ["US", "Colorado", "Jefferson County, CO", "Douglas County, CO", "Denver County, CO", "Boulder County, CO"]
coronaUS_can = covid19_can(US_locs, lastday, lmbd=lmbd, dbf=dbf, nsub=nsub)
```

```python
cvp.per_capita_US_plot(coronaUS_can, lastday, days_before=dbf)
```

```python
cvp.rate_US_plot(coronaUS_can, lastday, days_before=dbf)
```

```python
cvp.icu_US_plot(coronaUS_can, lastday, days_before=dbf)
```

```python
cvp.tests_vacc_US_plot(coronaUS_can, lastday, days_before=dbf)
```
