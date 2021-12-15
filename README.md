# COVID-19 time-series analysis

This purpose of this repository is to make available a set of time-series
analyses of COVID-19 data. Data sources are:

- https://github.com/CSSEGISandData/COVID-19 (COVID-19 data)
- https://covidactnow.org (COVID-19 US data)
- https://data.worldbank.org/indicator/sp.pop.totl (world population data)
- https://www.census.gov/data/tables/time-series/demo/popest/2010s-state-total.html (US population)

The python script `analyze_covid_time_series.py` is setup to produce an example
set of results. 

A jupyter notebook is implemented via
`analyze_covid_time_series_notebook.md` (jupytext required). A static snapshot
including figures is available at:

https://colab.research.google.com/drive/1rH7V7MzUNIXpS1R65beDz-klu-y1xJah?usp=sharing

The python notebook can be run "live" by clicking on this binder link (it may
take awhile to load, be patient):

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/jjstickel/covid-19_ts_analysis.git/master?filepath=analyze_covid_time_series_notebook.md)
