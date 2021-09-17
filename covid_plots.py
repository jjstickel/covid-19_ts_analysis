"""
module of plotting function-routines for covid data
"""

## TODO:
# - plot icu per-cap vs icu frac
# - plot vaccinations and tests together

from decimal import Decimal
from datetime import datetime
from matplotlib.pyplot import *
import matplotlib.dates as mdates


### plotting ###
ion()  ### does initiating ion() work here? seems too...
fw = 6
fh = 4
clr = ['C%g' % i for i in range(10)]
sbl = ["o", "s", "v", "d", "x", "^", "*"]
mew = 0.5
ms = 4
lw = 1.5
# def critlow_readable(corona):
#     # get inverse human readable form for t=0 criterium
#     critlow = corona["critlow"]
#     clinv_tup = Decimal("%.1g"%(1/critlow)).as_tuple()
#     clinv_dig = clinv_tup.digits[0]
#     clinv_exp = clinv_tup.exponent
#     corona["clinv_dig"] = clinv_dig
#     corona["clinv_exp"] = clinv_exp
#     return

#days_before = 75 # days before present day to show on many plots

# us yearly deaths
us_tot_d = 87. # per 10,000, average of 2017-2019
#IFR = 0.00725153*1e4 # average IFR, see covid-ifr.py
IFR = 0.0217777*1e4 # average IFR, see covid-ifr.py, updated since 3/21

#### global plotting functions #####

def total_global_plot(corona, N=1):
    # total cases
    countries = corona["locs"]
    nctry = len(countries)
    dates = corona["dates"]
    mult = corona["mult"]
    figure(N, figsize=(2*fw,fh))
    clf()
    subplot(121)
    for i in range(nctry):
        ctryd = corona[countries[i]]
        population = ctryd['population']
        plot(dates, ctryd['cnf'], sbl[i]+clr[i], mfc='none', mew=mew, ms=ms, label=ctryd["name"])
        #plot(dates, ctryd['cnf_expfit']*population/mult, '--'+clr[i], lw=lw)
        plot(dates, ctryd['cnf_pc_h']*population/mult, '-'+clr[i])
    ax = gca()
    ax.get_xaxis().set_major_formatter(mdates.DateFormatter("%m/%d"))
    setp(ax.get_xticklabels(), rotation=30, ha="right")
    cnf_max = max([corona[country]['cnf'].max() for country in countries])
    axis(xmin = datetime(year=2020, month=2, day = 15), ymin=0-cnf_max*0.1,
         ymax = cnf_max*1.1)
    xlabel("date [m/d]")
    ylabel("confirmed")
    #legend(loc='best')
    title("confirmed cases")
    subplot(122)
    for i in range(nctry):
        ctryd = corona[countries[i]]
        population = ctryd['population']
        plot(dates, ctryd['dth'], sbl[i]+clr[i], mfc='none', mew=mew, ms=ms, label=ctryd["name"])
        plot(dates, ctryd['dth_pc_h']*population/mult, '-'+clr[i])
    ax = gca()
    ax.get_xaxis().set_major_formatter(mdates.DateFormatter("%m/%d"))
    setp(ax.get_xticklabels(), rotation=30, ha="right")
    dth_max = max([corona[country]['dth'].max() for country in countries])
    axis(xmin = datetime(year=2020, month=2, day = 15), ymin=0-dth_max*0.1,
         ymax = dth_max*1.1)
    xlabel("date [m/d]")
    ylabel("deaths")
    legend(loc='best')
    title("deaths");
    return

def per_capita_global_plot(corona, lastday, N=1, savefigs=False, days_before=None):
    # per capita cases
    countries = corona["locs"]
    nctry = len(countries)
    dates = corona["dates"]
    mult = corona["mult"]
    days = corona["days"] #- ctryd["days"][-1]
    #    clinv_dig = corona["clinv_dig"]
    #    clinv_exp = corona["clinv_exp"]
    if days_before is not None:
        days_before = -days_before

    figure(N, figsize=(2*fw,fh))
    clf()
    subplot(121)
    for i in range(nctry):
        ctryd = corona[countries[i]]
        plot(days, ctryd["cnf_pc"], sbl[i]+clr[i], mfc='none', mew=mew, ms=ms)
        #plot(days, ctryd["cnf_expfit"], '--'+clr[i], lw=lw)
        plot(days, ctryd["cnf_pc_h"], '-'+clr[i], label=ctryd["name"])
    scaled_max = max([corona[country]['cnf_pc'].max() for country in countries])
    #axis(xmin=-5, ymin=0-scaled_max*0.1, ymax=scaled_max*1.1)
    axis(xmin = days_before)
    #xlabel('days since %i confirmed per $10^%i$' % (clinv_dig, clinv_exp))
    xlabel("days before %s" % dates[-1].date())
    ylabel("confirmed per $10^%i$" % np.log10(mult))
    #legend(loc='best')
    title("confirmed per capita")
    subplot(122)
    for i in range(nctry):
        ctryd = corona[countries[i]]
        plot(days, ctryd["dth_pc"], sbl[i]+clr[i], mfc='none', mew=mew, ms=ms)
        plot(days, ctryd["dth_pc_h"], '-'+clr[i], label=ctryd["name"])
#    plot(days, 50*np.ones(days.shape), '--k')#, label="herd immunity?")
    #scaled_max = max([corona[country]['dth_pc'].max() for country in countries])
    #axis(xmin=-5, ymin=0-scaled_max*0.1, ymax=scaled_max*1.1)
    axis(xmin = days_before)
#    annotate("herd immunity?", (0.5, 0.9), xycoords="axes fraction")
    #xlabel('days since %i confirmed per $10^%i$' % (clinv_dig, clinv_exp))
    xlabel("days before %s" % dates[-1].date())
    ylabel("deaths per $10^%i$" % np.log10(mult))
    #legend(loc='best')
    legend(loc='upper left')
    title("deaths per capita");
    if savefigs:  savefig("plots/per_capita_global.pdf", bbox_inches="tight")
    return

def rate_global_plot(corona, lastday, N=1, savefigs=False, days_before=None):
    # rate confirmed per capita
    countries = corona["locs"]
    nctry = len(countries)
    dates = corona["dates"]
    mult = corona["mult"]
    days = corona["days"]
#    clinv_dig = corona["clinv_dig"]
#    clinv_exp = corona["clinv_exp"]
    if days_before is not None:
        days_before = -days_before

    figure(N, figsize=(2*fw,fh))
    clf()
    subplot(121)
    for i in range(nctry):
        ctryd = corona[countries[i]]
        plot(days, ctryd["cnf_rate"], "-"+sbl[i]+clr[i], mfc='none', mew=mew, ms=ms,
             label=ctryd["name"])
    axis(xmin = days_before)
    xlabel("days before %s" % dates[-1].date())
    ylabel("rate [confirmed per $10^%i$ / day]" % np.log10(mult))
    #legend(loc="best")
    title("growth rate of confirmed cases")
    subplot(122)
    for i in range(nctry):
        ctryd = corona[countries[i]]
        plot(days, ctryd["dth_rate"], "-"+sbl[i]+clr[i], mfc='none', mew=mew, ms=ms,
             label=ctryd["name"])
    axis(xmin = days_before)#, ymax = 1.5e3/mult)
    xlabel("days before %s" % dates[-1].date())
    ylabel("rate [deaths per $10^%i$ / day]" % np.log10(mult))
    legend(loc="upper left")
    title("growth rate of deaths");
    if savefigs:  savefig("plots/rate_global.pdf", bbox_inches="tight")
    return

def active_CFR_global_plot(corona, lastday, N=1, savefigs=False, days_before=None):
    # active and CFR
    countries = corona["locs"]
    nctry = len(countries)
    dates = corona["dates"]
    days = corona["days"]
    mult = corona["mult"]
#    clinv_dig = corona["clinv_dig"]
#    clinv_exp = corona["clinv_exp"]
    if days_before is not None:
        days_before = -days_before

    figure(N, figsize=(2*fw,fh))
    clf()
    subplot(121)
    for i in range(nctry):
        ctryd = corona[countries[i]]
        #plot(days, ctryd["acv_pc"], sbl[i]+clr[i], mfc='none', mew=mew, ms=ms, label=ctryd["name"])
        plot(days, ctryd["acvest_pc"], "-"+sbl[i]+clr[i], mfc='none', mew=mew, ms=ms,
             label=ctryd["name"])
        #plot(days, ctryd["acvest_pc"], "-"+clr[i], lw=2)
    #scaled_max = max([corona[country]['cnf_pc'].max() for country in countries])
    axis(xmin=days_before)#, ymin=0-scaled_max*0.1, ymax=scaled_max*1.1)
    xlabel("days before %s" % dates[-1].date())
    ylabel("active per $10^%i$" % np.log10(mult))
    #legend(loc='best')
    title("active cases")
    subplot(122)
    for i in range(nctry):
        ctryd = corona[countries[i]]
        plot(days, ctryd["dth_pc_h"]/ctryd["cnf_pc_h"]*100, "-"+sbl[i]+clr[i], mfc='none',
             mew=mew, ms=ms, label=ctryd["name"])
    axis(xmin = days_before, ymin=0-0.05*15, ymax=15)
    xlabel("days before %s" % dates[-1].date())
    ylabel("case-fatality ratio [%]")
    legend(loc="best")
    title("case fatality ratio (CFR)");
    if savefigs:  savefig("plots/active_CFR_global.pdf", bbox_inches="tight")
    return

def exp_fit_confirmed_plot(corona, N=1):
    countries = corona["locs"]
    nctry = len(countries)
    dates = corona["dates"]
    mult = corona["mult"]
#    clinv_dig = corona["clinv_dig"]
#    clinv_exp = corona["clinv_exp"]
    figure(N, figsize=(2*fw,fh))
    clf()
    # exponential fit illustration
    subplot(121)
    for i in range(nctry):
        ctryd = corona[countries[i]]
        days = ctryd["days"]
        plot(days, ctryd["cnf_pc"], sbl[i]+clr[i], mfc='none', mew=mew, ms=ms, label=ctryd["name"])
        plot(days, ctryd["cnf_expfit"], '--'+clr[i], lw=lw)
        plot(days, ctryd["cnf_pc_h"], '-'+clr[i])
    scaled_max = max([corona[country]['cnf_pc'].max() for country in countries])
    axis(xmin=-5, ymin=0-scaled_max*0.1, ymax=scaled_max*1.1)
#    xlabel('days since %i confirmed per $10^%i$' % (clinv_dig, clinv_exp))
    ylabel("confirmed per $10^%i$" % np.log10(mult))
    #legend(loc='best')
    title("confirmed per capita")
    # log-scale confirmed per capita
    subplot(122)
    for i in range(nctry):
        ctryd = corona[countries[i]]
        semilogy(ctryd["days"], ctryd["cnf_pc"], sbl[i]+clr[i], mfc='none', mew=mew, ms=ms,
                 label=ctryd["name"])
        semilogy(ctryd["days"], ctryd["cnf_expfit"], '--'+clr[i], lw=lw)
        semilogy(ctryd["days"], ctryd["cnf_pc_h"], '-'+clr[i])
    axis(xmin=-5, ymin=1e-3)
#    xlabel('days since %i confirmed per $10^%i$' % (clinv_dig, clinv_exp))
    ylabel("confirmed per $10^%i$" % np.log10(mult))
    legend(loc='best')
    title("per capita confirmed, log scale");
    return

def confirmed_deaths_simul_global_plot(corona, lastday, N=1, days_before=None):
    # confirmed per-capita cases and deaths on one graph
    countries = corona["locs"]
    nctry = len(countries)
    dates = corona["dates"]
    mult = corona["mult"]
    days = corona["days"]
#    clinv_dig = corona["clinv_dig"]
#    clinv_exp = corona["clinv_exp"]
    if days_before is not None:
        days_before = -days_before

    figure(N, figsize=(fw,fh))
    clf()
    ax1 = subplot(111)
    ax2 = ax1.twinx()
    for i in range(nctry):
        ctryd = corona[countries[i]]
        ax1.plot(days, ctryd["cnf_pc"], "-"+sbl[i]+clr[i], mfc='none', mew=mew, ms=ms,
                 label=ctryd["name"])
        ax2.plot(days, ctryd["dth_pc"], ":"+sbl[i]+clr[i], mfc='none', mew=mew, ms=ms)
    cnf_max = max([corona[country]['cnf_pc'].max() for country in countries])
    dth_max = max([corona[country]['dth_pc'].max() for country in countries])
    ax1.axis(xmin=days_before, ymin=0-cnf_max*0.1, ymax=cnf_max*1.1)
    ax2.axis(xmin=days_before, ymin=0-dth_max*0.1, ymax=dth_max*1.1)
    ax1.set_xlabel("days before %s" % dates[-1].date())
    #ax1.set_xlabel('days since %i deaths per $10^%i$' % (clinv_dig, clinv_exp))
    ax1.set_ylabel("confirmed per $10^%i$" % np.log10(mult))
    ax2.set_ylabel("deaths per $10^%i$" % np.log10(mult))
    handles, labels = ax1.get_legend_handles_labels()
    legend(handles, labels, loc='best')
    title("per capita confirmed and deaths on one graph");
    return


##### US plotting functions ####
k = 0 # color shift so that US location colors are different from global colors

def per_capita_US_plot(corona, lastday, N=1, savefigs=False, days_before=None):
    # per capita cases
    locs = corona["locs"]
    nloc = len(locs)
    mult = corona["mult"]
    daysmaster = "days" in corona.keys()
    if daysmaster:
        days = corona["days"]
    if days_before is not None:
        days_before = -days_before
    figure(N, figsize=(2*fw,fh))
    clf()
    subplot(121)
    for i in range(nloc):
        locd = corona[locs[i]]
        if ~daysmaster:
            days = locd["days"]
        plot(days, locd["cnf_pc"], sbl[i]+clr[i], mfc='none', mew=mew, ms=ms)
        #plot(days, locd["positive"]/locd["population"]*100, sbl[i]+clr[i], mfc='none', mew=mew, ms=ms, label=locd["name"])
        plot(days, locd["cnf_pc_h"], '-'+clr[i], label=locd["name"])
    maxvals = [np.nanmax(corona[loc]['cnf_pc']) for loc in locs]
    scaled_max = max(maxvals)
    #axis(xmin=days_before, ymin=0-scaled_max*0.1, ymax=scaled_max*1.1)
    axis(xmin=days_before)
    xlabel("days before %s" % lastday.date())
    ylabel("confirmed per $10^%i$" % np.log10(mult))
    #ylabel("confirmed [%]")
    #legend(loc='best')
    title("confirmed per capita")
    # per capita deaths
    subplot(122)
    for i in range(nloc):
        locd = corona[locs[i]]
        if ~daysmaster:
            days = locd["days"]
        plot(days, locd["dth_pc"], sbl[i]+clr[i], mfc='none', mew=mew, ms=ms)
        plot(days, locd["dth_pc_h"], '-'+clr[i], label=locd["name"])
    maxvals = [np.nanmax(corona[loc]['dth_pc']) for loc in locs]
    scaled_max = max(maxvals)
    axis(xmin=days_before)
    xlabel("days before %s" % lastday.date())
    ylabel("deaths per $10^%i$" % np.log10(mult))
    #ylabel("deaths [%]")
    legend(loc='best')
    title("deaths per capita");
    if savefigs:  savefig("plots/per_capita_US.pdf", bbox_inches="tight")
    return

def rate_US_plot(corona, lastday, N=1, savefigs=False, days_before=None):
    # rate confirmed per capita
    locs = corona["locs"]
    nloc = len(locs)
    mult = corona["mult"]
    daysmaster = "days" in corona.keys()
    if daysmaster:
        days = corona["days"]
    if days_before is not None:
        days_before = -days_before

    figure(N, figsize=(2*fw,fh))
    clf()
    subplot(121)
    for i in range(nloc):
        locd = corona[locs[i]]
        if ~daysmaster:
            days = locd["days"] 
        plot(days, locd["cnf_rate"], "-"+sbl[i]+clr[i], mfc='none', mew=mew, ms=ms,
             label=locd["name"])
    axis(xmin = days_before)
    xlabel("days before %s" % lastday.date())
    ylabel("rate [confirmed per $10^%i$ / day]" % np.log10(mult))
    #ylabel("rate of increase in cases [% pop. per day]")
    #legend(loc="best")
    title("growth rate of confirmed cases")
    # rate deaths per capita
    subplot(122)
    for i in range(nloc):
        locd = corona[locs[i]]
        if ~daysmaster:
            days = locd["days"] 
        plot(days, locd["dth_rate"], "-"+sbl[i]+clr[i], mfc='none', mew=mew, ms=ms,
             label=locd["name"])
    #axis(xmin=days_before, ymax = 1.5e3/mult)
    axis(xmin=days_before)
    xlabel("days before %s" % lastday.date())
    ylabel("rate [deaths per $10^%i$ / day]" % np.log10(mult))
    #ylabel("rate of increase in deaths [% pop. per day]")
    legend(loc='best')
    title("growth rate of deaths");
    if savefigs:  savefig("plots/rate_US.pdf", bbox_inches="tight")
    return


def active_hosp_US_plot(corona, lastday, N=1, savefigs=False, days_before=None):
    # active cases per capita
    locs = corona["locs"]
    nloc = len(locs)
    mult = corona["mult"]
    figure(N, figsize=(2*fw,fh))
    clf()
    subplot(121)
    for i in range(nloc):
        locd = corona[locs[i]]
        days = locd["days"] 
        plot(days, locd["acvest_pc"], "-"+sbl[i]+clr[i], lw=lw, mfc='none', mew=mew, ms=ms,
             label=locd["name"])
    if days_before is not None:
        days_before = -days_before
    axis(xmin=days_before)
    xlabel("days before %s" % lastday.date())  
    ylabel("active cases per $10^%i$" % np.log10(mult))
    #ylabel("active cases [%]")
    #legend(loc='best')
    title("active cases")
    # hospitalizations per capita
    subplot(122)
    for i in range(nloc):
        locd = corona[locs[i]]
        days = locd["days"] 
        plot(days, locd["hosp_covid_pc"], "-"+sbl[i]+clr[i], lw=lw, mfc='none', mew=mew,
             ms=ms, label=locd["name"])
#        if capacity:  # no longer using this, switching to ICU data with covid act now data
#            plot(days, locd["hsp_cap_pc"]*np.ones(days.size), "--"+clr[i], lw=lw)
    axis(xmin=days_before)
    xlabel("days before %s" % lastday.date())
    ylabel("hospitalizations per $10^%i$" % np.log10(mult))
    #ylabel("hospitalizations [%]")
    legend(loc='upper left')
    title("current COVID-19 hospitalizations")
    if savefigs:  savefig("plots/active_hosp_US.pdf", bbox_inches="tight")
    return


### this is a bit messy; do a 2x2 subplot? JJS 8/30/21
def hosp_icu_US_plot(corona, lastday, N=1, savefigs=False, days_before=None):
    # COVID-19 hospitalizations and ICU per capita
    locs = corona["locs"]
    nloc = len(locs)
    mult = corona["mult"]
    figure(N, figsize=(2*fw,fh))
    clf()
    subplot(121)
    for i in range(nloc):
        locd = corona[locs[i]]
        days = locd["days"] 
        plot(days, locd["hosp_covid_pc"], "-"+sbl[i]+clr[i], lw=lw, mfc='none', mew=mew,
             ms=ms, label=locd["name"])
#        plot(days, locd["icu_covid_pc"], "--"+sbl[i]+clr[i], lw=lw, mfc='none', mew=mew,
#             ms=ms)
    if days_before is not None:
        days_before = -days_before
    axis(xmin=days_before)
    xlabel("days before %s" % lastday.date())  
    ylabel("hospital and icu beds per $10^%i$" % np.log10(mult))
    legend(loc='best')
    title("Current COVID-19 hospitalizations")
    # ICU bed usage, percent of total
    subplot(122)
    for i in range(nloc):
        locd = corona[locs[i]]
        days = locd["days"] 
        plot(days, locd["icu_covid_frac"]*100, "-"+sbl[i]+clr[i], lw=lw, mfc='none', mew=mew,
             ms=ms, label=locd["name"])
        plot(days, locd["icu_total_frac"]*100, "--"+sbl[i]+clr[i], lw=lw, mfc='none',
             mew=mew, ms=ms)
    axis(xmin=days_before, ymin=0, ymax=100)
    xlabel("days before %s" % lastday.date())
    ylabel("ICU bed usage (%)")
    annotate("Total ICU", (0.5, 0.85), xycoords="axes fraction", ha='center')
    annotate("COVID-19 ICU", (0.5, 0.25), xycoords="axes fraction", ha='center')
    #legend(loc='best')
    title("ICU bed usage")
    if savefigs:  savefig("plots/hosp_icu_US.pdf", bbox_inches="tight")
    return


# vaccinations; plot against deaths/hosp/icu?
def icu_vacc_US_plot(corona, lastday, N=1, savefigs=False, days_before=None):
    # tests fraction
    locs = corona["locs"]
    nloc = len(locs)
    mult = corona["mult"]
    figure(N, figsize=(2*fw,fh))
    clf()
    subplot(121)
    ymax = 100
    for i in range(nloc):
        locd = corona[locs[i]]
        days = locd["days"] 
        plot(days, locd["icu_covid_frac"]*100, "-"+sbl[i]+clr[i], lw=lw, mfc='none', mew=mew,
             ms=ms, label=locd["name"])
        plot(days, locd["icu_total_frac"]*100, "--"+sbl[i]+clr[i], lw=lw, mfc='none',
             mew=mew, ms=ms)
        # sometimes total icu usage is over 100%, e.g., Alabama
        ymax = max(ymax, np.nanmax(locd["icu_total_frac"])*100)
    if days_before is not None:
        days_before = -days_before
    ymax = min(120, ymax) # to counter some wild data, e.g. Anchorage Municipality, AK
    axis(xmin=days_before, ymin=0, ymax=ymax)
    xlabel("days before %s" % lastday.date())
    ylabel("ICU bed usage (%)")
    annotate("Total ICU", (0.5, 0.85), xycoords="axes fraction", ha='center')
    annotate("COVID-19 ICU", (0.5, 0.25), xycoords="axes fraction", ha='center')
    #legend(loc='best')
    title("ICU bed usage")
    subplot(122)
    for i in range(nloc):
        locd = corona[locs[i]]
        days = locd["days"] 
#        plot(days, locd["vacc_init_frac"]*100, "--"+sbl[i]+clr[i], lw=lw, mew=mew, ms=ms,
#             mfc='none')
        plot(days, locd["vacc_full_frac"]*100, "-"+sbl[i]+clr[i], lw=lw, mew=mew, ms=ms,
             mfc='none', label=locd["name"])
    #axis(xmin=days_before, ymin=0, ymax=100)
    axis(xmin=days_before)
    xlabel("days before %s" % lastday.date())
    ylabel("vaccinations, % of population")
    legend(loc='best')
    title("Completed vaccinations")
    if savefigs:  savefig("plots/icu_vacc_US.pdf", bbox_inches="tight")
    return


def deaths_persp_US_plot(corona, lastday, N=1, savefigs=False, days_before=None):
    if days_before is not None:
        days_before = -days_before
    locs = corona["locs"]
    nloc = len(locs)
    mult = corona["mult"]
    figure(N, figsize=(fw,fh))
    clf()
    for i in range(nloc):
        locd = corona[locs[i]]
        days = locd["days"] 
        plot(days, locd["dth_pc"], sbl[i]+clr[i], mfc='none', mew=mew, ms=ms)
        plot(days, locd["dth_pc_h"], '-'+clr[i], label=locd["name"])
    maxvals = [np.nanmax(corona[loc]['dth_pc']) for loc in locs]
    scaled_max = max(maxvals)
    plot(days, us_tot_d*np.ones(days.shape), '--k')#, label="total deaths 2018")
    annotate("total deaths 2018", (0.5, 0.42), xycoords="axes fraction")
    #plot(days, 0.75*IFR*np.ones(days.shape), ':k', lw=2)#, label="herd immunity?")
    fill_between(days, 0.5*IFR*np.ones(days.shape), 1.0*IFR*np.ones(days.shape),
                 facecolor="black", alpha=0.2)
    annotate("Full endemic penetration?\n (no vaccine)", (0.6, 0.7),
             xycoords="axes fraction", ha="center")
    axis(xmin=days_before, xmax=0)
    xlabel("days before %s" % lastday.date())
    ylabel("deaths per $10^%i$" % np.log10(mult))
    #ylabel("deaths [%]")
    legend(loc='upper left')
    title("deaths per capita");
    if savefigs:  savefig("plots/deaths_perspective_US.pdf", bbox_inches="tight")
    return


#### below plots are no longer interesting and are not maintained, JJS 8/30/21 ####

def hosp_cap_deaths_US_plot(corona, lastday, N=1, savefigs=False, days_before=None):
    locs = corona["locs"]
    nloc = len(locs)
    mult = corona["mult"]
    figure(N, figsize=(2*fw,fh))
    clf()
    # hospitalizations per capita with capacity -- FIXME:  switch to ICU plot
    subplot(121)
    for i in range(nloc):
        locd = corona[locs[i]]
        days = locd["days"] 
        plot(days, locd["hosp_covid_pc"], "-"+sbl[i]+clr[i], lw=lw, mfc='none', mew=mew,
             ms=ms, label=locd["name"])
        #plot(days, locd["hsp_cap_pc"]*np.ones(days.size), "--"+clr[i], lw=lw)
    if days_before is not None:
        days_before = -days_before
    axis(xmin=days_before)
    xlabel("days before %s" % lastday.date())
    ylabel("hospitalizations per $10^%i$" % np.log10(mult))
    #annotate("dashed lines = total capacity", [0.25,0.5], xycoords="axes fraction")
    #ylabel("hospitalizations [%]")
    #legend(loc='best')
    title("current hospitalizations")#\n(dashed lines = total capacity)")
    # deaths per capita with total US deaths
    subplot(122)
    for i in range(nloc):
        locd = corona[locs[i]]
        days = locd["days"] 
        plot(days, locd["dth_pc"], sbl[i]+clr[i], mfc='none', mew=mew, ms=ms)
        plot(days, locd["dth_pc_h"], '-'+clr[i], label=locd["name"])
    maxvals = [np.nanmax(corona[loc]['dth_pc']) for loc in locs]
    scaled_max = max(maxvals)
    plot(days, us_tot_d*np.ones(days.shape), '--k')#, label="total deaths 2018")
    annotate("total deaths 2018", (0.5, 0.9), xycoords="axes fraction")
    #plot(days, 0.75*IFR*np.ones(days.shape), ':k', lw=2)#, label="herd immunity?")
    fill_between(days, 0.5*IFR*np.ones(days.shape), 1.0*IFR*np.ones(days.shape),
                 facecolor="black", alpha=0.2)
    annotate("herd immunity? (no vaccine)", (0.4, 0.6), xycoords="axes fraction")
    axis(xmin=days_before, xmax=0)
    xlabel("days before %s" % lastday.date())
    ylabel("deaths per $10^%i$" % np.log10(mult))
    #ylabel("deaths [%]")
    legend(loc='upper left')
    title("deaths per capita");
    if savefigs:  savefig("plots/hosp_cap_deaths_US.pdf", bbox_inches="tight")
    return


def tests_CFR_US_plot(corona, lastday, N=1, savefigs=False, days_before=None):
    # tests fraction
    locs = corona["locs"]
    nloc = len(locs)
    mult = corona["mult"]
    figure(N, figsize=(2*fw,fh))
    clf()
    subplot(121)
    for i in range(nloc):
        locd = corona[locs[i]]
        days = locd["days"] 
        plot(days, locd["test_frac"]*100, "-"+sbl[i]+clr[i], lw=lw, mew=mew, ms=ms, mfc='none',
             label=locd["name"])
    if days_before is not None:
        days_before = -days_before
    axis(xmin=days_before)
    xlabel("days before %s" % lastday.date())  
    ylabel("tests, % of population")
    #legend(loc='best')
    title("total tests")
    # hospitalizations per capita
    subplot(122)
    for i in range(nloc):
        locd = corona[locs[i]]
        days = locd["days"] 
        plot(days, locd["cfr"]*100, "-"+sbl[i]+clr[i], lw=lw, mew=mew, ms=ms, mfc='none',
             label=locd["name"])
    axis(xmin=days_before)
    xlabel("days before %s" % lastday.date())
    ylabel("case fatality ratio [%]")
    legend(loc='best')
    title("case fatality ratio (CFR)")
    if savefigs:  savefig("plots/tests_CFR_US.pdf", bbox_inches="tight")
    return
