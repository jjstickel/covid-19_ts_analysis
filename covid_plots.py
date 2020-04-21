"""
module of plotting function-routines for covid data
"""

from decimal import Decimal
from datetime import datetime
from matplotlib.pyplot import *
import matplotlib.dates as mdates


### plotting ###
ion()  ### does initiating ion() work here? seems too...
fw = 6
fh = 4
clr = ['C%g' % i for i in range(10)]
sbl = ["o", "s", "v", "d", "x"]
savefigs = False

def critlow_readable(corona):
    # get inverse human readable form for t=0 criterium
    critlow = corona["critlow"]
    clinv_tup = Decimal("%.1g"%(1/critlow)).as_tuple()
    clinv_dig = clinv_tup.digits[0]
    clinv_exp = clinv_tup.exponent
    corona["clinv_dig"] = clinv_dig
    corona["clinv_exp"] = clinv_exp
    return


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
        plot(dates, ctryd['cnf'], sbl[i]+clr[i], mfc='none', mew=1.5, label=ctryd["name"])
        #plot(dates, ctryd['cnf_expfit']*population/mult, '--'+clr[i], lw=1.5)
        plot(dates, ctryd['cnf_pc_h']*population/mult, '-'+clr[i])
    ax = gca()
    ax.get_xaxis().set_major_formatter(mdates.DateFormatter("%m/%d"))
    setp(ax.get_xticklabels(), rotation=30, ha="right")
    cnf_max = max([corona[country]['cnf'].max() for country in countries])
    axis(xmin = datetime(year=2020, month=2, day = 15), ymin=0-cnf_max*0.1, ymax = cnf_max*1.1)
    xlabel("date [m/d]")
    ylabel("confirmed")
    #legend(loc='best')
    title("confirmed cases")
    subplot(122)
    for i in range(nctry):
        ctryd = corona[countries[i]]
        population = ctryd['population']
        plot(dates, ctryd['dth'], sbl[i]+clr[i], mfc='none', mew=1.5, label=ctryd["name"])
        plot(dates, ctryd['dth_pc_h']*population/mult, '-'+clr[i])
    ax = gca()
    ax.get_xaxis().set_major_formatter(mdates.DateFormatter("%m/%d"))
    setp(ax.get_xticklabels(), rotation=30, ha="right")
    dth_max = max([corona[country]['dth'].max() for country in countries])
    axis(xmin = datetime(year=2020, month=2, day = 15), ymin=0-dth_max*0.1, ymax = dth_max*1.1)
    xlabel("date [m/d]")
    ylabel("deaths")
    legend(loc='best')
    title("deaths");
    return

def per_capita_global_plot(corona, N=1):
    # per capita cases
    countries = corona["locs"]
    nctry = len(countries)
    dates = corona["dates"]
    mult = corona["mult"]
    clinv_dig = corona["clinv_dig"]
    clinv_exp = corona["clinv_exp"]
    figure(N, figsize=(2*fw,fh))
    clf()
    subplot(121)
    for i in range(nctry):
        ctryd = corona[countries[i]]
        days = ctryd["days"]
        plot(days, ctryd["cnf_pc"], sbl[i]+clr[i], mfc='none', mew=1.5, label=ctryd["name"])
        plot(days, ctryd["cnf_expfit"], '--'+clr[i], lw=1.5)
        plot(days, ctryd["cnf_pc_h"], '-'+clr[i])
    scaled_max = max([corona[country]['cnf_pc'].max() for country in countries])
    axis(xmin=-5, ymin=0-scaled_max*0.1, ymax=scaled_max*1.1)
    xlabel('days since %i confirmed per $10^%i$' % (clinv_dig, clinv_exp))
    ylabel("confirmed per $10^%i$" % np.log10(mult))
    #legend(loc='best')
    title("confirmed per capita")
    subplot(122)
    for i in range(nctry):
        ctryd = corona[countries[i]]
        days = ctryd["days"]
        plot(days, ctryd["dth_pc"], sbl[i]+clr[i], mfc='none', mew=1.5, label=ctryd["name"])
        plot(days, ctryd["dth_pc_h"], '-'+clr[i])
    scaled_max = max([corona[country]['dth_pc'].max() for country in countries])
    axis(xmin=-5, ymin=0-scaled_max*0.1, ymax=scaled_max*1.1)
    xlabel('days since %i confirmed per $10^%i$' % (clinv_dig, clinv_exp))
    ylabel("deaths per $10^%i$" % np.log10(mult))
    legend(loc='best')
    title("deaths per capita");
    return

def rate_global_plot(corona, N=1):
    # rate confirmed per capita
    countries = corona["locs"]
    nctry = len(countries)
    dates = corona["dates"]
    mult = corona["mult"]
    clinv_dig = corona["clinv_dig"]
    clinv_exp = corona["clinv_exp"]
    figure(N, figsize=(2*fw,fh))
    clf()
    subplot(121)
    for i in range(nctry):
        ctryd = corona[countries[i]]
        plot(ctryd["days"], ctryd["cnf_rate"], "-"+sbl[i]+clr[i], mfc='none', mew=1.5,
             label=ctryd["name"])
    axis(xmin = -5)
    xlabel('days since %i confirmed per $10^%i$' % (clinv_dig, clinv_exp))
    ylabel("rate [confirmed per $10^%i$ / day]" % np.log10(mult))
    #legend(loc="best")
    title("growth rate of confirmed cases")
    subplot(122)
    for i in range(nctry):
        ctryd = corona[countries[i]]
        plot(ctryd["days"], ctryd["dth_rate"], "-"+sbl[i]+clr[i], mfc='none', mew=1.5,
             label=ctryd["name"])
    axis(xmin = -5)
    xlabel('days since %i confirmed per $10^%i$' % (clinv_dig, clinv_exp))
    ylabel("rate [deaths per $10^%i$ / day]" % np.log10(mult))
    legend(loc="best")
    title("growth rate of deaths");
    return

def active_CFR_global_plot(corona, N=1):
    # active and CFR
    countries = corona["locs"]
    nctry = len(countries)
    dates = corona["dates"]
    mult = corona["mult"]
    clinv_dig = corona["clinv_dig"]
    clinv_exp = corona["clinv_exp"]
    figure(N, figsize=(2*fw,fh))
    clf()
    subplot(121)
    for i in range(nctry):
        ctryd = corona[countries[i]]
        days = ctryd["days"]
        plot(days, ctryd["acv_pc"], sbl[i]+clr[i], mfc='none', mew=1.5, label=ctryd["name"])
        plot(days, ctryd["acvest_pc"], "-"+clr[i], lw=2)
    #scaled_max = max([corona[country]['cnf_pc'].max() for country in countries])
    axis(xmin=-5)#, ymin=0-scaled_max*0.1, ymax=scaled_max*1.1)
    xlabel('days since %i confirmed per $10^%i$' % (clinv_dig, clinv_exp))
    ylabel("active per $10^%i$" % np.log10(mult))
    #legend(loc='best')
    title("active cases")
    subplot(122)
    for i in range(nctry):
        ctryd = corona[countries[i]]
        plot(ctryd["days"], ctryd["dth_pc_h"]/ctryd["cnf_pc_h"]*100, "-"+sbl[i]+clr[i], mfc='none',
             mew=1.5, label=ctryd["name"])
    axis(xmin = -5, ymax=15)
    xlabel('days since %i confirmed per $10^%i$' % (clinv_dig, clinv_exp))
    ylabel("case-fatality ratio [%]")
    legend(loc="best")
    title("case fatality ratio (CFR)");
    return

def log_confirmed_plot(corona, N=1):
    # log-scale confirmed per capita
    countries = corona["locs"]
    nctry = len(countries)
    dates = corona["dates"]
    mult = corona["mult"]
    clinv_dig = corona["clinv_dig"]
    clinv_exp = corona["clinv_exp"]
    figure(N, figsize=(fw,fh))
    clf()
    for i in range(nctry):
        ctryd = corona[countries[i]]
        semilogy(ctryd["days"], ctryd["cnf_pc"], sbl[i]+clr[i], mfc='none', mew=1.5,
                 label=ctryd["name"])
        semilogy(ctryd["days"], ctryd["cnf_expfit"], '--'+clr[i], lw=1.5)
        semilogy(ctryd["days"], ctryd["cnf_pc_h"], '-'+clr[i])
    axis(xmin=-5, ymin=1e-3)
    xlabel('days since %i confirmed per $10^%i$' % (clinv_dig, clinv_exp))
    ylabel("confirmed per $10^%i$" % np.log10(mult))
    legend(loc='best')
    title("per capita confirmed, log scale");
    return

def confirmed_deaths_simul_global_plot(corona, N=1):
    # confirmed per-capita cases and deaths on one graph
    countries = corona["locs"]
    nctry = len(countries)
    dates = corona["dates"]
    mult = corona["mult"]
    clinv_dig = corona["clinv_dig"]
    clinv_exp = corona["clinv_exp"]
    figure(N, figsize=(fw,fh))
    clf()
    ax1 = subplot(111)
    ax2 = ax1.twinx()
    for i in range(nctry):
        ctryd = corona[countries[i]]
        days = ctryd["days"]
        ax1.plot(days, ctryd["cnf_pc"], "-"+sbl[i]+clr[i], mfc='none', mew=1.5,
                 label=ctryd["name"])
        ax2.plot(days, ctryd["dth_pc"], ":"+sbl[i]+clr[i], mfc='none', mew=1.5, ms=3)
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
    title("per capita confirmed and deaths on one graph");
    return


##### US plotting functions ####
k = 4 # color shift so that US location colors are different from global colors

def per_capita_US_plot(coronaUS, corona, N=1):
    # per capita cases
    US_locs = coronaUS["locs"]
    nUSloc = len(US_locs)
    dates = corona["dates"]
    mult = corona["mult"]
    clinv_dig = corona["clinv_dig"]
    clinv_exp = corona["clinv_exp"]
    USd = corona["US"]
    figure(N, figsize=(2*fw,fh))
    clf()
    subplot(121)
    plot(USd["days"], USd["cnf_pc"], sbl[0]+clr[0], mfc='none', mew=1.5, label=USd["name"])
    plot(USd["days"], USd["cnf_pc_h"], "-"+clr[0], mfc='none', mew=1.5)
    for i in range(nUSloc):
        j = i+1
        locd = coronaUS[US_locs[i]]
        days = locd["days"]
        plot(days, locd["cnf_pc"], sbl[j]+clr[j+k], mfc='none', mew=1.5, label=locd["name"])
        plot(days, locd["cnf_pc_h"], '-'+clr[j+k])
    scaled_max = max([coronaUS[loc]['cnf_pc'].max() for loc in US_locs])
    axis(xmin=-5, ymin=0-scaled_max*0.1, ymax=scaled_max*1.1)
    xlabel('days since %i confirmed per $10^%i$' % (clinv_dig, clinv_exp))
    ylabel("confirmed per $10^%i$" % np.log10(mult))
    #legend(loc='best')
    title("confirmed per capita")
    subplot(122)
    plot(USd["days"], USd["dth_pc"], sbl[0]+clr[0], mfc='none', mew=1.5, label=USd["name"])
    plot(USd["days"], USd["dth_pc_h"], "-"+clr[0], mfc='none', mew=1.5)
    for i in range(nUSloc):
        j = i+1
        locd = coronaUS[US_locs[i]]
        days = locd["days"]
        plot(days, locd["dth_pc"], sbl[j]+clr[j+k], mfc='none', mew=1.5, label=locd["name"])
        plot(days, locd["dth_pc_h"], '-'+clr[j+k])
    scaled_max = max([coronaUS[loc]['dth_pc'].max() for loc in US_locs])
    axis(xmin=-5, ymin=0-scaled_max*0.1, ymax=scaled_max*1.1)
    xlabel('days since %i confirmed per $10^%i$' % (clinv_dig, clinv_exp))
    ylabel("deaths per $10^%i$" % np.log10(mult))
    legend(loc='best')
    title("deaths per capita");
    return

def rate_US_plot(coronaUS, corona, N=1):
    # rate confirmed per capita
    US_locs = coronaUS["locs"]
    nUSloc = len(US_locs)
    dates = corona["dates"]
    mult = corona["mult"]
    clinv_dig = corona["clinv_dig"]
    clinv_exp = corona["clinv_exp"]
    USd = corona["US"]
    N+=1; figure(N, figsize=(2*fw,fh))
    clf()
    subplot(121)
    plot(USd["days"], USd["cnf_rate"], "-"+sbl[0]+clr[0], mfc='none', mew=1.5, label=USd["name"])
    for i in range(nUSloc):
        j = i+1
        locd = coronaUS[US_locs[i]]
        plot(locd["days"], locd["cnf_rate"], "-"+sbl[j]+clr[j+k], mfc='none', mew=1.5,
             label=locd["name"])
    axis(xmin = -5)
    xlabel('days since %i confirmed per $10^%i$' % (clinv_dig, clinv_exp))
    ylabel("rate [confirmed per $10^%i$ / day]" % np.log10(mult))
    #legend(loc="best")
    title("growth rate of confirmed cases")
    subplot(122)
    plot(USd["days"], USd["dth_rate"], "-"+sbl[0]+clr[0], mfc='none', mew=1.5, label=USd["name"])
    for i in range(nUSloc):
        j=i+1
        locd = coronaUS[US_locs[i]]
        plot(locd["days"], locd["dth_rate"], "-"+sbl[j]+clr[j+k], mfc='none', mew=1.5,
             label=locd["name"])
    axis(xmin = -5)
    xlabel('days since %i confirmed per $10^%i$' % (clinv_dig, clinv_exp))
    ylabel("rate [deaths per $10^%i$ / day]" % np.log10(mult))
    legend(loc="best")
    title("growth rate of deaths");
    return

def active_CFR_US_plot(coronaUS, corona, N=1):
    # active and CFR
    US_locs = coronaUS["locs"]
    nUSloc = len(US_locs)
    dates = corona["dates"]
    mult = corona["mult"]
    clinv_dig = corona["clinv_dig"]
    clinv_exp = corona["clinv_exp"]
    USd = corona["US"]
    N+=1; figure(N, figsize=(2*fw,fh))
    clf()
    subplot(121)
    plot(USd["days"], USd["acvest_pc"], "-"+sbl[0]+clr[0], mfc='none', lw=2,
         label=USd["name"])
    for i in range(nUSloc):
        j=i+1
        locd = coronaUS[US_locs[i]]
        days = locd["days"]
        plot(days, locd["acvest_pc"], "-"+sbl[j]+clr[j+k], lw=2, mfc='none',
             label=locd["name"])
    axis(xmin=-5)#, ymin=0-scaled_max*0.1, ymax=scaled_max*1.1)
    xlabel('days since %i confirmed per $10^%i$' % (clinv_dig, clinv_exp))
    ylabel("active per $10^%i$" % np.log10(mult))
    #legend(loc='best')
    title("active cases")
    subplot(122)
    plot(USd["days"], USd["dth_pc_h"]/USd["cnf_pc_h"]*100, "-"+sbl[0]+clr[0], mfc='none',
         mew=1.5, label=USd["name"])
    for i in range(nUSloc):
        j=i+1
        locd = coronaUS[US_locs[i]]
        plot(locd["days"], locd["dth_pc_h"]/locd["cnf_pc_h"]*100, "-"+sbl[j]+clr[j+k],
             mfc='none', mew=1.5, label=locd["name"])
    axis(xmin = -5, ymax=15)
    xlabel('days since %i confirmed per $10^%i$' % (clinv_dig, clinv_exp))
    ylabel("case-fatality ratio [%]")
    legend(loc="best")
    title("case fatality ratio (CFR)");
    return
