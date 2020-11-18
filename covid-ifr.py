"""
Calculate overall IFR
data from:
https://www.cdc.gov/coronavirus/2019-ncov/hcp/planning-scenarios.html
https://www.census.gov/data/tables/2019/demo/age-and-sex/2019-age-sex-composition.html
"""

import numpy as np


age_ifr = np.array([0, 20, 50, 70])
ifr_frac = np.array([0.00003, 0.0002, 0.005, 0.054])

age_pop = np.arange(0,86,5)
pop_frac = np.array([6.1, 6.2, 6.4, 6.4, 6.6, 7.2, 6.8, 6.6, 6.0, 6.3, 6.3, 6.5, 6.3, 5.4, 4.4, 2.9, 1.9, 1.8])/100

pop_frac2 = []
idx1 = 0
for age in age_ifr[1:]:
    #print(age)
    idx2 = np.nonzero(age_pop<=age)[0][-1]
    pop_frac2.append(pop_frac[idx1:idx2].sum())
    idx1 = idx2
pop_frac2.append(pop_frac[idx1:].sum())
pop_frac2 = np.array(pop_frac2)

ifr_tot = np.sum(pop_frac2*ifr_frac)
print("estimated total population IFR = %g" % ifr_tot)

# estimated total population IFR = 0.00725153
