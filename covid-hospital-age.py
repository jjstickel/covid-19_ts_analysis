"""
analyze age-based hospitalization data from here:
https://gis.cdc.gov/grasp/COVIDNet/COVID19_3.html

Covers the states CA, CO, CT, GA, MD, MN, NM, NY, OR, TN, IA, MI, OH, and UT.
"""

import numpy as np
import pandas as pd


# must first manually delete the footer
data = pd.read_csv("../COVID-19Surveillance_All_Data.csv", skiprows=2)

# age classifications
ages = ["0-4 yr", "5-17 yr", "18-49 yr", "50-64 yr", "65+ yr"]

header = data.columns.values # useful for inspection

# "Entire Network", I presume COVID-NET = EIP + IHSP
data_etr = data[data["NETWORK"] == "COVID-NET"]

var = data_etr[data_etr["AGE CATEGORY"] == "0-4 yr"]
var2020 = var[var["MMWR-YEAR"]==2020]
