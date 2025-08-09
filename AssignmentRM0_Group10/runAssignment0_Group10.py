from utilities import read_mat_file_dates
from utilities import read_mat_file_rates
from utilities import convert_date_to_string
from utilities import zeroRates
from utilities import plot_bootstrap_results
import numpy as np
from bootstrap import bootstrap
import datetime as dt
import matplotlib.pyplot as plt



datesSet = read_mat_file_dates("datesSet.mat")
datesSetNew = {"settlement":convert_date_to_string(datesSet["settlement"]), "depos": convert_date_to_string(datesSet["depos"]),
               "futures_settlement": convert_date_to_string(datesSet["futures_settlement"]), "futures_expiry": convert_date_to_string(datesSet["futures_scadenza"]),
               "swaps": convert_date_to_string(datesSet["swaps"])}
ratesSet = read_mat_file_rates("ratesSet.mat")
print(ratesSet)
[dates, discounts] = bootstrap(datesSetNew,ratesSet)
discounts = [1] +  discounts
dates = datesSetNew["settlement"] + dates
print("dates:",dates)
discounts = np.array(discounts)
print("discounts:",discounts)
zeroRates = zeroRates(dates,discounts)*100
print("zeroRates:",zeroRates)

plot_bootstrap_results(dates, discounts)