"""
Mathematical Engineering - Financial Engineering, FY 2024-2025
Risk Management - Exercise 2: Corporate Bond Portfolio
"""

from typing import Union
import numpy as np
import pandas as pd
import datetime as dt
from ex1_utilities import (
    year_frac_act_x,
    business_date_offset,

)
from ex1_utilities import (get_discount_factor_by_zero_rates_linear_interp)
from scipy.optimize import fsolve


def bond_cash_flows(
    ref_date: Union[dt.date, pd.Timestamp],
    expiry: Union[dt.date, pd.Timestamp],
    coupon_rate: float,
    coupon_freq: int,
    notional: float = 1.0,
) -> pd.Series:
    """
    Calculate the cash flows of a bond.

    Parameters:
    ref_date (Union[dt.date, pd.Timestamp]): Reference date.
    expiry (Union[dt.date, pd.Timestamp]): Bond's expiry date.
    coupon_rate (float): Coupon rate.
    coupon_freq (int): Coupon frequency in payments per years.
    notional (float): Notional amount.

    Returns:
        pd.Series: Bond cash flows.
    """

    # Payment dates
    cash_flows_dates = []
    payment_dt = ref_date
    counter = 1
    while expiry > payment_dt:

        payment_dt = business_date_offset(
            ref_date, month_offset=(12 * counter) // coupon_freq
        )
        cash_flows_dates.append(payment_dt)
        counter += 1

        
    '''# Payment dates
    cash_flows_dates = []
    payment_dt = expiry
    counter = 1
    while payment_dt > ref_date:
        cash_flows_dates.append(payment_dt)

        payment_dt = business_date_offset(
            expiry, month_offset=(-12 * counter) // coupon_freq
        )
        counter += 1

    cash_flows_dates.sort()'''
    # Coupon payments
    cash_flows = pd.Series(
        data=[notional * coupon_rate ] * year_frac_act_x([ref_date]+cash_flows_dates[:-1], cash_flows_dates,6),
        index=cash_flows_dates,
    )

    # Notional payment
    cash_flows[expiry] += notional

    return cash_flows






def defaultable_bond_dirty_price_from_z_spread(
    ref_date: Union[dt.date, pd.Timestamp],
    expiry: Union[dt.date, pd.Timestamp],
    coupon_rate: float,
    coupon_freq: int,
    z_spread: float,
    dates: pd.Series,
    discount_factors: pd.Series,
    notional: float = 1.0,
) -> float:
    """
    Calculate the dirty price of a defaultable bond from the Z-spread.

    Parameters:
    ref_date (Union[dt.date, pd.Timestamp]): Reference date.
    expiry (Union[dt.date, pd.Timestamp]): Bond's expiry date.
    coupon_rate (float): Coupon rate.
    coupon_freq (int): Coupon frequency in payments a years.
    z_spread (float): Z-spread.
    discount_factors (pd.Series): Discount factors.
    notional (float): Notional amount.

    Returns:
        float: Dirty price of the bond.
    """

    # !!! COMPLETE THE FUNCTION !!!
    ref_date = pd.to_datetime(ref_date)
    cash_flows = bond_cash_flows(
        ref_date,
        expiry,
        coupon_rate,coupon_freq,notional)


    # we compute the discounts at the relevant dates for cash flows
    discounts_at_cash_flows = [
        get_discount_factor_by_zero_rates_linear_interp(ref_date, date, dates, discount_factors)[0]
        for date in cash_flows.index.tolist()]


    days_diff = np.array([(cash_flows.index[i] - ref_date).days for i in range(len(cash_flows))])

    #we compute the dirty price
    res = np.sum(np.array(cash_flows.values) * np.array(discounts_at_cash_flows) * np.exp(-z_spread*np.array(days_diff)/365))

    return res


def constant_intensity_to_be_solved(cash_flows: pd.Series,
                                    ref_date:Union[dt.date, pd.Timestamp],
                                    discounts_at_cash_flows: list[float],
                                    notional:float,
                                    dirty_price:float,
                                    recovery_rate:float,
                                    intensity:float) -> float:



    # Convert ref_date to a pandas Timestamp if it's a date object
    ref_date = pd.to_datetime(ref_date)

    # Convert the cash flows index to a series of days from the reference date
    days_diff = [(cash_flows.index[i] - ref_date).days/365 for i in range(len(cash_flows))]
    # First sum: Vectorized operation using element-wise multiplication
    first_sum = np.sum(np.array(cash_flows.values) * np.array(discounts_at_cash_flows) * np.exp(-intensity * days_diff))

    # Complete dates for second sum calculation
    complete_dates = [ref_date] + list(cash_flows.index)

    days_diff_complete = [(complete_dates[i] - ref_date).days/365 for i in range(len(complete_dates))]
    # Second sum: Vectorized calculation using the difference of exponential terms
    second_sum = np.sum(np.array(discounts_at_cash_flows) * (
                np.exp(-intensity * np.array(days_diff_complete[:-1])) - np.exp(
            -intensity * np.array(days_diff_complete[1:]))))



    # Calculate and return the final result
    return dirty_price -first_sum -recovery_rate * second_sum


def constant_intensity_from_bond (
    ref_date: Union[dt.date, pd.Timestamp],
    expiry: Union[dt.date, pd.Timestamp],
    coupon_rate: float,
    coupon_freq: int,
    dates: list[dt.datetime],
    discount_factors: list[float],
    dirty_price:float,
    recovery_rate:float,
    end: float,
    notional: float = 1.0,

)->float:
    #we compute the cash flow and their dates for the input bond
    cash_flows = bond_cash_flows(ref_date,expiry,coupon_rate,coupon_freq)
    cash_flows = cash_flows[end:]

    #we compute the discounts at the relevant dates for cash flows
    discounts_at_cash_flows = [
        get_discount_factor_by_zero_rates_linear_interp(ref_date, date, dates, discount_factors)[0]
        for date in cash_flows.index.tolist()]

    #we use fsolve to find lambda
    def to_be_minimized(intensity):
        return constant_intensity_to_be_solved(cash_flows,ref_date, discounts_at_cash_flows,notional, dirty_price,recovery_rate,intensity)

    return fsolve(to_be_minimized, 0.000001)

def piece_wise_constant_intensity_from_bonds(ref_date: Union[dt.date, pd.Timestamp],
    expiry: list[Union[dt.date, pd.Timestamp]],
    coupon_rate: list[float],
    coupon_freq: list[int],
    dates: list[dt.datetime],
    discount_factors: list[float],
    dirty_price: list[float],
    recovery_rate: float,

) ->list[float]:

  num_bond = len(expiry)

  intensity_vector = np.zeros(num_bond)
  computed_intensities = 0
  start = 0
  end = 0
  for i in range(0,num_bond):
   dirty_price_modifier = 0
   cash_flows = bond_cash_flows(ref_date, expiry[i], coupon_rate[i], coupon_freq[i])
   discounts_at_cash_flows = [
       get_discount_factor_by_zero_rates_linear_interp(ref_date, date, dates, discount_factors)[0]
       for date in cash_flows.index.tolist()]
   days_diff = np.array([(cash_flows.index[i] - ref_date).days/365 for i in range(len(cash_flows))])

   for j in range(1,computed_intensities+1):
     max_index = np.sum(np.array(cash_flows.index <expiry[j-1]))
     end = max_index
     dirty_price_modifier += sum(cash_flows.values[start:end]*discounts_at_cash_flows[start:end] \
     *np.exp(-intensity_vector[j-1]*days_diff[start:end]))
     days_diff_completo =[0]+days_diff
     dirty_price_modifier += recovery_rate*sum(discounts_at_cash_flows[start:end]*(np.exp(-intensity_vector[j-1]*days_diff_completo[start:end-1])
                                                                                   -np.exp(-intensity_vector[j-1]*days_diff_completo[start+1:end])))
     start = max_index

   intensity_vector[i] = constant_intensity_from_bond(ref_date,
   expiry[i],
   coupon_rate[i],
   coupon_freq[i],
   dates,
   discount_factors,
   dirty_price[i]-dirty_price_modifier,
   recovery_rate,end
   )

   computed_intensities +=1

  return intensity_vector
