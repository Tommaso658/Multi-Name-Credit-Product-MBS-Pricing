"""
Mathematical Engineering - Financial Engineering, FY 2024-2025
Risk Management - Exercise 1: Hedging a Swaption Portfolio
"""
import datetime
from enum import Enum
import numpy as np
import pandas as pd
import datetime as dt
import calendar
from scipy.stats import norm
from scipy.optimize import fsolve
from bootstrap import bootstrap
from utilities import yearfrac

from typing import Iterable, Union, List, Union, Tuple


class SwapType(Enum):
    """
    Types of swaptions.
    """

    RECEIVER = "receiver"
    PAYER = "payer"


def year_frac_act_x(start_dates: list[dt.datetime] , end_dates: list[dt.datetime], flag: int):
    """
    Compute the year fraction between two dates using the ACT/x convention.

    Parameters:
        t1 (dt.datetime): First date.
        t2 (dt.datetime): Second date.
        x (int): Number of days in a year.

    Returns:
        float: Year fraction between the two dates.
    """

    if len(start_dates) == 1:
        actual_days = np.array([(end - start_dates[0]).days for end in end_dates])
    else:
        actual_days = np.array([(end - start).days for start, end in zip(start_dates, end_dates)])

    if flag == 2:
        return actual_days / 360

    if flag == 3:
        return actual_days / 365

    if flag == 6:
        d1 = np.minimum(30, np.array([date.day for date in start_dates]))
        d2 = np.minimum(30, np.array([date.day for date in end_dates]))
        if len(start_dates) == 1:
            year_diff = np.array([end.year - start_dates[0].year for end in end_dates])
            month_diff = np.array([end.month - start_dates[0].month for end in end_dates])
        else:
            year_diff = np.array([end.year - start.year for start, end in zip(start_dates, end_dates)])
            month_diff = np.array([end.month - start.month for start, end in zip(start_dates, end_dates)])
        return (360 * year_diff + 30 * month_diff + (d2 - d1)) / 360


def from_discount_factors_to_zero_rates(
    dates: Union[List[float], pd.DatetimeIndex],
    discount_factors: Iterable[float],
) -> List[float]:
    """
    Compute the zero rates from the discount factors.

    Parameters:
        dates (Union[List[float], pd.DatetimeIndex]): List of year fractions or dates.
        discount_factors (Iterable[float]): List of discount factors.

    Returns:
        List[float]: List of zero rates.
    """

    effDates, effDf = dates, discount_factors
    if isinstance(effDates, pd.DatetimeIndex):
        effDates = [
            year_frac_act_x(effDates[i - 1], effDates[i], 3)
            for i in range(1, len(effDates))
        ]
        effDf = discount_factors[1:]

    return -np.log(np.array(effDf)) / np.array(effDates)


def get_discount_factor_by_zero_rates_linear_interp(
    reference_date: Union[dt.datetime, pd.Timestamp],
    interp_date: Union[dt.datetime, pd.Timestamp],
    dates: Union[List[dt.datetime], pd.DatetimeIndex],
    discount_factors: Iterable[float],
) -> float:
    """
    Given a list of discount factors, return the discount factor at a given date by linear
    interpolation.

    Parameters:
        reference_date (Union[dt.datetime, pd.Timestamp]): Reference date.
        interp_date (Union[dt.datetime, pd.Timestamp]): Date at which the discount factor is
            interpolated.
        dates (Union[List[dt.datetime], pd.DatetimeIndex]): List of dates.
        discount_factors (Iterable[float]): List of discount factors.

    Returns:
        float: Discount factor at the interpolated date.
    """

    if len(dates) != len(discount_factors):
        raise ValueError("Dates and discount factors must have the same length.")

    year_fractions = year_frac_act_x([reference_date], dates[1:], 3)
    zero_rates = from_discount_factors_to_zero_rates(
        year_fractions, discount_factors[1:]
    )
    inter_year_frac = year_frac_act_x([reference_date], [interp_date], 3)
    rate = np.interp(inter_year_frac, year_fractions, zero_rates)
    return np.exp(-inter_year_frac * rate)


def business_date_offset(
    base_date: Union[dt.date, pd.Timestamp],
    year_offset: int = 0,
    month_offset: int = 0,
    day_offset: int = 0,
) -> Union[dt.date, pd.Timestamp]:
    """
    Return the closest following business date to a reference date after applying the specified offset.

    Parameters:
        base_date (Union[dt.date, pd.Timestamp]): Reference date.
        year_offset (int): Number of years to add.
        month_offset (int): Number of months to add.
        day_offset (int): Number of days to add.

    Returns:
        Union[dt.date, pd.Timestamp]: Closest following business date to ref_date once the specified
            offset is applied.
    """

    # Adjust the year and month
    total_months = base_date.month + month_offset - 1
    year, month = divmod(total_months, 12)
    year += base_date.year + year_offset
    month += 1

    # Adjust the day and handle invalid days
    day = base_date.day
    try:
        adjusted_date = base_date.replace(
            year=year, month=month, day=day
        ) + dt.timedelta(days=day_offset)
    except ValueError:
        # Set to the last valid day of the adjusted month
        last_day_of_month = calendar.monthrange(year, month)[1]
        adjusted_date = base_date.replace(
            year=year, month=month, day=last_day_of_month
        ) + dt.timedelta(days=day_offset)

    # Adjust to the closest business day
    if adjusted_date.weekday() == 5:  # Saturday
        adjusted_date += dt.timedelta(days=2)
    elif adjusted_date.weekday() == 6:  # Sunday
        adjusted_date += dt.timedelta(days=1)

    return adjusted_date


def date_series(
    t0: Union[dt.date, pd.Timestamp], t1: Union[dt.date, pd.Timestamp], freq: int
) -> Union[List[dt.date], List[pd.Timestamp]]:
    """
    Return a list of dates from t0 to t1 inclusive with frequency freq, where freq is specified as
    the number of dates per year.
    """

    dates = [t0]
    print(dates[-1])
    while dates[-1] < t1:
        dates.append(business_date_offset(t0, month_offset=len(dates) * 12 // freq))

    if dates[-1] > t1:
        dates.pop()
    if dates[-1] != t1:
        dates.append(t1)

    return dates


def swaption_price_calculator(
    dates: list[dt.datetime],
    S0: float,
    strike: float,
    ref_date: Union[dt.date, pd.Timestamp],
    expiry: Union[dt.date, pd.Timestamp],
    underlying_expiry: Union[dt.date, pd.Timestamp],
    sigma_black: float,
    freq: int,
    discount_factors: pd.Series,
    swaption_type: SwapType = SwapType.RECEIVER,
    compute_delta: bool = False,
) -> Union[float, Tuple[float, float]]:
    """
    Return the swaption price defined by the input parameters.

    Parameters:
        S0 (float): Forward swap rate.
        strike (float): Swaption strike price.
        ref_date (Union[dt.date, pd.Timestamp]): Value date.
        expiry (Union[dt.date, pd.Timestamp]): Swaption expiry date.
        underlying_expiry (Union[dt.date, pd.Timestamp]): Underlying forward starting swap expiry.
        sigma_black (float): Swaption implied volatility.
        freq (int): Number of times a year the fixed leg pays the coupon.
        discount_factors (pd.Series): Discount factors.
        swaption_type (SwaptionType ["receiver", "payer"]): Swaption type, default to receiver.
        swaption_type (SwapType): Swaption type, default to receiver.

    Returns:
        Union[float, Tuple[float, float]]: Swaption price (and possibly delta).
    """
    #calculate fixed leg schedule
    fixed_leg_schedule = date_series(expiry, underlying_expiry, freq)
    #calculate B(0,tn), maybe useless
    #discount_tn = get_discount_factor_by_zero_rates_linear_interp(ref_date, expiry, dates, discount_factors)
    #compute BPV
    bpv = sum([year_frac_act_x([fixed_leg_schedule[i - 1]], [fixed_leg_schedule[i]],
                         6) * get_discount_factor_by_zero_rates_linear_interp(ref_date, fixed_leg_schedule[i], dates,
                         discount_factors) for i in range(1, len(fixed_leg_schedule))])
    T = (expiry - ref_date).days / 365.0

    # Calculate d1 and d2 using Blackformula
    d1 = (np.log(S0 / strike) + 0.5 * sigma_black ** 2 * T) / (sigma_black * np.sqrt(T))
    d2 = d1 - sigma_black * np.sqrt(T)

    # Calculate the price based on the swaption type
    if swaption_type == SwapType.RECEIVER:
        price = bpv*(-S0 * norm.cdf(-d1) +strike * norm.cdf(-d2))
        delta = -bpv * norm.cdf(-d1)
    else:
       price = bpv*(S0 * norm.cdf(d1) - strike * norm.cdf(d2))
       delta = bpv * norm.cdf(d1)
    if compute_delta:
        return [price,delta]
    else:
        return price



def swap_par_rate(dates: List[dt.datetime],
    fixed_leg_schedule: List[dt.datetime],
    discount_factors: pd.Series,
    fwd_start_date: dt.datetime | None = None,
) -> float:
    """
    Given a fixed leg payment schedule and the discount factors, return the swap par rate. If a
    forward start date is provided, a forward swap rate is returned.

    Parameters:
        dates: List[dt.datetime]: all dates corresponding to discounts
        fixed_leg_schedule (List[dt.datetime]): Fixed leg payment dates.
        discount_factors (pd.Series): Discount factors.
        fwd_start_date (dt.datetime | None): Forward start date, default to None.

    Returns:
        float: Swap par rate.
    """

    ### !!! MODIFY AS APPROPRIATE !!!
    discount_factor_t0 = 1 if fwd_start_date is None else get_discount_factor_by_zero_rates_linear_interp(
    dates[0],
    fwd_start_date,
    dates,
    discount_factors,
    )

    # !!! MODIFY AS APPROPRIATE !!!

    bpv = sum([year_frac_act_x([fixed_leg_schedule[i-1]],[fixed_leg_schedule[i]],6)*get_discount_factor_by_zero_rates_linear_interp(dates[0],fixed_leg_schedule[i],dates,
                                                               discount_factors) for i in range(1,len(fixed_leg_schedule))])

    discount_factor_tN = get_discount_factor_by_zero_rates_linear_interp(
    dates[0],
    fixed_leg_schedule[-1],
    dates,
    discount_factors,
    )

    return (discount_factor_t0 - discount_factor_tN) / bpv


def swap_mtm(
    dates,
    swap_rate: float,
    fixed_leg_schedule: List[dt.datetime],
    discount_factors: pd.Series,
    swap_type: SwapType = SwapType.PAYER,
) -> float:
    """
    Given a swap rate, a fixed leg payment schedule and the discount factors, return the swap
    mark-to-market.

    Parameters:
        swap_rate (float): Swap rate.
        fixed_leg_schedule (List[dt.datetime]): Fixed leg payment dates.
        discount_factors (pd.Series): Discount factors.
        swap_type (SwapType): Swap type, either 'payer' or 'receiver', default to 'payer'.

    Returns:
        float: Swap mark-to-market.
    """
    # Single curve framework, returns price and basis point value
    discount_factors = pd.Series(discount_factors)
    bpv = sum([year_frac_act_x([fixed_leg_schedule[i - 1]], [fixed_leg_schedule[i]],
                         6) * get_discount_factor_by_zero_rates_linear_interp(dates[0], fixed_leg_schedule[i], dates,
                         discount_factors) for i in range(1, len(fixed_leg_schedule))])
    P_term = get_discount_factor_by_zero_rates_linear_interp(
        dates[0],
        fixed_leg_schedule[-1],
        dates,
        discount_factors.values,
    )
    float_leg = 1.0 - P_term
    fixed_leg = swap_rate * bpv

    if swap_type == SwapType.RECEIVER:
        multiplier = -1
    elif swap_type == SwapType.PAYER:
        multiplier = 1
    else:
        raise ValueError("Unknown swap type.")

    return multiplier * (float_leg - fixed_leg)


def irs_proxy_duration(
    ref_date: dt.date,
    swap_rate: float,
    fixed_leg_payment_dates: List[dt.date],
    discount_factors: List[float],
    dates: List[dt.date],
) -> float:
    """
    Given the specifics of an interest rate swap (IRS), return its rate sensitivity calculated as
    the duration of a fixed coupon bond.

    Parameters:
        ref_date (dt.date): Reference date.
        swap_rate (float): Swap rate.
        fixed_leg_payment_dates (List[dt.date]): Fixed leg payment dates.
        discount_factors (pd.Series): Discount factors.

    Returns:
        (float): Swap duration.
    """

    #Compute discount factor
    discounts_interp = np.zeros(len(fixed_leg_payment_dates))
    for i in range(len(fixed_leg_payment_dates)):
        discounts_interp[i] = get_discount_factor_by_zero_rates_linear_interp(ref_date, fixed_leg_payment_dates[i], dates, discount_factors)
    discounts_interp = discounts_interp[1:]
    #Compute the year_frac
    yf = year_frac_act_x(fixed_leg_payment_dates[:-1],fixed_leg_payment_dates[1:],6) #check flag
    #Compute the coupon for every time i
    c_i = yf*swap_rate
    c_i[-1] += 1
    #Comupte Interbank Bonds at time 0 (must be 1)
    IB_bond_t0 = np.sum(c_i*discounts_interp)
    #print("IB_bond_t0:", IB_bond_t0)
    #Compute the Maculay Duration
    times = year_frac_act_x([ref_date],fixed_leg_payment_dates,6)
    times = times[1:]
    Sium = np.sum(c_i*discounts_interp*times)
    DU = Sium/IB_bond_t0

    return DU



