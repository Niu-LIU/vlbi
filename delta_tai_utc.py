#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: dat.py
"""
Created on Tue Dec 12 18:05:41 2017

@author: Neo(liuniu@smail.nju.edu.cn)


History
12/06/2018 : Now I use module "astropy.time" to calculate "TAI-UTC"
"""

import numpy as np
from astropy.time import Time


# -----------------------------  FUNCTIONS -----------------------------
def delta_tai_utc_calc0(mjd):
    '''For a given UTC date (MJD), calculate delta(AT) = TAI-UTC.
    '''

    # Dates and Delta(AT)s
    changes = np.array([
        [1960,  1,  1.4178180],
        [1961,  1,  1.4228180],
        [1961,  8,  1.3728180],
        [1962,  1,  1.8458580],
        [1963, 11,  1.9458580],
        [1964,  1,  3.2401300],
        [1964,  4,  3.3401300],
        [1964,  9,  3.4401300],
        [1965,  1,  3.5401300],
        [1965,  3,  3.6401300],
        [1965,  7,  3.7401300],
        [1965,  9,  3.8401300],
        [1966,  1,  4.3131700],
        [1968,  2,  4.2131700],
        [1972,  1, 10.0],
        [1972,  7, 11.0],
        [1973,  1, 12.0],
        [1974,  1, 13.0],
        [1975,  1, 14.0],
        [1976,  1, 15.0],
        [1977,  1, 16.0],
        [1978,  1, 17.0],
        [1979,  1, 18.0],
        [1980,  1, 19.0],
        [1981,  7, 20.0],
        [1982,  7, 21.0],
        [1983,  7, 22.0],
        [1985,  7, 23.0],
        [1988,  1, 24.0],
        [1990,  1, 25.0],
        [1991,  1, 26.0],
        [1992,  7, 27.0],
        [1993,  7, 28.0],
        [1994,  7, 29.0],
        [1996,  1, 30.0],
        [1997,  7, 31.0],
        [1999,  1, 32.0],
        [2006,  1, 33.0],
        [2009,  1, 34.0],
        [2012,  7, 35.0],
        [2015,  7, 36.0],
        [2017,  1, 37.0]
    ])

    year, mon, delta_tai_utc = np.transpose(changes)

    mjds = Time(["%d-%d-1" % (yi, mi)
                 for yi, mi in zip(year, mon)])

    if isinstance(mjd, float) or isinstance(mjd, int):

        index = np.searchsorted(mjds.mjd, mjd)
        if index:
            return delta_tai_utc[index - 1]
        else:
            print("Epoch too early!")

    else:
        delta = []

        for mjdi in mjd:

            index = np.searchsorted(mjds.mjd, mjdi)
            if index:
                delta.append(delta_tai_utc[index - 1])
            else:
                print("Epoch too early for %f!" % mjdi)

        return np.asarray(delta)


def delta_tai_utc_calc(mjd):
    '''For a given UTC date (MJD), calculate delta(AT) = TAI-UTC.

    Parameters
    ----------
    mjd : float, or array of float
        Modified Julian date

    Returns
    ----------
    dt : float, or array of float
        TAI-UTC in second
    '''

    t = Time(mjd, format='mjd', scale='utc')

    # t0 = Time("1972-01-01", scale='utc')

    # # If the date is earlier than 1972-01-01, I will use delta_tai_utc_calc0.
    # if mjd < t0.mjd:
    #     return delta_tai_utc_calc0

    # else:
    #     dt = Time(t.tai.iso) - Time(t.utc.iso)
    #     return dt.sec

    dt = Time(t.tai.iso) - Time(t.utc.iso)

    return dt.sec


if __name__ == '__main__':
    pass
# --------------------------------- END --------------------------------
