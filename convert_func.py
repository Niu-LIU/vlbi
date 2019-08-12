#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: convert_func.py
"""
Created on Thu Mar 21 10:58:04 2019

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np
from astropy.time import Time


__all__ = ["RA_conv", "DC_conv", "date2jyear", "date2mjd"]


# -----------------------------  FUNCTIONS -----------------------------
def RA_conv(RAstr):
    """Convert right ascension string of HH_MM_SS.ssssssss into float.

    Parameter
    ----------
    RAstr : str
        string of right ascension in format of HH_MM_SS.ssssssss

    Return
    ------
    RA : float
        right ascension in degree
    """

    hours, mins, secs = RAstr[4:].split("_")
    RA = (float(hours) + float(mins) / 60. + float(secs)/3600.) * 15
    return RA


def DC_conv(DCstr):
    """Convert declination string of +DD_AM_AS.ssssssss into float.

    Parameter
    ----------
    DCstr : str
        string of declination in format of +DD_AM_AS.ssssssss

    Return
    ------
    DC : float
        declination in degree
    """

    degs, ams, ass = DCstr[4:].split("_")

    # Determine the sign
    if DCstr[4] == "-":
        DC = float(degs) - float(ams) / 60. - float(ass) / 3600.
    else:
        DC = float(degs) + float(ams) / 60. + float(ass) / 3600.

    return DC


def date2jyear(s):
    """Convert a string of date into a float value in unit of year.

    Parameters
    ----------
    datafile : string
        date as "<YYYY>.<MM>.<DD>"

    Returns
    ----------
    year : float
        date in unit of year
    """

    if s == " " * 10:
        return 0

    t = Time(s.replace(".", "-"), scale="utc")

    return t.jyear


def date2mjd(s):
    """Convert a string of date into a float value in unit of year.

    Parameters
    ----------
    datafile : string
        date as "<YYYY>.<MM>.<DD>"

    Returns
    ----------
    year : float
        date in unit of year
    """

    if s == " " * 10:
        return 0

    t = Time(s.replace(".", "-"), scale="utc")

    return t.mjd


# --------------------------------- END --------------------------------
