#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 10:52:21 2017

@author: Neo

Retrieve the estimates of positions of global stations and the formal
uncertainties of these estimates from .sta file which is generated by
the program getpar.

"""

from astropy import units as u
from astropy.time import Time
from astropy.table import Table, Column
import numpy as np
import sys


__all__ = ["calc_sess_epoch_from_name", "read_rms"]


# ------------------------------  FUNCTIONS  ---------------------------
def calc_sess_epoch_from_name(sess_id):
    """Calculate the epoch from session identifiers.
    """

    # Mapping table between letter and number for month
    mon_str2num = {
        "JAN": "01",
        "FEB": "02",
        "MAR": "03",
        "APR": "04",
        "MAY": "05",
        "JUN": "06",
        "JUL": "07",
        "AUG": "08",
        "SEP": "09",
        "OCT": "10",
        "NOV": "11",
        "DEC": "12"
    }

    year = int(sess_id[1:3])
    mon_str = sess_id[3:6]
    date = sess_id[6:8]

    if year >= 79:
        year += 1900
    else:
        year += 2000

    epoch = Time(
        "{}-{}-{}".format(year, mon_str2num[mon_str], date), format="iso")

    return epoch.jyear


def read_rms(rmsfile):
    '''Retrieve the result from .rms file.

    Parameters
    ----------
    rmsfile : string
        name of data file

    Returns
    ----------
    dbname : array, string
       database name with leading dollar sign
    obsnum : array, int
        number of observations used
    epo : array, float
        epoch
    wrmsd : array, float
        overall wrms of postfit delay residuals, psec
    wrmsr : array, float
        overall wrms of postfit delay rate residuals, fsec/s
    '''

    rmstable = Table.read(rmsfile,
                          format="ascii.fixed_width_no_header",
                          data_start=3,
                          col_starts=[10, 22, 31, 48],
                          col_ends=[20, 28, 42, 58],
                          names=["sess_name", "obs_num", "delay_rms", "delay_rate_rms"])

    # Add unit information
    rmstable["delay_rms"].unit = u.second / 1e12
    rmstable["delay_rate_rms"].unit = 1 / 1e15

    # Calculate the epoch
    epo = [calc_sess_epoch_from_name(db) for db in rmstable["sess_name"]]

    epo_col = Column(np.array(epo) << u.yr, name="epoch")
    rmstable.add_column(epo_col, index=0)
    rmstable.sort(["epoch"])

    return rmstable


def main():
    print("Nothing to do :-)")


if __name__ == "__main__":
    main()
# ------------------------------ END -----------------------------------
