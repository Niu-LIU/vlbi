#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: check_eob_wrt_c04.py
"""
Created on Tue Dec 12 16:19:31 2017

@author: Neo(liuniu@smail.nju.edu.cn)

History
16/06/2018 : Now use .eops file to replace .eob file
"""

import numpy as np
import sys

# My modules
from read_eob import read_eops
from delta_tai_utc import delta_tai_utc_calc
from eops_diff_stats import eops_offset_stats


__all__ = ["eops_offset_calc", "check_eob_wrt_c04"]


# -----------------------------  FUNCTIONS -----------------------------
def eops_offset_calc(EOPS_file):
    '''Calculate the EOP and Nutation difference wrt. C04.

    Parameters
    ----------
    EOPS_file : string
        full path of .eops file

    Returns
    ----------
    None
    '''

    # A prioir EOP file based on C04 series
    APR_file = "%s.eops_c04_apr" % EOPS_file[:-5]

    # Output file
    DIF_file = "%s.eops_c04_dif" % EOPS_file[:-5]

    [db_name, mjd, obs_num, sess_len, rms,
     xp, yp, ut, dx, dy,
     xp_err, yp_err, ut_err, dx_err, dy_err,
     xp_yp_corr, xp_ut_corr, yp_ut_corr, dx_dy_corr,
     xpr, ypr, utr, xpr_err, ypr_err, utr_err] = read_eops(EOPS_file)

    mjd_apr, xp_apr, yp_apr, ut_apr, dx_apr, dy_apr = np.genfromtxt(
        APR_file, unpack=True)

    # print(mjd, mjd_apr)

    # if all(mjd != mjd_apr):
    #     print("These two files do not match!")
    #     exit()

    # for (mjdi, mjd_apri) in zip(mjd, mjd_apr):
    #     if mjdi != mjd_apri:
    #         print(mjdi, mjd_apri)

    # exit()

    # In C04, Uapr is UT1-TAI while it is UT1-UTC in .eops file,
    # so we need to calculate IAI - UTC and then convert Uapr to UT1-UTC
    # IAI - UTC
    delta_tai_utc = delta_tai_utc_calc(mjd)
    ut += delta_tai_utc * 15000  # unit: mas

    # compare
    fout = open(DIF_file, "w")
    print("# EOP and Nutation difference wrt C04:\n"
          "# Epoch  dXp  err  dYp  err  dUT  err  ddX  err  ddY  err\n"
          "# mjd    mas  mas  mas  mas  mas  mas  mas  mas  mas  mas",
          file=fout)

    # difference
    for (i, mjd) in enumerate(mjd_apr):
        print("%13.6f" % mjd, "  %+10.3f  %7.3f" * 5 %
              (xp[i] - xp_apr[i], xp_err[i], yp[i] - yp_apr[i], yp_err[i],
               ut[i] - ut_apr[i], ut_err[i], dx[i] - dx_apr[i], dx_err[i],
               dy[i] - dy_apr[i], dy_err[i]), file=fout)

    fout.close()


def check_eops_wrt_c04(EOPS_file):
    '''Check EOP and Nutation parameter wrt. C04 series.

    Parameters
    ----------
    EOPS_file : string
        full path of .eops file

    Returns
    ----------
    None
    '''

    # Calculate offset
    eops_offset_calc(EOPS_file)

    # Statistical result
    eops_offset_stats("%s.eops_c04_dif" % EOPS_file[:-5])


# -------------------- MAIN ----------------------------------
if len(sys.argv) == 2:
    check_eops_wrt_c04(sys.argv[1])
else:
    print("Input error!", sys.argv)
    exit()
# --------------------------------- END --------------------------------
