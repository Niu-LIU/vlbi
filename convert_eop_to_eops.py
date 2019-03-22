#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: convert_eop.py
"""
Created on Tue Jun 12 10:50:55 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np
import os
import sys


# -----------------------------  FUNCTIONS -----------------------------
def print_header(ofile, sln_label=""):
    """Print the header to the output file.

    Parameters
    ----------
    ofile : string
        output file name
    sln_label : string
        solution label
    """

    # Header
    print("# VLBI Earth Orientation Solution %s\n"
          "#\n"
          "# Columns  Units   Meaning\n"
          "#    1     days    Modified Julian date(TDT)\n"
          "#    2     arcsec  X pole coordinate\n"
          "#    3     arcsec  Y pole coordinate\n"
          "#    4     sec     UT1-TAI\n"
          "#    5     mas     Celestial pole offset dX wrt IAU 2006\n"
          "#    6     mas     Celestial pole offset dY wrt IAU 2006\n"
          "#    7     arcsec  Formal uncertainty of X pole coordinate\n"
          "#    8     arcsec  Formal uncertainty of Y pole coordinate\n"
          "#    9     sec     Formal uncertainty of UT1-TAI\n"
          "#   10     mas     Formal uncertainty of celestial pole offset dX\n"
          "#   11     mas     Formal uncertainty of celestial pole offset dY\n"
          "#   12     --      10-character database identifier\n"
          "#   13     --      Correlation between X and Y\n"
          "#   14     --      Correlation between X and UT1-TAI\n"
          "#   15     --      Correlation between Y and UT1-TAI\n"
          "#   16     --      Correlation between dX and dY\n"
          "#   17     --      Number of used observations in the session\n"
          "#   18     --      6-character session identifier\n"
          "#   19     hours   Session duration\n"
          "#   20     asc/day X pole coordinate rate\n"
          "#   21     asc/day Y pole coordinate rate\n"
          "#   22     ms      Excess of length of day\n"
          "#   23             Field not used\n"
          "#   24             Field not used\n"
          "#   25     asc/day Formal uncertainty of X pole coordinate rate\n"
          "#   26     asc/day Formal uncertainty of Y pole coordinate rate\n"
          "#   27     ms      Formal uncertainty of excess of length of day\n"
          "#   28             Field not used\n"
          "#   29             Field not used\n"
          "#   30     ps      Postfit rms delay\n"
          "#   31     --      Array structure\n"
          "#\n"
          "# Note: Some EOPs were not estimated for some sessions. \n"
          "# They are, however, left in this file for sake of completeness. \n"
          "# They can be identified by their sigma set to zero.\n"
          "#" % sln_label, file=ofile)


def fill_missing(s):
    """Filling the missing value as 0.000000.

    Parameters
    ----------
    s: string

    Returns
    ----------
    0.000000 or s
    """

    if s == "*" * len(s) or s == " " * len(s):
        return " " * (len(s)-2) + "-0"
    else:
        return s


def convert_eob_to_eops(eob_file, sln_label=""):
    """Convert the format of .eob file to the same used by Sebastien

    Parameters
    ----------
    eob_file : string
        full path of .eob file
    sln_label : string
        solution label

    Returns
    ----------
    None
    """

    # Try
    if not os.path.exists(eob_file):
        print("Couldn't find the file %s" % eob_file)
        sys.exit()

    # Output file
    ofile = open("%s.eops" % (eob_file[:-4]), "w")

    # Print header
    print_header(ofile, sln_label)

    # Read .eob file
    ifile = open(eob_file, "r")
    try:
        for line in ifile:
            if not len(line) or line.startswith('#'):
                continue

            # Get what we need
            mjd = line[2:14]

            # EOP
            xp = fill_missing(line[33:41])
            yp = fill_missing(line[42:50])
            ut = fill_missing(line[51:62])
            dX = fill_missing(line[63:71])
            dY = fill_missing(line[72:80])
            xp_err = fill_missing(line[109:117])
            yp_err = fill_missing(line[118:126])
            ut_err = fill_missing(line[127:136])
            dX_err = fill_missing(line[137:144])
            dY_err = fill_missing(line[145:152])

            # database name
            db = line[15:25]

            # correlation
            xp_yp_corr = line[181:187]
            xp_ut_corr = line[188:194]
            yp_ut_corr = line[195:201]
            dX_dY_corr = line[202:208]

            # number of used observation in the session
            obs_num = line[244:250]
            sess_code = fill_missing(line[26:32])
            sess_len = line[230:235]

            # Rate of X/Y-pole
            xpr = fill_missing(line[81:90])
            ypr = fill_missing(line[91:100])
            utr = fill_missing(line[101:108])
            xpr_err = fill_missing(line[153:162])
            ypr_err = fill_missing(line[163:172])
            utr_err = fill_missing(line[173:180])

            rms = line[236:243]
            net = line[264:328].strip("\n")

            print(mjd, xp, yp, ut, dX, dY,
                  xp_err, yp_err, ut_err, dX_err, dY_err,
                  db, xp_yp_corr, xp_ut_corr, yp_ut_corr, dX_dY_corr,
                  obs_num, sess_code, sess_len,
                  xpr, ypr, utr, "-0", "-0",
                  xpr_err, ypr_err, utr_err, "-0", "-0",
                  rms, net, file=ofile)

    finally:
        ifile.close()
        ofile.close()


# --------------------------------- MAIN -------------------------------
if __name__ == "__main__":
    if len(sys.argv) == 2:
        convert_eob_to_eops(sys.argv[1])
    elif len(sys.argv) == 3:
        convert_eob_to_eops(sys.argv[1], sys.argv[2])
    else:
        print("Input error!")
        exit()
# --------------------------------- END --------------------------------
