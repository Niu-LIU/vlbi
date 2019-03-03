#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: read_spm.py
"""
Created on Sun Jun 24 23:39:15 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""

from astropy.table import Table, join, unique, Column
import astropy.units as u
from astropy.units import cds
from astropy.coordinates import Angle
import numpy as np
import os
import sys
import time
# My module
from convert_func import RA_conv, DC_conv, date2mjd

__all__ = ["write_pmt", "read_pmt", "get_spm"]


# -----------------------------  FUNCTIONS -----------------------------
def print_pmt_header(ofile):
    """Print the header to the .pmt file.

    Parameters
    ----------
    ofile : file handle

    Returns
    ----------
    None
    """

    # Header
    print("#\n"
          "# Columns  Units   Meaning\n"
          "#    1     --      IVS designation\n"
          "#    2     mas/yr  Proper motion in right ascension\n"
          "#    3     mas/yr  Proper motion in declination\n"
          "#    4     mas/yr  Formal uncertainty of proper motion in right"
          " ascension\n"
          "#    5     mas/yr  Formal uncertainty of proper motion in "
          "declination\n"
          "#    6     --      Correlation between proper motions in right "
          "ascension and declination\n"
          "#\n"
          "# Note: Some proper motion component were not estimated for some "
          "sources. \n"
          "# They are however left in this file for sake of completeness. \n"
          "# They can be identified by their sigma set to 999.000000 or"
          " 206.264806.\n"
          # "# For correlation, "
          "#", file=ofile)


def check_info(pm_ra, pm_dec, pm_ra_err, pm_dec_err, corr,
               read_ra, read_dec, read_corr):
    """Check if the proper motion information for a source is complete.

    If no information provided in the spool file, these parameters will be set as some
    default value for further distinguish.

    Parameters
    ----------
    pm_ra/pm_dec : string
        Proper motion in right ascension (*cos(Dec))/Declination
    pm_ra_err/pm_dec_err : string
        Formal uncertainty of pm_ra/pm_dec
    corr : string
        Correlation between pm_ra and pm_dec
    read_ra/read_dec/read_corr : boolean
        Flag to marker if pm_ra/pm_dec/corr information is already read.

    Returns
    -------
    pm_ra/pm_dec : string
        Proper motion in right ascension (*cos(Dec))/Declination.
        If no pm_ra/pm_dec is provided, they will be set as "     0.000000".
    pm_ra_err/pm_dec_err : string
        Formal uncertainty of pm_ra/pm_dec.
        If no pm_ra_err/pm_dec_err is provided, they will be set as "999.000000".
    corr : string
        Correlation between pm_ra and pm_dec.
        If no corr is provided, they will be set as " 0.0000".
    """

    if not read_ra:
        # There is no information for R.A. component
        pm_ra = "     0.000000"
        pm_ra_err = "999.000000"
    elif pm_ra_err == "*" * 10:
        pm_ra_err = "999.000000"

    if not read_dec:
        # There is no information for Dec. component
        pm_dec = "     0.000000"
        pm_dec_err = "999.000000"
    elif pm_dec_err == "*" * 10:
        pm_dec_err = "999.000000"

    if not read_corr:
        # There is no information for correlation
        # (normally, this couldn't happen!)
        corr = " 1.9999"

    return [pm_ra, pm_dec, pm_ra_err, pm_dec_err, corr]


def write_pmt(spl_file, pmt_file):
    """Load the proper motion information and write into file.

    Parameters
    ----------
    spl_file : string
        spool file
    pmt_file : string
        source proper motion file
    """

    # output
    ofile = open(pmt_file, "w")

    # Print header
    print("#original file: ", spl_file, file=ofile)
    print_pmt_header(ofile)

    read_ra = False
    read_dec = False
    read_corr = False

    with open("pmt.temp", "r") as ftmp:
        line = ftmp.readline()
        isou = line[8:16]
        while len(line):
            if line[8:16] != isou:
                # Here we get a new source

                # Check the completeness of proper motion of the old source
                pm_ra, pm_dec, pm_ra_err, pm_dec_err, corr = check_info(
                    pm_ra, pm_dec, pm_ra_err, pm_dec_err, corr,
                    read_ra, read_dec, read_corr)
                # Print proper motion
                print(isou, pm_ra, pm_dec, pm_ra_err,
                      pm_dec_err, corr, file=ofile)

                # Renew source name and flags
                isou = line[8:16]
                read_ra = False
                read_dec = False
                read_corr = False

            if line[17:33] == "RT. ASC.velocity":
                pm_ra = line[33:46]
                pm_ra_err = line[71:81]
                read_ra = True
            elif line[17:33] == "DEC.    velocity":
                pm_dec = line[33:46]
                pm_dec_err = line[71:81]
                read_dec = True
            elif line[17:28] == "CORRELATION":
                corr = line[32:39]
                read_corr = True
            else:
                print("There is something wrong in line:\n", line)
                exit()

            line = ftmp.readline()

        # Check the completeness of the information of the last source
        pm_ra, pm_dec, pm_ra_err, pm_dec_err, corr = check_info(
            pm_ra, pm_dec, pm_ra_err, pm_dec_err, corr,
            read_ra, read_dec, read_corr)
        # Print proper motion of the last source
        print(isou, pm_ra, pm_dec, pm_ra_err, pm_dec_err, corr, file=ofile)

    ofile.close()


def retrieve_pmt(spl_file, pmt_file):
    """Read source proper motions from spool file.

    Parameters
    ----------
    spl_file : string
        spool file
    pmt_file : string
        output file file
    """

    tmpfile = "pmt.temp"

    # Get the related lines in spool file
    # To see if we are in area of proper motion
    find_pm = False

    # temp output
    fout = open(tmpfile, "w")

    # read file line by line
    with open(spl_file, "r") as fspl:
        while True:
            line = fspl.readline()
            if line[:6] == "1  APR" or not line:
                break

            if "RT. ASC.velocity" in line or "DEC.    velocity" in line:
                print(line, end="", file=fout)
                find_pm = True
            elif "CORRELATION" in line:
                if find_pm:
                    print(line, end="", file=fout)

    fout.close()

    # rewrite the proper motion
    write_pmt(spl_file, pmt_file)


# def read_pmt(ifile):
#     """Read the .spm file.

#     Parameters
#     ----------
#     ifile : string
#         the full path of .spm file

#     Returns
#     -------
#     sou : array_like of string
#         source name of IVS designation
#     pm_ra/pm_dec : array_like of float
#         Proper motion in right ascension (*cos(Dec))/Declination in mas/yr.
#     pm_ra_err/pm_dec_err : array_like of float
#         Formal uncertainty of pm_ra/pm_dec in mas/yr.
#     corr : array_like of float
#         Correlation between pm_ra and pm_dec
#     """

#     sou = np.genfromtxt(ifile, usecols=(0,), dtype=str)
#     pm_ra, pm_dec, pm_ra_err, pm_dec_err, corr = np.genfromtxt(
#         ifile, usecols=range(1, 6), unpack=True)

#     return sou, pm_ra, pm_dec, pm_ra_err, pm_dec_err, corr


def read_sou(sou_file):
    """Read radio source positions

    Parameters
    ----------
    sou_file : string
        the full path of .sou file

    Return
    ------
    t_sou : astropy.table object
    """

    sol_label = "quasarpm-180510a"
    # Table of radio source position
    sou_file = "%s/%s.sou" % (sol_label, sol_label)
    print("  Remove these sources with 0 observation used in Solve.")
    os.system("~/Astronomy/Works/201711_GDR2_ICRF3/progs/sou_elim %s"
              % sou_file)

    t_sou = Table.read(sou_file, format="ascii.fixed_width_no_header",
                       names=("ivs_name", "ra", "ra_err", "dec", "dec_err",
                              "ra_dec_corr", "obs_used", "obs_tot", "ses_used",
                              "ses_tot", "date_beg", "date_end"),
                       col_starts=(10, 24, 45, 61, 82, 98, 117,
                                   132, 150, 164, 181, 202),
                       col_ends=(18, 41, 55, 78, 92, 104, 122,
                                 139, 155, 170, 191, 212))

    # convert string into float for RA, Decl. and observing epoch
    Nsize = len(t_sou)
    ra = np.empty(shape=(Nsize,), dtype=float)
    dec = np.empty(shape=(Nsize,), dtype=float)
    date_beg = np.empty(shape=(Nsize,), dtype=float)
    date_end = np.empty(shape=(Nsize,), dtype=float)

    for i, (rai, deci, epbi, epei) in enumerate(
            zip(t_sou["ra"], t_sou["dec"],
                t_sou["date_beg"], t_sou["date_end"])):
        ra[i] = RA_conv(rai)
        dec[i] = DC_conv(deci)
        date_beg[i] = date2mjd(epbi)
        date_end[i] = date2mjd(epei)

    # replace original columns with new columns
    t_sou["ra"] = ra
    t_sou["dec"] = dec
    t_sou["date_beg"] = date_beg
    t_sou["date_end"] = date_end

    # unit
    t_sou["ra"].unit = u.deg
    t_sou["dec"].unit = u.deg
    t_sou["ra_err"].unit = u.mas
    t_sou["dec_err"].unit = u.mas
    t_sou["date_beg"].unit = cds.MJD
    t_sou["date_end"].unit = cds.MJD

    return t_sou


def read_pmt(pmt_file):
    """Read the .spm file.

    Parameters
    ----------
    pmt_file: string
        the full path of .spm file

    Returns
    -------
    sou: array_like of string
        source name of IVS designation
    pm_ra/pm_dec: array_like of float
        Proper motion in right ascension(*cos(Dec))/Declination in mas/yr.
    pm_ra_err/pm_dec_err: array_like of float
        Formal uncertainty of pm_ra/pm_dec in mas/yr.
    corr: array_like of float
        Correlation between pm_ra and pm_dec
    """

    # Table of source proper motion
    t_pmt = Table.read(pmt_file, format="ascii",
                       names=["ivs_name", "pmra", "pmdec", "pmra_err",
                              "pmdec_err", "pmra_pmdec_corr"])

    # unit
    t_pmt["pmra"].unit = u.mas/u.year
    t_pmt["pmdec"].unit = u.mas/u.year
    t_pmt["pmra_err"].unit = u.mas/u.year
    t_pmt["pmdec_err"].unit = u.mas/u.year

    return t_pmt


def get_spm(spl_file):
    """Read source position and proper motions from spool file.

    Parameter
    ----------
    spl_file: string
        spool file
    """

    # Read .sou file to get source positions
    sou_file = "%s.sou" % spl_file[:-4]
    t_sou = read_sou(sou_file)

    # Extract proper infotmation from spool file
    pmt_file = "%s.pmt" % spl_file[:-4]
    retrieve_pmt(spl_file, pmt_file)
    t_pmt = read_pmt(pmt_file)

    # Table of source name
    t_name = Table.read("/Users/Neo/Astronomy/Tools/Calc_Solve/mk5-opa/"
                        "save_files/source.names",
                        format="ascii.fixed_width_no_header",
                        names=["ivs_name", "iers_name"],
                        col_starts=(0, 22),
                        col_ends=(8, 30))

    t_name = unique(t_name, keys="ivs_name")

    # Join these three tables
    t_tmp = join(t_sou, t_pmt, keys="ivs_name", join_type="left")
    t_spm = join(t_tmp, t_name, keys="ivs_name", join_type="left")

    # Multipliy the formal error in R.A. by a factor of cos(Decl.)
    factor = np.cos(Angle(t_spm["dec"]).radian)
    t_spm["ra_err"] = t_spm["ra_err"] * factor
    t_spm["pmra"] = t_spm["pmra"] * factor
    t_spm["pmra_err"] = t_spm["pmra_err"] * factor

    # Change the order of columns
    # Copy the columns
    iers_name = t_spm["iers_name"]
    pmra = t_spm["pmra"]
    pmdec = t_spm["pmdec"]
    pmra_err = t_spm["pmra_err"]
    pmdec_err = t_spm["pmdec_err"]
    pmra_pmdec_corr = t_spm["pmra_pmdec_corr"]

    # Remove the existing columns
    t_spm.remove_column("iers_name")
    t_spm.remove_column("pmra")
    t_spm.remove_column("pmdec")
    t_spm.remove_column("pmra_err")
    t_spm.remove_column("pmdec_err")
    t_spm.remove_column("pmra_pmdec_corr")

    # New columns
    date_mean = Column((t_spm["date_beg"] + t_spm["date_end"]) / 2.,
                       name="date_mean")

    # insert these columns
    t_spm.add_columns([iers_name, pmra, pmra_err, pmdec, pmdec_err,
                       pmra_pmdec_corr, date_mean],
                      [1, 6, 6, 6, 6, 6, 10])

    # Fill the masked(missing) value
    t_spm["iers_name"].fill_value = "-" * 8
    t_spm["pmra"].fill_value = 0.00
    t_spm["pmra_err"].fill_value = 999
    t_spm["pmdec"].fill_value = 0.00
    t_spm["pmdec_err"].fill_value = 999
    t_spm["pmra_pmdec_corr"].fill_value = -1.99

    t_spm["iers_name"] = t_spm["iers_name"].filled()
    t_spm["pmra"] = t_spm["pmra"].filled()
    t_spm["pmra_err"] = t_spm["pmra_err"].filled()
    t_spm["pmdec"] = t_spm["pmdec"].filled()
    t_spm["pmdec_err"] = t_spm["pmdec_err"].filled()
    t_spm["pmra_pmdec_corr"] = t_spm["pmra_pmdec_corr"].filled()

    # Add comments
    t_spm.meta["comments"] = [
        "VLBI Celestial Reference Frame Solution quasarpm-180510a",
        " Columns  Units   Meaning",
        "    1     --      IVS designation",
        "    2     --      IERS designation",
        "    3     deg     Right ascension",
        "    4     mas     Formal uncertainty of the right ascension "
        "(*cos(Dec))",
        "    5     deg     Declination",
        "    6     mas     Formal uncertainty of the declination",
        "    7     --      Correlation between right ascension and "
        "declination",
        "    8     mas/yr  Proper motion in right ascension",
        "    9     mas/yr  Formal uncertainty of proper motion in right "
        "ascension (*cos(Dec))",
        "   10     mas/yr  Proper motion in declination",
        "   11     mas/yr  Formal uncertainty of proper motion in declination",
        "   12     --      Correlation between proper motions in right "
        "ascension and declination",
        "   13     --      Number of delays",
        "   14     --      Number of sessions",
        "   15     days    Average epoch of observation (MJD)",
        "   16     days    First epoch of observation (MJD)",
        "   17     days    Last epoch of observation (MJD)",
        " Note: Some proper motion component were not estimated for some"
        " sources.",
        " They are however left in this file for sake of completeness.",
        " They can be identified by their sigma set to 999.000000 or "
        "206.264806 (default value set in Solve).",
        " Created date: %s." % time.strftime("%d/%m/%Y", time.localtime())]

    # Write table into .spm file
    spm_file = "%s.spm" % spl_file[:-4]
    t_spm.write(spm_file, format="ascii.fixed_width_no_header",
                delimiter="  ",
                exclude_names=["obs_tot", "ses_tot"],
                formats={"ivs_name": "%8s", "iers_name": "%8s",
                         "ra": "%14.10f", "dec": "%+14.10f",
                         "ra_err": "%10.4f", "dec_err": "%10.4f",
                         "ra_dec_corr": "%+7.4f",
                         "pmra": "%+10.4f", "pmdec": "%+10.4f",
                         "pmra_err": "%10.4f", "pmdec_err": "%10.4f",
                         "pmra_pmdec_corr": "%+7.4f"},
                overwrite=True)

    print("Write results into %s : Done!" % spm_file)


# -------------------------------- MAIN --------------------------------
if __name__ == '__main__':
    get_spm(sys.argv[1])
# --------------------------------- END --------------------------------
