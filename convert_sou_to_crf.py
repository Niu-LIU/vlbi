#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: convert_sou_to_crf.py
"""
Created on Thu Apr 26 15:46:01 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""

from astropy.table import Table, join, unique, Column
import numpy as np
import sys
import time

# My module
from my_progs.vlbi.read_sou import read_sou

__all__ = ["write_into_crf", "convert_sou_to_crf"]


# -----------------------------  FUNCTIONS -----------------------------
def write_into_crf(t_sou, crf_file, sln_label=""):
    """
    Parameters
    ----------
    t_sou : astropy.table object
    crf_file : string
        the full path of .sou file

    """

    # Table of source name
    t_srcname = Table.read("aux_files/source.names",
                           format="ascii.fixed_width_no_header",
                           names=["ivs_name", "j2000_name", "iers_name"],
                           col_starts=(0, 10, 22),
                           col_ends=(8, 20, 30))

    t_srcname = unique(t_srcname, keys="ivs_name")

    # Join these three tables
    t_sou = join(t_sou, t_srcname, keys="ivs_name", join_type="left")

    # New columns
    mean_date = Column((t_sou["beg_date"] + t_sou["end_date"]) / 2.,
                       name="mean_date")

    # Change the order of columns
    # Copy the columns
    iers_name = t_sou["iers_name"]
    # Remove the existing columns
    t_sou.remove_column("iers_name")

    # Insert these columns
    t_sou.add_columns([iers_name, mean_date], [1, 10])

    # Fill the masked(missing) value
    t_sou["iers_name"] = t_sou["iers_name"].filled("-" * 8)

    # Add comments
    t_sou.meta["comments"] = [
        "VLBI Celestial Reference Frame Solution %s" % sln_label,
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
        "    8     --      Number of delays",
        "    9     --      Number of sessions",
        "   10     days    Average epoch of observation (MJD)",
        "   11     days    First epoch of observation (MJD)",
        "   12     days    Last epoch of observation (MJD)",
        " Created on %s." % time.strftime("%d/%m/%Y", time.localtime()),
        " "]

    t_sou.write(crf_file, format="ascii.fixed_width_no_header",
                delimiter="",
                exclude_names=["j2000_name", "total_obs",
                               "total_sess", "pos_err"],
                formats={"ivs_name": "%-8s", "iers_name": "%-8s",
                         "ra": "%14.10f", "dec": "%+14.10f",
                         "ra_err": "%10.4f", "dec_err": "%10.4f",
                         "ra_dec_corr": "%+7.4f"}, overwrite=True)


def convert_sou_to_crf(sou_file, sln_label=""):
    """Rewrite the source position into .crf file

    Parameters
    ----------
    sou_file : string
        the full path of .sou file

    Return
    ------
    t_sou : astropy.table object
    """

    # Read .sou file to get source positions
    t_sou = read_sou(sou_file)

    # write the source position into a .cat file
    crf_file = "%s.crf" % sou_file[:-4]

    write_into_crf(t_sou, crf_file, sln_label)
    print("%s \n --> \n%s \n Done!" % (sou_file, crf_file))


# -------------------------------- MAIN --------------------------------
if __name__ == "__main__":
    if len(sys.argv) == 2:
        convert_sou_to_crf(sys.argv[1])
    elif len(sys.argv) == 3:
        convert_sou_to_crf(sys.argv[1], sys.argv[2])
    else:
        print("Too many arguments!")
# -------------------------------- END --------------------------------
