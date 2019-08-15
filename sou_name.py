#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: souname_xmatch.py
"""
Created on Thu Apr 26 15:55:00 2018

@author: Neo(liuniu@smail.nju.edu.cn)

Read radio source name from data list. Obsolete!
"""

from astropy.io import fits
from astropy import table
from astropy.table import Table
import numpy as np
import os
import time

# My modules
from .get_dir import get_data_dir


__all__ = ["get_souname"]


# -----------------------------  FUNCTIONS -----------------------------
def rewrite_srcnames():
    """Read source names.

    Parameters
    ----------
    None

    Returns
    ----------
    None
    """

    datadir = get_datadir()

    ifile = open("{}/IVS_SrcNamesTable.txt".format(datadir), "r")
    ofile = open("{}/IVS_SrcNamesTable.csv".format(datadir), "w")

    print("#==================================\n"
          "# IVS source Name Translation Table\n"
          "#==================================\n"
          "# File format      : 26/04/2018\n"
          "# Date last revised: %s\n"
          "# Columns      Format             Content\n"
          "#   1          A8                 IVS source name\n"
          "#                                 Source name used "
          "now or in the past in IVS schedules and databases\n"
          "#   2          A16                J2000 source name long\n"
          "#                                 ICRF designations following"
          " therecommendations of the IAU Task Group on Designations\n"
          "#   3          A8                 IERS designations source name"
          " constructed from B1950 coordinates\n#"
          % time.strftime("%d/%m/%Y", time.localtime()), file=ofile)

    for line in ifile.readlines():
        if not len(line) or line.startswith("#"):
            continue

        if line[40] is "-" or line[40] is " ":
            print("%s,%s,%s" % (line[:8], line[10:26], line[:8]),
                  file=ofile)
        else:
            print("%s,%s,%s" % (line[:8], line[10:26], line[40:48]),
                  file=ofile)

    ofile.close()


def rewrite_sourcename(datafil=None):
    """Read source names from a file named "source.names" and rewrite it.

    Parameters
    ----------
    datafil : str
        full path of source.names. if None, default value will be used.

    """

    if datafil is None:
        datadir = get_datadir()
        datafil = "{}/source.names".format(datadir)

    # print(datafil)
    ifile = open(datafil, "r")

    ivs_name = []
    iers_name = []
    Typ = []

    for line in ifile.readlines():
        if not len(line) or line.startswith("#") or line.startswith("\n"):
            continue

        # print(line[0], len(line))
        ivs_name.append(line[:8].rstrip())
        iers_name.append(line[22:30])
        Typ.append(line[42])

    ifile.close()

    # convert list to array
    ivs_name = np.asarray(ivs_name)
    iers_name = np.asarray(iers_name)
    Typ = np.asarray(Typ)
    # icrf_name = np.empty_like(ivs_name, dtype="S16")

    ofile = open("{}/source.names.csv".format(datadir), "w")

    print("#==================================\n"
          "# IVS source Name Translation Table\n"
          "#==================================\n"
          "  # File format      : 26/04/2018\n"
          "  # Date last revised: %s\n"
          "# Columns      Format             Content\n"
          "#   1          A8                 IVS source name\n"
          "#                                 Source name used "
          "now or in the past in IVS schedules and databases\n"
          "#   2          A8                 IERS designations source name"
          " constructed from B1950 coordinates\n#"
          % time.strftime("%d/%m/%Y", time.localtime()), file=ofile)

    for ivs_namei, iers_namei in zip(ivs_name, iers_name):
        print("%s,%s" % (ivs_namei, iers_namei),
              file=ofile)

    ofile.close()


def get_souname():
    """Read source names.

    Parameters
    ----------
    Nonec

    Returns
    -------
    souname: table object
        different designations of radio source names
    """

    datadir = get_datadir()
    datafil = "{}/source.names.csv".format(datadir)

    rewrite_sourcename()

    # empty array to store data
    souname = Table.read(datafil, format="ascii.csv",
                         names=["ivs_name", "iers_name"])

    # Eliminate duplicate sources
    souname = table.unique(souname, keys="ivs_name")
    souname["iers_name"].fill_value = " " * 8

    return souname


if __name__ == "__main__":
    get_souname()
# --------------------------------- END --------------------------------
