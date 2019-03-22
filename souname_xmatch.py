#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: souname_xmatch.py
"""
Created on Thu Apr 26 15:55:00 2018

@author: Neo(liuniu@smail.nju.edu.cn)

Read radio source name from data list. Obsolete!
"""

from astropy.io import fits
import numpy as np
import os
import time
from cross_match import list_crossmatch


# -----------------------------  FUNCTIONS -----------------------------
def rewrite_srcnames():
    '''Read source names.

    Parameters
    ----------
    None

    Returns
    ----------
    None
    '''

    ifile = open("aux_files/IVS_SrcNamesTable.txt", "r")
    ofile = open("aux_files/IVS_SrcNamesTable.csv", "w")

    print("#==================================\n"
          "# IVS source Name Translation Table\n"
          "#==================================\n"
          " # File format      : 26/04/2018\n"
          " # Date last revised: %s\n"
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
        if not len(line) or line.startswith('#'):
            continue

        if line[40] is "-" or line[40] is " ":
            print("%s,%s,%s" % (line[:8], line[10:26], line[:8]),
                  file=ofile)
        else:
            print("%s,%s,%s" % (line[:8], line[10:26], line[40:48]),
                  file=ofile)

    ofile.close()


def rewrite_sourcename(datafile=None):
    '''Read source names from a file named 'source.names' and rewrite it.

    Parameters
    ----------
    datafile : str
        full path of source.names. if None, default value will be used.

    Returns
    -------
    IVS : array of string
        IVS source name
    ICRF : array of string
        J2000 source name long
    IERS : array of string
        IERS designations source name
    '''

    if datafile is None:
        datafile = ("aux_files/source.names")

    # print(datafile)
    ifile = open(datafile, "r")

    ivs_name = []
    iers_name = []
    Typ = []

    for line in ifile.readlines():
        if not len(line) or line.startswith('#') or line.startswith('\n'):
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
    # icrf_name = np.empty_like(ivs_name, dtype='S16')

    ofile = open("aux_files/source.names.csv", "w")

    print("#==================================\n"
          "# IVS source Name Translation Table\n"
          "#==================================\n"
          "  # File format      : 26/04/2018\n"
          "  # Date last revised: %s\n"
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

    for ivs_namei, iers_namei in zip(ivs_name, iers_name):
        print("%s,%s,%s" % (ivs_namei, " "*16, iers_namei),
              file=ofile)

    ofile.close()


def read_souname():
    '''Read source names.

    Parameters
    ----------
    Nonec

    Returns
    -------
    IVS : array of string
        IVS source name
    ICRF : array of string
        J2000 source name long
    IERS : array of string
        IERS designations source name
    '''

    datafile = "aux_files/source.names.csv"

    if not os.path.exists(datafile):
        print("Couldn't find the file" % datafile)

    # empty array to store data
    IVS, ICRF, IERS = np.genfromtxt(datafile, dtype=str,
                                    delimiter=",", unpack=True)
    return IVS, ICRF, IERS


def find_sou_designation(sou, nametype):
    """Find the other corresponding radio source designations.

    Parameters
    ----------
    sou : array_like
        radio source name
    nametype : string
        type of input source name, could be set as "IVS", "ICRF", or
        "IERS"

    Returns
    -------
    list1, list2, list3 : ndarrays
        IVS, ICRF, and IERS designation of input sources.
    """

    IVS, ICRF, IERS = read_souname()

    if nametype is "IVS" or "ivs":
        list0 = IVS
    elif nametype is "ICRF" or "icrf":
        list0 = ICRF
    elif nametype is "IERS" or "iers":
        list0 = IERS
    else:
        print("ERROR! Nametype can only be IVS, ICRF, or IERS!")
        exit()

    list1 = np.empty_like(sou)
    list2 = np.empty_like(sou)
    list3 = np.empty_like(sou)

    for i, soui in enumerate(sou):
        indarr = np.where(list0 == soui)[0]

        if indarr.size:
            j = indarr[0]
            list1[i] = IVS[j]
            list2[i] = ICRF[j]
            list3[i] = IERS[j]
        else:
            if nametype is "IVS" or "ivs":
                list1[i] = soui
                list2[i] = " " * 16
                list3[i] = " " * 8
            elif nametype is "ICRF" or "icrf":
                list1[i] = " " * 8
                list2[i] = soui
                list3[i] = " " * 8
            elif nametype is "IERS" or "iers":
                list1[i] = " " * 8
                list2[i] = " " * 16
                list3[i] = soui

    return list1, list2, list3
# --------------------------------- END --------------------------------
