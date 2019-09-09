#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: read_arc.py
"""
Created on Wed Aug 28 09:19:53 2019

@author: Neo(liuniu@smail.nju.edu.cn)
"""

from astropy.time import Time
import numpy as np
import re

# My modules
from .get_dir import get_data_dir

__all__ = ["get_ns_code", "sta_obs_hist"]


# -----------------------------  FUNCTIONS -----------------------------
def get_sess_and_ns(arcfile):
    """Extract session identifier and network from arc file.

    This maybe only works for the src file of Sebastien.

    Parameter
    ---------
    arcfile : string
        name and path of the arc file

    Return
    ------
    sessIDs : 1d array
        session identifiers
    networks : string of 1d array
        network consisting of network station code
    """

    sessIDs = []
    networks = []

    with open(arcfile, "r") as farc:
        for line in farc.readlines():

            line = line.strip().strip("\n")

            # The commented character is "*"
            if line.startswith("*"):
                continue

            line1, line2 = line.split("!")
            sessID = line1.split()[0]
            network = line2.split()[-1]

            sessIDs.append(sessID)
            networks.append(network)

    return sessIDs, networks


def sess_epoch_calc(sessID):
    """Calculate the epoch from session identifiers.

    Parameter
    ---------
    sessID : string
        session identifier

    Return
    ------
    epoch : 1d-array
        observing epoch in Julian year
    """

    # Mapping table between letter and number for month
    mstr2num = {
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

    year = int(sessID[1:3])
    monstr = sessID[3:6]
    date = sessID[6:8]

    if year >= 79:
        year += 1900
    else:
        year += 2000

    epoch = Time("{}-{}-{}".format(year, mstr2num[monstr], date), format="iso")

    return epoch.jyear


def get_ns_code(nsfile=None):
    """

    Parameter
    ---------
    nsfile : string
        name and direcotry of the network station code file

    Return
    ------
    ns_code_name : dict
        mapping table of network station and code
    """

    # network station code file
    if nsfile is None:
        datadir = get_data_dir()
        nsfile = "{}/ns-codes.txt".format(datadir)

    # Initialize the dict
    ns_code_name = {}

    with open(nsfile, "r") as nsfile:
        for line in nsfile.readlines():
            # The commented character is "*"
            if line.startswith("*"):
                continue

            stacode = line[1:3]
            staname = line[4:12].strip()

            ns_code_name[stacode] = staname

    return ns_code_name


def sta_obs_hist(arcfile):
    """Get the observing history of each station

    Parameter
    ---------
    arcfile : string
        name and path of the arc file

    Return
    ------
    staobs : dict object
        key is station name and value is an array of observing epoch
    """

    # Parse the arc file
    sessids, networks = get_sess_and_ns(arcfile)
    sess_epos = list(map(sess_epoch_calc, sessids))

    sta_code_name = get_ns_code()

    # Initialize an empty dict to store the observational epoch
    staobs = {}

    for code in sta_code_name.values():
        staobs[code] = []

    # A loop to parse each session
    for epochi, networki in zip(sess_epos, networks):

        # I do not know what "00" means
        if networki == "00":
            continue

        # Separate the network
        # All VLBA stations
        if networki == "Va":
            codes = re.findall(".{2}", "BrFdHnKpLaMkNlPtOvSc")
        else:
            codes = re.findall(".{2}", networki)

        for code in codes:
            sta_name = sta_code_name[code]
            staobs[sta_name].append(epochi)

    # Delete stations with no observation
    delkeys = []

    for key, value in staobs.items():
        if len(value):
            continue
        delkeys.append(key)

    sta_obs_all = staobs.copy()
    for key in delkeys:
        del staobs[key]

    return staobs


def sess_sta_count(arcfile):
    """Count the number of participated stations for VLBI sessions.

    Parameter
    ---------
    arcfile : string
        name and path of the arc file

    Returns
    -------
    obsepoch : 1d array
        observing epoch
    stanum : 1d-array
        station number
    """

    # Parse the arc file
    sessids, networks = get_sess_and_ns(arcfile)
    obsepoch = list(map(sess_epoch_calc, sessids))

    # Initialize an empty dict to store the observational epoch
    stanum = []

    # A loop to parse each session
    for networki in networks:
        # I do not know what "00" means
        if networki == "00":
            stanum.append(0)

        # Separate the network
        # All VLBA stations
        if networki == "Va":
            stanum.append(10)
        else:
            codes = re.findall(".{2}", networki)
            stanum.append(len(codes))

    return staobs

# --------------------------------- END --------------------------------
