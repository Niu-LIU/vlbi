#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: read_trf.py
"""
Created on Tue Mar 26 21:57:42 2019

@author: Neo(liuniu@smail.nju.edu.cn)

Read the estimate of the station position and velocity

"""

from astropy.table import join
import numpy as np
import os
import sys
from .read_sta import read_sta
from .read_vel import read_vel

__all__ = ["read_trf"]


# -----------------------------  FUNCTIONS -----------------------------
def read_trf(slnlabel):
    """Read the estimate of the station position and velocity

    Parameters
    ----------
    slnlabel : string
        solution label

    Returns
    -------
    t_crf : astropy.table object
    """

    # Full path of the .sta and .vel file
    stafile = "%s.sta" % slnlabel
    velfile = "%s.vel" % slnlabel

    # position information
    sta = read_sta(stafile)

    # velocity information
    vel = read_vel(velfile, get_corr=False)

    # Merge the two tables into one
    trf = join(sta, vel, keys="station", join_type="outer")

    return trf


# -------------------------------- MAIN --------------------------------
if __name__ == "__main__":
    if len(sys.argv) == 2:
        read_trf(sys.argv[1])
# --------------------------------- END --------------------------------
