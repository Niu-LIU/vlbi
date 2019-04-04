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
def read_crf(sln_label):
    """Read the estimate of the station position and velocity

    Parameters
    ----------
    sln_lable : string
        solution label

    Returns
    -------
    t_crf : astropy.table object
    """

    # Full path of the .sta and .vel file
    sta_file = "%s.sta" % sln_label
    vel_file = "%s.vel" % sln_lable

    # position information
    t_sta = read_sta(sta_file)

    # velocity information
    t_vel = read_vel(vel_file)

    # Merge the two tables into one
    t_crf = join(t_sta, t_vel, keys="station", join_type="outer")

    return crf


# -------------------------------- MAIN --------------------------------
if __name__ == "__main__":
    if len(sys.argv) == 2:
        read_crf(sys.argv[1])
# --------------------------------- END --------------------------------
