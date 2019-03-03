#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: read_itrf.py
"""
Created on Thu Jan 25 15:18:36 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np
from numpy import genfromtxt


__all__ = {"read_itrf", "find_stable_sta"}


# -----------------------------  FUNCTIONS -----------------------------
def read_itrf(datafile):
    '''Retrieve site positions/voelocities from ITRF file.

    Parameters
    ---------
    datafile : string
            path and name of ITRF file.

    Returns
    ---------
    sta : array, string
            Station name
    pvX : array, float
            Position/Velocity of X-component, meter
    pvY : array, float
            Position/Velocity of Y-component, meter
    pvZ : array, float
            Position/Velocity of Z-component, meter
    '''

    sta = []
    pvX = []
    pvY = []
    pvZ = []

    for line in open(datafile, 'r'):
        if line[0] != '$':
            sta.append(line[:12].strip())
            # pvX.append(float(line[12:27]))
            # pvY.append(float(line[27:43]))
            # pvZ.append(float(line[43:58]))

            pvX.append(float(line[12:28]))
            pvY.append(float(line[28:45]))
            pvZ.append(float(line[45:61]))

    sta = np.array(sta)
    pvX = np.array(pvX)
    pvY = np.array(pvY)
    pvZ = np.array(pvZ)

    return sta, pvX, pvY, pvZ


def find_stable_sta(sta, pvX, pvY, pvZ):
    '''Get the position/velocity for stable stations.

    'Stable' here means the stations which position is unique in ITRF.
    '''

    cond = np.ones_like(sta, dtype=bool)

    for i, stai in enumerate(sta):

        if cond[i]:

            ind = np.where(sta == stai)[0]
            if ind.size > 1:
                for j in ind:
                    cond[j] = False

    nsta = sta[cond]
    npvX = pvX[cond]
    npvY = pvY[cond]
    npvZ = pvZ[cond]

    return nsta, npvX, npvY, npvZ


def test_code():
    '''Test these codes.
    '''

    sta, pvX, pvY, pvZ = read_itrf(
        # "/home/oper/traitement/itrf2014.sit")  # vlbi2
        "/Users/Neo/Astronomy/Data/SOLVE/itrf2014_sitmod")  # My Mac

    nsta, npvX, npvY, npvZ = find_stable_sta(sta, pvX, pvY, pvZ)

    print(nsta, npvX, npvY, npvZ)

# --------------------------------- TEST ------------------------------
# test_code()
# --------------------------------- END -------------------------------
