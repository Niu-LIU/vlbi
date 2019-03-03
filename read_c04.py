#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: read_c04.py
"""
Created on Tue Dec 12 16:10:27 2017

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np


# -----------------------------  FUNCTIONS -----------------------------
def read_c04(C04_file):
    '''Fetch the C04 series.


    Parameters
    ----------
    C04_file : string
                path and name of C04 data file.

    Returns
    ----------
    mjd : array, float
        epoch in modified Julian date
    Xp : array, float
        Xp position of CIP in ITRS, mas
    Yp : array, float
        Yp position of CIP in ITRS, mas
    U : array, float
        UT1 - UTC, ms
    # LOD : array, float
    #     length of day, ms
    dX : array, float
        Xp component of CPO, mas
    dY : array, float
        Yp component of CPO, mas
    XpErr : array, float
        formal uncertainty of Xp, mas
    YpErr : array, float
        formal uncertainty of Yp, mas
    UErr : array, float
        formal uncertainty of U, ms
    # LODErr : array, float
    #     formal uncertainty of LOD, ms
    dXErr : array, float
        formal uncertainty of dX, mas
    dYErr : array, float
        formal uncertainty of dY, mas
    '''

    mjd, Xp, Yp, U = np.genfromtxt(C04_file, skip_header=14,
                                   usecols=np.arange(3, 7), unpack=True)
    dX, dY = np.genfromtxt(C04_file, skip_header=14,
                           usecols=(8, 9), unpack=True)
    XpErr, YpErr, UErr = np.genfromtxt(C04_file, skip_header=14,
                                       usecols=np.arange(10, 13),
                                       unpack=True)
    dXErr, dYErr = np.genfromtxt(C04_file, skip_header=14,
                                 usecols=(14, 15),
                                 unpack=True)

    # arc-sec -->  mas, second --> ms
    Xp, Yp, U = Xp * 1000, Yp * 1000, U * 1000
    XpErr, YpErr, UErr = XpErr * 1000, YpErr * 1000, UErr * 1000

    dX, dY = dX * 1000, dY * 1000
    dXErr, dYErr = dXErr * 1000, dYErr * 1000

    # print(mjd[0], Xp[0], Yp[0], U[0], LOD[0], dX[0], dY[0],
    #       XpErr[0], YpErr[0], UErr[0], LODErr[0],
    #       dXErr[0], dYErr[0])

    # return [mjd, Xp, Yp, U, LOD, dX, dY,
    #         XpErr, YpErr, UErr, LODErr,
    #         dXErr, dYErr]

    return mjd, Xp, Yp, U, dX, dY, XpErr, YpErr, UErr, dXErr, dYErr

# --------------------------------- END --------------------------------
