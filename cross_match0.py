#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: list_crossmatch.py
"""
Created on Fri Apr 27 00:24:29 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np
from numpy import cos, sqrt
from functools import reduce


__all__ = ["list_crossmatch", "pos_Xmatch"
           "nor_sep_calc", "postional_difference_calc", "pos_max_calc"]


# -----------------------------  FUNCTIONS -----------------------------
def list_crossmatch(X1, X2):
    '''Corssmatch between two list.

    Parameters
    ----------
    X1 : array_like
        first dataset, shape(N1, D)
    X2 : array_like
        second dataset, shape(N2, D)

    Returns
    -------
    '''

    # =======
    # Add some codes here to eliminate duplicate elemnets
    # =======

    com_list = []
    index1 = []
    index2 = []

    for i, x1 in enumerate(X1):
        indarr = np.where(X2 == x1)[0]

        if indarr:
            com_list.append(x1)
            index1.append(i)
            j = indarr[0]
            index2.append(j)

    com_list = np.asarray(com_list, dtype=str)
    index1 = np.asarray(index1, dtype=int)
    index2 = np.asarray(index2, dtype=int)

    return com_list, index1, index2


def position_taken(index, RA, RAc_err, Dec, Dec_err, cor):
    '''Extract the elements from array at specific index.

    Parameters
    ----------
    index :
        the indice corresponding to common sources
    RA / Dec :
        Right Ascension / Declination for all sources, degreees
    RAc_err / Dec_err :
        formal uncertainty of RA*cos(Dec) / Dec for all sources, micro-as.
    cor :
        correlation coeffient between RA and Dec for all sources.

    Returns
    ----------
    RAn / Decn :
        Right Ascension / Declination for common sources, degreees
    RAc_errn / Dec_err :
        formal uncertainty of RA*cos(Dec) / Dec for common sources, micro-as.
    corn :
        correlation coeffient between RA and Dec for common sources.
    '''

    RAn = np.take(RA, index)
    RAc_errn = np.take(RAc_err, index)
    Decn = np.take(Dec, index)
    Dec_errn = np.take(Dec_err, index)
    corn = np.take(cor, index)

    return RAn, RAc_errn, Decn, Dec_errn, corn


def pos_Xmatch(sou1, RA1, RAc_err1, Dec1, Dec_err1, cor1,
               sou2, RA2, RAc_err2, Dec2, Dec_err2, cor2):
    '''Crossmatch between Gaia and VLBI catalogs.

    Parameters
    ----------
    sou :
        source name (ICRF designation)
    RA / Dec :
        Right Ascension / Declination, degreees
    RAc_err / DC_err :
        formal uncertainty of RA*cos(Dec) / Dec, micro-as.
    cor :
        correlation coeffient between RA and Dec.

    Returns
    ----------
    soucom :
        name (ICRF designation) of common sources
    RAn / Decn :
        Right Ascension / Declination for common sources, degreees
    RAc_errn / Dec_err :
        formal uncertainty of RA*cos(Dec) / Dec for common sources, micro-as.
    cor :
        correlation coeffient between RA and Dec for common sources.
    '''

    soucom = []
    index1 = []
    index2 = []

    for i, soui in enumerate(sou1):
        indarr = np.where(sou2 == soui)[0]

        if indarr:
            soucom.append(soui)
            index1.append(i)
            j = indarr[0]
            index2.append(j)

    RA1n, RAc_err1n, Dec1n, Dec_err1n, cor1n = position_taken(
        index1, RA1, RAc_err1, Dec1, Dec_err1, cor1)
    RA2n, RAc_err2n, Dec2n, Dec_err2n, cor2n = position_taken(
        index2, RA2, RAc_err2, Dec2, Dec_err2, cor2)

    # list -> array
    soucom = np.asarray(soucom)

    return [soucom,
            RA1n, RAc_err1n, Dec1n, Dec_err1n, cor1n,
            RA2n, RAc_err2n, Dec2n, Dec_err2n, cor2n]


def nor_sep_calc(dRA, dRA_err, dDC, dDC_err, C):
    '''Calculate the normalized seperation.

    Parameters
    ----------
    dRA/dDC : Right Ascension / Declination differences in micro-as
    dRA_err/dDC_err : formal uncertainty of dRA*cos(Dec)/dDC in micro-as
    C : correlation coeffient between dRA*cos(Dec) and dDC.

    Returns
    ----------
    ang_sep : angular seperation, in micro-as
    X_a / X_d : normalized coordinate differences in RA / DC, unit-less
    X : Normalized separations, unit-less.
    '''

    # Angular seperations
    ang_sep = sqrt(dRA**2 + dDC**2)

    # Normalised coordinate difference
    X_a = dRA / dRA_err
    X_d = dDC / dDC_err

    # Normalised separation - Mignard's statistics (considering covariance)
    X = np.zeros_like(X_a)

    for i, (X_ai, X_di, Ci) in enumerate(zip(X_a, X_d, C)):
        if Ci == -1.:
            Ci = -0.99
        if Ci == 1.0:
            Ci = 0.99
        # print(Ci)

        wgt = np.linalg.inv(np.mat([[1, Ci],
                                    [Ci, 1]]))
        Xmat = np.mat([X_ai, X_di])
        X[i] = sqrt(reduce(np.dot, (Xmat, wgt, Xmat.T)))

    # # Normalised separation (normal way)
    X2 = sqrt(X_a**2 + X_d**2)

    # return ang_sep, X_a, X_d, X
    return ang_sep, X_a, X_d, X, X2


def postional_difference_calc(RA1, RA1_err, DC1, DC1_err, Cor1,
                              RA2, RA2_err, DC2, DC2_err, Cor2,
                              arccof=None):
    '''Calculate the normalized seperation between VLBI and Gaia positions.


    Parameters
    ----------
    RA / DC : Right Ascension / Declination, degress
    e_RA / e_DC : formal uncertainty of RA * cos(Dec) / DC, mas
    Cor : correlation coeffient between RA and DC.
    arccof : cos(Dec.)

    Note: suffix 'G' stands for GaiaDR1 and I for VLBI catalog.

    Returns
    ----------
    ang_sep : angular seperation in micro-as
    X_a / X_d : normalized seperation in RA / DC, unit-less
    X : Normalized separations, unit-less.
    '''

    if arccof is None:
        arccof = np.cos(np.deg2rad(DC1))

    # # deg -> uas
    # dRA = (RA1 - RA2) * 3.6e9 * arccof
    # dRA_err = sqrt(RA1_err**2 + RA2_err**2)
    # dDC = (DC1 - DC2) * 3.6e9
    # dDC_err = sqrt(DC1_err**2 + DC2_err**2)

    # deg -> mas
    dRA = (RA1 - RA2) * 3.6e6 * arccof
    dRA_err = sqrt(RA1_err**2 + RA2_err**2)
    dDC = (DC1 - DC2) * 3.6e6
    dDC_err = sqrt(DC1_err**2 + DC2_err**2)

    # Correlation coefficient of combined errors
    cov = RA1_err * DC1_err * Cor1 + RA2_err * DC2_err * Cor2
    corf = cov / (dRA_err * dDC_err)

    # Normalised separation
    ang_sep, X_a, X_d, X, X2 = nor_sep_calc(dRA, dRA_err,
                                            dDC, dDC_err, corf)

    # return ang_sep, X_a, X_d, X
    return dRA, dDC, dRA_err, dDC_err, cov, ang_sep, X_a, X_d, X, X2


def pos_max_calc(ra_err, dec_err, ra_dec_corr):
    """Calculate the semi-major axis of the dispersion ellipse.

    Parameters
    ----------
    ra_err/dec_err : formal uncertainty of RA/Dec, usually in micro-as
    ra_dec_corr : correlation coeffient between RA and Dec, unitless.

    Returns
    ----------
    sig_pos_max : semi-major axis of the dispersion ellipse for
                  characterising the positional uncertainty of a source;
                  same unit as ra_err/dec_err
    """

    sig_pos_max = sqrt(0.5 * (ra_err**2 + dec_err**2 +
                              sqrt((ra_err**2 - dec_err**2)**2 +
                                   (2*ra_err*dec_err*ra_dec_corr)**2)))

    return sig_pos_max


def overall_err_calc(ra_err, dec_err, ra_dec_corr):
    """Calculate the ovrall formal uncertainty.

    ovrall formal uncertainty = sqrt(RA_err^2+Dec_err^2+C*RA_err*Dec_err)

    Parameters
    ----------
    ra_err/dec_err : formal uncertainty of RA/Dec, usually in micro-as
    ra_dec_corr : correlation coeffient between RA and Dec, unitless.

    Returns
    ----------
    overall_err : ovrall formal uncertainty;
                  same unit as ra_err/dec_err
    """

    overall_err = sqrt(ra_err**2 + dec_err**2 + ra_err*dec_err*ra_dec_corr)

    return overall_err
# --------------------------------- END --------------------------------
