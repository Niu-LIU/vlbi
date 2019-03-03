#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 11:01:07 2017

@author: Neo

Retrieve the estimates of positions of global stations and the formal
uncertainties of these estimates from .crl file which is generated by
the program getpar.

.crl file contains off-diagonal coefficients of correlations between
the estimates of EOP at the same experiment. Correlations are ordered in
according the order of elements of a symmetric matrix in low-triangular
representation without the main diagonal. Order of rows/columns: X_pole,
X_pole rate, Y_pole, Y_pole rate, UT1, UT1 rate, Nutation in longitude,
Nutation in obliquity

   File contains lines of two types:
1) Comment. The first character is #. Header comment contain the full name of
   the spool file.

2) Correlations.

   Field   Format Units     Meaning
   1-8     A8     --        record type identifier:  CRL_LOC:
   11-20   A10    --        database name with leading dollar sign
   23-25   I3     --        database version number
   29-34   F6.4   d/l       correlation between X_pole rate and X_pole
   36-41   F6.4   d/l       correlation between Y_pole and X_pole
   43-48   F6.4   d/l       correlation between Y_pole and X_pole rate
   50-55   F6.4   d/l       correlation between Y_pole rate and X_pole
   57-62   F6.4   d/l       correlation between Y_pole rate and X_pole rate
   64-69   F6.4   d/l       correlation between Y_pole rate and Y_pole
   71-76   F6.4   d/l       correlation between UT1 and X_pole
   78-83   F6.4   d/l       correlation between UT1 and X_pole rate
   85-90   F6.4   d/l       correlation between UT1 and Y_pole
   92-97   F6.4   d/l       correlation between UT1 and Y_pole rate
   99-104  F6.4   d/l       correlation between UT1 rate and X-pole
   106-111 F6.4   d/l       correlation between UT1 rate and X-pole rate
   113-118 F6.4   d/l       correlation between UT1 rate and Y-pole
   120-125 F6.4   d/l       correlation between UT1 rate and Y-pole rate
   127-132 F6.4   d/l       correlation between UT1 rate and UT1
   134-139 F6.4   d/l       correlation between Nutation Psi and X_pole
   141-146 F6.4   d/l       correlation between Nutation Psi and X_pole rate
   148-153 F6.4   d/l       correlation between Nutation Psi and Y_pole
   155-160 F6.4   d/l       correlation between Nutation Psi and Y_pole rate
   162-167 F6.4   d/l       correlation between Nutation Psi and UT1
   169-174 F6.4   d/l       correlation between Nutation Psi and UT1 rate
   176-181 F6.4   d/l       correlation between Nutation Eps and X_pole
   183-188 F6.4   d/l       correlation between Nutation Eps and X_pole rate
   190-195 F6.4   d/l       correlation between Nutation Eps and Y_pole
   197-202 F6.4   d/l       correlation between Nutation Eps and Y_pole rate
   204-209 F6.4   d/l       correlation between Nutation Eps and UT1
   211-216 F6.4   d/l       correlation between Nutation Eps and UT1 rate
   218-223 F6.4   d/l       correlation between Nutation Eps and Nutation Psi

"""

import numpy as np


# ------------------------------  FUNCTIONS  ---------------------------
def read_crl(datafile):
    '''Retrieve the result from .crl file.

    Parameters
    ----------
    datafile : string
        name of data file

    Returns
    ----------
    dbname : array, string
       database name with leading dollar sign
    corXXR : array, float
        correlation between X and XR
    corXY : array, float
        correlation between X and Y
    corXYR : array, float
        correlation between X and YR
    corXU : array, float
        correlation between X and U
    corXUR : array, float
        correlation between X and UR
    corXP : array, float
        correlation between X and P
    corXE : array, float
        correlation between X and E
    corXRY : array, float
        correlation between XR and Y
    corXRYR : array, float
        correlation between XR and YR
    corXRU : array, float
        correlation between XR and U
    corXRUR : array, float
        correlation between XR and UR
    corXRP : array, float
        correlation between XR and P
    corXRE : array, float
        correlation between XR and E
    corYYR : array, float
        correlation between Y and YR
    corYU : array, float
        correlation between Y and U
    corYUR : array, float
        correlation between Y and UR
    corYP : array, float
        correlation between Y and P
    corYE : array, float
        correlation between Y and E
    corYRU : array, float
        correlation between YR and U
    corYRUR : array, float
        correlation between YR and UR
    corYRP : array, float
        correlation between YR and P
    corYRE : array, float
        correlation between YR and E
    corUUR : array, float
        correlation between U and UR
    corUP : array, float
        correlation between U and P
    corURE : array, float
        correlation between UR and E
    corURP : array, float
        correlation between UR and P
    corPE : array, float
        correlation between P and E
    '''

    dbname = np.genfromtxt(datafile, dtype=str, usecols=(1,))
    [corXXR,
     corXY, corXRY,
     corXYR, corXRYR, corYYR,
     corXU, corXRU, corYU, corYRU,
     corXUR, corXRUR, corYUR, corYRUR, corUUR,
     corXP, corXRP, corYP, corYRP, corUP, corURP,
     corXE, corXRE, corYE, corYRE, corUE, corURE, corPE] = np.genfromtxt(
        datafile, usecols=np.arange(3, 31), unpack=True,
        missing_values='*'*8,
        filling_values=0.)

    return [dbname,
            corXX,
            corXY, corXRY,
            corXYR, corXRYR, corYYR,
            corXU, corXRU, corYU, corYRU,
            corXUR, corXRUR, corYUR, corYRUR, corUUR,
            corXP, corXRP, corYP, corYRP, corUP, corURP,
            corXE, corXRE, corYE, corYRE, corUE, corURE, corPE]

# # ------------------------------  MAIN BODY  ---------------------------
# # Retrieve estimates.
# [dbname,
#  corXXR,
#  corXY, corXRY,
#  corXYR, corXRYR, corYYR,
#  corXU, corXRU, corYU, corYRU,
#  corXUR, corXRUR, corYUR, corYRUR, corUUR,
#  corXP, corXRP, corYP, corYRP, corUP, corURP,
#  corXE, corXRE, corYE, corYRE, corUE, corURE, corPE] = read_crl(
#     'result/test.crl')
#
# print(dbname[0])
# print(corXXR[0],
#       corXY[0], corXRY[0],
#       corXYR[0], corXRYR[0], corYYR[0],
#       corXU[0], corXRU[0], corYU[0], corYRU[0],
#       corXUR[0], corXRUR[0], corYUR[0], corYRUR[0], corUUR[0],
#       corXP[0], corXRP[0], corYP[0], corYRP[0], corUP[0], corURP[0],
#       corXE[0], corXRE[0], corYE[0], corYRE[0], corUE[0], corURE[0], corPE[0])
# ------------------------------ END -----------------------------------
