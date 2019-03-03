#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: write_solvesrc.py
"""
Created on Mon Jan 22 15:04:40 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np

__all__ = ['write_solve_src', 'write_nnrs']


# -----------------------------  FUNCTIONS -----------------------------
def deg2dms(angle):
    '''RA(deg) -> D, M, S
    '''

    sgn = np.sign(angle)
    angle = np.fabs(angle)

    D = angle // 1
    M = angle % 1 * 60 // 1
    S = angle % 1 * 60 % 1 * 60

    return sgn, D, M, S


def write_solve_src(ivs, RA, DC, errDC, comment, src_file):
    '''Write qso position used for SOLVE.

    Parameters
    ----------
    ivs : array, string
        IVS designation name of radio source
    RA : array, float
        Right ascension, degree
    DC : array, float
        Declination, degree
    errDC : array, float
        Formal uncertainty of DC, micro-arcsec
    comment : string
        comment
    src_file : string
        Output file

    Returns
    ----------
    None
    '''

    _, RAh, RAm, RAs = deg2dms(RA / 15.0)     # arc-sec -> second
    DCsgn, DCd, DCm, DCs = deg2dms(DC)

    sgn = np.where(DCsgn == -1, '-', ' ')

    # Loop for writing data
    lfmt = "    %8s  %02d %02d %11.8f   %s%02d %02d  %10.7f   %7.3f   %s"

    fopt = open(src_file, "w")

    for (souni, RAhi, RAmi, RAsi, gni, DCdi, DCmi, DCsi, errDCi) \
            in zip(ivs, RAh, RAm, RAs, sgn, DCd, DCm, DCs, errDC):

        print(lfmt % (souni, RAhi, RAmi, RAsi,
                      sgni, DCdi, DCmi, DCsi, errDCi, comment),
              file=fopt)

    fopt.close()


def write_nnrs(soun, list_file):
    '''Write the No_Net_Rotation constraint source list.

    Parameters
    ----------
    soun: array, string
        IVS/IERS designation name of radio source
    list_file: Output file

    Returns
    ----------
    None
    '''

    N = 8  # 8 source pre line

    fopt = open(list_file, "w")

    for i, souni in enumerate(soun):
        if not i % N:
            if i:
                print("\\\n     ", end="", file=fopt)
            else:
                print("     ", end="", file=fopt)

        print("%-9s" % souni, end="", file=fopt)

    fopt.close()


def main():
    print("Nothing to do at this moment!")


if __name__ == '__main__':
    main()

# --------------------------------- END -------------------------------
