#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: helmert_trans.py
"""
Created on Sun Dec 17 06:21:29 2017

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np
from numpy import concatenate
from functools import reduce
import sys


__all__ = ["trans_fitting", "helmert_trans"]


# -----------------------------  FUNCTIONS -----------------------------
def rad_to_as(rad):
    """Radian to arc-second.
    """

    return np.rad2deg(rad) * 3600


def rad_to_mas(rad):
    """Radian to milli-arc-second.
    """

    return np.rad2deg(rad) * 3600 * 1000


def vecmod_calc(x):
    """Calculate the module of an vector
    """
    return np.sqrt(np.sum(x ** 2))


def vecerr_calc(par, err):
    """Calculate formal error of vector.
    """
    return np.sqrt(np.dot(par**2, err**2))


def wgtmat_calc(dx_err, dy_err, dz_err,
                xy_cov, xz_cov, yz_cov):
    """Generate the weighted matrix.

    Parameters
    ----------
    dx_err, dy_err, dz_err : array of float
        formal uncertainty of X/Y/Z-component of station positions
        or velocities
    xy_cov/xz_cov/yz_cov : array of float
        covariance between dRA and dDE in uas^2, default is None

    Returns
    ----------
    wgt_mat : matrix
        weighted matrix used in the least squares fitting.
    """

    # Covariance matrix
    err = concatenate((dx_err, dy_err, dz_err), axis=0)
    cov_mat = np.diag(err**2)

    # Take the correlation into consideration.
    num = dx_err.size

    for i in range(num):
        # X-Y
        cov_mat[i, i+num] = xy_cov[i]
        cov_mat[i+num, i] = xy_cov[i]

        # X-Z
        cov_mat[i, i+num*2] = xz_cov[i]
        cov_mat[i+num*2, i] = xz_cov[i]

        # Y-Z
        cov_mat[i+num, i+num*2] = yz_cov[i]
        cov_mat[i+num*2, i+num] = yz_cov[i]

    # print(cov_mat)
    wgt_mat = np.linalg.inv(cov_mat)

    return wgt_mat


def jacmat_calc(x, y, z):
    '''Generate the Jacobian matrix.

    Parameters
    ----------
    x/y/z : X/Y/Z-component of station positions

    Returns
    ----------
    jac_mat/jac_mat_t : matrix
        Jacobian matrix and its transpose matrix
    '''
    # Jacobian matrix.

    # For \Delta_x
    # translation components
    parx_t1 = np.ones_like(x)
    parx_t2 = np.zeros_like(x)
    parx_t3 = np.zeros_like(x)
    # scale factor
    parx_d = x
    # rotation angles
    parx_r1 = np.zeros_like(x)
    parx_r2 = y
    parx_r3 = -z

    # For \Delta_y
    # translation components
    pary_t1 = np.zeros_like(y)
    pary_t2 = np.ones_like(y)
    pary_t3 = np.zeros_like(y)
    # scale factor
    pary_d = y
    # rotation angles
    pary_r1 = -z
    pary_r2 = np.zeros_like(y)
    pary_r3 = x

    # For \Delta_z
    # translation components
    parz_t1 = np.zeros_like(z)
    parz_t2 = np.zeros_like(z)
    parz_t3 = np.ones_like(z)
    # scale factor
    parz_d = z
    # rotation angles
    parz_r1 = y
    parz_r2 = -x
    parz_r3 = np.zeros_like(z)

    # (dx1, dy1, dz1, dx2, dy2, dz3, ...).
    part1 = concatenate((parx_t1, pary_t1, parz_t1), axis=0)
    part2 = concatenate((parx_t2, pary_t2, parz_t2), axis=0)
    part3 = concatenate((parx_t3, pary_t3, parz_t3), axis=0)
    pard = concatenate((parx_d, pary_d, parz_d), axis=0)
    parr1 = concatenate((parx_r1, pary_r1, parz_r1), axis=0)
    parr2 = concatenate((parx_r2, pary_r2, parz_r2), axis=0)
    parr3 = concatenate((parx_r3, pary_r3, parz_r3), axis=0)

    # Jacobian matrix.
    N = part1.size
    jac_mat = concatenate((part1.reshape(N, 1), part2.reshape(N, 1),
                           part3.reshape(N, 1), pard.reshape(N, 1),
                           parr1.reshape(N, 1), parr2.reshape(N, 1),
                           parr3.reshape(N, 1)), axis=1)

    jac_mat_t = np.transpose(jac_mat)

    return jac_mat, jac_mat_t


def trans_fitting(dx, dy, dz, dx_err, dy_err, dz_err,
                  xy_cov, xz_cov, yz_cov, x, y, z):

    # Jacobian matrix and its transpose.
    jac_mat, jac_mat_t = jacmat_calc(x, y, z)

    # Weighted matrix.
    wgt_mat = wgtmat_calc(dx_err, dy_err, dz_err, xy_cov, xz_cov, yz_cov)

    # Calculate matrix A and b of matrix equation:
    # A * w = b.
    A = reduce(np.dot, (jac_mat_t, wgt_mat, jac_mat))

    dpos = np.hstack((dx, dy, dz))
    b = reduce(np.dot, (jac_mat_t, wgt_mat, dpos))

    # Solve the equations.
    # w = (t1, t2, t3, d, r1, r2, r3)
    w = np.linalg.solve(A, b)

    # Covariance.
    cov = np.linalg.inv(A)
    sig = np.sqrt(cov.diagonal())

    # Correlation coefficient.
    corr_mat = np.array([cov[i, j] / sig[i] / sig[j]
                         for j in range(len(w))
                         for i in range(len(w))])
    corr_mat.resize((len(w), len(w)))

    return w, sig, corr_mat


def helmert_trans(dx, dy, dz, dx_err, dy_err, dz_err,
                  xy_cov, xz_cov, yz_cov, x, y, z,
                  data_type="p", fout=sys.stdout):
    """Hermert transformation.
    """

    # Name of estimated parameters.
    # xname = ['T_1', 'T_2', 'T_3', 'd', 'R_1', 'R_2', 'R_3']

    x, sig, corr_mat = trans_fitting(
        dx, dy, dz, dx_err, dy_err, dz_err,
        xy_cov, xz_cov, yz_cov, x, y, z)

    [tx,  ty,  tz,  d, rx,  ry,  rz] = x
    [tx_err,  ty_err,  tz_err,  d_err, rx_err,  ry_err,  rz_err] = sig

    # Calculate the modulus
    t = vecmod_calc(np.array([tx,  ty,  tz]))
    t_err = vecerr_calc(np.array([tx,  ty,  tz]),
                        np.array([tx_err,  ty_err,  tz_err]))

    r = vecmod_calc(np.array([rx,  ry,  rz]))
    r_err = vecerr_calc(np.array([rx,  ry,  rz]),
                        np.array([rx_err,  ry_err,  rz_err]))

    # convert the unit of scale factor into ppb
    d = d * 1e9
    d_err = d_err * 1e9

    # Convert the unit of rotation into arc-second.
    rx, ry, rz = rad_to_mas(rx), rad_to_mas(ry), rad_to_mas(rz)
    r, r_err = rad_to_mas(r), rad_to_mas(r_err)
    rx_err = rad_to_mas(rx_err)
    ry_err = rad_to_mas(ry_err)
    rz_err = rad_to_mas(rz_err)

    # For log file.
    if data_type == "p":
        print("#### Translation component (mm):\n",
              " %+8.3f +/- %8.3f |" * 3 % (tx, tx_err, ty, ty_err, tz, tz_err),
              "=> %8.3f +/- %8.3f" % (t, t_err), file=fout)
        print("#### Scale factor (ppb）:\n",
              " %7.3f +/- %7.3f" % (d, d_err), file=fout)
        print("#### Rotation component（mas）:\n",
              "  %+8.3f +/- %8.3f |" *
              3 % (rx, rx_err, ry, ry_err, rz, rz_err),
              "=> %+8.3f +/- %8.3f" % (r, r_err), file=fout)
        print("##   correlation coefficients are:\n", corr_mat, file=fout)
    else:
        print("#### Translation component (mm/yr):\n",
              " %+8.3f +/- %8.3f |" * 3 % (tx, tx_err, ty, ty_err, tz, tz_err),
              "=> %8.3f +/- %8.3f" % (t, t_err), file=fout)
        print("#### Scale factor (ppb/yr）:\n",
              " %7.3f +/- %7.3f" % (d, d_err), file=fout)
        print("#### Rotation component（mas/yr）:\n",
              "  %+8.3f +/- %8.3f |" *
              3 % (rx, rx_err, ry, ry_err, rz, rz_err),
              "=> %+8.3f +/- %8.3f" % (r, r_err), file=fout)
        print("##   correlation coefficients are:\n", corr_mat, file=fout)


# ------------------------------- MAIN --------------------------------
if __name__ == '__main__':
    pass
# --------------------------------- END --------------------------------
