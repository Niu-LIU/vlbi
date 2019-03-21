#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: pos_err.py
"""
Created on Fri Sep 21 15:36:35 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""

from numpy import sqrt, sin, cos


__all__ = ["pos_err_calc", "overall_error", "error_ellipse", "error_ellipse2"]


# -----------------------------  FUNCTIONS -----------------------------
def pos_err_calc(ra_err, dec_err, ra_dec_corr):
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


def overall_error(ra_err, dec_err, ra_dec_corr):
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


def error_ellipse(ra_err, dec_err, ra_dec_corr):
    """Calculate the ovrall formal uncertainty.

    ovrall formal uncertainty = sqrt(RA_err^2+Dec_err^2+C*RA_err*Dec_err)

    Parameters
    ----------
    ra_err/dec_err : formal uncertainty of RA/Dec, usually in micro-as
    ra_dec_corr : correlation coeffient between RA and Dec, unitless.

    Returns
    ----------
    M : semi-major axis of the error ellipse
    m : semi-minor axis of the error ellipse
    pa : the position angle of the major axis of the error ellipse
    """

    cov = ra_dec_corr * ra_err * dec_err
    cov_mat = np.array([[ra_err**2, cov], [cov, dec_err**2]])

    from numpy.linalg import eig

    eig_val, eig_vec = eig(cov_mat)

    M2 = np.max(eig_val)
    m2 = np.min(eig_val)

    M, m = np.sqrt(M2), np.sqrt(m2)

    index = np.where(eig_val == M2)[0][0]

    vec_M = eig_vec[:, index]

    nor_vec = np.array([0, 1])
    # nor_vec = np.array([1, 0])
    pa0 = np.rad2deg(np.arccos(np.dot(vec_M, nor_vec)))

    if vec_M[0] > 0:
        pa = 180 - pa0
    else:
        pa = pa0

    return M, m, pa


def error_ellipse2(ra_err, dec_err, ra_dec_corr):
    """Calculate the ovrall formal uncertainty.

    ovrall formal uncertainty = sqrt(RA_err^2+Dec_err^2+C*RA_err*Dec_err)

    Parameters
    ----------
    ra_err/dec_err : formal uncertainty of RA/Dec, usually in micro-as
    ra_dec_corr : correlation coeffient between RA and Dec, unitless.

    Returns
    ----------
    eema : semi-major axis of the error ellipse
    eena : semi-minor axis of the error ellipse
    pa : the position angle of the major axis of the error ellipse
    """

    a = ra_err
    b = dec_err
    r = ra_dec_corr
    rab = r * a * b

    # This angle is counter-clockwise angle measured from the X-axis
    theta = 0.5 * np.arctan2(2 * rab, (a**2-b**2))

    sx = a**2 * cos(theta)**2 + b**2 * sin(theta)**2 + \
        2 * rab * sin(theta) * cos(theta)
    sy = a**2 * sin(theta)**2 + b**2 * cos(theta)**2 - \
        2 * rab * sin(theta) * cos(theta)
    eema = sqrt(np.maximum(sx, sy))
    eena = sqrt(np.minimum(sx, sy))

    # For the position angle which is supposed be eastward from the north
    # where north is the Y-axis, i.e., clockwise angle measured from Y-axis.
    # (The RA direction is the east while the declination is the north.)
    if theta < 0:
        theta0 = 360 + theta
    if a >= b:
        pa = (360 - (theta0 - 90)) % 360
    else:
        pa = 360 - theta0

    return eema3, eena3, pa


def main():
    """Main function.
    """

    print("This code is just a module of functions.")


if __name__ == "__main__":
    main()
# --------------------------------- END --------------------------------
