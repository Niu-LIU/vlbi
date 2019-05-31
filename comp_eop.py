#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: eob_comp.py
import time
"""
Created on Wed Feb 27 14:59:41 2019

@author: Neo(liuniu@smail.nju.edu.cn)

Compare the EOP series and calculate the offsets (or differences).

"""

from astropy.table import Table, join, Column
import astropy.units as u
from astropy.units import cds
import numpy as np
import sys
import os
import time


__all__ = {"calc_eop_offset", "save_eop_offset_txt"}


# -----------------------------  FUNCTIONS -----------------------------
def root_sum_square(x, y):
    """Calculate the root-sum-square."""

    return np.sqrt(x**2 + y**2)


def calc_eop_offset(t_eop1,  t_eop2):
    """Calculate the EOP difference series between two solutions.

    Parameters
    ----------
    t_eop1, t_eop2 : astropy.table object
        EOP series from two solutions

    Return
    ------
    t_eop_offset : astropy.table object
        EOP difference series
    """

    # Copy the original tables and keep only the EOP information
    t_eop3 = Table(t_eop1)
    t_eop3.keep_columns(["epoch_pmr", "db_name",
                         "xp", "yp", "ut1_tai", "dX", "dY",
                         "xp_err", "yp_err", "dut1_err", "dX_err", "dY_err",
                         "epoch_nut"])

    t_eop4 = Table(t_eop2)
    t_eop4.keep_columns(["db_name",
                         "xp", "yp", "ut1_tai", "dX", "dY",
                         "xp_err", "yp_err", "dut1_err", "dX_err", "dY_err"])

    # Cross-match between two tables
    t_eop_com = join(t_eop3, t_eop4, keys="db_name")

    print("There are %d and %d points in series 1 and series 2, respectively,"
          "between which %d are common."
          % (len(t_eop1), len(t_eop2), len(t_eop_com)))

    # Calculate the offset and the uncertainties
    dxp = t_eop_com["xp_2"] - t_eop_com["xp_1"]
    dyp = t_eop_com["yp_2"] - t_eop_com["yp_1"]
    dut = t_eop_com["ut1_tai_2"] - t_eop_com["ut1_tai_1"]
    ddX = t_eop_com["dX_2"] - t_eop_com["dX_1"]
    ddY = t_eop_com["dY_2"] - t_eop_com["dY_1"]

    dxp_err = root_sum_square(t_eop_com["xp_err_1"], t_eop_com["xp_err_2"])
    dyp_err = root_sum_square(t_eop_com["yp_err_1"], t_eop_com["yp_err_2"])
    dut_err = root_sum_square(t_eop_com["dut1_err_1"], t_eop_com["dut1_err_2"])
    ddX_err = root_sum_square(t_eop_com["dX_err_1"], t_eop_com["dX_err_2"])
    ddY_err = root_sum_square(t_eop_com["dY_err_1"], t_eop_com["dY_err_2"])

    # Convert the unit
    # Time tag
    from astropy.time import Time
    t_pmr_mjd = Time(t_eop_com["epoch_pmr"], format="mjd")
    t_pmr = Column(t_pmr_mjd.jyear, unit=u.year)

    t_nut_mjd = Time(t_eop_com["epoch_nut"], format="mjd")
    t_nut = Column(t_nut_mjd.jyear, unit=u.year)

    # Polar motion
    dxp.convert_unit_to(u.uas)
    dxp_err.convert_unit_to(u.uas)
    dyp.convert_unit_to(u.uas)
    dyp_err.convert_unit_to(u.uas)

    # UT1-UTC (s -> us)
    dut.convert_unit_to(u.second / 1e6)
    dut_err.convert_unit_to(u.second / 1e6)

    # Nutation offset
    ddX.convert_unit_to(u.uas)
    ddX_err.convert_unit_to(u.uas)
    ddY.convert_unit_to(u.uas)
    ddY_err.convert_unit_to(u.uas)

    # Add these columns to the combined table.
    t_eop_offset = Table([t_eop_com["db_name"], t_pmr, t_nut,
                          dxp, dyp, dut, ddX, ddY,
                          dxp_err, dyp_err, dut_err, ddX_err, ddY_err],
                         names=["db_name", "epoch_pmr", "epoch_nut",
                                "dxp", "dyp", "dut", "ddX", "ddY",
                                "dxp_err", "dyp_err", "dut1_err",
                                "ddX_err", "ddY_err"])

    return t_eop_offset


def save_eop_offset_txt(t_eop_offset, offset_file):
    """Save the eop offset series into a text file.

    Parameters
    ----------
    t_eop_offset : Table object
        eop offset series
    offset_file : string
        name of file to store the data
    """

    # Header
    t_eop_offset.meta["comments"] = [
        " EOP Offset series",
        " Columns  Units   Meaning",
        "    1     day     Time Tag for polar motion and UT1 (MJD)",
        "    2     day     Time Tag for Nutation (MJD)",
        "    3     uas     offset of X pole coordinate",
        "    4     uas     offset of Y pole coordinate",
        "    5     us      offset of UT1",
        "    6     uas     offset of dX of Nutation offsets",
        "    7     uas     offset of dY of Nutation offsets",
        "    8     uas     formal uncertainty for offset of X pole coordinate",
        "    9     uas     formal uncertainty for offset of Y pole coordinate",
        "   10     us      formal uncertainty for offset of UT1",
        "   11     uas     formal uncertainty for offset of Nutation dX",
        "   12     uas     formal uncertainty for offset of Nutation dY",
        " Created date: %s." % time.strftime("%d/%m/%Y", time.localtime())]

    t_eop_offset.write(offset_file, format="ascii.fixed_width_no_header",
                       exclude_names=["db_name"],
                       formats={"epoch_pmr": "%14.6f", "epoch_nut": "%14.6f",
                                "dxp": "%+8.1f", "dyp": "%+8.1f",
                                "dut": "%+8.1f",
                                "ddX": "%+8.1f", "ddY": "%+8.1f",
                                "dxp_err": "%+8.1f", "dyp_err": "%+8.1f",
                                "dut1_err": "%+8.1f",
                                "ddX_err": "%+8.1f", "ddY_err": "%+8.1f", },
                       delimiter="", overwrite=True)


def read_eop_offset(eop_oft_file):
    """Read EOP offset data.

    Parameters
    ----------
    eop_oft_file : string
        EOP offset file

    Returns
    ----------
    t_eop_oft : astropy.table object
    """

    if not os.path.isfile(eop_oft_file):
        print("Couldn't find the file", eop_oft_file)
        sys.exit()

    t_eop_oft = Table.read(eop_oft_file, format="ascii",
                           names=["epoch_pmr", "epoch_nut",
                                  "dxp", "dyp", "dut", "ddX", "ddY",
                                  "dxp_err", "dyp_err", "dut_err",
                                  "ddX_err", "ddY_err"])

    # Add unit information
    t_eop_oft["epoch_pmr"].unit = cds.MJD
    t_eop_oft["epoch_nut"].unit = cds.MJD

    t_eop_oft["dxp"].unit = u.uas
    t_eop_oft["dyp"].unit = u.uas
    t_eop_oft["dxp_err"].unit = u.uas
    t_eop_oft["dyp_err"].unit = u.uas

    t_eop_oft["dut"].unit = u.second / 1e6
    t_eop_oft["dut_err"].unit = u.second / 1e6

    t_eop_oft["ddX"].unit = u.uas
    t_eop_oft["ddY"].unit = u.uas
    t_eop_oft["ddX_err"].unit = u.uas
    t_eop_oft["ddY_err"].unit = u.uas

    return t_eop_oft


if __name__ == '__main__':
    print('Nothing to do!')
    pass
# --------------------------------- END --------------------------------
