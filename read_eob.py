#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: read_eob.py
"""
Created on Thu Mar 21 10:16:00 2019

@author: Neo(liuniu@smail.nju.edu.cn)

Retrieve the estimates of X pole coordinate, Y pole coordinate, UT1-TAI
angle, UT1 rate, daily offsets of nutation angles as well as their formal
uncertainties and correlations from .eob file which is generated by
the program getpar.

.eob file contains series of the estimates of X pole coordinate,
Y pole coordinate, UT1-TAI angle, UT1 rate, daily offsets of nutation angles
as well as their formal uncertainties and correlations. Time tag and database
name is attached to each line. .EOB format is an extension of the IERS EOP
format.

   File contains lines of three types:
1) Comment. The first character is #. Header comments contain some information
   about solution.

2) Header. The first two symbols are blank. Header lines contain titles of the
   columns

3) Estimates.


  1    1-1    A1     ---     Usage flag
  2    3-14   F12.6  days    Modified Julian date of the TDT time tag for
                             pole coordinates and UT1
  3   16-25   A10    ---     Database name
  4   27-32   A6     ---     IVS session code (if available)
  5   34-41   F8.6   arcsec  The estimate of X pole coordinate
  6   43-50   F8.6   arcsec  The estimate of Y pole coordinate
  7   52-62   F11.7  sec     The UT1-TAI function
  8   64-71   F8.3   mas     Adjustment of the nutation in longitude angle with
                                        respect to IAU 1980 nutation expansion
  9   73-80   F8.3   mas     Adjustment of the nutation in obliquity angle with
                                        respect to IAU 1980 theory
 10   82-90   F9.6   asc/day The estimate of X pole rate
 11   92-100  F9.6   asc/day The estimate of Y pole rate
 12  102-108  F7.4   ms/day  The estimate of UT1 rate
 13  110-117  F8.6   arcsec  Formal uncertainty of X pole coordinate
 14  119-126  F8.6   arcsec  Formal uncertainty of Y pole coordinate
 15  128-136  F9.7   sec     Formal uncertainty of UT1-UTC function
 16  138-144  F7.3   mas     Formal uncertainty of nutation in longitude angle
 17  146-152  F7.3   mas     Formal uncertainty of nutation in obliquity angle
 18  154-162  F9.6   asc/day Formal uncertainty of X pole rate
 19  164-172  F9.6   asc/day Formal uncertainty of Y pole rate
 20  174-180  F7.4   asc/day Formal uncertainty of UT1 rate
 21  182-187  F6.4   --      Correlation between the estimates of X-pole
                                          positions and Y-pole position
 22  189-194  F6.4   --      Correlation between the estimates of X-pole
                                         positions and UT1-TAI angle
 23  196-201  F6.4   --      Correlation between the estimates of Y-pole
                                         positions and UT1-TAI angle
 24  203-208  F6.4   --      Correlation between the estimates of nutation in
                                         longitude and nutation in obliquity
 25  210-215  F6.4   --      Correlation between the estimates of X-pole
                                          positions and UT1 rate
 26  217-222  F6.4   --      Correlation between the estimates of Y-pole
                                         positions and UT1-TAI date
 27  224-229  F6.4   --      Correlation between the estimates of
                                         UT1-TAI angle UT1 rate
 28  231-235  F5.2   hours   Session duration
 29  237-243  F7.2   psec    Weighted root mean square of postfit residuals
 30  245-250  I6     --      Number of used observations in the session
 31  252-263  F12.6  days    Modified Julian date for nutation at TDT time
                             scale
 32  265-328  A64    --      The network configuration line. Consists of
                             two characters IVS station codes listed
                             in alphabetic order for stations that participated
                             in the experiment and supplied the data that have
                             been used in processing this experiment.

If the specific parameter was not estimated in this experiment, the field
for its value and formal uncertainty is replaced by filler: $$$$$$. The filler
takes entire field.
"""

import numpy as np


# -----------------------------  FUNCTIONS -----------------------------


# --------------------------------- END --------------------------------
