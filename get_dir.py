#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File Name: get_dir.py
"""
Created on Thu 15 Aug 2019 10:05:04 CST

@author: Neo (niu.liu@foxmail.com)
"""

import numpy as np
import os
from sys import platform as _platform

__all__ = ["get_home_dir", "get_data_dir"]

# -------------------- FUNCTIONS --------------------


def get_home_dir():
    """Get the path of home directory.

    Return
    ------
    homedir: string
        home directory
    """

    # Check the type of OS and get the home diretory path
    if _platform == "linux" or _platform == "linux2":
        # linux
        homedir = os.getenv("HOME")
    elif _platform == "darwin":
        # MAC OS X
        homedir = os.getenv("HOME")
    elif _platform == "win32":
        # Windows
        print("Not implemented yet")
        exit()
    elif _platform == "win64":
        # Windows 64-bit
        print("Not implemented yet")
        exit()
    else:
        print("Weird! What kind of OS do you use?")
        exit()

    return homedir


def get_data_dir():
    """Get the diretory of source name data

    Returns
    -------
    datadir: string
        diretory to put source name data

    """

    homedir = get_home_dir()

    datadir = "{}/.ipython/my_progs/vlbi/aux_files".format(homedir)

    return datadir
# ----------------------- END -----------------------
