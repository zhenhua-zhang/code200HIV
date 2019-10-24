# -*- UTF-8 -*-
"""Utilities for package `qtlmapping`.
TODO: 
    1. move all error class into a errors.py
"""

import os
import sys
import logging

#
## Loger
#
log_me = logging.getLogger()
log_me.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
fmtstr = "|{levelname: ^10}| {asctime} | {funcName: <20}| {message}"
fmt = logging.Formatter(fmtstr, style="{", datefmt="%Y-%m-%d %H:%M:%S")
ch.setFormatter(fmt)
log_me.addHandler(ch)


#
## A function to help setup working directories
#
def setup_wd(path, override=False):
    """Setup working directory at given path

    The structure of the output directory
    PATH
    |-- logs (Not decided yet)
    |-- PhenotypeLevelPlots/
    |-- QTLInformation/
    |-- LocusZoomPlots/
    |-- ManhattanPlots/
    |-- Preprocessing/
    |-- PathwayAnalysis/
    |-- TempDir/ (Will be removed after calling `clean_up()` function)
    |-- report.tsv (Not decided yet)

    Args:
        path (str): The path to create the working directory. Default: ./qtlmapping_report
    Raises:
    Returns:
    Todos:
    Notes:
    """
    if os.path.exists(path):
        if override:
            log_me.warn(
                "{} exists, but will use it directly.".format(path) + 
                " Start checking mandatory subdirectories"
            )
        else:
            log_me.error("{} exists, exiting ...".format(path))
            sys.exit()
    else:
        os.mkdir(path)

    compulsary_dirs = [
        "PhenotypeLevelPlots", "QTLInformation", "LocusZoomPlots",
        "ManhattanPlots", "Preprocessing", "PathwayAnalysis",
    ]

    for _dir in compulsary_dirs:
        sub_dir = os.path.join(path, _dir) 
        if override:
            if os.path.exists(sub_dir):
                log_me.info("{} exists ...".format(sub_dir))
            else:
                log_me.warn("{} doesn't exist, creating it...".format(sub_dir))
                os.mkdir(sub_dir)
        else:
            os.mkdir(sub_dir)


def clean_up():
    """A function to remove temporary files in `YOUR_WK_DIR/TempDir`.
    """
    log_me.info("Start cleaning up your sh*t ...")
    not_implemented()


def not_runnable():
    """A funnction to warn the user that the script is not runnable."""
    log_me.warn("Not a runnable script!!!!")


def not_implemented():
    """A funnction to warn the user that the function isn't implemented and exit."""
    log_me.warn("Not implemented yet!!!!")
    sys.exit(-1)


class UnknownCorrelationMethod(Exception):
    pass

class ArgumentDependencyError(Exception):
    pass


if __name__ == "__main__":
    not_runnable()
