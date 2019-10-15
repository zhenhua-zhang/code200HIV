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
def setup_wd(path="qtlmapping_report"):
    """Setup working directory at given path

    The structure of the output directory
    output_dir
    |-- logs
    |-- reports
    |   |-- GenotypeLevelPlots
    |   |-- QTLInformation
    |   |-- LocusZoomPlots
    |   |-- ManhattanPlots
    |   |-- Preprocessing
    |   |-- PathwayAnalysis
    |-- report.tsv

    TODO:
        1. report.tsv need a head line?
    """
    os.mkdir(path)
    os.mkdir(os.path.join(path, "reports"))
    compulsary_dirs = [
        "GenotypeLevelPlots", "QTLInformation", "LocusZoomPlots",
        "ManhattanPlots", "Preprocessing", "PathwayAnalysis",
    ]

    for _dir in compulsary_dirs:
        os.mkdir(os.path.join(path, _dir))


def not_runnable():
    """A funnction to warn the user that the script is not runnable."""
    log_me.warn("Not a runnable script!!!!", file=sys.stderr)


def not_implemented():
    """A funnction to warn the user that the function isn't implemented and exit."""
    log_me.warn("Not implemented yet!!!!", file=sys.stderr)
    sys.exit(-1)


class UnknownCorrelationMethod(Exception):
    pass

class ArgumentDependencyError(Exception):
    pass


if __name__ == "__main__":
    not_runnable()
