#!/usr/bin/env python3
#-*- coding: utf-8 -*-
# 这个脚本还没写完

import sys
import csv
import logging
import argparse

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
fmt = logging.Formatter(
    "|{levelname: ^10}| {asctime} | {funcName: <12}| {message}", style="{",
    datefmt="%Y-%m-%d %H:%M:%S"
)
ch.setFormatter(fmt)
logger.addHandler(ch)


def fetch_arg(lg=logger):
    """Parse arguments"""

    parser = argparse.ArgumentParser()
    group = parser.add_argument_group("Inputs")
    group.add_argument(
        "-i", "--input-file", dest="iptf", required=True, type=str,
        help="Input csv file including qurey genes."
    )
    group.add_argument(
        "-g", "--gene-databse", dest="gndb", required=True, type=str,
        help="Databse including gene informations. GFF/GTF/BED"
    )
    group.add_argument(
        "-s", "--snp-database", dest="snpdb", required=True, type=str,
        help="Database include SNPs."
    )

    group = parser.add_argument_group("Parameters")
    group.add_argument(
        "--qry-cols", dest="qrt_cols", default="gene_id", type=str,
        help=""
    )
    group.add_argument(
        "--snp-cols", dest="snp_cols", default=["snp_id", "chr", "pos"],
        type=str, nargs=3, help=""
    )
    group.add_argument(
        "--gene-cols", dest="gene_cols", 
        default=["gene_id", "start", "stop", "strand"], type=str, nargs=4,
        help="Input csv file including qurey genes."
    )
    group.add_argument(
        "-w", "--window-size", dest="win_size", default=5e5, type=int,
        help="Window size. Default: 5E5"
    )
    group.add_argument(
        "--upstream_size", dest="up_win_size", default=-1, type=int,
        help="Upstream window size, confilct with -w/--window-size. Default: -1"
    )
    group.add_argument(
        "--dwstream_size", dest="dw_win_size", default=-1, type=int,
        help="Downstream window size, confilct with -w/--window-size. Default: -1"
    )

    group = parser.add_argument_group("Outputs")
    group.add_argument(
        "-o", "--output-file", dest="optf", default="output.csv", type=str,
        help="Output file. Default: output.csv"
    )

    group = parser.add_argument_group("Misc")
    group.add_argument(
        "-v", "--verbose", dest="verbose", action="count", default=0,
        help="Running verbosely (accelerative)."
    )
    group.add_argument(
        "--debug", dest="debug", action="store_true",
        help="Entering DEBUG mode, verbose will be overwrited. Default: FALSE"
    )

    return parser.parse_args()


def parse_qry(csvfl, dlmt=",", qtch="\"", lg=logger):
    """Parse query"""

    lg.info("Parsing query file: {}".format(csvfl))

    with open(csvfl, newline='') as csvhdl:
        csvit = csv.reader(csvhdl, delimiter=dlmt, quotechar=qtch)

    return ""


def parse_snpdb(snpdb, zipped=False, lg=logger):
    """Parse SNP database"""
    lg.info("Parsing SNP database: {}".format(snpdb))

    return ""


def parse_gndb(gndb, zipped=False, lg=logger):
    """Parse gene database"""
    lg.info("Parsing gene database: {}".format(gndb))
    return ""


def parse_rfr(gndb, snpdb, lg=logger):
    """Parse reference databases"""
    lg.info("Parsing reference start")
    parse_gndb(gndb)
    parse_snpdb(snpdb)

    lg.info("Parsing reference DONE!")
    return ""


def fetch_snp(qry, rfr, optf, lg=logger):
    lg.info("Fetching your SNPs ...")
    return ""


def main(lg=logger):
    """Main function"""
    lg.info("STARTS")
    arguments = fetch_arg()

    debug = arguments.debug

    if(debug):
        verbose = 5
    else:
        verbose = arguments.verbose

    logger.setLevel(50 - verbose*10)

    iptf = arguments.iptf
    gndb = arguments.gndb
    snpdb = arguments.snpdb

    win_size = arguments.win_size
    up_win_size = arguments.up_win_size
    dw_win_size = arguments.dw_win_size

    optf = arguments.optf

    if(up_win_size == -1 or dw_win_size == -1):
        win = (int(win_size), int(win_size))
    else:
        win = (int(up_win_size), int(dw_win_size))

    lg.info("##### Command line arguments #####")
    lg.info("# Input file   : {}".format(iptf))
    lg.info("# SNP database : {}".format(snpdb))
    lg.info("# Gene database: {}".format(gndb))
    lg.info("# Window size  : {}".format(win))
    lg.info("# Output file  : {}".format(optf))
    lg.info("# Verbose      : {}".format(verbose))
    lg.info("# Debug        : {}".format(debug))
    lg.info("##################################")

    qry = parse_qry(iptf)
    rfr = parse_rfr(gndb, snpdb)
    fetch_snp(qry, rfr, optf)

    lg.info("DONE!")

if __name__ == "__main__":
    main()
