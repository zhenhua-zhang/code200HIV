#!/usr/bin/env python3
# -*- utf-8 -*-

"""A pipeline to do QTL mapping.

Notes:
    1. Currently the main mapping algorithm is backended by R package MatrixEQTL.
ToDos:
    - [ ] 1. Replace hard coded variable by config files.
    - [ ] 2. Update docstrings.
    - [ ] 3. A class to replace MatrixEQTL wrapper.
    - [ ] 4. Remove samples by intersecting genotypes, phenotypes and covariates files.
    - [ ] 5. More fancy analysis. Locus zoom, pathway enrichment.
    - [ ] 6. More efficient implementation, e.g parallel processing, SQL.
"""

#
## mapper.py
#

# Standard library
import math
import subprocess
from argparse import ArgumentParser

# Third party

# Customer modules
from covariate import Covariate
from phenotype import Phenotype
from genotype import Genotype
from report import Reporter
from utls import *

# black list
# "X1009", "X1012", "X1053", "X1068", "X1101", "X1109", "X1110", "X1125",
# "X1126", "X1128", "X1129", "X1142", "X1150", "X1165", "X1193", "X1214"
blist = [
    "X1009", "X1021", "X1053", "X1068", "X1109", "X1110", "X1125", "X1126",
    "X1128", "X1129", "X1142", "X1150", "X1161", "X1165", "X1193"
]

#
## Parsing commandline options
#
def get_arg():
    """Get CLI arguments for current script.

    Args:
    Raise:
    Returns:
        parser (argparse.ArigumentParser):
    """
    parser = ArgumentParser()
    parser.add_argument("--config", dest="config", type=str, help="Configuration file")
    sub_parser = parser.add_subparsers(dest="subcmd", help="Preprocessing")

    map_parser = sub_parser.add_parser("map", help="Mapping QTLs")
    _group = map_parser.add_argument_group("Input")
    _group.add_argument(
        "-w", "--wk-dir", dest="wk_dir", type=str, required=True, help="Working dir. Required"
    )

    _group.add_argument(
        "-g", "--genotype", dest="gntp_inpt", type=str, required=True, help="Genotype. Required"
    )

    _group.add_argument(
        "-p", "--phenotype", dest="pntp_inpt", type=str, required=True, help="Phenotype. Required"
    )

    _group.add_argument(
        "-c", "--covariate", dest="cvrt_inpt", type=str, required=True, help="Covariates. Required"
    )

    return parser

#
## Main entry of the QLT-mapping pipeline
#
def main():
    """The main entry of `qtlmapping` pipeline.
    Todos:
        1. Using configs or commandline arguments to replace hard coding parameters
        2. An argument for each method to decide whether a function will be executed or not.
    """
    parser = get_arg()
    args = parser.parse_args()

    wk_dir = args.wk_dir
    pntp_inpt = args.pntp_inpt
    cvrt_inpt = args.cvrt_inpt
    gntp_inpt = args.gntp_inpt

    # Setup working directories
    # setup_wd(wk_dir, True)

    # Preprocessing phenotypes
    run_pntp = False
    if run_pntp:
        pntp = Phenotype(wk_dir, pntp_inpt)
        pntp.load_files(dtfm_sig="pntp") \
            .remove_rows(blist) \
            .transform_by() \
            .check_norm("log10") \
            .pop_bad_phenos() \
            .calc_pca() \
            .calc_corr() \
            .mask_outliers(pad=4) \
            .check_norm("mask_outliers") \
            .transpose_matrix() \
            .dump_pntp() \
            .into_report(dump_pw_corr=True)

    # Preprocessing covariates
    run_cvrt = False
    if run_cvrt:
        cvrt = Covariate(wk_dir, cvrt_inpt)
        cvrt.load_files(sep=",", dtfm_sig="cvrt") \
            .keep_cols(["age", "gender"]) \
            .remove_rows(blist) \
            .transpose_matrix() \
            .dump_cvrt()

    # Preprocessing genotypes
    run_gntp = False
    if run_gntp:
        gntp = Genotype(
            wk_dir=wk_dir,
            gntp_path=gntp_inpt, gntp_ptrn="chr*dosage.gz",
            info_path=gntp_inpt, info_ptrn="chr*variantInfo.gz"
        )
        gntp.load_gntp() \
            .load_gtif(sep=" ") \
            .remove_cols(blist) \
            .mask_maf() \
            .dump_gntp()

    # QTL mapping wrapper of MatrixEQTL
    pntp_proc = os.path.join(wk_dir, "Preprocessing", "phenotypes.tsv")
    cvrt_proc = os.path.join(wk_dir, "Preprocessing", "covariates.tsv")
    gntp_proc = os.path.join(wk_dir, "Preprocessing", "genotypes.tsv")
    qtls_otpt = os.path.join(wk_dir, "QTLInformation", "qtls.tsv")

    run_map = False
    if run_map:
        subprocess.run(
            [
                "module load R/3.5.1-foss-2015b-bare &&" +
                "Rscript ./mewrapper.R " +
                "-p {} -c {} -g {} -o {}".format(pntp_proc, cvrt_proc, gntp_proc, qtls_otpt)
            ], shell=True, check=True
        )

    # Generate mapping report
    rpot = Reporter(
        wk_dir=wk_dir, qtls=qtls_otpt,
        gntp=gntp_inpt, gntp_ptrn="chr*dosage.gz",
        gtif=gntp_inpt, gtif_ptrn="chr*variantInfo.gz",
        pntp=pntp_inpt
    )
    # rpot.plot_qq().plot_manhattan()
    rpot.plot_pntp_per_gntp(5e-7)


if __name__ == "__main__":
    log_me.info("Job start ...")
    main()
    clean_up()
    log_me.info("Job done!")
