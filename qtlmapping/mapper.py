#!/usr/bin/env python3
# -*- utf-8 -*-

#
## mapper.py
#

# Standard library
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
    sub_parser = parser.add_subparsers(help="Preprocessing")
    map_parser = sub_parser.add_parser("map")
    group = map_parser.add_argument_group("Input")
    group.add_argument("-w", "--wk-dir", dest="wk_dir", type=str, required=True, help="Working dir. Required")
    group.add_argument("-g", "--genotype", dest="gntp", type=str, required=True, help="Genotype. Required")
    group.add_argument("-p", "--phenotype", dest="pntp", type=str, required=True, help="Phenotype. Required")
    group.add_argument("-c", "--covariate", dest="cvrt", type=str, required=True, help="Covariates. Required")

    return parser


#
## Preprocessing the phenotypes
#
def proc_phenotypes(wk_dir=None, pntps=None, rblist=None):
    """Preporcessing the phenotypes.

    Args:
        pntps (str): None. File including measurements of phenotypes
        wk_dir (str): None. Path of working directory
    Raises:
    Returns:
        phtp (Phenotype):
    """
    pntp = Phenotype(pntps, wk_dir)
    pntp.load_files() \
        .remove_rows(rblist) \
        .transform_by() \
        .check_norm() \
        .pop_bad_phenos() \
        .calc_pca() \
        .calc_corr() \
        .mask_outliers() \
        .transpose_matrix() \
        .into_report(dump_pw_corr=True)

    return pntp


#
## Preprocessing the genotypes
#
def proc_genotypes(wk_dir=None, gntps=None, gntp_ptrn=None, infos=None,
                   info_ptrn=None, cblist=None):
    """Preprocessing genotypes.
    """
    gntp = Genotype(gntps, gntp_ptrn=gntp_ptrn, info_path=infos, info_ptrn=info_ptrn)
    gntp_proc_dir = os.path.join(wk_dir, "Preprocessing")
    gntp.load_gntp() \
        .load_gntp_info(sep=" ") \
        .concat_gntp_info() \
        .remove_cols(cblist)
#         .mask_maf() \
#         .dump_wkdtfm(gntp_proc_dir, flnm="genotypes.tsv.gz")
    
    return gntp


#
## Preprocessing the covariates
#
def proc_covariates(wk_dir=None, cvrts=None, rblist=None):
    """Preprocessing covariates.
    """
    cvrt = Covariate(cvrts)
    cvrt_proc_dir = os.path.join(wk_dir, "Preprocessing")
    cvrt.load_files(sep=",") \
        .keep_cols(["age", "gender"]) \
        .remove_rows(rblist) \
        .transpose_matrix()
        # .dump_wkdtfm(cvrt_proc_dir, flnm="covariates.tsv")
    
    return cvrt


#
## Subprocesses to run the R package MatrixEQTL
#
def xtnl_mewrapper_r(pntp_path, cvrt_path, gntp_path, opt, mewp='./mewrapper.R'):
    """Run the MatrixEQTL R package using a subprocess"""
    subprocess.run(
        [
            "module load R/3.5.1-foss-2015b-bare"
        ], shell=True
    )

    subprocess.run(
        [
            "$(which Rscript) {} -p {} -c {} -g {} -o {}".format(
                mewp, pntp_path, cvrt_path, gntp_path, opt
            )
        ], shell=True
    )


#
## Dump the QTL mapping report
#
def dump_report(wk_dir, gntp, pntp, annos=None):
    """Generator QTL mapping report
    """
    rpot = Reporter(wk_dir)
    rpot.load_files() \
        .plot_qq() \
        .plot_manhattan(gntp) \
        .plot_pntp_per_gntp(gntp, pntp, gwth=5e-6)

    return rpot

#
## Main entry of the QLT-mapping pipeline
#
def main():
    parser = get_arg()
    args = parser.parse_args()

    wk_dir = args.wk_dir
    pntp = args.pntp
    cvrt = args.cvrt
    gntp = args.gntp

    # setup_wd(wk_dir, override=True)
    # cvrt_obj = proc_covariates(wk_dir=wk_dir, cvrts=cvrt, rblist=blist)
    pntp_obj = proc_phenotypes(wk_dir=wk_dir, pntps=pntp, rblist=blist)

    gntp_obj = proc_genotypes(
        wk_dir=wk_dir, gntps=gntp, gntp_ptrn="chr22_*dosage.gz", infos=gntp,
        info_ptrn="chr22.variantInfo.gz", cblist=blist
    )


    pntp_path = os.path.join(wk_dir, "Preprocessing", "phenotypes.tsv")
    cvrt_path = os.path.join(wk_dir, "Preprocessing", "covariates.tsv")
    gntp_path = os.path.join(wk_dir, "Preprocessing", "genotypes.tsv.gz")
    qtl_opt = os.path.join(wk_dir, "QTLInformation", "qtls.tsv")

    # xtnl_mewrapper_r(pntp_path, cvrt_path, gntp_path, qtl_opt)

    rpot_obj = dump_report(wk_dir, gntp_obj, pntp_obj)

if __name__ == "__main__":
    log_me.info("Job start ...")
    main()
    clean_up()
    log_me.info("Job done!")
