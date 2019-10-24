#!/usr/bin/env python3
# -*- UTF-8 -*-
import os
import sys

import numpy as np
import pandas as pd

from utls import setup_wd
from genotype import Genotype
from phenotype import Phenotype
from covariate import Covariate


def simulater(pntp, cvrt, seed=3541):
    """Make simulated phenotype measurements and covariates(e.g age, gender)"""
    np.random.seed(seed)

    n_pntp = int(np.random.randint(5, 10))
    pntp_names = ["pheno_{}".format(x) for x in range(n_pntp)]

    n_samples = int(np.random.randint(150, 250))
    sample_names = ["sample_{}".format(x) for x in range(n_samples)]

    pntp_msmt = 10 + np.random.standard_normal((n_samples, n_pntp))
    pntp_dtfm = pd.DataFrame(pntp_msmt, index=sample_names, columns=pntp_names)
    pntp_dtfm.to_csv(pntp, sep="\t")
    
    n_cvrt = int(np.random.randint(3, 5))
    cvrt_names = ["cvrt_{}".format(x) for x in range(n_cvrt)]

    covariates_cont = 5 * np.around(np.random.standard_gamma(0.5, (n_cvrt-1, n_samples)), 2)
    covariates_disr = np.random.binomial(1, 0.6, (1, n_samples))
    covariates = np.concatenate((covariates_cont, covariates_disr), axis=0)

    cvrt_dtfm = pd.DataFrame(covariates, index=cvrt_names, columns=sample_names)
    cvrt_dtfm.to_csv(cvrt, sep="\t")


# pntp_file = "phenotypes.tsv"
# cvrt_file = "covariates.tsv"
# simulater(pntp_file, cvrt_file)

# if not os.path.exists("qtlmapping_report"):
#     setup_wd()

# phnt = Phenotype("./", "qtlmapping_report")
# phnt.load_files("./", ptrn="test0.tsv") \
#     .load_files("./", ptrn="test1.tsv") \
#     .load_files("./", ptrn="test2.tsv") \
#     .remove_cols("pheno_0") \
#     .remove_rows(["sample_1", "sample_111"]) \
#     .check_norm() \
#     .transform_by(target=["pheno_1", "pheno_2"]) \
#     .pop_bad_phenos(threshold=0.1) \
#     .calc_corr() \
#     .calc_pca(pad=2) \
#     .mask_outliers(pad=2) \
#     .transpose_matrix() \
#     .into_report(dump_pw_corr=True)

# target = [
#     "rs892665", "rs6111385", "rs6045990", "rs6046657", "rs1836444",
#     "rs6036082", "rs62190471", "rs34383360", "rs1935386", "rs6051504"
# ]
# gntp = Genotype(gntp_path="../test/genotypes", info_path="../test/genotypes")
# gntp.load_gntp(sep="\t") \
#     .load_gntp_info(sep=" ") \
#     .mask_maf(0.05)

dtfm = pd.DataFrame(dict(a=range(10), b=))
