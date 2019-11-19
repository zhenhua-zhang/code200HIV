# -*- UTF-8 -*-

#
## genotype.py
#

"""A module use to read genotype files. However it isn't used too often as they
regularly well formated to MatrixEQTL package.

Notes:
    1. Split dosages into files by chromosome, and name the file following the
      pattern: chr[1-22X]_dosage.gz and chr[1-22X]_variantInfo.gz. By this,
      the tool will infer the chromsome by itself, which make life much easier.
TODOs:
"""

import os
import copy
import pandas as pd

from preproc import PreProc
from utls import log_me, not_runnable

class Genotype(PreProc):
    """A genotype class.

    Attributes:
    Methods:
    InherenteFrom:
    Examples:
    Todos:
        1. Plan to hide load_files() from PreProc class, instead, use load_gntp() and load_gntp_info() as wrappers of load_files() to load genotype dosage and genotype information files.
    Notes:
    """

    def __init__(self, wk_dir, gntp_path, gntp_ptrn="*dosage.gz", info_path=None,
                 info_ptrn="*variantInfo.gz"):
        """Construct an Genotype instance.

        Args:
            - path (None, str): Genotypes path that could be a normal-, zipped-, gzipped-file or directory including phenotypes. Default: None
            - info_path (None, str): The path of information of genotypes, i.e.  genotypes in bases. Default: None
            - ptns (str): Patterns used to find genotype files.
            - info_ptns (str): Patterns used to find information of genotypes.  Default: *variantInfo.gz
        Raises:
        Returns:
        Todos:
            1. A regular expression could be more professional? Default: *dosage.gz
        """
        super(Genotype, self).__init__(gntp_path, gntp_ptrn)

        self.wk_dir = wk_dir
        self.info_path = info_path  # The path to the variant infomation
        self.info_ptrn = info_ptrn  # The pattern of the file names of variant information
        self.rare_vars = pd.DataFrame()
        self.gntp_dtfm = pd.DataFrame()

    def load_gntp(self, **kwargs):
        """Load genotypes.

        This is an convenient API of `PreProc.load_files()` to load genotype files

        Args:
            - **kwargs:  Any supported arguments for PreProc.load_files()
        Raises:
        Rturns:
            - self: the instance itself
        """
        if "dtfm_sig" not in kwargs:
            kwargs["dtfm_sig"] = "gntp"

        self.load_files(**kwargs)
        return self

    def load_gtif(self, **kwargs):
        """Load genotype infomations.

        This is an convenient API of `PreProc.load_files()` to load information of genotype files

        Args:
            - **kwargs (): Any supported arguments for PreProc.load_files()
        Raises:
        Rturns:
            - self: the instance itself
        """
        if "dtfm_sig" not in kwargs:
            kwargs["dtfm_sig"] = "gtif"

        self.load_files(self.info_path, self.info_ptrn, **kwargs)
        return self

    @staticmethod
    def encoder(record, trmx):
        """A static method to translate dosage into genotypes.

        Args:
            record (pandas.Series): Each row or column in given pandas.DataFrame.
            trmx (pandas.DataFrame): Translation matrix, index / column names should be consistent with record.
        Currently, only called by Genotype.encode()
        """
        rec_id = record.name
        rec_dct = trmx.loc[rec_id, :].to_dict()
        rec_trmx = {
            0.0: rec_dct["EffectAllele"] * 2,
            1.0: rec_dct["EffectAllele"] + rec_dct["AlternativeAllele"],
            2.0: rec_dct["AlternativeAllele"] * 2
        }

        return record.replace(rec_trmx)

    def encode(self, target, dosage_dtfm=None, info_dtfm=None):
        """Encode dosage into genotypes.

        Args:
            - target (str, list, tuple): SNP will be encoded. Required
            - dosage_dtfm (pandas.DataFrame, None): Dataframe including dosage. Default: None
                If it's None dataframe with signature 'gntp' will be used.
            - info_dtfm (pandas.DataFrame, None): Dataframe including infomation of dosage. Default: None
                If it's None dataframe with signature 'info' will be used, and the info_dtfm should be correspoinding to dosage_dtfm.
        Raises:
        Returns:
            self: the instance itself.
        Todos:
        """

        log_me.info("Encode dosage into genotypes")
        target = self._wrap_into_list(target)

        if dosage_dtfm is None:
            dosage_dtfm = self.get_dtfm('gntp')

        if info_dtfm is None:
            info_dtfm = self.get_dtfm('gtif')

        target_inf = info_dtfm.loc[target, :]
        target_dsg = dosage_dtfm.loc[target, :]
        self.gntp_dtfm = target_dsg \
            .round() \
            .apply(self.encoder, axis=1, trmx=target_inf)

        return self

    def mask_maf(self, threshold: float = 0.05):
        """Mask out SNPs with MAF lower than given threshold.

        Args:
            threshold (float): [0.05] Genotypes with MAF < threshold will be removed
        Raises:
        Returns:
            self: the instance itself
        """
        log_me.info("Mask low MAF SNPs, the threshold is {}".format(threshold))

        _, n_samples = self.wk_dtfm.shape
        self.wk_dtfm["AAFreq"] = self.wk_dtfm.round().sum(axis=1) / (2.0 * n_samples)

        is_normal = self.wk_dtfm.loc[:, "AAFreq"].between(threshold, 1 - threshold)

        self.dtfm_bin["rare_vars"] = copy.deepcopy(self.wk_dtfm.loc[~is_normal, :])
        self.wk_dtfm = self.wk_dtfm.loc[is_normal, :]

        del self.wk_dtfm["AAFreq"]

        return self

    def dump_gntp(self, opt_dir=None, opt_fnm="genotypes.tsv"):
        """Dump processed genotypes into disk."""
        if opt_dir is None:
            opt_dir = os.path.join(self.wk_dir, "Preprocessing")

        log_me.info("Dump processed genotypes to {} in {}".format(opt_fnm, opt_dir))
        self.dump_wkdtfm(opt_dir, opt_fnm)

        return self

    def dump_rare_vars(self, opt_dir=None, opt_fnm="genotypes_rareVarsMAFle0.05.tsv"):
        """Dump processed genotypes into disk."""
        if opt_dir is None:
            opt_dir = os.path.join(self.wk_dir, "Preprocessing")

        log_me.info("Dump rare variants to {} in {}".format(opt_fnm, opt_dir))
        self.set_wkdf("rare_vars")
        self.dump_wkdtfm(opt_dir, opt_fnm)

        return self


if __name__ == "__main__":
    not_runnable()
