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

    def __init__(self, gntp_path, gntp_ptrn="*dosage.gz", info_path=None,
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
        super().__init__(gntp_path, gntp_ptrn)

        self.info_path = info_path  # The path to the variant infomation
        self.info_ptrn = info_ptrn  # The pattern of the file names of variant information
        self.rare_vars = pd.DataFrame()
        self.gntp_dtfm = pd.DataFrame()
        self.gntp_info_dtfm = pd.DataFrame()

    def load_gntp(self, **kwargs):
        """Load genotypes.

        This is an convenient API of `PreProc.load_files()` to load genotype files

        Args:
            - **kwargs:  Any supported arguments for PreProc.load_files()
        Raises:
        Rturns:
            - self: the instance itself
        """
        super().load_files(dtfm_sig="gntp", **kwargs)
        return self

    def load_gntp_info(self, **kwargs):
        """Load genotype infomations.

        This is an convenient API of `PreProc.load_files()` to load information of genotype files

        Args:
            - **kwargs (): Any supported arguments for PreProc.load_files()
        Raises:
        Rturns:
            - self: the instance itself
        """
        super().load_files(self.info_path, self.info_ptrn, dtfm_sig="info", **kwargs)
        return self

    @staticmethod
    def _encode(record, trmx):
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
        target = self._wrap_into_list(target)

        if dosage_dtfm is None:
            dosage_dtfm = self.get_dtfm(0, 'gntp')

        if info_dtfm is None:
            info_dtfm = self.get_dtfm(1, 'info')

        target_inf = info_dtfm.loc[target, :]
        target_dsg = dosage_dtfm.loc[target, :]
        self.gntp_dtfm = target_dsg \
            .round() \
            .apply(self._encode, axis=1, trmx=target_inf)

        return self

    @staticmethod
    def _check_maf(record: pd.Series, threshold: float =0.05):
        """Check the minor allele frequency
        Currently only called by Genotype.mask_maf().

        Args:
            record (pandas.Series): A container include tricode-genotype to be estimated. Required
            threshold (float): Threshold of MAF. Default: 0.05
        Returns:
            True / False: True if the tested MAF is smaller than threshold; False if not.
        """
        total_allele_n = 2 * len(record)
        alt_allele_n = record.sum()

        if threshold <= float(alt_allele_n) / total_allele_n < 1 - threshold:
            return False

        return True

    def mask_maf(self, threshold: float =0.05):
        """Mask out SNPs with MAF lower than given threshold.

        Args:
            threshold (float): [0.05] Genotypes with MAF < threshold will be removed
        Raises:
        Returns:
            self: the instance itself
        """
        is_rare = self.wk_dtfm \
            .round() \
            .apply(self._check_maf, axis=1, threshold=threshold)

        self.rare_vars = self.wk_dtfm.loc[is_rare, :]
        self.wk_dtfm = self.wk_dtfm.loc[~is_rare, :]

        return self

    def concat_gntp_info(self):
        """Itegrate genotypes and information of genotype.

        """
        # FIXME: Duplicated function with encode

        dosage_dtfm = self.get_dtfm(0, "gntp")
        info_dtfm = self.get_dtfm(1, "info")
        self.gntp_info_dtfm = pd.concat([dosage_dtfm, info_dtfm], axis=1)

        return self


if __name__ == "__main__":
    not_runnable()
