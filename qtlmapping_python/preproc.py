# -*- UTF-8 -*-
"""A module to preprocessing QTL mapping data"""

#
## preproc.py
#

import os
import glob
from os.path import join as pjoin
from os.path import isdir as pisdir
from os.path import isfile as pisfile

import pandas as pd

# Customer scripts in current package
from utls import log_me

class PreProc:
    """A class to preprocess raw data.

    This class is inheriten by phenotype.Phenotype, genotype.Genotype,
    dovariate.Covariate in qtlmapping package.

    Attributes:
    Methods:
    Todos:
        1. Complete the documents for each methods and the class
        2. A method to summarize the working dataframe
        3. A container to record operations self.operations
    """
    def __init__(self, inpt_path: str = "./", inpt_ptrn: str = "*.tsv"):
        """A class to preprocess data.

        Args:
            inpt_path (str, optional): The path of input file of dir. Default: ./
            inpt_ptrn (str, optional): Pattern to search input file. Default: *.tsv
        """
        self.inpt_path = inpt_path
        self.inpt_ptrn = inpt_ptrn
        self.dtfm_bin = {}         # A bin of loaded files
        self.dtfm_act = {}         # Tracks whether the dataframe is activated
        self.wk_dtfm = pd.DataFrame()
        self.rmd_cols = pd.DataFrame()
        self.rmd_rows = pd.DataFrame()
        self.mltp_index = False
        self.trps_n = 0
        self.pw_corr = pd.DataFrame()

    def load_files(self, path: str = None, ptrn: str = None, dtfm_sig: str = None,
                   as_wkdf: bool = False, mtidx: bool = False, **kwargs):
        """Load the given file(s) into memory.

        The function will load files under given pattern (ptrn) from given path (path). It can be used more than once to handle multiple data-loading events. By default, it'll set the firstly loaded dataframe as activated one (aka self.wk_dtfm), no matter how many times it's called. However, you can specify the activated dataframe by set as_wkdf=True when calling the method, which will deactivate the last activated one and you get a warning message.  In case you don't use as_wkdf=True at all, the activated dataframe is either the first loaded dataframe or the one that specified by calling set_wkdf() which requires the signature(given by dtfm_sig) of dataframe you want to work on.

        Args:
            path (str, optional): Path to the files will be loaded. Default: None
            ptrn (str, optional): In given path, the pattern of files to be loaded. Default: None
            dtfm_sig (str, optional): The signature for new added DataFrame.  Default: None
            as_wkdf (bool, optional): Whether the current loaded files will be the working DataFrame. Default: False
            mtidx (bool, optional): Using multiple index or not. Default: False
            **kwargs: Any arguments for pandas.read_csv().
        Raises:
            TypeError:
            KeyError:
            ValueError:
        Returns:
            self: the instance itself.
        Notes:
            1. Prone to path when self.inpt_path is given; prone to prtn when self.inpt_ptrn is given at the same time
            2. For kwargs, there're some default value: sep("\t"), index_col(0), header(0).
        Todos:
        """
        if not isinstance(dtfm_sig, (str, int)):
            raise TypeError("dtfm_sig should be a str or int!!!")

        if dtfm_sig in self.dtfm_bin:
            raise KeyError("DataFrame signature already exists!!!")

        if path is None:
            if self.inpt_path is None:
                raise ValueError("When path is None, self.inpt_path should be supplied")
            path = self.inpt_path

        if "sep" not in kwargs:
            kwargs["sep"] = "\t"

        if "index_col" not in kwargs:
            kwargs["index_col"] = 0

        if "header" not in kwargs:
            kwargs["header"] = 0

        if pisdir(path):
            if ptrn is None:
                if self.inpt_ptrn is None:
                    raise ValueError(
                        "When path is a dir, either ptrn or self.inpt_ptrn should be given"
                    )
                ptrn = self.inpt_ptrn

            frames, frame_names = [], []
            files = glob.glob(pjoin(path, ptrn))
            log_me.info("Read files following patterns({}) from dir ({}). ".format(ptrn, path))

            for _file in files:
                if pisfile(_file):
                    _, flnm = os.path.split(_file)
                    dtfmnm, _ = os.path.splitext(flnm)
                    frames.append(pd.read_csv(_file, **kwargs))
                    if mtidx:
                        frame_names.append(dtfmnm)
                        self.mltp_index = True
                else:
                    log_me.warn("{} is an direcotry, skipping...".format(_file))

            if frame_names:
                rw_dtfm = pd.concat(frames, keys=frame_names)
            else:
                rw_dtfm = pd.concat(frames)
        elif pisfile(path):
            log_me.info("Read file {}".format(path))
            rw_dtfm = pd.read_csv(path, **kwargs)
        else:
            raise TypeError("Only direcotries or general file is valid")

        activated = (len(self.dtfm_bin) == 0 or as_wkdf)
        if activated:
            self.wk_dtfm = rw_dtfm
            self.dtfm_act = {sig: sig == dtfm_sig for sig in self.dtfm_act.keys()}

        self.dtfm_act[dtfm_sig] = activated
        self.dtfm_bin[dtfm_sig] = rw_dtfm

        return self

    def set_wkdf(self, dtfm_sig=None):
        """Set the given DataFrame as the actived one which will be opertated.

        The method will set the given DataFrame as the self.wk_dtfm. However, the method options
        could be both None.

        Args:
            dtfm_sig (str; optional): The signature of the DataFrame to be set as the active
                DataFrame.
        """
        if not isinstance(dtfm_sig, (str, int)):
            raise TypeError("dtfm_sig should be a str or int!!!")

        self.wk_dtfm = self.dtfm_bin.get(dtfm_sig)
        self.dtfm_act = {sig: sig == dtfm_sig for sig in self.dtfm_act.keys()}

        return self

    def remove_rows(self, rows=None):
        """Remove given rows.

        Args:
            rows (None, list, tuple; optional): [None] Rows to be removed.
        Return:
            self: the instance itself
        Notes:
            1. When removing the targeted rows, non-existing rows will be ingnored
        """
        if rows:
            if isinstance(rows, str):
                rows = [rows]
            log_me.info("Will remove {} rows: {} etc.".format(len(rows), ", ".join(rows[:5])))
            self.rmd_rows = self.wk_dtfm.reindex(index=rows)
            self.wk_dtfm.drop(index=rows, inplace=True, errors="ignore")
        else:
            log_me.info("rows is empty or None, skipping remove_rows()")

        return self

    def remove_cols(self, cols=None):
        """Remove given columns.

        Args:
            cols (list, tuple; optional): Columns to be removed. Default: None
        Return:
            self: the instance itself
        Notes:
            1. When removing the targeted columns, non-existing columns will be ingnored
        """
        if cols:
            if isinstance(cols, str):
                cols = [cols]
            log_me.info("Will remove {} columns: {} etc.".format(len(cols), ", ".join(cols[:5])))
            self.rmd_cols = self.wk_dtfm.reindex(columns=cols)
            self.wk_dtfm.drop(columns=cols, inplace=True, errors="ignore")
        else:
            log_me.info("cols is empty or None, skipping remove_cols()")

        return self

    def keep_rows(self, rows=None):
        """Only keep given rows.

        Args:
            rows (list, tuple; optional): Rows to be removed. Default: None
        Return:
            self: the instance itself
        """
        if rows:
            if isinstance(rows, str):
                rows = [rows]
            log_me.info("Will keep {} row(s): {} ...".format(len(rows), ", ".join(rows[:5])))
            self.wk_dtfm = self.wk_dtfm.loc[rows, :]
        else:
            log_me.info("rows is empty or None, skipping keep_rows()")

        return self

    def keep_cols(self, cols=None):
        """Only keep given columns.

        Args:
            cols (list, tuple; optional): Columns to be removed. Default: None
        Return:
            self: the instance itself
        """
        if cols:
            if isinstance(cols, str):
                cols = [cols]
            log_me.info("Will keep {} columns: {} ...".format(len(cols), ", ".join(cols[:5])))
            self.wk_dtfm = self.wk_dtfm.loc[:, cols]
        else:
            log_me.info("cols is empty or None, skipping keep_cols()")

        return self

    def transpose_matrix(self):
        """Transpose the matrix.

        Times of transposition will be recorded into self.tprs_counts

        Returns:
            self: the instance itself
        """
        if self.trps_n % 2 == 1:
            log_me.info("Transpose wk_dtfm, while the old wk_dtfm is transformed: {}")

        self.trps_n += 1
        self.wk_dtfm = self.wk_dtfm.transpose()

        return self

    def get_dtfm(self, dtfm_sig=None):
        """Return the pre-processed phenotype measurements as a pandas.DataFrame.

        Args:
            dtfm_sig (str; optional): Default: None
        Returns:
            pandas.DataFrame: the working dataframe (self.wk_dtfm) of current instance.
        """
        if dtfm_sig is None:
            log_me.info("Get self.wk_dtfm")
            return self.wk_dtfm
        log_me.info("Get DataFrame with signature {}".format(dtfm_sig))
        return self.dtfm_bin.get(dtfm_sig)

    def dump_wkdtfm(self, opt_dir, opt_fnm="preprocessed.tsv", **kwargs):
        """Dump the preprocessed dataframe into an given path.

        Args:
            opt_dir (str; required): The path where the working dataframe should be dumpped to.
            opt_fnm (str; optional): The name of output working dataframe. Default: preprocessed.tsv
            **kwargs: Any keyword parameters for pandas.DataFrame.to_csv()
        Raises:
            TypeError:
        Returns:
            self: the instance itself
        """
        if not pisdir(opt_dir):
            raise TypeError("Argument opt should be a dir!")

        if not os.path.exists(opt_dir):
            os.makedirs(opt_dir)

        if "sep" not in kwargs:
            kwargs["sep"] = "\t"

        opt_path = pjoin(opt_dir, opt_fnm)
        self.wk_dtfm.to_csv(opt_path, index=True, header=True, **kwargs)

        return self

    def calc_corr(self, target=None, method="spearman"):
        """Caluculate the correlation between phenotypes (Pairwise)

        Args:
            target (list, tuple; optional): All the combination or subset of all phenotypes, the
                later can be specified by and list.
            method (str, callable; optional): The method use to calculate the correlation.
                Options: [spearman, pearson, kendall]. Read more at pandas.DataFrame.corr().
        Raises:
            TypeError: When the target augument given value type beyound None, tuple, and list.
        Returns:
            self: The instance itself
        """
        if method in ["spearman", "pearson", "kendall"]:
            log_me.info("Will use {}".format(method))
        elif callable(method):
            log_me.info("User custom callable {}".format(method))
        else:
            raise TypeError("Unknown correlation methods")

        if target is None:
            self.pw_corr = self.wk_dtfm.corr(method=method)
        elif isinstance(target, (tuple, list)):
            self.pw_corr = self.wk_dtfm.loc[:, target].corr(method=method)
        else:
            raise TypeError(
                "target should be a list or tuple, otherwise leave target as None using all"
            )

        return self

    def concat_dtfm(self, signatures=None, **kwargs):
        """Concatenate multiple dataframe in dtfm_bin

        Args:
            signatures (list, str; optional): Signatures of datafram will be concatenated.

        Returns:
            pandas.DataFrame
        """
        if signatures is None:
            dtfm_list = self.dtfm_bin.values()
        else:
            dtfm_list = [dtfm for key, dtfm in self.dtfm_bin.items() if key in signatures]

        return pd.concat(dtfm_list, **kwargs)

    @staticmethod
    def _wrap_into_list(target):
        """An useful method to check whether the target is legal

        Args:
            target (str, list, tuple):
        Raises:
            TypeError: The type of parameter target should be str, list or tuple
        returns:
            target (list, tuple)
        """
        if isinstance(target, str):
            target = [target]
        elif not isinstance(target, (list, tuple)):
            raise TypeError("target should be str, list or tuple")

        return target

