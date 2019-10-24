# -*- UTF-8 -*-

#
## preproc.py
#

import os
import sys
import copy
import glob
from os.path import join as pjoin
from os.path import isdir as pisdir
from os.path import isfile as pisfile

import pandas as pd

# Customer scripts in current package
from utls import log_me
from utls import ArgumentDependencyError

class PreProc(object):
    """A class to preprocess raw data.

    This class is inheriten by phenotype.Phenotype, genotype.Genotype,
    dovariate.Covariate in qtlmapping package.

    Attributes:
    Methods:
    Todos:
        1. Complete the documents for each methods and the class
        2. A method to summarize the working dataframe
        3. A container to record operations `self.operations`
    """
    def __init__(self, inpt_path="./", inpt_ptrn ="*.tsv"):
        """A class to preprocess data.

        Args:
            inpt_path (str): The path of input file of dir. Required.
            inpt_ptrn (str, optional): Pattern to search input file. Default: *.tsv
        """
        self.inpt_path = inpt_path
        self.inpt_ptrn = inpt_ptrn
        self.dtfm_bin = {}         # A bin of loaded files
        self.dtfm_idx = []         # Index for each load_files operation
        self.dtfm_sig = []
        self.dtfm_act = []         # Tracks whether the dataframe is activated
        self.wk_dtfm = pd.DataFrame()
        self.wk_dtfm_id = None
        self.rmd_cols = pd.DataFrame()
        self.rmd_rows = pd.DataFrame()
        self.mltp_index = False
        self.trps_counts = 0

    def load_files(self, path=None, ptrn=None, as_wkdf=False, mtidx=False,
                   dtfm_sig=None, **kwargs):
        """Load the given file into memory.

        The function will load files under given pattern (`ptrn`) from given
        path (`path`). It can be used more than once to handle multiple
        data-loading events. By default, it'll set the firstly loaded dataframe
        as activated one (aka `self.wk_dtfm`), no matter how many times it's
        called. However, you can specify the activated dataframe by set
        `as_wkdf=True` when calling the method, which will deactivate the
        last activated one and you get a warning message.  In case you don't
        use `as_wkdf=True` at all, the activated dataframe is either the first
        loaded dataframe or the one that specified by calling `set_wkdf()`
        which requires the `index` of dataframe you want to work on.

        Args:
            `path` (`str`, optional): Path to the files will be loaded. Default: `None`
            `ptrn` (`str`, optional): In given path, the pattern of files to be loaded. Default: `None`
            `as_wkdf` (`bool`, optional): Whether the current loaded files will be the working DataFrame. Default: `False`
            `mtidx` (`bool`, optional): Using multiple index or not. Default: `False`
            `**kwargs`: Any arguments for `pandas.read_csv()`.
        Raises: 
            `ValueError`:
        Returns:
            self: the instance itself
        Todos:
            -[x] Deal with the conflicts of input file patterns between `__init__()` and `load_files()`
        """
        if path is None and self.inpt_path is None:
            raise ValueError("When `path` is None, `self.inpt_path` should be supplied")
        elif path is None:
            path = self.inpt_path

        if ptrn is None and self.inpt_ptrn is None:
            raise ValueError("When `ptrn` is None, `self.inpt_ptrn` should be supplied")
        elif ptrn is None:
            ptrn = self.inpt_ptrn
        
        if "sep" not in kwargs:
            kwargs["sep"] = "\t"

        if "index_col" not in kwargs:
            kwargs["index_col"] = 0

        if "header" not in kwargs:
            kwargs["header"] = 0

        if pisdir(path):
            if ptrn is None:
                raise ValueError("if `path` is a dir, `ptrn` is required")

            frames, frame_names = [], []
            files = glob.glob(pjoin(path, ptrn))
            log_me.info("Read from directory: {}".format(path))
            log_me.info("Will find files following patterns: {}".format(ptrn))

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
            rw_dtfm = pd.read_csv(path, **kwargs)
        else:
            raise TypeError("Only direcotries or general file is valid")

        n_idx = len(self.dtfm_idx)
        activated = (n_idx == 0 or as_wkdf == True)

        if activated == True:
            self.wk_dtfm = rw_dtfm
            self.wk_dtfm_id = (n_idx, dtfm_sig)
            self.dtfm_act = [False] * n_idx

        self.dtfm_idx.append(n_idx)
        self.dtfm_sig.append(dtfm_sig)
        self.dtfm_act.append(activated)
        self.dtfm_bin[(n_idx, dtfm_sig)] = rw_dtfm

        return self

    def set_wkdf(self, index=None, dtfm_sig=None):
        """Set the given DataFrame as the actived one which will be opertated.

        The method will set the given DataFrame as the `self.wk_dtfm`. However,
        the method options could be both None.

        Args:
            index (int; optional): The index of the DataFrame to be set as the active DataFrame.
            dtfm_sig (str; optional): The signature of the DataFrame to be set as the active DataFrame.
        """
        if index is None:
            index = self.dtfm_act.index(True)

        if dtfm_sig is None:
            dtfm_sig = self.dtfm_sig[index]

        self.wk_dtfm = self.dtfm_bin.get((index, dtfm_sig), None)
        self.wk_dtfm_id = (index, dtfm_sig)

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
            log_me.info("Will remove {} rows: {} etc.".format(len(rows), ", ".join(rows[: min(5, len(rows))])))
            self.rmd_rows = self.wk_dtfm.reindex(index=rows)
            self.wk_dtfm.drop(index=rows, inplace=True, errors="ignore")
        else:
            log_me.info("`rows` is empty or None, skipping `remove_rows()`")

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
            log_me.info("Will remove {} columns: {} etc.".format(len(cols), ", ".join(cols[: min(5, len(cols))])))
            self.rmd_cols = self.wk_dtfm.reindex(columns=cols)
            self.wk_dtfm.drop(columns=cols, inplace=True, errors="ignore")
        else:
            log_me.info("`cols` is empty or None, skipping `remove_cols()`")

        return self

    def keep_rows(self, rows=None):
        """Only keep given rows.

        Args:
            rows (None, list, tuple; optional): [None] Rows to be removed.
        Return:
            self: the instance itself
        """
        if rows:
            if isinstance(rows, str):
                rows = [rows]
            log_me.info("Will keep {} row(s): {} ...".format(len(rows), ", ".join(rows[: min(5, len(rows))])))
            self.kpt_rows = rows
            self.wk_dtfm = self.wk_dtfm.loc[rows, :]
        else:
            log_me.info("`rows` is empty or None, skipping `keep_rows()`")

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
            log_me.info("Will keep {} columns: {} ...".format(len(cols), ", ".join(cols[: min(5, len(cols))])))
            self.kpt_cols = cols
            self.wk_dtfm = self.wk_dtfm.loc[:, cols]
        else:
            log_me.info("`cols` is empty or None, skipping `keep_cols()`")

        return self

    def transpose_matrix(self):
        """Transpose the matrix.
        Times of transposition will be recorded into self.tprs_counts

        Return:
            self: the instance itself
        """
        log_me.info(
            "Will transpose the dataframe, " +
            "while the old is transformed: {}".format(self.trps_counts % 2 == 1)
        )

        self.trps_counts += 1
        self.wk_dtfm = self.wk_dtfm.transpose()

        return self

    def get_dtfm(self, index=None, dtfm_sig=None):
        """Return the pre-processed phenotype measurements as a pandas.DataFrame.

        Args:
            index (int; optional): Default: None
            dtfm_sig (str; optional): Default: None
        Return:
            pandas.DataFrame: the working dataframe (self.wk_dtfm) of current instance.
        """
        if index is None:
            return self.wk_dtfm
        elif dtfm_sig is None:
            return self.dtfm_bin.get((index, None), None)
        else:
            return self.dtfm_bin.get((index, dtfm_sig), None)

    def dump_wkdtfm(self, opt, flnm="preprocessed.tsv", discard_wkdtfm=False, reset_wkdtfm=False, **kwargs):
        """Dump the preprocessed dataframe into an given path.

        Args:
            opt (str): required. The path where the working dataframe should be dumpped to.
            flnm (str; optional): [preprocessed.tsv] The name of output working dataframe.
            discard_wkdtfm (bool; optional): Pop out wk_dtfm after dumpping it to the disk which helps reduce resource usage. Default: False
            **kwargs: Any keyword parameters for `pandas.DataFrame.to_csv()`
        Raises:
            TypeError:
        Returns:
            self: the instance itself
        """
        if not pisdir(opt):
            raise TypeError("Argument `opt` should be a dir!")

        if not os.path.exists(opt):
            os.makedirs(opt)

        if "sep" not in kwargs:
            kwargs["sep"] = "\t"

        opt_path = pjoin(opt, flnm)
        if self.wk_dtfm.index.name is None:
            self.wk_dtfm.index.name = "id"
        self.wk_dtfm.to_csv(opt_path, index=True, header=True, **kwargs)

        if discard_wkdtfm:
            self.dtfm_bin.pop(self.wk_dtfm_id)

        if reset_wkdtfm:
            rest_keys = list(self.dtfm_bin.keys())
            if len(rest_keys) != 0:
                self.wk_dtfm_id = rest_keys[0]
                self.wk_dtfm = self.dtfm_bin.get(self.wk_dtfm_id)
            else:
                log_me.warn("The working dataframe bin is empty now!")

        return self

    def calc_corr(self, target=None, method="spearman", dist_mthd="ward"):
        """Caluculate the correlation between phenotypes (Pairwise).

        Args:
            target (list, tuple, None): [None] All the combination or subset of all phenotypes, the later can be specified by and list
            method (str, callable): ['spearman'] The method use to calculate the correlation. Options: [spearman, pearson, kendall]. For more info, refer to `pandas.DataFrame.coor()`
            dist_mthd (str): ['ward'] Method used to calculate the distance between correlations.
        Raises:
            TypeError: When the target augument given value type beyound None, tuple, and list
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
            self.pw_corr = self.wk_dtfm.corr(method=method)  # Pairwise
        elif isinstance(target, (tuple, list)):
            self.pw_corr = self.wk_dtfm.loc[:, target].corr(method=method)
        else:
            raise TypeError("`target` should be a list or tuple, otherwise leave target as None using all")
        
        return self

    @staticmethod
    def _wrap_into_list(target):
        """An useful method to check whether the target is legal

        Args:
            target (str, list, tuple):
        Raises:
            TypeError: The type of parameter `target` should be str, list or tuple
        returns:
            target (list, tuple)
        """
        if isinstance(target, str):
            target = [target]
        elif not isinstance(target, (list, tuple)):
            raise TypeError("`target` should be str, list or tuple")

        return target

    # def summarize(self):
    #     """"""
