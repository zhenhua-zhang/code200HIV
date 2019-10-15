# -*- UTF-8 -*-

#
## Preprocessing super class. 
#
import os
import sys
import copy
import glob
import pandas as pd
from utls import log_me
from utls import ArgumentDependencyError

class PreProc(object):
    """A class to preprocess raw data
    """
    def __init__(self, input_path, input_ptrn="*.tsv"):
        """A class to preprocess data.

        Args:
            input_path (str): required The path of input file of dir
            input_ptrn (str): *.tsv Pattern to search input file
        """
        self.input_path = input_path
        self.input_ptrn = input_ptrn
        self.rw_dtfm = None
        self.wk_dtfm = None
        self.rmd_cols = None
        self.rmd_rows = None
        self.mltp_index = False
        self.trps_counts = 0
    
    def load_files(self, keep_rwdt=False, ptrn=None, mtidx=False, **kwargs):
        """Load the given file into memory.

        Args:
            keep_raw_data (bool): [False] If it's True it will deep copy the raw_dataframe in the instance, whch is memory-cost.
            ptrn (str, None): [None]
            mtidx (bool): [False] Using multiple index or not.
            **kwargs: Any arguments for `pandas.read_csv()`.

        Raises: None

        Returns:
            self: the instance itself
        """

        if os.path.isdir(self.input_path):
            if self.input_ptrn is None:
                log_me.error("if `input_path` is a dir, `input_ptrn` is required")
                raise ArgumentDependencyError

            frames, frame_names = [], []
            files = glob.glob(os.path.join(self.input_path + self.input_ptrn))
            log_me.info("Read from directory: {}".format(self.input_path))

            if len(files) > 1:
                log_me.info("Files following pattern: {}".format(self.input_ptrn))
            else:
                log_me.info("But only one file following patern: {}".format(self.input_ptrn))

            for _file in files:
                if os.path.isfile(_file):
                    _, flnm = os.path.split(_file)
                    dtfmnm, _ = os.path.splitext(flnm)
                    frames.append(pd.read_csv(_file, header=0, index_col=0, sep="\t" **kwargs))
                    if mtidx:
                        frame_names.append(dtfmnm)
                        self.mltp_index = True
                else:
                    log_me.warn("{} is an direcotry, skipping...".format(_file))

            self.rw_dtfm = pd.concat(frames, keys=frame_names)

        elif os.path.isfile(self.input_path):
            self.rw_dtfm = pd.read_csv(self.input_path, header=0, index_col=0,  sep="\t", **kwargs)
        else:
            raise TypeError("Only direcotries or general file is valide")

        if  keep_rwdt:
            log_me.info("Will keep raw data frame in memory")
            self.wk_dtfm = copy.deepcopy(self.rw_dtfm)
        else:
            log_me.info("Will NOT keep raw data frame")
            self.wk_dtfm = self.rw_dtfm
        
        return self

    def remove_rows(self, rows=None):
        """Remove given rows.

        Args:
            rows (None, list, tuple): [None] Rows to be removed.

        Raises: None

        Return:
            self: the instance itself
        """

        if rows:
            if isinstance(rows, str):
                rows = [rows]
            log_me.info("Will remove: {} ...".format(", ".join(rows[: min(5, len(rows))])))
            self.rmd_rows = rows
            self.wk_dtfm = self.wk_dtfm.drop(index=rows)

        return self
    
    def remove_cols(self, cols=None):
        """Remove given columns.

        Args:
            cols (None, list, tuple): [None] Columns to be removed.

        Raises: None

        Return:
            self: the instance itself
        """

        if cols:
            if isinstance(cols, str):
                cols = [cols]
            log_me.info("Will remove: {} ...".format(", ".join(cols[: min(5, len(cols))])))
            self.rmd_cols = cols
            self.wk_dtfm = self.wk_dtfm.drop(columns=cols)

        return self

    def transpose_matrix(self):
        """Transpose the matrix.
        Times of transposition will be recorded into self.tprs_counts

        Args: None

        Raises: None

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

    def get_wkdtfm(self):
        """Return the pre-processed phenotype measurements as a pandas.DataFrame.

        Args: None

        Raises: None

        Return:
            pandas.DataFrame: the working dataframe (self.wk_dtfm) of current instance.
        """
        return self.wk_dtfm
    
    def dump_wkdtfm(self, opt, flnm="preprocessed_phenotypes.tsv"):
        """Dump the preprocessed dataframe into an given path.

        Args:
            opt (str): required. The path where the working dataframe should be dumpped to.
            flnm (str): [preprocessed_phenotypes.tsv] The name of output working dataframe
        Raises: None
        Returns:
            self: the instance itself
        """
        if not os.path.isdir(opt):
            raise TypeError("Argument `opt` should be a dir!")
        elif not os.path.exists(opt):
            os.makedirs(opt)

        opt_path = os.path.join(opt, flnm)
        if os.path.exists(opt_path):
            log_me.warn("File: {} exists, will overwrite it!".format(opt_path))
        self.wk_dtfm.to_csv(opt_path, sep="\t", index=True, header=True)

        return self

    # def summarize(self):
    #     """"""