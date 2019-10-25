# -*- UTF-8 -*-
"""A class to preprocessing covariates
"""

#
## covariat.py
#
import os

from preproc import PreProc

class Covariate(PreProc):
    """A class handling covariates.

    Currently, just a wrapper of `PrePorc` as it can cover all need fucntions.
    """

    def __init__(self, wk_dir, inpt_path, inpt_ptrn="*.tsv"):
        """"""
        self.wk_dir = wk_dir
        super(Covariate, self).__init__(inpt_path, inpt_ptrn)

    def dump_cvrt(self, opt_dir=None, opt_fnm="covariates.tsv"):
        """Dump processed covariates into disk.
        Args:
            opt_dir (str; optional): The output directory.
            opt_fnm (str; optional): The output file name.
        Returns:
            self
        """

        if opt_dir is None:
            opt_dir = os.path.join(self.wk_dir, "Preprocessing")

        self.dump_wkdtfm(opt_dir, opt_fnm)

        return self
