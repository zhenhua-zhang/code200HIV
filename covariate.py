# -*- UTF-8 -*-
"""A class to preprocessing covariates
"""

#
## covariat.py
#

from preproc import PreProc

class Covariate(PreProc):
    """A class handling covariates.

    Currently, just a wrapper of `PrePorc` as it can cover all need fucntions.
    """

    def __init__(self, inpt_path, inpt_ptrn="*.tsv"):
        """"""
        super(Covariate, self).__init__(inpt_path, inpt_ptrn)
