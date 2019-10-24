# -*- UTF-8 -*-

#
## py_matrix_qtl.py
#

"""A module of Python mirror to impliment R MatrixEQTL package

The module includes two Class: `MatrixQTL` and `SlicedData`. The `MatrixQTL` is
the main class, while the `SlicedData` is the main data structure.
"""

from preproc import PreProc
from utls import not_runnable, not_implemented

class MatrixQTL(PreProc):
    """A class to do QTL mapping forked from MatrixEQTL

    Attributes:
    Methods:
    InherenteFrom:
    InherentedBy:
    Todos:
    Notes:
    References:
    """

    def __init__(self, gntp, pntp, cvrt):
        """Create a MatrixQTL instance"""
        not_implemented()

    def set_params(self, **kwargs):
        """Set given params"""
        not_implemented()

    def map():
        """Main function to run the mapping"""
        not_implemented()


class SlicedData:
    """A Python mirror of SlicedData object in MatrixEQTL

    Attributes:
    Methods:
    InherenteFrom:
    InherentedBy:
    Todos:
    Notes:
    References:
    """

    def __init__(self):
        """Create a SlicedData instance"""
        not_implemented()


if __name__ == "__main__":
    not_runnable()
