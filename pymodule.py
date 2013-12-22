import psi4
import re
import os
import inputparser
import math
import warnings
from driver import *
from wrappers import *
from molutil import *
import p4util
from psiexceptions import *


def run_psi4_gamess_mos(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    mos can be called via :py:func:`~driver.energy`. For post-scf plugins.

    >>> energy('psi4_gamess_mos')

    """
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

    # Your plugin's psi4 run sequence goes here
    psi4.set_global_option('BASIS', 'sto-3g')
    psi4.set_local_option('PSI4_GAMESS_MOS', 'PRINT', 1)
    scf_helper(name, **kwargs)
    returnvalue = psi4.plugin('psi4_gamess_mos.so')
    psi4.set_variable('CURRENT ENERGY', returnvalue)


# Integration with driver routines
procedures['energy']['psi4_gamess_mos'] = run_psi4_gamess_mos


def exampleFN():
    # Your Python code goes here
    pass
