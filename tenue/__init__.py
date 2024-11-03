from tenue.cook import (
    readbias,
    makebias,
    readdark,
    makedark,
    readflat,
    readmask,
    makeflatandmask,
)
from tenue.object import makeobject

import sys

if not sys.warnoptions:
    import warnings

    warnings.simplefilter("ignore")
