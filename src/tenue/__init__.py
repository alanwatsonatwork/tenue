from tenue.cook import (
    readbias,
    makebias,
    readdark,
    makedark,
    readflat,
    makeflat,
    usefakebias,
    usefakedark,
    usefakeflat,
)
from tenue.object import makeobject

import sys

if not sys.warnoptions:
    import warnings

    warnings.simplefilter("ignore", UserWarning)
