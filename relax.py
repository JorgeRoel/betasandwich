import sys
from filters import check_long_hairpins, check_short_hairpins
from protocols import relax
from pyrosetta import *

pyrosetta.init(silent=True)

relax(sys.argv[1])
