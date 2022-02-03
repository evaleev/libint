from __future__ import print_function
import math
import re
import sys

# first arg is the input file, if missing read sys.stdin
try:
    instr = open(sys.argv[1]) if len(sys.argv) > 1 else sys.stdin
except Exception as e:
    sys.exit(2)


eref = -74.942080057696
tol = 1e-11

returncode = 1
for line in instr:
    print(line,end="")
    x = re.match(r'\*\* Hartree-Fock energy =\s*([-\d.]+)', line)
    if x:
        efound_str = x.group(1)
        if (math.fabs(eref - float(efound_str)) < tol):
            returncode = 0
            print("Hartree-Fock energy check: passed")

if len(sys.argv) > 1:
    instr.close()

sys.exit(returncode)
