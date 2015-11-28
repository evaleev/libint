import sys, re, math, os

def pat_numbers(n):
    result = ''
    for i in range(n):
        result += '\s*([+-e\d.]+)'
    return result

path_to_libfeatures = sys.argv[1]
if os.path.exists(path_to_libfeatures): 
    # define boolean constants used in MakeVars.features
    no = False
    yes = True
    execfile(path_to_libfeatures)

eref = -76.003354058456
etol = 1e-11
eok = False

muref = [-0.263282355181168, -0.091203683186213, -0.105312942071265]
mutol = 1e-10
muok = False

Qref = [-7.01719638095415, 0.0246243454339316, -0.57620129865407, \
        -8.04716668958708, 0.638204907503245,  -6.35134519301054]
Qtol = 1e-10
Qok = False

F1ref = [-5.4356955590748,  -1.88298017661221, -2.17427822361541, \
          3.47022732552768, -2.96798871234053,  2.59189820377569, \
          1.96546823354713,  4.85096888895273, -0.417619980160289]
F1tol = 1e-10
F1ok = False if LIBINT_ONEBODY_DERIV>0 else True

FPref = [ 0.355265310364197,  0.123067513554652,  0.142106124142559, \
         -0.224258642051479,  0.180741977858114, -0.16430503581098,  \
         -0.131006668312719, -0.303809491412766,  0.0221989116684212]
FPtol = 1e-10
FPok = False if LIBINT_ONEBODY_DERIV>0 else True

for line in sys.stdin:
    match1 = re.match('\*\* Hartree-Fock energy =' + pat_numbers(1), line)
    match2 = re.match('\*\* edipole =' + pat_numbers(3), line)
    match3 = re.match('\*\* equadrupole =' + pat_numbers(6), line)
    match4 = re.match('\*\* 1-body forces =' + pat_numbers(9), line)
    match5 = re.match('\*\* Pulay forces =' + pat_numbers(9), line)
    if match1:
        efound = float(match1.group(1)) 
        if (math.fabs(eref - efound) < etol):
            eok = True
            print "HF energy check: passed"
        else:
            print "HF energy check: failed\n", line
    elif match2:
        for i in range(3):
            mufound = float(match2.group(i+1))
            if (math.fabs(muref[i] - mufound) < mutol):
                muok = True
            else:
                print "electric dipole moment check: failed\n", line
        if (muok): print "electric dipole moment check: passed"
    elif match3:
        for i in range(6):
            Qfound = float(match3.group(i+1))
            if (math.fabs(Qref[i] - Qfound) < Qtol):
                Qok = True
            else:
                print "electric quadrupole moment check: failed\n", line
        if (Qok): print "electric quadrupole moment check: passed"
    elif match4:
        for i in range(9):
            F1found = float(match4.group(i+1))
            if (math.fabs(F1ref[i] - F1found) < F1tol):
                F1ok = True
            else:
                print "1-body force check: failed\n", line
        if (F1ok): print "1-body force check: passed"
    elif match5:
        for i in range(9):
            FPfound = float(match5.group(i+1))
            if (math.fabs(FPref[i] - FPfound) < FPtol):
                FPok = True
            else:
                print "Pulay force check: failed\n", line
        if (FPok): print "Pulay force check: passed"
    else:
        print(line),

ok = eok and muok and Qok and F1ok and FPok
if not ok: sys.exit(1)


