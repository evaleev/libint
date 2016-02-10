from __future__ import print_function
import sys, re, math, os

def pat_numbers(n):
    result = ''
    for i in range(n):
        result += '\s*([+-e\d.]+)'
    return result

def validate(label, data, refdata, tolerance, textline):
    ok = True
    ndata = len(refdata)
    for i in range(ndata):
        datum = float(data[i])
        refdatum = refdata[i]
        if (math.fabs(refdatum - datum) > tolerance):
            ok = False
            print(label, "check: failed\nreference:", refdata, "\nactual:", textline)
            break
    if (ok): print(label, "check: passed")
    return ok

path_to_libfeatures = sys.argv[1]
if os.path.exists(path_to_libfeatures): 
    # define boolean constants used in MakeVars.features
    no = False
    yes = True
    execfile(path_to_libfeatures)

eref = [-76.003354058456]
etol = 5e-12

muref = [-0.263282355852899, -0.0912036834147694, -0.105312942341114]
mutol = 1e-9

Qref = [-7.01719638259305, 0.0246243455752736, -0.576201298992812,
        -8.04716669233508, 0.638204908066147,  -6.35134519438538]
Qtol = 1e-8

F1ref = [-5.43569555885951, -1.88298017648151, -2.17427822354382, \
          3.47022732532752, -2.96798871197376,  2.59189820356784, \
          1.96546823353199,  4.85096888845526, -0.417619980024024 ]
F1tol = 1e-9
F1ok = False if LIBINT_ONEBODY_DERIV>0 else True

FPref = [ 0.35526531023214,   0.123067513497363,  0.14210612409289, \
         -0.224258641960886,  0.180741977756425, -0.164305035736127, \
         -0.131006668271254, -0.303809491253788,  0.0221989116432378]
FPtol = 1e-9
FPok = False if LIBINT_ONEBODY_DERIV>0 else True

for line in sys.stdin:
    match1 = re.match('\*\* Hartree-Fock energy =' + pat_numbers(1), line)
    match2 = re.match('\*\* edipole =' + pat_numbers(3), line)
    match3 = re.match('\*\* equadrupole =' + pat_numbers(6), line)
    match4 = re.match('\*\* 1-body forces =' + pat_numbers(9), line)
    match5 = re.match('\*\* Pulay forces =' + pat_numbers(9), line)
    if match1:
        eok = validate("HF energy", match1.groups(), eref, etol, line)
    elif match2:
        muok = validate("electric dipole moment", match2.groups(), muref, mutol, line)
    elif match3:
        Qok = validate("electric quadrupole moment", match3.groups(), Qref, Qtol, line)
    elif match4:
        F1ok = validate("1-body force", match4.groups(), F1ref, F1tol, line)
    elif match5:
        FPok = validate("Pulay force", match5.groups(), FPref, FPtol, line)
    else:
        print(line,end="")

ok = eok and muok and Qok and F1ok and FPok
if not ok: sys.exit(1)


