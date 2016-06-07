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

F1ref = [-5.43569555885951, -1.88298017648151, -2.17427822354382, \
          3.47022732532752, -2.96798871197376,  2.59189820356784, \
          1.96546823353199,  4.85096888845526, -0.417619980024024 ]
F1tol = 1e-9
F1ok = False if LIBINT_ONEBODY_DERIV>0 else True

F2ref = [ 2.95100851676078,  1.02225933691908,  1.18040340670247, \
         -1.89116113658198,  1.64868682622021, -1.42151545975317, \
         -1.0598473801788 , -2.67094616313929,  0.241112053050697 ]
F2tol = 1e-9
F2ok = False if LIBINT_ERI_DERIV>0 else True

N1ref = [ 2.01500148345464 ,  0.698016989334019,  0.806000593381856, \
         -1.28187870177    ,  1.07670120714609 , -0.951756216776298, \
         -0.733122781684634, -1.77471819648011 ,  0.145755623394443 ]
N1tol = 1e-10
N1ok = False if LIBINT_ONEBODY_DERIV>0 and LIBINT_ERI_DERIV>0 else True

Fref = [ -0.114420248668043 , -0.039636336821566 , -0.0457680994672225, \
          0.0729288451097068, -0.0618587006295552,  0.0543214912849466, \
          0.0414914035583673,  0.101495037451119 , -0.00855339181771989 ]
Ftol = 1e-9
Fok = False if LIBINT_ONEBODY_DERIV>0 and LIBINT_ERI_DERIV>0 else True

for line in sys.stdin:
    match1 = re.match('\*\* Hartree-Fock energy =' + pat_numbers(1), line)
    match2 = re.match('\*\* edipole =' + pat_numbers(3), line)
    match3 = re.match('\*\* equadrupole =' + pat_numbers(6), line)
    match4 = re.match('\*\* 1-body forces =' + pat_numbers(9), line)
    match5 = re.match('\*\* Pulay forces =' + pat_numbers(9), line)
    match6 = re.match('\*\* 2-body forces =' + pat_numbers(9), line)
    match7 = re.match('\*\* nuclear repulsion forces =' + pat_numbers(9), line)
    match8 = re.match('\*\* Hartree-Fock forces =' + pat_numbers(9), line)
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
    elif match6:
        F2ok = validate("2-body force", match6.groups(), F2ref, F2tol, line)
    elif match7:
        N1ok = validate("nuclear repulsion force", match7.groups(), N1ref, N1tol, line)
    elif match8:
        Fok = validate("HF force", match8.groups(), Fref, Ftol, line)
    else:
        print(line,end="")

ok = eok and muok and Qok and F1ok and FPok and F2ok and N1ok and Fok
if not ok: sys.exit(1)


