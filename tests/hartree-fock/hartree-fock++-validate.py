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
        -8.04716669233508, 0.638204908066147,  -6.35134519438538  ]
Qtol = 1e-8

F1ref = [-5.43569555885951, -1.88298017648151, -2.17427822354382,
          3.47022732532752, -2.96798871197376,  2.59189820356784,
          1.96546823353199,  4.85096888845526, -0.417619980024024 ]
F1tol = 1e-9
F1ok = False if LIBINT_ONEBODY_DERIV>0 else True

FPref = [ 0.35526531023214,   0.123067513497363,  0.14210612409289,
         -0.224258641960886,  0.180741977756425, -0.164305035736127,
         -0.131006668271254, -0.303809491253788,  0.0221989116432378 ]
FPtol = 1e-9
FPok = False if LIBINT_ONEBODY_DERIV>0 else True

F2ref = [ 2.95100851676078,  1.02225933691908,  1.18040340670247,
         -1.89116113658198,  1.64868682622021, -1.42151545975317,
         -1.0598473801788 , -2.67094616313929,  0.241112053050697 ]
F2tol = 1e-9
F2ok = False if LIBINT_ERI_DERIV>0 else True

N1ref = [ 2.01500148345464 ,  0.698016989334019,  0.806000593381856,
         -1.28187870177    ,  1.07670120714609 , -0.951756216776298,
         -0.733122781684634, -1.77471819648011 ,  0.145755623394443  ]
N1tol = 1e-10
N1ok = False if LIBINT_ONEBODY_DERIV>0 and LIBINT_ERI_DERIV>0 else True

Fref = [ -0.114420248668043 , -0.039636336821566 , -0.0457680994672225,
          0.0729288451097068, -0.0618587006295552,  0.0543214912849466,
          0.0414914035583673,  0.101495037451119 , -0.00855339181771989 ]
Ftol = 1e-9
Fok = False if LIBINT_ONEBODY_DERIV>0 and LIBINT_ERI_DERIV>0 else True

H1ref = [3.54760842854909, 0.0561425790505338, -1.8408803442899, -0.887601426469258, -2.17549130380983, 1.98344832770147, -2.66000700208114, 2.1193487247593, -0.142567983411562, 0.355853985911328, 2.00292709865962, -2.17850638724121, -1.77309201100161, -1.61619662190203, 2.12236380819067, 1.41723802509053, -0.386730476757594, 5.69831976206994, 1.98431870731686, -1.61468908018631, -2.14019765079036, -0.143438363026958, -0.388238018473309, -3.55812211127938, 0.9042631343743, 2.18335388767572, -1.98924643245079, -0.0166617079050424, -0.00484750043451748, 0.00492772513392813, 1.7532880067497, 1.63118422317771, -0.00786258386589356, 0.0198040042519126, -0.0164951429913965, 2.15669050309239, 0.00579810474932402, -0.0149876012756808, -0.016492852302032, 2.67666870998618, -2.11450122432478, 0.137640258277634, -1.43704202934245, 0.40322561974899, 3.57461496358141 ]
H1tol = 1e-9
H1ok = False if LIBINT_ONEBODY_DERIV>1 else True

HPref = [ -0.373657119720633, -0.00133454389437712, 0.0913757616940999, 0.144293176043723, 0.103810096192085, -0.0967303414264739, 0.22936394367691, -0.102475552297708, 0.00535457973237397, -0.221600106506835, -0.0970957473316287, 0.103810096192085, 0.187363744122523, 0.0780172520214974, -0.102475552297708, 0.0342363623844577, 0.0190784953101312, -0.481921135055919, -0.0967303414264739, 0.0780172520214974, 0.206932260474631, 0.00535457973237397, 0.0190784953101312, 0.274988874581143, -0.141868433967373, -0.102526043934845, 0.0963037328217154, -0.00242474207634997, -0.0012840522572401, 0.000426608604758467, -0.191255819978244, -0.0755194237428591, -0.0012840522572401, 0.00389207585572005, -0.00249782827863833, -0.205330014322054, 0.000426608604758466, -0.00249782827863833, -0.00160224615257662, -0.22693920160056, 0.103759604554948, -0.00578118833713244, -0.0381284382401776, -0.0165806670314929, -0.273386628428567 ]
HPtol = 1e-9
HPok = False if LIBINT_ONEBODY_DERIV>1 else True

for line in sys.stdin:
    match1 = re.match('\*\* Hartree-Fock energy =' + pat_numbers(1), line)
    match2 = re.match('\*\* edipole =' + pat_numbers(3), line)
    match3 = re.match('\*\* equadrupole =' + pat_numbers(6), line)
    match4 = re.match('\*\* 1-body forces =' + pat_numbers(9), line)
    match5 = re.match('\*\* Pulay forces =' + pat_numbers(9), line)
    match6 = re.match('\*\* 2-body forces =' + pat_numbers(9), line)
    match7 = re.match('\*\* nuclear repulsion forces =' + pat_numbers(9), line)
    match8 = re.match('\*\* Hartree-Fock forces =' + pat_numbers(9), line)
    match9 = re.match('\*\* 1-body hessian =' + pat_numbers(45), line)
    match10 = re.match('\*\* Pulay hessian =' + pat_numbers(45), line)
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
    elif match9:
        H1ok = validate("1-body hessian", match9.groups(), H1ref, H1tol, line)
    elif match10:
        HPok = validate("Pulay hessian", match10.groups(), HPref, HPtol, line)
    else:
        print(line,end="")

ok = eok and muok and Qok and F1ok and FPok and F2ok and N1ok and Fok and H1ok and HPok
if not ok: sys.exit(1)


