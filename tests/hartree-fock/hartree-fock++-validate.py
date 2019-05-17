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

assert(len(sys.argv) > 1)
# first arg is the path to MakeVars.features
path_to_libfeatures = sys.argv[1]
# second arg is the input file, if missing read sys.stdin
try:
    instr = open(sys.argv[2]) if len(sys.argv) > 2 else sys.stdin
except Exception as e:
    sys.exit(2)

if os.path.exists(path_to_libfeatures):
    # define boolean constants used in MakeVars.features
    no = False
    yes = True
    exec(open(path_to_libfeatures).read())

eref = [-76.003354058454]
etol = 5e-12

muref = [-0.263282355852899, -0.0912036834147694, -0.105312942341114]
mutol = 1e-9

Qref = [-7.01719638259305, 0.0246243455752736, -0.576201298992812,
        -8.04716669233508, 0.638204908066147,  -6.35134519438538  ]
Qtol = 1e-8

smuref = [-muref[1]/2, muref[2], -muref[0]/2]
smutol = 1e-9

sQref = [Qref[1]/4, -Qref[4]/2, (2*Qref[5] - Qref[0] - Qref[3])/4,
         -Qref[2]/2, (Qref[0] - Qref[3])/8  ]
sQtol = 1e-8

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

FNref = [ 2.01500148345464 ,  0.698016989334019,  0.806000593381856,
         -1.28187870177    ,  1.07670120714609 , -0.951756216776298,
         -0.733122781684634, -1.77471819648011 ,  0.145755623394443  ]
FNtol = 1e-10
FNok = False if LIBINT_ONEBODY_DERIV>0 and LIBINT_ERI_DERIV>0 else True

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

H2ref = [-2.97603189186708, -0.0420345748265952, 0.898032292940286, 1.06445776378604, 1.08388172283809, -0.9676177714966, 1.91157412808118, -1.04184714801149, 0.0695854785563141, -1.35473830767443, -1.00051785464708, 1.01230086802415, 1.43977388170385, 0.77581343019648, -0.97026629319755, -0.0850355740292769, 0.224704424450602, -4.00998703183294, -0.946954158598755, 0.811603857603426, 1.6661469701995, 0.0489218656584689, 0.188913997043655, 2.34384006163361, -1.07351307133863, -1.05201061709615, 0.957147912326479, 0.00905530755259405, 0.039709749072005, -0.0101937537277215, -1.43185263083151, -0.799836089967502, -0.0318711057419372, -0.00792125087234215, -0.0117677676359248, -1.67096353593361, 0.0104698591701241, 0.0240226597710215, 0.00481656573411065, -1.92062943563377, 1.00213739893949, -0.0593917248285927, 0.0929568249016185, -0.212936656814677, -2.34865662736772]
H2tol = 1e-9
H2ok = False if LIBINT_ERI_DERIV>1 else True

HNref = [-0.09513423158289, 0.00935257129329314, 0.840399911385928, -0.357531998623518, 0.977636861928377, -0.90631889299092, 0.452666230206408, -0.98698943322167, 0.0659189816049918, 1.20846783402573, -0.858405093048817, 0.977636861928378, 0.124944488934069, 0.709863445831558, -0.986989433221671, -1.3334123229598, 0.148541647217259, -1.11333360244284, -0.90631889299092, 0.709863445831558, 0.23258750968945, 0.0659189816049918, 0.148541647217259, 0.880746092753391, 0.331757531648002, -0.991491488723426, 0.911651541330682, 0.0257744669755165, 0.0138546267950485, -0.00533264833976203, -0.0813945274926787, -0.737572699421655, 0.0138546267950485, -0.0435499614413899, 0.027709253590097, -0.250363004155323, -0.00533264833976203, 0.027709253590097, 0.0177754944658734, -0.478440697181925, 0.973134806426622, -0.0605863332652298, 1.37696228440119, -0.176250900807356, -0.898521587219264]
HNtol = 1e-10
HNok = False if LIBINT_ONEBODY_DERIV>1 and LIBINT_ERI_DERIV>1 else True

Href = [0.10278518545586, 0.022126031631778, -0.0110723782691211, -0.0363824853068044, -0.0101626228447786, 0.0127813217832919, -0.0664027001502164, -0.0119634087869938, -0.00170894351417136, -0.0120165941959616, 0.04690840364466, -0.0847585610957569, -0.0210098962611897, -0.0525024938577335, 0.062632529463985, 0.0330264904576947, 0.0055940902130715, 0.0930779928063925, 0.034315314298125, -0.0152045247322665, -0.0345309104656565, -0.0232429360290052, -0.0317038789123941, -0.0585470823405211, 0.0206391607592142, 0.0373257379155841, -0.0241432459680846, 0.0157433245475938, 0.0474328231801721, -0.0101720683300368, 0.0487850284761392, 0.0182560100458433, -0.0271631150708055, -0.0277751322149512, -0.00305148531357672, 0.0300339487209513, 0.0113619241847971, 0.0342464838118885, 0.00449696174470318, 0.050659375602623, -0.0354694143931759, 0.0118810118442084, -0.00525135824274559, -0.00254260489949476, 0.0540501205958149]
Htol = 1e-9
Hok = False if LIBINT_ONEBODY_DERIV>1 and LIBINT_ERI_DERIV>1 else True

for line in instr:
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
    match11 = re.match('\*\* 2-body hessian =' + pat_numbers(45), line)
    match12 = re.match('\*\* nuclear repulsion hessian =' + pat_numbers(45), line)
    match13 = re.match('\*\* Hartree-Fock hessian =' + pat_numbers(45), line)
    match14 = re.match('\*\* sph edipole =' + pat_numbers(3), line)
    match15 = re.match('\*\* sph equadrupole =' + pat_numbers(5), line)
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
        FNok = validate("nuclear repulsion force", match7.groups(), FNref, FNtol, line)
    elif match8:
        Fok = validate("HF force", match8.groups(), Fref, Ftol, line)
    elif match9:
        H1ok = validate("1-body hessian", match9.groups(), H1ref, H1tol, line)
    elif match10:
        HPok = validate("Pulay hessian", match10.groups(), HPref, HPtol, line)
    elif match11:
        H2ok = validate("2-body hessian", match11.groups(), H2ref, H2tol, line)
    elif match12:
        HNok = validate("nuclear repulsion hessian", match12.groups(), HNref, HNtol, line)
    elif match13:
        Hok = validate("HF hessian", match13.groups(), Href, Htol, line)
    elif match14:
        smuok = validate("spherical electric dipole moment", match14.groups(), smuref, smutol, line)
    elif match15:
        sQok = validate("spherical electric quadrupole moment", match15.groups(), sQref, sQtol, line)
    else:
        print(line,end="")

if len(sys.argv) > 2:
    instr.close()
ok = eok and muok and Qok and F1ok and FPok and F2ok and FNok and Fok and H1ok and HPok and H2ok and HNok and Hok and smuok and sQok
if not ok: sys.exit(1)


