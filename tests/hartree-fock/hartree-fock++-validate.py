from __future__ import print_function
import math
import os
import re
import sys

def pat_numbers(n):
    result = ''
    for i in range(n):
        result += r'\s*([+-e\d.]+)'
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
else:
    sys.exit(2)

# ignore test if ERI_MAX_AM  is too low
if LIBINT_ERI_MAX_AM<2:
    sys.exit(0)

eref = [-76.003354058439]
etol = 5e-12

muref = [-0.263282355191187, -0.0912036831854911, -0.105312942076475]
mutol = 1e-9

Qref = [-7.01719638074453, 0.024624345507934, -0.576201299024159,
        -8.04716669024124, 0.638204907990955, -6.35134519242917 ]
Qtol = 1e-8

# * Note that the sphemultipole order is fixed at generator build time (not influenced by
#   libint2::solid_harmonics_ordering())
# * As of https://github.com/evaleev/libint/commit/cdbb9f3, sphemultipole is hard-coded at Standard
#   (and not influenced by LIBINT_SHGSHELL_ORDERING)
smuref = [-muref[1]/2, muref[2], -muref[0]/2]  # [ 0.046, -0.105,  0.132]
smutol = 1e-9

sQref = [Qref[1]/4, -Qref[4]/2, (2*Qref[5] - Qref[0] - Qref[3])/4,
         -Qref[2]/2, (Qref[0] - Qref[3])/8  ]
sQtol = 1e-8
print("Checking sphemultipole Standard\n");

# Reference values from LIBINT_SHGSHELL_ORDERING=Gaussian. These are no longer accessible after cdbb9f3.
#   print("Checking sphemultipole Gaussian\n");
#   smuref = [ muref[2],   -muref[0]/2, -muref[1]/2]  # [-0.105,  0.132,  0.046]
#   sQref = [(2*Qref[5] - Qref[0] - Qref[3])/4, -Qref[2]/2, -Qref[4]/2, (Qref[0] - Qref[3])/8, Qref[1]/4]

F1ref = [-5.43569555903312, -1.88298017654395, -2.17427822361352,
         3.47022732536532, -2.96798871167808, 2.59189820350226,
         1.9654682336678, 4.85096888822203, -0.417619979888738 ]
F1tol = 1e-9
F1ok = False if LIBINT_ONEBODY_DERIV>0 else True

FPref = [ 0.355265310323155, 0.123067513529209, 0.142106124129298,
          -0.224258642015539, 0.180741977786854, -0.164305035772324,
          -0.131006668307617, -0.303809491316068, 0.0221989116430271 ]
FPtol = 1e-9
FPok = False if LIBINT_ONEBODY_DERIV>0 else True

F2ref = [ 2.95100851670393, 1.02225933689958, 1.18040340668181,
          -1.89116113654409, 1.64868682617684, -1.42151545972338,
          -1.05984738015984, -2.67094616307642, 0.241112053041571 ]
F2tol = 1e-9
F2ok = False if LIBINT_ERI_DERIV>0 else True

FNref = [ 2.01500148332517 ,  0.698016989289171,  0.806000593330069,
         -1.28187870168764 ,  1.07670120707691 , -0.951756216715147,
         -0.733122781637531, -1.77471819636609 ,  0.145755623385078 ]
FNtol = 1e-10
FNok = False if LIBINT_ONEBODY_DERIV>0 and LIBINT_ERI_DERIV>0 else True

Fref = [ -0.114420248680859, -0.0396363368259882, -0.0457680994723476,
         0.0729288451180514, -0.0618587006374771, 0.0543214912914121,
         0.0414914035628126, 0.101495037463456, -0.00855339181906306 ]
Ftol = 1e-9
Fok = False if LIBINT_ONEBODY_DERIV>0 and LIBINT_ERI_DERIV>0 else True

H1ref = [ 3.54760842822125, 0.0561425790590885, -1.84088034418252, -0.887601426360438, -2.17549130368644, 1.98344832758524, -2.66000700186171, 2.11934872462735, -0.142567983402711, 0.355853985722548, 2.00292709856151, -2.17850638711514, -1.77309201081378, -1.61619662181537, 2.12236380805605, 1.41723802509125, -0.386730476746138, 5.69831976160228, 1.98431870719983, -1.61468908010102, -2.14019765060108, -0.143438363017304, -0.38823801846049, -3.55812211100127, 0.904263134266358, 2.18335388755013, -1.98924643233292, -0.0166617079059194, -0.00484750043498585, 0.00492772513309295, 1.75328800656853, 1.63118422308913, -0.00786258386369225, 0.0198040042452449, -0.0164951429881064, 2.15669050290334, 0.00579810474768649, -0.0149876012737541, -0.0164928523022598, 2.67666870976763, -2.11450122419236, 0.137640258269618, -1.43704202933649, 0.403225619734244, 3.57461496330354 ]
H1tol = 1e-9
H1ok = False if LIBINT_ONEBODY_DERIV>1 else True

HPref = [ -0.373657119664361, -0.00133454389521286, 0.0913757616853391, 0.144293176019879, 0.103810096183052, -0.0967303414175233, 0.229363943644482, -0.102475552287839, 0.00535457973218418, -0.221600106461453, -0.0970957473238979, 0.103810096183052, 0.187363744093183, 0.0780172520149254, -0.102475552287839, 0.0342363623684173, 0.0190784953089725, -0.481921134987926, -0.0967303414175233, 0.0780172520149254, 0.206932260444167, 0.00535457973218418, 0.0190784953089725, 0.27498887454376, -0.141868433944075, -0.102526043925858, 0.096303732812767, -0.00242474207580482, -0.00128405225719423, 0.000426608604756398, -0.191255819949221, -0.0755194237363812, -0.00128405225719423, 0.00389207585603787, -0.00249782827854418, -0.205330014292075, 0.000426608604756397, -0.00249782827854418, -0.00160224615209187, -0.226939201568677, 0.103759604545033, -0.00578118833694058, -0.038128438224455, -0.0165806670304283, -0.273386628391668 ]
HPtol = 1e-9
HPok = False if LIBINT_ONEBODY_DERIV>1 else True

H2ref = [ -2.97603189152785, -0.04203457482293, 0.898032292915088, 1.06445776362688, 1.08388172281076, -0.967617771471307, 1.91157412790108, -1.04184714798783, 0.0695854785562199, -1.3547383073865, -1.00051785461617, 1.01230086799738, 1.43977388154118, 0.775813430173996, -0.970266293174445, -0.0850355741545849, 0.224704424442172, -4.00998703146547, -0.946954158573641, 0.811603857580693, 1.66614697002347, 0.0489218656585546, 0.188913997035478, 2.34384006144209, -1.07351307118258, -1.0520106170714, 0.957147912300498, 0.00905530755570372, 0.0397097490740273, -0.0101937537268573, -1.43185263066289, -0.799836089947826, -0.0318711057393559, -0.00792125087828533, -0.0117677676328675, -1.67096353575741, 0.0104698591708088, 0.024022659773826, 0.004816565733943, -1.92062943545679, 1.0021373989138, -0.0593917248293631, 0.0929568250328693, -0.212936656809304, -2.34865662717603 ]
H2tol = 1e-9
H2ok = False if LIBINT_ERI_DERIV>1 else True

HNref = [ -0.0951342315737215, 0.00935257129239186, 0.840399911304933, -0.357531998589061, 0.977636861834156, -0.906318892903572, 0.452666230162782, -0.986989433126548, 0.0659189815986388, 1.20846783390926, -0.858405092966087, 0.977636861834156, 0.124944488922027, 0.709863445763144, -0.986989433126549, -1.33341232283129, 0.148541647202943, -1.11333360233554, -0.906318892903572, 0.709863445763144, 0.232587509667034, 0.0659189815986388, 0.148541647202943, 0.880746092668508, 0.331757531616028, -0.99149148862787, 0.91165154124282, 0.0257744669730324, 0.0138546267937132, -0.00533264833924809, -0.0813945274848341, -0.737572699350571, 0.0138546267937132, -0.0435499614371927, 0.0277092535874264, -0.250363004131194, -0.00533264833924809, 0.0277092535874264, 0.0177754944641603, -0.478440697135815, 0.973134806332835, -0.0605863332593907, 1.37696228426848, -0.176250900790369, -0.898521587132668 ]
HNtol = 1e-10
HNok = False if LIBINT_ONEBODY_DERIV>1 and LIBINT_ERI_DERIV>1 else True

Href = [ 0.102785185455314, 0.0221260316333375, -0.0110723782771648, -0.0363824853027409, -0.0101626228584681, 0.012781321792833, -0.0664027001533642, -0.0119634087748708, -0.00170894351566768, -0.0120165942161412, 0.0469084036553551, -0.0847585611005566, -0.0210098962573883, -0.0525024938633059, 0.0626325294672223, 0.0330264904737867, 0.0055940902079491, 0.0930779928133423, 0.0343153143050929, -0.0152045247422558, -0.0345309104664071, -0.0232429360279265, -0.0317038789130969, -0.0585470823469183, 0.0206391607557284, 0.0373257379249965, -0.0241432459768363, 0.015743324547012, 0.0474328231755604, -0.010172068328256, 0.0487850284715841, 0.0182560100543473, -0.0271631150665291, -0.0277751322141953, -0.00305148531209163, 0.0300339487226562, 0.0113619241840036, 0.0342464838089542, 0.00449696174375165, 0.0506593756063487, -0.0354694144006902, 0.0118810118439231, -0.005251358259593, -0.00254260489585742, 0.0540501206031722 ]
Htol = 1e-9
Hok = False if LIBINT_ONEBODY_DERIV>1 and LIBINT_ERI_DERIV>1 else True

for line in instr:
    match1 = re.match(r'\*\* Hartree-Fock energy =' + pat_numbers(1), line)
    match2 = re.match(r'\*\* edipole =' + pat_numbers(3), line)
    match3 = re.match(r'\*\* equadrupole =' + pat_numbers(6), line)
    match4 = re.match(r'\*\* 1-body forces =' + pat_numbers(9), line)
    match5 = re.match(r'\*\* Pulay forces =' + pat_numbers(9), line)
    match6 = re.match(r'\*\* 2-body forces =' + pat_numbers(9), line)
    match7 = re.match(r'\*\* nuclear repulsion forces =' + pat_numbers(9), line)
    match8 = re.match(r'\*\* Hartree-Fock forces =' + pat_numbers(9), line)
    match9 = re.match(r'\*\* 1-body hessian =' + pat_numbers(45), line)
    match10 = re.match(r'\*\* Pulay hessian =' + pat_numbers(45), line)
    match11 = re.match(r'\*\* 2-body hessian =' + pat_numbers(45), line)
    match12 = re.match(r'\*\* nuclear repulsion hessian =' + pat_numbers(45), line)
    match13 = re.match(r'\*\* Hartree-Fock hessian =' + pat_numbers(45), line)
    match14 = re.match(r'\*\* sph edipole =' + pat_numbers(3), line)
    match15 = re.match(r'\*\* sph equadrupole =' + pat_numbers(5), line)
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


