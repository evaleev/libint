# <<<  User edit

# * For each integrals class, supply a list with max AM for each enabled derivative level.
# * User is responsible for "distributing" defaults across classes and derivs.
#   e.g., `--enable-eri2=2 --with-max-am=4` (`--with-eri2-max-am` not specified) becomes `eri2_max_am = [4, 4, 4]`
orderings = "ss"
# max am:
# * only used for paired on eri3.
# * must be same length as eri3_max_am.
# * if `--with-max-am` not given to ./configure, duplicate `eri3_max_am` as `max_am` here.
# * this script will NOT use this as default for other integrals classes.
max_am = [6, 5]
multipole = [10]
onebody = [6, 5, 4]
eri_max_am = [5, 4]
eri3_max_am = [6, 5]
eri3_pure_sh = False
eri2_max_am = [6, 5]
eri2_pure_sh = False
g12_max_am = [4, 4]

# >>>  End user edit

amstr = "SPDFGHIKLMNOQRTUVWXYZ"
components = [orderings]

# multipole
for deriv in range(len(multipole)):
    for am in range(multipole[deriv], 1, -1):  # don't enumerate s, p
        centers = amstr[am].lower() * 2
        comp = f"multipole_{centers}_d{deriv}"
        components.append(comp)

# onebody
for deriv in range(len(onebody)):
    for am in range(onebody[deriv], 1, -1):
        centers = amstr[am].lower() * 2
        comp = f"onebody_{centers}_d{deriv}"
        components.append(comp)

# eri (4-center)
for deriv in range(len(eri_max_am)):
    for am in range(eri_max_am[deriv], 1, -1):
        centers = amstr[am].lower() * 4
        comp = f"eri_{centers}_d{deriv}"
        components.append(comp)

# eri3
no_pure_sh = []
for deriv in range(len(eri3_max_am)):
    for am_fitting in range(eri3_max_am[deriv], 1, -1):
        for am_paired in range(max_am[deriv], 1, -1):
            if am_fitting >= am_paired:
                centers = amstr[am_paired].lower() * 2 + amstr[am_fitting].upper()
                comp = f"eri_{centers}_d{deriv}"
                #print(deriv, am_fitting, am_paired, centers, comp)
                components.append(comp)
                no_pure_sh.append(comp.lower())
if not eri3_pure_sh:
    components.extend(no_pure_sh)

# eri2
no_pure_sh = []
for deriv in range(len(eri2_max_am)):
    for am in range(eri2_max_am[deriv], 1, -1):
        centers = amstr[am].upper() * 2
        comp = f"eri_{centers}_d{deriv}"
        components.append(comp)
        no_pure_sh.append(comp.lower())
if not eri2_pure_sh:
    components.extend(no_pure_sh)

# g12
for deriv in range(len(g12_max_am)):
    for am in range(g12_max_am[deriv], 1, -1):
        centers = amstr[am].lower() * 4
        comp = f"g12_{centers}_d{deriv}"
        components.append(comp)

for comp in components:
    print(comp)
components = ";".join(components)
print(components)


# An example

#./configure \
#    --enable-eri=1 \
#    --enable-eri3=1 \
#    --enable-eri2=1 \
#    --enable-1body=2 \
#    --enable-g12=1 \
#    --disable-1body-property-derivs \
#    --with-multipole-max-order=10 \
#    --with-g12-max-am=4 \
#    --with-eri-max-am=5,4 \
#    --with-eri3-max-am=6,5 \
#    --with-eri2-max-am=6,5 \
#    --with-max-am=6,5

## script headmatter
#orderings = "ss"
#max_am = [6, 5]
#multipole = [10]
#onebody = [6, 5, 4]
#eri_max_am = [5, 4]
#eri3_max_am = [6, 5]
#eri3_pure_sh = False
#eri2_max_am = [6, 5]
#eri2_pure_sh = False
#g12_max_am = [4, 4]

#ans = "ss;multipole_nn_d0;multipole_mm_d0;multipole_ll_d0;multipole_kk_d0;multipole_ii_d0;multipole_hh_d0;multipole_gg_d0;multipole_ff_d0;multipole_dd_d0;onebody_ii_d0;onebody_hh_d0;onebody_gg_d0;onebody_ff_d0;onebody_dd_d0;onebody_hh_d1;onebody_gg_d1;onebody_ff_d1;onebody_dd_d1;onebody_gg_d2;onebody_ff_d2;onebody_dd_d2;eri_hhhh_d0;eri_gggg_d0;eri_ffff_d0;eri_dddd_d0;eri_gggg_d1;eri_ffff_d1;eri_dddd_d1;eri_iiI_d0;eri_hhI_d0;eri_hhH_d0;eri_ggI_d0;eri_ggH_d0;eri_ggG_d0;eri_ffI_d0;eri_ffH_d0;eri_ffG_d0;eri_ffF_d0;eri_ddI_d0;eri_ddH_d0;eri_ddG_d0;eri_ddF_d0;eri_ddD_d0;eri_hhH_d1;eri_ggH_d1;eri_ggG_d1;eri_ffH_d1;eri_ffG_d1;eri_ffF_d1;eri_ddH_d1;eri_ddG_d1;eri_ddF_d1;eri_ddD_d1;eri_iii_d0;eri_hhi_d0;eri_hhh_d0;eri_ggi_d0;eri_ggh_d0;eri_ggg_d0;eri_ffi_d0;eri_ffh_d0;eri_ffg_d0;eri_fff_d0;eri_ddi_d0;eri_ddh_d0;eri_ddg_d0;eri_ddf_d0;eri_ddd_d0;eri_hhh_d1;eri_ggh_d1;eri_ggg_d1;eri_ffh_d1;eri_ffg_d1;eri_fff_d1;eri_ddh_d1;eri_ddg_d1;eri_ddf_d1;eri_ddd_d1;eri_II_d0;eri_HH_d0;eri_GG_d0;eri_FF_d0;eri_DD_d0;eri_HH_d1;eri_GG_d1;eri_FF_d1;eri_DD_d1;eri_ii_d0;eri_hh_d0;eri_gg_d0;eri_ff_d0;eri_dd_d0;eri_hh_d1;eri_gg_d1;eri_ff_d1;eri_dd_d1;g12_gggg_d0;g12_ffff_d0;g12_dddd_d0;g12_gggg_d1;g12_ffff_d1;g12_dddd_d1"

# Another example

#./configure \
#    --with-max-am=2,2 \
#    --with-eri-max-am=2,2 \
#    --with-eri3-max-am=3,2 \
#    --enable-eri=1 \
#    --enable-eri3=1 \
#    --enable-1body=1 \
#    --disable-1body-property-derivs \
#    --with-multipole-max-order=2 \
#    --enable-eri3-pure-sh

## script headmatter
#orderings = "ss"
#max_am = [2, 2]
#multipole = [2]
#onebody = [2, 2]
#eri_max_am = [2, 2]
#eri3_max_am = [3, 2]
#eri3_pure_sh = True
#eri2_max_am = []
#eri2_pure_sh = False
#g12_max_am = []

#ans = "ss;multipole_dd_d0;onebody_dd_d0;onebody_dd_d1;eri_dddd_d0;eri_dddd_d1;eri_ddF_d0;eri_ddD_d0;eri_ddD_d1"

# Check examples

#ans = ans.split(";")
#for idx, comp in enumerate(components.split(";")):
#    print(comp, ans[idx], comp == ans[idx])
#ans = ";".join(ans)
#
#assert components == ans

