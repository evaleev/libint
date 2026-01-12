/*
 *  Copyright (C) 2004-2026 Edward F. Valeev
 *
 *  This file is part of Libint library.
 *
 *  Libint library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Libint library.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef _libint2_src_lib_libint_shell_h_
#define _libint2_src_lib_libint_shell_h_

#include <libint2/util/cxxstd.h>
#if LIBINT2_CPLUSPLUS_STD < 2011
#error "libint2/shell.h requires C++11 support"
#endif

#include <libint2.h>
#include <libint2/config.h>
#include <libint2/util/small_vector.h>

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <type_traits>
#include <vector>

namespace libint2 {

namespace math {

/// @brief fac[k] = k!, `0<=k<=20`
static constexpr std::array<int64_t, 21> fac = {{1LL,
                                                 1LL,
                                                 2LL,
                                                 6LL,
                                                 24LL,
                                                 120LL,
                                                 720LL,
                                                 5040LL,
                                                 40320LL,
                                                 362880LL,
                                                 3628800LL,
                                                 39916800LL,
                                                 479001600LL,
                                                 6227020800LL,
                                                 87178291200LL,
                                                 1307674368000LL,
                                                 20922789888000LL,
                                                 355687428096000LL,
                                                 6402373705728000LL,
                                                 121645100408832000LL,
                                                 2432902008176640000LL}};

/// @brief fac_ld[k] = static_cast<long double>(k!), `0<=k<=170`
/// @note This is an inexact representation for most values, so should only be
/// used when `fac` does not
///       suffice, i.e. for `k>=21`.
/// @note This should be sufficient for computing with `double` reals
static constexpr std::array<long double, 171> fac_ld = {
    {1.l,
     1.l,
     2.l,
     6.l,
     24.l,
     120.l,
     720.l,
     5040.l,
     40320.l,
     362880.l,
     3.6288e6l,
     3.99168e7l,
     4.790016e8l,
     6.2270208e9l,
     8.71782912e10l,
     1.307674368e12l,
     2.0922789888e13l,
     3.55687428096e14l,
     6.402373705728e15l,
     1.21645100408832e17l,
     2.43290200817664e18l,
     5.109094217170944e19l,
     1.12400072777760768e21l,
     2.585201673888497664e22l,
     6.2044840173323943936e23l,
     1.5511210043330985984e25l,
     4.03291461126605635584e26l,
     1.0888869450418352160768e28l,
     3.04888344611713860501504e29l,
     8.841761993739701954543616e30l,
     2.6525285981219105863630848e32l,
     8.22283865417792281772556288e33l,
     2.6313083693369353016721801216e35l,
     8.68331761881188649551819440128e36l,
     2.9523279903960414084761860964352e38l,
     1.03331479663861449296666513375232e40l,
     3.719933267899012174679994481508352e41l,
     1.37637530912263450463159795815809024e43l,
     5.23022617466601111760007224100074291e44l,
     2.03978820811974433586402817399028974e46l,
     8.15915283247897734345611269596115894e47l,
     3.34525266131638071081700620534407517e49l,
     1.40500611775287989854314260624451157e51l,
     6.04152630633738356373551320685139975e52l,
     2.65827157478844876804362581101461589e54l,
     1.19622220865480194561963161495657715e56l,
     5.50262215981208894985030542880025489e57l,
     2.5862324151116818064296435515361198e59l,
     1.2413915592536072670862289047373375e61l,
     6.08281864034267560872252163321295377e62l,
     3.04140932017133780436126081660647688e64l,
     1.55111875328738228022424301646930321e66l,
     8.0658175170943878571660636856403767e67l,
     4.27488328406002556429801375338939965e69l,
     2.30843697339241380472092742683027581e71l,
     1.2696403353658275925965100847566517e73l,
     7.1099858780486345185404564746372495e74l,
     4.05269195048772167556806019054323221e76l,
     2.35056133128287857182947491051507468e78l,
     1.38683118545689835737939019720389406e80l,
     8.32098711274139014427634118322336438e81l,
     5.07580213877224798800856812176625227e83l,
     3.14699732603879375256531223549507641e85l,
     1.98260831540444006411614670836189814e87l,
     1.26886932185884164103433389335161481e89l,
     8.24765059208247066672317030678549625e90l,
     5.44344939077443064003729240247842753e92l,
     3.64711109181886852882498590966054644e94l,
     2.48003554243683059960099041856917158e96l,
     1.71122452428141311372468338881272839e98l,
     1.19785716699698917960727837216890987e100l,
     8.5047858856786231752116764423992601e101l,
     6.12344583768860868615240703852746727e103l,
     4.47011546151268434089125713812505111e105l,
     3.30788544151938641225953028221253782e107l,
     2.48091408113953980919464771165940337e109l,
     1.88549470166605025498793226086114656e111l,
     1.45183092028285869634070784086308285e113l,
     1.13242811782062978314575211587320462e115l,
     8.94618213078297528685144171539831652e116l,
     7.15694570462638022948115337231865322e118l,
     5.79712602074736798587973423157810911e120l,
     4.75364333701284174842138206989404947e122l,
     3.94552396972065865118974711801206106e124l,
     3.31424013456535326699938757913013129e126l,
     2.81710411438055027694947944226061159e128l,
     2.42270953836727323817655232034412597e130l,
     2.1077572983795277172136005186993896e132l,
     1.85482642257398439114796845645546284e134l,
     1.65079551609084610812169192624536193e136l,
     1.48571596448176149730952273362082574e138l,
     1.35200152767840296255166568759495142e140l,
     1.24384140546413072554753243258735531e142l,
     1.15677250708164157475920516230624044e144l,
     1.08736615665674308027365285256786601e146l,
     1.03299784882390592625997020993947271e148l,
     9.91677934870949689209571401541893801e149l,
     9.61927596824821198533284259495636987e151l,
     9.42689044888324774562618574305724247e153l,
     9.33262154439441526816992388562667005e155l,
     9.33262154439441526816992388562667005e157l,
     9.42594775983835942085162312448293675e159l,
     9.61446671503512660926865558697259548e161l,
     9.90290071648618040754671525458177335e163l,
     1.02990167451456276238485838647650443e166l,
     1.08139675824029090050410130580032965e168l,
     1.14628056373470835453434738414834943e170l,
     1.22652020319613793935175170103873389e172l,
     1.3246418194518289744998918371218326e174l,
     1.44385958320249358220488210246279753e176l,
     1.58824554152274294042537031270907729e178l,
     1.76295255109024466387216104710707579e180l,
     1.97450685722107402353682037275992488e182l,
     2.23119274865981364659660702121871512e184l,
     2.54355973347218755712013200418933523e186l,
     2.92509369349301569068815180481773552e188l,
     3.3931086844518982011982560935885732e190l,
     3.96993716080872089540195962949863065e192l,
     4.68452584975429065657431236280838416e194l,
     5.57458576120760588132343171174197716e196l,
     6.68950291344912705758811805409037259e198l,
     8.09429852527344373968162284544935083e200l,
     9.87504420083360136241157987144820801e202l,
     1.21463043670253296757662432418812959e205l,
     1.50614174151114087979501416199328069e207l,
     1.88267717688892609974376770249160086e209l,
     2.37217324288004688567714730513941708e211l,
     3.01266001845765954480997707752705969e213l,
     3.85620482362580421735677065923463641e215l,
     4.97450422247728744039023415041268096e217l,
     6.46685548922047367250730439553648525e219l,
     8.47158069087882051098456875815279568e221l,
     1.11824865119600430744996307607616903e224l,
     1.48727070609068572890845089118130481e226l,
     1.99294274616151887673732419418294845e228l,
     2.6904727073180504835953876621469804e230l,
     3.65904288195254865768972722051989335e232l,
     5.01288874827499166103492629211225388e234l,
     6.91778647261948849222819828311491036e236l,
     9.6157231969410890041971956135297254e238l,
     1.34620124757175246058760738589416156e241l,
     1.89814375907617096942852641411076779e243l,
     2.69536413788816277658850750803729027e245l,
     3.85437071718007277052156573649332508e247l,
     5.55029383273930478955105466055038812e249l,
     8.04792605747199194484902925779806277e251l,
     1.17499720439091082394795827163851716e254l,
     1.72724589045463891120349865930862023e256l,
     2.55632391787286558858117801577675794e258l,
     3.80892263763056972698595524350736934e260l,
     5.713383956445854590478932865261054e262l,
     8.62720977423324043162318862654419154e264l,
     1.31133588568345254560672467123471711e267l,
     2.00634390509568239477828874698911719e269l,
     3.08976961384735088795856467036324047e271l,
     4.78914290146339387633577523906302272e273l,
     7.47106292628289444708380937293831545e275l,
     1.17295687942641442819215807155131553e278l,
     1.85327186949373479654360975305107853e280l,
     2.94670227249503832650433950735121486e282l,
     4.71472363599206132240694321176194378e284l,
     7.59070505394721872907517857093672949e286l,
     1.22969421873944943411017892849175018e289l,
     2.00440157654530257759959165344155279e291l,
     3.28721858553429622726333031164414657e293l,
     5.42391066613158877498449501421284184e295l,
     9.00369170577843736647426172359331746e297l,
     1.50361651486499904020120170784008402e300l,
     2.52607574497319838753801886917134115e302l,
     4.26906800900470527493925188889956654e304l,
     7.25741561530799896739672821112926311e306l}};

/// @brief computes factorial as an integer

/// @param k non-negative integer
/// @pre `k<fac.size()`
/// @return factorial of \p k
template <typename Int, typename = typename std::enable_if<
                            std::is_integral<Int>::value>::type>
constexpr int64_t fac_int(Int k) {
#if LIBINT2_CPLUSPLUS_STD >= 2014
  assert(!std::is_signed<Int>::value || k >= 0);
  assert(k < fac.size());
#endif
  return fac[k];
}

/// @brief computes factorial as a real number

/// @param k non-negative integer
/// @pre `k<fac_ld.size()`
/// @return factorial of \p k , uses exact representation when possible
template <
    typename Real, typename Int,
    typename = typename std::enable_if<std::is_integral<Int>::value>::type>
constexpr Real fac_real(Int k) {
#if LIBINT2_CPLUSPLUS_STD >= 2014
  assert(!std::is_signed<Int>::value || k >= 0);
  assert(k < fac.size() || k < fac_ld.size());
  // for extended-precision types should do better ...
  assert(k < fac.size() || (std::is_same<Real, double>::value ||
                            std::is_same<Real, long double>::value));
#endif

  // use exact representation when available
  return (k < fac.size()) ? static_cast<Real>(fac[k])
                          : static_cast<Real>(fac_ld[k]);
}

/// `df_Kminus1[k] = (k-1)!!`, `0<=k<=34`
static constexpr std::array<int64_t, 35> df_Kminus1 = {{1LL,
                                                        1LL,
                                                        1LL,
                                                        2LL,
                                                        3LL,
                                                        8LL,
                                                        15LL,
                                                        48LL,
                                                        105LL,
                                                        384LL,
                                                        945LL,
                                                        3840LL,
                                                        10395LL,
                                                        46080LL,
                                                        135135LL,
                                                        645120LL,
                                                        2027025LL,
                                                        10321920LL,
                                                        34459425LL,
                                                        185794560LL,
                                                        654729075LL,
                                                        3715891200LL,
                                                        13749310575LL,
                                                        81749606400LL,
                                                        316234143225LL,
                                                        1961990553600LL,
                                                        7905853580625LL,
                                                        51011754393600LL,
                                                        213458046676875LL,
                                                        1428329123020800LL,
                                                        6190283353629375LL,
                                                        42849873690624000LL,
                                                        191898783962510625LL,
                                                        1371195958099968000LL,
                                                        6332659870762850625LL}};

/// @brief fac_ld[k] = static_cast<long double>(k!), `0<=k<=170`
/// @note This is an inexact representation for most values, so should only be
/// used when `fac` does not
///       suffice, i.e. for `k>=21`.
/// @note This should be sufficient for computing with `double` reals
static constexpr std::array<long double, 302> df_Kminus1_ld = {
    {1.l,
     1.l,
     1.l,
     2.l,
     3.l,
     8.l,
     15.l,
     48.l,
     105.l,
     384.l,
     945.l,
     3840.l,
     10395.l,
     46080.l,
     135135.l,
     645120.l,
     2.027025e6l,
     1.032192e7l,
     3.4459425e7l,
     1.8579456e8l,
     6.54729075e8l,
     3.7158912e9l,
     1.3749310575e10l,
     8.17496064e10l,
     3.16234143225e11l,
     1.9619905536e12l,
     7.905853580625e12l,
     5.10117543936e13l,
     2.13458046676875e14l,
     1.4283291230208e15l,
     6.190283353629375e15l,
     4.2849873690624e16l,
     1.91898783962510625e17l,
     1.371195958099968e18l,
     6.332659870762850625e18l,
     4.6620662575398912e19l,
     2.21643095476699771875e20l,
     1.678343852714360832e21l,
     8.200794532637891559375e21l,
     6.3777066403145711616e22l,
     3.19830986772877770815625e23l,
     2.55108265612582846464e24l,
     1.3113070457687988603440625e25l,
     1.0714547155728479551488e26l,
     5.63862029680583509947946875e26l,
     4.71440074852053100265472e27l,
     2.5373791335626257947657609375e28l,
     2.1686243443194442612211712e29l,
     1.192568192774434123539907640625e30l,
     1.040939685273333245386162176e31l,
     5.8435841445947272053455474390625e31l,
     5.20469842636666622693081088e32l,
     2.980227913743310874726229193921875e33l,
     2.7064431817106664380040216576e34l,
     1.57952079428395476360490147277859375e35l,
     1.461479318123759876522171695104e36l,
     8.68736436856175119982695810028226562e36l,
     8.1842841814930553085241614925824e37l,
     4.95179769008019818390136611716089141e38l,
     4.746884825265972078944013665697792e39l,
     2.92156063714731692850180600912492593e40l,
     2.8481308951595832473664081994186752e41l,
     1.78215198865986332638610166556620482e42l,
     1.76584115499894161336717308363957862e43l,
     1.12275575285571389562324404930670903e44l,
     1.13013833919932263255499077352933032e45l,
     7.29791239356214032155108632049360873e45l,
     7.45891303871552937486293910529358011e46l,
     4.88960130368663401543922783473071785e47l,
     5.07206086632655997490679859159963447e48l,
     3.37382489954377747065306720596419531e49l,
     3.55044260642859198243475901411974413e50l,
     2.39541567867608200416367771623457867e51l,
     2.55631867662858622735302649016621577e52l,
     1.74865344543353986303948473285124243e53l,
     1.89167582070515380824123960272299967e54l,
     1.31149008407515489727961354963843182e55l,
     1.43767362373591689426334209806947975e56l,
     1.0098473647378692709053024332215925e57l,
     1.12138542651401517752540683649419421e58l,
     7.97779418142916724015188922245058078e58l,
     8.97108341211212142020325469195355365e59l,
     6.46201328695762546452303027018497043e60l,
     7.35628839793193956456666884740191399e61l,
     5.36347102817482913555411512425352546e62l,
     6.17928225426282923423600183181760775e63l,
     4.55895037394860476522099785561549664e64l,
     5.31418273866603314144296157536314267e65l,
     3.96628682533528614574226813438548208e66l,
     4.67648081002610916446980618631956555e67l,
     3.52999527454840466971061863960307905e68l,
     4.20883272902349824802282556768760899e69l,
     3.21229569983904824943666296203880193e70l,
     3.87212611070161838818099952227260027e71l,
     2.9874350008503148719760965546960858e72l,
     3.63979854405952128489013955093624426e73l,
     2.83806325080779912837729172696128151e74l,
     3.49420660229714043349453396889879449e75l,
     2.75292135328356515452597297515244306e76l,
     3.4243224702511976248246432895208186e77l,
     2.72539213975072950298071324540091863e78l,
     3.4243224702511976248246432895208186e79l,
     2.75264606114823679801052037785492782e80l,
     3.49280891965622157732113615531123497e81l,
     2.83522544298268390195083598919057565e82l,
     3.63252127644247044041398160152368437e83l,
     2.97698671513181809704837778865010444e84l,
     3.85047255302901866683882049761510543e85l,
     3.18537578519104536384176423385561175e86l,
     4.15851035727134016018592613742431386e87l,
     3.4720596058582394465875230149026168e88l,
     4.57436139299847417620451875116674525e89l,
     3.85398616250264578571215054654190465e90l,
     5.12328476015829107734906100130675468e91l,
     4.35500436362798973785473011759235226e92l,
     5.84054462658045182817792954148970034e93l,
     5.0082550181721881985329396352312051e94l,
     6.77503176683332412068639826812805239e95l,
     5.85965837126146019228353937322050996e96l,
     7.99453748486332246240994995639110182e97l,
     6.97299346180113762881741185413240686e98l,
     9.59344498183598695489193994766932219e99l,
     8.4373220887793765308690683435002123e100l,
     1.17040028778399040849681667361565731e102l,
     1.03779061691986331329689540625052611e103l,
     1.45129635685214810653605267528341506e104l,
     1.29723827114982914162111925781315764e105l,
     1.82863340963370661423542637085710298e106l,
     1.6474926043602830098588214574227102e107l,
     2.34065076433114446622134575469709181e108l,
     2.12526545962476508271787968007529616e109l,
     3.04284599363048780608774948110621935e110l,
     2.78409775210844225836042238089863797e111l,
     4.01655671159224390403582931506020954e112l,
     3.7028500103042282036193617665951885e113l,
     5.38218599353360683140801128218068079e114l,
     4.99884751391070807488613838490350448e115l,
     7.31977295120570529071489534376572587e116l,
     6.84842109405767006259400958731780114e117l,
     1.01012866726638733011865555743967017e119l,
     9.51930532074016138700567332637174358e119l,
     1.41418013417294226216611778041553824e121l,
     1.34222205022436275556779993901841585e122l,
     2.0081357905255780122758872481900643e123l,
     1.91937753182083874046195391279633466e124l,
     2.89171553835683233767727763739369259e125l,
     2.78309742114021617366983317355468525e126l,
     4.22190468600097521300882535059479118e127l,
     4.09115320907611777529465476512538732e128l,
     6.24841893528144331525306151888029095e129l,
     6.09581828152341548518903560003682711e130l,
     9.37262840292216497287959227832043642e131l,
     9.20468560510035738263544375605560894e132l,
     1.42463951724416907587769802630470634e134l,
     1.40831689758035467954322289467650817e135l,
     2.19394485655602037685165496050924776e136l,
     2.18289119124954975329199548674858766e137l,
     3.4225539762273917878885817383944265e138l,
     3.42713917026179311266843291419528263e139l,
     5.40763528243927902486395914666319387e140l,
     5.44915128071625104914280833357049938e141l,
     8.6522164519028464397823346346611102e142l,
     8.773133561953164189119921417048504e143l,
     1.40165906520826112324473821081509985e145l,
     1.43002077059836576282654719097890615e146l,
     2.29872086694154824212137066573676376e147l,
     2.35953427148730350866380286511519515e148l,
     3.81587663912297008192147530512302784e149l,
     3.9404222333837968594685507847423759e150l,
     6.41067275372658973762807851260668677e151l,
     6.65931357441861669250185082621461527e152l,
     1.08981436813352025539677334714313675e154l,
     1.13874262122558345441781649128269921e155l,
     1.87448071318965483928245015708619521e156l,
     1.97002473472025937614282252991906964e157l,
     3.26159644094999942035146327332997967e158l,
     3.44754328576045390824993942735837186e159l,
     5.74040973607199897981857536106076421e160l,
     6.1021516157960034176023927864243182e161l,
     1.02179293302081581840770641426881603e163l,
     1.09228513922748461175082830876995296e164l,
     1.83922727943746847313387154568386885e165l,
     1.97703610200174714726899923887361485e166l,
     3.34739364857619262110364621314464131e167l,
     3.61797606666319727950226860713871518e168l,
     6.15920431338019442283070903218614002e169l,
     6.69325572332691496707919692320662308e170l,
     1.14561200228871616264651187998662204e172l,
     1.25163882026213309884380982463963852e173l,
     2.15375056430278638577544233437484944e174l,
     2.3655973702954315568148005685689168e175l,
     4.09212607217529413297334043531221394e176l,
     4.51829097726427427351626908596663108e177l,
     7.85688205857656473530881363579945076e178l,
     8.72030158612004934788639933591559799e179l,
     1.52423511936385355864990984534509345e181l,
     1.70045880929340962283784787050354161e182l,
     2.98750083395315297495382329687638316e183l,
     3.34990385430801695699056030489197697e184l,
     5.91525165122724289040857012781523865e185l,
     6.66630867007295374441121500673503416e186l,
     1.18305033024544857808171402556304773e188l,
     1.33992804268466370262665421635374187e189l,
     2.38976166709580612772506233163735642e190l,
     2.72005392664986731633210805919809599e191l,
     4.87511380087544450055912715654020709e192l,
     5.57611054963222799848082152135609678e193l,
     1.00427344298034156711518019424728266e195l,
     1.15425488377387119568553005492071203e196l,
     2.08888876139911045959957480403434793e197l,
     2.41239270708739079898275781478428815e198l,
     4.38666639893813196515910708847213066e199l,
     5.090148611954394585853618989194848e200l,
     9.299732765748839766137307027560917e201l,
     1.08420165434628604678682084469850262e203l,
     1.99014281187025170995338370389803624e204l,
     2.33103355684451500059166481610178064e205l,
     4.29870847363974369349930880041975827e206l,
     5.05834281835259755128391265094086399e207l,
     9.37118447253464125182849318491507304e208l,
     1.10777707721921886373117687055604921e210l,
     2.06166058395762107540226850068131607e211l,
     2.44818734065447368884590088392886876e212l,
     4.57688649638591878739303607151252167e213l,
     5.45945776965947632612635897116137734e214l,
     1.02522257519044580837604008001880485e216l,
     1.2283779981733821733784307685113099e217l,
     2.31700301993040752692985058084249897e218l,
     2.78841805585357753356903784452067348e219l,
     5.28276688544132916140005932432089765e220l,
     6.38547734790469255187309666395234226e221l,
     1.21503638365150570712201364459380646e223l,
     1.47504526736598397948268532937299106e224l,
     2.81888441007149324052307165545763099e225l,
     3.43685547296274267219465681743906917e226l,
     6.59618951956729418282398767377085651e227l,
     8.07661036146244527965744352098181256e228l,
     1.55670072661788142714646109100992214e230l,
     1.91415665566659953127881411447268958e231l,
     3.70494772935055779660857739660361469e232l,
     4.57483440704317287975636573358972809e233l,
     8.89187455044133871186058575184867524e234l,
     1.10253509209740466402128414179512447e236l,
     2.15183364120680396827026175194737941e237l,
     2.67916027379669333357172046456215246e238l,
     5.25047408454460168257943867475160576e239l,
     6.56394267080189866725071513817727353e240l,
     1.29161662479797201391454191398889502e242l,
     1.62129383968806897081092663912978656e243l,
     3.20320922949897059450806394669245964e244l,
     4.03702166082329173731920733143316854e245l,
     8.0080230737474264862701598667311491e246l,
     1.0132924368666462260671210401897253e248l,
     2.01802181458435147454008028641624957e249l,
     2.56362986527261495194981623168000502e250l,
     5.12577540904425274533180392749727392e251l,
     6.53725615644516812747203139078401279e252l,
     1.31219850471532870280494180543930212e254l,
     1.68007483220640820876031206743149129e255l,
     3.38547214216554805323674985803339948e256l,
     4.35139381541459726068920825464756243e257l,
     8.80222756963042493841554963088683864e258l,
     1.1357137858232098850398833544630138e260l,
     2.30618362324317133386487400329235172e261l,
     2.98692725671504199765489322223772628e262l,
     6.08832476536197232140326736869180855e263l,
     7.91535723029486129378546703892997465e264l,
     1.61949438758628463749326912007202107e266l,
     2.11340038048872796544071969939430323e267l,
     4.34024495873124282848196124179301648e268l,
     5.68504702351467822703553599137067569e269l,
     1.17186613885743556369012953528411445e271l,
     1.54064774337247779952663025366145311e272l,
     3.1874758976922247332371523359727913e273l,
     4.205968339406864392707700592495767e274l,
     8.73368395967669576906979740056544817e275l,
     1.15664129333688770799461766293633592e277l,
     2.41049677287076803226326408255606369e278l,
     3.20389638254317895114509092633365051e279l,
     6.70118102858073512969187414950585707e280l,
     8.93887090729546927369480368447088492e281l,
     1.87633068800260583631372476186163998e283l,
     2.51182272495002686590823983533631866e284l,
     5.29125254016734845840470382844982474e285l,
     7.10845831160857603052031873400178182e286l,
     1.50271572140752696218693588727975023e288l,
     2.02591061880844416869829083919050782e289l,
     4.29776696322552711185463663762008565e290l,
     5.81436347598023476416409470847675744e291l,
     1.23775688540895180821413535163458467e293l,
     1.6803510445582878468434233707497829e294l,
     3.58949496768596024382099251974029553e295l,
     4.88982153966461763431436200888186824e296l,
     1.0481325305643003911957298157641663e298l,
     1.43271771112173296685410806860238739e299l,
     3.08150963985904315011544565834664891e300l,
     4.22651724780911225221961880237704281e301l,
     9.12126853398276772434171914870608078e302l,
     1.25527562259930633890922678430598171e304l,
     2.71813802312686478185383230631441207e305l,
     3.75327411157192595333858808507488533e306l,
     8.15441406938059434556149691894323621e307l}};

/// @brief computes double factorial as an integer

/// @param k non-negative integer
/// @pre `k<df_Kminus1.size()`
/// @return double factorial of `k-1`; N.B.: `-1!! = 1`
template <typename Int, typename = typename std::enable_if<
                            std::is_integral<Int>::value>::type>
constexpr int64_t df_Kminus1_int(Int k) {
#if LIBINT2_CPLUSPLUS_STD >= 2014
  assert(!std::is_signed<Int>::value || k >= 0);
  assert(k < df_Kminus1.size());
#endif
  return df_Kminus1[k];
}

/// @brief computes double factorial as a real number

/// @param k non-negative integer
/// @pre `k<df_Kminus1.size()`
/// @return double factorial of `k-1` (note that `-1!!=1`), uses exact
/// representation when possible
template <
    typename Real, typename Int,
    typename = typename std::enable_if<std::is_integral<Int>::value>::type>
constexpr Real df_Kminus1_real(Int k) {
#if LIBINT2_CPLUSPLUS_STD >= 2014
  assert(!std::is_signed<Int>::value || k >= 0);
  assert(k < df_Kminus1.size() || k < df_Kminus1_ld.size());
  // for extended-precision types should do better ...
  assert(k < df_Kminus1.size() || (std::is_same<Real, double>::value ||
                                   std::is_same<Real, long double>::value));
#endif
  // use exact representation when available
  return (k < df_Kminus1.size()) ? static_cast<Real>(df_Kminus1[k])
                                 : static_cast<Real>(df_Kminus1_ld[k]);
}

/// @brief computes binomial coefficient as an integer

/// @param i non-negative integer
/// @param j non-negative integer
/// @pre `i>=j`
/// return `i! / (j! * (i-j)!)`
template <typename Int, typename = typename std::enable_if<
                            std::is_integral<Int>::value>::type>
constexpr int64_t bc_int(Int i, Int j) {
#if LIBINT2_CPLUSPLUS_STD >= 2014
  assert(i >= j);
#endif
  return fac_int(i) / (fac_int(j) * fac_int(i - j));
}

/// @brief computes binomial coefficient as a real number

/// @param i non-negative integer
/// @param j non-negative integer
/// @pre `i>=j`
/// @return binomial coefficient, (\p i , \p j), uses exact representation when
/// possible
template <
    typename Real, typename Int,
    typename = typename std::enable_if<std::is_integral<Int>::value>::type>
constexpr Real bc_real(Int i, Int j) {
#if LIBINT2_CPLUSPLUS_STD >= 2014
  assert(!std::is_signed<Int>::value || (i >= 0 && j >= 0));
  assert(i >= j);
#endif
  return fac_real<Real>(i) / (fac_real<Real>(j) * fac_real<Real>(i - j));
}

}  // namespace math

/// generally-contracted Solid-Harmonic/Cartesion Gaussian Shell

/** A simple-to-use Gaussian shell. Here's an example of how to create an s+p
 * shell of the STO-3G basis on the oxygen atom located at the origin. \verbatim
 *  auto s = Shell{
 *                  {5.033151300, 1.169596100, 0.380389000},
 *                  {
 *                    {0, false, {-0.09996723, 0.39951283, 0.70011547}},
 *                    {1, false, {0.15591627, 0.60768372, 0.39195739}}
 *                  },
 *                   {{0.0, 0.0, 0.0}}
 *                };
 *  \endverbatim
 *  \note The contraction coefficients correspond to <em>unity-normalized</em>
 * primitives. EMSL Gaussian Basis Set Database, as well as basis set libraries
 * embedded into most quantum chemistry programs use this convention. However,
 * coefficients are automatically converted internally to refer to
 *  normalization-free primitives before computing integrals (see \c
 * Shell::renorm() ).
 */
struct Shell {
  typedef double real_t;

  /// contracted Gaussian = angular momentum + sph/cart flag + contraction
  /// coefficients
  struct Contraction {
    int l;
    bool pure;
    svector<real_t> coeff;

    bool operator==(const Contraction& other) const {
      return &other == this ||
             (l == other.l && pure == other.pure && coeff == other.coeff);
    }
    bool operator!=(const Contraction& other) const {
      return not this->operator==(other);
    }
    size_t cartesian_size() const { return (l + 1) * (l + 2) / 2; }
    size_t size() const { return pure ? (2 * l + 1) : cartesian_size(); }
  };

  svector<real_t> alpha;         //!< exponents
  svector<Contraction> contr;    //!< contractions
  std::array<real_t, 3> O;       //!< origin
  svector<real_t> max_ln_coeff;  //!< maximum ln of (absolute) contraction
                                 //!< coefficient for each primitive

  Shell() = default;
  Shell(const Shell&) = default;
  // intel does not support "move ctor = default"
  Shell(Shell&& other) noexcept
      : alpha(std::move(other.alpha)),
        contr(std::move(other.contr)),
        O(std::move(other.O)),
        max_ln_coeff(std::move(other.max_ln_coeff)) {}
  Shell& operator=(const Shell&) = default;
  // intel does not support "move asgnmt = default"
  Shell& operator=(Shell&& other) noexcept {
    alpha = std::move(other.alpha);
    contr = std::move(other.contr);
    O = std::move(other.O);
    max_ln_coeff = std::move(other.max_ln_coeff);
    return *this;
  }
  /// @param embed_normalization_into_coefficients if true, will embed
  /// normalization factors into coefficients, else will use the coefficients in
  /// @p _contr as given
  Shell(svector<real_t> _alpha, svector<Contraction> _contr,
        std::array<real_t, 3> _O,
        bool embed_normalization_into_coefficients = true)
      : alpha(std::move(_alpha)), contr(std::move(_contr)), O(std::move(_O)) {
    // embed normalization factors into contraction coefficients
    if (embed_normalization_into_coefficients)
      renorm();
    else {
      update_max_ln_coeff();
    }
  }

  Shell& move(std::array<real_t, 3> new_origin) {
    O = std::move(new_origin);
    return *this;
  }

  size_t cartesian_size() const {
    size_t s = 0;
    for (const auto& c : contr) {
      s += c.cartesian_size();
    }
    return s;
  }
  size_t size() const {
    size_t s = 0;
    for (const auto& c : contr) {
      s += c.size();
    }
    return s;
  }

  size_t ncontr() const { return contr.size(); }
  size_t nprim() const { return alpha.size(); }

  bool operator==(const Shell& other) const {
    return &other == this ||
           (O == other.O && alpha == other.alpha && contr == other.contr);
  }
  bool operator!=(const Shell& other) const {
    return not this->operator==(other);
  }

  /// @param l angular momentum quantum number
  /// @return (lower-case) letter symbol corresponding to @p l ; e.g., `s` for
  /// `l=0`, `p` for `l=1`, etc.
  /// @throw std::invalid_argument if \c l is greater than 20
  static char am_symbol(size_t l) {
    assert(l <= sizeof(LIBINT_AM2SYMBOL) - 1);
    static char lsymb[] = LIBINT_AM2SYMBOL;
    return lsymb[l];
  }

  /// inverse of am_symbol()
  /// @param am_symbol letter symbol denoting orbital angular momentum @p l ;
  /// e.g., `s` for `l=0`, `p` for `l=1`, etc.
  /// @note this function is case insensitive, i.e. `am_symbol_to_l('s') ==
  /// am_symbol_to_l('S')`
  /// @return angular momentum quantum number
  /// @sa am_symbol()
  static unsigned short am_symbol_to_l(char am_symbol) {
    const char AM_SYMBOL = ::toupper(am_symbol);
    switch (AM_SYMBOL) {
      case 'S':
        return 0;
      case 'P':
        return 1;
      case 'D':
        return 2;
      case 'F':
        return 3;
      case 'G':
        return 4;
      case 'H':
        return 5;
      case 'I':
        return 6;
      case 'K':
        return 7;
      case 'L':
        return 8;
      case 'M':
        return 9;
      case 'N':
        return 10;
      case 'O':
        return 11;
      case 'Q':
        return 12;
      case 'R':
        return 13;
      case 'T':
        return 14;
      case 'U':
        return 15;
      case 'V':
        return 16;
      case 'W':
        return 17;
      case 'X':
        return 18;
      case 'Y':
        return 19;
      case 'Z':
        return 20;
      default:
        throw std::invalid_argument{"invalid angular momentum label"};
    }
  }

  struct defaultable_boolean {
    typedef enum { false_value = 0, true_value = 1, default_value = 2 } value_t;
    defaultable_boolean() : value_(default_value) {}
    defaultable_boolean(bool v) : value_(static_cast<value_t>(v ? 1 : 0)) {}
    bool is_default() const { return value_ == default_value; }
    operator bool() const {
      assert(value_ != default_value);
      return value_ == true_value;
    }

   private:
    value_t value_;
  };

  /// sets and/or reports whether the auto-renormalization to unity is set
  /// if called without arguments, returns the current value of the flag
  /// otherwise, will set the flag to \c flag
  /// \note by default, shells WILL be re-normalized to unity
  static bool do_enforce_unit_normalization(
      defaultable_boolean flag = defaultable_boolean()) {
    static bool result{true};
    if (not flag.is_default()) {
      result = flag;
    }
    return result;
  }

  /// @return "unit" Shell, with exponent=0. and coefficient=1., located at the
  /// origin
  static const Shell& unit() {
    static const Shell unitshell{make_unit()};
    return unitshell;
  }

  /// @return the coefficient of primitive \c p in contraction \c c assuming
  /// unit normalized primitive
  ///         (coeff contains coefficients of normalization-free primitives; @sa
  ///         Shell::renorm() )
  real_t coeff_normalized(size_t c, size_t p) const {
    const auto alpha = this->alpha.at(p);
    assert(alpha >= 0.0);
    const auto l = contr.at(c).l;

    using libint2::math::df_Kminus1_real;
    const auto sqrt_Pi_cubed = real_t{5.56832799683170784528481798212};
    const auto two_alpha = 2 * alpha;
    const auto two_alpha_to_am32 = pow(two_alpha, l + 1) * sqrt(two_alpha);
    const auto one_over_N =
        sqrt((sqrt_Pi_cubed * df_Kminus1_real<real_t>(2 * l)) /
             (pow(2, l) * two_alpha_to_am32));
    return contr.at(c).coeff[p] * one_over_N;
  }

  /// extract primitive shell

  /// @param p the index of the primitive to extract
  /// @param unit_normalized whether to produce unit-normalized primitive; set
  /// to false to produce normalization-free primitive
  /// @return a primitive Shell
  Shell extract_primitive(size_t p, bool unit_normalized = true) const {
    assert(p < nprim());
    svector<Contraction> prim_contr;
    prim_contr.reserve(ncontr());
    for (auto&& c : contr) {
      prim_contr.emplace_back(Contraction{c.l, c.pure, {1.}});
    }
    return Shell({alpha[p]}, prim_contr, O, unit_normalized);
  }

 private:
  // this makes a unit shell
  struct make_unit {};
  Shell(make_unit)
      : alpha{0.0},              // exponent = 0
        contr{{0, false, {1}}},  // contraction coefficient = 1
        O{{0, 0, 0}},            // placed at origin
        max_ln_coeff{0} {}

  /// embeds normalization constants into contraction coefficients. Do this
  /// before computing integrals. \warning Must be done only once. \note this is
  /// now private
  void renorm() {
    using libint2::math::df_Kminus1_real;
    using std::pow;
    const auto sqrt_Pi_cubed = real_t{5.56832799683170784528481798212};
    const auto np = nprim();
    for (auto& c : contr) {
      for (auto p = 0ul; p != np; ++p) {
        assert(alpha[p] >= 0);
        if (alpha[p] != 0) {
          const auto two_alpha = 2 * alpha[p];
          const auto two_alpha_to_am32 =
              pow(two_alpha, c.l + 1) * sqrt(two_alpha);
          const auto normalization_factor =
              sqrt(pow(2, c.l) * two_alpha_to_am32 /
                   (sqrt_Pi_cubed * df_Kminus1_real<real_t>(2 * c.l)));

          c.coeff[p] *= normalization_factor;
        }
      }

      // need to force normalization to unity?
      if (do_enforce_unit_normalization()) {
        // compute the self-overlap of the , scale coefficients by its inverse
        // square root
        double norm{0};
        for (auto p = 0ul; p != np; ++p) {
          for (decltype(p) q = 0ul; q <= p; ++q) {
            auto gamma = alpha[p] + alpha[q];
            norm += (p == q ? 1 : 2) * df_Kminus1_real<real_t>(2 * c.l) *
                    sqrt_Pi_cubed * c.coeff[p] * c.coeff[q] /
                    (pow(2, c.l) * pow(gamma, c.l + 1) * sqrt(gamma));
          }
        }
        auto normalization_factor = 1 / sqrt(norm);
        for (auto p = 0ul; p != np; ++p) {
          c.coeff[p] *= normalization_factor;
        }
      }
    }

    update_max_ln_coeff();
  }

  void update_max_ln_coeff() {
    // update max log coefficients
    max_ln_coeff.resize(nprim());
    for (auto p = 0ul; p != nprim(); ++p) {
      real_t max_ln_c = -std::numeric_limits<real_t>::max();
      for (auto& c : contr) {
        max_ln_c = std::max(max_ln_c, std::log(std::abs(c.coeff[p])));
      }
      max_ln_coeff[p] = max_ln_c;
    }
  }
};

inline std::ostream& operator<<(std::ostream& os, const Shell& sh) {
  os << "Shell:( O={" << sh.O[0] << "," << sh.O[1] << "," << sh.O[2] << "}"
     << std::endl;
  os << "  ";
  for (const auto& c : sh.contr) {
    os << " {l=" << c.l << ",sph=" << c.pure << "}";
  }
  os << std::endl;

  for (auto i = 0ul; i < sh.alpha.size(); ++i) {
    os << "  " << sh.alpha[i];
    for (const auto& c : sh.contr) {
      os << " " << c.coeff.at(i);
    }
    os << std::endl;
  }

  return os;
}

// clang-format off
  /// @brief describes method for primitive screening used by ShellPair and Engine
  ///
  /// @note *Rationale*. Since Engine::compute2() can compute primitive data on-the-fly it needs cheap methods for screening primitives.
  ///       - ScreeningMethod::Original was the fast approach that works well for uncontracted integrals over spherical (L=0) Gaussians, but can introduce significant errors in integrals over nonspherical and contracted Gaussians.
  ///       - ScreeningMethod::Conservative addresses the weaknesses of the original method and works well across the board.
  ///       - ScreeningMethod::Schwarz and ScreeningMethod::SchwarzInf should be used when it is possible to precompute the shell pair data _and_ the data can be computed for a specific operator type. These approaches provide rigorous guarantees of precision.
  enum class ScreeningMethod {
    /// standard screening method:
    /// - omit primitive pair if \f$ {\rm ln_scr} \equiv \max|c_a| \max|c_b| \exp(- \alpha_a \alpha_b |AB|^2 / (\alpha_a + \alpha_b) ) < \epsilon \f$
    /// - omit primitive 2-body integral if `spbra.primpairs[pb].ln_scr + spket.primpairs[pk].ln_scr < log(ε)` or `scale * c_a * c_b * c_c * c_d * spbrapp.K * spketpp.K / sqrt(spbrapp.gamma + spketpp.gamma) < ε`
    Original = 0x0001,
    /// conservative screening:
    /// - omit primitive pair if \f$ {\rm ln_scr} \equiv (k_a k_b) \max|c_a| \max|c_b| \exp(- \rho |AB|^2 ) \sqrt{2 \pi^{5/2}} \max( (\max_i|PA_i|)^{L_a} (\max_i|PB_i|)^{L_b}, L_a! L_b! / (\alpha_a + \alpha_b)^{L_a+L_b} ) / (\alpha_a + \alpha_b) < \epsilon \f$
    /// - omit primitive 2-body integral if `spbra.primpairs[pb].ln_scr + spket.primpairs[pk].ln_scr < log(ε)` or `scale * c_a * c_b * c_c * c_d * max(1, spbrapp.nonsph_screen_fac * spketpp.nonsph_screen_fac) * spbrapp.K * spketpp.K / sqrt(spbrapp.gamma + spketpp.gamma) < ε / (npbra * npket)`
    Conservative = 0x0010,
    /// Schwarz screening method using Frobenius norm:
    /// - omit primitive pair if \f$ {\rm ln_scr} \equiv (k_a k_b) \max|c_a| \max|c_b| K_{ab}  < \epsilon \f$, where \f$ K_{ab} \equiv \sqrt{||(ab|ab)||_2} \f$
    /// - omit primitive 2-body integral if `spbra.primpairs[pb].ln_scr + spket.primpairs[pk].ln_scr < log(ε)`
    Schwarz = 0x0100,
    /// Schwarz screening method using infinity norm:
    /// - omit primitive pair if \f$ {\rm ln_scr} \equiv (k_a k_b)  \max|c_a| \max|c_b| K_{ab}  < \epsilon \f$, where \f$ K_{ab} \equiv \sqrt{\max{|(ab|ab)|}} \f$
    /// - omit primitive 2-body integral if `spbra.primpairs[pb].ln_scr + spket.primpairs[pk].ln_scr < log(ε)`
    SchwarzInf = 0x1000,
    Invalid = 0x0000
  };
// clang-format on

namespace detail {
inline ScreeningMethod& default_screening_method_accessor() {
  static ScreeningMethod default_screening_method = ScreeningMethod::Original;
  return default_screening_method;
}
}  // namespace detail

inline ScreeningMethod default_screening_method() {
  return detail::default_screening_method_accessor();
}

inline void default_screening_method(ScreeningMethod screening_method) {
  detail::default_screening_method_accessor() = screening_method;
}

/// ShellPair contains pre-computed shell-pair data, primitive pairs are
/// screened to finite precision
struct ShellPair {
  typedef Shell::real_t real_t;

  // clang-format off
      /// PrimPairData contains pre-computed primitive pair data
      struct PrimPairData {
          real_t P[3]; //!< \f$ (\alpha_1 \vec{A} + \alpha_2 \vec{B})/(\alpha_1 + \alpha_2) \f$
          real_t K;    //!< \f$ \sqrt{2} \pi^{5/4} \exp(-|\vec{A}-\vec{B}|^2 \alpha_a * \alpha_b / (\alpha_a + \alpha_b)) / (\alpha_a + \alpha_b)  \f$
          real_t one_over_gamma; //!< \f$ 1 / (\alpha_a + \alpha_b)  \f$
          real_t nonsph_screen_fac; //!< used only when `screening_method_==ScreeningMethod::Conservative`: approximate upper bound for the modulation of the integrals due to nonspherical bra: \f$  \max( (\max_i|PA_i|)^{L_a} (\max_i|PB_i|)^{L_b}, L_a! L_b! / (\alpha_a + \alpha_b)^{L_a+L_b} ) \f$
          real_t ln_scr; //!< natural log of the primitive pair screening factor (see ScreeningMethod )
          int p1;  //!< first primitive index
          int p2;  //!< second primitive index
      };
  // clang-format on

  std::vector<PrimPairData> primpairs;
  real_t AB[3];
  real_t ln_prec = std::numeric_limits<real_t>::lowest();
  ScreeningMethod screening_method_ = ScreeningMethod::Invalid;

  ShellPair() : primpairs() {
    for (int i = 0; i != 3; ++i) AB[i] = 0.;
  }

  ShellPair(size_t max_nprim) : primpairs() {
    primpairs.reserve(max_nprim * max_nprim);
    for (int i = 0; i != 3; ++i) AB[i] = 0.;
  }
  template <typename Real>
  ShellPair(const Shell& s1, const Shell& s2, Real ln_prec,
            ScreeningMethod screening_method = default_screening_method()) {
    init(s1, s2, ln_prec, screening_method);
  }
  template <typename Real, typename SchwarzFactorEvaluator>
  ShellPair(const Shell& s1, const Shell& s2, Real ln_prec,
            ScreeningMethod screening_method,
            SchwarzFactorEvaluator&& schwarz_factor_evaluator) {
    init(s1, s2, ln_prec, screening_method,
         std::forward<SchwarzFactorEvaluator>(schwarz_factor_evaluator));
  }

  /// makes this equivalent to a default-initialized ShellPair, however the
  /// memory allocated in primpairs is not released
  void reset() {
    primpairs.clear();
    for (int i = 0; i != 3; ++i) AB[i] = 0.;
    ln_prec = std::numeric_limits<real_t>::lowest();
    screening_method_ = ScreeningMethod::Invalid;
  }

  void resize(std::size_t max_nprim) {
    const auto max_nprim2 = max_nprim * max_nprim;
    if (max_nprim * max_nprim > primpairs.size()) primpairs.resize(max_nprim2);
  }

  /// initializes the shell pair data using original or conservative screening
  /// methods
  template <typename Real>
  void init(const Shell& s1, const Shell& s2, Real ln_prec,
            ScreeningMethod screening_method = ScreeningMethod::Original) {
    assert(screening_method == ScreeningMethod::Original ||
           screening_method == ScreeningMethod::Conservative);

    using std::log;

    primpairs.clear();

    const auto& A = s1.O;
    const auto& B = s2.O;
    real_t AB2 = 0.;
    for (int i = 0; i != 3; ++i) {
      AB[i] = A[i] - B[i];
      AB2 += AB[i] * AB[i];
    }

    auto max_l = [](const Shell& s) {
      using std::begin;
      using std::end;
      return std::max_element(
                 begin(s.contr), end(s.contr),
                 [](const Shell::Contraction& c1,
                    const Shell::Contraction& c2) { return c1.l < c2.l; })
          ->l;
    };
    const auto max_l1 = max_l(s1);
    const auto max_l2 = max_l(s2);

    const auto nprim1 = s1.alpha.size();
    const auto nprim2 = s2.alpha.size();
    size_t c = 0;
    for (size_t p1 = 0; p1 != nprim1; ++p1) {
      for (size_t p2 = 0; p2 != nprim2; ++p2) {
        const auto& a1 = s1.alpha[p1];
        const auto& a2 = s2.alpha[p2];
        const auto gamma = a1 + a2;
        const auto oogamma = 1 / gamma;

        const auto rho = a1 * a2 * oogamma;
        const auto minus_rho_times_AB2 = -rho * AB2;
        real_t ln_screen_fac =
            minus_rho_times_AB2 + s1.max_ln_coeff[p1] + s2.max_ln_coeff[p2];
        if (screening_method == ScreeningMethod::Original &&
            ln_screen_fac < ln_prec)
          continue;

        real_t P[3];
        if (AB2 == 0.) {  // this buys a bit more precision
          P[0] = A[0];
          P[1] = A[1];
          P[2] = A[2];
        } else {
          P[0] = (a1 * A[0] + a2 * B[0]) * oogamma;
          P[1] = (a1 * A[1] + a2 * B[1]) * oogamma;
          P[2] = (a1 * A[2] + a2 * B[2]) * oogamma;
        }

        // conservative screening:
        // - partitions the error among all primitive pairs (use \epsilon /
        // nprim to screen, instead of \epsilon itself), and
        // - accounts for the proper spherical gaussian prefactor in the
        // integrals (namely, adds extra \sqrt{2 \pi^{52}}/\gamma_{ab} factor)
        // - accounts for the nonspherical gaussians ... namely
        //   magnitude of primitive (ab|00) integral for nonzero L differs from
        //   that of (00|00) by the magnitude of:
        //   - (max_i|PA_i|)^La (max_i|PB_i|)^Lb when A-B separation is large,
        //   or
        //   - La! Lb! / gammap^(La+Lb) when the separation is small
        real_t nonspherical_scr_factor = 0;
        if (screening_method == ScreeningMethod::Conservative) {
          const auto maxabs_PA_i_to_l1 = std::pow(
              std::max(std::max(std::abs(P[0] - A[0]), std::abs(P[1] - A[1])),
                       std::abs(P[2] - A[2])),
              max_l1);
          const auto maxabs_PB_i_to_l2 = std::pow(
              std::max(std::max(std::abs(P[0] - B[0]), std::abs(P[1] - B[1])),
                       std::abs(P[2] - B[2])),
              max_l2);
          const auto fac_l1_fac_l2_oogamma_to_l =
              math::fac_real<Real>(max_l1) * math::fac_real<Real>(max_l2) *
              std::pow(oogamma, max_l1 + max_l2);
          nonspherical_scr_factor =
              std::max(maxabs_PA_i_to_l1 * maxabs_PB_i_to_l2,
                       fac_l1_fac_l2_oogamma_to_l);
          const auto ln_nonspherical_scr_factor =
              log(std::max(nonspherical_scr_factor, static_cast<real_t>(1)));

          constexpr decltype(rho) ln_sqrt_two_times_M_PI_to_1pt25 =
              1.777485947591722872387900;  // \ln(\sqrt{2} (\pi)^{5/4})
          const auto ln_spherical_scr_extra_factor =
              ln_sqrt_two_times_M_PI_to_1pt25 + log(oogamma);
          const auto ln_nprim = log(nprim1 * nprim2);
          ln_screen_fac += ln_spherical_scr_extra_factor +
                           ln_nonspherical_scr_factor + ln_nprim;
          if (ln_screen_fac < ln_prec) continue;
        }

        primpairs.resize(c + 1);
        PrimPairData& p = primpairs[c];
        p.ln_scr = ln_screen_fac;
        p.p1 = p1;
        p.p2 = p2;
        constexpr decltype(rho) sqrt_two_times_M_PI_to_1pt25 =
            5.9149671727956128778;  // \sqrt{2} (\pi)^{5/4}
        p.K = sqrt_two_times_M_PI_to_1pt25 * exp(minus_rho_times_AB2) * oogamma;
        p.P[0] = P[0];
        p.P[1] = P[1];
        p.P[2] = P[2];
        p.nonsph_screen_fac = nonspherical_scr_factor;
        p.one_over_gamma = oogamma;

        ++c;
      }
    }

    this->ln_prec = ln_prec;
    this->screening_method_ = screening_method;
  }

  /// initializes the shell pair data using Schwarz screening methods
  template <typename Real, typename SchwarzFactorEvaluator>
  void init(const Shell& s1, const Shell& s2, Real ln_prec,
            ScreeningMethod screening_method,
            SchwarzFactorEvaluator&& schwarz_factor_evaluator) {
    assert(screening_method == ScreeningMethod::Schwarz ||
           screening_method == ScreeningMethod::SchwarzInf);

    using std::log;

    primpairs.clear();

    const auto& A = s1.O;
    const auto& B = s2.O;
    real_t AB2 = 0.;
    for (int i = 0; i != 3; ++i) {
      AB[i] = A[i] - B[i];
      AB2 += AB[i] * AB[i];
    }

    const auto nprim1 = s1.alpha.size();
    const auto nprim2 = s2.alpha.size();
    const auto nprim12 = nprim1 * nprim2;
    size_t c = 0;
    for (size_t p1 = 0; p1 != nprim1; ++p1) {
      for (size_t p2 = 0; p2 != nprim2; ++p2) {
        const auto ln_screen_fac =
            log(nprim12 * schwarz_factor_evaluator(s1, p1, s2, p2)) +
            s1.max_ln_coeff[p1] + s2.max_ln_coeff[p2];
        if (ln_screen_fac < ln_prec) continue;

        const auto& a1 = s1.alpha[p1];
        const auto& a2 = s2.alpha[p2];
        const auto gamma = a1 + a2;
        const auto oogamma = 1 / gamma;

        const auto rho = a1 * a2 * oogamma;
        const auto minus_rho_times_AB2 = -rho * AB2;

        real_t P[3];
        if (AB2 == 0.) {  // this buys a bit more precision
          P[0] = A[0];
          P[1] = A[1];
          P[2] = A[2];
        } else {
          P[0] = (a1 * A[0] + a2 * B[0]) * oogamma;
          P[1] = (a1 * A[1] + a2 * B[1]) * oogamma;
          P[2] = (a1 * A[2] + a2 * B[2]) * oogamma;
        }

        primpairs.resize(c + 1);
        PrimPairData& p = primpairs[c];
        p.ln_scr = ln_screen_fac;
        p.p1 = p1;
        p.p2 = p2;
        constexpr decltype(rho) sqrt_two_times_M_PI_to_1pt25 =
            5.9149671727956128778;  // \sqrt{2} (\pi)^{5/4}
        p.K = sqrt_two_times_M_PI_to_1pt25 * exp(minus_rho_times_AB2) * oogamma;
        p.P[0] = P[0];
        p.P[1] = P[1];
        p.P[2] = P[2];
        p.nonsph_screen_fac = 0;
        p.one_over_gamma = oogamma;

        ++c;
      }
    }

    this->ln_prec = ln_prec;
    this->screening_method_ = screening_method;
  }
};

}  // namespace libint2

#endif /* _libint2_src_lib_libint_shell_h_ */
