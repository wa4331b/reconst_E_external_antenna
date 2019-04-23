/**********************************************************************
 *
 * GL.cpp
 *
 * Copyright (C) 2014 Idesbald Van den Bosch
 *
 * This file is part of Puma-EM.
 * 
 * Puma-EM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Puma-EM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Puma-EM.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Suggestions/bugs : <vandenbosch.idesbald@gmail.com>
 *
 **********************************************************************/

#include <cstdlib>
#include <iostream>

using namespace std;

#include "GL.hpp"

const double XGL_1[1] = {0.0};
const double WGL_1[1] = {2.0};

const double XGL_2[2] = {-0.577350269189625764509148780502,
                          0.577350269189625764509148780502};
const double WGL_2[2] = {1.0, 1.0};

const double XGL_3[3] = {-0.774596669241483377035853079956,
                          0.0,
                          0.774596669241483377035853079956};
const double WGL_3[3] = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};

const double XGL_4[4] = {-0.861136311594052575223946488893,
                         -0.339981043584856264802665759103,
                          0.339981043584856264802665759103,
                          0.861136311594052575223946488893};
const double WGL_4[4] = {0.347854845137453857373063949222,
                         0.652145154862546142626936050778,
                         0.652145154862546142626936050778,
                         0.347854845137453857373063949222};

const double XGL_5[5] = {-0.906179845938663992797626878299,
                         -0.538469310105683091036314420700,
                          0.0,
                          0.538469310105683091036314420700,
                          0.906179845938663992797626878299};
const double WGL_5[5] = {0.236926885056189087514264040720,
                         0.478628670499366468041291514836,
                         0.568888888888888888888888888889,
                         0.478628670499366468041291514836,
                         0.236926885056189087514264040720};

const double XGL_6[6] = {-0.932469514203152027812301554494,
                         -0.661209386466264513661399595020,
                         -0.238619186083196908630501721681,
                          0.238619186083196908630501721681,
                          0.661209386466264513661399595020,
                          0.932469514203152027812301554494};
const double WGL_6[6] = {0.171324492379170345040296142173,
                         0.360761573048138607569833513838,
                         0.467913934572691047389870343990,
                         0.467913934572691047389870343990,
                         0.360761573048138607569833513838,
                         0.171324492379170345040296142173};

const double XGL_7[7] = {-0.949107912342758524526189684048,
                         -0.741531185599394439863864773281,
                         -0.405845151377397166906606412077,
                          0.0,
                          0.405845151377397166906606412077,
                          0.741531185599394439863864773281,
                          0.949107912342758524526189684048};
const double WGL_7[7] = {0.129484966168869693270611432679,
                         0.279705391489276667901467771424,
                         0.381830050505118944950369775489,
                         0.417959183673469387755102040816,
                         0.381830050505118944950369775489,
                         0.279705391489276667901467771424,
                         0.129484966168869693270611432679};

const double XGL_8[8] = {-0.960289856497536231683560868569,
                         -0.796666477413626739591553936476,
                         -0.525532409916328985817739049189,
                         -0.183434642495649804939476142360,
                          0.183434642495649804939476142360,
                          0.525532409916328985817739049189,
                          0.796666477413626739591553936476,
                          0.960289856497536231683560868569};
const double WGL_8[8] = {0.101228536290376259152531354310,
                         0.222381034453374470544355994426,
                         0.313706645877887287337962201987,
                         0.362683783378361982965150449277,
                         0.362683783378361982965150449277,
                         0.313706645877887287337962201987,
                         0.222381034453374470544355994426,
                         0.101228536290376259152531354310};

const double XGL_9[9] = {-0.968160239507626089835576202904,
          -0.836031107326635794299429788070,
          -0.613371432700590397308702039341,
          -0.324253423403808929038538014643,
            0.0,
            0.324253423403808929038538014643,
            0.613371432700590397308702039341,
            0.836031107326635794299429788070,
            0.968160239507626089835576202904};
const double WGL_9[9] = {0.812743883615744119718921581105E-01,
          0.180648160694857404058472031243,
          0.260610696402935462318742869419,
          0.312347077040002840068630406584,
          0.330239355001259763164525069287,
          0.312347077040002840068630406584,
          0.260610696402935462318742869419,
          0.180648160694857404058472031243,
          0.812743883615744119718921581105E-01};

const double XGL_10[10] = { -0.973906528517171720077964012084,
          -0.865063366688984510732096688423,
          -0.679409568299024406234327365115,
          -0.433395394129247190799265943166,
          -0.148874338981631210884826001130,
            0.148874338981631210884826001130,
            0.433395394129247190799265943166,
            0.679409568299024406234327365115,
            0.865063366688984510732096688423,
            0.973906528517171720077964012084};
const double WGL_10[10] = {0.666713443086881375935688098933E-01,
          0.149451349150580593145776339658,
          0.219086362515982043995534934228,
          0.269266719309996355091226921569,
          0.295524224714752870173892994651,
          0.295524224714752870173892994651,
          0.269266719309996355091226921569,
          0.219086362515982043995534934228,
          0.149451349150580593145776339658,
          0.666713443086881375935688098933E-01};

const double XGL_11[11] = {-0.978228658146056992803938001123,
          -0.887062599768095299075157769304,
          -0.730152005574049324093416252031,
          -0.519096129206811815925725669459,
          -0.269543155952344972331531985401,
            0.0,
            0.269543155952344972331531985401,
            0.519096129206811815925725669459,
            0.730152005574049324093416252031,
            0.887062599768095299075157769304,
            0.978228658146056992803938001123};
const double WGL_11[11] = {0.556685671161736664827537204425E-01,
          0.125580369464904624634694299224,
          0.186290210927734251426097641432,
          0.233193764591990479918523704843,
          0.262804544510246662180688869891,
          0.272925086777900630714483528336,
          0.262804544510246662180688869891,
          0.233193764591990479918523704843,
          0.186290210927734251426097641432,
          0.125580369464904624634694299224,
          0.556685671161736664827537204425E-01};

const double XGL_12[12] = { -0.981560634246719250690549090149,
          -0.904117256370474856678465866119,
          -0.769902674194304687036893833213,
          -0.587317954286617447296702418941,
          -0.367831498998180193752691536644,
          -0.125233408511468915472441369464,
            0.125233408511468915472441369464,
            0.367831498998180193752691536644,
            0.587317954286617447296702418941,
            0.769902674194304687036893833213,
            0.904117256370474856678465866119,
            0.981560634246719250690549090149};
const double WGL_12[12] = {0.471753363865118271946159614850E-01,
          0.106939325995318430960254718194,
          0.160078328543346226334652529543,
          0.203167426723065921749064455810,
          0.233492536538354808760849898925,
          0.249147045813402785000562436043,
          0.249147045813402785000562436043,
          0.233492536538354808760849898925,
          0.203167426723065921749064455810,
          0.160078328543346226334652529543,
          0.106939325995318430960254718194,
          0.471753363865118271946159614850E-01};

const double XGL_13[13] = {-0.984183054718588149472829448807,
          -0.917598399222977965206547836501,
          -0.801578090733309912794206489583,
          -0.642349339440340220643984606996,
          -0.448492751036446852877912852128,
          -0.230458315955134794065528121098,
            0.0,
            0.230458315955134794065528121098,
            0.448492751036446852877912852128,
            0.642349339440340220643984606996,
            0.801578090733309912794206489583,
            0.917598399222977965206547836501,
            0.984183054718588149472829448807};
const double WGL_13[13] = {0.404840047653158795200215922010E-01,
          0.921214998377284479144217759538E-01,
          0.138873510219787238463601776869,
          0.178145980761945738280046691996,
          0.207816047536888502312523219306,
          0.226283180262897238412090186040,
          0.232551553230873910194589515269,
          0.226283180262897238412090186040,
          0.207816047536888502312523219306,
          0.178145980761945738280046691996,
          0.138873510219787238463601776869,
          0.921214998377284479144217759538E-01,
          0.404840047653158795200215922010E-01};

const double XGL_14[14] = {-0.986283808696812338841597266704,
          -0.928434883663573517336391139378,
          -0.827201315069764993189794742650,
          -0.687292904811685470148019803019,
          -0.515248636358154091965290718551,
          -0.319112368927889760435671824168,
          -0.108054948707343662066244650220,
            0.108054948707343662066244650220,
            0.319112368927889760435671824168,
            0.515248636358154091965290718551,
            0.687292904811685470148019803019,
            0.827201315069764993189794742650,
            0.928434883663573517336391139378,
            0.986283808696812338841597266704};
const double WGL_14[14] = {0.351194603317518630318328761382E-01,
          0.801580871597602098056332770629E-01,
          0.121518570687903184689414809072,
          0.157203167158193534569601938624,
          0.185538397477937813741716590125,
          0.205198463721295603965924065661,
          0.215263853463157790195876443316,
          0.215263853463157790195876443316,
          0.205198463721295603965924065661,
          0.185538397477937813741716590125,
          0.157203167158193534569601938624,
          0.121518570687903184689414809072,
          0.801580871597602098056332770629E-01,
          0.351194603317518630318328761382E-01};

const double XGL_15[15] = {-0.987992518020485428489565718587,
          -0.937273392400705904307758947710,
          -0.848206583410427216200648320774,
          -0.724417731360170047416186054614,
          -0.570972172608538847537226737254,
          -0.394151347077563369897207370981,
          -0.201194093997434522300628303395,
           0.0,
            0.201194093997434522300628303395,
            0.394151347077563369897207370981,
            0.570972172608538847537226737254,
            0.724417731360170047416186054614,
            0.848206583410427216200648320774,
            0.937273392400705904307758947710,
            0.987992518020485428489565718587};
const double WGL_15[15] = {0.307532419961172683546283935772E-01,
          0.703660474881081247092674164507E-01,
          0.107159220467171935011869546686,
          0.139570677926154314447804794511,
          0.166269205816993933553200860481,
          0.186161000015562211026800561866,
          0.198431485327111576456118326444,
          0.202578241925561272880620199968,
          0.198431485327111576456118326444,
          0.186161000015562211026800561866,
          0.166269205816993933553200860481,
          0.139570677926154314447804794511,
          0.107159220467171935011869546686,
          0.703660474881081247092674164507E-01,
          0.307532419961172683546283935772E-01};

const double XGL_16[16] = {-0.989400934991649932596154173450,
          -0.944575023073232576077988415535,
          -0.865631202387831743880467897712,
          -0.755404408355003033895101194847,
          -0.617876244402643748446671764049,
          -0.458016777657227386342419442984,
          -0.281603550779258913230460501460,
          -0.950125098376374401853193354250E-01,
           0.950125098376374401853193354250E-01,
            0.281603550779258913230460501460,
            0.458016777657227386342419442984,
            0.617876244402643748446671764049,
            0.755404408355003033895101194847,
            0.865631202387831743880467897712,
            0.944575023073232576077988415535,
            0.989400934991649932596154173450};
const double WGL_16[16] = {0.271524594117540948517805724560E-01,
          0.622535239386478928628438369944E-01,
          0.951585116824927848099251076022E-01,
          0.124628971255533872052476282192,
          0.149595988816576732081501730547,
          0.169156519395002538189312079030,
          0.182603415044923588866763667969,
          0.189450610455068496285396723208,
          0.189450610455068496285396723208,
          0.182603415044923588866763667969,
          0.169156519395002538189312079030,
          0.149595988816576732081501730547,
          0.124628971255533872052476282192,
          0.951585116824927848099251076022E-01,
          0.622535239386478928628438369944E-01,
          0.271524594117540948517805724560E-01};

const double XGL_17[17] = {-0.990575475314417335675434019941,
    -0.950675521768767761222716957896,
    -0.880239153726985902122955694488,
    -0.781514003896801406925230055520,
    -0.657671159216690765850302216643,
    -0.512690537086476967886246568630,
    -0.351231763453876315297185517095,
    -0.178484181495847855850677493654,
     0.0,
      0.178484181495847855850677493654,
      0.351231763453876315297185517095,
      0.512690537086476967886246568630,
      0.657671159216690765850302216643,
      0.781514003896801406925230055520,
      0.880239153726985902122955694488,
      0.950675521768767761222716957896,
      0.990575475314417335675434019941};
const double WGL_17[17] = {0.241483028685479319601100262876E-01,
    0.554595293739872011294401653582E-01,
    0.850361483171791808835353701911E-01,
    0.111883847193403971094788385626,
    0.135136368468525473286319981702,
    0.154045761076810288081431594802,
    0.168004102156450044509970663788,
    0.176562705366992646325270990113,
    0.179446470356206525458265644262,
    0.176562705366992646325270990113,
    0.168004102156450044509970663788,
    0.154045761076810288081431594802,
    0.135136368468525473286319981702,
    0.111883847193403971094788385626,
    0.850361483171791808835353701911E-01,
    0.554595293739872011294401653582E-01,
    0.241483028685479319601100262876E-01};

const double XGL_18[18] = {-0.991565168420930946730016004706,
    -0.955823949571397755181195892930,
    -0.892602466497555739206060591127,
    -0.803704958972523115682417455015,
    -0.691687043060353207874891081289,
    -0.559770831073947534607871548525,
    -0.411751161462842646035931793833,
    -0.251886225691505509588972854878,
    -0.847750130417353012422618529358E-01,
     0.847750130417353012422618529358E-01,
      0.251886225691505509588972854878,
      0.411751161462842646035931793833,
      0.559770831073947534607871548525,
      0.691687043060353207874891081289,
      0.803704958972523115682417455015,
      0.892602466497555739206060591127,
      0.955823949571397755181195892930,
      0.991565168420930946730016004706};
const double WGL_18[18] = {0.216160135264833103133427102665E-01,
    0.497145488949697964533349462026E-01,
    0.764257302548890565291296776166E-01,
    0.100942044106287165562813984925,
    0.122555206711478460184519126800,
    0.140642914670650651204731303752,
    0.154684675126265244925418003836,
    0.164276483745832722986053776466,
    0.169142382963143591840656470135,
    0.169142382963143591840656470135,
    0.164276483745832722986053776466,
    0.154684675126265244925418003836,
    0.140642914670650651204731303752,
    0.122555206711478460184519126800,
    0.100942044106287165562813984925,
    0.764257302548890565291296776166E-01,
    0.497145488949697964533349462026E-01,
    0.216160135264833103133427102665E-01};

const double XGL_19[19] = {-0.992406843843584403189017670253,
    -0.960208152134830030852778840688,
    -0.903155903614817901642660928532,
    -0.822714656537142824978922486713,
    -0.720966177335229378617095860824,
    -0.600545304661681023469638164946,
    -0.464570741375960945717267148104,
    -0.316564099963629831990117328850,
    -0.160358645640225375868096115741,
     0.0,
      0.160358645640225375868096115741,
      0.316564099963629831990117328850,
      0.464570741375960945717267148104,
      0.600545304661681023469638164946,
      0.720966177335229378617095860824,
      0.822714656537142824978922486713,
      0.903155903614817901642660928532,
      0.960208152134830030852778840688,
      0.992406843843584403189017670253};
const double WGL_19[19] = {0.194617882297264770363120414644E-01,
    0.448142267656996003328381574020E-01,
    0.690445427376412265807082580060E-01,
    0.914900216224499994644620941238E-01,
    0.111566645547333994716023901682,
    0.128753962539336227675515784857,
    0.142606702173606611775746109442,
    0.152766042065859666778855400898,
    0.158968843393954347649956439465,
    0.161054449848783695979163625321,
    0.158968843393954347649956439465,
    0.152766042065859666778855400898,
    0.142606702173606611775746109442,
    0.128753962539336227675515784857,
    0.111566645547333994716023901682,
    0.914900216224499994644620941238E-01,
    0.690445427376412265807082580060E-01,
    0.448142267656996003328381574020E-01,
    0.194617882297264770363120414644E-01};

const double XGL_20[20] = {-0.993128599185094924786122388471,
    -0.963971927277913791267666131197,
    -0.912234428251325905867752441203,
    -0.839116971822218823394529061702,
    -0.746331906460150792614305070356,
    -0.636053680726515025452836696226,
    -0.510867001950827098004364050955,
    -0.373706088715419560672548177025,
    -0.227785851141645078080496195369,
    -0.765265211334973337546404093988E-01,
      0.765265211334973337546404093988E-01,
      0.227785851141645078080496195369,
      0.373706088715419560672548177025,
      0.510867001950827098004364050955,
      0.636053680726515025452836696226,
      0.746331906460150792614305070356,
      0.839116971822218823394529061702,
      0.912234428251325905867752441203,
      0.963971927277913791267666131197,
      0.993128599185094924786122388471};
const double WGL_20[20] = {0.176140071391521183118619623519E-01,
    0.406014298003869413310399522749E-01,
    0.626720483341090635695065351870E-01,
    0.832767415767047487247581432220E-01,
    0.101930119817240435036750135480,
    0.118194531961518417312377377711,
    0.131688638449176626898494499748,
    0.142096109318382051329298325067,
    0.149172986472603746787828737002,
    0.152753387130725850698084331955,
    0.152753387130725850698084331955,
    0.149172986472603746787828737002,
    0.142096109318382051329298325067,
    0.131688638449176626898494499748,
    0.118194531961518417312377377711,
    0.101930119817240435036750135480,
    0.832767415767047487247581432220E-01,
    0.626720483341090635695065351870E-01,
    0.406014298003869413310399522749E-01,
    0.176140071391521183118619623519E-01};

const double XGL_32[32] = {-0.997263861849481563544981128665,
    -0.985611511545268335400175044631,
    -0.964762255587506430773811928118,
    -0.934906075937739689170919134835,
    -0.896321155766052123965307243719,
    -0.849367613732569970133693004968,
    -0.794483795967942406963097298970,
    -0.732182118740289680387426665091,
    -0.663044266930215200975115168663,
    -0.587715757240762329040745476402,
    -0.506899908932229390023747474378,
    -0.421351276130635345364119436172,
    -0.331868602282127649779916805730,
    -0.239287362252137074544603209166,
    -0.144471961582796493485186373599,
    -0.483076656877383162348125704405E-01,
      0.483076656877383162348125704405E-01,
      0.144471961582796493485186373599,
      0.239287362252137074544603209166,
      0.331868602282127649779916805730,
      0.421351276130635345364119436172,
      0.506899908932229390023747474378,
      0.587715757240762329040745476402,
      0.663044266930215200975115168663,
      0.732182118740289680387426665091,
      0.794483795967942406963097298970,
      0.849367613732569970133693004968,
      0.896321155766052123965307243719,
      0.934906075937739689170919134835,
      0.964762255587506430773811928118,
      0.985611511545268335400175044631,
      0.997263861849481563544981128665};
const double WGL_32[32] = {0.701861000947009660040706373885E-02,
    0.162743947309056706051705622064E-01,
    0.253920653092620594557525897892E-01,
    0.342738629130214331026877322524E-01,
    0.428358980222266806568786466061E-01,
    0.509980592623761761961632446895E-01,
    0.586840934785355471452836373002E-01,
    0.658222227763618468376500637069E-01,
    0.723457941088485062253993564785E-01,
    0.781938957870703064717409188283E-01,
    0.833119242269467552221990746043E-01,
    0.876520930044038111427714627518E-01,
    0.911738786957638847128685771116E-01,
    0.938443990808045656391802376681E-01,
    0.956387200792748594190820022041E-01,
    0.965400885147278005667648300636E-01,
    0.965400885147278005667648300636E-01,
    0.956387200792748594190820022041E-01,
    0.938443990808045656391802376681E-01,
    0.911738786957638847128685771116E-01,
    0.876520930044038111427714627518E-01,
    0.833119242269467552221990746043E-01,
    0.781938957870703064717409188283E-01,
    0.723457941088485062253993564785E-01,
    0.658222227763618468376500637069E-01,
    0.586840934785355471452836373002E-01,
    0.509980592623761761961632446895E-01,
    0.428358980222266806568786466061E-01,
    0.342738629130214331026877322524E-01,
    0.253920653092620594557525897892E-01,
    0.162743947309056706051705622064E-01,
    0.701861000947009660040706373885E-02};

const double XGL_64[64] = {-0.999305041735772139456905624346,
    -0.996340116771955279346924500676,
    -0.991013371476744320739382383443,
    -0.983336253884625956931299302157,
    -0.973326827789910963741853507352,
    -0.961008799652053718918614121897,
    -0.946411374858402816062481491347,
    -0.929569172131939575821490154559,
    -0.910522137078502805756380668008,
    -0.889315445995114105853404038273,
    -0.865999398154092819760783385070,
    -0.840629296252580362751691544696,
    -0.813265315122797559741923338086,
    -0.783972358943341407610220525214,
    -0.752819907260531896611863774886,
    -0.719881850171610826848940217832,
    -0.685236313054233242563558371031,
    -0.648965471254657339857761231993,
    -0.611155355172393250248852971019,
    -0.571895646202634034283878116659,
    -0.531279464019894545658013903544,
    -0.489403145707052957478526307022,
    -0.446366017253464087984947714759,
    -0.402270157963991603695766771260,
    -0.357220158337668115950442615046,
    -0.311322871990210956157512698560,
    -0.264687162208767416373964172510,
    -0.217423643740007084149648748989,
    -0.169644420423992818037313629748,
    -0.121462819296120554470376463492,
    -0.729931217877990394495429419403E-01,
    -0.243502926634244325089558428537E-01,
     0.243502926634244325089558428537E-01,
      0.729931217877990394495429419403E-01,
      0.121462819296120554470376463492,
      0.169644420423992818037313629748,
      0.217423643740007084149648748989,
      0.264687162208767416373964172510,
      0.311322871990210956157512698560,
      0.357220158337668115950442615046,
      0.402270157963991603695766771260,
      0.446366017253464087984947714759,
      0.489403145707052957478526307022,
      0.531279464019894545658013903544,
      0.571895646202634034283878116659,
      0.611155355172393250248852971019,
      0.648965471254657339857761231993,
      0.685236313054233242563558371031,
      0.719881850171610826848940217832,
      0.752819907260531896611863774886,
      0.783972358943341407610220525214,
      0.813265315122797559741923338086,
      0.840629296252580362751691544696,
      0.865999398154092819760783385070,
      0.889315445995114105853404038273,
      0.910522137078502805756380668008,
      0.929569172131939575821490154559,
      0.946411374858402816062481491347,
      0.961008799652053718918614121897,
      0.973326827789910963741853507352,
      0.983336253884625956931299302157,
      0.991013371476744320739382383443,
      0.996340116771955279346924500676,
      0.999305041735772139456905624346};
const double WGL_64[64] = {0.178328072169643294729607914497E-02,
    0.414703326056246763528753572855E-02,
    0.650445796897836285611736039998E-02,
    0.884675982636394772303091465973E-02,
    0.111681394601311288185904930192E-01,
    0.134630478967186425980607666860E-01,
    0.157260304760247193219659952975E-01,
    0.179517157756973430850453020011E-01,
    0.201348231535302093723403167285E-01,
    0.222701738083832541592983303842E-01,
    0.243527025687108733381775504091E-01,
    0.263774697150546586716917926252E-01,
    0.283396726142594832275113052002E-01,
    0.302346570724024788679740598195E-01,
    0.320579283548515535854675043479E-01,
    0.338051618371416093915654821107E-01,
    0.354722132568823838106931467152E-01,
    0.370551285402400460404151018096E-01,
    0.385501531786156291289624969468E-01,
    0.399537411327203413866569261283E-01,
    0.412625632426235286101562974736E-01,
    0.424735151236535890073397679088E-01,
    0.435837245293234533768278609737E-01,
    0.445905581637565630601347100309E-01,
    0.454916279274181444797709969713E-01,
    0.462847965813144172959532492323E-01,
    0.469681828162100173253262857546E-01,
    0.475401657148303086622822069442E-01,
    0.479993885964583077281261798713E-01,
    0.483447622348029571697695271580E-01,
    0.485754674415034269347990667840E-01,
    0.486909570091397203833653907347E-01,
    0.486909570091397203833653907347E-01,
    0.485754674415034269347990667840E-01,
    0.483447622348029571697695271580E-01,
    0.479993885964583077281261798713E-01,
    0.475401657148303086622822069442E-01,
    0.469681828162100173253262857546E-01,
    0.462847965813144172959532492323E-01,
    0.454916279274181444797709969713E-01,
    0.445905581637565630601347100309E-01,
    0.435837245293234533768278609737E-01,
    0.424735151236535890073397679088E-01,
    0.412625632426235286101562974736E-01,
    0.399537411327203413866569261283E-01,
    0.385501531786156291289624969468E-01,
    0.370551285402400460404151018096E-01,
    0.354722132568823838106931467152E-01,
    0.338051618371416093915654821107E-01,
    0.320579283548515535854675043479E-01,
    0.302346570724024788679740598195E-01,
    0.283396726142594832275113052002E-01,
    0.263774697150546586716917926252E-01,
    0.243527025687108733381775504091E-01,
    0.222701738083832541592983303842E-01,
    0.201348231535302093723403167285E-01,
    0.179517157756973430850453020011E-01,
    0.157260304760247193219659952975E-01,
    0.134630478967186425980607666860E-01,
    0.111681394601311288185904930192E-01,
    0.884675982636394772303091465973E-02,
    0.650445796897836285611736039998E-02,
    0.414703326056246763528753572855E-02,
    0.178328072169643294729607914497E-02};

void Gauss_Legendre (const double * & XGL, const double * & WGL, const int N_points) {
  int N = N_points;
  if (N_points < 1) {
    std::cout << "Gauss_Legendre(): N_points < 1!" << endl;
    exit(EXIT_FAILURE);
  }
  else if (N_points > 20) {
    cout << "Warning: Gauss_Legendre(): N_points > 20! Next choices: 32 or 64..." << endl;
    if (N_points <= 32) N = 32;
    else if (N_points <= 64) N = 64;
    else {
      cout << "Gauss_Legendre(): N_points > 64!!" << endl;
      exit(EXIT_FAILURE);
    }
  }
  switch(N) {
    case 1: XGL = XGL_1; WGL = WGL_1; break;
    case 2: XGL = XGL_2; WGL = WGL_2; break;
    case 3: XGL = XGL_3; WGL = WGL_3; break;
    case 4: XGL = XGL_4; WGL = WGL_4; break;
    case 5: XGL = XGL_5; WGL = WGL_5; break;
    case 6: XGL = XGL_6; WGL = WGL_6; break;
    case 7: XGL = XGL_7; WGL = WGL_7; break;
    case 8: XGL = XGL_8; WGL = WGL_8; break;
    case 9: XGL = XGL_9; WGL = WGL_9; break;
    case 10: XGL = XGL_10; WGL = WGL_10; break;
    case 11: XGL = XGL_11; WGL = WGL_11; break;
    case 12: XGL = XGL_12; WGL = WGL_12; break;
    case 13: XGL = XGL_13; WGL = WGL_13; break;
    case 14: XGL = XGL_14; WGL = WGL_14; break;
    case 15: XGL = XGL_15; WGL = WGL_15; break;
    case 16: XGL = XGL_16; WGL = WGL_16; break;
    case 17: XGL = XGL_17; WGL = WGL_17; break;
    case 18: XGL = XGL_18; WGL = WGL_18; break;
    case 19: XGL = XGL_19; WGL = WGL_19; break;
    case 20: XGL = XGL_20; WGL = WGL_20; break;
    case 32: XGL = XGL_32; WGL = WGL_32; break;
    case 64: XGL = XGL_64; WGL = WGL_64; break;
    default: break;
  }
}
