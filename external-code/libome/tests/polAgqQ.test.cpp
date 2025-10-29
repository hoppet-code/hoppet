/******************************************************************************
 * This file is part of libome                                                *
 * Copyright (C) 2025 Arnd Behring, Kay Schoenwald                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 ******************************************************************************/

#include <limits>
#include <gtest/gtest.h>
#include <ome/integration_engine_gsl.h>
#include <ome/mellin.h>
#include <ome/polAgqQ.h>

using namespace ome;

const double eps = std::numeric_limits<double>::epsilon();
const double testAs = 0.25;
const double testLM = -5.;
const double testNF = 3.;



TEST(polAgqQXspace,FullPowersOf4RegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto& eval = polAgqQ_reg;

  {
    testVal  = eval(testAs,testLM,testNF,9.31322574615478515625e-10);
    refVal   = -5864.5142543212087412144437886;
    errBound = 5864.5142543212087412144437886 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullPowersOf4, reg, x = 4^-15 = 9.31322574615478515625e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.7252902984619140625e-9);
    refVal   = -3125.04674983650805342296463118;
    errBound = 3125.04674983650805342296463118 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullPowersOf4, reg, x = 4^-14 = 3.7252902984619140625e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.490116119384765625e-8);
    refVal   = -1162.88394591719136040640158731;
    errBound = 1162.88394591719136040640158731 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullPowersOf4, reg, x = 4^-13 = 1.490116119384765625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,5.9604644775390625e-8);
    refVal   = 161.476433990065407314628160465;
    errBound = 161.476433990065407314628160465 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullPowersOf4, reg, x = 4^-12 = 5.9604644775390625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,2.384185791015625e-7);
    refVal   = 975.134411184787801504518441785;
    errBound = 975.134411184787801504518441785 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullPowersOf4, reg, x = 4^-11 = 2.384185791015625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,9.5367431640625e-7);
    refVal   = 1392.7878850744338955843824548;
    errBound = 1392.7878850744338955843824548 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullPowersOf4, reg, x = 4^-10 = 9.5367431640625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.814697265625e-6);
    refVal   = 1516.73266527659648318024397577;
    errBound = 1516.73266527659648318024397577 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullPowersOf4, reg, x = 4^-9 = 3.814697265625e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0000152587890625);
    refVal   = 1436.86144201094908568454902057;
    errBound = 1436.86144201094908568454902057 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullPowersOf4, reg, x = 4^-8 = 0.0000152587890625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00006103515625);
    refVal   = 1230.65755535976387962158965494;
    errBound = 1230.65755535976387962158965494 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullPowersOf4, reg, x = 4^-7 = 0.00006103515625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.000244140625);
    refVal   = 963.168670014583159977730677508;
    errBound = 963.168670014583159977730677508 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullPowersOf4, reg, x = 4^-6 = 0.000244140625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0009765625);
    refVal   = 686.913712505116993727973216291;
    errBound = 686.913712505116993727973216291 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullPowersOf4, reg, x = 4^-5 = 0.0009765625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00390625);
    refVal   = 441.596768861523741448862247604;
    errBound = 441.596768861523741448862247604 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullPowersOf4, reg, x = 4^-4 = 0.00390625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.015625);
    refVal   = 253.348019781124276314971278223;
    errBound = 253.348019781124276314971278223 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullPowersOf4, reg, x = 4^-3 = 0.015625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0625);
    refVal   = 133.07279798645902377019309066;
    errBound = 133.07279798645902377019309066 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullPowersOf4, reg, x = 4^-2 = 0.0625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9375);
    refVal   = 54.3014786886121466209010309635;
    errBound = 54.3014786886121466209010309635 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullPowersOf4, reg, x = 1-4^-2 = 0.9375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.984375);
    refVal   = 56.3376274624626634801937297186;
    errBound = 56.3376274624626634801937297186 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullPowersOf4, reg, x = 1-4^-3 = 0.984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99609375);
    refVal   = 53.137160355368286273136755414;
    errBound = 53.137160355368286273136755414 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullPowersOf4, reg, x = 1-4^-4 = 0.99609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9990234375);
    refVal   = 40.9419020999756532974699534511;
    errBound = 40.9419020999756532974699534511 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullPowersOf4, reg, x = 1-4^-5 = 0.9990234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999755859375);
    refVal   = 16.4936821067731044144988291174;
    errBound = 16.4936821067731044144988291174 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullPowersOf4, reg, x = 1-4^-6 = 0.999755859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99993896484375);
    refVal   = -23.7917571613185626821020315381;
    errBound = 23.7917571613185626821020315381 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullPowersOf4, reg, x = 1-4^-7 = 0.99993896484375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999847412109375);
    refVal   = -84.0352152287090219278269991027;
    errBound = 84.0352152287090219278269991027 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullPowersOf4, reg, x = 1-4^-8 = 0.9999847412109375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999996185302734375);
    refVal   = -168.932286811602588563093679449;
    errBound = 168.932286811602588563093679449 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullPowersOf4, reg, x = 1-4^-9 = 0.999996185302734375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999904632568359375);
    refVal   = -283.7549013975768698030234838;
    errBound = 283.7549013975768698030234838 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullPowersOf4, reg, x = 1-4^-10 = 0.99999904632568359375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999997615814208984375);
    refVal   = -434.348392956537857596546628632;
    errBound = 434.348392956537857596546628632 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullPowersOf4, reg, x = 1-4^-11 = 0.9999997615814208984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999940395355224609375);
    refVal   = -627.129546832176406041938378165;
    errBound = 627.129546832176406041938378165 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullPowersOf4, reg, x = 1-4^-12 = 0.999999940395355224609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999998509883880615234375);
    refVal   = -869.085689400231211517321756791;
    errBound = 869.085689400231211517321756791 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullPowersOf4, reg, x = 1-4^-13 = 0.99999998509883880615234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999962747097015380859375);
    refVal   = -1167.77432037836458951995436224;
    errBound = 1167.77432037836458951995436224 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullPowersOf4, reg, x = 1-4^-14 = 0.9999999962747097015380859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999068677425384521484375);
    refVal   = -1531.32297575340645293996717891;
    errBound = 1531.32297575340645293996717891 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullPowersOf4, reg, x = 1-4^-15 = 0.999999999068677425384521484375";
  }
}



TEST(polAgqQXspace,FullNormalDoublesRegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto& eval = polAgqQ_reg;

  {
    testVal  = eval(testAs,testLM,testNF,1.e-10);
    refVal   = -12295.2932998604657138311494396;
    errBound = 12295.2932998604657138311494396 * std::max(eps/(1.-1.e-10),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullNormalDoubles, reg, x = 10^-10 = 1.e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-9);
    refVal   = -5702.64114785747835997444773944;
    errBound = 5702.64114785747835997444773944 * std::max(eps/(1.-1.e-9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullNormalDoubles, reg, x = 10^-9 = 1.e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-8);
    refVal   = -1656.16661025819369445329791964;
    errBound = 1656.16661025819369445329791964 * std::max(eps/(1.-1.e-8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullNormalDoubles, reg, x = 10^-8 = 1.e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-7);
    refVal   = 518.499982321423315250145097263;
    errBound = 518.499982321423315250145097263 * std::max(eps/(1.-1.e-7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullNormalDoubles, reg, x = 10^-7 = 1.e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-6);
    refVal   = 1401.33480077093327057734004037;
    errBound = 1401.33480077093327057734004037 * std::max(eps/(1.-1.e-6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullNormalDoubles, reg, x = 10^-6 = 1.e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-5);
    refVal   = 1477.92123079632944361880799917;
    errBound = 1477.92123079632944361880799917 * std::max(eps/(1.-1.e-5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullNormalDoubles, reg, x = 10^-5 = 1.e-5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0001);
    refVal   = 1139.4290898075323485028016202;
    errBound = 1139.4290898075323485028016202 * std::max(eps/(1.-0.0001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullNormalDoubles, reg, x = 10^-4 = 0.0001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.001);
    refVal   = 682.362452266163818357719600591;
    errBound = 682.362452266163818357719600591 * std::max(eps/(1.-0.001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullNormalDoubles, reg, x = 10^-3 = 0.001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.01);
    refVal   = 306.723560378491719133575541367;
    errBound = 306.723560378491719133575541367 * std::max(eps/(1.-0.01),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullNormalDoubles, reg, x = 10^-2 = 0.01";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.05);
    refVal   = 147.895509156925415611010424994;
    errBound = 147.895509156925415611010424994 * std::max(eps/(1.-0.05),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullNormalDoubles, reg, x = 0.05 = 0.05";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.1);
    refVal   = 107.076513522237357530241052406;
    errBound = 107.076513522237357530241052406 * std::max(eps/(1.-0.1),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullNormalDoubles, reg, x = 0.1 = 0.1";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.15);
    refVal   = 89.8367754709730851110121580964;
    errBound = 89.8367754709730851110121580964 * std::max(eps/(1.-0.15),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullNormalDoubles, reg, x = 0.15 = 0.15";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.2);
    refVal   = 80.1415817381059121838898042215;
    errBound = 80.1415817381059121838898042215 * std::max(eps/(1.-0.2),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullNormalDoubles, reg, x = 0.2 = 0.2";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.25);
    refVal   = 73.8663463028922428001820071232;
    errBound = 73.8663463028922428001820071232 * std::max(eps/(1.-0.25),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullNormalDoubles, reg, x = 0.25 = 0.25";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.3);
    refVal   = 69.4321047732760340441941511182;
    errBound = 69.4321047732760340441941511182 * std::max(eps/(1.-0.3),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullNormalDoubles, reg, x = 0.3 = 0.3";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.35);
    refVal   = 66.0990862294612236273236266188;
    errBound = 66.0990862294612236273236266188 * std::max(eps/(1.-0.35),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullNormalDoubles, reg, x = 0.35 = 0.35";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.4);
    refVal   = 63.4749886562811730238842098138;
    errBound = 63.4749886562811730238842098138 * std::max(eps/(1.-0.4),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullNormalDoubles, reg, x = 0.4 = 0.4";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.45);
    refVal   = 61.3341680544947518943787789744;
    errBound = 61.3341680544947518943787789744 * std::max(eps/(1.-0.45),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullNormalDoubles, reg, x = 0.45 = 0.45";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.5);
    refVal   = 59.5400164553177613699271831444;
    errBound = 59.5400164553177613699271831444 * std::max(eps/(1.-0.5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullNormalDoubles, reg, x = 0.5 = 0.5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.55);
    refVal   = 58.0080125021197086787370151817;
    errBound = 58.0080125021197086787370151817 * std::max(eps/(1.-0.55),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullNormalDoubles, reg, x = 0.55 = 0.55";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.6);
    refVal   = 56.6870603349778285566078352629;
    errBound = 56.6870603349778285566078352629 * std::max(eps/(1.-0.6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullNormalDoubles, reg, x = 0.6 = 0.6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.65);
    refVal   = 55.55017510252496346525002805;
    errBound = 55.55017510252496346525002805 * std::max(eps/(1.-0.65),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullNormalDoubles, reg, x = 0.65 = 0.65";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.7);
    refVal   = 54.5909044859470685221982246629;
    errBound = 54.5909044859470685221982246629 * std::max(eps/(1.-0.7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullNormalDoubles, reg, x = 0.7 = 0.7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.75);
    refVal   = 53.8246209732169010367379256209;
    errBound = 53.8246209732169010367379256209 * std::max(eps/(1.-0.75),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullNormalDoubles, reg, x = 0.75 = 0.75";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.8);
    refVal   = 53.2963284053901378830033164373;
    errBound = 53.2963284053901378830033164373 * std::max(eps/(1.-0.8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullNormalDoubles, reg, x = 0.8 = 0.8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.85);
    refVal   = 53.1014319655211460134947237931;
    errBound = 53.1014319655211460134947237931 * std::max(eps/(1.-0.85),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullNormalDoubles, reg, x = 0.85 = 0.85";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9);
    refVal   = 53.4401718313850076803343985295;
    errBound = 53.4401718313850076803343985295 * std::max(eps/(1.-0.9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullNormalDoubles, reg, x = 0.9 = 0.9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.95);
    refVal   = 54.7683881290352012100537972805;
    errBound = 54.7683881290352012100537972805 * std::max(eps/(1.-0.95),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullNormalDoubles, reg, x = 0.95 = 0.95";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99);
    refVal   = 56.0864002000572412331476455484;
    errBound = 56.0864002000572412331476455484 * std::max(eps/(1.-0.99),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullNormalDoubles, reg, x = 1-10^-2 = 0.99";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999);
    refVal   = 41.24401935348072909876250687;
    errBound = 41.24401935348072909876250687 * std::max(eps/(1.-0.999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullNormalDoubles, reg, x = 1-10^-3 = 0.999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999);
    refVal   = -7.38183169342911794985678099919;
    errBound = 7.38183169342911794985678099919 * std::max(eps/(1.-0.9999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullNormalDoubles, reg, x = 1-10^-4 = 0.9999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999);
    refVal   = -107.069732216635878089170433553;
    errBound = 107.069732216635878089170433553 * std::max(eps/(1.-0.99999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullNormalDoubles, reg, x = 1-10^-5 = 0.99999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999);
    refVal   = -279.270035010546256637791435863;
    errBound = 279.270035010546256637791435863 * std::max(eps/(1.-0.999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullNormalDoubles, reg, x = 1-10^-6 = 0.999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999);
    refVal   = -549.808320031093197978294956775;
    errBound = 549.808320031093197978294956775 * std::max(eps/(1.-0.9999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullNormalDoubles, reg, x = 1-10^-7 = 0.9999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999999);
    refVal   = -948.8621604022100903280776822;
    errBound = 948.8621604022100903280776822 * std::max(eps/(1.-0.99999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullNormalDoubles, reg, x = 1-10^-8 = 0.99999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999);
    refVal   = -1510.95019198122283604307030278;
    errBound = 1510.95019198122283604307030278 * std::max(eps/(1.-0.999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullNormalDoubles, reg, x = 1-10^-9 = 0.999999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999999);
    refVal   = -2274.92959317430684348602023851;
    errBound = 2274.92959317430684348602023851 * std::max(eps/(1.-0.9999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, FullNormalDoubles, reg, x = 1-10^-10 = 0.9999999999";
  }
}



TEST(polAgqQXspace,TruncAs2PowersOf4RegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto eval = polAgqQ_reg.truncate(2);

  {
    testVal  = eval(testAs,testLM,testNF,9.31322574615478515625e-10);
    refVal   = 4.61728394943300956561241620867;
    errBound = 4.61728394943300956561241620867 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2PowersOf4, reg, x = 4^-15 = 9.31322574615478515625e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.7252902984619140625e-9);
    refVal   = 4.61728394588018641059781299119;
    errBound = 4.61728394588018641059781299119 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2PowersOf4, reg, x = 4^-14 = 3.7252902984619140625e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.490116119384765625e-8);
    refVal   = 4.61728393166889379053940064904;
    errBound = 4.61728393166889379053940064904 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2PowersOf4, reg, x = 4^-13 = 1.490116119384765625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,5.9604644775390625e-8);
    refVal   = 4.61728387482372331030578505698;
    errBound = 4.61728387482372331030578505698 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2PowersOf4, reg, x = 4^-12 = 5.9604644775390625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,2.384185791015625e-7);
    refVal   = 4.617283647443041389373484388;
    errBound = 4.617283647443041389373484388 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2PowersOf4, reg, x = 4^-11 = 2.384185791015625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,9.5367431640625e-7);
    refVal   = 4.61728273792031370578263057448;
    errBound = 4.61728273792031370578263057448 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2PowersOf4, reg, x = 4^-10 = 9.5367431640625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.814697265625e-6);
    refVal   = 4.61727909982940298027357087263;
    errBound = 4.61727909982940298027357087263 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2PowersOf4, reg, x = 4^-9 = 3.814697265625e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0000152587890625);
    refVal   = 4.61726454746576064492334708225;
    errBound = 4.61726454746576064492334708225 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2PowersOf4, reg, x = 4^-8 = 0.0000152587890625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00006103515625);
    refVal   = 4.61720633801122757328599014306;
    errBound = 4.61720633801122757328599014306 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2PowersOf4, reg, x = 4^-7 = 0.00006103515625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.000244140625);
    refVal   = 4.61697350019541702750161834139;
    errBound = 4.61697350019541702750161834139 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2PowersOf4, reg, x = 4^-6 = 0.000244140625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0009765625);
    refVal   = 4.61604214908088818879713982257;
    errBound = 4.61604214908088818879713982257 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2PowersOf4, reg, x = 4^-5 = 0.0009765625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00390625);
    refVal   = 4.61231675417175069421594693764;
    errBound = 4.61231675417175069421594693764 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2PowersOf4, reg, x = 4^-4 = 0.00390625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.015625);
    refVal   = 4.59741579380119943604803859767;
    errBound = 4.59741579380119943604803859767 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2PowersOf4, reg, x = 4^-3 = 0.015625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0625);
    refVal   = 4.5378537890749454776697437006;
    errBound = 4.5378537890749454776697437006 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2PowersOf4, reg, x = 4^-2 = 0.0625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9375);
    refVal   = 4.7006417491449966698290571176;
    errBound = 4.7006417491449966698290571176 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2PowersOf4, reg, x = 1-4^-2 = 0.9375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.984375);
    refVal   = 6.22522861949026531551228592284;
    errBound = 6.22522861949026531551228592284 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2PowersOf4, reg, x = 1-4^-3 = 0.984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99609375);
    refVal   = 8.15891776208110510632974224262;
    errBound = 8.15891776208110510632974224262 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2PowersOf4, reg, x = 1-4^-4 = 0.99609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9990234375);
    refVal   = 10.3724635201790262574495242021;
    errBound = 10.3724635201790262574495242021 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2PowersOf4, reg, x = 1-4^-5 = 0.9990234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999755859375);
    refVal   = 12.8209511398910947418951318607;
    errBound = 12.8209511398910947418951318607 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2PowersOf4, reg, x = 1-4^-6 = 0.999755859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99993896484375);
    refVal   = 15.4896492035196957875462267648;
    errBound = 15.4896492035196957875462267648 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2PowersOf4, reg, x = 1-4^-7 = 0.99993896484375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999847412109375);
    refVal   = 18.3739112007600045424826503462;
    errBound = 18.3739112007600045424826503462 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2PowersOf4, reg, x = 1-4^-8 = 0.9999847412109375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999996185302734375);
    refVal   = 21.4723125966721280433653393273;
    errBound = 21.4723125966721280433653393273 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2PowersOf4, reg, x = 1-4^-9 = 0.999996185302734375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999904632568359375);
    refVal   = 24.7844260323317759966106456931;
    errBound = 24.7844260323317759966106456931 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2PowersOf4, reg, x = 1-4^-10 = 0.99999904632568359375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999997615814208984375);
    refVal   = 28.3101254871234477243456553193;
    errBound = 28.3101254871234477243456553193 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2PowersOf4, reg, x = 1-4^-11 = 0.9999997615814208984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999940395355224609375);
    refVal   = 32.0493743170253863206604437516;
    errBound = 32.0493743170253863206604437516 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2PowersOf4, reg, x = 1-4^-12 = 0.999999940395355224609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999998509883880615234375);
    refVal   = 36.0021619904034884473864945475;
    errBound = 36.0021619904034884473864945475 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2PowersOf4, reg, x = 1-4^-13 = 0.99999998509883880615234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999962747097015380859375);
    refVal   = 40.1684855102141388226637669272;
    errBound = 40.1684855102141388226637669272 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2PowersOf4, reg, x = 1-4^-14 = 0.9999999962747097015380859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999068677425384521484375);
    refVal   = 44.5483440307931804093489845689;
    errBound = 44.5483440307931804093489845689 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2PowersOf4, reg, x = 1-4^-15 = 0.999999999068677425384521484375";
  }
}



TEST(polAgqQXspace,TruncAs2NormalDoublesRegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto eval = polAgqQ_reg.truncate(2);

  {
    testVal  = eval(testAs,testLM,testNF,1.e-10);
    refVal   = 4.61728395049012345679012345679;
    errBound = 4.61728395049012345679012345679 * std::max(eps/(1.-1.e-10),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2NormalDoubles, reg, x = 10^-10 = 1.e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-9);
    refVal   = 4.61728394934567901234567901252;
    errBound = 4.61728394934567901234567901252 * std::max(eps/(1.-1.e-9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2NormalDoubles, reg, x = 10^-9 = 1.e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-8);
    refVal   = 4.61728393790123456790123474074;
    errBound = 4.61728393790123456790123474074 * std::max(eps/(1.-1.e-8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2NormalDoubles, reg, x = 10^-8 = 1.e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-7);
    refVal   = 4.61728382345679012345696296298;
    errBound = 4.61728382345679012345696296298 * std::max(eps/(1.-1.e-7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2NormalDoubles, reg, x = 10^-7 = 1.e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-6);
    refVal   = 4.61728267901234567918518536728;
    errBound = 4.61728267901234567918518536728 * std::max(eps/(1.-1.e-6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2NormalDoubles, reg, x = 10^-6 = 1.e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-5);
    refVal   = 4.61727123456790140740922841201;
    errBound = 4.61727123456790140740922841201 * std::max(eps/(1.-1.e-5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2NormalDoubles, reg, x = 10^-5 = 1.e-5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0001);
    refVal   = 4.61715679012362964784120077162;
    errBound = 4.61715679012362964784120077162 * std::max(eps/(1.-0.0001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2NormalDoubles, reg, x = 10^-4 = 0.0001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.001);
    refVal   = 4.61601234585203412021618994331;
    errBound = 4.61601234585203412021618994331 * std::max(eps/(1.-0.001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2NormalDoubles, reg, x = 10^-3 = 0.001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.01);
    refVal   = 4.60456807591216191110458128969;
    errBound = 4.60456807591216191110458128969 * std::max(eps/(1.-0.01),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2NormalDoubles, reg, x = 10^-2 = 0.01";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.05);
    refVal   = 4.5537265022368824340734588473;
    errBound = 4.5537265022368824340734588473 * std::max(eps/(1.-0.05),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2NormalDoubles, reg, x = 0.05 = 0.05";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.1);
    refVal   = 4.49031637039527560713538209828;
    errBound = 4.49031637039527560713538209828 * std::max(eps/(1.-0.1),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2NormalDoubles, reg, x = 0.1 = 0.1";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.15);
    refVal   = 4.42723363382083587193555992282;
    errBound = 4.42723363382083587193555992282 * std::max(eps/(1.-0.15),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2NormalDoubles, reg, x = 0.15 = 0.15";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.2);
    refVal   = 4.36470333613255787835907196541;
    errBound = 4.36470333613255787835907196541 * std::max(eps/(1.-0.2),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2NormalDoubles, reg, x = 0.2 = 0.2";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.25);
    refVal   = 4.30300892177808931020793920638;
    errBound = 4.30300892177808931020793920638 * std::max(eps/(1.-0.25),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2NormalDoubles, reg, x = 0.25 = 0.25";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.3);
    refVal   = 4.24251073846190364331713676506;
    errBound = 4.24251073846190364331713676506 * std::max(eps/(1.-0.3),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2NormalDoubles, reg, x = 0.3 = 0.3";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.35);
    refVal   = 4.1836720042795251100605392741;
    errBound = 4.1836720042795251100605392741 * std::max(eps/(1.-0.35),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2NormalDoubles, reg, x = 0.35 = 0.35";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.4);
    refVal   = 4.12709604786965627549307627695;
    errBound = 4.12709604786965627549307627695 * std::max(eps/(1.-0.4),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2NormalDoubles, reg, x = 0.4 = 0.4";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.45);
    refVal   = 4.07358106141533380093308352096;
    errBound = 4.07358106141533380093308352096 * std::max(eps/(1.-0.45),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2NormalDoubles, reg, x = 0.45 = 0.45";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.5);
    refVal   = 4.02420296644421252729255117049;
    errBound = 4.02420296644421252729255117049 * std::max(eps/(1.-0.5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2NormalDoubles, reg, x = 0.5 = 0.5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.55);
    refVal   = 3.98044517509312558920755269379;
    errBound = 3.98044517509312558920755269379 * std::max(eps/(1.-0.55),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2NormalDoubles, reg, x = 0.55 = 0.55";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.6);
    refVal   = 3.94441023407143937169294940008;
    errBound = 3.94441023407143937169294940008 * std::max(eps/(1.-0.6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2NormalDoubles, reg, x = 0.6 = 0.6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.65);
    refVal   = 3.91918256710291119580898662975;
    errBound = 3.91918256710291119580898662975 * std::max(eps/(1.-0.65),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2NormalDoubles, reg, x = 0.65 = 0.65";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.7);
    refVal   = 3.90948977863511229883165024577;
    errBound = 3.90948977863511229883165024577 * std::max(eps/(1.-0.7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2NormalDoubles, reg, x = 0.7 = 0.7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.75);
    refVal   = 3.92300757769368870809254220628;
    errBound = 3.92300757769368870809254220628 * std::max(eps/(1.-0.75),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2NormalDoubles, reg, x = 0.75 = 0.75";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.8);
    refVal   = 3.97322168804105016995807403725;
    errBound = 3.97322168804105016995807403725 * std::max(eps/(1.-0.8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2NormalDoubles, reg, x = 0.8 = 0.8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.85);
    refVal   = 4.08671580348580889448617730769;
    errBound = 4.08671580348580889448617730769 * std::max(eps/(1.-0.85),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2NormalDoubles, reg, x = 0.85 = 0.85";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9);
    refVal   = 4.32657183933661878166662874387;
    errBound = 4.32657183933661878166662874387 * std::max(eps/(1.-0.9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2NormalDoubles, reg, x = 0.9 = 0.9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.95);
    refVal   = 4.90819964941440292807470513688;
    errBound = 4.90819964941440292807470513688 * std::max(eps/(1.-0.95),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2NormalDoubles, reg, x = 0.95 = 0.95";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99);
    refVal   = 6.81196394901942041243189231177;
    errBound = 6.81196394901942041243189231177 * std::max(eps/(1.-0.99),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2NormalDoubles, reg, x = 1-10^-2 = 0.99";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999);
    refVal   = 10.3325514507710507678649917672;
    errBound = 10.3325514507710507678649917672 * std::max(eps/(1.-0.999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2NormalDoubles, reg, x = 1-10^-3 = 0.999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999);
    refVal   = 14.5143595269692985793945546872;
    errBound = 14.5143595269692985793945546872 * std::max(eps/(1.-0.9999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2NormalDoubles, reg, x = 1-10^-4 = 0.9999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999);
    refVal   = 19.2957077649091804716682324662;
    errBound = 19.2957077649091804716682324662 * std::max(eps/(1.-0.99999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2NormalDoubles, reg, x = 1-10^-5 = 0.99999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999);
    refVal   = 24.6675703764789731320624136122;
    errBound = 24.6675703764789731320624136122 * std::max(eps/(1.-0.999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2NormalDoubles, reg, x = 1-10^-6 = 0.999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999);
    refVal   = 30.6287160628086868567820137636;
    errBound = 30.6287160628086868567820137636 * std::max(eps/(1.-0.9999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2NormalDoubles, reg, x = 1-10^-7 = 0.9999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999999);
    refVal   = 37.1789845274679714506918636079;
    errBound = 37.1789845274679714506918636079 * std::max(eps/(1.-0.99999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2NormalDoubles, reg, x = 1-10^-8 = 0.99999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999);
    refVal   = 44.3183555947057098230526412707;
    errBound = 44.3183555947057098230526412707 * std::max(eps/(1.-0.999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2NormalDoubles, reg, x = 1-10^-9 = 0.999999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999999);
    refVal   = 52.0468267893906073627928101848;
    errBound = 52.0468267893906073627928101848 * std::max(eps/(1.-0.9999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAgqQ, TruncAs2NormalDoubles, reg, x = 1-10^-10 = 0.9999999999";
  }
}



TEST(polAgqQNspace,FullMoments)
{
  double refVal = 0., errBound = 0.;
  auto& rpd = polAgqQ;
  integration_engine_gsl engine;
  auto mom = make_mellin_moment(rpd, engine, testAs, testLM, testNF);
  std::tuple<integration_status, double, double> res;

  {
    SCOPED_TRACE("polAgqQ Full Mellin moment N=3");
    refVal = 18.6631553335352302607848775339;
    errBound = 18.6631553335352302607848775339 * 3.e-11;
    res = mom.integrate(3, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAgqQ Full Mellin moment N=5");
    refVal = 10.8934663752502961030327661892;
    errBound = 10.8934663752502961030327661892 * 3.e-11;
    res = mom.integrate(5, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAgqQ Full Mellin moment N=7");
    refVal = 7.72882343371172501741001475838;
    errBound = 7.72882343371172501741001475838 * 3.e-11;
    res = mom.integrate(7, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAgqQ Full Mellin moment N=9");
    refVal = 6.00008396005396427186055102338;
    errBound = 6.00008396005396427186055102338 * 3.e-11;
    res = mom.integrate(9, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

}


TEST(polAgqQNspace,TruncAs2Moments)
{
  double refVal = 0., errBound = 0.;
  auto rpd = polAgqQ.truncate(2);
  integration_engine_gsl engine;
  auto mom = make_mellin_moment(rpd, engine, testAs, testLM, testNF);
  std::tuple<integration_status, double, double> res;

  {
    SCOPED_TRACE("polAgqQ TruncAs2 Mellin moment N=3");
    refVal = 1.46701388888888888888888888889;
    errBound = 1.46701388888888888888888888889 * 3.e-11;
    res = mom.integrate(3, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAgqQ TruncAs2 Mellin moment N=5");
    refVal = 0.923171810699588477366255144033;
    errBound = 0.923171810699588477366255144033 * 3.e-11;
    res = mom.integrate(5, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAgqQ TruncAs2 Mellin moment N=7");
    refVal = 0.688398905693229672821509556203;
    errBound = 0.688398905693229672821509556203 * 3.e-11;
    res = mom.integrate(7, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAgqQ TruncAs2 Mellin moment N=9");
    refVal = 0.555660581116305690085819421504;
    errBound = 0.555660581116305690085819421504 * 3.e-11;
    res = mom.integrate(9, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

}
