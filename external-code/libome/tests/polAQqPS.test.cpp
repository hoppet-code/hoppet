/******************************************************************************
 * This file is part of libome                                                *
 * Copyright (C) 2025 Arnd Behring, Kay Schoenwald                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 ******************************************************************************/

#include <limits>
#include <gtest/gtest.h>
#include <ome/integration_engine_gsl.h>
#include <ome/mellin.h>
#include <ome/polAQqPS.h>

using namespace ome;

const double eps = std::numeric_limits<double>::epsilon();
const double testAs = 0.25;
const double testLM = -5.;
const double testNF = 3.;



TEST(polAQqPSXspace,FullPowersOf4RegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto& eval = polAQqPS_reg;

  {
    testVal  = eval(testAs,testLM,testNF,9.31322574615478515625e-10);
    refVal   = -14550.9779556345755310463050614;
    errBound = 14550.9779556345755310463050614 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullPowersOf4, reg, x = 4^-15 = 9.31322574615478515625e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.7252902984619140625e-9);
    refVal   = -13812.1953914756785783268086642;
    errBound = 13812.1953914756785783268086642 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullPowersOf4, reg, x = 4^-14 = 3.7252902984619140625e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.490116119384765625e-8);
    refVal   = -12017.8460942575964759283671119;
    errBound = 12017.8460942575964759283671119 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullPowersOf4, reg, x = 4^-13 = 1.490116119384765625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,5.9604644775390625e-8);
    refVal   = -9751.42790633823500743761503504;
    errBound = 9751.42790633823500743761503504 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullPowersOf4, reg, x = 4^-12 = 5.9604644775390625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,2.384185791015625e-7);
    refVal   = -7427.38966963621865387218589931;
    errBound = 7427.38966963621865387218589931 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullPowersOf4, reg, x = 4^-11 = 2.384185791015625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,9.5367431640625e-7);
    refVal   = -5315.59304057740049035860360784;
    errBound = 5315.59304057740049035860360784 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullPowersOf4, reg, x = 4^-10 = 9.5367431640625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.814697265625e-6);
    refVal   = -3565.77252154140963047568033758;
    errBound = 3565.77252154140963047568033758 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullPowersOf4, reg, x = 4^-9 = 3.814697265625e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0000152587890625);
    refVal   = -2231.99007364738097325866121262;
    errBound = 2231.99007364738097325866121262 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullPowersOf4, reg, x = 4^-8 = 0.0000152587890625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00006103515625);
    refVal   = -1297.07304286278653166765685949;
    errBound = 1297.07304286278653166765685949 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullPowersOf4, reg, x = 4^-7 = 0.00006103515625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.000244140625);
    refVal   = -696.997764938323050254402810672;
    errBound = 696.997764938323050254402810672 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullPowersOf4, reg, x = 4^-6 = 0.000244140625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0009765625);
    refVal   = -345.081068756373088443379449233;
    errBound = 345.081068756373088443379449233 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullPowersOf4, reg, x = 4^-5 = 0.0009765625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00390625);
    refVal   = -155.443419982509668452009519007;
    errBound = 155.443419982509668452009519007 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullPowersOf4, reg, x = 4^-4 = 0.00390625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.015625);
    refVal   = -63.655191951623183245197621411;
    errBound = 63.655191951623183245197621411 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullPowersOf4, reg, x = 4^-3 = 0.015625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0625);
    refVal   = -36.9717931674720967981249668138;
    errBound = 36.9717931674720967981249668138 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullPowersOf4, reg, x = 4^-2 = 0.0625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9375);
    refVal   = -5.91327709777633647816710181427;
    errBound = 5.91327709777633647816710181427 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullPowersOf4, reg, x = 1-4^-2 = 0.9375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.984375);
    refVal   = -1.22249483416429423621575432351;
    errBound = 1.22249483416429423621575432351 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullPowersOf4, reg, x = 1-4^-3 = 0.984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99609375);
    refVal   = -0.21218232083968820941000913462;
    errBound = 0.21218232083968820941000913462 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullPowersOf4, reg, x = 1-4^-4 = 0.99609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9990234375);
    refVal   = -0.0220888003641932255748993332029;
    errBound = 0.0220888003641932255748993332029 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullPowersOf4, reg, x = 1-4^-5 = 0.9990234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999755859375);
    refVal   = 0.00416319656298259802573147998768;
    errBound = 0.00416319656298259802573147998768 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullPowersOf4, reg, x = 1-4^-6 = 0.999755859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99993896484375);
    refVal   = 0.00392985861882642026861464095169;
    errBound = 0.00392985861882642026861464095169 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullPowersOf4, reg, x = 1-4^-7 = 0.99993896484375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999847412109375);
    refVal   = 0.00180732695850124044340991642679;
    errBound = 0.00180732695850124044340991642679 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullPowersOf4, reg, x = 1-4^-8 = 0.9999847412109375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999996185302734375);
    refVal   = 0.000677563171526391371583399763772;
    errBound = 0.000677563171526391371583399763772 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullPowersOf4, reg, x = 1-4^-9 = 0.999996185302734375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999904632568359375);
    refVal   = 0.000228520991401414467332000147153;
    errBound = 0.000228520991401414467332000147153 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullPowersOf4, reg, x = 1-4^-10 = 0.99999904632568359375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999997615814208984375);
    refVal   = 0.0000718789323272358368775234913458;
    errBound = 0.0000718789323272358368775234913458 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullPowersOf4, reg, x = 1-4^-11 = 0.9999997615814208984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999940395355224609375);
    refVal   = 0.0000214306392730148195231224839987;
    errBound = 0.0000214306392730148195231224839987 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullPowersOf4, reg, x = 1-4^-12 = 0.999999940395355224609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999998509883880615234375);
    refVal   = 6.10167996822887079872573392194e-6;
    errBound = 6.10167996822887079872573392194e-6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullPowersOf4, reg, x = 1-4^-13 = 0.99999998509883880615234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999962747097015380859375);
    refVal   = 1.66241553699523718770990669527e-6;
    errBound = 1.66241553699523718770990669527e-6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullPowersOf4, reg, x = 1-4^-14 = 0.9999999962747097015380859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999068677425384521484375);
    refVal   = 4.3228644609618155961889709983e-7;
    errBound = 4.3228644609618155961889709983e-7 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullPowersOf4, reg, x = 1-4^-15 = 0.999999999068677425384521484375";
  }
}



TEST(polAQqPSXspace,FullNormalDoublesRegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto& eval = polAQqPS_reg;

  {
    testVal  = eval(testAs,testLM,testNF,1.e-10);
    refVal   = -11476.3569222839619080870993685;
    errBound = 11476.3569222839619080870993685 * std::max(eps/(1.-1.e-10),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullNormalDoubles, reg, x = 10^-10 = 1.e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-9);
    refVal   = -14550.1813532716533368029585726;
    errBound = 14550.1813532716533368029585726 * std::max(eps/(1.-1.e-9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullNormalDoubles, reg, x = 10^-9 = 1.e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-8);
    refVal   = -12604.7443256780746682577305643;
    errBound = 12604.7443256780746682577305643 * std::max(eps/(1.-1.e-8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullNormalDoubles, reg, x = 10^-8 = 1.e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-7);
    refVal   = -8873.37197946933305448892692526;
    errBound = 8873.37197946933305448892692526 * std::max(eps/(1.-1.e-7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullNormalDoubles, reg, x = 10^-7 = 1.e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-6);
    refVal   = -5249.19302937878963370852838324;
    errBound = 5249.19302937878963370852838324 * std::max(eps/(1.-1.e-6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullNormalDoubles, reg, x = 10^-6 = 1.e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-5);
    refVal   = -2594.36769175909856271699359614;
    errBound = 2594.36769175909856271699359614 * std::max(eps/(1.-1.e-5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullNormalDoubles, reg, x = 10^-5 = 1.e-5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0001);
    refVal   = -1049.2201858008046775438042814;
    errBound = 1049.2201858008046775438042814 * std::max(eps/(1.-0.0001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullNormalDoubles, reg, x = 10^-4 = 0.0001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.001);
    refVal   = -340.696652105762594640057561391;
    errBound = 340.696652105762594640057561391 * std::max(eps/(1.-0.001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullNormalDoubles, reg, x = 10^-3 = 0.001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.01);
    refVal   = -85.1341267864517043984357032428;
    errBound = 85.1341267864517043984357032428 * std::max(eps/(1.-0.01),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullNormalDoubles, reg, x = 10^-2 = 0.01";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.05);
    refVal   = -37.5472968537352940165079656137;
    errBound = 37.5472968537352940165079656137 * std::max(eps/(1.-0.05),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullNormalDoubles, reg, x = 0.05 = 0.05";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.1);
    refVal   = -39.5424646153036818080807254033;
    errBound = 39.5424646153036818080807254033 * std::max(eps/(1.-0.1),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullNormalDoubles, reg, x = 0.1 = 0.1";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.15);
    refVal   = -44.6063825357561530654738361958;
    errBound = 44.6063825357561530654738361958 * std::max(eps/(1.-0.15),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullNormalDoubles, reg, x = 0.15 = 0.15";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.2);
    refVal   = -48.5611540959001906502158750065;
    errBound = 48.5611540959001906502158750065 * std::max(eps/(1.-0.2),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullNormalDoubles, reg, x = 0.2 = 0.2";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.25);
    refVal   = -50.9932421445315976662193458295;
    errBound = 50.9932421445315976662193458295 * std::max(eps/(1.-0.25),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullNormalDoubles, reg, x = 0.25 = 0.25";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.3);
    refVal   = -52.0217083410804945272387692177;
    errBound = 52.0217083410804945272387692177 * std::max(eps/(1.-0.3),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullNormalDoubles, reg, x = 0.3 = 0.3";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.35);
    refVal   = -51.8473048532975148594258894995;
    errBound = 51.8473048532975148594258894995 * std::max(eps/(1.-0.35),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullNormalDoubles, reg, x = 0.35 = 0.35";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.4);
    refVal   = -50.6623123181647806341144192365;
    errBound = 50.6623123181647806341144192365 * std::max(eps/(1.-0.4),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullNormalDoubles, reg, x = 0.4 = 0.4";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.45);
    refVal   = -48.6330256485727144574720898597;
    errBound = 48.6330256485727144574720898597 * std::max(eps/(1.-0.45),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullNormalDoubles, reg, x = 0.45 = 0.45";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.5);
    refVal   = -45.899522397045080095115919971;
    errBound = 45.899522397045080095115919971 * std::max(eps/(1.-0.5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullNormalDoubles, reg, x = 0.5 = 0.5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.55);
    refVal   = -42.5793896604565544599733832024;
    errBound = 42.5793896604565544599733832024 * std::max(eps/(1.-0.55),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullNormalDoubles, reg, x = 0.55 = 0.55";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.6);
    refVal   = -38.7718949300146020456578614548;
    errBound = 38.7718949300146020456578614548 * std::max(eps/(1.-0.6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullNormalDoubles, reg, x = 0.6 = 0.6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.65);
    refVal   = -34.5617387235110338681063230666;
    errBound = 34.5617387235110338681063230666 * std::max(eps/(1.-0.65),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullNormalDoubles, reg, x = 0.65 = 0.65";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.7);
    refVal   = -30.0223159109290793233209714975;
    errBound = 30.0223159109290793233209714975 * std::max(eps/(1.-0.7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullNormalDoubles, reg, x = 0.7 = 0.7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.75);
    refVal   = -25.218700173945311444131242133;
    errBound = 25.218700173945311444131242133 * std::max(eps/(1.-0.75),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullNormalDoubles, reg, x = 0.75 = 0.75";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.8);
    refVal   = -20.210812730146510510454642218;
    errBound = 20.210812730146510510454642218 * std::max(eps/(1.-0.8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullNormalDoubles, reg, x = 0.8 = 0.8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.85);
    refVal   = -15.0578449520193298772730152722;
    errBound = 15.0578449520193298772730152722 * std::max(eps/(1.-0.85),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullNormalDoubles, reg, x = 0.85 = 0.85";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9);
    refVal   = -9.82731790460264349484403427753;
    errBound = 9.82731790460264349484403427753 * std::max(eps/(1.-0.9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullNormalDoubles, reg, x = 0.9 = 0.9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.95);
    refVal   = -4.62533288002915454920725852522;
    errBound = 4.62533288002915454920725852522 * std::max(eps/(1.-0.95),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullNormalDoubles, reg, x = 0.95 = 0.95";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99);
    refVal   = -0.713646238014366027781186031711;
    errBound = 0.713646238014366027781186031711 * std::max(eps/(1.-0.99),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullNormalDoubles, reg, x = 1-10^-2 = 0.99";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999);
    refVal   = -0.0232283087372029495808061602361;
    errBound = 0.0232283087372029495808061602361 * std::max(eps/(1.-0.999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullNormalDoubles, reg, x = 1-10^-3 = 0.999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999);
    refVal   = 0.00466953401524600547159236972589;
    errBound = 0.00466953401524600547159236972589 * std::max(eps/(1.-0.9999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullNormalDoubles, reg, x = 1-10^-4 = 0.9999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999);
    refVal   = 0.00136032507454343475622471940041;
    errBound = 0.00136032507454343475622471940041 * std::max(eps/(1.-0.99999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullNormalDoubles, reg, x = 1-10^-5 = 0.99999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999);
    refVal   = 0.00023748368127672885836334414251;
    errBound = 0.00023748368127672885836334414251 * std::max(eps/(1.-0.999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullNormalDoubles, reg, x = 1-10^-6 = 0.999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999);
    refVal   = 0.0000338578103568962267428930344744;
    errBound = 0.0000338578103568962267428930344744 * std::max(eps/(1.-0.9999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullNormalDoubles, reg, x = 1-10^-7 = 0.9999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999999);
    refVal   = 4.21643402428458319064302343402e-6;
    errBound = 4.21643402428458319064302343402e-6 * std::max(eps/(1.-0.99999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullNormalDoubles, reg, x = 1-10^-8 = 0.99999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999);
    refVal   = 4.63802073845580104276018115926e-7;
    errBound = 4.63802073845580104276018115926e-7 * std::max(eps/(1.-0.999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullNormalDoubles, reg, x = 1-10^-9 = 0.999999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999999);
    refVal   = 4.3697108561469609789528489475e-8;
    errBound = 4.3697108561469609789528489475e-8 * std::max(eps/(1.-0.9999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, FullNormalDoubles, reg, x = 1-10^-10 = 0.9999999999";
  }
}



TEST(polAQqPSXspace,TruncAs2PowersOf4RegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto eval = polAQqPS_reg.truncate(2);

  {
    testVal  = eval(testAs,testLM,testNF,9.31322574615478515625e-10);
    refVal   = -59.3990812086629697196511738494;
    errBound = 59.3990812086629697196511738494 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2PowersOf4, reg, x = 4^-15 = 9.31322574615478515625e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.7252902984619140625e-9);
    refVal   = -69.8181238146933988877175854757;
    errBound = 69.8181238146933988877175854757 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2PowersOf4, reg, x = 4^-14 = 3.7252902984619140625e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.490116119384765625e-8);
    refVal   = -74.5305887651230232307099778816;
    errBound = 74.5305887651230232307099778816 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2PowersOf4, reg, x = 4^-13 = 1.490116119384765625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,5.9604644775390625e-8);
    refVal   = -74.4245433449897062191689439992;
    errBound = 74.4245433449897062191689439992 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2PowersOf4, reg, x = 4^-12 = 5.9604644775390625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,2.384185791015625e-7);
    refVal   = -70.3880590840178478437425290815;
    errBound = 70.3880590840178478437425290815 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2PowersOf4, reg, x = 4^-11 = 2.384185791015625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,9.5367431640625e-7);
    refVal   = -63.3092219249125391919162406075;
    errBound = 63.3092219249125391919162406075 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2PowersOf4, reg, x = 4^-10 = 9.5367431640625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.814697265625e-6);
    refVal   = -54.076162481772035924257230728;
    errBound = 54.076162481772035924257230728 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2PowersOf4, reg, x = 4^-9 = 3.814697265625e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0000152587890625);
    refVal   = -43.5771315498569752655623249415;
    errBound = 43.5771315498569752655623249415 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2PowersOf4, reg, x = 4^-8 = 0.0000152587890625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00006103515625);
    refVal   = -32.7006174528134434844385150827;
    errBound = 32.7006174528134434844385150827 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2PowersOf4, reg, x = 4^-7 = 0.00006103515625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.000244140625);
    refVal   = -22.3351177340981417646302483661;
    errBound = 22.3351177340981417646302483661 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2PowersOf4, reg, x = 4^-6 = 0.000244140625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0009765625);
    refVal   = -13.3658052881549061440826086477;
    errBound = 13.3658052881549061440826086477 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2PowersOf4, reg, x = 4^-5 = 0.0009765625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00390625);
    refVal   = -6.65381533603946500732591491233;
    errBound = 6.65381533603946500732591491233 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2PowersOf4, reg, x = 4^-4 = 0.00390625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.015625);
    refVal   = -2.93641412485751121217061259639;
    errBound = 2.93641412485751121217061259639 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2PowersOf4, reg, x = 4^-3 = 0.015625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0625);
    refVal   = -2.42333024996566615272937173056;
    errBound = 2.42333024996566615272937173056 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2PowersOf4, reg, x = 4^-2 = 0.0625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9375);
    refVal   = -0.397877362825530297555952532611;
    errBound = 0.397877362825530297555952532611 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2PowersOf4, reg, x = 1-4^-2 = 0.9375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.984375);
    refVal   = -0.0996676207259408022276561510504;
    errBound = 0.0996676207259408022276561510504 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2PowersOf4, reg, x = 1-4^-3 = 0.984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99609375);
    refVal   = -0.0249256807458307260515480786516;
    errBound = 0.0249256807458307260515480786516 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2PowersOf4, reg, x = 1-4^-4 = 0.99609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9990234375);
    refVal   = -0.00623191476581169835748846576146;
    errBound = 0.00623191476581169835748846576146 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2PowersOf4, reg, x = 1-4^-5 = 0.9990234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999755859375);
    refVal   = -0.00155800877035294503499974172078;
    errBound = 0.00155800877035294503499974172078 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2PowersOf4, reg, x = 1-4^-6 = 0.999755859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99993896484375);
    refVal   = -0.000389504059552401467117168397412;
    errBound = 0.000389504059552401467117168397412 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2PowersOf4, reg, x = 1-4^-7 = 0.99993896484375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999847412109375);
    refVal   = -0.0000973761313708970220963451465086;
    errBound = 0.0000973761313708970220963451465086 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2PowersOf4, reg, x = 1-4^-8 = 0.9999847412109375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999996185302734375);
    refVal   = -0.0000243440401197361256407695508288;
    errBound = 0.0000243440401197361256407695508288 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2PowersOf4, reg, x = 1-4^-9 = 0.999996185302734375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999904632568359375);
    refVal   = -6.08601048469785489616597970888e-6;
    errBound = 6.08601048469785489616597970888e-6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2PowersOf4, reg, x = 1-4^-10 = 0.99999904632568359375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999997615814208984375);
    refVal   = -1.5215026495964305382290640492e-6;
    errBound = 1.5215026495964305382290640492e-6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2PowersOf4, reg, x = 1-4^-11 = 0.9999997615814208984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999940395355224609375);
    refVal   = -3.8037566417546849557685426782e-7;
    errBound = 3.8037566417546849557685426782e-7 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2PowersOf4, reg, x = 1-4^-12 = 0.999999940395355224609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999998509883880615234375);
    refVal   = -9.50939161548894891945224663268e-8;
    errBound = 9.50939161548894891945224663268e-8 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2PowersOf4, reg, x = 1-4^-13 = 0.99999998509883880615234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999962747097015380859375);
    refVal   = -2.37734790456612671843783296221e-8;
    errBound = 2.37734790456612671843783296221e-8 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2PowersOf4, reg, x = 1-4^-14 = 0.9999999962747097015380859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999068677425384521484375);
    refVal   = -5.94336976184899768043004161607e-9;
    errBound = 5.94336976184899768043004161607e-9 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2PowersOf4, reg, x = 1-4^-15 = 0.999999999068677425384521484375";
  }
}



TEST(polAQqPSXspace,TruncAs2NormalDoublesRegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto eval = polAQqPS_reg.truncate(2);

  {
    testVal  = eval(testAs,testLM,testNF,1.e-10);
    refVal   = -28.3984612644582784449692932653;
    errBound = 28.3984612644582784449692932653 * std::max(eps/(1.-1.e-10),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2NormalDoubles, reg, x = 10^-10 = 1.e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-9);
    refVal   = -60.0867920948464965970917321012;
    errBound = 60.0867920948464965970917321012 * std::max(eps/(1.-1.e-9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2NormalDoubles, reg, x = 10^-9 = 1.e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-8);
    refVal   = -73.7075520885662476831273100175;
    errBound = 73.7075520885662476831273100175 * std::max(eps/(1.-1.e-8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2NormalDoubles, reg, x = 10^-8 = 1.e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-7);
    refVal   = -73.3301028712697833375741070199;
    errBound = 73.3301028712697833375741070199 * std::max(eps/(1.-1.e-7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2NormalDoubles, reg, x = 10^-7 = 1.e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-6);
    refVal   = -63.0238403617535335329540080481;
    errBound = 63.0238403617535335329540080481 * std::max(eps/(1.-1.e-6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2NormalDoubles, reg, x = 10^-6 = 1.e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-5);
    refVal   = -46.8583937278439029257594364377;
    errBound = 46.8583937278439029257594364377 * std::max(eps/(1.-1.e-5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2NormalDoubles, reg, x = 10^-5 = 1.e-5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0001);
    refVal   = -28.9043971116741630855670763838;
    errBound = 28.9043971116741630855670763838 * std::max(eps/(1.-0.0001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2NormalDoubles, reg, x = 10^-4 = 0.0001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.001);
    refVal   = -13.2296223998274471171222666628;
    errBound = 13.2296223998274471171222666628 * std::max(eps/(1.-0.001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2NormalDoubles, reg, x = 10^-3 = 0.001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.01);
    refVal   = -3.77444693155044997366315088525;
    errBound = 3.77444693155044997366315088525 * std::max(eps/(1.-0.01),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2NormalDoubles, reg, x = 10^-2 = 0.01";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.05);
    refVal   = -2.31924069708402567642905047546;
    errBound = 2.31924069708402567642905047546 * std::max(eps/(1.-0.05),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2NormalDoubles, reg, x = 0.05 = 0.05";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.1);
    refVal   = -2.78970249994209802590754441145;
    errBound = 2.78970249994209802590754441145 * std::max(eps/(1.-0.1),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2NormalDoubles, reg, x = 0.1 = 0.1";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.15);
    refVal   = -3.16500946613673234317587876731;
    errBound = 3.16500946613673234317587876731 * std::max(eps/(1.-0.15),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2NormalDoubles, reg, x = 0.15 = 0.15";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.2);
    refVal   = -3.38264565633502664500340957029;
    errBound = 3.38264565633502664500340957029 * std::max(eps/(1.-0.2),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2NormalDoubles, reg, x = 0.2 = 0.2";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.25);
    refVal   = -3.4756552202451714115802930902;
    errBound = 3.4756552202451714115802930902 * std::max(eps/(1.-0.25),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2NormalDoubles, reg, x = 0.25 = 0.25";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.3);
    refVal   = -3.47443700279242618189019750141;
    errBound = 3.47443700279242618189019750141 * std::max(eps/(1.-0.3),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2NormalDoubles, reg, x = 0.3 = 0.3";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.35);
    refVal   = -3.40150106736440269790946329967;
    errBound = 3.40150106736440269790946329967 * std::max(eps/(1.-0.35),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2NormalDoubles, reg, x = 0.35 = 0.35";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.4);
    refVal   = -3.27325680027028688031620218518;
    errBound = 3.27325680027028688031620218518 * std::max(eps/(1.-0.4),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2NormalDoubles, reg, x = 0.4 = 0.4";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.45);
    refVal   = -3.10184673511955548308572892198;
    errBound = 3.10184673511955548308572892198 * std::max(eps/(1.-0.45),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2NormalDoubles, reg, x = 0.45 = 0.45";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.5);
    refVal   = -2.8964391762470998897499572542;
    errBound = 2.8964391762470998897499572542 * std::max(eps/(1.-0.5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2NormalDoubles, reg, x = 0.5 = 0.5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.55);
    refVal   = -2.66409373695044917262350635363;
    errBound = 2.66409373695044917262350635363 * std::max(eps/(1.-0.55),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2NormalDoubles, reg, x = 0.55 = 0.55";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.6);
    refVal   = -2.41034311116220807233409255793;
    errBound = 2.41034311116220807233409255793 * std::max(eps/(1.-0.6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2NormalDoubles, reg, x = 0.6 = 0.6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.65);
    refVal   = -2.13959115111208454364150837002;
    errBound = 2.13959115111208454364150837002 * std::max(eps/(1.-0.65),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2NormalDoubles, reg, x = 0.65 = 0.65";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.7);
    refVal   = -1.85539104897020818448578588742;
    errBound = 1.85539104897020818448578588742 * std::max(eps/(1.-0.7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2NormalDoubles, reg, x = 0.7 = 0.7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.75);
    refVal   = -1.56064382621114576767932623466;
    errBound = 1.56064382621114576767932623466 * std::max(eps/(1.-0.75),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2NormalDoubles, reg, x = 0.75 = 0.75";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.8);
    refVal   = -1.25774274698059216987910727669;
    errBound = 1.25774274698059216987910727669 * std::max(eps/(1.-0.8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2NormalDoubles, reg, x = 0.8 = 0.8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.85);
    refVal   = -0.948680277825833730919262117472;
    errBound = 0.948680277825833730919262117472 * std::max(eps/(1.-0.85),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2NormalDoubles, reg, x = 0.85 = 0.85";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9);
    refVal   = -0.635128599494371993226121409051;
    errBound = 0.635128599494371993226121409051 * std::max(eps/(1.-0.9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2NormalDoubles, reg, x = 0.9 = 0.9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.95);
    refVal   = -0.31850110475889604151463097659;
    errBound = 0.31850110475889604151463097659 * std::max(eps/(1.-0.95),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2NormalDoubles, reg, x = 0.95 = 0.95";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99);
    refVal   = -0.0637985044898329143144934370971;
    errBound = 0.0637985044898329143144934370971 * std::max(eps/(1.-0.99),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2NormalDoubles, reg, x = 1-10^-2 = 0.99";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999);
    refVal   = -0.00638147675522234277627290380143;
    errBound = 0.00638147675522234277627290380143 * std::max(eps/(1.-0.999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2NormalDoubles, reg, x = 1-10^-3 = 0.999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999);
    refVal   = -0.000638162800962074830532715846037;
    errBound = 0.000638162800962074830532715846037 * std::max(eps/(1.-0.9999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2NormalDoubles, reg, x = 1-10^-4 = 0.9999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999);
    refVal   = -0.0000638164302215613081104471864501;
    errBound = 0.0000638164302215613081104471864501 * std::max(eps/(1.-0.99999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2NormalDoubles, reg, x = 1-10^-5 = 0.99999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999);
    refVal   = -6.38164452228147604074968040199e-6;
    errBound = 6.38164452228147604074968040199e-6 * std::max(eps/(1.-0.999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2NormalDoubles, reg, x = 1-10^-6 = 0.999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999);
    refVal   = -6.38164467228272948445199708329e-7;
    errBound = 6.38164467228272948445199708329e-7 * std::max(eps/(1.-0.9999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2NormalDoubles, reg, x = 1-10^-7 = 0.9999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999999);
    refVal   = -6.38164468728274201888042557909e-8;
    errBound = 6.38164468728274201888042557909e-8 * std::max(eps/(1.-0.99999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2NormalDoubles, reg, x = 1-10^-8 = 0.99999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999);
    refVal   = -6.38164468878274214422470126939e-9;
    errBound = 6.38164468878274214422470126939e-9 * std::max(eps/(1.-0.999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2NormalDoubles, reg, x = 1-10^-9 = 0.999999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999999);
    refVal   = -6.3816446889327421454781440177e-10;
    errBound = 6.3816446889327421454781440177e-10 * std::max(eps/(1.-0.9999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPS, TruncAs2NormalDoubles, reg, x = 1-10^-10 = 0.9999999999";
  }
}



TEST(polAQqPSNspace,FullMoments)
{
  double refVal = 0., errBound = 0.;
  auto& rpd = polAQqPS;
  integration_engine_gsl engine;
  auto mom = make_mellin_moment(rpd, engine, testAs, testLM, testNF);
  std::tuple<integration_status, double, double> res;

  {
    SCOPED_TRACE("polAQqPS Full Mellin moment N=3");
    refVal = -7.68267731234743584216572285594;
    errBound = 7.68267731234743584216572285594 * 3.e-11;
    res = mom.integrate(3, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAQqPS Full Mellin moment N=5");
    refVal = -3.22885196133530965051579013891;
    errBound = 3.22885196133530965051579013891 * 3.e-11;
    res = mom.integrate(5, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAQqPS Full Mellin moment N=7");
    refVal = -1.74878750092145554306292312844;
    errBound = 1.74878750092145554306292312844 * 3.e-11;
    res = mom.integrate(7, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

}


TEST(polAQqPSNspace,TruncAs2Moments)
{
  double refVal = 0., errBound = 0.;
  auto rpd = polAQqPS.truncate(2);
  integration_engine_gsl engine;
  auto mom = make_mellin_moment(rpd, engine, testAs, testLM, testNF);
  std::tuple<integration_status, double, double> res;

  {
    SCOPED_TRACE("polAQqPS TruncAs2 Mellin moment N=3");
    refVal = -0.487316743827160493827160493827;
    errBound = 0.487316743827160493827160493827 * 3.e-11;
    res = mom.integrate(3, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAQqPS TruncAs2 Mellin moment N=5");
    refVal = -0.204126337448559670781893004115;
    errBound = 0.204126337448559670781893004115 * 3.e-11;
    res = mom.integrate(5, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAQqPS TruncAs2 Mellin moment N=7");
    refVal = -0.111178948062005414410662224073;
    errBound = 0.111178948062005414410662224073 * 3.e-11;
    res = mom.integrate(7, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

}
