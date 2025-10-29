/******************************************************************************
 * This file is part of libome                                                *
 * Copyright (C) 2025 Arnd Behring, Kay Schoenwald                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 ******************************************************************************/

#include <limits>
#include <gtest/gtest.h>
#include <ome/integration_engine_gsl.h>
#include <ome/mellin.h>
#include <ome/polAQg.h>

using namespace ome;

const double eps = std::numeric_limits<double>::epsilon();
const double testAs = 0.25;
const double testLM = -5.;
const double testNF = 3.;



TEST(polAQgXspace,FullPowersOf4RegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto& eval = polAQg_reg;

  {
    testVal  = eval(testAs,testLM,testNF,9.31322574615478515625e-10);
    refVal   = -40651.0480198341236973719714047;
    errBound = 40651.0480198341236973719714047 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullPowersOf4, reg, x = 4^-15 = 9.31322574615478515625e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.7252902984619140625e-9);
    refVal   = -41471.8204206167248654495347957;
    errBound = 41471.8204206167248654495347957 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullPowersOf4, reg, x = 4^-14 = 3.7252902984619140625e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.490116119384765625e-8);
    refVal   = -37600.7869156532349703280419558;
    errBound = 37600.7869156532349703280419558 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullPowersOf4, reg, x = 4^-13 = 1.490116119384765625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,5.9604644775390625e-8);
    refVal   = -31293.4282619730105093551797204;
    errBound = 31293.4282619730105093551797204 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullPowersOf4, reg, x = 4^-12 = 5.9604644775390625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,2.384185791015625e-7);
    refVal   = -24190.6539790880121121693432289;
    errBound = 24190.6539790880121121693432289 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullPowersOf4, reg, x = 4^-11 = 2.384185791015625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,9.5367431640625e-7);
    refVal   = -17403.6860980893429239287162684;
    errBound = 17403.6860980893429239287162684 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullPowersOf4, reg, x = 4^-10 = 9.5367431640625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.814697265625e-6);
    refVal   = -11598.9405255505529380667091439;
    errBound = 11598.9405255505529380667091439 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullPowersOf4, reg, x = 4^-9 = 3.814697265625e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0000152587890625);
    refVal   = -7082.8994955168171096242374714;
    errBound = 7082.8994955168171096242374714 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullPowersOf4, reg, x = 4^-8 = 0.0000152587890625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00006103515625);
    refVal   = -3886.95073088169778503145260141;
    errBound = 3886.95073088169778503145260141 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullPowersOf4, reg, x = 4^-7 = 0.00006103515625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.000244140625);
    refVal   = -1852.09612412182715943170971072;
    errBound = 1852.09612412182715943170971072 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullPowersOf4, reg, x = 4^-6 = 0.000244140625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0009765625);
    refVal   = -713.123159974766489161034777255;
    errBound = 713.123159974766489161034777255 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullPowersOf4, reg, x = 4^-5 = 0.0009765625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00390625);
    refVal   = -180.533813509245094970604962079;
    errBound = 180.533813509245094970604962079 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullPowersOf4, reg, x = 4^-4 = 0.00390625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.015625);
    refVal   = -13.4207613794804493630949646136;
    errBound = 13.4207613794804493630949646136 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullPowersOf4, reg, x = 4^-3 = 0.015625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0625);
    refVal   = -58.3680183989056710610646077198;
    errBound = 58.3680183989056710610646077198 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullPowersOf4, reg, x = 4^-2 = 0.0625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9375);
    refVal   = 379.725627427044366602252604329;
    errBound = 379.725627427044366602252604329 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullPowersOf4, reg, x = 1-4^-2 = 0.9375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.984375);
    refVal   = 460.55490862080132843760770475;
    errBound = 460.55490862080132843760770475 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullPowersOf4, reg, x = 1-4^-3 = 0.984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99609375);
    refVal   = 416.938369865732367015553789973;
    errBound = 416.938369865732367015553789973 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullPowersOf4, reg, x = 1-4^-4 = 0.99609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9990234375);
    refVal   = 307.157639746335080222046357354;
    errBound = 307.157639746335080222046357354 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullPowersOf4, reg, x = 1-4^-5 = 0.9990234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999755859375);
    refVal   = 171.964306111688153407016699828;
    errBound = 171.964306111688153407016699828 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullPowersOf4, reg, x = 1-4^-6 = 0.999755859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99993896484375);
    refVal   = 29.0106836770145756359667744066;
    errBound = 29.0106836770145756359667744066 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullPowersOf4, reg, x = 1-4^-7 = 0.99993896484375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999847412109375);
    refVal   = -140.30018795667430366415419705;
    errBound = 140.30018795667430366415419705 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullPowersOf4, reg, x = 1-4^-8 = 0.9999847412109375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999996185302734375);
    refVal   = -407.281094377249157794858419049;
    errBound = 407.281094377249157794858419049 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullPowersOf4, reg, x = 1-4^-9 = 0.999996185302734375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999904632568359375);
    refVal   = -913.428261759586287647281260091;
    errBound = 913.428261759586287647281260091 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullPowersOf4, reg, x = 1-4^-10 = 0.99999904632568359375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999997615814208984375);
    refVal   = -1888.13663000932301864049042185;
    errBound = 1888.13663000932301864049042185 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullPowersOf4, reg, x = 1-4^-11 = 0.9999997615814208984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999940395355224609375);
    refVal   = -3666.46757177140560627865095349;
    errBound = 3666.46757177140560627865095349 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullPowersOf4, reg, x = 1-4^-12 = 0.999999940395355224609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999998509883880615234375);
    refVal   = -6706.92598807459477001812417781;
    errBound = 6706.92598807459477001812417781 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullPowersOf4, reg, x = 1-4^-13 = 0.99999998509883880615234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999962747097015380859375);
    refVal   = -11609.2386340612418901950914271;
    errBound = 11609.2386340612418901950914271 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullPowersOf4, reg, x = 1-4^-14 = 0.9999999962747097015380859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999068677425384521484375);
    refVal   = -19132.1324333927122887534808086;
    errBound = 19132.1324333927122887534808086 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullPowersOf4, reg, x = 1-4^-15 = 0.999999999068677425384521484375";
  }
}



TEST(polAQgXspace,FullNormalDoublesRegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto& eval = polAQg_reg;

  {
    testVal  = eval(testAs,testLM,testNF,1.e-10);
    refVal   = -21708.7825510519154726743582026;
    errBound = 21708.7825510519154726743582026 * std::max(eps/(1.-1.e-10),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullNormalDoubles, reg, x = 10^-10 = 1.e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-9);
    refVal   = -40851.0059887552659186014263776;
    errBound = 40851.0059887552659186014263776 * std::max(eps/(1.-1.e-9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullNormalDoubles, reg, x = 10^-9 = 1.e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-8);
    refVal   = -39051.0739001223026782714104383;
    errBound = 39051.0739001223026782714104383 * std::max(eps/(1.-1.e-8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullNormalDoubles, reg, x = 10^-8 = 1.e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-7);
    refVal   = -28665.1901524030737694709637648;
    errBound = 28665.1901524030737694709637648 * std::max(eps/(1.-1.e-7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullNormalDoubles, reg, x = 10^-7 = 1.e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-6);
    refVal   = -17186.1872066508824547300457091;
    errBound = 17186.1872066508824547300457091 * std::max(eps/(1.-1.e-6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullNormalDoubles, reg, x = 10^-6 = 1.e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-5);
    refVal   = -8316.57979705564308746077942461;
    errBound = 8316.57979705564308746077942461 * std::max(eps/(1.-1.e-5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullNormalDoubles, reg, x = 10^-5 = 1.e-5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0001);
    refVal   = -3041.44924172939047206329879733;
    errBound = 3041.44924172939047206329879733 * std::max(eps/(1.-0.0001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullNormalDoubles, reg, x = 10^-4 = 0.0001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.001);
    refVal   = -699.705578367939114569993801423;
    errBound = 699.705578367939114569993801423 * std::max(eps/(1.-0.001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullNormalDoubles, reg, x = 10^-3 = 0.001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.01);
    refVal   = -38.4229269205927832904105918395;
    errBound = 38.4229269205927832904105918395 * std::max(eps/(1.-0.01),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullNormalDoubles, reg, x = 10^-2 = 0.01";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.05);
    refVal   = -41.3236367219493127452868411549;
    errBound = 41.3236367219493127452868411549 * std::max(eps/(1.-0.05),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullNormalDoubles, reg, x = 0.05 = 0.05";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.1);
    refVal   = -101.158051805773854972284151978;
    errBound = 101.158051805773854972284151978 * std::max(eps/(1.-0.1),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullNormalDoubles, reg, x = 0.1 = 0.1";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.15);
    refVal   = -139.160104604561259132854292399;
    errBound = 139.160104604561259132854292399 * std::max(eps/(1.-0.15),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullNormalDoubles, reg, x = 0.15 = 0.15";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.2);
    refVal   = -160.703960183129106485406370383;
    errBound = 160.703960183129106485406370383 * std::max(eps/(1.-0.2),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullNormalDoubles, reg, x = 0.2 = 0.2";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.25);
    refVal   = -170.056596512741042629185082107;
    errBound = 170.056596512741042629185082107 * std::max(eps/(1.-0.25),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullNormalDoubles, reg, x = 0.25 = 0.25";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.3);
    refVal   = -169.911194307055003456016573361;
    errBound = 169.911194307055003456016573361 * std::max(eps/(1.-0.3),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullNormalDoubles, reg, x = 0.3 = 0.3";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.35);
    refVal   = -162.015778572384455974001065923;
    errBound = 162.015778572384455974001065923 * std::max(eps/(1.-0.35),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullNormalDoubles, reg, x = 0.35 = 0.35";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.4);
    refVal   = -147.53422673010897871899639903;
    errBound = 147.53422673010897871899639903 * std::max(eps/(1.-0.4),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullNormalDoubles, reg, x = 0.4 = 0.4";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.45);
    refVal   = -127.24096303566441810720119726;
    errBound = 127.24096303566441810720119726 * std::max(eps/(1.-0.45),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullNormalDoubles, reg, x = 0.45 = 0.45";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.5);
    refVal   = -101.627687211753577043583465207;
    errBound = 101.627687211753577043583465207 * std::max(eps/(1.-0.5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullNormalDoubles, reg, x = 0.5 = 0.5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.55);
    refVal   = -70.9605394319641577816455634374;
    errBound = 70.9605394319641577816455634374 * std::max(eps/(1.-0.55),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullNormalDoubles, reg, x = 0.55 = 0.55";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.6);
    refVal   = -35.3055379051012492068668894908;
    errBound = 35.3055379051012492068668894908 * std::max(eps/(1.-0.6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullNormalDoubles, reg, x = 0.6 = 0.6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.65);
    refVal   = 5.47071505920533851962852550195;
    errBound = 5.47071505920533851962852550195 * std::max(eps/(1.-0.65),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullNormalDoubles, reg, x = 0.65 = 0.65";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.7);
    refVal   = 51.7252917307341931304116130636;
    errBound = 51.7252917307341931304116130636 * std::max(eps/(1.-0.7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullNormalDoubles, reg, x = 0.7 = 0.7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.75);
    refVal   = 104.095587910077569648000070387;
    errBound = 104.095587910077569648000070387 * std::max(eps/(1.-0.75),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullNormalDoubles, reg, x = 0.75 = 0.75";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.8);
    refVal   = 163.60570694385662057965999417;
    errBound = 163.60570694385662057965999417 * std::max(eps/(1.-0.8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullNormalDoubles, reg, x = 0.8 = 0.8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.85);
    refVal   = 231.851370054606330379050875246;
    errBound = 231.851370054606330379050875246 * std::max(eps/(1.-0.85),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullNormalDoubles, reg, x = 0.85 = 0.85";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9);
    refVal   = 311.245129901422356904646887611;
    errBound = 311.245129901422356904646887611 * std::max(eps/(1.-0.9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullNormalDoubles, reg, x = 0.9 = 0.9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.95);
    refVal   = 404.083523760977497379235968401;
    errBound = 404.083523760977497379235968401 * std::max(eps/(1.-0.95),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullNormalDoubles, reg, x = 0.95 = 0.95";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99);
    refVal   = 456.93443873983794977652333914;
    errBound = 456.93443873983794977652333914 * std::max(eps/(1.-0.99),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullNormalDoubles, reg, x = 1-10^-2 = 0.99";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999);
    refVal   = 309.338836551054913353309117766;
    errBound = 309.338836551054913353309117766 * std::max(eps/(1.-0.999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullNormalDoubles, reg, x = 1-10^-3 = 0.999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999);
    refVal   = 81.0304302457009211747902191841;
    errBound = 81.0304302457009211747902191841 * std::max(eps/(1.-0.9999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullNormalDoubles, reg, x = 1-10^-4 = 0.9999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999);
    refVal   = -206.342177032883330361094398683;
    errBound = 206.342177032883330361094398683 * std::max(eps/(1.-0.99999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullNormalDoubles, reg, x = 1-10^-5 = 0.99999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999);
    refVal   = -889.950435241849837587483415461;
    errBound = 889.950435241849837587483415461 * std::max(eps/(1.-0.999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullNormalDoubles, reg, x = 1-10^-6 = 0.999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999);
    refVal   = -2882.51533174194653897442127092;
    errBound = 2882.51533174194653897442127092 * std::max(eps/(1.-0.9999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullNormalDoubles, reg, x = 1-10^-7 = 0.9999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999999);
    refVal   = -7896.37936992379969610309800945;
    errBound = 7896.37936992379969610309800945 * std::max(eps/(1.-0.99999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullNormalDoubles, reg, x = 1-10^-8 = 0.99999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999);
    refVal   = -18668.2181475539241458037250125;
    errBound = 18668.2181475539241458037250125 * std::max(eps/(1.-0.999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullNormalDoubles, reg, x = 1-10^-9 = 0.999999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999999);
    refVal   = -39183.7854646474923900045776552;
    errBound = 39183.7854646474923900045776552 * std::max(eps/(1.-0.9999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, FullNormalDoubles, reg, x = 1-10^-10 = 0.9999999999";
  }
}



TEST(polAQgXspace,TruncAs1PowersOf4RegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto eval = polAQg_reg.truncate(1);

  {
    testVal  = eval(testAs,testLM,testNF,9.31322574615478515625e-10);
    refVal   = -2.49999999534338712692260742188;
    errBound = 2.49999999534338712692260742188 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1PowersOf4, reg, x = 4^-15 = 9.31322574615478515625e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.7252902984619140625e-9);
    refVal   = -2.4999999813735485076904296875;
    errBound = 2.4999999813735485076904296875 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1PowersOf4, reg, x = 4^-14 = 3.7252902984619140625e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.490116119384765625e-8);
    refVal   = -2.49999992549419403076171875;
    errBound = 2.49999992549419403076171875 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1PowersOf4, reg, x = 4^-13 = 1.490116119384765625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,5.9604644775390625e-8);
    refVal   = -2.499999701976776123046875;
    errBound = 2.499999701976776123046875 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1PowersOf4, reg, x = 4^-12 = 5.9604644775390625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,2.384185791015625e-7);
    refVal   = -2.4999988079071044921875;
    errBound = 2.4999988079071044921875 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1PowersOf4, reg, x = 4^-11 = 2.384185791015625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,9.5367431640625e-7);
    refVal   = -2.49999523162841796875;
    errBound = 2.49999523162841796875 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1PowersOf4, reg, x = 4^-10 = 9.5367431640625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.814697265625e-6);
    refVal   = -2.499980926513671875;
    errBound = 2.499980926513671875 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1PowersOf4, reg, x = 4^-9 = 3.814697265625e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0000152587890625);
    refVal   = -2.4999237060546875;
    errBound = 2.4999237060546875 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1PowersOf4, reg, x = 4^-8 = 0.0000152587890625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00006103515625);
    refVal   = -2.49969482421875;
    errBound = 2.49969482421875 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1PowersOf4, reg, x = 4^-7 = 0.00006103515625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.000244140625);
    refVal   = -2.498779296875;
    errBound = 2.498779296875 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1PowersOf4, reg, x = 4^-6 = 0.000244140625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0009765625);
    refVal   = -2.4951171875;
    errBound = 2.4951171875 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1PowersOf4, reg, x = 4^-5 = 0.0009765625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00390625);
    refVal   = -2.48046875;
    errBound = 2.48046875 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1PowersOf4, reg, x = 4^-4 = 0.00390625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.015625);
    refVal   = -2.421875;
    errBound = 2.421875 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1PowersOf4, reg, x = 4^-3 = 0.015625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0625);
    refVal   = -2.1875;
    errBound = 2.1875 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1PowersOf4, reg, x = 4^-2 = 0.0625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9375);
    refVal   = 2.1875;
    errBound = 2.1875 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1PowersOf4, reg, x = 1-4^-2 = 0.9375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.984375);
    refVal   = 2.421875;
    errBound = 2.421875 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1PowersOf4, reg, x = 1-4^-3 = 0.984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99609375);
    refVal   = 2.48046875;
    errBound = 2.48046875 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1PowersOf4, reg, x = 1-4^-4 = 0.99609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9990234375);
    refVal   = 2.4951171875;
    errBound = 2.4951171875 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1PowersOf4, reg, x = 1-4^-5 = 0.9990234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999755859375);
    refVal   = 2.498779296875;
    errBound = 2.498779296875 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1PowersOf4, reg, x = 1-4^-6 = 0.999755859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99993896484375);
    refVal   = 2.49969482421875;
    errBound = 2.49969482421875 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1PowersOf4, reg, x = 1-4^-7 = 0.99993896484375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999847412109375);
    refVal   = 2.4999237060546875;
    errBound = 2.4999237060546875 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1PowersOf4, reg, x = 1-4^-8 = 0.9999847412109375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999996185302734375);
    refVal   = 2.499980926513671875;
    errBound = 2.499980926513671875 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1PowersOf4, reg, x = 1-4^-9 = 0.999996185302734375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999904632568359375);
    refVal   = 2.49999523162841796875;
    errBound = 2.49999523162841796875 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1PowersOf4, reg, x = 1-4^-10 = 0.99999904632568359375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999997615814208984375);
    refVal   = 2.4999988079071044921875;
    errBound = 2.4999988079071044921875 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1PowersOf4, reg, x = 1-4^-11 = 0.9999997615814208984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999940395355224609375);
    refVal   = 2.499999701976776123046875;
    errBound = 2.499999701976776123046875 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1PowersOf4, reg, x = 1-4^-12 = 0.999999940395355224609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999998509883880615234375);
    refVal   = 2.49999992549419403076171875;
    errBound = 2.49999992549419403076171875 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1PowersOf4, reg, x = 1-4^-13 = 0.99999998509883880615234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999962747097015380859375);
    refVal   = 2.4999999813735485076904296875;
    errBound = 2.4999999813735485076904296875 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1PowersOf4, reg, x = 1-4^-14 = 0.9999999962747097015380859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999068677425384521484375);
    refVal   = 2.49999999534338712692260742188;
    errBound = 2.49999999534338712692260742188 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1PowersOf4, reg, x = 1-4^-15 = 0.999999999068677425384521484375";
  }
}



TEST(polAQgXspace,TruncAs1NormalDoublesRegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto eval = polAQg_reg.truncate(1);

  {
    testVal  = eval(testAs,testLM,testNF,1.e-10);
    refVal   = -2.4999999995;
    errBound = 2.4999999995 * std::max(eps/(1.-1.e-10),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1NormalDoubles, reg, x = 10^-10 = 1.e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-9);
    refVal   = -2.499999995;
    errBound = 2.499999995 * std::max(eps/(1.-1.e-9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1NormalDoubles, reg, x = 10^-9 = 1.e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-8);
    refVal   = -2.49999995;
    errBound = 2.49999995 * std::max(eps/(1.-1.e-8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1NormalDoubles, reg, x = 10^-8 = 1.e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-7);
    refVal   = -2.4999995;
    errBound = 2.4999995 * std::max(eps/(1.-1.e-7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1NormalDoubles, reg, x = 10^-7 = 1.e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-6);
    refVal   = -2.499995;
    errBound = 2.499995 * std::max(eps/(1.-1.e-6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1NormalDoubles, reg, x = 10^-6 = 1.e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-5);
    refVal   = -2.49995;
    errBound = 2.49995 * std::max(eps/(1.-1.e-5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1NormalDoubles, reg, x = 10^-5 = 1.e-5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0001);
    refVal   = -2.4995;
    errBound = 2.4995 * std::max(eps/(1.-0.0001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1NormalDoubles, reg, x = 10^-4 = 0.0001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.001);
    refVal   = -2.495;
    errBound = 2.495 * std::max(eps/(1.-0.001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1NormalDoubles, reg, x = 10^-3 = 0.001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.01);
    refVal   = -2.45;
    errBound = 2.45 * std::max(eps/(1.-0.01),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1NormalDoubles, reg, x = 10^-2 = 0.01";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.05);
    refVal   = -2.25;
    errBound = 2.25 * std::max(eps/(1.-0.05),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1NormalDoubles, reg, x = 0.05 = 0.05";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.1);
    refVal   = -2.;
    errBound = 2. * std::max(eps/(1.-0.1),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1NormalDoubles, reg, x = 0.1 = 0.1";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.15);
    refVal   = -1.75;
    errBound = 1.75 * std::max(eps/(1.-0.15),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1NormalDoubles, reg, x = 0.15 = 0.15";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.2);
    refVal   = -1.5;
    errBound = 1.5 * std::max(eps/(1.-0.2),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1NormalDoubles, reg, x = 0.2 = 0.2";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.25);
    refVal   = -1.25;
    errBound = 1.25 * std::max(eps/(1.-0.25),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1NormalDoubles, reg, x = 0.25 = 0.25";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.3);
    refVal   = -1.;
    errBound = 1. * std::max(eps/(1.-0.3),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1NormalDoubles, reg, x = 0.3 = 0.3";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.35);
    refVal   = -0.75;
    errBound = 0.75 * std::max(eps/(1.-0.35),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1NormalDoubles, reg, x = 0.35 = 0.35";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.4);
    refVal   = -0.5;
    errBound = 0.5 * std::max(eps/(1.-0.4),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1NormalDoubles, reg, x = 0.4 = 0.4";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.45);
    refVal   = -0.25;
    errBound = 0.25 * std::max(eps/(1.-0.45),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1NormalDoubles, reg, x = 0.45 = 0.45";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.5);
    refVal   = 0;
    errBound = 0 * std::max(eps/(1.-0.5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1NormalDoubles, reg, x = 0.5 = 0.5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.55);
    refVal   = 0.25;
    errBound = 0.25 * std::max(eps/(1.-0.55),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1NormalDoubles, reg, x = 0.55 = 0.55";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.6);
    refVal   = 0.5;
    errBound = 0.5 * std::max(eps/(1.-0.6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1NormalDoubles, reg, x = 0.6 = 0.6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.65);
    refVal   = 0.75;
    errBound = 0.75 * std::max(eps/(1.-0.65),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1NormalDoubles, reg, x = 0.65 = 0.65";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.7);
    refVal   = 1.;
    errBound = 1. * std::max(eps/(1.-0.7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1NormalDoubles, reg, x = 0.7 = 0.7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.75);
    refVal   = 1.25;
    errBound = 1.25 * std::max(eps/(1.-0.75),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1NormalDoubles, reg, x = 0.75 = 0.75";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.8);
    refVal   = 1.5;
    errBound = 1.5 * std::max(eps/(1.-0.8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1NormalDoubles, reg, x = 0.8 = 0.8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.85);
    refVal   = 1.75;
    errBound = 1.75 * std::max(eps/(1.-0.85),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1NormalDoubles, reg, x = 0.85 = 0.85";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9);
    refVal   = 2.;
    errBound = 2. * std::max(eps/(1.-0.9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1NormalDoubles, reg, x = 0.9 = 0.9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.95);
    refVal   = 2.25;
    errBound = 2.25 * std::max(eps/(1.-0.95),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1NormalDoubles, reg, x = 0.95 = 0.95";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99);
    refVal   = 2.45;
    errBound = 2.45 * std::max(eps/(1.-0.99),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1NormalDoubles, reg, x = 1-10^-2 = 0.99";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999);
    refVal   = 2.495;
    errBound = 2.495 * std::max(eps/(1.-0.999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1NormalDoubles, reg, x = 1-10^-3 = 0.999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999);
    refVal   = 2.4995;
    errBound = 2.4995 * std::max(eps/(1.-0.9999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1NormalDoubles, reg, x = 1-10^-4 = 0.9999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999);
    refVal   = 2.49995;
    errBound = 2.49995 * std::max(eps/(1.-0.99999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1NormalDoubles, reg, x = 1-10^-5 = 0.99999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999);
    refVal   = 2.499995;
    errBound = 2.499995 * std::max(eps/(1.-0.999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1NormalDoubles, reg, x = 1-10^-6 = 0.999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999);
    refVal   = 2.4999995;
    errBound = 2.4999995 * std::max(eps/(1.-0.9999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1NormalDoubles, reg, x = 1-10^-7 = 0.9999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999999);
    refVal   = 2.49999995;
    errBound = 2.49999995 * std::max(eps/(1.-0.99999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1NormalDoubles, reg, x = 1-10^-8 = 0.99999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999);
    refVal   = 2.499999995;
    errBound = 2.499999995 * std::max(eps/(1.-0.999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1NormalDoubles, reg, x = 1-10^-9 = 0.999999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999999);
    refVal   = 2.4999999995;
    errBound = 2.4999999995 * std::max(eps/(1.-0.9999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs1NormalDoubles, reg, x = 1-10^-10 = 0.9999999999";
  }
}



TEST(polAQgXspace,TruncAs2PowersOf4RegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto eval = polAQg_reg.truncate(2);

  {
    testVal  = eval(testAs,testLM,testNF,9.31322574615478515625e-10);
    refVal   = -268.577319978029352221066082097;
    errBound = 268.577319978029352221066082097 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2PowersOf4, reg, x = 4^-15 = 9.31322574615478515625e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.7252902984619140625e-9);
    refVal   = -277.380668339246718151156700838;
    errBound = 277.380668339246718151156700838 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2PowersOf4, reg, x = 4^-14 = 3.7252902984619140625e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.490116119384765625e-8);
    refVal   = -272.572891559139625438049264763;
    errBound = 272.572891559139625438049264763 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2PowersOf4, reg, x = 4^-13 = 1.490116119384765625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,5.9604644775390625e-8);
    refVal   = -256.596175509571716737199503766;
    errBound = 256.596175509571716737199503766 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2PowersOf4, reg, x = 4^-12 = 5.9604644775390625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,2.384185791015625e-7);
    refVal   = -231.892719309123044090683194432;
    errBound = 231.892719309123044090683194432 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2PowersOf4, reg, x = 4^-11 = 2.384185791015625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,9.5367431640625e-7);
    refVal   = -200.904764665480342855082095473;
    errBound = 200.904764665480342855082095473 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2PowersOf4, reg, x = 4^-10 = 9.5367431640625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.814697265625e-6);
    refVal   = -166.074674085604130645572486437;
    errBound = 166.074674085604130645572486437 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2PowersOf4, reg, x = 4^-9 = 3.814697265625e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0000152587890625);
    refVal   = -129.845078173834009227237836993;
    errBound = 129.845078173834009227237836993 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2PowersOf4, reg, x = 4^-8 = 0.0000152587890625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00006103515625);
    refVal   = -94.6588118943817963874674517063;
    errBound = 94.6588118943817963874674517063 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2PowersOf4, reg, x = 4^-7 = 0.00006103515625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.000244140625);
    refVal   = -62.956262895918789051236548;
    errBound = 62.956262895918789051236548 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2PowersOf4, reg, x = 4^-6 = 0.000244140625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0009765625);
    refVal   = -37.1567326302197993053814709359;
    errBound = 37.1567326302197993053814709359 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2PowersOf4, reg, x = 4^-5 = 0.0009765625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00390625);
    refVal   = -19.5605439171586284734761051431;
    errBound = 19.5605439171586284734761051431 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2PowersOf4, reg, x = 4^-4 = 0.00390625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.015625);
    refVal   = -11.9072888555976279350536326671;
    errBound = 11.9072888555976279350536326671 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2PowersOf4, reg, x = 4^-3 = 0.015625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0625);
    refVal   = -13.589717840111394255172493168;
    errBound = 13.589717840111394255172493168 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2PowersOf4, reg, x = 4^-2 = 0.0625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9375);
    refVal   = 28.1624000832133889802259868509;
    errBound = 28.1624000832133889802259868509 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2PowersOf4, reg, x = 1-4^-2 = 0.9375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.984375);
    refVal   = 29.3822359471102586912707628552;
    errBound = 29.3822359471102586912707628552 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2PowersOf4, reg, x = 1-4^-3 = 0.984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99609375);
    refVal   = 23.2173496643018839668565967152;
    errBound = 23.2173496643018839668565967152 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2PowersOf4, reg, x = 1-4^-4 = 0.99609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9990234375);
    refVal   = 12.8024787468191045739985706114;
    errBound = 12.8024787468191045739985706114 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2PowersOf4, reg, x = 1-4^-5 = 0.9990234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999755859375);
    refVal   = -0.219355810713932759931666128652;
    errBound = 0.219355810713932759931666128652 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2PowersOf4, reg, x = 1-4^-6 = 0.999755859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99993896484375);
    refVal   = -14.6162759501893970275520214601;
    errBound = 14.6162759501893970275520214601 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2PowersOf4, reg, x = 1-4^-7 = 0.99993896484375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999847412109375);
    refVal   = -29.254089323029100445893499638;
    errBound = 29.254089323029100445893499638 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2PowersOf4, reg, x = 1-4^-8 = 0.9999847412109375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999996185302734375);
    refVal   = -43.0188974711533031978522694375;
    errBound = 43.0188974711533031978522694375 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2PowersOf4, reg, x = 1-4^-9 = 0.999996185302734375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999904632568359375);
    refVal   = -54.8003151776154838301905184168;
    errBound = 54.8003151776154838301905184168 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2PowersOf4, reg, x = 1-4^-10 = 0.99999904632568359375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999997615814208984375);
    refVal   = -63.4883574293064121391431458917;
    errBound = 63.4883574293064121391431458917 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2PowersOf4, reg, x = 1-4^-11 = 0.9999997615814208984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999940395355224609375);
    refVal   = -67.9730085343871963729299670186;
    errBound = 67.9730085343871963729299670186 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2PowersOf4, reg, x = 1-4^-12 = 0.999999940395355224609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999998509883880615234375);
    refVal   = -67.1442123343487741967642971738;
    errBound = 67.1442123343487741967642971738 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2PowersOf4, reg, x = 1-4^-13 = 0.99999998509883880615234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999962747097015380859375);
    refVal   = -59.8918949964987155052102270829;
    errBound = 59.8918949964987155052102270829 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2PowersOf4, reg, x = 1-4^-14 = 0.9999999962747097015380859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999068677425384521484375);
    refVal   = -45.1059767081641960799897239696;
    errBound = 45.1059767081641960799897239696 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2PowersOf4, reg, x = 1-4^-15 = 0.999999999068677425384521484375";
  }
}



TEST(polAQgXspace,TruncAs2NormalDoublesRegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto eval = polAQg_reg.truncate(2);

  {
    testVal  = eval(testAs,testLM,testNF,1.e-10);
    refVal   = -219.648128592756308796864079824;
    errBound = 219.648128592756308796864079824 * std::max(eps/(1.-1.e-10),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2NormalDoubles, reg, x = 10^-10 = 1.e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-9);
    refVal   = -269.399117513294107422029297948;
    errBound = 269.399117513294107422029297948 * std::max(eps/(1.-1.e-9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2NormalDoubles, reg, x = 10^-9 = 1.e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-8);
    refVal   = -275.208009942799953851665261567;
    errBound = 275.208009942799953851665261567 * std::max(eps/(1.-1.e-8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2NormalDoubles, reg, x = 10^-8 = 1.e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-7);
    refVal   = -248.265552802056355517592012477;
    errBound = 248.265552802056355517592012477 * std::max(eps/(1.-1.e-7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2NormalDoubles, reg, x = 10^-7 = 1.e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-6);
    refVal   = -199.762597755575586859711417579;
    errBound = 199.762597755575586859711417579 * std::max(eps/(1.-1.e-6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2NormalDoubles, reg, x = 10^-6 = 1.e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-5);
    refVal   = -140.890619552032974536655124794;
    errBound = 140.890619552032974536655124794 * std::max(eps/(1.-1.e-5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2NormalDoubles, reg, x = 10^-5 = 1.e-5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0001);
    refVal   = -82.8425035030529701224997453483;
    errBound = 82.8425035030529701224997453483 * std::max(eps/(1.-0.0001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2NormalDoubles, reg, x = 10^-4 = 0.0001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.001);
    refVal   = -36.7805412079531507639547177931;
    errBound = 36.7805412079531507639547177931 * std::max(eps/(1.-0.001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2NormalDoubles, reg, x = 10^-3 = 0.001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.01);
    refVal   = -13.2378402733070595315971107401;
    errBound = 13.2378402733070595315971107401 * std::max(eps/(1.-0.01),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2NormalDoubles, reg, x = 10^-2 = 0.01";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.05);
    refVal   = -12.8883134015802743907714330962;
    errBound = 12.8883134015802743907714330962 * std::max(eps/(1.-0.05),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2NormalDoubles, reg, x = 0.05 = 0.05";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.1);
    refVal   = -15.1314462873937055311662019907;
    errBound = 15.1314462873937055311662019907 * std::max(eps/(1.-0.1),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2NormalDoubles, reg, x = 0.1 = 0.1";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.15);
    refVal   = -15.9751329561440742507098971252;
    errBound = 15.9751329561440742507098971252 * std::max(eps/(1.-0.15),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2NormalDoubles, reg, x = 0.15 = 0.15";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.2);
    refVal   = -15.8234200575319721314540185312;
    errBound = 15.8234200575319721314540185312 * std::max(eps/(1.-0.2),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2NormalDoubles, reg, x = 0.2 = 0.2";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.25);
    refVal   = -14.9692422212437361254608394254;
    errBound = 14.9692422212437361254608394254 * std::max(eps/(1.-0.25),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2NormalDoubles, reg, x = 0.25 = 0.25";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.3);
    refVal   = -13.5902938390980799858718325082;
    errBound = 13.5902938390980799858718325082 * std::max(eps/(1.-0.3),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2NormalDoubles, reg, x = 0.3 = 0.3";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.35);
    refVal   = -11.7983882405508167846752529842;
    errBound = 11.7983882405508167846752529842 * std::max(eps/(1.-0.35),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2NormalDoubles, reg, x = 0.35 = 0.35";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.4);
    refVal   = -9.66653419879803405086734860169;
    errBound = 9.66653419879803405086734860169 * std::max(eps/(1.-0.4),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2NormalDoubles, reg, x = 0.4 = 0.4";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.45);
    refVal   = -7.24339729687010771280786272692;
    errBound = 7.24339729687010771280786272692 * std::max(eps/(1.-0.45),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2NormalDoubles, reg, x = 0.45 = 0.45";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.5);
    refVal   = -4.56138898834694203746166427746;
    errBound = 4.56138898834694203746166427746 * std::max(eps/(1.-0.5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2NormalDoubles, reg, x = 0.5 = 0.5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.55);
    refVal   = -1.64140516193688042026350557476;
    errBound = 1.64140516193688042026350557476 * std::max(eps/(1.-0.55),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2NormalDoubles, reg, x = 0.55 = 0.55";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.6);
    refVal   = 1.50428705471523194305028256832;
    errBound = 1.50428705471523194305028256832 * std::max(eps/(1.-0.6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2NormalDoubles, reg, x = 0.6 = 0.6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.65);
    refVal   = 4.87019740639888828102407734781;
    errBound = 4.87019740639888828102407734781 * std::max(eps/(1.-0.65),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2NormalDoubles, reg, x = 0.65 = 0.65";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.7);
    refVal   = 8.45622612587937083475102597094;
    errBound = 8.45622612587937083475102597094 * std::max(eps/(1.-0.7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2NormalDoubles, reg, x = 0.7 = 0.7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.75);
    refVal   = 12.2660364074760804940601960741;
    errBound = 12.2660364074760804940601960741 * std::max(eps/(1.-0.75),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2NormalDoubles, reg, x = 0.75 = 0.75";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.8);
    refVal   = 16.3034859944280049634773316682;
    errBound = 16.3034859944280049634773316682 * std::max(eps/(1.-0.8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2NormalDoubles, reg, x = 0.8 = 0.8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.85);
    refVal   = 20.5610520803132559412629591057;
    errBound = 20.5610520803132559412629591057 * std::max(eps/(1.-0.85),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2NormalDoubles, reg, x = 0.85 = 0.85";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9);
    refVal   = 24.970825784519866077030802692;
    errBound = 24.970825784519866077030802692 * std::max(eps/(1.-0.9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2NormalDoubles, reg, x = 0.9 = 0.9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.95);
    refVal   = 29.0844969127567032058507042949;
    errBound = 29.0844969127567032058507042949 * std::max(eps/(1.-0.95),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2NormalDoubles, reg, x = 0.95 = 0.95";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99);
    refVal   = 28.0017915718054738354261443316;
    errBound = 28.0017915718054738354261443316 * std::max(eps/(1.-0.99),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2NormalDoubles, reg, x = 1-10^-2 = 0.99";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999);
    refVal   = 13.0063746719023947581880189116;
    errBound = 13.0063746719023947581880189116 * std::max(eps/(1.-0.999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2NormalDoubles, reg, x = 1-10^-3 = 0.999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999);
    refVal   = -9.40326274426414337455815301544;
    errBound = 9.40326274426414337455815301544 * std::max(eps/(1.-0.9999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2NormalDoubles, reg, x = 1-10^-4 = 0.9999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999);
    refVal   = -33.5935681051980943556364149631;
    errBound = 33.5935681051980943556364149631 * std::max(eps/(1.-0.99999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2NormalDoubles, reg, x = 1-10^-5 = 0.99999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999);
    refVal   = -54.4419934148818120588515335188;
    errBound = 54.4419934148818120588515335188 * std::max(eps/(1.-0.999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2NormalDoubles, reg, x = 1-10^-6 = 0.999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999);
    refVal   = -66.8611702105102290610309738964;
    errBound = 66.8611702105102290610309738964 * std::max(eps/(1.-0.9999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2NormalDoubles, reg, x = 1-10^-7 = 0.9999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999999);
    refVal   = -65.76465390091176598964621635;
    errBound = 65.76465390091176598964621635 * std::max(eps/(1.-0.99999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2NormalDoubles, reg, x = 1-10^-8 = 0.99999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999);
    refVal   = -46.0658004494490074032563223599;
    errBound = 46.0658004494490074032563223599 * std::max(eps/(1.-0.999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2NormalDoubles, reg, x = 1-10^-9 = 0.999999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999999);
    refVal   = -2.67792013814202834192605025602;
    errBound = 2.67792013814202834192605025602 * std::max(eps/(1.-0.9999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQg, TruncAs2NormalDoubles, reg, x = 1-10^-10 = 0.9999999999";
  }
}



TEST(polAQgNspace,FullMoments)
{
  double refVal = 0., errBound = 0.;
  auto& rpd = polAQg;
  integration_engine_gsl engine;
  auto mom = make_mellin_moment(rpd, engine, testAs, testLM, testNF);
  std::tuple<integration_status, double, double> res;

  {
    SCOPED_TRACE("polAQg Full Mellin moment N=3");
    refVal = 51.3560206671198806249251011596;
    errBound = 51.3560206671198806249251011596 * 1.e-9;
    res = mom.integrate(3, 0., 1.e-9);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAQg Full Mellin moment N=5");
    refVal = 47.4909719813215682887281840246;
    errBound = 47.4909719813215682887281840246 * 1.e-9;
    res = mom.integrate(5, 0., 1.e-9);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAQg Full Mellin moment N=7");
    refVal = 40.9926412761823619129362416508;
    errBound = 40.9926412761823619129362416508 * 1.e-9;
    res = mom.integrate(7, 0., 1.e-9);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

}


TEST(polAQgNspace,TruncAs1Moments)
{
  double refVal = 0., errBound = 0.;
  auto rpd = polAQg.truncate(1);
  integration_engine_gsl engine;
  auto mom = make_mellin_moment(rpd, engine, testAs, testLM, testNF);
  std::tuple<integration_status, double, double> res;

  {
    SCOPED_TRACE("polAQg TruncAs1 Mellin moment N=3");
    refVal = 0.416666666666666666666666666667;
    errBound = 0.416666666666666666666666666667 * 1.e-9;
    res = mom.integrate(3, 0., 1.e-9);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAQg TruncAs1 Mellin moment N=5");
    refVal = 0.333333333333333333333333333333;
    errBound = 0.333333333333333333333333333333 * 1.e-9;
    res = mom.integrate(5, 0., 1.e-9);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAQg TruncAs1 Mellin moment N=7");
    refVal = 0.267857142857142857142857142857;
    errBound = 0.267857142857142857142857142857 * 1.e-9;
    res = mom.integrate(7, 0., 1.e-9);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

}


TEST(polAQgNspace,TruncAs2Moments)
{
  double refVal = 0., errBound = 0.;
  auto rpd = polAQg.truncate(2);
  integration_engine_gsl engine;
  auto mom = make_mellin_moment(rpd, engine, testAs, testLM, testNF);
  std::tuple<integration_status, double, double> res;

  {
    SCOPED_TRACE("polAQg TruncAs2 Mellin moment N=3");
    refVal = 4.44966081532921810699588477366;
    errBound = 4.44966081532921810699588477366 * 1.e-9;
    res = mom.integrate(3, 0., 1.e-9);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAQg TruncAs2 Mellin moment N=5");
    refVal = 3.79404722222222222222222222222;
    errBound = 3.79404722222222222222222222222 * 1.e-9;
    res = mom.integrate(5, 0., 1.e-9);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAQg TruncAs2 Mellin moment N=7");
    refVal = 3.12941521954084685650526272489;
    errBound = 3.12941521954084685650526272489 * 1.e-9;
    res = mom.integrate(7, 0., 1.e-9);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

}
