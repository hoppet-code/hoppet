/******************************************************************************
 * This file is part of libome                                                *
 * Copyright (C) 2025 Arnd Behring, Kay Schoenwald                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 ******************************************************************************/

#include <limits>
#include <gtest/gtest.h>
#include <ome/integration_engine_gsl.h>
#include <ome/mellin.h>
#include <ome/AQg.h>

using namespace ome;

const double eps = std::numeric_limits<double>::epsilon();
const double testAs = 0.25;
const double testLM = -5.;
const double testNF = 3.;



TEST(AQgXspace,FullPowersOf4RegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto& eval = AQg_reg;

  {
    testVal  = eval(testAs,testLM,testNF,9.31322574615478515625e-10);
    refVal   = 6.07631635691092879708538178023e11;
    errBound = 6.07631635691092879708538178023e11 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullPowersOf4, reg, x = 4^-15 = 9.31322574615478515625e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.7252902984619140625e-9);
    refVal   = 1.38403386030002196518754010239e11;
    errBound = 1.38403386030002196518754010239e11 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullPowersOf4, reg, x = 4^-14 = 3.7252902984619140625e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.490116119384765625e-8);
    refVal   = 3.12247120727017365463248837449e10;
    errBound = 3.12247120727017365463248837449e10 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullPowersOf4, reg, x = 4^-13 = 1.490116119384765625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,5.9604644775390625e-8);
    refVal   = 6.96214302240950048589384058537e9;
    errBound = 6.96214302240950048589384058537e9 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullPowersOf4, reg, x = 4^-12 = 5.9604644775390625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,2.384185791015625e-7);
    refVal   = 1.52952690748449615363924419804e9;
    errBound = 1.52952690748449615363924419804e9 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullPowersOf4, reg, x = 4^-11 = 2.384185791015625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,9.5367431640625e-7);
    refVal   = 3.29630031016696155910945976084e8;
    errBound = 3.29630031016696155910945976084e8 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullPowersOf4, reg, x = 4^-10 = 9.5367431640625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.814697265625e-6);
    refVal   = 6.92203129502850852469556118051e7;
    errBound = 6.92203129502850852469556118051e7 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullPowersOf4, reg, x = 4^-9 = 3.814697265625e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0000152587890625);
    refVal   = 1.40090092749370776310785896765e7;
    errBound = 1.40090092749370776310785896765e7 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullPowersOf4, reg, x = 4^-8 = 0.0000152587890625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00006103515625);
    refVal   = 2.6788764440810881303709871579e6;
    errBound = 2.6788764440810881303709871579e6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullPowersOf4, reg, x = 4^-7 = 0.00006103515625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.000244140625);
    refVal   = 464409.823863053462927957465671;
    errBound = 464409.823863053462927957465671 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullPowersOf4, reg, x = 4^-6 = 0.000244140625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0009765625);
    refVal   = 65216.3522979666494418657024416;
    errBound = 65216.3522979666494418657024416 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullPowersOf4, reg, x = 4^-5 = 0.0009765625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00390625);
    refVal   = 3949.68568571336139404881981445;
    errBound = 3949.68568571336139404881981445 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullPowersOf4, reg, x = 4^-4 = 0.00390625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.015625);
    refVal   = -1794.75517625650155736063067342;
    errBound = 1794.75517625650155736063067342 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullPowersOf4, reg, x = 4^-3 = 0.015625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0625);
    refVal   = -892.005815468712993011917097751;
    errBound = 892.005815468712993011917097751 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullPowersOf4, reg, x = 4^-2 = 0.0625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9375);
    refVal   = 384.346555385233451623714468211;
    errBound = 384.346555385233451623714468211 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullPowersOf4, reg, x = 1-4^-2 = 0.9375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.984375);
    refVal   = 460.817956794887586909664127309;
    errBound = 460.817956794887586909664127309 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullPowersOf4, reg, x = 1-4^-3 = 0.984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99609375);
    refVal   = 416.949973997442022424074187943;
    errBound = 416.949973997442022424074187943 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullPowersOf4, reg, x = 1-4^-4 = 0.99609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9990234375);
    refVal   = 307.157950667584481152867756493;
    errBound = 307.157950667584481152867756493 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullPowersOf4, reg, x = 1-4^-5 = 0.9990234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999755859375);
    refVal   = 171.964295761817018607989219838;
    errBound = 171.964295761817018607989219838 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullPowersOf4, reg, x = 1-4^-6 = 0.999755859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99993896484375);
    refVal   = 29.0106808833009061514040194239;
    errBound = 29.0106808833009061514040194239 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullPowersOf4, reg, x = 1-4^-7 = 0.99993896484375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999847412109375);
    refVal   = -140.3001883081381902386164391276;
    errBound = 140.3001883081381902386164391276 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullPowersOf4, reg, x = 1-4^-8 = 0.9999847412109375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999996185302734375);
    refVal   = -407.281094416377159584889464642;
    errBound = 407.281094416377159584889464642 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullPowersOf4, reg, x = 1-4^-9 = 0.999996185302734375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999904632568359375);
    refVal   = -913.428261763851286119466691645;
    errBound = 913.428261763851286119466691645 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullPowersOf4, reg, x = 1-4^-10 = 0.99999904632568359375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999997615814208984375);
    refVal   = -1888.13663000978420459376795852;
    errBound = 1888.13663000978420459376795852 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullPowersOf4, reg, x = 1-4^-11 = 0.9999997615814208984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999940395355224609375);
    refVal   = -3666.46757177145462861406403409;
    errBound = 3666.46757177145462861406403409 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullPowersOf4, reg, x = 1-4^-12 = 0.999999940395355224609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999998509883880615234375);
    refVal   = -6706.92598807459984471457184819;
    errBound = 6706.92598807459984471457184819 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullPowersOf4, reg, x = 1-4^-13 = 0.99999998509883880615234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999962747097015380859375);
    refVal   = -11609.238634061242399311787843;
    errBound = 11609.238634061242399311787843 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullPowersOf4, reg, x = 1-4^-14 = 0.9999999962747097015380859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999068677425384521484375);
    refVal   = -19132.132433392712338069755134;
    errBound = 19132.132433392712338069755134 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullPowersOf4, reg, x = 1-4^-15 = 0.999999999068677425384521484375";
  }
}



TEST(AQgXspace,FullNormalDoublesRegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto& eval = AQg_reg;

  {
    testVal  = eval(testAs,testLM,testNF,1.e-10);
    refVal   = 6.46879231932386526745861353674e12;
    errBound = 6.46879231932386526745861353674e12 * std::max(eps/(1.-1.e-10),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullNormalDoubles, reg, x = 10^-10 = 1.e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-9);
    refVal   = 5.63319060163524021851148708532e11;
    errBound = 5.63319060163524021851148708532e11 * std::max(eps/(1.-1.e-9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullNormalDoubles, reg, x = 10^-9 = 1.e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-8);
    refVal   = 4.79758802897807383214982031078e10;
    errBound = 4.79758802897807383214982031078e10 * std::max(eps/(1.-1.e-8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullNormalDoubles, reg, x = 10^-8 = 1.e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-7);
    refVal   = 3.96198348463796593754535174398e9;
    errBound = 3.96198348463796593754535174398e9 * std::max(eps/(1.-1.e-7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullNormalDoubles, reg, x = 10^-7 = 1.e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-6);
    refVal   = 3.12638398846795554115194040911e8;
    errBound = 3.12638398846795554115194040911e8 * std::max(eps/(1.-1.e-6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullNormalDoubles, reg, x = 10^-6 = 1.e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-5);
    refVal   = 2.29088939628048290080019504989e7;
    errBound = 2.29088939628048290080019504989e7 * std::max(eps/(1.-1.e-5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullNormalDoubles, reg, x = 10^-5 = 1.e-5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0001);
    refVal   = 1.45630622196767787209126368634e6;
    errBound = 1.45630622196767787209126368634e6 * std::max(eps/(1.-0.0001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullNormalDoubles, reg, x = 10^-4 = 0.0001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.001);
    refVal   = 62845.4446951065555619347114162;
    errBound = 62845.4446951065555619347114162 * std::max(eps/(1.-0.001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullNormalDoubles, reg, x = 10^-3 = 0.001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.01);
    refVal   = -1498.48148289126393383003066007;
    errBound = 1498.48148289126393383003066007 * std::max(eps/(1.-0.01),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullNormalDoubles, reg, x = 10^-2 = 0.01";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.05);
    refVal   = -1066.13371752555342365436472312;
    errBound = 1066.13371752555342365436472312 * std::max(eps/(1.-0.05),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullNormalDoubles, reg, x = 0.05 = 0.05";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.1);
    refVal   = -567.68167582039480514532097545;
    errBound = 567.68167582039480514532097545 * std::max(eps/(1.-0.1),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullNormalDoubles, reg, x = 0.1 = 0.1";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.15);
    refVal   = -344.281970233330701552219147323;
    errBound = 344.281970233330701552219147323 * std::max(eps/(1.-0.15),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullNormalDoubles, reg, x = 0.15 = 0.15";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.2);
    refVal   = -218.021754540676991737462811037;
    errBound = 218.021754540676991737462811037 * std::max(eps/(1.-0.2),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullNormalDoubles, reg, x = 0.2 = 0.2";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.25);
    refVal   = -137.554685102874615888659253396;
    errBound = 137.554685102874615888659253396 * std::max(eps/(1.-0.25),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullNormalDoubles, reg, x = 0.25 = 0.25";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.3);
    refVal   = -82.3298274799967202828988818219;
    errBound = 82.3298274799967202828988818219 * std::max(eps/(1.-0.3),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullNormalDoubles, reg, x = 0.3 = 0.3";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.35);
    refVal   = -42.2601898421571268619679960316;
    errBound = 42.2601898421571268619679960316 * std::max(eps/(1.-0.35),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullNormalDoubles, reg, x = 0.35 = 0.35";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.4);
    refVal   = -11.6122813903340581830510412802;
    errBound = 11.6122813903340581830510412802 * std::max(eps/(1.-0.4),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullNormalDoubles, reg, x = 0.4 = 0.4";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.45);
    refVal   = 13.2862671311020975210235579735;
    errBound = 13.2862671311020975210235579735 * std::max(eps/(1.-0.45),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullNormalDoubles, reg, x = 0.45 = 0.45";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.5);
    refVal   = 35.0502958826255689243032965156;
    errBound = 35.0502958826255689243032965156 * std::max(eps/(1.-0.5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullNormalDoubles, reg, x = 0.5 = 0.5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.55);
    refVal   = 55.7334919014703752836221924243;
    errBound = 55.7334919014703752836221924243 * std::max(eps/(1.-0.55),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullNormalDoubles, reg, x = 0.55 = 0.55";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.6);
    refVal   = 77.1039228846050193373831877955;
    errBound = 77.1039228846050193373831877955 * std::max(eps/(1.-0.6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullNormalDoubles, reg, x = 0.6 = 0.6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.65);
    refVal   = 100.819400492566687115459649859;
    errBound = 100.819400492566687115459649859 * std::max(eps/(1.-0.65),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullNormalDoubles, reg, x = 0.65 = 0.65";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.7);
    refVal   = 128.5619819162091303377505793317;
    errBound = 128.5619819162091303377505793317 * std::max(eps/(1.-0.7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullNormalDoubles, reg, x = 0.7 = 0.7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.75);
    refVal   = 162.168791292895660555077555452;
    errBound = 162.168791292895660555077555452 * std::max(eps/(1.-0.75),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullNormalDoubles, reg, x = 0.75 = 0.75";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.8);
    refVal   = 203.792888769135024707421429032;
    errBound = 203.792888769135024707421429032 * std::max(eps/(1.-0.8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullNormalDoubles, reg, x = 0.8 = 0.8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.85);
    refVal   = 256.13321035257629555994338157;
    errBound = 256.13321035257629555994338157 * std::max(eps/(1.-0.85),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullNormalDoubles, reg, x = 0.85 = 0.85";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9);
    refVal   = 322.723215463771376867218056767;
    errBound = 322.723215463771376867218056767 * std::max(eps/(1.-0.9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullNormalDoubles, reg, x = 0.9 = 0.9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.95);
    refVal   = 407.046722480520275535268389533;
    errBound = 407.046722480520275535268389533 * std::max(eps/(1.-0.95),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullNormalDoubles, reg, x = 0.95 = 0.95";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99);
    refVal   = 457.033763573451576686167545172;
    errBound = 457.033763573451576686167545172 * std::max(eps/(1.-0.99),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullNormalDoubles, reg, x = 1-10^-2 = 0.99";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999);
    refVal   = 309.339170632094793499773841642;
    errBound = 309.339170632094793499773841642 * std::max(eps/(1.-0.999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullNormalDoubles, reg, x = 1-10^-3 = 0.999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999);
    refVal   = 81.0304249295600728288441775597;
    errBound = 81.0304249295600728288441775597 * std::max(eps/(1.-0.9999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullNormalDoubles, reg, x = 1-10^-4 = 0.9999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999);
    refVal   = -206.342177213841176309414351552;
    errBound = 206.342177213841176309414351552 * std::max(eps/(1.-0.99999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullNormalDoubles, reg, x = 1-10^-5 = 0.99999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999);
    refVal   = -889.950435246451335824762022382;
    errBound = 889.950435246451335824762022382 * std::max(eps/(1.-0.999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullNormalDoubles, reg, x = 1-10^-6 = 0.999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999);
    refVal   = -2882.51533174206001838401190582;
    errBound = 2882.51533174206001838401190582 * std::max(eps/(1.-0.9999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullNormalDoubles, reg, x = 1-10^-7 = 0.9999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999999);
    refVal   = -7896.37936992380232364282193417;
    errBound = 7896.37936992380232364282193417 * std::max(eps/(1.-0.99999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullNormalDoubles, reg, x = 1-10^-8 = 0.99999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999);
    refVal   = -18668.2181475539242014546048099;
    errBound = 18668.2181475539242014546048099 * std::max(eps/(1.-0.999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullNormalDoubles, reg, x = 1-10^-9 = 0.999999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999999);
    refVal   = -39183.7854646474923908923015858;
    errBound = 39183.7854646474923908923015858 * std::max(eps/(1.-0.9999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, FullNormalDoubles, reg, x = 1-10^-10 = 0.9999999999";
  }
}



TEST(AQgXspace,TruncAs1PowersOf4RegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto eval = AQg_reg.truncate(1);

  {
    testVal  = eval(testAs,testLM,testNF,9.31322574615478515625e-10);
    refVal   = 2.49999999534338713125941611182;
    errBound = 2.49999999534338713125941611182 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1PowersOf4, reg, x = 4^-15 = 9.31322574615478515625e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.7252902984619140625e-9);
    refVal   = 2.49999998137354857707936872657;
    errBound = 2.49999998137354857707936872657 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1PowersOf4, reg, x = 4^-14 = 3.7252902984619140625e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.490116119384765625e-8);
    refVal   = 2.49999992549419514098474337516;
    errBound = 2.49999992549419514098474337516 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1PowersOf4, reg, x = 4^-13 = 1.490116119384765625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,5.9604644775390625e-8);
    refVal   = 2.4999997019767938866152690025;
    errBound = 2.4999997019767938866152690025 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1PowersOf4, reg, x = 4^-12 = 5.9604644775390625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,2.384185791015625e-7);
    refVal   = 2.49999880790738870928180404007;
    errBound = 2.49999880790738870928180404007 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1PowersOf4, reg, x = 4^-11 = 2.384185791015625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,9.5367431640625e-7);
    refVal   = 2.49999523163296544225886464119;
    errBound = 2.49999523163296544225886464119 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1PowersOf4, reg, x = 4^-10 = 9.5367431640625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.814697265625e-6);
    refVal   = 2.49998092658643145114183425903;
    errBound = 2.49998092658643145114183425903 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1PowersOf4, reg, x = 4^-9 = 3.814697265625e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0000152587890625);
    refVal   = 2.49992370721884071826934814453;
    errBound = 2.49992370721884071826934814453 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1PowersOf4, reg, x = 4^-8 = 0.0000152587890625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00006103515625);
    refVal   = 2.4996948428452014923095703125;
    errBound = 2.4996948428452014923095703125 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1PowersOf4, reg, x = 4^-7 = 0.00006103515625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.000244140625);
    refVal   = 2.498779594898223876953125;
    errBound = 2.498779594898223876953125 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1PowersOf4, reg, x = 4^-6 = 0.000244140625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0009765625);
    refVal   = 2.49512195587158203125;
    errBound = 2.49512195587158203125 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1PowersOf4, reg, x = 4^-5 = 0.0009765625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00390625);
    refVal   = 2.4805450439453125;
    errBound = 2.4805450439453125 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1PowersOf4, reg, x = 4^-4 = 0.00390625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.015625);
    refVal   = 2.423095703125;
    errBound = 2.423095703125 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1PowersOf4, reg, x = 4^-3 = 0.015625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0625);
    refVal   = 2.20703125;
    errBound = 2.20703125 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1PowersOf4, reg, x = 4^-2 = 0.0625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9375);
    refVal   = 2.20703125;
    errBound = 2.20703125 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1PowersOf4, reg, x = 1-4^-2 = 0.9375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.984375);
    refVal   = 2.423095703125;
    errBound = 2.423095703125 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1PowersOf4, reg, x = 1-4^-3 = 0.984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99609375);
    refVal   = 2.4805450439453125;
    errBound = 2.4805450439453125 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1PowersOf4, reg, x = 1-4^-4 = 0.99609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9990234375);
    refVal   = 2.49512195587158203125;
    errBound = 2.49512195587158203125 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1PowersOf4, reg, x = 1-4^-5 = 0.9990234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999755859375);
    refVal   = 2.498779594898223876953125;
    errBound = 2.498779594898223876953125 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1PowersOf4, reg, x = 1-4^-6 = 0.999755859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99993896484375);
    refVal   = 2.4996948428452014923095703125;
    errBound = 2.4996948428452014923095703125 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1PowersOf4, reg, x = 1-4^-7 = 0.99993896484375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999847412109375);
    refVal   = 2.49992370721884071826934814453;
    errBound = 2.49992370721884071826934814453 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1PowersOf4, reg, x = 1-4^-8 = 0.9999847412109375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999996185302734375);
    refVal   = 2.49998092658643145114183425903;
    errBound = 2.49998092658643145114183425903 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1PowersOf4, reg, x = 1-4^-9 = 0.999996185302734375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999904632568359375);
    refVal   = 2.49999523163296544225886464119;
    errBound = 2.49999523163296544225886464119 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1PowersOf4, reg, x = 1-4^-10 = 0.99999904632568359375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999997615814208984375);
    refVal   = 2.49999880790738870928180404007;
    errBound = 2.49999880790738870928180404007 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1PowersOf4, reg, x = 1-4^-11 = 0.9999997615814208984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999940395355224609375);
    refVal   = 2.4999997019767938866152690025;
    errBound = 2.4999997019767938866152690025 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1PowersOf4, reg, x = 1-4^-12 = 0.999999940395355224609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999998509883880615234375);
    refVal   = 2.49999992549419514098474337516;
    errBound = 2.49999992549419514098474337516 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1PowersOf4, reg, x = 1-4^-13 = 0.99999998509883880615234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999962747097015380859375);
    refVal   = 2.49999998137354857707936872657;
    errBound = 2.49999998137354857707936872657 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1PowersOf4, reg, x = 1-4^-14 = 0.9999999962747097015380859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999068677425384521484375);
    refVal   = 2.49999999534338713125941611182;
    errBound = 2.49999999534338713125941611182 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1PowersOf4, reg, x = 1-4^-15 = 0.999999999068677425384521484375";
  }
}



TEST(AQgXspace,TruncAs1NormalDoublesRegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto eval = AQg_reg.truncate(1);

  {
    testVal  = eval(testAs,testLM,testNF,1.e-10);
    refVal   = 2.49999999950000000005;
    errBound = 2.49999999950000000005 * std::max(eps/(1.-1.e-10),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1NormalDoubles, reg, x = 10^-10 = 1.e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-9);
    refVal   = 2.499999995000000005;
    errBound = 2.499999995000000005 * std::max(eps/(1.-1.e-9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1NormalDoubles, reg, x = 10^-9 = 1.e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-8);
    refVal   = 2.4999999500000005;
    errBound = 2.4999999500000005 * std::max(eps/(1.-1.e-8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1NormalDoubles, reg, x = 10^-8 = 1.e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-7);
    refVal   = 2.49999950000005;
    errBound = 2.49999950000005 * std::max(eps/(1.-1.e-7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1NormalDoubles, reg, x = 10^-7 = 1.e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-6);
    refVal   = 2.499995000005;
    errBound = 2.499995000005 * std::max(eps/(1.-1.e-6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1NormalDoubles, reg, x = 10^-6 = 1.e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-5);
    refVal   = 2.4999500005;
    errBound = 2.4999500005 * std::max(eps/(1.-1.e-5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1NormalDoubles, reg, x = 10^-5 = 1.e-5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0001);
    refVal   = 2.49950005;
    errBound = 2.49950005 * std::max(eps/(1.-0.0001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1NormalDoubles, reg, x = 10^-4 = 0.0001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.001);
    refVal   = 2.495005;
    errBound = 2.495005 * std::max(eps/(1.-0.001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1NormalDoubles, reg, x = 10^-3 = 0.001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.01);
    refVal   = 2.4505;
    errBound = 2.4505 * std::max(eps/(1.-0.01),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1NormalDoubles, reg, x = 10^-2 = 0.01";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.05);
    refVal   = 2.2625;
    errBound = 2.2625 * std::max(eps/(1.-0.05),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1NormalDoubles, reg, x = 0.05 = 0.05";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.1);
    refVal   = 2.05;
    errBound = 2.05 * std::max(eps/(1.-0.1),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1NormalDoubles, reg, x = 0.1 = 0.1";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.15);
    refVal   = 1.8625;
    errBound = 1.8625 * std::max(eps/(1.-0.15),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1NormalDoubles, reg, x = 0.15 = 0.15";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.2);
    refVal   = 1.7;
    errBound = 1.7 * std::max(eps/(1.-0.2),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1NormalDoubles, reg, x = 0.2 = 0.2";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.25);
    refVal   = 1.5625;
    errBound = 1.5625 * std::max(eps/(1.-0.25),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1NormalDoubles, reg, x = 0.25 = 0.25";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.3);
    refVal   = 1.45;
    errBound = 1.45 * std::max(eps/(1.-0.3),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1NormalDoubles, reg, x = 0.3 = 0.3";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.35);
    refVal   = 1.3625;
    errBound = 1.3625 * std::max(eps/(1.-0.35),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1NormalDoubles, reg, x = 0.35 = 0.35";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.4);
    refVal   = 1.3;
    errBound = 1.3 * std::max(eps/(1.-0.4),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1NormalDoubles, reg, x = 0.4 = 0.4";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.45);
    refVal   = 1.2625;
    errBound = 1.2625 * std::max(eps/(1.-0.45),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1NormalDoubles, reg, x = 0.45 = 0.45";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.5);
    refVal   = 1.25;
    errBound = 1.25 * std::max(eps/(1.-0.5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1NormalDoubles, reg, x = 0.5 = 0.5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.55);
    refVal   = 1.2625;
    errBound = 1.2625 * std::max(eps/(1.-0.55),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1NormalDoubles, reg, x = 0.55 = 0.55";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.6);
    refVal   = 1.3;
    errBound = 1.3 * std::max(eps/(1.-0.6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1NormalDoubles, reg, x = 0.6 = 0.6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.65);
    refVal   = 1.3625;
    errBound = 1.3625 * std::max(eps/(1.-0.65),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1NormalDoubles, reg, x = 0.65 = 0.65";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.7);
    refVal   = 1.45;
    errBound = 1.45 * std::max(eps/(1.-0.7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1NormalDoubles, reg, x = 0.7 = 0.7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.75);
    refVal   = 1.5625;
    errBound = 1.5625 * std::max(eps/(1.-0.75),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1NormalDoubles, reg, x = 0.75 = 0.75";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.8);
    refVal   = 1.7;
    errBound = 1.7 * std::max(eps/(1.-0.8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1NormalDoubles, reg, x = 0.8 = 0.8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.85);
    refVal   = 1.8625;
    errBound = 1.8625 * std::max(eps/(1.-0.85),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1NormalDoubles, reg, x = 0.85 = 0.85";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9);
    refVal   = 2.05;
    errBound = 2.05 * std::max(eps/(1.-0.9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1NormalDoubles, reg, x = 0.9 = 0.9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.95);
    refVal   = 2.2625;
    errBound = 2.2625 * std::max(eps/(1.-0.95),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1NormalDoubles, reg, x = 0.95 = 0.95";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99);
    refVal   = 2.4505;
    errBound = 2.4505 * std::max(eps/(1.-0.99),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1NormalDoubles, reg, x = 1-10^-2 = 0.99";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999);
    refVal   = 2.495005;
    errBound = 2.495005 * std::max(eps/(1.-0.999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1NormalDoubles, reg, x = 1-10^-3 = 0.999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999);
    refVal   = 2.49950005;
    errBound = 2.49950005 * std::max(eps/(1.-0.9999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1NormalDoubles, reg, x = 1-10^-4 = 0.9999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999);
    refVal   = 2.4999500005;
    errBound = 2.4999500005 * std::max(eps/(1.-0.99999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1NormalDoubles, reg, x = 1-10^-5 = 0.99999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999);
    refVal   = 2.499995000005;
    errBound = 2.499995000005 * std::max(eps/(1.-0.999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1NormalDoubles, reg, x = 1-10^-6 = 0.999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999);
    refVal   = 2.49999950000005;
    errBound = 2.49999950000005 * std::max(eps/(1.-0.9999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1NormalDoubles, reg, x = 1-10^-7 = 0.9999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999999);
    refVal   = 2.4999999500000005;
    errBound = 2.4999999500000005 * std::max(eps/(1.-0.99999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1NormalDoubles, reg, x = 1-10^-8 = 0.99999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999);
    refVal   = 2.499999995000000005;
    errBound = 2.499999995000000005 * std::max(eps/(1.-0.999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1NormalDoubles, reg, x = 1-10^-9 = 0.999999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999999);
    refVal   = 2.49999999950000000005;
    errBound = 2.49999999950000000005 * std::max(eps/(1.-0.9999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs1NormalDoubles, reg, x = 1-10^-10 = 0.9999999999";
  }
}



TEST(AQgXspace,TruncAs2PowersOf4RegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto eval = AQg_reg.truncate(2);

  {
    testVal  = eval(testAs,testLM,testNF,9.31322574615478515625e-10);
    refVal   = -6.14418928387247514331748915598e9;
    errBound = 6.14418928387247514331748915598e9 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2PowersOf4, reg, x = 4^-15 = 9.31322574615478515625e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.7252902984619140625e-9);
    refVal   = -1.53604732329012400832425818381e9;
    errBound = 1.53604732329012400832425818381e9 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2PowersOf4, reg, x = 4^-14 = 3.7252902984619140625e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.490116119384765625e-8);
    refVal   = -3.84011847600205050218748455126e8;
    errBound = 3.84011847600205050218748455126e8 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2PowersOf4, reg, x = 4^-13 = 1.490116119384765625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,5.9604644775390625e-8);
    refVal   = -9.60029863568370851859786889939e7;
    errBound = 9.60029863568370851859786889939e7 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2PowersOf4, reg, x = 4^-12 = 5.9604644775390625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,2.384185791015625e-7);
    refVal   = -2.40007731141497792471013814016e7;
    errBound = 2.40007731141497792471013814016e7 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2PowersOf4, reg, x = 4^-11 = 2.384185791015625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,9.5367431640625e-7);
    refVal   = -6.00021742630808481689996408807e6;
    errBound = 6.00021742630808481689996408807e6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2PowersOf4, reg, x = 4^-10 = 9.5367431640625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.814697265625e-6);
    refVal   = -1.50007284757340591308658155557e6;
    errBound = 1.50007284757340591308658155557e6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2PowersOf4, reg, x = 4^-9 = 3.814697265625e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0000152587890625);
    refVal   = -375028.93238877165374908531045;
    errBound = 375028.93238877165374908531045 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2PowersOf4, reg, x = 4^-8 = 0.0000152587890625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00006103515625);
    refVal   = -93759.2350994929586321728928923;
    errBound = 93759.2350994929586321728928923 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2PowersOf4, reg, x = 4^-7 = 0.00006103515625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.000244140625);
    refVal   = -23433.306708578983269920239834;
    errBound = 23433.306708578983269920239834 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2PowersOf4, reg, x = 4^-6 = 0.000244140625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0009765625);
    refVal   = -5844.67592094030442992091801865;
    errBound = 5844.67592094030442992091801865 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2PowersOf4, reg, x = 4^-5 = 0.0009765625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00390625);
    refVal   = -1442.76402483873223464768167071;
    errBound = 1442.76402483873223464768167071 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2PowersOf4, reg, x = 4^-4 = 0.00390625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.015625);
    refVal   = -340.58427876557697119391979231;
    errBound = 340.58427876557697119391979231 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2PowersOf4, reg, x = 4^-3 = 0.015625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0625);
    refVal   = -66.0769126367088995291825843574;
    errBound = 66.0769126367088995291825843574 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2PowersOf4, reg, x = 4^-2 = 0.0625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9375);
    refVal   = 28.4629709443407041251895186472;
    errBound = 28.4629709443407041251895186472 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2PowersOf4, reg, x = 1-4^-2 = 0.9375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.984375);
    refVal   = 29.3988370665728829331778801139;
    errBound = 29.3988370665728829331778801139 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2PowersOf4, reg, x = 1-4^-3 = 0.984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99609375);
    refVal   = 23.2181086139172541817930507787;
    errBound = 23.2181086139172541817930507787 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2PowersOf4, reg, x = 1-4^-4 = 0.99609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9990234375);
    refVal   = 12.80250308007301914000623760548;
    errBound = 12.80250308007301914000623760548 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2PowersOf4, reg, x = 1-4^-5 = 0.9990234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999755859375);
    refVal   = -0.219355938171967756391615502983;
    errBound = 0.219355938171967756391615502983 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2PowersOf4, reg, x = 1-4^-6 = 0.999755859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99993896484375);
    refVal   = -14.6162760653605389328796959466;
    errBound = 14.6162760653605389328796959466 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2PowersOf4, reg, x = 1-4^-7 = 0.99993896484375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999847412109375);
    refVal   = -29.2540893366669810652682300385;
    errBound = 29.2540893366669810652682300385 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2PowersOf4, reg, x = 1-4^-8 = 0.9999847412109375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999996185302734375);
    refVal   = -43.0188974723594654624336355795;
    errBound = 43.0188974723594654624336355795 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2PowersOf4, reg, x = 1-4^-9 = 0.999996185302734375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999904632568359375);
    refVal   = -54.8003151777079176483285775;
    errBound = 54.8003151777079176483285775 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2PowersOf4, reg, x = 1-4^-10 = 0.99999904632568359375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999997615814208984375);
    refVal   = -63.4883574293128121026934032793;
    errBound = 63.4883574293128121026934032793 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2PowersOf4, reg, x = 1-4^-11 = 0.9999997615814208984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999940395355224609375);
    refVal   = -67.9730085343875997425479572344;
    errBound = 67.9730085343875997425479572344 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2PowersOf4, reg, x = 1-4^-12 = 0.999999940395355224609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999998509883880615234375);
    refVal   = -67.1442123343487969028642898389;
    errBound = 67.1442123343487969028642898389 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2PowersOf4, reg, x = 1-4^-13 = 0.99999998509883880615234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999962747097015380859375);
    refVal   = -59.8918949964987165672963209795;
    errBound = 59.8918949964987165672963209795 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2PowersOf4, reg, x = 1-4^-14 = 0.9999999962747097015380859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999068677425384521484375);
    refVal   = -45.1059767081641961095969812509;
    errBound = 45.1059767081641961095969812509 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2PowersOf4, reg, x = 1-4^-15 = 0.999999999068677425384521484375";
  }
}



TEST(AQgXspace,TruncAs2NormalDoublesRegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto eval = AQg_reg.truncate(2);

  {
    testVal  = eval(testAs,testLM,testNF,1.e-10);
    refVal   = -5.72222220977524469766928686488e10;
    errBound = 5.72222220977524469766928686488e10 * std::max(eps/(1.-1.e-10),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2NormalDoubles, reg, x = 10^-10 = 1.e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-9);
    refVal   = -5.722222181916790520678952863e9;
    errBound = 5.722222181916790520678952863e9 * std::max(eps/(1.-1.e-9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2NormalDoubles, reg, x = 10^-9 = 1.e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-8);
    refVal   = -5.72222231370111115090173980345e8;
    errBound = 5.72222231370111115090173980345e8 * std::max(eps/(1.-1.e-8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2NormalDoubles, reg, x = 10^-8 = 1.e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-7);
    refVal   = -5.72222532337944831842181862046e7;
    errBound = 5.72222532337944831842181862046e7 * std::max(eps/(1.-1.e-7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2NormalDoubles, reg, x = 10^-7 = 1.e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-6);
    refVal   = -5.72225462931213096901893427447e6;
    errBound = 5.72225462931213096901893427447e6 * std::max(eps/(1.-1.e-6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2NormalDoubles, reg, x = 10^-6 = 1.e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-5);
    refVal   = -572242.678625563843528361618877;
    errBound = 572242.678625563843528361618877 * std::max(eps/(1.-1.e-5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2NormalDoubles, reg, x = 10^-5 = 1.e-5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0001);
    refVal   = -57224.504053434965256287007706;
    errBound = 57224.504053434965256287007706 * std::max(eps/(1.-0.0001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2NormalDoubles, reg, x = 10^-4 = 0.0001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.001);
    refVal   = -5707.1920403666855672235009316;
    errBound = 5707.1920403666855672235009316 * std::max(eps/(1.-0.001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2NormalDoubles, reg, x = 10^-3 = 0.001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.01);
    refVal   = -547.281508164520849433014235602;
    errBound = 547.281508164520849433014235602 * std::max(eps/(1.-0.01),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2NormalDoubles, reg, x = 10^-2 = 0.01";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.05);
    refVal   = -88.7357576529195408883884560983;
    errBound = 88.7357576529195408883884560983 * std::max(eps/(1.-0.05),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2NormalDoubles, reg, x = 0.05 = 0.05";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.1);
    refVal   = -32.4520718577767695519165912229;
    errBound = 32.4520718577767695519165912229 * std::max(eps/(1.-0.1),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2NormalDoubles, reg, x = 0.1 = 0.1";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.15);
    refVal   = -14.2642231631371965135989272849;
    errBound = 14.2642231631371965135989272849 * std::max(eps/(1.-0.15),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2NormalDoubles, reg, x = 0.15 = 0.15";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.2);
    refVal   = -5.55775152106242571895060670763;
    errBound = 5.55775152106242571895060670763 * std::max(eps/(1.-0.2),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2NormalDoubles, reg, x = 0.2 = 0.2";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.25);
    refVal   = -0.614901849564000643729273697693;
    errBound = 0.614901849564000643729273697693 * std::max(eps/(1.-0.25),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2NormalDoubles, reg, x = 0.25 = 0.25";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.3);
    refVal   = 2.48892522612888928253052895919;
    errBound = 2.48892522612888928253052895919 * std::max(eps/(1.-0.3),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2NormalDoubles, reg, x = 0.3 = 0.3";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.35);
    refVal   = 4.60021178820168106249077911068;
    errBound = 4.60021178820168106249077911068 * std::max(eps/(1.-0.35),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2NormalDoubles, reg, x = 0.35 = 0.35";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.4);
    refVal   = 6.16275961667800946141337988692;
    errBound = 6.16275961667800946141337988692 * std::max(eps/(1.-0.4),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2NormalDoubles, reg, x = 0.4 = 0.4";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.45);
    refVal   = 7.44200074715267850724237468578;
    errBound = 7.44200074715267850724237468578 * std::max(eps/(1.-0.45),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2NormalDoubles, reg, x = 0.45 = 0.45";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.5);
    refVal   = 8.61459576259538526682894362735;
    errBound = 8.61459576259538526682894362735 * std::max(eps/(1.-0.5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2NormalDoubles, reg, x = 0.5 = 0.5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.55);
    refVal   = 9.80926269077433789371418629677;
    errBound = 9.80926269077433789371418629677 * std::max(eps/(1.-0.55),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2NormalDoubles, reg, x = 0.55 = 0.55";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.6);
    refVal   = 11.1273620249890748796840075517;
    errBound = 11.1273620249890748796840075517 * std::max(eps/(1.-0.6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2NormalDoubles, reg, x = 0.6 = 0.6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.65);
    refVal   = 12.65413054383250890323783782334;
    errBound = 12.65413054383250890323783782334 * std::max(eps/(1.-0.65),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2NormalDoubles, reg, x = 0.65 = 0.65";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.7);
    refVal   = 14.4650927822982595569302418895;
    errBound = 14.4650927822982595569302418895 * std::max(eps/(1.-0.7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2NormalDoubles, reg, x = 0.7 = 0.7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.75);
    refVal   = 16.6293673986411701977835482829;
    errBound = 16.6293673986411701977835482829 * std::max(eps/(1.-0.75),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2NormalDoubles, reg, x = 0.75 = 0.75";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.8);
    refVal   = 19.2094088069687244289978165192;
    errBound = 19.2094088069687244289978165192 * std::max(eps/(1.-0.8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2NormalDoubles, reg, x = 0.8 = 0.8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.85);
    refVal   = 22.2517905846407684899888749936;
    errBound = 22.2517905846407684899888749936 * std::max(eps/(1.-0.85),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2NormalDoubles, reg, x = 0.85 = 0.85";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9);
    refVal   = 25.7399600862304706112079569078;
    errBound = 25.7399600862304706112079569078 * std::max(eps/(1.-0.9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2NormalDoubles, reg, x = 0.9 = 0.9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.95);
    refVal   = 29.2753277691405533065530754853;
    errBound = 29.2753277691405533065530754853 * std::max(eps/(1.-0.95),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2NormalDoubles, reg, x = 0.95 = 0.95";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99);
    refVal   = 28.0080862255227553233282075876;
    errBound = 28.0080862255227553233282075876 * std::max(eps/(1.-0.99),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2NormalDoubles, reg, x = 1-10^-2 = 0.99";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999);
    refVal   = 13.00640063716497596337895173695;
    errBound = 13.00640063716497596337895173695 * std::max(eps/(1.-0.999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2NormalDoubles, reg, x = 1-10^-3 = 0.999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999);
    refVal   = -9.40326295106120627303974566481;
    errBound = 9.40326295106120627303974566481 * std::max(eps/(1.-0.9999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2NormalDoubles, reg, x = 1-10^-4 = 0.9999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999);
    refVal   = -33.5935681118423067595206883401;
    errBound = 33.5935681118423067595206883401 * std::max(eps/(1.-0.99999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2NormalDoubles, reg, x = 1-10^-5 = 0.99999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999);
    refVal   = -54.4419934149829187647596053617;
    errBound = 54.4419934149829187647596053617 * std::max(eps/(1.-0.999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2NormalDoubles, reg, x = 1-10^-6 = 0.999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999);
    refVal   = -66.8611702105113740175127970445;
    errBound = 66.8611702105113740175127970445 * std::max(eps/(1.-0.9999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2NormalDoubles, reg, x = 1-10^-7 = 0.9999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999999);
    refVal   = -65.7646539009117756331582395491;
    errBound = 65.7646539009117756331582395491 * std::max(eps/(1.-0.99999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2NormalDoubles, reg, x = 1-10^-8 = 0.99999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999);
    refVal   = -46.0658004494490074400080106041;
    errBound = 46.0658004494490074400080106041 * std::max(eps/(1.-0.999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2NormalDoubles, reg, x = 1-10^-9 = 0.999999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999999);
    refVal   = -2.67792013814202834117876981324;
    errBound = 2.67792013814202834117876981324 * std::max(eps/(1.-0.9999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQg, TruncAs2NormalDoubles, reg, x = 1-10^-10 = 0.9999999999";
  }
}



TEST(AQgNspace,FullMoments)
{
  double refVal = 0., errBound = 0.;
  auto& rpd = AQg;
  integration_engine_gsl engine;
  auto mom = make_mellin_moment(rpd, engine, testAs, testLM, testNF);
  std::tuple<integration_status, double, double> res;

  {
    SCOPED_TRACE("AQg Full Mellin moment N=2");
    refVal = 68.3552495306151645846642628695;
    errBound = 68.3552495306151645846642628695 * 1.e-9;
    res = mom.integrate(2, 0., 1.e-9);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AQg Full Mellin moment N=4");
    refVal = 61.1193765214428097015862781084;
    errBound = 61.1193765214428097015862781084 * 1.e-9;
    res = mom.integrate(4, 0., 1.e-9);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AQg Full Mellin moment N=6");
    refVal = 48.9231135876527279876057809573;
    errBound = 48.9231135876527279876057809573 * 1.e-9;
    res = mom.integrate(6, 0., 1.e-9);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AQg Full Mellin moment N=8");
    refVal = 40.5653152099417635925154297063;
    errBound = 40.5653152099417635925154297063 * 1.e-9;
    res = mom.integrate(8, 0., 1.e-9);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

}


TEST(AQgNspace,TruncAs1Moments)
{
  double refVal = 0., errBound = 0.;
  auto rpd = AQg.truncate(1);
  integration_engine_gsl engine;
  auto mom = make_mellin_moment(rpd, engine, testAs, testLM, testNF);
  std::tuple<integration_status, double, double> res;

  {
    SCOPED_TRACE("AQg TruncAs1 Mellin moment N=2");
    refVal = 0.833333333333333333333333333333;
    errBound = 0.833333333333333333333333333333 * 1.e-9;
    res = mom.integrate(2, 0., 1.e-9);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AQg TruncAs1 Mellin moment N=4");
    refVal = 0.458333333333333333333333333333;
    errBound = 0.458333333333333333333333333333 * 1.e-9;
    res = mom.integrate(4, 0., 1.e-9);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AQg TruncAs1 Mellin moment N=6");
    refVal = 0.327380952380952380952380952381;
    errBound = 0.327380952380952380952380952381 * 1.e-9;
    res = mom.integrate(6, 0., 1.e-9);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AQg TruncAs1 Mellin moment N=8");
    refVal = 0.256944444444444444444444444444;
    errBound = 0.256944444444444444444444444444 * 1.e-9;
    res = mom.integrate(8, 0., 1.e-9);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

}


TEST(AQgNspace,TruncAs2Moments)
{
  double refVal = 0., errBound = 0.;
  auto rpd = AQg.truncate(2);
  integration_engine_gsl engine;
  auto mom = make_mellin_moment(rpd, engine, testAs, testLM, testNF);
  std::tuple<integration_status, double, double> res;

  {
    SCOPED_TRACE("AQg TruncAs2 Mellin moment N=2");
    refVal = 6.99601337448559670781893004115;
    errBound = 6.99601337448559670781893004115 * 1.e-9;
    res = mom.integrate(2, 0., 1.e-9);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AQg TruncAs2 Mellin moment N=4");
    refVal = 5.10096838348765432098765432099;
    errBound = 5.10096838348765432098765432099 * 1.e-9;
    res = mom.integrate(4, 0., 1.e-9);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AQg TruncAs2 Mellin moment N=6");
    refVal = 3.8192172743176068613386397643;
    errBound = 3.8192172743176068613386397643 * 1.e-9;
    res = mom.integrate(6, 0., 1.e-9);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AQg TruncAs2 Mellin moment N=8");
    refVal = 3.04737355771343864668257182465;
    errBound = 3.04737355771343864668257182465 * 1.e-9;
    res = mom.integrate(8, 0., 1.e-9);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

}
