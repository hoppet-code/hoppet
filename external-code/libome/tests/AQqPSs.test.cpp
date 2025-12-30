/******************************************************************************
 * This file is part of libome                                                *
 * Copyright (C) 2025 Arnd Behring, Kay Schoenwald                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 ******************************************************************************/

#include <limits>
#include <gtest/gtest.h>
#include <ome/integration_engine_gsl.h>
#include <ome/mellin.h>
#include <ome/AQqPSs.h>

using namespace ome;

const double eps = std::numeric_limits<double>::epsilon();
const double testAs = 0.25;
const double testLM = -5.;
const double testNF = 3.;



TEST(AQqPSsXspace,FullPowersOf4RegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto& eval = AQqPSs_reg;

  {
    testVal  = eval(testAs,testLM,testNF,9.31322574615478515625e-10);
    refVal   = 3486.23705272201552917965513077;
    errBound = 3486.23705272201552917965513077 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullPowersOf4, reg, x = 4^-15 = 9.31322574615478515625e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.7252902984619140625e-9);
    refVal   = 3656.55227887637293069948788493;
    errBound = 3656.55227887637293069948788493 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullPowersOf4, reg, x = 4^-14 = 3.7252902984619140625e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.490116119384765625e-8);
    refVal   = 3479.04677742695656471759969968;
    errBound = 3479.04677742695656471759969968 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullPowersOf4, reg, x = 4^-13 = 1.490116119384765625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,5.9604644775390625e-8);
    refVal   = 3085.22417731407490977867582891;
    errBound = 3085.22417731407490977867582891 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullPowersOf4, reg, x = 4^-12 = 5.9604644775390625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,2.384185791015625e-7);
    refVal   = 2578.84307938487106125081556457;
    errBound = 2578.84307938487106125081556457 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullPowersOf4, reg, x = 4^-11 = 2.384185791015625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,9.5367431640625e-7);
    refVal   = 2038.76146506583109262150821504;
    errBound = 2038.76146506583109262150821504 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullPowersOf4, reg, x = 4^-10 = 9.5367431640625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.814697265625e-6);
    refVal   = 1521.7808899895146166541436794;
    errBound = 1521.7808899895146166541436794 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullPowersOf4, reg, x = 4^-9 = 3.814697265625e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0000152587890625);
    refVal   = 1065.48994648597961349832581306;
    errBound = 1065.48994648597961349832581306 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullPowersOf4, reg, x = 4^-8 = 0.0000152587890625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00006103515625);
    refVal   = 691.105348283077668373579789183;
    errBound = 691.105348283077668373579789183 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullPowersOf4, reg, x = 4^-7 = 0.00006103515625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.000244140625);
    refVal   = 406.305720704232917758556738499;
    errBound = 406.305720704232917758556738499 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullPowersOf4, reg, x = 4^-6 = 0.000244140625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0009765625);
    refVal   = 208.04483103120483467436420626;
    errBound = 208.04483103120483467436420626 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullPowersOf4, reg, x = 4^-5 = 0.0009765625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00390625);
    refVal   = 85.3145208152306968419345070931;
    errBound = 85.3145208152306968419345070931 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullPowersOf4, reg, x = 4^-4 = 0.00390625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.015625);
    refVal   = 21.816158517285804434193711371;
    errBound = 21.816158517285804434193711371 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullPowersOf4, reg, x = 4^-3 = 0.015625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0625);
    refVal   = -1.42143670137779381366114343548;
    errBound = 1.42143670137779381366114343548 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullPowersOf4, reg, x = 4^-2 = 0.0625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9375);
    refVal   = -0.0960419417730695687376165257576;
    errBound = 0.0960419417730695687376165257576 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullPowersOf4, reg, x = 1-4^-2 = 0.9375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.984375);
    refVal   = -0.0244237573942589820176777543394;
    errBound = 0.0244237573942589820176777543394 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullPowersOf4, reg, x = 1-4^-3 = 0.984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99609375);
    refVal   = -0.00614317534366446828677888019875;
    errBound = 0.00614317534366446828677888019875 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullPowersOf4, reg, x = 1-4^-4 = 0.99609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9990234375);
    refVal   = -0.00151634753296332606439043097879;
    errBound = 0.00151634753296332606439043097879 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullPowersOf4, reg, x = 1-4^-5 = 0.9990234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999755859375);
    refVal   = -0.000367523372200125314324005405923;
    errBound = 0.000367523372200125314324005405923 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullPowersOf4, reg, x = 1-4^-6 = 0.999755859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99993896484375);
    refVal   = -0.0000873825880528582509430742042541;
    errBound = 0.0000873825880528582509430742042541 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullPowersOf4, reg, x = 1-4^-7 = 0.99993896484375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999847412109375);
    refVal   = -0.0000203250454987846963799512367525;
    errBound = 0.0000203250454987846963799512367525 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullPowersOf4, reg, x = 1-4^-8 = 0.9999847412109375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999996185302734375);
    refVal   = -4.60254612804722303236226426319e-6;
    errBound = 4.60254612804722303236226426319e-6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullPowersOf4, reg, x = 1-4^-9 = 0.999996185302734375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999904632568359375);
    refVal   = -1.00634911996580551011445714114e-6;
    errBound = 1.00634911996580551011445714114e-6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullPowersOf4, reg, x = 1-4^-10 = 0.99999904632568359375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999997615814208984375);
    refVal   = -2.09365612631151778018690432447e-7;
    errBound = 2.09365612631151778018690432447e-7 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullPowersOf4, reg, x = 1-4^-11 = 0.9999997615814208984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999940395355224609375);
    refVal   = -4.02486972113381504276700135573e-8;
    errBound = 4.02486972113381504276700135573e-8 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullPowersOf4, reg, x = 1-4^-12 = 0.999999940395355224609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999998509883880615234375);
    refVal   = -6.65468697443533219362889673418e-9;
    errBound = 6.65468697443533219362889673418e-9 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullPowersOf4, reg, x = 1-4^-13 = 0.99999998509883880615234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999962747097015380859375);
    refVal   = -7.15722987889583312903394059273e-10;
    errBound = 7.15722987889583312903394059273e-10 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullPowersOf4, reg, x = 1-4^-14 = 0.9999999962747097015380859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999068677425384521484375);
    refVal   = 8.20756191764054124325567359804e-11;
    errBound = 8.20756191764054124325567359804e-11 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullPowersOf4, reg, x = 1-4^-15 = 0.999999999068677425384521484375";
  }
}



TEST(AQqPSsXspace,FullNormalDoublesRegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto& eval = AQqPSs_reg;

  {
    testVal  = eval(testAs,testLM,testNF,1.e-10);
    refVal   = 2059.14873437031496417241448767;
    errBound = 2059.14873437031496417241448767 * std::max(eps/(1.-1.e-10),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullNormalDoubles, reg, x = 10^-10 = 1.e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-9);
    refVal   = 3505.87496483559602629250743471;
    errBound = 3505.87496483559602629250743471 * std::max(eps/(1.-1.e-9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullNormalDoubles, reg, x = 10^-9 = 1.e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-8);
    refVal   = 3557.51966867236951227492325591;
    errBound = 3557.51966867236951227492325591 * std::max(eps/(1.-1.e-8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullNormalDoubles, reg, x = 10^-8 = 1.e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-7);
    refVal   = 2904.64725093152031867669198617;
    errBound = 2904.64725093152031867669198617 * std::max(eps/(1.-1.e-7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullNormalDoubles, reg, x = 10^-7 = 1.e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-6);
    refVal   = 2020.42627115720954576474884947;
    errBound = 2020.42627115720954576474884947 * std::max(eps/(1.-1.e-6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullNormalDoubles, reg, x = 10^-6 = 1.e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-5);
    refVal   = 1196.58635344536514429080046329;
    errBound = 1196.58635344536514429080046329 * std::max(eps/(1.-1.e-5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullNormalDoubles, reg, x = 10^-5 = 1.e-5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0001);
    refVal   = 579.362052714239210493881952328;
    errBound = 579.362052714239210493881952328 * std::max(eps/(1.-0.0001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullNormalDoubles, reg, x = 10^-4 = 0.0001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.001);
    refVal   = 205.35062289419745161907332294;
    errBound = 205.35062289419745161907332294 * std::max(eps/(1.-0.001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullNormalDoubles, reg, x = 10^-3 = 0.001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.01);
    refVal   = 36.926360219130694892410694679;
    errBound = 36.926360219130694892410694679 * std::max(eps/(1.-0.01),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullNormalDoubles, reg, x = 10^-2 = 0.01";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.05);
    refVal   = 0.401898767041831628726687109079;
    errBound = 0.401898767041831628726687109079 * std::max(eps/(1.-0.05),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullNormalDoubles, reg, x = 0.05 = 0.05";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.1);
    refVal   = -3.55491814060246691944104238624;
    errBound = 3.55491814060246691944104238624 * std::max(eps/(1.-0.1),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullNormalDoubles, reg, x = 0.1 = 0.1";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.15);
    refVal   = -3.95657368162054781758182678459;
    errBound = 3.95657368162054781758182678459 * std::max(eps/(1.-0.15),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullNormalDoubles, reg, x = 0.15 = 0.15";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.2);
    refVal   = -3.67170761289455942240568848069;
    errBound = 3.67170761289455942240568848069 * std::max(eps/(1.-0.2),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullNormalDoubles, reg, x = 0.2 = 0.2";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.25);
    refVal   = -3.22802247786239797824579078343;
    errBound = 3.22802247786239797824579078343 * std::max(eps/(1.-0.25),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullNormalDoubles, reg, x = 0.25 = 0.25";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.3);
    refVal   = -2.77071319268462933614761984569;
    errBound = 2.77071319268462933614761984569 * std::max(eps/(1.-0.3),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullNormalDoubles, reg, x = 0.3 = 0.3";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.35);
    refVal   = -2.3456265065859783362366947807;
    errBound = 2.3456265065859783362366947807 * std::max(eps/(1.-0.35),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullNormalDoubles, reg, x = 0.35 = 0.35";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.4);
    refVal   = -1.96641078398449258860880278127;
    errBound = 1.96641078398449258860880278127 * std::max(eps/(1.-0.4),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullNormalDoubles, reg, x = 0.4 = 0.4";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.45);
    refVal   = -1.63498902387747561531673825035;
    errBound = 1.63498902387747561531673825035 * std::max(eps/(1.-0.45),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullNormalDoubles, reg, x = 0.45 = 0.45";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.5);
    refVal   = -1.34881621425471613654390759769;
    errBound = 1.34881621425471613654390759769 * std::max(eps/(1.-0.5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullNormalDoubles, reg, x = 0.5 = 0.5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.55);
    refVal   = -1.10370226579823535227744678105;
    errBound = 1.10370226579823535227744678105 * std::max(eps/(1.-0.55),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullNormalDoubles, reg, x = 0.55 = 0.55";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.6);
    refVal   = -0.894965095215363438725143610922;
    errBound = 0.894965095215363438725143610922 * std::max(eps/(1.-0.6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullNormalDoubles, reg, x = 0.6 = 0.6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.65);
    refVal   = -0.717901807532403680623623829585;
    errBound = 0.717901807532403680623623829585 * std::max(eps/(1.-0.65),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullNormalDoubles, reg, x = 0.65 = 0.65";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.7);
    refVal   = -0.567964554410041843440562017573;
    errBound = 0.567964554410041843440562017573 * std::max(eps/(1.-0.7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullNormalDoubles, reg, x = 0.7 = 0.7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.75);
    refVal   = -0.440801739841483120788042110649;
    errBound = 0.440801739841483120788042110649 * std::max(eps/(1.-0.75),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullNormalDoubles, reg, x = 0.75 = 0.75";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.8);
    refVal   = -0.332232217527925274311430445765;
    errBound = 0.332232217527925274311430445765 * std::max(eps/(1.-0.8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullNormalDoubles, reg, x = 0.8 = 0.8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.85);
    refVal   = -0.238176218723743584762457641756;
    errBound = 0.238176218723743584762457641756 * std::max(eps/(1.-0.85),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullNormalDoubles, reg, x = 0.85 = 0.85";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9);
    refVal   = -0.154535365075644246085897295606;
    errBound = 0.154535365075644246085897295606 * std::max(eps/(1.-0.9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullNormalDoubles, reg, x = 0.9 = 0.9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.95);
    refVal   = -0.0769521068054244942801449115057;
    errBound = 0.0769521068054244942801449115057 * std::max(eps/(1.-0.95),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullNormalDoubles, reg, x = 0.95 = 0.95";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99);
    refVal   = -0.0156944612724226609376929944652;
    errBound = 0.0156944612724226609376929944652 * std::max(eps/(1.-0.99),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullNormalDoubles, reg, x = 1-10^-2 = 0.99";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999);
    refVal   = -0.00155331620638077004069559366481;
    errBound = 0.00155331620638077004069559366481 * std::max(eps/(1.-0.999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullNormalDoubles, reg, x = 1-10^-3 = 0.999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999);
    refVal   = -0.000146091075709798894812044132522;
    errBound = 0.000146091075709798894812044132522 * std::max(eps/(1.-0.9999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullNormalDoubles, reg, x = 1-10^-4 = 0.9999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999);
    refVal   = -0.0000129650480492973207275146107826;
    errBound = 0.0000129650480492973207275146107826 * std::max(eps/(1.-0.99999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullNormalDoubles, reg, x = 1-10^-5 = 0.99999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999);
    refVal   = -1.06083645773287777342150773521e-6;
    errBound = 1.06083645773287777342150773521e-6 * std::max(eps/(1.-0.999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullNormalDoubles, reg, x = 1-10^-6 = 0.999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999);
    refVal   = -7.5400379847687458228953232959e-8;
    errBound = 7.5400379847687458228953232959e-8 * std::max(eps/(1.-0.9999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullNormalDoubles, reg, x = 1-10^-7 = 0.9999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999999);
    refVal   = -3.7601897759656575531668912318e-9;
    errBound = 3.7601897759656575531668912318e-9 * std::max(eps/(1.-0.99999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullNormalDoubles, reg, x = 1-10^-8 = 0.99999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999);
    refVal   = 7.31165566152502410534788071481e-11;
    errBound = 7.31165566152502410534788071481e-11 * std::max(eps/(1.-0.999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullNormalDoubles, reg, x = 1-10^-9 = 0.999999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999999);
    refVal   = 5.93402636014515517361614468235e-11;
    errBound = 5.93402636014515517361614468235e-11 * std::max(eps/(1.-0.9999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPSs, FullNormalDoubles, reg, x = 1-10^-10 = 0.9999999999";
  }
}



TEST(AQqPSsNspace,FullMoments)
{
  double refVal = 0., errBound = 0.;
  auto& rpd = AQqPSs;
  integration_engine_gsl engine;
  auto mom = make_mellin_moment(rpd, engine, testAs, testLM, testNF);
  std::tuple<integration_status, double, double> res;

  {
    SCOPED_TRACE("AQqPSs Full Mellin moment N=3");
    refVal = -0.203066980704295765918208800682;
    errBound = 0.203066980704295765918208800682 * 3.e-11;
    res = mom.integrate(3, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AQqPSs Full Mellin moment N=5");
    refVal = -0.0666094084540538534010648362852;
    errBound = 0.0666094084540538534010648362852 * 3.e-11;
    res = mom.integrate(5, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AQqPSs Full Mellin moment N=7");
    refVal = -0.0321997482513328566195183392847;
    errBound = 0.0321997482513328566195183392847 * 3.e-11;
    res = mom.integrate(7, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AQqPSs Full Mellin moment N=9");
    refVal = -0.0189417536645648204346631769842;
    errBound = 0.0189417536645648204346631769842 * 3.e-11;
    res = mom.integrate(9, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

}
