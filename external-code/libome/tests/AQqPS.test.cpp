/******************************************************************************
 * This file is part of libome                                                *
 * Copyright (C) 2025 Arnd Behring, Kay Schoenwald                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 ******************************************************************************/

#include <limits>
#include <gtest/gtest.h>
#include <ome/integration_engine_gsl.h>
#include <ome/mellin.h>
#include <ome/AQqPS.h>

using namespace ome;

const double eps = std::numeric_limits<double>::epsilon();
const double testAs = 0.25;
const double testLM = -5.;
const double testNF = 3.;



TEST(AQqPSXspace,FullPowersOf4RegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto& eval = AQqPS_reg;

  {
    testVal  = eval(testAs,testLM,testNF,9.31322574615478515625e-10);
    refVal   = 2.71385701515686759164116171416e11;
    errBound = 2.71385701515686759164116171416e11 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullPowersOf4, reg, x = 4^-15 = 9.31322574615478515625e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.7252902984619140625e-9);
    refVal   = 6.18444173017536489608775334601e10;
    errBound = 6.18444173017536489608775334601e10 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullPowersOf4, reg, x = 4^-14 = 3.7252902984619140625e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.490116119384765625e-8);
    refVal   = 1.3960601816105550880240348884e10;
    errBound = 1.3960601816105550880240348884e10 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullPowersOf4, reg, x = 4^-13 = 1.490116119384765625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,5.9604644775390625e-8);
    refVal   = 3.11502509157568402873069057516e9;
    errBound = 3.11502509157568402873069057516e9 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullPowersOf4, reg, x = 4^-12 = 5.9604644775390625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,2.384185791015625e-7);
    refVal   = 6.84975526435125406265832286497e8;
    errBound = 6.84975526435125406265832286497e8 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullPowersOf4, reg, x = 4^-11 = 2.384185791015625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,9.5367431640625e-7);
    refVal   = 1.47799358808982048092346588171e8;
    errBound = 1.47799358808982048092346588171e8 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullPowersOf4, reg, x = 4^-10 = 9.5367431640625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.814697265625e-6);
    refVal   = 3.10893044315240882393331370976e7;
    errBound = 3.10893044315240882393331370976e7 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullPowersOf4, reg, x = 4^-9 = 3.814697265625e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0000152587890625);
    refVal   = 6.30766476776470287484978190107e6;
    errBound = 6.30766476776470287484978190107e6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullPowersOf4, reg, x = 4^-8 = 0.0000152587890625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00006103515625);
    refVal   = 1.21109760451343802763121016049e6;
    errBound = 1.21109760451343802763121016049e6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullPowersOf4, reg, x = 4^-7 = 0.00006103515625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.000244140625);
    refVal   = 211564.470800420706289318165785;
    errBound = 211564.470800420706289318165785 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullPowersOf4, reg, x = 4^-6 = 0.000244140625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0009765625);
    refVal   = 30262.5561682196259990832457285;
    errBound = 30262.5561682196259990832457285 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullPowersOf4, reg, x = 4^-5 = 0.0009765625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00390625);
    refVal   = 2037.78142389761271712358324824;
    errBound = 2037.78142389761271712358324824 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullPowersOf4, reg, x = 4^-4 = 0.00390625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.015625);
    refVal   = -773.256840074113909316653215157;
    errBound = 773.256840074113909316653215157 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullPowersOf4, reg, x = 4^-3 = 0.015625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0625);
    refVal   = -442.145967140465503706720765883;
    errBound = 442.145967140465503706720765883 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullPowersOf4, reg, x = 4^-2 = 0.0625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9375);
    refVal   = -4.42629021000634080573005433033;
    errBound = 4.42629021000634080573005433033 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullPowersOf4, reg, x = 1-4^-2 = 0.9375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.984375);
    refVal   = -0.846141609783329533399158380965;
    errBound = 0.846141609783329533399158380965 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullPowersOf4, reg, x = 1-4^-3 = 0.984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99609375);
    refVal   = -0.1248614255154691595046134237003;
    errBound = 0.1248614255154691595046134237003 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullPowersOf4, reg, x = 1-4^-4 = 0.99609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9990234375);
    refVal   = -0.00328118869839263875792184369731;
    errBound = 0.00328118869839263875792184369731 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullPowersOf4, reg, x = 1-4^-5 = 0.9990234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999755859375);
    refVal   = 0.00785297062433552698644320908546;
    errBound = 0.00785297062433552698644320908546 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullPowersOf4, reg, x = 1-4^-6 = 0.999755859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99993896484375);
    refVal   = 0.00455125018627822894073531101025;
    errBound = 0.00455125018627822894073531101025 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullPowersOf4, reg, x = 1-4^-7 = 0.99993896484375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999847412109375);
    refVal   = 0.00187924679867558753245349867238;
    errBound = 0.00187924679867558753245349867238 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullPowersOf4, reg, x = 1-4^-8 = 0.9999847412109375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999996185302734375);
    refVal   = 0.000673591344372366876931615100425;
    errBound = 0.000673591344372366876931615100425 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullPowersOf4, reg, x = 1-4^-9 = 0.999996185302734375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999904632568359375);
    refVal   = 0.000222002037335632192788431017081;
    errBound = 0.000222002037335632192788431017081 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullPowersOf4, reg, x = 1-4^-10 = 0.99999904632568359375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999997615814208984375);
    refVal   = 0.0000689170229373789256455471567594;
    errBound = 0.0000689170229373789256455471567594 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullPowersOf4, reg, x = 1-4^-11 = 0.9999997615814208984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999940395355224609375);
    refVal   = 0.0000203841567205084135920553031147;
    errBound = 0.0000203841567205084135920553031147 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullPowersOf4, reg, x = 1-4^-12 = 0.999999940395355224609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999998509883880615234375);
    refVal   = 5.77399342322360264861932314276e-6;
    errBound = 5.77399342322360264861932314276e-6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullPowersOf4, reg, x = 1-4^-13 = 0.99999998509883880615234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999962747097015380859375);
    refVal   = 1.56750525391745770784390600046e-6;
    errBound = 1.56750525391745770784390600046e-6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullPowersOf4, reg, x = 1-4^-14 = 0.9999999962747097015380859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999068677425384521484375);
    refVal   = 4.0642341488971257579715841135e-7;
    errBound = 4.0642341488971257579715841135e-7 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullPowersOf4, reg, x = 1-4^-15 = 0.999999999068677425384521484375";
  }
}



TEST(AQqPSXspace,FullNormalDoublesRegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto& eval = AQqPS_reg;

  {
    testVal  = eval(testAs,testLM,testNF,1.e-10);
    refVal   = 2.88737926298890757288466568624e12;
    errBound = 2.88737926298890757288466568624e12 * std::max(eps/(1.-1.e-10),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullNormalDoubles, reg, x = 10^-10 = 1.e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-9);
    refVal   = 2.51600075261730014297228433584e11;
    errBound = 2.51600075261730014297228433584e11 * std::max(eps/(1.-1.e-9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullNormalDoubles, reg, x = 10^-9 = 1.e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-8);
    refVal   = 2.14462205493710998073282758433e10;
    errBound = 2.14462205493710998073282758433e10 * std::max(eps/(1.-1.e-8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullNormalDoubles, reg, x = 10^-8 = 1.e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-7);
    refVal   = 1.77324372693659428322046368556e9;
    errBound = 1.77324372693659428322046368556e9 * std::max(eps/(1.-1.e-7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullNormalDoubles, reg, x = 10^-7 = 1.e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-6);
    refVal   = 1.40187467490248758402253223961e8;
    errBound = 1.40187467490248758402253223961e8 * std::max(eps/(1.-1.e-6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullNormalDoubles, reg, x = 10^-6 = 1.e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-5);
    refVal   = 1.03058418820706654567037739368e7;
    errBound = 1.03058418820706654567037739368e7 * std::max(eps/(1.-1.e-5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullNormalDoubles, reg, x = 10^-5 = 1.e-5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0001);
    refVal   = 659786.684980947234309442994079;
    errBound = 659786.684980947234309442994079 * std::max(eps/(1.-0.0001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullNormalDoubles, reg, x = 10^-4 = 0.0001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.001);
    refVal   = 29178.0749554457439940272792459;
    errBound = 29178.0749554457439940272792459 * std::max(eps/(1.-0.001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullNormalDoubles, reg, x = 10^-3 = 0.001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.01);
    refVal   = -592.247132539866852431227179464;
    errBound = 592.247132539866852431227179464 * std::max(eps/(1.-0.01),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullNormalDoubles, reg, x = 10^-2 = 0.01";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.05);
    refVal   = -512.951956996595622575720427512;
    errBound = 512.951956996595622575720427512 * std::max(eps/(1.-0.05),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullNormalDoubles, reg, x = 0.05 = 0.05";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.1);
    refVal   = -308.964737327568383465605933051;
    errBound = 308.964737327568383465605933051 * std::max(eps/(1.-0.1),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullNormalDoubles, reg, x = 0.1 = 0.1";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.15);
    refVal   = -216.655798976572424612557503065;
    errBound = 216.655798976572424612557503065 * std::max(eps/(1.-0.15),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullNormalDoubles, reg, x = 0.15 = 0.15";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.2);
    refVal   = -164.035442242361660313472648484;
    errBound = 164.035442242361660313472648484 * std::max(eps/(1.-0.2),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullNormalDoubles, reg, x = 0.2 = 0.2";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.25);
    refVal   = -129.920215441476281120574212085;
    errBound = 129.920215441476281120574212085 * std::max(eps/(1.-0.25),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullNormalDoubles, reg, x = 0.25 = 0.25";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.3);
    refVal   = -105.909911002571665922474405673;
    errBound = 105.909911002571665922474405673 * std::max(eps/(1.-0.3),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullNormalDoubles, reg, x = 0.3 = 0.3";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.35);
    refVal   = -88.0030266596137485790000790001;
    errBound = 88.0030266596137485790000790001 * std::max(eps/(1.-0.35),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullNormalDoubles, reg, x = 0.35 = 0.35";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.4);
    refVal   = -74.0473618090319688799534433636;
    errBound = 74.0473618090319688799534433636 * std::max(eps/(1.-0.4),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullNormalDoubles, reg, x = 0.4 = 0.4";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.45);
    refVal   = -62.7803717343560345227889899991;
    errBound = 62.7803717343560345227889899991 * std::max(eps/(1.-0.45),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullNormalDoubles, reg, x = 0.45 = 0.45";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.5);
    refVal   = -53.4111590352759380297240652135;
    errBound = 53.4111590352759380297240652135 * std::max(eps/(1.-0.5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullNormalDoubles, reg, x = 0.5 = 0.5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.55);
    refVal   = -45.418021493996974162369003844;
    errBound = 45.418021493996974162369003844 * std::max(eps/(1.-0.55),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullNormalDoubles, reg, x = 0.55 = 0.55";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.6);
    refVal   = -38.4422005129060459680192286671;
    errBound = 38.4422005129060459680192286671 * std::max(eps/(1.-0.6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullNormalDoubles, reg, x = 0.6 = 0.6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.65);
    refVal   = -32.2285585076793237136124747393;
    errBound = 32.2285585076793237136124747393 * std::max(eps/(1.-0.65),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullNormalDoubles, reg, x = 0.65 = 0.65";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.7);
    refVal   = -26.5908991505905871203157248243;
    errBound = 26.5908991505905871203157248243 * std::max(eps/(1.-0.7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullNormalDoubles, reg, x = 0.7 = 0.7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.75);
    refVal   = -21.3911684126213277502976633743;
    errBound = 21.3911684126213277502976633743 * std::max(eps/(1.-0.75),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullNormalDoubles, reg, x = 0.75 = 0.75";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.8);
    refVal   = -16.5272506100638183170028435512;
    errBound = 16.5272506100638183170028435512 * std::max(eps/(1.-0.8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullNormalDoubles, reg, x = 0.8 = 0.8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.85);
    refVal   = -11.92730617088141199672441700543;
    errBound = 11.92730617088141199672441700543 * std::max(eps/(1.-0.85),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullNormalDoubles, reg, x = 0.85 = 0.85";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9);
    refVal   = -7.55235656007184082201368892365;
    errBound = 7.55235656007184082201368892365 * std::max(eps/(1.-0.9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullNormalDoubles, reg, x = 0.9 = 0.9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.95);
    refVal   = -3.42282117674308679281933388634;
    errBound = 3.42282117674308679281933388634 * std::max(eps/(1.-0.95),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullNormalDoubles, reg, x = 0.95 = 0.95";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99);
    refVal   = -0.476624352663855956577735478686;
    errBound = 0.476624352663855956577735478686 * std::max(eps/(1.-0.99),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullNormalDoubles, reg, x = 1-10^-2 = 0.99";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999);
    refVal   = -0.00390674600924077426702700773085;
    errBound = 0.00390674600924077426702700773085 * std::max(eps/(1.-0.999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullNormalDoubles, reg, x = 1-10^-3 = 0.999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999);
    refVal   = 0.00587071509289112332791584471489;
    errBound = 0.00587071509289112332791584471489 * std::max(eps/(1.-0.9999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullNormalDoubles, reg, x = 1-10^-4 = 0.9999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999);
    refVal   = 0.001390107597031348997944448834245;
    errBound = 0.001390107597031348997944448834245 * std::max(eps/(1.-0.99999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullNormalDoubles, reg, x = 1-10^-5 = 0.99999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999);
    refVal   = 0.000230844309924262939736675864178;
    errBound = 0.000230844309924262939736675864178 * std::max(eps/(1.-0.999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullNormalDoubles, reg, x = 1-10^-6 = 0.999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999);
    refVal   = 0.0000322868577287302020002292689042;
    errBound = 0.0000322868577287302020002292689042 * std::max(eps/(1.-0.9999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullNormalDoubles, reg, x = 1-10^-7 = 0.9999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999999);
    refVal   = 3.9854170324192378564245263039e-6;
    errBound = 3.9854170324192378564245263039e-6 * std::max(eps/(1.-0.99999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullNormalDoubles, reg, x = 1-10^-8 = 0.99999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999);
    refVal   = 4.36116576329611674822372754252e-7;
    errBound = 4.36116576329611674822372754252e-7 * std::max(eps/(1.-0.999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullNormalDoubles, reg, x = 1-10^-9 = 0.999999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999999);
    refVal   = 4.0864054307351879149683013165e-8;
    errBound = 4.0864054307351879149683013165e-8 * std::max(eps/(1.-0.9999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, FullNormalDoubles, reg, x = 1-10^-10 = 0.9999999999";
  }
}



TEST(AQqPSXspace,TruncAs2PowersOf4RegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto eval = AQqPS_reg.truncate(2);

  {
    testVal  = eval(testAs,testLM,testNF,9.31322574615478515625e-10);
    refVal   = -2.73075085759923740929447876764e9;
    errBound = 2.73075085759923740929447876764e9 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2PowersOf4, reg, x = 4^-15 = 9.31322574615478515625e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.7252902984619140625e-9);
    refVal   = -6.82687763714267322017310017118e8;
    errBound = 6.82687763714267322017310017118e8 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2PowersOf4, reg, x = 4^-14 = 3.7252902984619140625e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.490116119384765625e-8);
    refVal   = -1.70671995037671282978441737946e8;
    errBound = 1.70671995037671282978441737946e8 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2PowersOf4, reg, x = 4^-13 = 1.490116119384765625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,5.9604644775390625e-8);
    refVal   = -4.26680537908505025509561191424e7;
    errBound = 4.26680537908505025509561191424e7 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2PowersOf4, reg, x = 4^-12 = 5.9604644775390625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,2.384185791015625e-7);
    refVal   = -1.0667066195212400323850578714e7;
    errBound = 1.0667066195212400323850578714e7 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2PowersOf4, reg, x = 4^-11 = 2.384185791015625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,9.5367431640625e-7);
    refVal   = -2.66681447218697482402758671092e6;
    errBound = 2.66681447218697482402758671092e6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2PowersOf4, reg, x = 4^-10 = 9.5367431640625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.814697265625e-6);
    refVal   = -666744.843282018028104357452422;
    errBound = 666744.843282018028104357452422 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2PowersOf4, reg, x = 4^-9 = 3.814697265625e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0000152587890625);
    refVal   = -166719.530256732784469226563468;
    errBound = 166719.530256732784469226563468 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2PowersOf4, reg, x = 4^-8 = 0.0000152587890625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00006103515625);
    refVal   = -41704.7556153880388561672154811;
    errBound = 41704.7556153880388561672154811 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2PowersOf4, reg, x = 4^-7 = 0.00006103515625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.000244140625);
    refVal   = -10442.7437830101373565613414705;
    errBound = 10442.7437830101373565613414705 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2PowersOf4, reg, x = 4^-6 = 0.000244140625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0009765625);
    refVal   = -2619.72283986132822341968293366;
    errBound = 2619.72283986132822341968293366 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2PowersOf4, reg, x = 4^-5 = 0.0009765625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00390625);
    refVal   = -657.92138019003881227456031795;
    errBound = 657.92138019003881227456031795 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2PowersOf4, reg, x = 4^-4 = 0.00390625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.015625);
    refVal   = -163.526883133583988149005306978;
    errBound = 163.526883133583988149005306978 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2PowersOf4, reg, x = 4^-3 = 0.015625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0625);
    refVal   = -38.4695051614245119010817195546;
    errBound = 38.4695051614245119010817195546 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2PowersOf4, reg, x = 4^-2 = 0.0625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9375);
    refVal   = -0.28441416593558800666261457494;
    errBound = 0.28441416593558800666261457494 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2PowersOf4, reg, x = 1-4^-2 = 0.9375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.984375);
    refVal   = -0.0710386562938892088027443897225;
    errBound = 0.0710386562938892088027443897225 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2PowersOf4, reg, x = 1-4^-3 = 0.984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99609375);
    refVal   = -0.0177644836510657961296073800316;
    errBound = 0.0177644836510657961296073800316 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2PowersOf4, reg, x = 1-4^-4 = 0.99609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9990234375);
    refVal   = -0.00444155425532261867392595662161;
    errBound = 0.00444155425532261867392595662161 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2PowersOf4, reg, x = 1-4^-5 = 0.9990234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999755859375);
    refVal   = -0.00111041768812273395625187851773;
    errBound = 0.00111041768812273395625187851773 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2PowersOf4, reg, x = 1-4^-6 = 0.999755859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99993896484375);
    refVal   = -0.000277606274087730046440618508022;
    errBound = 0.000277606274087730046440618508022 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2PowersOf4, reg, x = 1-4^-7 = 0.99993896484375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999847412109375);
    refVal   = -0.0000694016847718391216034818248114;
    errBound = 0.0000694016847718391216034818248114 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2PowersOf4, reg, x = 1-4^-8 = 0.9999847412109375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999996185302734375);
    refVal   = -0.0000173504284663328751076139651337;
    errBound = 0.0000173504284663328751076139651337 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2PowersOf4, reg, x = 1-4^-9 = 0.999996185302734375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999904632568359375);
    refVal   = -4.33760757129018691092945943116e-6;
    errBound = 4.33760757129018691092945943116e-6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2PowersOf4, reg, x = 1-4^-10 = 0.99999904632568359375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999997615814208984375);
    refVal   = -1.08440192124362517905288429035e-6;
    errBound = 1.08440192124362517905288429035e-6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2PowersOf4, reg, x = 1-4^-11 = 0.9999997615814208984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999940395355224609375);
    refVal   = -2.71100482087253275120852009467e-7;
    errBound = 2.71100482087253275120852009467e-7 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2PowersOf4, reg, x = 1-4^-12 = 0.999999940395355224609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999998509883880615234375);
    refVal   = -6.77751206328354671952094449088e-8;
    errBound = 6.77751206328354671952094449088e-8 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2PowersOf4, reg, x = 1-4^-13 = 0.99999998509883880615234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999962747097015380859375);
    refVal   = -1.69437801651477582957171867635e-8;
    errBound = 1.69437801651477582957171867635e-8 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2PowersOf4, reg, x = 1-4^-14 = 0.9999999962747097015380859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999068677425384521484375);
    refVal   = -4.23594504172062040531424250148e-9;
    errBound = 4.23594504172062040531424250148e-9 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2PowersOf4, reg, x = 1-4^-15 = 0.999999999068677425384521484375";
  }
}



TEST(AQqPSXspace,TruncAs2NormalDoublesRegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto eval = AQqPS_reg.truncate(2);

  {
    testVal  = eval(testAs,testLM,testNF,1.e-10);
    refVal   = -2.54320987718894031101145602826e10;
    errBound = 2.54320987718894031101145602826e10 * std::max(eps/(1.-1.e-10),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2NormalDoubles, reg, x = 10^-10 = 1.e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-9);
    refVal   = -2.54320992341914125203291162277e9;
    errBound = 2.54320992341914125203291162277e9 * std::max(eps/(1.-1.e-9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2NormalDoubles, reg, x = 10^-9 = 1.e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-8);
    refVal   = -2.54321055114009332472295297621e8;
    errBound = 2.54321055114009332472295297621e8 * std::max(eps/(1.-1.e-8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2NormalDoubles, reg, x = 10^-8 = 1.e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-7);
    refVal   = -2.54321710433706742657419105961e7;
    errBound = 2.54321710433706742657419105961e7 * std::max(eps/(1.-1.e-7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2NormalDoubles, reg, x = 10^-7 = 1.e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-6);
    refVal   = -2.54327527664023787125621637603e6;
    errBound = 2.54327527664023787125621637603e6 * std::max(eps/(1.-1.e-6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2NormalDoubles, reg, x = 10^-6 = 1.e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-5);
    refVal   = -254371.883645423598454705924501;
    errBound = 254371.883645423598454705924501 * std::max(eps/(1.-1.e-5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2NormalDoubles, reg, x = 10^-5 = 1.e-5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0001);
    refVal   = -25464.9369521361850641145438301;
    errBound = 25464.9369521361850641145438301 * std::max(eps/(1.-0.0001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2NormalDoubles, reg, x = 10^-4 = 0.0001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.001);
    refVal   = -2558.52186767614735255136853727;
    errBound = 2558.52186767614735255136853727 * std::max(eps/(1.-0.001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2NormalDoubles, reg, x = 10^-3 = 0.001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.01);
    refVal   = -256.724434112796770084081670136;
    errBound = 256.724434112796770084081670136 * std::max(eps/(1.-0.01),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2NormalDoubles, reg, x = 10^-2 = 0.01";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.05);
    refVal   = -48.9078512208815899859448735308;
    errBound = 48.9078512208815899859448735308 * std::max(eps/(1.-0.05),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2NormalDoubles, reg, x = 0.05 = 0.05";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.1);
    refVal   = -22.8930872297648058883893355704;
    errBound = 22.8930872297648058883893355704 * std::max(eps/(1.-0.1),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2NormalDoubles, reg, x = 0.1 = 0.1";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.15);
    refVal   = -14.3669862816390225672719931519;
    errBound = 14.3669862816390225672719931519 * std::max(eps/(1.-0.15),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2NormalDoubles, reg, x = 0.15 = 0.15";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.2);
    refVal   = -10.1880495866015431374017173232;
    errBound = 10.1880495866015431374017173232 * std::max(eps/(1.-0.2),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2NormalDoubles, reg, x = 0.2 = 0.2";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.25);
    refVal   = -7.72668898096843758143817098164;
    errBound = 7.72668898096843758143817098164 * std::max(eps/(1.-0.25),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2NormalDoubles, reg, x = 0.25 = 0.25";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.3);
    refVal   = -6.11049327874046892082381924689;
    errBound = 6.11049327874046892082381924689 * std::max(eps/(1.-0.3),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2NormalDoubles, reg, x = 0.3 = 0.3";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.35);
    refVal   = -4.96830265274403074337901843211;
    errBound = 4.96830265274403074337901843211 * std::max(eps/(1.-0.35),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2NormalDoubles, reg, x = 0.35 = 0.35";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.4);
    refVal   = -4.11611162989638330261430552401;
    errBound = 4.11611162989638330261430552401 * std::max(eps/(1.-0.4),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2NormalDoubles, reg, x = 0.4 = 0.4";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.45);
    refVal   = -3.45257176693040391587682522813;
    errBound = 3.45257176693040391587682522813 * std::max(eps/(1.-0.45),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2NormalDoubles, reg, x = 0.45 = 0.45";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.5);
    refVal   = -2.91734359098792785988671321303;
    errBound = 2.91734359098792785988671321303 * std::max(eps/(1.-0.5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2NormalDoubles, reg, x = 0.5 = 0.5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.55);
    refVal   = -2.47225918818199503414092128709;
    errBound = 2.47225918818199503414092128709 * std::max(eps/(1.-0.55),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2NormalDoubles, reg, x = 0.55 = 0.55";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.6);
    refVal   = -2.09196027763492453368301192062;
    errBound = 2.09196027763492453368301192062 * std::max(eps/(1.-0.6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2NormalDoubles, reg, x = 0.6 = 0.6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.65);
    refVal   = -1.75889211205602323897623259085;
    errBound = 1.75889211205602323897623259085 * std::max(eps/(1.-0.65),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2NormalDoubles, reg, x = 0.65 = 0.65";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.7);
    refVal   = -1.46046482020242964132989363019;
    errBound = 1.46046482020242964132989363019 * std::max(eps/(1.-0.7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2NormalDoubles, reg, x = 0.7 = 0.7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.75);
    refVal   = -1.18736434645573117666735526113;
    errBound = 1.18736434645573117666735526113 * std::max(eps/(1.-0.75),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2NormalDoubles, reg, x = 0.75 = 0.75";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.8);
    refVal   = -0.93250607685227720988901246537;
    errBound = 0.93250607685227720988901246537 * std::max(eps/(1.-0.8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2NormalDoubles, reg, x = 0.8 = 0.8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.85);
    refVal   = -0.690363996393134263929669610266;
    errBound = 0.690363996393134263929669610266 * std::max(eps/(1.-0.85),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2NormalDoubles, reg, x = 0.85 = 0.85";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9);
    refVal   = -0.456527700426807299742447235963;
    errBound = 0.456527700426807299742447235963 * std::max(eps/(1.-0.9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2NormalDoubles, reg, x = 0.9 = 0.9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.95);
    refVal   = -0.227402228487268293507606910589;
    errBound = 0.227402228487268293507606910589 * std::max(eps/(1.-0.95),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2NormalDoubles, reg, x = 0.95 = 0.95";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99);
    refVal   = -0.0454695743018762281966931905853;
    errBound = 0.0454695743018762281966931905853 * std::max(eps/(1.-0.99),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2NormalDoubles, reg, x = 1-10^-2 = 0.99";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999);
    refVal   = -0.0045481477951493994370184128976;
    errBound = 0.0045481477951493994370184128976 * std::max(eps/(1.-0.999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2NormalDoubles, reg, x = 1-10^-3 = 0.999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999);
    refVal   = -0.000454829471999040887848664122157;
    errBound = 0.000454829471999040887848664122157 * std::max(eps/(1.-0.9999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2NormalDoubles, reg, x = 1-10^-4 = 0.9999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999);
    refVal   = -0.0000454830968925979783401698064308;
    errBound = 0.0000454830968925979783401698064308 * std::max(eps/(1.-0.99999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2NormalDoubles, reg, x = 1-10^-5 = 0.99999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999);
    refVal   = -4.54831118895251268139936447453e-6;
    errBound = 4.54831118895251268139936447453e-6 * std::max(eps/(1.-0.999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2NormalDoubles, reg, x = 1-10^-6 = 0.999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999);
    refVal   = -4.54831133894943985082891415752e-7;
    errBound = 4.54831133894943985082891415752e-7 * std::max(eps/(1.-0.9999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2NormalDoubles, reg, x = 1-10^-7 = 0.9999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999999);
    refVal   = -4.5483113539494091225441651701e-8;
    errBound = 4.5483113539494091225441651701e-8 * std::max(eps/(1.-0.99999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2NormalDoubles, reg, x = 1-10^-8 = 0.99999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999);
    refVal   = -4.54831135544940881526133863572e-9;
    errBound = 4.54831135544940881526133863572e-9 * std::max(eps/(1.-0.999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2NormalDoubles, reg, x = 1-10^-9 = 0.999999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999999);
    refVal   = -4.54831135559940881218851039133e-10;
    errBound = 4.54831135559940881218851039133e-10 * std::max(eps/(1.-0.9999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AQqPS, TruncAs2NormalDoubles, reg, x = 1-10^-10 = 0.9999999999";
  }
}



TEST(AQqPSNspace,FullMoments)
{
  double refVal = 0., errBound = 0.;
  auto& rpd = AQqPS;
  integration_engine_gsl engine;
  auto mom = make_mellin_moment(rpd, engine, testAs, testLM, testNF);
  std::tuple<integration_status, double, double> res;

  {
    SCOPED_TRACE("AQqPS Full Mellin moment N=2");
    refVal = -21.9987970416932404460059955113;
    errBound = 21.9987970416932404460059955113 * 3.e-11;
    res = mom.integrate(2, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AQqPS Full Mellin moment N=4");
    refVal = -4.77399482861080891205508206879;
    errBound = 4.77399482861080891205508206879 * 3.e-11;
    res = mom.integrate(4, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AQqPS Full Mellin moment N=6");
    refVal = -2.05321796632440898501494923712;
    errBound = 2.05321796632440898501494923712 * 3.e-11;
    res = mom.integrate(6, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AQqPS Full Mellin moment N=8");
    refVal = -1.13668634294632596834371488214;
    errBound = 1.13668634294632596834371488214 * 3.e-11;
    res = mom.integrate(8, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

}


TEST(AQqPSNspace,TruncAs2Moments)
{
  double refVal = 0., errBound = 0.;
  auto rpd = AQqPS.truncate(2);
  integration_engine_gsl engine;
  auto mom = make_mellin_moment(rpd, engine, testAs, testLM, testNF);
  std::tuple<integration_status, double, double> res;

  {
    SCOPED_TRACE("AQqPS TruncAs2 Mellin moment N=2");
    refVal = -1.4012345679012345679012345679;
    errBound = 1.4012345679012345679012345679 * 3.e-11;
    res = mom.integrate(2, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AQqPS TruncAs2 Mellin moment N=4");
    refVal = -0.269832484567901234567901234568;
    errBound = 0.269832484567901234567901234568 * 3.e-11;
    res = mom.integrate(4, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AQqPS TruncAs2 Mellin moment N=6");
    refVal = -0.117241578654641498826791991677;
    errBound = 0.117241578654641498826791991677 * 3.e-11;
    res = mom.integrate(6, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AQqPS TruncAs2 Mellin moment N=8");
    refVal = -0.0660837779589451072087336469633;
    errBound = 0.0660837779589451072087336469633 * 3.e-11;
    res = mom.integrate(8, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

}
