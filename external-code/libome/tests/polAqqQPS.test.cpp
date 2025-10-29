/******************************************************************************
 * This file is part of libome                                                *
 * Copyright (C) 2025 Arnd Behring, Kay Schoenwald                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 ******************************************************************************/

#include <limits>
#include <gtest/gtest.h>
#include <ome/integration_engine_gsl.h>
#include <ome/mellin.h>
#include <ome/polAqqQPS.h>

using namespace ome;

const double eps = std::numeric_limits<double>::epsilon();
const double testAs = 0.25;
const double testLM = -5.;
const double testNF = 3.;



TEST(polAqqQPSXspace,FullPowersOf4RegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto& eval = polAqqQPS_reg;

  {
    testVal  = eval(testAs,testLM,testNF,9.31322574615478515625e-10);
    refVal   = -4702.95666659521128536240820155;
    errBound = 4702.95666659521128536240820155 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullPowersOf4, reg, x = 4^-15 = 9.31322574615478515625e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.7252902984619140625e-9);
    refVal   = -3884.26933339131506925944643446;
    errBound = 3884.26933339131506925944643446 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullPowersOf4, reg, x = 4^-14 = 3.7252902984619140625e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.490116119384765625e-8);
    refVal   = -3172.94134433191067199762700536;
    errBound = 3172.94134433191067199762700536 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullPowersOf4, reg, x = 4^-13 = 1.490116119384765625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,5.9604644775390625e-8);
    refVal   = -2559.53508950704257629731663141;
    errBound = 2559.53508950704257629731663141 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullPowersOf4, reg, x = 4^-12 = 5.9604644775390625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,2.384185791015625e-7);
    refVal   = -2035.02345035965887089636906741;
    errBound = 2035.02345035965887089636906741 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullPowersOf4, reg, x = 4^-11 = 2.384185791015625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,9.5367431640625e-7);
    refVal   = -1590.79004685727988749594445853;
    errBound = 1590.79004685727988749594445853 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullPowersOf4, reg, x = 4^-10 = 9.5367431640625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.814697265625e-6);
    refVal   = -1218.62998120724126895399650187;
    errBound = 1218.62998120724126895399650187 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullPowersOf4, reg, x = 4^-9 = 3.814697265625e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0000152587890625);
    refVal   = -910.752020171011256714333243046;
    errBound = 910.752020171011256714333243046 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullPowersOf4, reg, x = 4^-8 = 0.0000152587890625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00006103515625);
    refVal   = -659.78480642213907526930562572;
    errBound = 659.78480642213907526930562572 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullPowersOf4, reg, x = 4^-7 = 0.00006103515625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.000244140625);
    refVal   = -458.793877029308123685113234067;
    errBound = 458.793877029308123685113234067 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullPowersOf4, reg, x = 4^-6 = 0.000244140625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0009765625);
    refVal   = -301.32604528638193839093135071;
    errBound = 301.32604528638193839093135071 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullPowersOf4, reg, x = 4^-5 = 0.0009765625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00390625);
    refVal   = -181.517551201153881062427436742;
    errBound = 181.517551201153881062427436742 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullPowersOf4, reg, x = 4^-4 = 0.00390625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.015625);
    refVal   = -94.3319597330301033834260344888;
    errBound = 94.3319597330301033834260344888 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullPowersOf4, reg, x = 4^-3 = 0.015625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0625);
    refVal   = -35.9952200343124762861852214318;
    errBound = 35.9952200343124762861852214318 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullPowersOf4, reg, x = 4^-2 = 0.0625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9375);
    refVal   = 1.00541342787457013490216005375;
    errBound = 1.00541342787457013490216005375 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullPowersOf4, reg, x = 1-4^-2 = 0.9375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.984375);
    refVal   = 0.431251498252852453002391146623;
    errBound = 0.431251498252852453002391146623 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullPowersOf4, reg, x = 1-4^-3 = 0.984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99609375);
    refVal   = 0.170743580030682888433255825441;
    errBound = 0.170743580030682888433255825441 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullPowersOf4, reg, x = 1-4^-4 = 0.99609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9990234375);
    refVal   = 0.0636689389160134136082754102319;
    errBound = 0.0636689389160134136082754102319 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullPowersOf4, reg, x = 1-4^-5 = 0.9990234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999755859375);
    refVal   = 0.0226571326542703407590017168651;
    errBound = 0.0226571326542703407590017168651 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullPowersOf4, reg, x = 1-4^-6 = 0.999755859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99993896484375);
    refVal   = 0.00776798067275199002905563220042;
    errBound = 0.00776798067275199002905563220042 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullPowersOf4, reg, x = 1-4^-7 = 0.99993896484375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999847412109375);
    refVal   = 0.00258389306472717252945935894924;
    errBound = 0.00258389306472717252945935894924 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullPowersOf4, reg, x = 1-4^-8 = 0.9999847412109375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999996185302734375);
    refVal   = 0.000838264213990267670007788447864;
    errBound = 0.000838264213990267670007788447864 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullPowersOf4, reg, x = 1-4^-9 = 0.999996185302734375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999904632568359375);
    refVal   = 0.000266298680289548958641089972013;
    errBound = 0.000266298680289548958641089972013 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullPowersOf4, reg, x = 1-4^-10 = 0.99999904632568359375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999997615814208984375);
    refVal   = 0.0000830992420436494315687803599776;
    errBound = 0.0000830992420436494315687803599776 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullPowersOf4, reg, x = 1-4^-11 = 0.9999997615814208984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999940395355224609375);
    refVal   = 0.0000255354180245240949688184037349;
    errBound = 0.0000255354180245240949688184037349 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullPowersOf4, reg, x = 1-4^-12 = 0.999999940395355224609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999998509883880615234375);
    refVal   = 7.7424001903362100688097868537e-6;
    errBound = 7.7424001903362100688097868537e-6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullPowersOf4, reg, x = 1-4^-13 = 0.99999998509883880615234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999962747097015380859375);
    refVal   = 2.32009183898638998670835048517e-6;
    errBound = 2.32009183898638998670835048517e-6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullPowersOf4, reg, x = 1-4^-14 = 0.9999999962747097015380859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999068677425384521484375);
    refVal   = 6.88048979910055375605364091871e-7;
    errBound = 6.88048979910055375605364091871e-7 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullPowersOf4, reg, x = 1-4^-15 = 0.999999999068677425384521484375";
  }
}



TEST(polAqqQPSXspace,FullNormalDoublesRegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto& eval = polAqqQPS_reg;

  {
    testVal  = eval(testAs,testLM,testNF,1.e-10);
    refVal   = -6271.27983456148402236618070676;
    errBound = 6271.27983456148402236618070676 * std::max(eps/(1.-1.e-10),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullNormalDoubles, reg, x = 10^-10 = 1.e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-9);
    refVal   = -4658.17104865177990106340178128;
    errBound = 4658.17104865177990106340178128 * std::max(eps/(1.-1.e-9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullNormalDoubles, reg, x = 10^-9 = 1.e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-8);
    refVal   = -3367.15810434698142348009528103;
    errBound = 3367.15810434698142348009528103 * std::max(eps/(1.-1.e-8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullNormalDoubles, reg, x = 10^-8 = 1.e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-7);
    refVal   = -2353.83275606021552993210557187;
    errBound = 2353.83275606021552993210557187 * std::max(eps/(1.-1.e-7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullNormalDoubles, reg, x = 10^-7 = 1.e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-6);
    refVal   = -1576.91103248722437846924384995;
    errBound = 1576.91103248722437846924384995 * std::max(eps/(1.-1.e-6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullNormalDoubles, reg, x = 10^-6 = 1.e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-5);
    refVal   = -998.238239843366910694548027747;
    errBound = 998.238239843366910694548027747 * std::max(eps/(1.-1.e-5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullNormalDoubles, reg, x = 10^-5 = 1.e-5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0001);
    refVal   = -582.81861980160413399954412398;
    errBound = 582.81861980160413399954412398 * std::max(eps/(1.-0.0001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullNormalDoubles, reg, x = 10^-4 = 0.0001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.001);
    refVal   = -298.975466332338462611522771224;
    errBound = 298.975466332338462611522771224 * std::max(eps/(1.-0.001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullNormalDoubles, reg, x = 10^-3 = 0.001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.01);
    refVal   = -119.10072086554307808924848154;
    errBound = 119.10072086554307808924848154 * std::max(eps/(1.-0.01),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullNormalDoubles, reg, x = 10^-2 = 0.01";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.05);
    refVal   = -43.5438238758275180349162540923;
    errBound = 43.5438238758275180349162540923 * std::max(eps/(1.-0.05),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullNormalDoubles, reg, x = 0.05 = 0.05";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.1);
    refVal   = -22.3526618903710881531332699834;
    errBound = 22.3526618903710881531332699834 * std::max(eps/(1.-0.1),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullNormalDoubles, reg, x = 0.1 = 0.1";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.15);
    refVal   = -13.0243730243280827491558464665;
    errBound = 13.0243730243280827491558464665 * std::max(eps/(1.-0.15),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullNormalDoubles, reg, x = 0.15 = 0.15";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.2);
    refVal   = -7.76963741786407613440845112337;
    errBound = 7.76963741786407613440845112337 * std::max(eps/(1.-0.2),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullNormalDoubles, reg, x = 0.2 = 0.2";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.25);
    refVal   = -4.47006755344809194513646571131;
    errBound = 4.47006755344809194513646571131 * std::max(eps/(1.-0.25),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullNormalDoubles, reg, x = 0.25 = 0.25";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.3);
    refVal   = -2.27530153879249694231025562235;
    errBound = 2.27530153879249694231025562235 * std::max(eps/(1.-0.3),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullNormalDoubles, reg, x = 0.3 = 0.3";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.35);
    refVal   = -0.769472784276148984905030377773;
    errBound = 0.769472784276148984905030377773 * std::max(eps/(1.-0.35),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullNormalDoubles, reg, x = 0.35 = 0.35";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.4);
    refVal   = 0.277223932671857719614091676307;
    errBound = 0.277223932671857719614091676307 * std::max(eps/(1.-0.4),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullNormalDoubles, reg, x = 0.4 = 0.4";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.45);
    refVal   = 1.00279838897181238814020387219;
    errBound = 1.00279838897181238814020387219 * std::max(eps/(1.-0.45),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullNormalDoubles, reg, x = 0.45 = 0.45";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.5);
    refVal   = 1.49522532309716582590719811747;
    errBound = 1.49522532309716582590719811747 * std::max(eps/(1.-0.5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullNormalDoubles, reg, x = 0.5 = 0.5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.55);
    refVal   = 1.81306355169848553838464284782;
    errBound = 1.81306355169848553838464284782 * std::max(eps/(1.-0.55),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullNormalDoubles, reg, x = 0.55 = 0.55";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.6);
    refVal   = 1.99646858288444395175817873164;
    errBound = 1.99646858288444395175817873164 * std::max(eps/(1.-0.6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullNormalDoubles, reg, x = 0.6 = 0.6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.65);
    refVal   = 2.07340045614604552185322686773;
    errBound = 2.07340045614604552185322686773 * std::max(eps/(1.-0.65),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullNormalDoubles, reg, x = 0.65 = 0.65";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.7);
    refVal   = 2.0632007078016893221571502742;
    errBound = 2.0632007078016893221571502742 * std::max(eps/(1.-0.7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullNormalDoubles, reg, x = 0.7 = 0.7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.75);
    refVal   = 1.97854185796427831873526697426;
    errBound = 1.97854185796427831873526697426 * std::max(eps/(1.-0.75),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullNormalDoubles, reg, x = 0.75 = 0.75";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.8);
    refVal   = 1.82609717925190692031805219242;
    errBound = 1.82609717925190692031805219242 * std::max(eps/(1.-0.8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullNormalDoubles, reg, x = 0.8 = 0.8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.85);
    refVal   = 1.60561248273470638937616889201;
    errBound = 1.60561248273470638937616889201 * std::max(eps/(1.-0.85),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullNormalDoubles, reg, x = 0.85 = 0.85";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9);
    refVal   = 1.3053817757742156115518444699;
    errBound = 1.3053817757742156115518444699 * std::max(eps/(1.-0.9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullNormalDoubles, reg, x = 0.9 = 0.9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.95);
    refVal   = 0.883320850776577964034150143304;
    errBound = 0.883320850776577964034150143304 * std::max(eps/(1.-0.95),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullNormalDoubles, reg, x = 0.95 = 0.95";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99);
    refVal   = 0.322455509688682580710251884004;
    errBound = 0.322455509688682580710251884004 * std::max(eps/(1.-0.99),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullNormalDoubles, reg, x = 1-10^-2 = 0.99";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999);
    refVal   = 0.0647800603558420821469480042175;
    errBound = 0.0647800603558420821469480042175 * std::max(eps/(1.-0.999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullNormalDoubles, reg, x = 1-10^-3 = 0.999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999);
    refVal   = 0.0114162339022702941123185548191;
    errBound = 0.0114162339022702941123185548191 * std::max(eps/(1.-0.9999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullNormalDoubles, reg, x = 1-10^-4 = 0.9999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999);
    refVal   = 0.00183785569552584311290593502414;
    errBound = 0.00183785569552584311290593502414 * std::max(eps/(1.-0.99999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullNormalDoubles, reg, x = 1-10^-5 = 0.99999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999);
    refVal   = 0.000277040918720688364130162289562;
    errBound = 0.000277040918720688364130162289562 * std::max(eps/(1.-0.999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullNormalDoubles, reg, x = 1-10^-6 = 0.999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999);
    refVal   = 0.0000397319518090177418495575895881;
    errBound = 0.0000397319518090177418495575895881 * std::max(eps/(1.-0.9999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullNormalDoubles, reg, x = 1-10^-7 = 0.9999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999999);
    refVal   = 5.48012498656899546661620642735e-6;
    errBound = 5.48012498656899546661620642735e-6 * std::max(eps/(1.-0.99999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullNormalDoubles, reg, x = 1-10^-8 = 0.99999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999);
    refVal   = 7.32510991408135370431506142263e-7;
    errBound = 7.32510991408135370431506142263e-7 * std::max(eps/(1.-0.999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullNormalDoubles, reg, x = 1-10^-9 = 0.999999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999999);
    refVal   = 9.54206127331042407563686953315e-8;
    errBound = 9.54206127331042407563686953315e-8 * std::max(eps/(1.-0.9999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQPS, FullNormalDoubles, reg, x = 1-10^-10 = 0.9999999999";
  }
}



TEST(polAqqQPSNspace,FullMoments)
{
  double refVal = 0., errBound = 0.;
  auto& rpd = polAqqQPS;
  integration_engine_gsl engine;
  auto mom = make_mellin_moment(rpd, engine, testAs, testLM, testNF);
  std::tuple<integration_status, double, double> res;

  {
    SCOPED_TRACE("polAqqQPS Full Mellin moment N=3");
    refVal = 0.384961282733549434963183505214;
    errBound = 0.384961282733549434963183505214 * 3.e-11;
    res = mom.integrate(3, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAqqQPS Full Mellin moment N=5");
    refVal = 0.264771441815979026607447682672;
    errBound = 0.264771441815979026607447682672 * 3.e-11;
    res = mom.integrate(5, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAqqQPS Full Mellin moment N=7");
    refVal = 0.174637951875937731531411329203;
    errBound = 0.174637951875937731531411329203 * 3.e-11;
    res = mom.integrate(7, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAqqQPS Full Mellin moment N=9");
    refVal = 0.12360315529395906598053582999;
    errBound = 0.12360315529395906598053582999 * 3.e-11;
    res = mom.integrate(9, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAqqQPS Full Mellin moment N=11");
    refVal = 0.0925816953054485219376788338187;
    errBound = 0.0925816953054485219376788338187 * 3.e-11;
    res = mom.integrate(11, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

}
