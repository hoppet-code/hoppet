/******************************************************************************
 * This file is part of libome                                                *
 * Copyright (C) 2025 Arnd Behring, Kay Schoenwald                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 ******************************************************************************/

#include <limits>
#include <gtest/gtest.h>
#include <ome/integration_engine_gsl.h>
#include <ome/mellin.h>
#include <ome/AqqQPS.h>

using namespace ome;

const double eps = std::numeric_limits<double>::epsilon();
const double testAs = 0.25;
const double testLM = -5.;
const double testNF = 3.;



TEST(AqqQPSXspace,FullPowersOf4RegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto& eval = AqqQPS_reg;

  {
    testVal  = eval(testAs,testLM,testNF,9.31322574615478515625e-10);
    refVal   = 1.41063550790701309630467295285e10;
    errBound = 1.41063550790701309630467295285e10 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullPowersOf4, reg, x = 4^-15 = 9.31322574615478515625e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.7252902984619140625e-9);
    refVal   = 3.52658685575897270466605919847e9;
    errBound = 3.52658685575897270466605919847e9 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullPowersOf4, reg, x = 4^-14 = 3.7252902984619140625e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.490116119384765625e-8);
    refVal   = 8.81645173891556427899541072818e8;
    errBound = 8.81645173891556427899541072818e8 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullPowersOf4, reg, x = 4^-13 = 1.490116119384765625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,5.9604644775390625e-8);
    refVal   = 2.20410070265037907946089443412e8;
    errBound = 2.20410070265037907946089443412e8 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullPowersOf4, reg, x = 4^-12 = 5.9604644775390625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,2.384185791015625e-7);
    refVal   = 5.51015600804338995273622140261e7;
    errBound = 5.51015600804338995273622140261e7 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullPowersOf4, reg, x = 4^-11 = 2.384185791015625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,9.5367431640625e-7);
    refVal   = 1.37746528317486737709515364662e7;
    errBound = 1.37746528317486737709515364662e7 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullPowersOf4, reg, x = 4^-10 = 9.5367431640625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.814697265625e-6);
    refVal   = 3.44310627787524353509206898791e6;
    errBound = 3.44310627787524353509206898791e6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullPowersOf4, reg, x = 4^-9 = 3.814697265625e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0000152587890625);
    refVal   = 860364.934507150221927768845758;
    errBound = 860364.934507150221927768845758 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullPowersOf4, reg, x = 4^-8 = 0.0000152587890625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00006103515625);
    refVal   = 214794.694264937765750565719293;
    errBound = 214794.694264937765750565719293 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullPowersOf4, reg, x = 4^-7 = 0.00006103515625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.000244140625);
    refVal   = 53491.4745047636494330399810293;
    errBound = 53491.4745047636494330399810293 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullPowersOf4, reg, x = 4^-6 = 0.000244140625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0009765625);
    refVal   = 13233.3624770775977448262283666;
    errBound = 13233.3624770775977448262283666 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullPowersOf4, reg, x = 4^-5 = 0.0009765625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00390625);
    refVal   = 3218.61823603683479947196847506;
    errBound = 3218.61823603683479947196847506 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullPowersOf4, reg, x = 4^-4 = 0.00390625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.015625);
    refVal   = 750.132574681304448315845820349;
    errBound = 750.132574681304448315845820349 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullPowersOf4, reg, x = 4^-3 = 0.015625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0625);
    refVal   = 156.599903312514476259703470914;
    errBound = 156.599903312514476259703470914 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullPowersOf4, reg, x = 4^-2 = 0.0625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9375);
    refVal   = 1.06720576278165102916813268185;
    errBound = 1.06720576278165102916813268185 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullPowersOf4, reg, x = 1-4^-2 = 0.9375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.984375);
    refVal   = 0.44335396715543417083034907851;
    errBound = 0.44335396715543417083034907851 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullPowersOf4, reg, x = 1-4^-3 = 0.984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99609375);
    refVal   = 0.173056950686560363316270562907;
    errBound = 0.173056950686560363316270562907 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullPowersOf4, reg, x = 1-4^-4 = 0.99609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9990234375);
    refVal   = 0.0639774079201381208989901232223;
    errBound = 0.0639774079201381208989901232223 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullPowersOf4, reg, x = 1-4^-5 = 0.9990234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999755859375);
    refVal   = 0.022640999057537964585414781222;
    errBound = 0.022640999057537964585414781222 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullPowersOf4, reg, x = 1-4^-6 = 0.999755859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99993896484375);
    refVal   = 0.00773412425816817205814877227051;
    errBound = 0.00773412425816817205814877227051 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullPowersOf4, reg, x = 1-4^-7 = 0.99993896484375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999847412109375);
    refVal   = 0.00256634420048336919010939408707;
    errBound = 0.00256634420048336919010939408707 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullPowersOf4, reg, x = 1-4^-8 = 0.9999847412109375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999996185302734375);
    refVal   = 0.000831198525482569942468373965986;
    errBound = 0.000831198525482569942468373965986 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullPowersOf4, reg, x = 1-4^-9 = 0.999996185302734375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999904632568359375);
    refVal   = 0.000263760818840419925725584975931;
    errBound = 0.000263760818840419925725584975931 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullPowersOf4, reg, x = 1-4^-10 = 0.99999904632568359375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999997615814208984375);
    refVal   = 0.0000822464615352829550745303861325;
    errBound = 0.0000822464615352829550745303861325 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullPowersOf4, reg, x = 1-4^-11 = 0.9999997615814208984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999940395355224609375);
    refVal   = 0.0000252612802817717987106133371735;
    errBound = 0.0000252612802817717987106133371735 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullPowersOf4, reg, x = 1-4^-12 = 0.999999940395355224609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999998509883880615234375);
    refVal   = 7.65703914344305742977144190395e-6;
    errBound = 7.65703914344305742977144190395e-6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullPowersOf4, reg, x = 1-4^-13 = 0.99999998509883880615234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999962747097015380859375);
    refVal   = 2.29414718513925867981502427334e-6;
    errBound = 2.29414718513925867981502427334e-6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullPowersOf4, reg, x = 1-4^-14 = 0.9999999962747097015380859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999068677425384521484375);
    refVal   = 6.80312283586666561919627471812e-7;
    errBound = 6.80312283586666561919627471812e-7 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullPowersOf4, reg, x = 1-4^-15 = 0.999999999068677425384521484375";
  }
}



TEST(AqqQPSXspace,FullNormalDoublesRegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto& eval = AqqQPS_reg;

  {
    testVal  = eval(testAs,testLM,testNF,1.e-10);
    refVal   = 1.31375696145868182054634539093e11;
    errBound = 1.31375696145868182054634539093e11 * std::max(eps/(1.-1.e-10),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullNormalDoubles, reg, x = 10^-10 = 1.e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-9);
    refVal   = 1.31375667329197743358227988698e10;
    errBound = 1.31375667329197743358227988698e10 * std::max(eps/(1.-1.e-9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullNormalDoubles, reg, x = 10^-9 = 1.e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-8);
    refVal   = 1.31375463443418844646861972373e9;
    errBound = 1.31375463443418844646861972373e9 * std::max(eps/(1.-1.e-8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullNormalDoubles, reg, x = 10^-8 = 1.e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-7);
    refVal   = 1.3137406910092429259494694422e8;
    errBound = 1.3137406910092429259494694422e8 * std::max(eps/(1.-1.e-7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullNormalDoubles, reg, x = 10^-7 = 1.e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-6);
    refVal   = 1.31364919700784297806808303087e7;
    errBound = 1.31364919700784297806808303087e7 * std::max(eps/(1.-1.e-6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullNormalDoubles, reg, x = 10^-6 = 1.e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-5);
    refVal   = 1.31307891347092888364989622452e6;
    errBound = 1.31307891347092888364989622452e6 * std::max(eps/(1.-1.e-5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullNormalDoubles, reg, x = 10^-5 = 1.e-5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0001);
    refVal   = 130975.055561612312914297250594;
    errBound = 130975.055561612312914297250594 * std::max(eps/(1.-0.0001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullNormalDoubles, reg, x = 10^-4 = 0.0001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.001);
    refVal   = 12919.5435916486702552379763827;
    errBound = 12919.5435916486702552379763827 * std::max(eps/(1.-0.001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullNormalDoubles, reg, x = 10^-3 = 0.001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.01);
    refVal   = 1207.75518283762074083140984558;
    errBound = 1207.75518283762074083140984558 * std::max(eps/(1.-0.01),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullNormalDoubles, reg, x = 10^-2 = 0.01";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.05);
    refVal   = 204.155776357254703554082850272;
    errBound = 204.155776357254703554082850272 * std::max(eps/(1.-0.05),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullNormalDoubles, reg, x = 0.05 = 0.05";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.1);
    refVal   = 87.2408688156872525476936053661;
    errBound = 87.2408688156872525476936053661 * std::max(eps/(1.-0.1),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullNormalDoubles, reg, x = 0.1 = 0.1";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.15);
    refVal   = 50.5978261837396571284553056977;
    errBound = 50.5978261837396571284553056977 * std::max(eps/(1.-0.15),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullNormalDoubles, reg, x = 0.15 = 0.15";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.2);
    refVal   = 33.2922621541876996148487802053;
    errBound = 33.2922621541876996148487802053 * std::max(eps/(1.-0.2),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullNormalDoubles, reg, x = 0.2 = 0.2";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.25);
    refVal   = 23.486234784281615627102825854;
    errBound = 23.486234784281615627102825854 * std::max(eps/(1.-0.25),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullNormalDoubles, reg, x = 0.25 = 0.25";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.3);
    refVal   = 17.3254390601734715214525422754;
    errBound = 17.3254390601734715214525422754 * std::max(eps/(1.-0.3),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullNormalDoubles, reg, x = 0.3 = 0.3";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.35);
    refVal   = 13.1912855467511210916124289714;
    errBound = 13.1912855467511210916124289714 * std::max(eps/(1.-0.35),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullNormalDoubles, reg, x = 0.35 = 0.35";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.4);
    refVal   = 10.2887693571585834025154623042;
    errBound = 10.2887693571585834025154623042 * std::max(eps/(1.-0.4),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullNormalDoubles, reg, x = 0.4 = 0.4";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.45);
    refVal   = 8.18294746972746694469070551126;
    errBound = 8.18294746972746694469070551126 * std::max(eps/(1.-0.45),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullNormalDoubles, reg, x = 0.45 = 0.45";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.5);
    refVal   = 6.61609320946221467665535163018;
    errBound = 6.61609320946221467665535163018 * std::max(eps/(1.-0.5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullNormalDoubles, reg, x = 0.5 = 0.5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.55);
    refVal   = 5.42561621542884564356659450762;
    errBound = 5.42561621542884564356659450762 * std::max(eps/(1.-0.55),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullNormalDoubles, reg, x = 0.55 = 0.55";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.6);
    refVal   = 4.50344759021697299799628269277;
    errBound = 4.50344759021697299799628269277 * std::max(eps/(1.-0.6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullNormalDoubles, reg, x = 0.6 = 0.6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.65);
    refVal   = 3.77430432276167230248063693508;
    errBound = 3.77430432276167230248063693508 * std::max(eps/(1.-0.65),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullNormalDoubles, reg, x = 0.65 = 0.65";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.7);
    refVal   = 3.1832114404749293707204848522;
    errBound = 3.1832114404749293707204848522 * std::max(eps/(1.-0.7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullNormalDoubles, reg, x = 0.7 = 0.7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.75);
    refVal   = 2.68775290258963401426425668174;
    errBound = 2.68775290258963401426425668174 * std::max(eps/(1.-0.75),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullNormalDoubles, reg, x = 0.75 = 0.75";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.8);
    refVal   = 2.25261449643185423236338204807;
    errBound = 2.25261449643185423236338204807 * std::max(eps/(1.-0.8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullNormalDoubles, reg, x = 0.8 = 0.8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.85);
    refVal   = 1.84460651506539911927551680471;
    errBound = 1.84460651506539911927551680471 * std::max(eps/(1.-0.85),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullNormalDoubles, reg, x = 0.85 = 0.85";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9);
    refVal   = 1.42530405166203142888658208185;
    errBound = 1.42530405166203142888658208185 * std::max(eps/(1.-0.9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullNormalDoubles, reg, x = 0.9 = 0.9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.95);
    refVal   = 0.929929661134234321228170583717;
    errBound = 0.929929661134234321228170583717 * std::max(eps/(1.-0.95),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullNormalDoubles, reg, x = 0.95 = 0.95";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99);
    refVal   = 0.329692970935589827477807570862;
    errBound = 0.329692970935589827477807570862 * std::max(eps/(1.-0.99),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullNormalDoubles, reg, x = 1-10^-2 = 0.99";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999);
    refVal   = 0.0651015414775625029500248734525;
    errBound = 0.0651015414775625029500248734525 * std::max(eps/(1.-0.999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullNormalDoubles, reg, x = 1-10^-3 = 0.999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999);
    refVal   = 0.0113793893616488395828643990834;
    errBound = 0.0113793893616488395828643990834 * std::max(eps/(1.-0.9999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullNormalDoubles, reg, x = 1-10^-4 = 0.9999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999);
    refVal   = 0.00182432771274337255292258140542;
    errBound = 0.00182432771274337255292258140542 * std::max(eps/(1.-0.99999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullNormalDoubles, reg, x = 1-10^-5 = 0.99999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999);
    refVal   = 0.000274409219705682745049354868681;
    errBound = 0.000274409219705682745049354868681 * std::max(eps/(1.-0.999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullNormalDoubles, reg, x = 1-10^-6 = 0.999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999);
    refVal   = 0.0000393114368474520411163329112263;
    errBound = 0.0000393114368474520411163329112263 * std::max(eps/(1.-0.9999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullNormalDoubles, reg, x = 1-10^-7 = 0.9999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999999);
    refVal   = 5.41939348546017400573294179703e-6;
    errBound = 5.41939348546017400573294179703e-6 * std::max(eps/(1.-0.99999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullNormalDoubles, reg, x = 1-10^-8 = 0.99999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999);
    refVal   = 7.24275290907007827093176928143e-7;
    errBound = 7.24275290907007827093176928143e-7 * std::max(eps/(1.-0.999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullNormalDoubles, reg, x = 1-10^-9 = 0.999999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999999);
    refVal   = 9.43513326544642974611852053737e-8;
    errBound = 9.43513326544642974611852053737e-8 * std::max(eps/(1.-0.9999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQPS, FullNormalDoubles, reg, x = 1-10^-10 = 0.9999999999";
  }
}



TEST(AqqQPSNspace,FullMoments)
{
  double refVal = 0., errBound = 0.;
  auto& rpd = AqqQPS;
  integration_engine_gsl engine;
  auto mom = make_mellin_moment(rpd, engine, testAs, testLM, testNF);
  std::tuple<integration_status, double, double> res;

  {
    SCOPED_TRACE("AqqQPS Full Mellin moment N=2");
    refVal = 4.2062933617580908693610822523;
    errBound = 4.2062933617580908693610822523 * 3.e-11;
    res = mom.integrate(2, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AqqQPS Full Mellin moment N=4");
    refVal = 0.679491035060561449752180563088;
    errBound = 0.679491035060561449752180563088 * 3.e-11;
    res = mom.integrate(4, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AqqQPS Full Mellin moment N=6");
    refVal = 0.305561977081704354671689330447;
    errBound = 0.305561977081704354671689330447 * 3.e-11;
    res = mom.integrate(6, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AqqQPS Full Mellin moment N=8");
    refVal = 0.18338742329786524282373118559;
    errBound = 0.18338742329786524282373118559 * 3.e-11;
    res = mom.integrate(8, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

}
