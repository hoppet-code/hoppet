/******************************************************************************
 * This file is part of libome                                                *
 * Copyright (C) 2025 Arnd Behring, Kay Schoenwald                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 ******************************************************************************/

#include <limits>
#include <gtest/gtest.h>
#include <ome/integration_engine_gsl.h>
#include <ome/mellin.h>
#include <ome/polAqgQ.h>

using namespace ome;

const double eps = std::numeric_limits<double>::epsilon();
const double testAs = 0.25;
const double testLM = -5.;
const double testNF = 3.;



TEST(polAqgQXspace,FullPowersOf4RegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto& eval = polAqgQ_reg;

  {
    testVal  = eval(testAs,testLM,testNF,9.31322574615478515625e-10);
    refVal   = 85983.963609840719402926960169;
    errBound = 85983.963609840719402926960169 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullPowersOf4, reg, x = 4^-15 = 9.31322574615478515625e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.7252902984619140625e-9);
    refVal   = 61147.1099164621832970253460334;
    errBound = 61147.1099164621832970253460334 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullPowersOf4, reg, x = 4^-14 = 3.7252902984619140625e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.490116119384765625e-8);
    refVal   = 42297.6202150377890676339635687;
    errBound = 42297.6202150377890676339635687 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullPowersOf4, reg, x = 4^-13 = 1.490116119384765625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,5.9604644775390625e-8);
    refVal   = 28304.8854255093830984284679704;
    errBound = 28304.8854255093830984284679704 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullPowersOf4, reg, x = 4^-12 = 5.9604644775390625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,2.384185791015625e-7);
    refVal   = 18187.212492485574716822344172;
    errBound = 18187.212492485574716822344172 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullPowersOf4, reg, x = 4^-11 = 2.384185791015625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,9.5367431640625e-7);
    refVal   = 11101.574892491987877539208412;
    errBound = 11101.574892491987877539208412 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullPowersOf4, reg, x = 4^-10 = 9.5367431640625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.814697265625e-6);
    refVal   = 6333.34916175536057324923707827;
    errBound = 6333.34916175536057324923707827 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullPowersOf4, reg, x = 4^-9 = 3.814697265625e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0000152587890625);
    refVal   = 3286.0196665393992092489250591;
    errBound = 3286.0196665393992092489250591 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullPowersOf4, reg, x = 4^-8 = 0.0000152587890625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00006103515625);
    refVal   = 1470.81762363990534817600529799;
    errBound = 1470.81762363990534817600529799 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullPowersOf4, reg, x = 4^-7 = 0.00006103515625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.000244140625);
    refVal   = 496.241943944046969636123371475;
    errBound = 496.241943944046969636123371475 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullPowersOf4, reg, x = 4^-6 = 0.000244140625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0009765625);
    refVal   = 57.4164104248711519908727937968;
    errBound = 57.4164104248711519908727937968 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullPowersOf4, reg, x = 4^-5 = 0.0009765625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00390625);
    refVal   = -74.6491675733096256880445957687;
    errBound = 74.6491675733096256880445957687 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullPowersOf4, reg, x = 4^-4 = 0.00390625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.015625);
    refVal   = -63.4674267984405642506810625141;
    errBound = 63.4674267984405642506810625141 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullPowersOf4, reg, x = 4^-3 = 0.015625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0625);
    refVal   = -16.0941049750992541404468517856;
    errBound = 16.0941049750992541404468517856 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullPowersOf4, reg, x = 4^-2 = 0.0625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9375);
    refVal   = -14.9343205339335993467288161592;
    errBound = 14.9343205339335993467288161592 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullPowersOf4, reg, x = 1-4^-2 = 0.9375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.984375);
    refVal   = -39.3726956322273503180571162862;
    errBound = 39.3726956322273503180571162862 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullPowersOf4, reg, x = 1-4^-3 = 0.984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99609375);
    refVal   = -74.4411019745222138435004729823;
    errBound = 74.4411019745222138435004729823 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullPowersOf4, reg, x = 1-4^-4 = 0.99609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9990234375);
    refVal   = -123.901831827239980801953606125;
    errBound = 123.901831827239980801953606125 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullPowersOf4, reg, x = 1-4^-5 = 0.9990234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999755859375);
    refVal   = -191.469930785662438298835754022;
    errBound = 191.469930785662438298835754022 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullPowersOf4, reg, x = 1-4^-6 = 0.999755859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99993896484375);
    refVal   = -280.90778035446262131135728375;
    errBound = 280.90778035446262131135728375 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullPowersOf4, reg, x = 1-4^-7 = 0.99993896484375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999847412109375);
    refVal   = -396.123838024035154757623306336;
    errBound = 396.123838024035154757623306336 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullPowersOf4, reg, x = 1-4^-8 = 0.9999847412109375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999996185302734375);
    refVal   = -541.23303784950248473599791149;
    errBound = 541.23303784950248473599791149 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullPowersOf4, reg, x = 1-4^-9 = 0.999996185302734375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999904632568359375);
    refVal   = -720.586258834415238105135269058;
    errBound = 720.586258834415238105135269058 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullPowersOf4, reg, x = 1-4^-10 = 0.99999904632568359375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999997615814208984375);
    refVal   = -938.782982799857731744598553002;
    errBound = 938.782982799857731744598553002 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullPowersOf4, reg, x = 1-4^-11 = 0.9999997615814208984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999940395355224609375);
    refVal   = -1200.67629953992039974226314165;
    errBound = 1200.67629953992039974226314165 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullPowersOf4, reg, x = 1-4^-12 = 0.999999940395355224609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999998509883880615234375);
    refVal   = -1511.37477310488000116352987288;
    errBound = 1511.37477310488000116352987288 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullPowersOf4, reg, x = 1-4^-13 = 0.99999998509883880615234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999962747097015380859375);
    refVal   = -1876.24310763554643116880971005;
    errBound = 1876.24310763554643116880971005 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullPowersOf4, reg, x = 1-4^-14 = 0.9999999962747097015380859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999068677425384521484375);
    refVal   = -2300.90237683626726093427908397;
    errBound = 2300.90237683626726093427908397 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullPowersOf4, reg, x = 1-4^-15 = 0.999999999068677425384521484375";
  }
}



TEST(polAqgQXspace,FullNormalDoublesRegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto& eval = polAqgQ_reg;

  {
    testVal  = eval(testAs,testLM,testNF,1.e-10);
    refVal   = 141861.889535542969284585455058;
    errBound = 141861.889535542969284585455058 * std::max(eps/(1.-1.e-10),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullNormalDoubles, reg, x = 10^-10 = 1.e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-9);
    refVal   = 84543.7778265445756756882292078;
    errBound = 84543.7778265445756756882292078 * std::max(eps/(1.-1.e-9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullNormalDoubles, reg, x = 10^-9 = 1.e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-8);
    refVal   = 47176.3595343951425401456196008;
    errBound = 47176.3595343951425401456196008 * std::max(eps/(1.-1.e-8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullNormalDoubles, reg, x = 10^-8 = 1.e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-7);
    refVal   = 24123.2417965836411346399335854;
    errBound = 24123.2417965836411346399335854 * std::max(eps/(1.-1.e-7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullNormalDoubles, reg, x = 10^-7 = 1.e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-6);
    refVal   = 10903.8537808397395991788421767;
    errBound = 10903.8537808397395991788421767 * std::max(eps/(1.-1.e-6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullNormalDoubles, reg, x = 10^-6 = 1.e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-5);
    refVal   = 4063.83111335000749637334155927;
    errBound = 4063.83111335000749637334155927 * std::max(eps/(1.-1.e-5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullNormalDoubles, reg, x = 10^-5 = 1.e-5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0001);
    refVal   = 1044.88704324968536367511206457;
    errBound = 1044.88704324968536367511206457 * std::max(eps/(1.-0.0001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullNormalDoubles, reg, x = 10^-4 = 0.0001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.001);
    refVal   = 53.1308429721627730970124929605;
    errBound = 53.1308429721627730970124929605 * std::max(eps/(1.-0.001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullNormalDoubles, reg, x = 10^-3 = 0.001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.01);
    refVal   = -75.1085963831326329677490133525;
    errBound = 75.1085963831326329677490133525 * std::max(eps/(1.-0.01),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullNormalDoubles, reg, x = 10^-2 = 0.01";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.05);
    refVal   = -23.2345297264486274202212266883;
    errBound = 23.2345297264486274202212266883 * std::max(eps/(1.-0.05),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullNormalDoubles, reg, x = 0.05 = 0.05";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.1);
    refVal   = -3.27720038494824537802627158297;
    errBound = 3.27720038494824537802627158297 * std::max(eps/(1.-0.1),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullNormalDoubles, reg, x = 0.1 = 0.1";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.15);
    refVal   = 4.77181987126501350212886763301;
    errBound = 4.77181987126501350212886763301 * std::max(eps/(1.-0.15),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullNormalDoubles, reg, x = 0.15 = 0.15";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.2);
    refVal   = 8.59629169180838156947876282325;
    errBound = 8.59629169180838156947876282325 * std::max(eps/(1.-0.2),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullNormalDoubles, reg, x = 0.2 = 0.2";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.25);
    refVal   = 10.452727907627936702336231476;
    errBound = 10.452727907627936702336231476 * std::max(eps/(1.-0.25),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullNormalDoubles, reg, x = 0.25 = 0.25";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.3);
    refVal   = 11.2460755073695596390319908843;
    errBound = 11.2460755073695596390319908843 * std::max(eps/(1.-0.3),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullNormalDoubles, reg, x = 0.3 = 0.3";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.35);
    refVal   = 11.4041431908610452459964726383;
    errBound = 11.4041431908610452459964726383 * std::max(eps/(1.-0.35),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullNormalDoubles, reg, x = 0.35 = 0.35";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.4);
    refVal   = 11.1492674181827818210958838911;
    errBound = 11.1492674181827818210958838911 * std::max(eps/(1.-0.4),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullNormalDoubles, reg, x = 0.4 = 0.4";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.45);
    refVal   = 10.6023111536867835247807694944;
    errBound = 10.6023111536867835247807694944 * std::max(eps/(1.-0.45),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullNormalDoubles, reg, x = 0.45 = 0.45";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.5);
    refVal   = 9.82795232960019276822898276663;
    errBound = 9.82795232960019276822898276663 * std::max(eps/(1.-0.5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullNormalDoubles, reg, x = 0.5 = 0.5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.55);
    refVal   = 8.85570763649434832432682762308;
    errBound = 8.85570763649434832432682762308 * std::max(eps/(1.-0.55),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullNormalDoubles, reg, x = 0.55 = 0.55";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.6);
    refVal   = 7.68919216464029534819931518507;
    errBound = 7.68919216464029534819931518507 * std::max(eps/(1.-0.6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullNormalDoubles, reg, x = 0.6 = 0.6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.65);
    refVal   = 6.30815356942624088317864779482;
    errBound = 6.30815356942624088317864779482 * std::max(eps/(1.-0.65),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullNormalDoubles, reg, x = 0.65 = 0.65";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.7);
    refVal   = 4.66384349559683680170599832864;
    errBound = 4.66384349559683680170599832864 * std::max(eps/(1.-0.7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullNormalDoubles, reg, x = 0.7 = 0.7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.75);
    refVal   = 2.66437068287244920608339070629;
    errBound = 2.66437068287244920608339070629 * std::max(eps/(1.-0.75),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullNormalDoubles, reg, x = 0.75 = 0.75";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.8);
    refVal   = 0.138325547156137073549973930555;
    errBound = 0.138325547156137073549973930555 * std::max(eps/(1.-0.8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullNormalDoubles, reg, x = 0.8 = 0.8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.85);
    refVal   = -3.26269655150377934328733578201;
    errBound = 3.26269655150377934328733578201 * std::max(eps/(1.-0.85),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullNormalDoubles, reg, x = 0.85 = 0.85";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9);
    refVal   = -8.38984649207653939384371608647;
    errBound = 8.38984649207653939384371608647 * std::max(eps/(1.-0.9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullNormalDoubles, reg, x = 0.9 = 0.9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.95);
    refVal   = -18.3074938437820318233899283319;
    errBound = 18.3074938437820318233899283319 * std::max(eps/(1.-0.95),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullNormalDoubles, reg, x = 0.95 = 0.95";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99);
    refVal   = -49.3195963025339874908001569915;
    errBound = 49.3195963025339874908001569915 * std::max(eps/(1.-0.99),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullNormalDoubles, reg, x = 1-10^-2 = 0.99";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999);
    refVal   = -122.914016222627472062201070247;
    errBound = 122.914016222627472062201070247 * std::max(eps/(1.-0.999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullNormalDoubles, reg, x = 1-10^-3 = 0.999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999);
    refVal   = -246.306034008438484577555205762;
    errBound = 246.306034008438484577555205762 * std::max(eps/(1.-0.9999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullNormalDoubles, reg, x = 1-10^-4 = 0.9999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999);
    refVal   = -436.993073593771785142481485759;
    errBound = 436.993073593771785142481485759 * std::max(eps/(1.-0.99999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullNormalDoubles, reg, x = 1-10^-5 = 0.99999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999);
    refVal   = -713.834673354334964547035625858;
    errBound = 713.834673354334964547035625858 * std::max(eps/(1.-0.999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullNormalDoubles, reg, x = 1-10^-6 = 0.999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999);
    refVal   = -1097.49466421115010330912942297;
    errBound = 1097.49466421115010330912942297 * std::max(eps/(1.-0.9999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullNormalDoubles, reg, x = 1-10^-7 = 0.9999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999999);
    refVal   = -1610.55888565636578835702015292;
    errBound = 1610.55888565636578835702015292 * std::max(eps/(1.-0.99999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullNormalDoubles, reg, x = 1-10^-8 = 0.99999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999);
    refVal   = -2277.55987429198864632708476339;
    errBound = 2277.55987429198864632708476339 * std::max(eps/(1.-0.999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullNormalDoubles, reg, x = 1-10^-9 = 0.999999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999999);
    refVal   = -3124.98137087834241167707204727;
    errBound = 3124.98137087834241167707204727 * std::max(eps/(1.-0.9999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqgQ, FullNormalDoubles, reg, x = 1-10^-10 = 0.9999999999";
  }
}



TEST(polAqgQNspace,FullMoments)
{
  double refVal = 0., errBound = 0.;
  auto& rpd = polAqgQ;
  integration_engine_gsl engine;
  auto mom = make_mellin_moment(rpd, engine, testAs, testLM, testNF);
  std::tuple<integration_status, double, double> res;

  {
    SCOPED_TRACE("polAqgQ Full Mellin moment N=3");
    refVal = -1.60784142722128662055464770777;
    errBound = 1.60784142722128662055464770777 * 3.e-11;
    res = mom.integrate(3, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAqgQ Full Mellin moment N=5");
    refVal = -2.16641818102603542534083795309;
    errBound = 2.16641818102603542534083795309 * 3.e-11;
    res = mom.integrate(5, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAqgQ Full Mellin moment N=7");
    refVal = -2.20252140600818713960951388479;
    errBound = 2.20252140600818713960951388479 * 3.e-11;
    res = mom.integrate(7, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAqgQ Full Mellin moment N=9");
    refVal = -2.12659159958437711618050718182;
    errBound = 2.12659159958437711618050718182 * 3.e-11;
    res = mom.integrate(9, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAqgQ Full Mellin moment N=11");
    refVal = -2.02701727328366679876707040098;
    errBound = 2.02701727328366679876707040098 * 3.e-11;
    res = mom.integrate(11, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

}
