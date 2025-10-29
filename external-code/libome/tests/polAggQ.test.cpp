/******************************************************************************
 * This file is part of libome                                                *
 * Copyright (C) 2025 Arnd Behring, Kay Schoenwald                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 ******************************************************************************/

#include <limits>
#include <gtest/gtest.h>
#include <ome/integration_engine_gsl.h>
#include <ome/mellin.h>
#include <ome/polAggQ.h>

using namespace ome;

const double eps = std::numeric_limits<double>::epsilon();
const double testAs = 0.25;
const double testLM = -5.;
const double testNF = 3.;



TEST(polAggQXspace,FullPowersOf4RegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto& eval = polAggQ_reg;

  {
    testVal  = eval(testAs,testLM,testNF,9.31322574615478515625e-10);
    refVal   = -301352.223282196061961918550699;
    errBound = 301352.223282196061961918550699 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, reg, x = 4^-15 = 9.31322574615478515625e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.7252902984619140625e-9);
    refVal   = -215464.331827915517619740901978;
    errBound = 215464.331827915517619740901978 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, reg, x = 4^-14 = 3.7252902984619140625e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.490116119384765625e-8);
    refVal   = -149963.378454188811347633798094;
    errBound = 149963.378454188811347633798094 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, reg, x = 4^-13 = 1.490116119384765625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,5.9604644775390625e-8);
    refVal   = -101095.8032837299095895967577742;
    errBound = 101095.8032837299095895967577742 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, reg, x = 4^-12 = 5.9604644775390625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,2.384185791015625e-7);
    refVal   = -65578.5012850802051136822343977;
    errBound = 65578.5012850802051136822343977 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, reg, x = 4^-11 = 2.384185791015625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,9.5367431640625e-7);
    refVal   = -40569.3724116672592501759464275;
    errBound = 40569.3724116672592501759464275 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, reg, x = 4^-10 = 9.5367431640625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.814697265625e-6);
    refVal   = -23637.8596521242058503954491547;
    errBound = 23637.8596521242058503954491547 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, reg, x = 4^-9 = 3.814697265625e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0000152587890625);
    refVal   = -12735.46504446737433441489509918;
    errBound = 12735.46504446737433441489509918 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, reg, x = 4^-8 = 0.0000152587890625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00006103515625);
    refVal   = -6166.23980413308851917532068289;
    errBound = 6166.23980413308851917532068289 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, reg, x = 4^-7 = 0.00006103515625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.000244140625);
    refVal   = -2557.29921948932938057195731228;
    errBound = 2557.29921948932938057195731228 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, reg, x = 4^-6 = 0.000244140625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0009765625);
    refVal   = -829.632492214494714859431099269;
    errBound = 829.632492214494714859431099269 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, reg, x = 4^-5 = 0.0009765625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00390625);
    refVal   = -170.123071911824095106760138108;
    errBound = 170.123071911824095106760138108 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, reg, x = 4^-4 = 0.00390625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.015625);
    refVal   = -7.06449048108370461002458154404;
    errBound = 7.06449048108370461002458154404 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, reg, x = 4^-3 = 0.015625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0625);
    refVal   = 7.29724839852870475005158833589;
    errBound = 7.29724839852870475005158833589 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, reg, x = 4^-2 = 0.0625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9375);
    refVal   = -75.8170571727567236245409206211;
    errBound = 75.8170571727567236245409206211 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, reg, x = 1-4^-2 = 0.9375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.984375);
    refVal   = -104.2383520051001876369851908983;
    errBound = 104.2383520051001876369851908983 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, reg, x = 1-4^-3 = 0.984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99609375);
    refVal   = -131.0864142948087427522285746732;
    errBound = 131.0864142948087427522285746732 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, reg, x = 1-4^-4 = 0.99609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9990234375);
    refVal   = -160.908481843380413752088960859;
    errBound = 160.908481843380413752088960859 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, reg, x = 1-4^-5 = 0.9990234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999755859375);
    refVal   = -195.935519500196451215094612751;
    errBound = 195.935519500196451215094612751 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, reg, x = 1-4^-6 = 0.999755859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99993896484375);
    refVal   = -237.721479060201146643824328667;
    errBound = 237.721479060201146643824328667 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, reg, x = 1-4^-7 = 0.99993896484375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999847412109375);
    refVal   = -287.647031251753212672015256001;
    errBound = 287.647031251753212672015256001 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, reg, x = 1-4^-8 = 0.9999847412109375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999996185302734375);
    refVal   = -347.0531603653328256917377772;
    errBound = 347.0531603653328256917377772 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, reg, x = 1-4^-9 = 0.999996185302734375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999904632568359375);
    refVal   = -417.272987430386343719882714076;
    errBound = 417.272987430386343719882714076 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, reg, x = 1-4^-10 = 0.99999904632568359375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999997615814208984375);
    refVal   = -499.638472855707618994082181058;
    errBound = 499.638472855707618994082181058 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, reg, x = 1-4^-11 = 0.9999997615814208984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999940395355224609375);
    refVal   = -595.48156018272751290940278768;
    errBound = 595.48156018272751290940278768 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, reg, x = 1-4^-12 = 0.999999940395355224609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999998509883880615234375);
    refVal   = -706.134274229220303937324560118;
    errBound = 706.134274229220303937324560118 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, reg, x = 1-4^-13 = 0.99999998509883880615234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999962747097015380859375);
    refVal   = -832.9286853991835475451279636;
    errBound = 832.9286853991835475451279636 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, reg, x = 1-4^-14 = 0.9999999962747097015380859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999068677425384521484375);
    refVal   = -977.196882633215357569814403361;
    errBound = 977.196882633215357569814403361 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, reg, x = 1-4^-15 = 0.999999999068677425384521484375";
  }
}

TEST(polAggQXspace,FullPowersOf4PlusEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;
  auto& eval = polAggQ_plus;
  {
    testVal  = eval(testAs,testLM,testNF,9.31322574615478515625e-10);
    refVal   = 42.2734346371860400952092060316;
    errBound = 42.2734346371860400952092060316 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, plus, x = 4^-15 = 9.31322574615478515625e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.7252902984619140625e-9);
    refVal   = 42.2734347552966524876352738601;
    errBound = 42.2734347552966524876352738601 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, plus, x = 4^-14 = 3.7252902984619140625e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.490116119384765625e-8);
    refVal   = 42.2734352277391086572844208845;
    errBound = 42.2734352277391086572844208845 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, plus, x = 4^-13 = 1.490116119384765625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,5.9604644775390625e-8);
    refVal   = 42.273437117509038935005216199;
    errBound = 42.273437117509038935005216199 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, plus, x = 4^-12 = 5.9604644775390625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,2.384185791015625e-7);
    refVal   = 42.273444676590449632272247516;
    errBound = 42.273444676590449632272247516 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, plus, x = 4^-11 = 2.384185791015625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,9.5367431640625e-7);
    refVal   = 42.2734749129431258288602113757;
    errBound = 42.2734749129431258288602113757 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, plus, x = 4^-10 = 9.5367431640625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.814697265625e-6);
    refVal   = 42.2735958587863667597428656719;
    errBound = 42.2735958587863667597428656719 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, plus, x = 4^-9 = 3.814697265625e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0000152587890625);
    refVal   = 42.2740796490800127470026475769;
    errBound = 42.2740796490800127470026475769 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, plus, x = 4^-8 = 0.0000152587890625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00006103515625);
    refVal   = 42.2760149209921661993349586302;
    errBound = 42.2760149209921661993349586302 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, plus, x = 4^-7 = 0.00006103515625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.000244140625);
    refVal   = 42.2837577808678057902139577069;
    errBound = 42.2837577808678057902139577069 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, plus, x = 4^-6 = 0.000244140625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0009765625);
    refVal   = 42.3147576032877968501774576759;
    errBound = 42.3147576032877968501774576759 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, plus, x = 4^-5 = 0.0009765625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00390625);
    refVal   = 42.439212772709231546795626669;
    errBound = 42.439212772709231546795626669 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, plus, x = 4^-4 = 0.00390625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.015625);
    refVal   = 42.944441496193865255686050796;
    errBound = 42.944441496193865255686050796 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, plus, x = 4^-3 = 0.015625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0625);
    refVal   = 45.0916635710035585184703533359;
    errBound = 45.0916635710035585184703533359 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, plus, x = 4^-2 = 0.0625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9375);
    refVal   = 676.374953565053377777055300038;
    errBound = 676.374953565053377777055300038 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, plus, x = 1-4^-2 = 0.9375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.984375);
    refVal   = 2705.49981426021351110822120015;
    errBound = 2705.49981426021351110822120015 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, plus, x = 1-4^-3 = 0.984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99609375);
    refVal   = 10821.9992570408540444328848006;
    errBound = 10821.9992570408540444328848006 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, plus, x = 1-4^-4 = 0.99609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9990234375);
    refVal   = 43287.9970281634161777315392024;
    errBound = 43287.9970281634161777315392024 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, plus, x = 1-4^-5 = 0.9990234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999755859375);
    refVal   = 173151.98811265366471092615681;
    errBound = 173151.98811265366471092615681 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, plus, x = 1-4^-6 = 0.999755859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99993896484375);
    refVal   = 692607.952450614658843704627239;
    errBound = 692607.952450614658843704627239 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, plus, x = 1-4^-7 = 0.99993896484375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999847412109375);
    refVal   = 2.77043180980245863537481850895e6;
    errBound = 2.77043180980245863537481850895e6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, plus, x = 1-4^-8 = 0.9999847412109375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999996185302734375);
    refVal   = 1.10817272392098345414992740358e7;
    errBound = 1.10817272392098345414992740358e7 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, plus, x = 1-4^-9 = 0.999996185302734375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999904632568359375);
    refVal   = 4.43269089568393381659970961433e7;
    errBound = 4.43269089568393381659970961433e7 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, plus, x = 1-4^-10 = 0.99999904632568359375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999997615814208984375);
    refVal   = 1.77307635827357352663988384573e8;
    errBound = 1.77307635827357352663988384573e8 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, plus, x = 1-4^-11 = 0.9999997615814208984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999940395355224609375);
    refVal   = 7.09230543309429410655953538292e8;
    errBound = 7.09230543309429410655953538292e8 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, plus, x = 1-4^-12 = 0.999999940395355224609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999998509883880615234375);
    refVal   = 2.83692217323771764262381415317e9;
    errBound = 2.83692217323771764262381415317e9 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, plus, x = 1-4^-13 = 0.99999998509883880615234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999962747097015380859375);
    refVal   = 1.13476886929508705704952566127e10;
    errBound = 1.13476886929508705704952566127e10 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, plus, x = 1-4^-14 = 0.9999999962747097015380859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999068677425384521484375);
    refVal   = 4.53907547718034822819810264507e10;
    errBound = 4.53907547718034822819810264507e10 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullPowersOf4, plus, x = 1-4^-15 = 0.999999999068677425384521484375";
  }
}



TEST(polAggQXspace,FullNormalDoublesRegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto& eval = polAggQ_reg;

  {
    testVal  = eval(testAs,testLM,testNF,1.e-10);
    refVal   = -493371.267114392590467328436807;
    errBound = 493371.267114392590467328436807 * std::max(eps/(1.-1.e-10),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, reg, x = 10^-10 = 1.e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-9);
    refVal   = -296382.723762630380504101892328;
    errBound = 296382.723762630380504101892328 * std::max(eps/(1.-1.e-9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, reg, x = 10^-9 = 1.e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-8);
    refVal   = -166948.461217354831941387508249;
    errBound = 166948.461217354831941387508249 * std::max(eps/(1.-1.e-8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, reg, x = 10^-8 = 1.e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-7);
    refVal   = -86438.8337512770400890550284728;
    errBound = 86438.8337512770400890550284728 * std::max(eps/(1.-1.e-7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, reg, x = 10^-7 = 1.e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-6);
    refVal   = -39869.2901898580944182458292707;
    errBound = 39869.2901898580944182458292707 * std::max(eps/(1.-1.e-6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, reg, x = 10^-6 = 1.e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-5);
    refVal   = -15528.0458284036774969724228205;
    errBound = 15528.0458284036774969724228205 * std::max(eps/(1.-1.e-5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, reg, x = 10^-5 = 1.e-5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0001);
    refVal   = -4603.42694328933703437389066498;
    errBound = 4603.42694328933703437389066498 * std::max(eps/(1.-0.0001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, reg, x = 10^-4 = 0.0001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.001);
    refVal   = -811.32744086799807271638296696;
    errBound = 811.32744086799807271638296696 * std::max(eps/(1.-0.001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, reg, x = 10^-3 = 0.001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.01);
    refVal   = -30.9916298733120080065009748804;
    errBound = 30.9916298733120080065009748804 * std::max(eps/(1.-0.01),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, reg, x = 10^-2 = 0.01";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.05);
    refVal   = 6.96166039982300424784612199591;
    errBound = 6.96166039982300424784612199591 * std::max(eps/(1.-0.05),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, reg, x = 0.05 = 0.05";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.1);
    refVal   = 8.30579452697699212846350596618;
    errBound = 8.30579452697699212846350596618 * std::max(eps/(1.-0.1),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, reg, x = 0.1 = 0.1";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.15);
    refVal   = 9.37949340213819374861042952031;
    errBound = 9.37949340213819374861042952031 * std::max(eps/(1.-0.15),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, reg, x = 0.15 = 0.15";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.2);
    refVal   = 9.54449892260017978879601324649;
    errBound = 9.54449892260017978879601324649 * std::max(eps/(1.-0.2),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, reg, x = 0.2 = 0.2";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.25);
    refVal   = 8.73714879221702516879120033346;
    errBound = 8.73714879221702516879120033346 * std::max(eps/(1.-0.25),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, reg, x = 0.25 = 0.25";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.3);
    refVal   = 7.06673038723516507942681567943;
    errBound = 7.06673038723516507942681567943 * std::max(eps/(1.-0.3),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, reg, x = 0.3 = 0.3";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.35);
    refVal   = 4.6542540509292334025212474369;
    errBound = 4.6542540509292334025212474369 * std::max(eps/(1.-0.35),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, reg, x = 0.35 = 0.35";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.4);
    refVal   = 1.5971929564798132126246913212;
    errBound = 1.5971929564798132126246913212 * std::max(eps/(1.-0.4),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, reg, x = 0.4 = 0.4";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.45);
    refVal   = -2.03505731193340194377382693811;
    errBound = 2.03505731193340194377382693811 * std::max(eps/(1.-0.45),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, reg, x = 0.45 = 0.45";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.5);
    refVal   = -6.19922159217635218084781215718;
    errBound = 6.19922159217635218084781215718 * std::max(eps/(1.-0.5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, reg, x = 0.5 = 0.5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.55);
    refVal   = -10.8760669573299942218039304156;
    errBound = 10.8760669573299942218039304156 * std::max(eps/(1.-0.55),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, reg, x = 0.55 = 0.55";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.6);
    refVal   = -16.070157065326325622731284441;
    errBound = 16.070157065326325622731284441 * std::max(eps/(1.-0.6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, reg, x = 0.6 = 0.6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.65);
    refVal   = -21.8126242333173780536305285737;
    errBound = 21.8126242333173780536305285737 * std::max(eps/(1.-0.65),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, reg, x = 0.65 = 0.65";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.7);
    refVal   = -28.1688931509767087556804795917;
    errBound = 28.1688931509767087556804795917 * std::max(eps/(1.-0.7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, reg, x = 0.7 = 0.7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.75);
    refVal   = -35.2558619173097685033401472618;
    errBound = 35.2558619173097685033401472618 * std::max(eps/(1.-0.75),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, reg, x = 0.75 = 0.75";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.8);
    refVal   = -43.2804659681407129625680132543;
    errBound = 43.2804659681407129625680132543 * std::max(eps/(1.-0.8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, reg, x = 0.8 = 0.8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.85);
    refVal   = -52.636612527567525848874358992;
    errBound = 52.636612527567525848874358992 * std::max(eps/(1.-0.85),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, reg, x = 0.85 = 0.85";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9);
    refVal   = -64.2085635337221734567143698848;
    errBound = 64.2085635337221734567143698848 * std::max(eps/(1.-0.9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, reg, x = 0.9 = 0.9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.95);
    refVal   = -80.8321572107792691214242545055;
    errBound = 80.8321572107792691214242545055 * std::max(eps/(1.-0.95),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, reg, x = 0.95 = 0.95";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99);
    refVal   = -112.7567683080966970168246604819;
    errBound = 112.7567683080966970168246604819 * std::max(eps/(1.-0.99),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, reg, x = 1-10^-2 = 0.99";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999);
    refVal   = -160.359456783466444404315729131;
    errBound = 160.359456783466444404315729131 * std::max(eps/(1.-0.999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, reg, x = 1-10^-3 = 0.999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999);
    refVal   = -221.976547364453515707628731544;
    errBound = 221.976547364453515707628731544 * std::max(eps/(1.-0.9999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, reg, x = 1-10^-4 = 0.9999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999);
    refVal   = -304.689245532046187721122289812;
    errBound = 304.689245532046187721122289812 * std::max(eps/(1.-0.99999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, reg, x = 1-10^-5 = 0.99999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999);
    refVal   = -414.677275589994284908604354947;
    errBound = 414.677275589994284908604354947 * std::max(eps/(1.-0.999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, reg, x = 1-10^-6 = 0.999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999);
    refVal   = -558.047044536384721278347163754;
    errBound = 558.047044536384721278347163754 * std::max(eps/(1.-0.9999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, reg, x = 1-10^-7 = 0.9999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999999);
    refVal   = -740.902016351564835172870669732;
    errBound = 740.902016351564835172870669732 * std::max(eps/(1.-0.99999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, reg, x = 1-10^-8 = 0.99999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999);
    refVal   = -969.346063782600670312727337305;
    errBound = 969.346063782600670312727337305 * std::max(eps/(1.-0.999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, reg, x = 1-10^-9 = 0.999999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999999);
    refVal   = -1249.483193982832126072198105272;
    errBound = 1249.483193982832126072198105272 * std::max(eps/(1.-0.9999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, reg, x = 1-10^-10 = 0.9999999999";
  }
}

TEST(polAggQXspace,FullNormalDoublesPlusEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;
  auto& eval = polAggQ_plus;
  {
    testVal  = eval(testAs,testLM,testNF,1.e-10);
    refVal   = 42.2734346020431795712702742095;
    errBound = 42.2734346020431795712702742095 * std::max(eps/(1.-1.e-10),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, plus, x = 10^-10 = 1.e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-9);
    refVal   = 42.2734346400892707511552270035;
    errBound = 42.2734346400892707511552270035 * std::max(eps/(1.-1.e-9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, plus, x = 10^-9 = 1.e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-8);
    refVal   = 42.273435020550186316567819418;
    errBound = 42.273435020550186316567819418 * std::max(eps/(1.-1.e-8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, plus, x = 10^-8 = 1.e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-7);
    refVal   = 42.2734388251597186270378189561;
    errBound = 42.2734388251597186270378189561 * std::max(eps/(1.-1.e-7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, plus, x = 10^-7 = 1.e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-6);
    refVal   = 42.2734768712927074037733600257;
    errBound = 42.2734768712927074037733600257 * std::max(eps/(1.-1.e-6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, plus, x = 10^-6 = 1.e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-5);
    refVal   = 42.2738573363892000030659869122;
    errBound = 42.2738573363892000030659869122 * std::max(eps/(1.-1.e-5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, plus, x = 10^-5 = 1.e-5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0001);
    refVal   = 42.2776623640522413351994762;
    errBound = 42.2776623640522413351994762 * std::max(eps/(1.-0.0001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, plus, x = 10^-4 = 0.0001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.001);
    refVal   = 42.3157503481640001111771333857;
    errBound = 42.3157503481640001111771333857 * std::max(eps/(1.-0.001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, plus, x = 10^-3 = 0.001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.01);
    refVal   = 42.700438987692763748551470962;
    errBound = 42.700438987692763748551470962 * std::max(eps/(1.-0.01),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, plus, x = 10^-2 = 0.01";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.05);
    refVal   = 44.4983522082271959063852171077;
    errBound = 44.4983522082271959063852171077 * std::max(eps/(1.-0.05),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, plus, x = 0.05 = 0.05";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.1);
    refVal   = 46.9704828864620401234066180582;
    errBound = 46.9704828864620401234066180582 * std::max(eps/(1.-0.1),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, plus, x = 0.1 = 0.1";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.15);
    refVal   = 49.7334524680186307189011250028;
    errBound = 49.7334524680186307189011250028 * std::max(eps/(1.-0.15),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, plus, x = 0.15 = 0.15";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.2);
    refVal   = 52.8417932472697951388324453155;
    errBound = 52.8417932472697951388324453155 * std::max(eps/(1.-0.2),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, plus, x = 0.2 = 0.2";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.25);
    refVal   = 56.3645794637544481480879416698;
    errBound = 56.3645794637544481480879416698 * std::max(eps/(1.-0.25),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, plus, x = 0.25 = 0.25";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.3);
    refVal   = 60.3906208540226230158085089319;
    errBound = 60.3906208540226230158085089319 * std::max(eps/(1.-0.3),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, plus, x = 0.3 = 0.3";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.35);
    refVal   = 65.0360532274089786324091634652;
    errBound = 65.0360532274089786324091634652 * std::max(eps/(1.-0.35),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, plus, x = 0.35 = 0.35";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.4);
    refVal   = 70.4557243296930601851099270873;
    errBound = 70.4557243296930601851099270873 * std::max(eps/(1.-0.4),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, plus, x = 0.4 = 0.4";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.45);
    refVal   = 76.8607901778469747473926477316;
    errBound = 76.8607901778469747473926477316 * std::max(eps/(1.-0.45),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, plus, x = 0.45 = 0.45";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.5);
    refVal   = 84.5468691956316722221319125047;
    errBound = 84.5468691956316722221319125047 * std::max(eps/(1.-0.5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, plus, x = 0.5 = 0.5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.55);
    refVal   = 93.9409657729240802468132361164;
    errBound = 93.9409657729240802468132361164 * std::max(eps/(1.-0.55),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, plus, x = 0.55 = 0.55";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.6);
    refVal   = 105.683586494539590277664890631;
    errBound = 105.683586494539590277664890631 * std::max(eps/(1.-0.6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, plus, x = 0.6 = 0.6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.65);
    refVal   = 120.781241708045246031617017864;
    errBound = 120.781241708045246031617017864 * std::max(eps/(1.-0.65),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, plus, x = 0.65 = 0.65";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.7);
    refVal   = 140.911448659386120370219854175;
    errBound = 140.911448659386120370219854175 * std::max(eps/(1.-0.7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, plus, x = 0.7 = 0.7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.75);
    refVal   = 169.093738391263344444263825009;
    errBound = 169.093738391263344444263825009 * std::max(eps/(1.-0.75),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, plus, x = 0.75 = 0.75";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.8);
    refVal   = 211.367172989079180555329781262;
    errBound = 211.367172989079180555329781262 * std::max(eps/(1.-0.8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, plus, x = 0.8 = 0.8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.85);
    refVal   = 281.822897318772240740439708349;
    errBound = 281.822897318772240740439708349 * std::max(eps/(1.-0.85),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, plus, x = 0.85 = 0.85";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9);
    refVal   = 422.734345978158361110659562524;
    errBound = 422.734345978158361110659562524 * std::max(eps/(1.-0.9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, plus, x = 0.9 = 0.9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.95);
    refVal   = 845.468691956316722221319125047;
    errBound = 845.468691956316722221319125047 * std::max(eps/(1.-0.95),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, plus, x = 0.95 = 0.95";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99);
    refVal   = 4227.34345978158361110659562524;
    errBound = 4227.34345978158361110659562524 * std::max(eps/(1.-0.99),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, plus, x = 1-10^-2 = 0.99";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999);
    refVal   = 42273.4345978158361110659562524;
    errBound = 42273.4345978158361110659562524 * std::max(eps/(1.-0.999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, plus, x = 1-10^-3 = 0.999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999);
    refVal   = 422734.345978158361110659562524;
    errBound = 422734.345978158361110659562524 * std::max(eps/(1.-0.9999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, plus, x = 1-10^-4 = 0.9999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999);
    refVal   = 4.22734345978158361110659562524e6;
    errBound = 4.22734345978158361110659562524e6 * std::max(eps/(1.-0.99999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, plus, x = 1-10^-5 = 0.99999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999);
    refVal   = 4.22734345978158361110659562524e7;
    errBound = 4.22734345978158361110659562524e7 * std::max(eps/(1.-0.999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, plus, x = 1-10^-6 = 0.999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999);
    refVal   = 4.22734345978158361110659562524e8;
    errBound = 4.22734345978158361110659562524e8 * std::max(eps/(1.-0.9999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, plus, x = 1-10^-7 = 0.9999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999999);
    refVal   = 4.22734345978158361110659562524e9;
    errBound = 4.22734345978158361110659562524e9 * std::max(eps/(1.-0.99999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, plus, x = 1-10^-8 = 0.99999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999);
    refVal   = 4.22734345978158361110659562524e10;
    errBound = 4.22734345978158361110659562524e10 * std::max(eps/(1.-0.999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, plus, x = 1-10^-9 = 0.999999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999999);
    refVal   = 4.22734345978158361110659562524e11;
    errBound = 4.22734345978158361110659562524e11 * std::max(eps/(1.-0.9999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, FullNormalDoubles, plus, x = 1-10^-10 = 0.9999999999";
  }
}



TEST(polAggQXspace,FullDeltaEvaluation)
{
  auto& eval = polAggQ_delta;
  EXPECT_DOUBLE_EQ(eval(testAs,testLM,testNF),-13.005021117466894651792463195)
    << "^-- polAggQ, Full, delta";
}







TEST(polAggQXspace,TruncAs0NormalDoublesDeltaEvaluation)
{
  auto eval = polAggQ_delta.truncate(0);
  EXPECT_DOUBLE_EQ(eval(testAs,testLM,testNF),1.)
    << "^-- polAggQ, TruncAs0NormalDoubles, delta";
}







TEST(polAggQXspace,TruncAs1NormalDoublesDeltaEvaluation)
{
  auto eval = polAggQ_delta.truncate(1);
  EXPECT_DOUBLE_EQ(eval(testAs,testLM,testNF),0.166666666666666666666666666667)
    << "^-- polAggQ, TruncAs1NormalDoubles, delta";
}



TEST(polAggQXspace,TruncAs2PowersOf4RegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto eval = polAggQ_reg.truncate(2);

  {
    testVal  = eval(testAs,testLM,testNF,9.31322574615478515625e-10);
    refVal   = -972.505843468715920657078969741;
    errBound = 972.505843468715920657078969741 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, reg, x = 4^-15 = 9.31322574615478515625e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.7252902984619140625e-9);
    refVal   = -816.527339309269861087969985555;
    errBound = 816.527339309269861087969985555 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, reg, x = 4^-14 = 3.7252902984619140625e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.490116119384765625e-8);
    refVal   = -677.305840134214086934209845759;
    errBound = 677.305840134214086934209845759 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, reg, x = 4^-13 = 1.490116119384765625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,5.9604644775390625e-8);
    refVal   = -553.953297835093970684040661732;
    errBound = 553.953297835093970684040661732 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, reg, x = 4^-12 = 5.9604644775390625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,2.384185791015625e-7);
    refVal   = -445.581704608118546780563990328;
    errBound = 445.581704608118546780563990328 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, reg, x = 4^-11 = 2.384185791015625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,9.5367431640625e-7);
    refVal   = -351.303182196383705805108663127;
    errBound = 351.303182196383705805108663127 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, reg, x = 4^-10 = 9.5367431640625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.814697265625e-6);
    refVal   = -270.23026067192585294869238167;
    errBound = 270.23026067192585294869238167 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, reg, x = 4^-9 = 3.814697265625e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0000152587890625);
    refVal   = -201.476727014403080440681237576;
    errBound = 201.476727014403080440681237576 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, reg, x = 4^-8 = 0.0000152587890625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00006103515625);
    refVal   = -144.16012515615495602089111541;
    errBound = 144.16012515615495602089111541 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, reg, x = 4^-7 = 0.00006103515625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.000244140625);
    refVal   = -97.4088153446575654847696879831;
    errBound = 97.4088153446575654847696879831 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, reg, x = 4^-6 = 0.000244140625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0009765625);
    refVal   = -60.3807875253508291743453768088;
    errBound = 60.3807875253508291743453768088 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, reg, x = 4^-5 = 0.0009765625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00390625);
    refVal   = -32.3097516243980376991514786291;
    errBound = 32.3097516243980376991514786291 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, reg, x = 4^-4 = 0.00390625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.015625);
    refVal   = -12.6036754092522980508425973416;
    errBound = 12.6036754092522980508425973416 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, reg, x = 4^-3 = 0.015625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0625);
    refVal   = -1.00562893252612422036401115463;
    errBound = 1.00562893252612422036401115463 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, reg, x = 4^-2 = 0.0625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9375);
    refVal   = -2.11481189432922646483727742957;
    errBound = 2.11481189432922646483727742957 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, reg, x = 1-4^-2 = 0.9375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.984375);
    refVal   = -2.24499446589499602512972129335;
    errBound = 2.24499446589499602512972129335 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, reg, x = 1-4^-3 = 0.984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99609375);
    refVal   = -2.14463696184354599697038573064;
    errBound = 2.14463696184354599697038573064 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, reg, x = 1-4^-4 = 0.99609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9990234375);
    refVal   = -1.98901837553779974174050116609;
    errBound = 1.98901837553779974174050116609 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, reg, x = 1-4^-5 = 0.9990234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999755859375);
    refVal   = -1.82001819137133657502696481497;
    errBound = 1.82001819137133657502696481497 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, reg, x = 1-4^-6 = 0.999755859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99993896484375);
    refVal   = -1.64777109916164039167905431744;
    errBound = 1.64777109916164039167905431744 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, reg, x = 1-4^-7 = 0.99993896484375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999847412109375);
    refVal   = -1.47473628356226093271057663703;
    errBound = 1.47473628356226093271057663703 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, reg, x = 1-4^-8 = 0.9999847412109375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999996185302734375);
    refVal   = -1.30151049933795614052950254834;
    errBound = 1.30151049933795614052950254834 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, reg, x = 1-4^-9 = 0.999996185302734375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999904632568359375);
    refVal   = -1.12823846109567536593103512823;
    errBound = 1.12823846109567536593103512823 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, reg, x = 1-4^-10 = 0.99999904632568359375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999997615814208984375);
    refVal   = -0.954955231232383687632147313173;
    errBound = 0.954955231232383687632147313173 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, reg, x = 1-4^-11 = 0.9999997615814208984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999940395355224609375);
    refVal   = -0.781669296425267584609321856981;
    errBound = 0.781669296425267584609321856981 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, reg, x = 1-4^-12 = 0.999999940395355224609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999998509883880615234375);
    refVal   = -0.608382708621962297389414833422;
    errBound = 0.608382708621962297389414833422 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, reg, x = 1-4^-13 = 0.99999998509883880615234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999962747097015380859375);
    refVal   = -0.435095963379514560052213535526;
    errBound = 0.435095963379514560052213535526 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, reg, x = 1-4^-14 = 0.9999999962747097015380859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999068677425384521484375);
    refVal   = -0.261809180229755116383718269117;
    errBound = 0.261809180229755116383718269117 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, reg, x = 1-4^-15 = 0.999999999068677425384521484375";
  }
}

TEST(polAggQXspace,TruncAs2PowersOf4PlusEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;
  auto eval = polAggQ_plus.truncate(2);
  {
    testVal  = eval(testAs,testLM,testNF,9.31322574615478515625e-10);
    refVal   = 2.86111111377572847985368183904;
    errBound = 2.86111111377572847985368183904 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, plus, x = 4^-15 = 9.31322574615478515625e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.7252902984619140625e-9);
    refVal   = 2.86111112176958061586081383243;
    errBound = 2.86111112176958061586081383243 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, plus, x = 4^-14 = 3.7252902984619140625e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.490116119384765625e-8);
    refVal   = 2.8611111537449896065806456063;
    errBound = 2.8611111537449896065806456063 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, plus, x = 4^-13 = 1.490116119384765625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,5.9604644775390625e-8);
    refVal   = 2.86111128164663271652125284865;
    errBound = 2.86111128164663271652125284865 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, plus, x = 4^-12 = 5.9604644775390625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,2.384185791015625e-7);
    refVal   = 2.86111179325331950929100205154;
    errBound = 2.86111179325331950929100205154 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, plus, x = 4^-11 = 2.384185791015625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,9.5367431640625e-7);
    refVal   = 2.86111383968189633020474877281;
    errBound = 2.86111383968189633020474877281 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, plus, x = 4^-10 = 9.5367431640625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.814697265625e-6);
    refVal   = 2.86112202542547812114422704826;
    errBound = 2.86112202542547812114422704826 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, plus, x = 4^-9 = 3.814697265625e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0000152587890625);
    refVal   = 2.86115476886820443698447818384;
    errBound = 2.86115476886820443698447818384 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, plus, x = 4^-8 = 0.0000152587890625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00006103515625);
    refVal   = 2.86128575013394643499019986843;
    errBound = 2.86128575013394643499019986843 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, plus, x = 4^-7 = 0.00006103515625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.000244140625);
    refVal   = 2.86180979514312847646180979514;
    errBound = 2.86180979514312847646180979514 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, plus, x = 4^-6 = 0.000244140625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0009765625);
    refVal   = 2.86390789616596068208971434778;
    errBound = 2.86390789616596068208971434778 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, plus, x = 4^-5 = 0.0009765625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00390625);
    refVal   = 2.87233115468409586056644880174;
    errBound = 2.87233115468409586056644880174 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, plus, x = 4^-4 = 0.00390625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.015625);
    refVal   = 2.90652557319223985890652557319;
    errBound = 2.90652557319223985890652557319 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, plus, x = 4^-3 = 0.015625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0625);
    refVal   = 3.05185185185185185185185185185;
    errBound = 3.05185185185185185185185185185 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, plus, x = 4^-2 = 0.0625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9375);
    refVal   = 45.7777777777777777777777777778;
    errBound = 45.7777777777777777777777777778 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, plus, x = 1-4^-2 = 0.9375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.984375);
    refVal   = 183.111111111111111111111111111;
    errBound = 183.111111111111111111111111111 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, plus, x = 1-4^-3 = 0.984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99609375);
    refVal   = 732.444444444444444444444444444;
    errBound = 732.444444444444444444444444444 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, plus, x = 1-4^-4 = 0.99609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9990234375);
    refVal   = 2929.77777777777777777777777778;
    errBound = 2929.77777777777777777777777778 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, plus, x = 1-4^-5 = 0.9990234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999755859375);
    refVal   = 11719.1111111111111111111111111;
    errBound = 11719.1111111111111111111111111 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, plus, x = 1-4^-6 = 0.999755859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99993896484375);
    refVal   = 46876.4444444444444444444444444;
    errBound = 46876.4444444444444444444444444 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, plus, x = 1-4^-7 = 0.99993896484375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999847412109375);
    refVal   = 187505.777777777777777777777778;
    errBound = 187505.777777777777777777777778 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, plus, x = 1-4^-8 = 0.9999847412109375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999996185302734375);
    refVal   = 750023.111111111111111111111111;
    errBound = 750023.111111111111111111111111 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, plus, x = 1-4^-9 = 0.999996185302734375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999904632568359375);
    refVal   = 3.00009244444444444444444444444e6;
    errBound = 3.00009244444444444444444444444e6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, plus, x = 1-4^-10 = 0.99999904632568359375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999997615814208984375);
    refVal   = 1.20003697777777777777777777778e7;
    errBound = 1.20003697777777777777777777778e7 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, plus, x = 1-4^-11 = 0.9999997615814208984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999940395355224609375);
    refVal   = 4.80014791111111111111111111111e7;
    errBound = 4.80014791111111111111111111111e7 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, plus, x = 1-4^-12 = 0.999999940395355224609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999998509883880615234375);
    refVal   = 1.92005916444444444444444444444e8;
    errBound = 1.92005916444444444444444444444e8 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, plus, x = 1-4^-13 = 0.99999998509883880615234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999962747097015380859375);
    refVal   = 7.68023665777777777777777777778e8;
    errBound = 7.68023665777777777777777777778e8 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, plus, x = 1-4^-14 = 0.9999999962747097015380859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999068677425384521484375);
    refVal   = 3.07209466311111111111111111111e9;
    errBound = 3.07209466311111111111111111111e9 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2PowersOf4, plus, x = 1-4^-15 = 0.999999999068677425384521484375";
  }
}



TEST(polAggQXspace,TruncAs2NormalDoublesRegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto eval = polAggQ_reg.truncate(2);

  {
    testVal  = eval(testAs,testLM,testNF,1.e-10);
    refVal   = -1261.01396967866258400082969246;
    errBound = 1261.01396967866258400082969246 * std::max(eps/(1.-1.e-10),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, reg, x = 10^-10 = 1.e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-9);
    refVal   = -964.078481354025644214813237118;
    errBound = 964.078481354025644214813237118 * std::max(eps/(1.-1.e-9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, reg, x = 10^-9 = 1.e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-8);
    refVal   = -715.696486480093250566764290404;
    errBound = 715.696486480093250566764290404 * std::max(eps/(1.-1.e-8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, reg, x = 10^-8 = 1.e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-7);
    refVal   = -511.798679879947914600626578669;
    errBound = 511.798679879947914600626578669 * std::max(eps/(1.-1.e-7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, reg, x = 10^-7 = 1.e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-6);
    refVal   = -348.31607750424821394789251344;
    errBound = 348.31607750424821394789251344 * std::max(eps/(1.-1.e-6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, reg, x = 10^-6 = 1.e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-5);
    refVal   = -221.181877826727856415417283993;
    errBound = 221.181877826727856415417283993 * std::max(eps/(1.-1.e-5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, reg, x = 10^-5 = 1.e-5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0001);
    refVal   = -126.343163851201623865773566172;
    errBound = 126.343163851201623865773566172 * std::max(eps/(1.-0.0001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, reg, x = 10^-4 = 0.0001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.001);
    refVal   = -59.827255144599992829383595975;
    errBound = 59.827255144599992829383595975 * std::max(eps/(1.-0.001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, reg, x = 10^-3 = 0.001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.01);
    refVal   = -18.0590243513931451094983367318;
    errBound = 18.0590243513931451094983367318 * std::max(eps/(1.-0.01),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, reg, x = 10^-2 = 0.01";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.05);
    refVal   = -2.32219960070189274186255374867;
    errBound = 2.32219960070189274186255374867 * std::max(eps/(1.-0.05),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, reg, x = 0.05 = 0.05";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.1);
    refVal   = 1.06443317305285840965508823195;
    errBound = 1.06443317305285840965508823195 * std::max(eps/(1.-0.1),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, reg, x = 0.1 = 0.1";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.15);
    refVal   = 2.07453233939688527484711287793;
    errBound = 2.07453233939688527484711287793 * std::max(eps/(1.-0.15),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, reg, x = 0.15 = 0.15";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.2);
    refVal   = 2.35067059676380132013063294656;
    errBound = 2.35067059676380132013063294656 * std::max(eps/(1.-0.2),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, reg, x = 0.2 = 0.2";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.25);
    refVal   = 2.31212745611316217022930516173;
    errBound = 2.31212745611316217022930516173 * std::max(eps/(1.-0.25),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, reg, x = 0.25 = 0.25";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.3);
    refVal   = 2.11720593855296094877270960213;
    errBound = 2.11720593855296094877270960213 * std::max(eps/(1.-0.3),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, reg, x = 0.3 = 0.3";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.35);
    refVal   = 1.8386405691431567318709146995;
    errBound = 1.8386405691431567318709146995 * std::max(eps/(1.-0.35),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, reg, x = 0.35 = 0.35";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.4);
    refVal   = 1.51411424221024572329091023285;
    errBound = 1.51411424221024572329091023285 * std::max(eps/(1.-0.4),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, reg, x = 0.4 = 0.4";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.45);
    refVal   = 1.16480430141339010418168349102;
    errBound = 1.16480430141339010418168349102 * std::max(eps/(1.-0.45),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, reg, x = 0.45 = 0.45";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.5);
    refVal   = 0.803341573197952542405809919374;
    errBound = 0.803341573197952542405809919374 * std::max(eps/(1.-0.5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, reg, x = 0.5 = 0.5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.55);
    refVal   = 0.437628663913182922686094986075;
    errBound = 0.437628663913182922686094986075 * std::max(eps/(1.-0.55),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, reg, x = 0.55 = 0.55";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.6);
    refVal   = 0.0728396077103877238975076288508;
    errBound = 0.0728396077103877238975076288508 * std::max(eps/(1.-0.6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, reg, x = 0.6 = 0.6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.65);
    refVal   = -0.287444381959330564132445073872;
    errBound = 0.287444381959330564132445073872 * std::max(eps/(1.-0.65),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, reg, x = 0.65 = 0.65";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.7);
    refVal   = -0.640520543378644465308583674734;
    errBound = 0.640520543378644465308583674734 * std::max(eps/(1.-0.7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, reg, x = 0.7 = 0.7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.75);
    refVal   = -0.984026903199855884389213646719;
    errBound = 0.984026903199855884389213646719 * std::max(eps/(1.-0.75),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, reg, x = 0.75 = 0.75";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.8);
    refVal   = -1.31537699217552558188981670168;
    errBound = 1.31537699217552558188981670168 * std::max(eps/(1.-0.8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, reg, x = 0.8 = 0.8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.85);
    refVal   = -1.63081239341568478864533824065;
    errBound = 1.63081239341568478864533824065 * std::max(eps/(1.-0.85),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, reg, x = 0.85 = 0.85";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9);
    refVal   = -1.92287829741947773559544486194;
    errBound = 1.92287829741947773559544486194 * std::max(eps/(1.-0.9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, reg, x = 0.9 = 0.9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.95);
    refVal   = -2.16908618816595511716143601257;
    errBound = 2.16908618816595511716143601257 * std::max(eps/(1.-0.95),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, reg, x = 0.95 = 0.95";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99);
    refVal   = -2.22449269963815887524830359084;
    errBound = 2.22449269963815887524830359084 * std::max(eps/(1.-0.99),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, reg, x = 1-10^-2 = 0.99";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999);
    refVal   = -1.99184413579436654730452268892;
    errBound = 1.99184413579436654730452268892 * std::max(eps/(1.-0.999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, reg, x = 1-10^-3 = 0.999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999);
    refVal   = -1.70926770852138383632181685146;
    errBound = 1.70926770852138383632181685146 * std::max(eps/(1.-0.9999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, reg, x = 1-10^-4 = 0.9999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999);
    refVal   = -1.42194322179319141722188466938;
    errBound = 1.42194322179319141722188466938 * std::max(eps/(1.-0.99999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, reg, x = 1-10^-5 = 0.99999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999);
    refVal   = -1.1341673585808618712777141829;
    errBound = 1.1341673585808618712777141829 * std::max(eps/(1.-0.999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, reg, x = 1-10^-6 = 0.999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999);
    refVal   = -0.846348690245281427352011773201;
    errBound = 0.846348690245281427352011773201 * std::max(eps/(1.-0.9999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, reg, x = 1-10^-7 = 0.9999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999999);
    refVal   = -0.55852597454569400356135414244;
    errBound = 0.55852597454569400356135414244 * std::max(eps/(1.-0.99999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, reg, x = 1-10^-8 = 0.99999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999);
    refVal   = -0.270702877423495722674449896757;
    errBound = 0.270702877423495722674449896757 * std::max(eps/(1.-0.999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, reg, x = 1-10^-9 = 0.999999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999999);
    refVal   = 0.0171202555095950795224697189644;
    errBound = 0.0171202555095950795224697189644 * std::max(eps/(1.-0.9999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, reg, x = 1-10^-10 = 0.9999999999";
  }
}

TEST(polAggQXspace,TruncAs2NormalDoublesPlusEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;
  auto eval = polAggQ_plus.truncate(2);
  {
    testVal  = eval(testAs,testLM,testNF,1.e-10);
    refVal   = 2.86111111139722222225083333334;
    errBound = 2.86111111139722222225083333334 * std::max(eps/(1.-1.e-10),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, plus, x = 10^-10 = 1.e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-9);
    refVal   = 2.86111111397222222508333333619;
    errBound = 2.86111111397222222508333333619 * std::max(eps/(1.-1.e-9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, plus, x = 10^-9 = 1.e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-8);
    refVal   = 2.86111113972222250833333619444;
    errBound = 2.86111113972222250833333619444 * std::max(eps/(1.-1.e-8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, plus, x = 10^-8 = 1.e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-7);
    refVal   = 2.86111139722225083333619444473;
    errBound = 2.86111139722225083333619444473 * std::max(eps/(1.-1.e-7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, plus, x = 10^-7 = 1.e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-6);
    refVal   = 2.86111397222508333619444730556;
    errBound = 2.86111397222508333619444730556 * std::max(eps/(1.-1.e-6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, plus, x = 10^-6 = 1.e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-5);
    refVal   = 2.86113972250833619447305584167;
    errBound = 2.86113972250833619447305584167 * std::max(eps/(1.-1.e-5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, plus, x = 10^-5 = 1.e-5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0001);
    refVal   = 2.86139725083619473058416952806;
    errBound = 2.86139725083619473058416952806 * std::max(eps/(1.-0.0001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, plus, x = 10^-4 = 0.0001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.001);
    refVal   = 2.86397508619730841953064175286;
    errBound = 2.86397508619730841953064175286 * std::max(eps/(1.-0.001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, plus, x = 10^-3 = 0.001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.01);
    refVal   = 2.89001122334455667789001122334;
    errBound = 2.89001122334455667789001122334 * std::max(eps/(1.-0.01),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, plus, x = 10^-2 = 0.01";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.05);
    refVal   = 3.01169590643274853801169590643;
    errBound = 3.01169590643274853801169590643 * std::max(eps/(1.-0.05),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, plus, x = 0.05 = 0.05";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.1);
    refVal   = 3.17901234567901234567901234568;
    errBound = 3.17901234567901234567901234568 * std::max(eps/(1.-0.1),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, plus, x = 0.1 = 0.1";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.15);
    refVal   = 3.36601307189542483660130718954;
    errBound = 3.36601307189542483660130718954 * std::max(eps/(1.-0.15),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, plus, x = 0.15 = 0.15";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.2);
    refVal   = 3.57638888888888888888888888889;
    errBound = 3.57638888888888888888888888889 * std::max(eps/(1.-0.2),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, plus, x = 0.2 = 0.2";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.25);
    refVal   = 3.81481481481481481481481481481;
    errBound = 3.81481481481481481481481481481 * std::max(eps/(1.-0.25),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, plus, x = 0.25 = 0.25";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.3);
    refVal   = 4.0873015873015873015873015873;
    errBound = 4.0873015873015873015873015873 * std::max(eps/(1.-0.3),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, plus, x = 0.3 = 0.3";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.35);
    refVal   = 4.40170940170940170940170940171;
    errBound = 4.40170940170940170940170940171 * std::max(eps/(1.-0.35),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, plus, x = 0.35 = 0.35";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.4);
    refVal   = 4.76851851851851851851851851852;
    errBound = 4.76851851851851851851851851852 * std::max(eps/(1.-0.4),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, plus, x = 0.4 = 0.4";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.45);
    refVal   = 5.20202020202020202020202020202;
    errBound = 5.20202020202020202020202020202 * std::max(eps/(1.-0.45),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, plus, x = 0.45 = 0.45";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.5);
    refVal   = 5.72222222222222222222222222222;
    errBound = 5.72222222222222222222222222222 * std::max(eps/(1.-0.5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, plus, x = 0.5 = 0.5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.55);
    refVal   = 6.35802469135802469135802469136;
    errBound = 6.35802469135802469135802469136 * std::max(eps/(1.-0.55),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, plus, x = 0.55 = 0.55";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.6);
    refVal   = 7.15277777777777777777777777778;
    errBound = 7.15277777777777777777777777778 * std::max(eps/(1.-0.6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, plus, x = 0.6 = 0.6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.65);
    refVal   = 8.1746031746031746031746031746;
    errBound = 8.1746031746031746031746031746 * std::max(eps/(1.-0.65),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, plus, x = 0.65 = 0.65";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.7);
    refVal   = 9.53703703703703703703703703704;
    errBound = 9.53703703703703703703703703704 * std::max(eps/(1.-0.7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, plus, x = 0.7 = 0.7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.75);
    refVal   = 11.4444444444444444444444444444;
    errBound = 11.4444444444444444444444444444 * std::max(eps/(1.-0.75),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, plus, x = 0.75 = 0.75";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.8);
    refVal   = 14.3055555555555555555555555556;
    errBound = 14.3055555555555555555555555556 * std::max(eps/(1.-0.8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, plus, x = 0.8 = 0.8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.85);
    refVal   = 19.0740740740740740740740740741;
    errBound = 19.0740740740740740740740740741 * std::max(eps/(1.-0.85),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, plus, x = 0.85 = 0.85";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9);
    refVal   = 28.6111111111111111111111111111;
    errBound = 28.6111111111111111111111111111 * std::max(eps/(1.-0.9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, plus, x = 0.9 = 0.9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.95);
    refVal   = 57.2222222222222222222222222222;
    errBound = 57.2222222222222222222222222222 * std::max(eps/(1.-0.95),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, plus, x = 0.95 = 0.95";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99);
    refVal   = 286.111111111111111111111111111;
    errBound = 286.111111111111111111111111111 * std::max(eps/(1.-0.99),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, plus, x = 1-10^-2 = 0.99";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999);
    refVal   = 2861.11111111111111111111111111;
    errBound = 2861.11111111111111111111111111 * std::max(eps/(1.-0.999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, plus, x = 1-10^-3 = 0.999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999);
    refVal   = 28611.1111111111111111111111111;
    errBound = 28611.1111111111111111111111111 * std::max(eps/(1.-0.9999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, plus, x = 1-10^-4 = 0.9999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999);
    refVal   = 286111.111111111111111111111111;
    errBound = 286111.111111111111111111111111 * std::max(eps/(1.-0.99999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, plus, x = 1-10^-5 = 0.99999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999);
    refVal   = 2.86111111111111111111111111111e6;
    errBound = 2.86111111111111111111111111111e6 * std::max(eps/(1.-0.999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, plus, x = 1-10^-6 = 0.999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999);
    refVal   = 2.86111111111111111111111111111e7;
    errBound = 2.86111111111111111111111111111e7 * std::max(eps/(1.-0.9999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, plus, x = 1-10^-7 = 0.9999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999999);
    refVal   = 2.86111111111111111111111111111e8;
    errBound = 2.86111111111111111111111111111e8 * std::max(eps/(1.-0.99999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, plus, x = 1-10^-8 = 0.99999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999);
    refVal   = 2.86111111111111111111111111111e9;
    errBound = 2.86111111111111111111111111111e9 * std::max(eps/(1.-0.999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, plus, x = 1-10^-9 = 0.999999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999999);
    refVal   = 2.86111111111111111111111111111e10;
    errBound = 2.86111111111111111111111111111e10 * std::max(eps/(1.-0.9999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAggQ, TruncAs2NormalDoubles, plus, x = 1-10^-10 = 0.9999999999";
  }
}



TEST(polAggQXspace,TruncAs2NormalDoublesDeltaEvaluation)
{
  auto eval = polAggQ_delta.truncate(2);
  EXPECT_DOUBLE_EQ(eval(testAs,testLM,testNF),-2.99305555555555555555555555556)
    << "^-- polAggQ, TruncAs2NormalDoubles, delta";
}



TEST(polAggQNspace,FullMoments)
{
  double refVal = 0., errBound = 0.;
  auto& rpd = polAggQ;
  integration_engine_gsl engine;
  auto mom = make_mellin_moment(rpd, engine, testAs, testLM, testNF);
  std::tuple<integration_status, double, double> res;

  {
    SCOPED_TRACE("polAggQ Full Mellin moment N=3");
    refVal = -91.4049757113646348459689114366;
    errBound = 91.4049757113646348459689114366 * 3.e-11;
    res = mom.integrate(3, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAggQ Full Mellin moment N=5");
    refVal = -112.836305694573515954828305599;
    errBound = 112.836305694573515954828305599 * 3.e-11;
    res = mom.integrate(5, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAggQ Full Mellin moment N=7");
    refVal = -126.227565121148249696850037437;
    errBound = 126.227565121148249696850037437 * 3.e-11;
    res = mom.integrate(7, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

}


TEST(polAggQNspace,TruncAs0Moments)
{
  double refVal = 0., errBound = 0.;
  auto rpd = polAggQ.truncate(0);
  integration_engine_gsl engine;
  auto mom = make_mellin_moment(rpd, engine, testAs, testLM, testNF);
  std::tuple<integration_status, double, double> res;

  {
    SCOPED_TRACE("polAggQ TruncAs0 Mellin moment N=3");
    refVal = 1.;
    errBound = 1. * 3.e-11;
    res = mom.integrate(3, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAggQ TruncAs0 Mellin moment N=5");
    refVal = 1.;
    errBound = 1. * 3.e-11;
    res = mom.integrate(5, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAggQ TruncAs0 Mellin moment N=7");
    refVal = 1.;
    errBound = 1. * 3.e-11;
    res = mom.integrate(7, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

}


TEST(polAggQNspace,TruncAs1Moments)
{
  double refVal = 0., errBound = 0.;
  auto rpd = polAggQ.truncate(1);
  integration_engine_gsl engine;
  auto mom = make_mellin_moment(rpd, engine, testAs, testLM, testNF);
  std::tuple<integration_status, double, double> res;

  {
    SCOPED_TRACE("polAggQ TruncAs1 Mellin moment N=3");
    refVal = 0.166666666666666666666666666667;
    errBound = 0.166666666666666666666666666667 * 3.e-11;
    res = mom.integrate(3, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAggQ TruncAs1 Mellin moment N=5");
    refVal = 0.166666666666666666666666666667;
    errBound = 0.166666666666666666666666666667 * 3.e-11;
    res = mom.integrate(5, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAggQ TruncAs1 Mellin moment N=7");
    refVal = 0.166666666666666666666666666667;
    errBound = 0.166666666666666666666666666667 * 3.e-11;
    res = mom.integrate(7, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

}


TEST(polAggQNspace,TruncAs2Moments)
{
  double refVal = 0., errBound = 0.;
  auto rpd = polAggQ.truncate(2);
  integration_engine_gsl engine;
  auto mom = make_mellin_moment(rpd, engine, testAs, testLM, testNF);
  std::tuple<integration_status, double, double> res;

  {
    SCOPED_TRACE("polAggQ TruncAs2 Mellin moment N=3");
    refVal = -7.58802726337448559670781893004;
    errBound = 7.58802726337448559670781893004 * 3.e-11;
    res = mom.integrate(3, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAggQ TruncAs2 Mellin moment N=5");
    refVal = -9.24149886831275720164609053498;
    errBound = 9.24149886831275720164609053498 * 3.e-11;
    res = mom.integrate(5, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAggQ TruncAs2 Mellin moment N=7");
    refVal = -10.2450827079609479846360312833;
    errBound = 10.2450827079609479846360312833 * 3.e-11;
    res = mom.integrate(7, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

}
