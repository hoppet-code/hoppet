/******************************************************************************
 * This file is part of libome                                                *
 * Copyright (C) 2025 Arnd Behring, Kay Schoenwald                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 ******************************************************************************/

#include <limits>
#include <gtest/gtest.h>
#include <ome/integration_engine_gsl.h>
#include <ome/mellin.h>
#include <ome/AqqQNSOdd.h>

using namespace ome;

const double eps = std::numeric_limits<double>::epsilon();
const double testAs = 0.25;
const double testLM = -5.;
const double testNF = 3.;



TEST(AqqQNSOddXspace,FullPowersOf4RegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto& eval = AqqQNSOdd_reg;

  {
    testVal  = eval(testAs,testLM,testNF,9.31322574615478515625e-10);
    refVal   = 2903.39912854132603833359239114;
    errBound = 2903.39912854132603833359239114 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, reg, x = 4^-15 = 9.31322574615478515625e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.7252902984619140625e-9);
    refVal   = 2277.55932067424925042648430908;
    errBound = 2277.55932067424925042648430908 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, reg, x = 4^-14 = 3.7252902984619140625e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.490116119384765625e-8);
    refVal   = 1758.94199318133288314000508386;
    errBound = 1758.94199318133288314000508386 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, reg, x = 4^-13 = 1.490116119384765625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,5.9604644775390625e-8);
    refVal   = 1334.54842951154439741676594844;
    errBound = 1334.54842951154439741676594844 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, reg, x = 4^-12 = 5.9604644775390625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,2.384185791015625e-7);
    refVal   = 992.200715973558465842279637697;
    errBound = 992.200715973558465842279637697 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, reg, x = 4^-11 = 2.384185791015625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,9.5367431640625e-7);
    refVal   = 720.541843464974854585965068413;
    errBound = 720.541843464974854585965068413 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, reg, x = 4^-10 = 9.5367431640625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.814697265625e-6);
    refVal   = 509.0359781921611953114765338;
    errBound = 509.0359781921611953114765338 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, reg, x = 4^-9 = 3.814697265625e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0000152587890625);
    refVal   = 347.969144115015315870803137225;
    errBound = 347.969144115015315870803137225 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, reg, x = 4^-8 = 0.0000152587890625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00006103515625);
    refVal   = 228.450823341025795686974093328;
    errBound = 228.450823341025795686974093328 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, reg, x = 4^-7 = 0.00006103515625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.000244140625);
    refVal   = 142.417348552595966601201245589;
    errBound = 142.417348552595966601201245589 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, reg, x = 4^-6 = 0.000244140625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0009765625);
    refVal   = 82.6381059138399407898782346791;
    errBound = 82.6381059138399407898782346791 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, reg, x = 4^-5 = 0.0009765625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00390625);
    refVal   = 42.7247087363808077000289006015;
    errBound = 42.7247087363808077000289006015 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, reg, x = 4^-4 = 0.00390625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.015625);
    refVal   = 17.1426688516013205182074678324;
    errBound = 17.1426688516013205182074678324 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, reg, x = 4^-3 = 0.015625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0625);
    refVal   = 1.23623809115101078305667407077;
    errBound = 1.23623809115101078305667407077 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, reg, x = 4^-2 = 0.0625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9375);
    refVal   = -20.5676464195395071440499736196;
    errBound = 20.5676464195395071440499736196 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, reg, x = 1-4^-2 = 0.9375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.984375);
    refVal   = -25.7414153711535213547436846194;
    errBound = 25.7414153711535213547436846194 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, reg, x = 1-4^-3 = 0.984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99609375);
    refVal   = -31.4794389957750058629460685294;
    errBound = 31.4794389957750058629460685294 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, reg, x = 1-4^-4 = 0.99609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9990234375);
    refVal   = -38.0053647469800817310505419477;
    errBound = 38.0053647469800817310505419477 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, reg, x = 1-4^-5 = 0.9990234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999755859375);
    refVal   = -45.5736350210459188333431030083;
    errBound = 45.5736350210459188333431030083 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, reg, x = 1-4^-6 = 0.999755859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99993896484375);
    refVal   = -54.4467110515578175278629519678;
    errBound = 54.4467110515578175278629519678 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, reg, x = 1-4^-7 = 0.99993896484375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999847412109375);
    refVal   = -64.8882063063985092953180191294;
    errBound = 64.8882063063985092953180191294 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, reg, x = 1-4^-8 = 0.9999847412109375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999996185302734375);
    refVal   = -77.1616126420689057354583907682;
    errBound = 77.1616126420689057354583907682 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, reg, x = 1-4^-9 = 0.999996185302734375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999904632568359375);
    refVal   = -91.5302337900869068079450244269;
    errBound = 91.5302337900869068079450244269 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, reg, x = 1-4^-10 = 0.99999904632568359375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999997615814208984375);
    refVal   = -108.2572708437582300212548431737;
    errBound = 108.2572708437582300212548431737 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, reg, x = 1-4^-11 = 0.9999997615814208984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999940395355224609375);
    refVal   = -127.6058806627417640401968136727;
    errBound = 127.6058806627417640401968136727 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, reg, x = 1-4^-12 = 0.999999940395355224609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999998509883880615234375);
    refVal   = -149.839203077892859915042399361;
    errBound = 149.839203077892859915042399361 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, reg, x = 1-4^-13 = 0.99999998509883880615234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999962747097015380859375);
    refVal   = -175.220371798818152471152734319;
    errBound = 175.220371798818152471152734319 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, reg, x = 1-4^-14 = 0.9999999962747097015380859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999068677425384521484375);
    refVal   = -204.012518436007929437924463479;
    errBound = 204.012518436007929437924463479 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, reg, x = 1-4^-15 = 0.999999999068677425384521484375";
  }
}

TEST(AqqQNSOddXspace,FullPowersOf4PlusEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;
  auto& eval = AqqQNSOdd_plus;
  {
    testVal  = eval(testAs,testLM,testNF,9.31322574615478515625e-10);
    refVal   = 19.8478639549625839237424922507;
    errBound = 19.8478639549625839237424922507 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, plus, x = 4^-15 = 9.31322574615478515625e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.7252902984619140625e-9);
    refVal   = 19.8478640104168754077863490609;
    errBound = 19.8478640104168754077863490609 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, plus, x = 4^-14 = 3.7252902984619140625e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.490116119384765625e-8);
    refVal   = 19.8478642322340444427118335806;
    errBound = 19.8478642322340444427118335806 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, plus, x = 4^-13 = 1.490116119384765625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,5.9604644775390625e-8);
    refVal   = 19.8478651195027701624175971435;
    errBound = 19.8478651195027701624175971435 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, plus, x = 4^-12 = 5.9604644775390625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,2.384185791015625e-7);
    refVal   = 19.8478686685784663214880366929;
    errBound = 19.8478686685784663214880366929 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, plus, x = 4^-11 = 2.384185791015625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,9.5367431640625e-7);
    refVal   = 19.847882864893943453643334231;
    errBound = 19.847882864893943453643334231 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, plus, x = 4^-10 = 9.5367431640625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.814697265625e-6);
    refVal   = 19.8479396503589326788280243906;
    errBound = 19.8479396503589326788280243906 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, plus, x = 4^-9 = 3.814697265625e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0000152587890625);
    refVal   = 19.8481667954682295308843167689;
    errBound = 19.8481667954682295308843167689 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, plus, x = 4^-8 = 0.0000152587890625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00006103515625);
    refVal   = 19.8490754278979799534045000832;
    errBound = 19.8490754278979799534045000832 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, plus, x = 4^-7 = 0.00006103515625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.000244140625);
    refVal   = 19.8527107896979612684142811272;
    errBound = 19.8527107896979612684142811272 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, plus, x = 4^-6 = 0.000244140625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0009765625);
    refVal   = 19.8672655630041914452972827996;
    errBound = 19.8672655630041914452972827996 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, plus, x = 4^-5 = 0.0009765625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00390625);
    refVal   = 19.9256986970130273024893336313;
    errBound = 19.9256986970130273024893336313 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, plus, x = 4^-4 = 0.00390625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.015625);
    refVal   = 20.1629093957869919132332542698;
    errBound = 20.1629093957869919132332542698 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, plus, x = 4^-3 = 0.015625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0625);
    refVal   = 21.1710548655763415088949169833;
    errBound = 21.1710548655763415088949169833 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, plus, x = 4^-2 = 0.0625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9375);
    refVal   = 317.565822983645122633423754749;
    errBound = 317.565822983645122633423754749 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, plus, x = 1-4^-2 = 0.9375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.984375);
    refVal   = 1270.263291934580490533695019;
    errBound = 1270.263291934580490533695019 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, plus, x = 1-4^-3 = 0.984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99609375);
    refVal   = 5081.05316773832196213478007599;
    errBound = 5081.05316773832196213478007599 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, plus, x = 1-4^-4 = 0.99609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9990234375);
    refVal   = 20324.212670953287848539120304;
    errBound = 20324.212670953287848539120304 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, plus, x = 1-4^-5 = 0.9990234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999755859375);
    refVal   = 81296.8506838131513941564812158;
    errBound = 81296.8506838131513941564812158 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, plus, x = 1-4^-6 = 0.999755859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99993896484375);
    refVal   = 325187.402735252605576625924863;
    errBound = 325187.402735252605576625924863 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, plus, x = 1-4^-7 = 0.99993896484375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999847412109375);
    refVal   = 1.30074961094101042230650369945e6;
    errBound = 1.30074961094101042230650369945e6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, plus, x = 1-4^-8 = 0.9999847412109375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999996185302734375);
    refVal   = 5.20299844376404168922601479781e6;
    errBound = 5.20299844376404168922601479781e6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, plus, x = 1-4^-9 = 0.999996185302734375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999904632568359375);
    refVal   = 2.08119937750561667569040591912e7;
    errBound = 2.08119937750561667569040591912e7 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, plus, x = 1-4^-10 = 0.99999904632568359375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999997615814208984375);
    refVal   = 8.3247975100224667027616236765e7;
    errBound = 8.3247975100224667027616236765e7 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, plus, x = 1-4^-11 = 0.9999997615814208984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999940395355224609375);
    refVal   = 3.3299190040089866811046494706e8;
    errBound = 3.3299190040089866811046494706e8 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, plus, x = 1-4^-12 = 0.999999940395355224609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999998509883880615234375);
    refVal   = 1.33196760160359467244185978824e9;
    errBound = 1.33196760160359467244185978824e9 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, plus, x = 1-4^-13 = 0.99999998509883880615234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999962747097015380859375);
    refVal   = 5.32787040641437868976743915296e9;
    errBound = 5.32787040641437868976743915296e9 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, plus, x = 1-4^-14 = 0.9999999962747097015380859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999068677425384521484375);
    refVal   = 2.13114816256575147590697566118e10;
    errBound = 2.13114816256575147590697566118e10 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullPowersOf4, plus, x = 1-4^-15 = 0.999999999068677425384521484375";
  }
}



TEST(AqqQNSOddXspace,FullNormalDoublesRegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto& eval = AqqQNSOdd_reg;

  {
    testVal  = eval(testAs,testLM,testNF,1.e-10);
    refVal   = 4171.21399354055368668779268519;
    errBound = 4171.21399354055368668779268519 * std::max(eps/(1.-1.e-10),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, reg, x = 10^-10 = 1.e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-9);
    refVal   = 2868.45337007264971132169774003;
    errBound = 2868.45337007264971132169774003 * std::max(eps/(1.-1.e-9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, reg, x = 10^-9 = 1.e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-8);
    refVal   = 1897.94352891113413474192046475;
    errBound = 1897.94352891113413474192046475 * std::max(eps/(1.-1.e-8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, reg, x = 10^-8 = 1.e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-7);
    refVal   = 1197.79583483047081801022127129;
    errBound = 1197.79583483047081801022127129 * std::max(eps/(1.-1.e-7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, reg, x = 10^-7 = 1.e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-6);
    refVal   = 712.368772898317869010708816664;
    errBound = 712.368772898317869010708816664 * std::max(eps/(1.-1.e-6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, reg, x = 10^-6 = 1.e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-5);
    refVal   = 392.269815423841784944976818346;
    errBound = 392.269815423841784944976818346 * std::max(eps/(1.-1.e-5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, reg, x = 10^-5 = 1.e-5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0001);
    refVal   = 194.363982236635406886863303647;
    errBound = 194.363982236635406886863303647 * std::max(eps/(1.-0.0001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, reg, x = 10^-4 = 0.0001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.001);
    refVal   = 81.8052394393326426541235667784;
    errBound = 81.8052394393326426541235667784 * std::max(eps/(1.-0.001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, reg, x = 10^-3 = 0.001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.01);
    refVal   = 24.1165115027787853353279442653;
    errBound = 24.1165115027787853353279442653 * std::max(eps/(1.-0.01),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, reg, x = 10^-2 = 0.01";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.05);
    refVal   = 3.30890121703068430854414856296;
    errBound = 3.30890121703068430854414856296 * std::max(eps/(1.-0.05),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, reg, x = 0.05 = 0.05";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.1);
    refVal   = -2.64779113163787358652265882781;
    errBound = 2.64779113163787358652265882781 * std::max(eps/(1.-0.1),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, reg, x = 0.1 = 0.1";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.15);
    refVal   = -5.55052599638825291633857410323;
    errBound = 5.55052599638825291633857410323 * std::max(eps/(1.-0.15),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, reg, x = 0.15 = 0.15";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.2);
    refVal   = -7.40390561521698636512059202837;
    errBound = 7.40390561521698636512059202837 * std::max(eps/(1.-0.2),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, reg, x = 0.2 = 0.2";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.25);
    refVal   = -8.7463459808140937283833769944;
    errBound = 8.7463459808140937283833769944 * std::max(eps/(1.-0.25),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, reg, x = 0.25 = 0.25";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.3);
    refVal   = -9.79561044867764414986540450641;
    errBound = 9.79561044867764414986540450641 * std::max(eps/(1.-0.3),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, reg, x = 0.3 = 0.3";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.35);
    refVal   = -10.6605952576011577538268406022;
    errBound = 10.6605952576011577538268406022 * std::max(eps/(1.-0.35),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, reg, x = 0.35 = 0.35";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.4);
    refVal   = -11.4039510193060619091726791845;
    errBound = 11.4039510193060619091726791845 * std::max(eps/(1.-0.4),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, reg, x = 0.4 = 0.4";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.45);
    refVal   = -12.0659151371789954763226083221;
    errBound = 12.0659151371789954763226083221 * std::max(eps/(1.-0.45),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, reg, x = 0.45 = 0.45";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.5);
    refVal   = -12.6750864156939643451842303613;
    errBound = 12.6750864156939643451842303613 * std::max(eps/(1.-0.5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, reg, x = 0.5 = 0.5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.55);
    refVal   = -13.25407525365028954431777913637;
    errBound = 13.25407525365028954431777913637 * std::max(eps/(1.-0.55),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, reg, x = 0.55 = 0.55";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.6);
    refVal   = -13.82302851608448007961445921745;
    errBound = 13.82302851608448007961445921745 * std::max(eps/(1.-0.6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, reg, x = 0.6 = 0.6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.65);
    refVal   = -14.4024522829480081103051209236;
    errBound = 14.4024522829480081103051209236 * std::max(eps/(1.-0.65),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, reg, x = 0.65 = 0.65";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.7);
    refVal   = -15.016345521963826093858542078;
    errBound = 15.016345521963826093858542078 * std::max(eps/(1.-0.7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, reg, x = 0.7 = 0.7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.75);
    refVal   = -15.6969771965170546069472478931;
    errBound = 15.6969771965170546069472478931 * std::max(eps/(1.-0.75),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, reg, x = 0.75 = 0.75";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.8);
    refVal   = -16.4942070805109232025013435174;
    errBound = 16.4942070805109232025013435174 * std::max(eps/(1.-0.8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, reg, x = 0.8 = 0.8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.85);
    refVal   = -17.497947791731271498149373258;
    errBound = 17.497947791731271498149373258 * std::max(eps/(1.-0.85),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, reg, x = 0.85 = 0.85";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9);
    refVal   = -18.9079442394228220141994118981;
    errBound = 18.9079442394228220141994118981 * std::max(eps/(1.-0.9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, reg, x = 0.9 = 0.9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.95);
    refVal   = -21.370040255524012286612827796;
    errBound = 21.370040255524012286612827796 * std::max(eps/(1.-0.95),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, reg, x = 0.95 = 0.95";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99);
    refVal   = -27.5173148761878782376964250776;
    errBound = 27.5173148761878782376964250776 * std::max(eps/(1.-0.99),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, reg, x = 1-10^-2 = 0.99";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999);
    refVal   = -37.8857000524258674484656762278;
    errBound = 37.8857000524258674484656762278 * std::max(eps/(1.-0.999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, reg, x = 1-10^-3 = 0.999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999);
    refVal   = -51.1204597796980893194303979581;
    errBound = 51.1204597796980893194303979581 * std::max(eps/(1.-0.9999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, reg, x = 1-10^-4 = 0.9999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999);
    refVal   = -68.4231568096435914195898626223;
    errBound = 68.4231568096435914195898626223 * std::max(eps/(1.-0.99999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, reg, x = 1-10^-5 = 0.99999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999);
    refVal   = -91.0011332217618289498403785141;
    errBound = 91.0011332217618289498403785141 * std::max(eps/(1.-0.999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, reg, x = 1-10^-6 = 0.999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999);
    refVal   = -120.0606511645336770101492909084;
    errBound = 120.0606511645336770101492909084 * std::max(eps/(1.-0.9999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, reg, x = 1-10^-7 = 0.9999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999999);
    refVal   = -156.807555960290542000615305058;
    errBound = 156.807555960290542000615305058 * std::max(eps/(1.-0.99999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, reg, x = 1-10^-8 = 0.99999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999);
    refVal   = -202.447602274376662346641936703;
    errBound = 202.447602274376662346641936703 * std::max(eps/(1.-0.999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, reg, x = 1-10^-9 = 0.999999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999999);
    refVal   = -258.186528604711704204407175938;
    errBound = 258.186528604711704204407175938 * std::max(eps/(1.-0.9999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, reg, x = 1-10^-10 = 0.9999999999";
  }
}

TEST(AqqQNSOddXspace,FullNormalDoublesPlusEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;
  auto& eval = AqqQNSOdd_plus;
  {
    testVal  = eval(testAs,testLM,testNF,1.e-10);
    refVal   = 19.8478639384626065584352453277;
    errBound = 19.8478639384626065584352453277 * std::max(eps/(1.-1.e-10),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, plus, x = 10^-10 = 1.e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-9);
    refVal   = 19.8478639563256841209146687927;
    errBound = 19.8478639563256841209146687927 * std::max(eps/(1.-1.e-9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, plus, x = 10^-9 = 1.e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-8);
    refVal   = 19.8478641349564615141535998134;
    errBound = 19.8478641349564615141535998134 * std::max(eps/(1.-1.e-8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, plus, x = 10^-8 = 1.e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-7);
    refVal   = 19.8478659212644122910302137749;
    errBound = 19.8478659212644122910302137749 * std::max(eps/(1.-1.e-7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, plus, x = 10^-7 = 1.e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-6);
    refVal   = 19.8478837843616045261935108653;
    errBound = 19.8478837843616045261935108653 * std::max(eps/(1.-1.e-6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, plus, x = 10^-6 = 1.e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-5);
    refVal   = 19.8480624171019911845008296801;
    errBound = 19.8480624171019911845008296801 * std::max(eps/(1.-1.e-5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, plus, x = 10^-5 = 1.e-5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0001);
    refVal   = 19.8498489213699571603050151733;
    errBound = 19.8498489213699571603050151733 * std::max(eps/(1.-0.0001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, plus, x = 10^-4 = 0.0001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.001);
    refVal   = 19.8677316681459661307197043762;
    errBound = 19.8677316681459661307197043762 * std::max(eps/(1.-0.001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, plus, x = 10^-3 = 0.001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.01);
    refVal   = 20.0483474105836567319080653251;
    errBound = 20.0483474105836567319080653251 * std::max(eps/(1.-0.01),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, plus, x = 10^-2 = 0.01";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.05);
    refVal   = 20.8924883541871791206199838651;
    errBound = 20.8924883541871791206199838651 * std::max(eps/(1.-0.05),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, plus, x = 0.05 = 0.05";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.1);
    refVal   = 22.0531821516420224050988718576;
    errBound = 22.0531821516420224050988718576 * std::max(eps/(1.-0.1),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, plus, x = 0.1 = 0.1";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.15);
    refVal   = 23.3504281605621413701046878492;
    errBound = 23.3504281605621413701046878492 * std::max(eps/(1.-0.15),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, plus, x = 0.15 = 0.15";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.2);
    refVal   = 24.8098299205972752057362308398;
    errBound = 24.8098299205972752057362308398 * std::max(eps/(1.-0.2),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, plus, x = 0.2 = 0.2";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.25);
    refVal   = 26.4638185819704268861186462291;
    errBound = 26.4638185819704268861186462291 * std::max(eps/(1.-0.25),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, plus, x = 0.25 = 0.25";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.3);
    refVal   = 28.3540913378254573779842638169;
    errBound = 28.3540913378254573779842638169 * std::max(eps/(1.-0.3),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, plus, x = 0.3 = 0.3";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.35);
    refVal   = 30.5351752868889540993676687259;
    errBound = 30.5351752868889540993676687259 * std::max(eps/(1.-0.35),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, plus, x = 0.35 = 0.35";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.4);
    refVal   = 33.0797732274630336076483077864;
    errBound = 33.0797732274630336076483077864 * std::max(eps/(1.-0.4),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, plus, x = 0.4 = 0.4";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.45);
    refVal   = 36.0870253390505821174345175851;
    errBound = 36.0870253390505821174345175851 * std::max(eps/(1.-0.45),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, plus, x = 0.45 = 0.45";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.5);
    refVal   = 39.6957278729556403291779693437;
    errBound = 39.6957278729556403291779693437 * std::max(eps/(1.-0.5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, plus, x = 0.5 = 0.5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.55);
    refVal   = 44.1063643032840448101977437152;
    errBound = 44.1063643032840448101977437152 * std::max(eps/(1.-0.55),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, plus, x = 0.55 = 0.55";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.6);
    refVal   = 49.6196598411945504114724616796;
    errBound = 49.6196598411945504114724616796 * std::max(eps/(1.-0.6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, plus, x = 0.6 = 0.6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.65);
    refVal   = 56.7081826756509147559685276338;
    errBound = 56.7081826756509147559685276338 * std::max(eps/(1.-0.65),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, plus, x = 0.65 = 0.65";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.7);
    refVal   = 66.1595464549260672152966155728;
    errBound = 66.1595464549260672152966155728 * std::max(eps/(1.-0.7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, plus, x = 0.7 = 0.7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.75);
    refVal   = 79.3914557459112806583559386873;
    errBound = 79.3914557459112806583559386873 * std::max(eps/(1.-0.75),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, plus, x = 0.75 = 0.75";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.8);
    refVal   = 99.2393196823891008229449233592;
    errBound = 99.2393196823891008229449233592 * std::max(eps/(1.-0.8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, plus, x = 0.8 = 0.8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.85);
    refVal   = 132.319092909852134430593231146;
    errBound = 132.319092909852134430593231146 * std::max(eps/(1.-0.85),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, plus, x = 0.85 = 0.85";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9);
    refVal   = 198.478639364778201645889846718;
    errBound = 198.478639364778201645889846718 * std::max(eps/(1.-0.9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, plus, x = 0.9 = 0.9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.95);
    refVal   = 396.957278729556403291779693437;
    errBound = 396.957278729556403291779693437 * std::max(eps/(1.-0.95),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, plus, x = 0.95 = 0.95";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99);
    refVal   = 1984.78639364778201645889846718;
    errBound = 1984.78639364778201645889846718 * std::max(eps/(1.-0.99),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, plus, x = 1-10^-2 = 0.99";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999);
    refVal   = 19847.8639364778201645889846718;
    errBound = 19847.8639364778201645889846718 * std::max(eps/(1.-0.999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, plus, x = 1-10^-3 = 0.999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999);
    refVal   = 198478.639364778201645889846718;
    errBound = 198478.639364778201645889846718 * std::max(eps/(1.-0.9999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, plus, x = 1-10^-4 = 0.9999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999);
    refVal   = 1.98478639364778201645889846718e6;
    errBound = 1.98478639364778201645889846718e6 * std::max(eps/(1.-0.99999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, plus, x = 1-10^-5 = 0.99999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999);
    refVal   = 1.98478639364778201645889846718e7;
    errBound = 1.98478639364778201645889846718e7 * std::max(eps/(1.-0.999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, plus, x = 1-10^-6 = 0.999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999);
    refVal   = 1.98478639364778201645889846718e8;
    errBound = 1.98478639364778201645889846718e8 * std::max(eps/(1.-0.9999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, plus, x = 1-10^-7 = 0.9999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999999);
    refVal   = 1.98478639364778201645889846718e9;
    errBound = 1.98478639364778201645889846718e9 * std::max(eps/(1.-0.99999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, plus, x = 1-10^-8 = 0.99999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999);
    refVal   = 1.98478639364778201645889846718e10;
    errBound = 1.98478639364778201645889846718e10 * std::max(eps/(1.-0.999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, plus, x = 1-10^-9 = 0.999999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999999);
    refVal   = 1.98478639364778201645889846718e11;
    errBound = 1.98478639364778201645889846718e11 * std::max(eps/(1.-0.9999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, FullNormalDoubles, plus, x = 1-10^-10 = 0.9999999999";
  }
}



TEST(AqqQNSOddXspace,FullDeltaEvaluation)
{
  auto& eval = AqqQNSOdd_delta;
  EXPECT_DOUBLE_EQ(eval(testAs,testLM,testNF),12.1501748418929015300593526677)
    << "^-- AqqQNSOdd, Full, delta";
}







TEST(AqqQNSOddXspace,TruncAs0NormalDoublesDeltaEvaluation)
{
  auto eval = AqqQNSOdd_delta.truncate(0);
  EXPECT_DOUBLE_EQ(eval(testAs,testLM,testNF),1.)
    << "^-- AqqQNSOdd, TruncAs0NormalDoubles, delta";
}







TEST(AqqQNSOddXspace,TruncAs1NormalDoublesDeltaEvaluation)
{
  auto eval = AqqQNSOdd_delta.truncate(1);
  EXPECT_DOUBLE_EQ(eval(testAs,testLM,testNF),1.)
    << "^-- AqqQNSOdd, TruncAs1NormalDoubles, delta";
}



TEST(AqqQNSOddXspace,TruncAs2PowersOf4RegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto eval = AqqQNSOdd_reg.truncate(2);

  {
    testVal  = eval(testAs,testLM,testNF,9.31322574615478515625e-10);
    refVal   = 17.8217061035486357547375500369;
    errBound = 17.8217061035486357547375500369 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, reg, x = 4^-15 = 9.31322574615478515625e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.7252902984619140625e-9);
    refVal   = 15.7858095460317299789869394768;
    errBound = 15.7858095460317299789869394768 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, reg, x = 4^-14 = 3.7252902984619140625e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.490116119384765625e-8);
    refVal   = 13.856680476843573151637347734;
    errBound = 13.856680476843573151637347734 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, reg, x = 4^-13 = 1.490116119384765625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,5.9604644775390625e-8);
    refVal   = 12.0343192826499716016331048158;
    errBound = 12.0343192826499716016331048158 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, reg, x = 4^-12 = 5.9604644775390625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,2.384185791015625e-7);
    refVal   = 10.3187273173085245697509487745;
    errBound = 10.3187273173085245697509487745 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, reg, x = 4^-11 = 2.384185791015625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,9.5367431640625e-7);
    refVal   = 8.70990926799161824778922135492;
    errBound = 8.70990926799161824778922135492 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, reg, x = 4^-10 = 9.5367431640625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.814697265625e-6);
    refVal   = 7.20788114229739815587154947429;
    errBound = 7.20788114229739815587154947429 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, reg, x = 4^-9 = 3.814697265625e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0000152587890625);
    refVal   = 5.81269669510812129899603649672;
    errBound = 5.81269669510812129899603649672 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, reg, x = 4^-8 = 0.0000152587890625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00006103515625);
    refVal   = 4.52453261341313405656388449925;
    errBound = 4.52453261341313405656388449925 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, reg, x = 4^-7 = 0.00006103515625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.000244140625);
    refVal   = 3.34395356834756400818569948679;
    errBound = 3.34395356834756400818569948679 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, reg, x = 4^-6 = 0.000244140625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0009765625);
    refVal   = 2.27269703009075740014577399795;
    errBound = 2.27269703009075740014577399795 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, reg, x = 4^-5 = 0.0009765625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00390625);
    refVal   = 1.31583931756388828143992875279;
    errBound = 1.31583931756388828143992875279 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, reg, x = 4^-4 = 0.00390625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.015625);
    refVal   = 0.48726127075627485257793501304;
    errBound = 0.48726127075627485257793501304 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, reg, x = 4^-3 = 0.015625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0625);
    refVal   = -0.177135651603830542552739162964;
    errBound = 0.177135651603830542552739162964 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, reg, x = 4^-2 = 0.0625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9375);
    refVal   = -0.38499958453489546857703568513;
    errBound = 0.38499958453489546857703568513 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, reg, x = 1-4^-2 = 0.9375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.984375);
    refVal   = -0.3556611804860762279683890258;
    errBound = 0.3556611804860762279683890258 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, reg, x = 1-4^-3 = 0.984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99609375);
    refVal   = -0.34818370702466657859671556732;
    errBound = 0.34818370702466657859671556732 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, reg, x = 1-4^-4 = 0.99609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9990234375);
    refVal   = -0.346305752864286253367270487748;
    errBound = 0.346305752864286253367270487748 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, reg, x = 1-4^-5 = 0.9990234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999755859375);
    refVal   = -0.345835732823234028838374314547;
    errBound = 0.345835732823234028838374314547 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, reg, x = 1-4^-6 = 0.999755859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99993896484375);
    refVal   = -0.345718194673063083996800636007;
    errBound = 0.345718194673063083996800636007 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, reg, x = 1-4^-7 = 0.99993896484375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999847412109375);
    refVal   = -0.345688808065505242198483338781;
    errBound = 0.345688808065505242198483338781 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, reg, x = 1-4^-8 = 0.9999847412109375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999996185302734375);
    refVal   = -0.345681461284259028495949349392;
    errBound = 0.345681461284259028495949349392 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, reg, x = 1-4^-9 = 0.999996185302734375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999904632568359375);
    refVal   = -0.345679624580862977792671971337;
    errBound = 0.345679624580862977792671971337 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, reg, x = 1-4^-10 = 0.99999904632568359375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999997615814208984375);
    refVal   = -0.345679165404508688721165380759;
    errBound = 0.345679165404508688721165380759 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, reg, x = 1-4^-11 = 0.9999997615814208984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999940395355224609375);
    refVal   = -0.345679050610388536751747507152;
    errBound = 0.345679050610388536751747507152 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, reg, x = 1-4^-12 = 0.999999940395355224609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999998509883880615234375);
    refVal   = -0.345679021911856525029190290443;
    errBound = 0.345679021911856525029190290443 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, reg, x = 1-4^-13 = 0.99999998509883880615234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999962747097015380859375);
    refVal   = -0.345679014737223398740431182894;
    errBound = 0.345679014737223398740431182894 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, reg, x = 1-4^-14 = 0.9999999962747097015380859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999068677425384521484375);
    refVal   = -0.34567901294356510945835919749;
    errBound = 0.34567901294356510945835919749 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, reg, x = 1-4^-15 = 0.999999999068677425384521484375";
  }
}

TEST(AqqQNSOddXspace,TruncAs2PowersOf4PlusEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;
  auto eval = AqqQNSOdd_plus.truncate(2);
  {
    testVal  = eval(testAs,testLM,testNF,9.31322574615478515625e-10);
    refVal   = 1.27160493945587932437941415068;
    errBound = 1.27160493945587932437941415068 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, plus, x = 4^-15 = 9.31322574615478515625e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.7252902984619140625e-9);
    refVal   = 1.27160494300870249593813948108;
    errBound = 1.27160494300870249593813948108 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, plus, x = 4^-14 = 3.7252902984619140625e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.490116119384765625e-8);
    refVal   = 1.27160495721999538070250915835;
    errBound = 1.27160495721999538070250915835 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, plus, x = 4^-13 = 1.490116119384765625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,5.9604644775390625e-8);
    refVal   = 1.27160501406517009623166793273;
    errBound = 1.27160501406517009623166793273 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, plus, x = 4^-12 = 5.9604644775390625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,2.384185791015625e-7);
    refVal   = 1.27160524144591978190711202291;
    errBound = 1.27160524144591978190711202291 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, plus, x = 4^-11 = 2.384185791015625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,9.5367431640625e-7);
    refVal   = 1.2716061509697317023132216768;
    errBound = 1.2716061509697317023132216768 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, plus, x = 4^-10 = 9.5367431640625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.814697265625e-6);
    refVal   = 1.27160978907799027606410091034;
    errBound = 1.27160978907799027606410091034 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, plus, x = 4^-9 = 3.814697265625e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0000152587890625);
    refVal   = 1.27162434171920197199310141504;
    errBound = 1.27162434171920197199310141504 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, plus, x = 4^-8 = 0.0000152587890625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00006103515625);
    refVal   = 1.27168255561508730444008883041;
    errBound = 1.27168255561508730444008883041 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, plus, x = 4^-7 = 0.00006103515625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.000244140625);
    refVal   = 1.27191546450805710064969324229;
    errBound = 1.27191546450805710064969324229 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, plus, x = 4^-6 = 0.000244140625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0009765625);
    refVal   = 1.27284795385153808092876193235;
    errBound = 1.27284795385153808092876193235 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, plus, x = 4^-5 = 0.0009765625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00390625);
    refVal   = 1.27659162430404260469619946744;
    errBound = 1.27659162430404260469619946744 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, plus, x = 4^-4 = 0.00390625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.015625);
    refVal   = 1.2917891436409954928473446992;
    errBound = 1.2917891436409954928473446992 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, plus, x = 4^-3 = 0.015625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0625);
    refVal   = 1.35637860082304526748971193416;
    errBound = 1.35637860082304526748971193416 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, plus, x = 4^-2 = 0.0625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9375);
    refVal   = 20.3456790123456790123456790123;
    errBound = 20.3456790123456790123456790123 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, plus, x = 1-4^-2 = 0.9375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.984375);
    refVal   = 81.3827160493827160493827160494;
    errBound = 81.3827160493827160493827160494 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, plus, x = 1-4^-3 = 0.984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99609375);
    refVal   = 325.530864197530864197530864198;
    errBound = 325.530864197530864197530864198 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, plus, x = 1-4^-4 = 0.99609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9990234375);
    refVal   = 1302.12345679012345679012345679;
    errBound = 1302.12345679012345679012345679 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, plus, x = 1-4^-5 = 0.9990234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999755859375);
    refVal   = 5208.49382716049382716049382716;
    errBound = 5208.49382716049382716049382716 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, plus, x = 1-4^-6 = 0.999755859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99993896484375);
    refVal   = 20833.9753086419753086419753086;
    errBound = 20833.9753086419753086419753086 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, plus, x = 1-4^-7 = 0.99993896484375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999847412109375);
    refVal   = 83335.9012345679012345679012346;
    errBound = 83335.9012345679012345679012346 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, plus, x = 1-4^-8 = 0.9999847412109375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999996185302734375);
    refVal   = 333343.604938271604938271604938;
    errBound = 333343.604938271604938271604938 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, plus, x = 1-4^-9 = 0.999996185302734375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999904632568359375);
    refVal   = 1.33337441975308641975308641975e6;
    errBound = 1.33337441975308641975308641975e6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, plus, x = 1-4^-10 = 0.99999904632568359375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999997615814208984375);
    refVal   = 5.33349767901234567901234567901e6;
    errBound = 5.33349767901234567901234567901e6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, plus, x = 1-4^-11 = 0.9999997615814208984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999940395355224609375);
    refVal   = 2.1333990716049382716049382716e7;
    errBound = 2.1333990716049382716049382716e7 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, plus, x = 1-4^-12 = 0.999999940395355224609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999998509883880615234375);
    refVal   = 8.53359628641975308641975308642e7;
    errBound = 8.53359628641975308641975308642e7 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, plus, x = 1-4^-13 = 0.99999998509883880615234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999962747097015380859375);
    refVal   = 3.41343851456790123456790123457e8;
    errBound = 3.41343851456790123456790123457e8 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, plus, x = 1-4^-14 = 0.9999999962747097015380859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999068677425384521484375);
    refVal   = 1.36537540582716049382716049383e9;
    errBound = 1.36537540582716049382716049383e9 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2PowersOf4, plus, x = 1-4^-15 = 0.999999999068677425384521484375";
  }
}



TEST(AqqQNSOddXspace,TruncAs2NormalDoublesRegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto eval = AqqQNSOdd_reg.truncate(2);

  {
    testVal  = eval(testAs,testLM,testNF,1.e-10);
    refVal   = 21.3230102048050060083039198921;
    errBound = 21.3230102048050060083039198921 * std::max(eps/(1.-1.e-10),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, reg, x = 10^-10 = 1.e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-9);
    refVal   = 17.7146173942785623569748329156;
    errBound = 17.7146173942785623569748329156 * std::max(eps/(1.-1.e-9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, reg, x = 10^-9 = 1.e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-8);
    refVal   = 14.4007746358257017944876439667;
    errBound = 14.4007746358257017944876439667 * std::max(eps/(1.-1.e-8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, reg, x = 10^-8 = 1.e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-7);
    refVal   = 11.3814830621623545539010845103;
    errBound = 11.3814830621623545539010845103 * std::max(eps/(1.-1.e-7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, reg, x = 10^-7 = 1.e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-6);
    refVal   = 8.65675171268660515546140388389;
    errBound = 8.65675171268660515546140388389 * std::max(eps/(1.-1.e-6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, reg, x = 10^-6 = 1.e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-5);
    refVal   = 6.22665025234672232098097860793;
    errBound = 6.22665025234672232098097860793 * std::max(eps/(1.-1.e-5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, reg, x = 10^-5 = 1.e-5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0001);
    refVal   = 4.09168960766253686199028926281;
    errBound = 4.09168960766253686199028926281 * std::max(eps/(1.-0.0001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, reg, x = 10^-4 = 0.0001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.001);
    refVal   = 2.25534283688209607561653854713;
    errBound = 2.25534283688209607561653854713 * std::max(eps/(1.-0.001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, reg, x = 10^-3 = 0.001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.01);
    refVal   = 0.738488095420545443932353103742;
    errBound = 0.738488095420545443932353103742 * std::max(eps/(1.-0.01),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, reg, x = 10^-2 = 0.01";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.05);
    refVal   = -0.0840366145940149601139494039936;
    errBound = 0.0840366145940149601139494039936 * std::max(eps/(1.-0.05),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, reg, x = 0.05 = 0.05";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.1);
    refVal   = -0.351397413919272335738111721405;
    errBound = 0.351397413919272335738111721405 * std::max(eps/(1.-0.1),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, reg, x = 0.1 = 0.1";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.15);
    refVal   = -0.473358476701255463509266244707;
    errBound = 0.473358476701255463509266244707 * std::max(eps/(1.-0.15),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, reg, x = 0.15 = 0.15";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.2);
    refVal   = -0.540138954765073360749215426109;
    errBound = 0.540138954765073360749215426109 * std::max(eps/(1.-0.2),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, reg, x = 0.2 = 0.2";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.25);
    refVal   = -0.578207692291456844284520537334;
    errBound = 0.578207692291456844284520537334 * std::max(eps/(1.-0.25),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, reg, x = 0.25 = 0.25";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.3);
    refVal   = -0.598802628849339143129449818016;
    errBound = 0.598802628849339143129449818016 * std::max(eps/(1.-0.3),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, reg, x = 0.3 = 0.3";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.35);
    refVal   = -0.607691174761393407886918353833;
    errBound = 0.607691174761393407886918353833 * std::max(eps/(1.-0.35),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, reg, x = 0.35 = 0.35";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.4);
    refVal   = -0.608205938556899136474781397016;
    errBound = 0.608205938556899136474781397016 * std::max(eps/(1.-0.4),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, reg, x = 0.4 = 0.4";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.45);
    refVal   = -0.602437372263944781429130689028;
    errBound = 0.602437372263944781429130689028 * std::max(eps/(1.-0.45),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, reg, x = 0.45 = 0.45";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.5);
    refVal   = -0.591778591749814791614315820361;
    errBound = 0.591778591749814791614315820361 * std::max(eps/(1.-0.5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, reg, x = 0.5 = 0.5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.55);
    refVal   = -0.577202330444190885352981519144;
    errBound = 0.577202330444190885352981519144 * std::max(eps/(1.-0.55),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, reg, x = 0.55 = 0.55";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.6);
    refVal   = -0.559413489610987856045682073043;
    errBound = 0.559413489610987856045682073043 * std::max(eps/(1.-0.6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, reg, x = 0.6 = 0.6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.65);
    refVal   = -0.538938605094666549378838062134;
    errBound = 0.538938605094666549378838062134 * std::max(eps/(1.-0.65),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, reg, x = 0.65 = 0.65";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.7);
    refVal   = -0.516181035683198740514234446871;
    errBound = 0.516181035683198740514234446871 * std::max(eps/(1.-0.7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, reg, x = 0.7 = 0.7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.75);
    refVal   = -0.491456453381019672887120740997;
    errBound = 0.491456453381019672887120740997 * std::max(eps/(1.-0.75),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, reg, x = 0.75 = 0.75";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.8);
    refVal   = -0.465016477683309411427140319021;
    errBound = 0.465016477683309411427140319021 * std::max(eps/(1.-0.8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, reg, x = 0.8 = 0.8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.85);
    refVal   = -0.43706489002058434918027274715;
    errBound = 0.43706489002058434918027274715 * std::max(eps/(1.-0.85),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, reg, x = 0.85 = 0.85";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9);
    refVal   = -0.407769047200208012914042224406;
    errBound = 0.407769047200208012914042224406 * std::max(eps/(1.-0.9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, reg, x = 0.9 = 0.9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.95);
    refVal   = -0.377268097411200815373475106491;
    errBound = 0.377268097411200815373475106491 * std::max(eps/(1.-0.95),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, reg, x = 0.95 = 0.95";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99);
    refVal   = -0.352078889482817115652856466876;
    errBound = 0.352078889482817115652856466876 * std::max(eps/(1.-0.99),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, reg, x = 1-10^-2 = 0.99";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999);
    refVal   = -0.346320790001448995456553329598;
    errBound = 0.346320790001448995456553329598 * std::max(eps/(1.-0.999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, reg, x = 1-10^-3 = 0.999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999);
    refVal   = -0.345743207901112644906621845497;
    errBound = 0.345743207901112644906621845497 * std::max(eps/(1.-0.9999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, reg, x = 1-10^-4 = 0.9999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999);
    refVal   = -0.345685432079012223764490732886;
    errBound = 0.345685432079012223764490732886 * std::max(eps/(1.-0.99999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, reg, x = 1-10^-5 = 0.99999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999);
    refVal   = -0.345679654320790123334876449074;
    errBound = 0.345679654320790123334876449074 * std::max(eps/(1.-0.999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, reg, x = 1-10^-6 = 0.999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999);
    refVal   = -0.345679076543207901234445987645;
    errBound = 0.345679076543207901234445987645 * std::max(eps/(1.-0.9999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, reg, x = 1-10^-7 = 0.9999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999999);
    refVal   = -0.345679018765432079012345557099;
    errBound = 0.345679018765432079012345557099 * std::max(eps/(1.-0.99999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, reg, x = 1-10^-8 = 0.99999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999);
    refVal   = -0.345679012987654320790123456668;
    errBound = 0.345679012987654320790123456668 * std::max(eps/(1.-0.999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, reg, x = 1-10^-9 = 0.999999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999999);
    refVal   = -0.345679012409876543207901234568;
    errBound = 0.345679012409876543207901234568 * std::max(eps/(1.-0.9999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, reg, x = 1-10^-10 = 0.9999999999";
  }
}

TEST(AqqQNSOddXspace,TruncAs2NormalDoublesPlusEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;
  auto eval = AqqQNSOdd_plus.truncate(2);
  {
    testVal  = eval(testAs,testLM,testNF,1.e-10);
    refVal   = 1.27160493839876543211148148148;
    errBound = 1.27160493839876543211148148148 * std::max(eps/(1.-1.e-10),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, plus, x = 10^-10 = 1.e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-9);
    refVal   = 1.27160493954320987781481481609;
    errBound = 1.27160493954320987781481481609 * std::max(eps/(1.-1.e-9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, plus, x = 10^-9 = 1.e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-8);
    refVal   = 1.27160495098765444814814941975;
    errBound = 1.27160495098765444814814941975 * std::max(eps/(1.-1.e-8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, plus, x = 10^-8 = 1.e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-7);
    refVal   = 1.27160506543211148148275308655;
    errBound = 1.27160506543211148148275308655 * std::max(eps/(1.-1.e-7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, plus, x = 10^-7 = 1.e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-6);
    refVal   = 1.27160620987781481608642102469;
    errBound = 1.27160620987781481608642102469 * std::max(eps/(1.-1.e-6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, plus, x = 10^-6 = 1.e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-5);
    refVal   = 1.2716176544481494197658025963;
    errBound = 1.2716176544481494197658025963 * std::max(eps/(1.-1.e-5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, plus, x = 10^-5 = 1.e-5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0001);
    refVal   = 1.2717321114827532135929642347;
    errBound = 1.2717321114827532135929642347 * std::max(eps/(1.-0.0001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, plus, x = 10^-4 = 0.0001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.001);
    refVal   = 1.27287781608769263090250744572;
    errBound = 1.27287781608769263090250744572 * std::max(eps/(1.-0.001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, plus, x = 10^-3 = 0.001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.01);
    refVal   = 1.28444943259758074572889387704;
    errBound = 1.28444943259758074572889387704 * std::max(eps/(1.-0.01),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, plus, x = 10^-2 = 0.01";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.05);
    refVal   = 1.33853151397011046133853151397;
    errBound = 1.33853151397011046133853151397 * std::max(eps/(1.-0.05),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, plus, x = 0.05 = 0.05";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.1);
    refVal   = 1.41289437585733882030178326475;
    errBound = 1.41289437585733882030178326475 * std::max(eps/(1.-0.1),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, plus, x = 0.1 = 0.1";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.15);
    refVal   = 1.49600580973129992737835875091;
    errBound = 1.49600580973129992737835875091 * std::max(eps/(1.-0.15),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, plus, x = 0.15 = 0.15";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.2);
    refVal   = 1.58950617283950617283950617284;
    errBound = 1.58950617283950617283950617284 * std::max(eps/(1.-0.2),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, plus, x = 0.2 = 0.2";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.25);
    refVal   = 1.6954732510288065843621399177;
    errBound = 1.6954732510288065843621399177 * std::max(eps/(1.-0.25),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, plus, x = 0.25 = 0.25";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.3);
    refVal   = 1.81657848324514991181657848325;
    errBound = 1.81657848324514991181657848325 * std::max(eps/(1.-0.3),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, plus, x = 0.3 = 0.3";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.35);
    refVal   = 1.95631528964862298195631528965;
    errBound = 1.95631528964862298195631528965 * std::max(eps/(1.-0.35),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, plus, x = 0.35 = 0.35";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.4);
    refVal   = 2.11934156378600823045267489712;
    errBound = 2.11934156378600823045267489712 * std::max(eps/(1.-0.4),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, plus, x = 0.4 = 0.4";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.45);
    refVal   = 2.31200897867564534231200897868;
    errBound = 2.31200897867564534231200897868 * std::max(eps/(1.-0.45),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, plus, x = 0.45 = 0.45";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.5);
    refVal   = 2.54320987654320987654320987654;
    errBound = 2.54320987654320987654320987654 * std::max(eps/(1.-0.5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, plus, x = 0.5 = 0.5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.55);
    refVal   = 2.82578875171467764060356652949;
    errBound = 2.82578875171467764060356652949 * std::max(eps/(1.-0.55),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, plus, x = 0.55 = 0.55";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.6);
    refVal   = 3.17901234567901234567901234568;
    errBound = 3.17901234567901234567901234568 * std::max(eps/(1.-0.6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, plus, x = 0.6 = 0.6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.65);
    refVal   = 3.63315696649029982363315696649;
    errBound = 3.63315696649029982363315696649 * std::max(eps/(1.-0.65),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, plus, x = 0.65 = 0.65";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.7);
    refVal   = 4.23868312757201646090534979424;
    errBound = 4.23868312757201646090534979424 * std::max(eps/(1.-0.7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, plus, x = 0.7 = 0.7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.75);
    refVal   = 5.08641975308641975308641975309;
    errBound = 5.08641975308641975308641975309 * std::max(eps/(1.-0.75),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, plus, x = 0.75 = 0.75";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.8);
    refVal   = 6.35802469135802469135802469136;
    errBound = 6.35802469135802469135802469136 * std::max(eps/(1.-0.8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, plus, x = 0.8 = 0.8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.85);
    refVal   = 8.47736625514403292181069958848;
    errBound = 8.47736625514403292181069958848 * std::max(eps/(1.-0.85),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, plus, x = 0.85 = 0.85";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9);
    refVal   = 12.7160493827160493827160493827;
    errBound = 12.7160493827160493827160493827 * std::max(eps/(1.-0.9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, plus, x = 0.9 = 0.9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.95);
    refVal   = 25.4320987654320987654320987654;
    errBound = 25.4320987654320987654320987654 * std::max(eps/(1.-0.95),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, plus, x = 0.95 = 0.95";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99);
    refVal   = 127.160493827160493827160493827;
    errBound = 127.160493827160493827160493827 * std::max(eps/(1.-0.99),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, plus, x = 1-10^-2 = 0.99";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999);
    refVal   = 1271.60493827160493827160493827;
    errBound = 1271.60493827160493827160493827 * std::max(eps/(1.-0.999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, plus, x = 1-10^-3 = 0.999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999);
    refVal   = 12716.0493827160493827160493827;
    errBound = 12716.0493827160493827160493827 * std::max(eps/(1.-0.9999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, plus, x = 1-10^-4 = 0.9999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999);
    refVal   = 127160.493827160493827160493827;
    errBound = 127160.493827160493827160493827 * std::max(eps/(1.-0.99999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, plus, x = 1-10^-5 = 0.99999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999);
    refVal   = 1.27160493827160493827160493827e6;
    errBound = 1.27160493827160493827160493827e6 * std::max(eps/(1.-0.999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, plus, x = 1-10^-6 = 0.999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999);
    refVal   = 1.27160493827160493827160493827e7;
    errBound = 1.27160493827160493827160493827e7 * std::max(eps/(1.-0.9999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, plus, x = 1-10^-7 = 0.9999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999999);
    refVal   = 1.27160493827160493827160493827e8;
    errBound = 1.27160493827160493827160493827e8 * std::max(eps/(1.-0.99999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, plus, x = 1-10^-8 = 0.99999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999);
    refVal   = 1.27160493827160493827160493827e9;
    errBound = 1.27160493827160493827160493827e9 * std::max(eps/(1.-0.999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, plus, x = 1-10^-9 = 0.999999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999999);
    refVal   = 1.27160493827160493827160493827e10;
    errBound = 1.27160493827160493827160493827e10 * std::max(eps/(1.-0.9999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSOdd, TruncAs2NormalDoubles, plus, x = 1-10^-10 = 0.9999999999";
  }
}



TEST(AqqQNSOddXspace,TruncAs2NormalDoublesDeltaEvaluation)
{
  auto eval = AqqQNSOdd_delta.truncate(2);
  EXPECT_DOUBLE_EQ(eval(testAs,testLM,testNF),1.45677694886353911970334838331)
    << "^-- AqqQNSOdd, TruncAs2NormalDoubles, delta";
}



TEST(AqqQNSOddNspace,FullMoments)
{
  double refVal = 0., errBound = 0.;
  auto& rpd = AqqQNSOdd;
  integration_engine_gsl engine;
  auto mom = make_mellin_moment(rpd, engine, testAs, testLM, testNF);
  std::tuple<integration_status, double, double> res;

  {
    SCOPED_TRACE("AqqQNSOdd Full Mellin moment N=3");
    refVal = -23.3078491429586047800120406239;
    errBound = 23.3078491429586047800120406239 * 3.e-11;
    res = mom.integrate(3, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AqqQNSOdd Full Mellin moment N=5");
    refVal = -32.9762130808769540631005999065;
    errBound = 32.9762130808769540631005999065 * 3.e-11;
    res = mom.integrate(5, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AqqQNSOdd Full Mellin moment N=7");
    refVal = -39.3417230373840618703749129601;
    errBound = 39.3417230373840618703749129601 * 3.e-11;
    res = mom.integrate(7, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AqqQNSOdd Full Mellin moment N=9");
    refVal = -44.1192067454349801977051537392;
    errBound = 44.1192067454349801977051537392 * 3.e-11;
    res = mom.integrate(9, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

}


TEST(AqqQNSOddNspace,TruncAs0Moments)
{
  double refVal = 0., errBound = 0.;
  auto rpd = AqqQNSOdd.truncate(0);
  integration_engine_gsl engine;
  auto mom = make_mellin_moment(rpd, engine, testAs, testLM, testNF);
  std::tuple<integration_status, double, double> res;

  {
    SCOPED_TRACE("AqqQNSOdd TruncAs0 Mellin moment N=3");
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
    SCOPED_TRACE("AqqQNSOdd TruncAs0 Mellin moment N=5");
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
    SCOPED_TRACE("AqqQNSOdd TruncAs0 Mellin moment N=7");
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

  {
    SCOPED_TRACE("AqqQNSOdd TruncAs0 Mellin moment N=9");
    refVal = 1.;
    errBound = 1. * 3.e-11;
    res = mom.integrate(9, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

}


TEST(AqqQNSOddNspace,TruncAs1Moments)
{
  double refVal = 0., errBound = 0.;
  auto rpd = AqqQNSOdd.truncate(1);
  integration_engine_gsl engine;
  auto mom = make_mellin_moment(rpd, engine, testAs, testLM, testNF);
  std::tuple<integration_status, double, double> res;

  {
    SCOPED_TRACE("AqqQNSOdd TruncAs1 Mellin moment N=3");
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
    SCOPED_TRACE("AqqQNSOdd TruncAs1 Mellin moment N=5");
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
    SCOPED_TRACE("AqqQNSOdd TruncAs1 Mellin moment N=7");
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

  {
    SCOPED_TRACE("AqqQNSOdd TruncAs1 Mellin moment N=9");
    refVal = 1.;
    errBound = 1. * 3.e-11;
    res = mom.integrate(9, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

}


TEST(AqqQNSOddNspace,TruncAs2Moments)
{
  double refVal = 0., errBound = 0.;
  auto rpd = AqqQNSOdd.truncate(2);
  integration_engine_gsl engine;
  auto mom = make_mellin_moment(rpd, engine, testAs, testLM, testNF);
  std::tuple<integration_status, double, double> res;

  {
    SCOPED_TRACE("AqqQNSOdd TruncAs2 Mellin moment N=3");
    refVal = -0.608506944444444444444444444444;
    errBound = 0.608506944444444444444444444444 * 3.e-11;
    res = mom.integrate(3, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AqqQNSOdd TruncAs2 Mellin moment N=5");
    refVal = -1.28033899176954732510288065844;
    errBound = 1.28033899176954732510288065844 * 3.e-11;
    res = mom.integrate(5, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AqqQNSOdd TruncAs2 Mellin moment N=7");
    refVal = -1.71850696028521637932068771071;
    errBound = 1.71850696028521637932068771071 * 3.e-11;
    res = mom.integrate(7, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AqqQNSOdd TruncAs2 Mellin moment N=9");
    refVal = -2.04432326280298192828637105113;
    errBound = 2.04432326280298192828637105113 * 3.e-11;
    res = mom.integrate(9, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

}
