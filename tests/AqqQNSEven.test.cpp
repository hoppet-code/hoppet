/******************************************************************************
 * This file is part of libome                                                *
 * Copyright (C) 2025 Arnd Behring, Kay Schoenwald                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 ******************************************************************************/

#include <limits>
#include <gtest/gtest.h>
#include <ome/integration_engine_gsl.h>
#include <ome/mellin.h>
#include <ome/AqqQNSEven.h>

using namespace ome;

const double eps = std::numeric_limits<double>::epsilon();
const double testAs = 0.25;
const double testLM = -5.;
const double testNF = 3.;



TEST(AqqQNSEvenXspace,FullPowersOf4RegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto& eval = AqqQNSEven_reg;

  {
    testVal  = eval(testAs,testLM,testNF,9.31322574615478515625e-10);
    refVal   = 2163.01257897213413594095782306;
    errBound = 2163.01257897213413594095782306 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, reg, x = 4^-15 = 9.31322574615478515625e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.7252902984619140625e-9);
    refVal   = 1694.3241182457240241235335396;
    errBound = 1694.3241182457240241235335396 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, reg, x = 4^-14 = 3.7252902984619140625e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.490116119384765625e-8);
    refVal   = 1306.69699693582389247892965236;
    errBound = 1306.69699693582389247892965236 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, reg, x = 4^-13 = 1.490116119384765625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,5.9604644775390625e-8);
    refVal   = 990.164864636953800181033774287;
    errBound = 990.164864636953800181033774287 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, reg, x = 4^-12 = 5.9604644775390625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,2.384185791015625e-7);
    refVal   = 735.399802763619240187111269966;
    errBound = 735.399802763619240187111269966 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, reg, x = 4^-11 = 2.384185791015625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,9.5367431640625e-7);
    refVal   = 533.712457515297956113756451468;
    errBound = 533.712457515297956113756451468 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, reg, x = 4^-10 = 9.5367431640625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.814697265625e-6);
    refVal   = 377.052392861420274601188957891;
    errBound = 377.052392861420274601188957891 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, reg, x = 4^-9 = 3.814697265625e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0000152587890625);
    refVal   = 258.00897498880900526229725398;
    errBound = 258.00897498880900526229725398 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, reg, x = 4^-8 = 0.0000152587890625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00006103515625);
    refVal   = 169.813416919913244925714868905;
    errBound = 169.813416919913244925714868905 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, reg, x = 4^-7 = 0.00006103515625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.000244140625);
    refVal   = 106.342960250369724136402795407;
    errBound = 106.342960250369724136402795407 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, reg, x = 4^-6 = 0.000244140625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0009765625);
    refVal   = 62.1277576737131427761892565522;
    errBound = 62.1277576737131427761892565522 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, reg, x = 4^-5 = 0.0009765625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00390625);
    refVal   = 32.3573609026138544277952086662;
    errBound = 32.3573609026138544277952086662 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, reg, x = 4^-4 = 0.00390625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.015625);
    refVal   = 12.8737856842499793216701865;
    errBound = 12.8737856842499793216701865 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, reg, x = 4^-3 = 0.015625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0625);
    refVal   = 0.137089640189828775255456845247;
    errBound = 0.137089640189828775255456845247 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, reg, x = 4^-2 = 0.0625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9375);
    refVal   = -20.5606144911050942416354711787;
    errBound = 20.5606144911050942416354711787 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, reg, x = 1-4^-2 = 0.9375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.984375);
    refVal   = -25.7393913115021134405354237444;
    errBound = 25.7393913115021134405354237444 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, reg, x = 1-4^-3 = 0.984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99609375);
    refVal   = -31.4788661492384389522093125103;
    errBound = 31.4788661492384389522093125103 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, reg, x = 1-4^-4 = 0.99609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9990234375);
    refVal   = -38.0052048221206458849255701642;
    errBound = 38.0052048221206458849255701642 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, reg, x = 1-4^-5 = 0.9990234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999755859375);
    refVal   = -45.5735908614332167370004037345;
    errBound = 45.5735908614332167370004037345 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, reg, x = 1-4^-6 = 0.999755859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99993896484375);
    refVal   = -54.4466989670536096670982061209;
    errBound = 54.4466989670536096670982061209 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, reg, x = 1-4^-7 = 0.99993896484375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999847412109375);
    refVal   = -64.888203024122172429529717644;
    errBound = 64.888203024122172429529717644 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, reg, x = 1-4^-8 = 0.9999847412109375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999996185302734375);
    refVal   = -77.1616117562122498393108523892;
    errBound = 77.1616117562122498393108523892 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, reg, x = 1-4^-9 = 0.999996185302734375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999904632568359375);
    refVal   = -91.5302335523008499063243382;
    errBound = 91.5302335523008499063243382 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, reg, x = 1-4^-10 = 0.99999904632568359375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999997615814208984375);
    refVal   = -108.2572707802312425638253576743;
    errBound = 108.2572707802312425638253576743 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, reg, x = 1-4^-11 = 0.9999997615814208984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999940395355224609375);
    refVal   = -127.6058806458398988678312222311;
    errBound = 127.6058806458398988678312222311 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, reg, x = 1-4^-12 = 0.999999940395355224609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999998509883880615234375);
    refVal   = -149.839203073412364044948910898;
    errBound = 149.839203073412364044948910898 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, reg, x = 1-4^-13 = 0.99999998509883880615234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999962747097015380859375);
    refVal   = -175.220371797634271109378838964;
    errBound = 175.220371797634271109378838964 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, reg, x = 1-4^-14 = 0.9999999962747097015380859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999068677425384521484375);
    refVal   = -204.01251843569601974891835882;
    errBound = 204.01251843569601974891835882 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, reg, x = 1-4^-15 = 0.999999999068677425384521484375";
  }
}

TEST(AqqQNSEvenXspace,FullPowersOf4PlusEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;
  auto& eval = AqqQNSEven_plus;
  {
    testVal  = eval(testAs,testLM,testNF,9.31322574615478515625e-10);
    refVal   = 19.8478639549625839237424922507;
    errBound = 19.8478639549625839237424922507 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, plus, x = 4^-15 = 9.31322574615478515625e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.7252902984619140625e-9);
    refVal   = 19.8478640104168754077863490609;
    errBound = 19.8478640104168754077863490609 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, plus, x = 4^-14 = 3.7252902984619140625e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.490116119384765625e-8);
    refVal   = 19.8478642322340444427118335806;
    errBound = 19.8478642322340444427118335806 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, plus, x = 4^-13 = 1.490116119384765625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,5.9604644775390625e-8);
    refVal   = 19.8478651195027701624175971435;
    errBound = 19.8478651195027701624175971435 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, plus, x = 4^-12 = 5.9604644775390625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,2.384185791015625e-7);
    refVal   = 19.8478686685784663214880366929;
    errBound = 19.8478686685784663214880366929 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, plus, x = 4^-11 = 2.384185791015625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,9.5367431640625e-7);
    refVal   = 19.847882864893943453643334231;
    errBound = 19.847882864893943453643334231 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, plus, x = 4^-10 = 9.5367431640625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.814697265625e-6);
    refVal   = 19.8479396503589326788280243906;
    errBound = 19.8479396503589326788280243906 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, plus, x = 4^-9 = 3.814697265625e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0000152587890625);
    refVal   = 19.8481667954682295308843167689;
    errBound = 19.8481667954682295308843167689 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, plus, x = 4^-8 = 0.0000152587890625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00006103515625);
    refVal   = 19.8490754278979799534045000832;
    errBound = 19.8490754278979799534045000832 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, plus, x = 4^-7 = 0.00006103515625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.000244140625);
    refVal   = 19.8527107896979612684142811272;
    errBound = 19.8527107896979612684142811272 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, plus, x = 4^-6 = 0.000244140625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0009765625);
    refVal   = 19.8672655630041914452972827996;
    errBound = 19.8672655630041914452972827996 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, plus, x = 4^-5 = 0.0009765625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00390625);
    refVal   = 19.9256986970130273024893336313;
    errBound = 19.9256986970130273024893336313 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, plus, x = 4^-4 = 0.00390625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.015625);
    refVal   = 20.1629093957869919132332542698;
    errBound = 20.1629093957869919132332542698 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, plus, x = 4^-3 = 0.015625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0625);
    refVal   = 21.1710548655763415088949169833;
    errBound = 21.1710548655763415088949169833 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, plus, x = 4^-2 = 0.0625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9375);
    refVal   = 317.565822983645122633423754749;
    errBound = 317.565822983645122633423754749 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, plus, x = 1-4^-2 = 0.9375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.984375);
    refVal   = 1270.263291934580490533695019;
    errBound = 1270.263291934580490533695019 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, plus, x = 1-4^-3 = 0.984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99609375);
    refVal   = 5081.05316773832196213478007599;
    errBound = 5081.05316773832196213478007599 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, plus, x = 1-4^-4 = 0.99609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9990234375);
    refVal   = 20324.212670953287848539120304;
    errBound = 20324.212670953287848539120304 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, plus, x = 1-4^-5 = 0.9990234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999755859375);
    refVal   = 81296.8506838131513941564812158;
    errBound = 81296.8506838131513941564812158 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, plus, x = 1-4^-6 = 0.999755859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99993896484375);
    refVal   = 325187.402735252605576625924863;
    errBound = 325187.402735252605576625924863 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, plus, x = 1-4^-7 = 0.99993896484375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999847412109375);
    refVal   = 1.30074961094101042230650369945e6;
    errBound = 1.30074961094101042230650369945e6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, plus, x = 1-4^-8 = 0.9999847412109375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999996185302734375);
    refVal   = 5.20299844376404168922601479781e6;
    errBound = 5.20299844376404168922601479781e6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, plus, x = 1-4^-9 = 0.999996185302734375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999904632568359375);
    refVal   = 2.08119937750561667569040591912e7;
    errBound = 2.08119937750561667569040591912e7 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, plus, x = 1-4^-10 = 0.99999904632568359375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999997615814208984375);
    refVal   = 8.3247975100224667027616236765e7;
    errBound = 8.3247975100224667027616236765e7 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, plus, x = 1-4^-11 = 0.9999997615814208984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999940395355224609375);
    refVal   = 3.3299190040089866811046494706e8;
    errBound = 3.3299190040089866811046494706e8 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, plus, x = 1-4^-12 = 0.999999940395355224609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999998509883880615234375);
    refVal   = 1.33196760160359467244185978824e9;
    errBound = 1.33196760160359467244185978824e9 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, plus, x = 1-4^-13 = 0.99999998509883880615234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999962747097015380859375);
    refVal   = 5.32787040641437868976743915296e9;
    errBound = 5.32787040641437868976743915296e9 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, plus, x = 1-4^-14 = 0.9999999962747097015380859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999068677425384521484375);
    refVal   = 2.13114816256575147590697566118e10;
    errBound = 2.13114816256575147590697566118e10 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullPowersOf4, plus, x = 1-4^-15 = 0.999999999068677425384521484375";
  }
}



TEST(AqqQNSEvenXspace,FullNormalDoublesRegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto& eval = AqqQNSEven_reg;

  {
    testVal  = eval(testAs,testLM,testNF,1.e-10);
    refVal   = 3114.73001916455105898966354378;
    errBound = 3114.73001916455105898966354378 * std::max(eps/(1.-1.e-10),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, reg, x = 10^-10 = 1.e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-9);
    refVal   = 2136.81935296373636488670387411;
    errBound = 2136.81935296373636488670387411 * std::max(eps/(1.-1.e-9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, reg, x = 10^-9 = 1.e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-8);
    refVal   = 1410.510902299631024359251334;
    errBound = 1410.510902299631024359251334 * std::max(eps/(1.-1.e-8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, reg, x = 10^-8 = 1.e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-7);
    refVal   = 888.327796272382111359670238662;
    errBound = 888.327796272382111359670238662 * std::max(eps/(1.-1.e-7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, reg, x = 10^-7 = 1.e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-6);
    refVal   = 527.652261077826462907439640711;
    errBound = 527.652261077826462907439640711 * std::max(eps/(1.-1.e-6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, reg, x = 10^-6 = 1.e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-5);
    refVal   = 290.728054567630464576456663041;
    errBound = 290.728054567630464576456663041 * std::max(eps/(1.-1.e-5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, reg, x = 10^-5 = 1.e-5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0001);
    refVal   = 144.671489555516015724766514984;
    errBound = 144.671489555516015724766514984 * std::max(eps/(1.-0.0001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, reg, x = 10^-4 = 0.0001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.001);
    refVal   = 61.5097668290227279613307664954;
    errBound = 61.5097668290227279613307664954 * std::max(eps/(1.-0.001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, reg, x = 10^-3 = 0.001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.01);
    refVal   = 18.251447613607098523063637119;
    errBound = 18.251447613607098523063637119 * std::max(eps/(1.-0.01),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, reg, x = 10^-2 = 0.01";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.05);
    refVal   = 1.86058928218824556446284086578;
    errBound = 1.86058928218824556446284086578 * std::max(eps/(1.-0.05),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, reg, x = 0.05 = 0.05";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.1);
    refVal   = -3.18596704259147533517278785729;
    errBound = 3.18596704259147533517278785729 * std::max(eps/(1.-0.1),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, reg, x = 0.1 = 0.1";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.15);
    refVal   = -5.7761070736252844631768481898;
    errBound = 5.7761070736252844631768481898 * std::max(eps/(1.-0.15),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, reg, x = 0.15 = 0.15";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.2);
    refVal   = -7.49051656312764463605575424831;
    errBound = 7.49051656312764463605575424831 * std::max(eps/(1.-0.2),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, reg, x = 0.2 = 0.2";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.25);
    refVal   = -8.76483678585161549027725291887;
    errBound = 8.76483678585161549027725291887 * std::max(eps/(1.-0.25),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, reg, x = 0.25 = 0.25";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.3);
    refVal   = -9.77968997618572581476766717749;
    errBound = 9.77968997618572581476766717749 * std::max(eps/(1.-0.3),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, reg, x = 0.3 = 0.3";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.35);
    refVal   = -10.6277008478908758924459622356;
    errBound = 10.6277008478908758924459622356 * std::max(eps/(1.-0.35),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, reg, x = 0.35 = 0.35";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.4);
    refVal   = -11.3635457895064387324337238184;
    errBound = 11.3635457895064387324337238184 * std::max(eps/(1.-0.4),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, reg, x = 0.4 = 0.4";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.45);
    refVal   = -12.0232690806995261345957109834;
    errBound = 12.0232690806995261345957109834 * std::max(eps/(1.-0.45),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, reg, x = 0.45 = 0.45";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.5);
    refVal   = -12.6331637229158185956132973196;
    errBound = 12.6331637229158185956132973196 * std::max(eps/(1.-0.5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, reg, x = 0.5 = 0.5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.55);
    refVal   = -13.21453841095870756026858325408;
    errBound = 13.21453841095870756026858325408 * std::max(eps/(1.-0.55),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, reg, x = 0.55 = 0.55";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.6);
    refVal   = -13.78679758362243204247098758267;
    errBound = 13.78679758362243204247098758267 * std::max(eps/(1.-0.6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, reg, x = 0.6 = 0.6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.65);
    refVal   = -14.3700270244950286892936920001;
    errBound = 14.3700270244950286892936920001 * std::max(eps/(1.-0.65),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, reg, x = 0.65 = 0.65";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.7);
    refVal   = -14.9879957275802406179687634022;
    errBound = 14.9879957275802406179687634022 * std::max(eps/(1.-0.7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, reg, x = 0.7 = 0.7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.75);
    refVal   = -15.6728577781169907734330470834;
    errBound = 15.6728577781169907734330470834 * std::max(eps/(1.-0.75),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, reg, x = 0.75 = 0.75";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.8);
    refVal   = -16.4744301811104800026969892763;
    errBound = 16.4744301811104800026969892763 * std::max(eps/(1.-0.8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, reg, x = 0.8 = 0.8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.85);
    refVal   = -17.4826322599955600627165407843;
    errBound = 17.4826322599955600627165407843 * std::max(eps/(1.-0.85),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, reg, x = 0.85 = 0.85";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9);
    refVal   = -18.8972597634602817291367015833;
    errBound = 18.8972597634602817291367015833 * std::max(eps/(1.-0.9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, reg, x = 0.9 = 0.9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.95);
    refVal   = -21.3642785871994788022081869879;
    errBound = 21.3642785871994788022081869879 * std::max(eps/(1.-0.95),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, reg, x = 0.95 = 0.95";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99);
    refVal   = -27.5159644171247246331130729548;
    errBound = 27.5159644171247246331130729548 * std::max(eps/(1.-0.99),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, reg, x = 1-10^-2 = 0.99";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999);
    refVal   = -37.8855365821650376186509007741;
    errBound = 37.8855365821650376186509007741 * std::max(eps/(1.-0.999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, reg, x = 1-10^-3 = 0.999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999);
    refVal   = -51.1204405899774297217903320948;
    errBound = 51.1204405899774297217903320948 * std::max(eps/(1.-0.9999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, reg, x = 1-10^-4 = 0.9999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999);
    refVal   = -68.4231546064017644280394632081;
    errBound = 68.4231546064017644280394632081 * std::max(eps/(1.-0.99999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, reg, x = 1-10^-5 = 0.99999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999);
    refVal   = -91.001132973010669797972833609;
    errBound = 91.001132973010669797972833609 * std::max(eps/(1.-0.999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, reg, x = 1-10^-6 = 0.999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999);
    refVal   = -120.0606511368158634492953481972;
    errBound = 120.0606511368158634492953481972 * std::max(eps/(1.-0.9999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, reg, x = 1-10^-7 = 0.9999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999999);
    refVal   = -156.807555957234490879962749595;
    errBound = 156.807555957234490879962749595 * std::max(eps/(1.-0.99999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, reg, x = 1-10^-8 = 0.99999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999);
    refVal   = -202.447602274042630258119964549;
    errBound = 202.447602274042630258119964549 * std::max(eps/(1.-0.999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, reg, x = 1-10^-9 = 0.999999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999999);
    refVal   = -258.186528604675458297909307061;
    errBound = 258.186528604675458297909307061 * std::max(eps/(1.-0.9999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, reg, x = 1-10^-10 = 0.9999999999";
  }
}

TEST(AqqQNSEvenXspace,FullNormalDoublesPlusEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;
  auto& eval = AqqQNSEven_plus;
  {
    testVal  = eval(testAs,testLM,testNF,1.e-10);
    refVal   = 19.8478639384626065584352453277;
    errBound = 19.8478639384626065584352453277 * std::max(eps/(1.-1.e-10),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, plus, x = 10^-10 = 1.e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-9);
    refVal   = 19.8478639563256841209146687927;
    errBound = 19.8478639563256841209146687927 * std::max(eps/(1.-1.e-9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, plus, x = 10^-9 = 1.e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-8);
    refVal   = 19.8478641349564615141535998134;
    errBound = 19.8478641349564615141535998134 * std::max(eps/(1.-1.e-8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, plus, x = 10^-8 = 1.e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-7);
    refVal   = 19.8478659212644122910302137749;
    errBound = 19.8478659212644122910302137749 * std::max(eps/(1.-1.e-7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, plus, x = 10^-7 = 1.e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-6);
    refVal   = 19.8478837843616045261935108653;
    errBound = 19.8478837843616045261935108653 * std::max(eps/(1.-1.e-6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, plus, x = 10^-6 = 1.e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-5);
    refVal   = 19.8480624171019911845008296801;
    errBound = 19.8480624171019911845008296801 * std::max(eps/(1.-1.e-5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, plus, x = 10^-5 = 1.e-5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0001);
    refVal   = 19.8498489213699571603050151733;
    errBound = 19.8498489213699571603050151733 * std::max(eps/(1.-0.0001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, plus, x = 10^-4 = 0.0001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.001);
    refVal   = 19.8677316681459661307197043762;
    errBound = 19.8677316681459661307197043762 * std::max(eps/(1.-0.001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, plus, x = 10^-3 = 0.001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.01);
    refVal   = 20.0483474105836567319080653251;
    errBound = 20.0483474105836567319080653251 * std::max(eps/(1.-0.01),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, plus, x = 10^-2 = 0.01";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.05);
    refVal   = 20.8924883541871791206199838651;
    errBound = 20.8924883541871791206199838651 * std::max(eps/(1.-0.05),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, plus, x = 0.05 = 0.05";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.1);
    refVal   = 22.0531821516420224050988718576;
    errBound = 22.0531821516420224050988718576 * std::max(eps/(1.-0.1),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, plus, x = 0.1 = 0.1";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.15);
    refVal   = 23.3504281605621413701046878492;
    errBound = 23.3504281605621413701046878492 * std::max(eps/(1.-0.15),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, plus, x = 0.15 = 0.15";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.2);
    refVal   = 24.8098299205972752057362308398;
    errBound = 24.8098299205972752057362308398 * std::max(eps/(1.-0.2),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, plus, x = 0.2 = 0.2";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.25);
    refVal   = 26.4638185819704268861186462291;
    errBound = 26.4638185819704268861186462291 * std::max(eps/(1.-0.25),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, plus, x = 0.25 = 0.25";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.3);
    refVal   = 28.3540913378254573779842638169;
    errBound = 28.3540913378254573779842638169 * std::max(eps/(1.-0.3),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, plus, x = 0.3 = 0.3";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.35);
    refVal   = 30.5351752868889540993676687259;
    errBound = 30.5351752868889540993676687259 * std::max(eps/(1.-0.35),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, plus, x = 0.35 = 0.35";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.4);
    refVal   = 33.0797732274630336076483077864;
    errBound = 33.0797732274630336076483077864 * std::max(eps/(1.-0.4),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, plus, x = 0.4 = 0.4";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.45);
    refVal   = 36.0870253390505821174345175851;
    errBound = 36.0870253390505821174345175851 * std::max(eps/(1.-0.45),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, plus, x = 0.45 = 0.45";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.5);
    refVal   = 39.6957278729556403291779693437;
    errBound = 39.6957278729556403291779693437 * std::max(eps/(1.-0.5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, plus, x = 0.5 = 0.5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.55);
    refVal   = 44.1063643032840448101977437152;
    errBound = 44.1063643032840448101977437152 * std::max(eps/(1.-0.55),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, plus, x = 0.55 = 0.55";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.6);
    refVal   = 49.6196598411945504114724616796;
    errBound = 49.6196598411945504114724616796 * std::max(eps/(1.-0.6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, plus, x = 0.6 = 0.6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.65);
    refVal   = 56.7081826756509147559685276338;
    errBound = 56.7081826756509147559685276338 * std::max(eps/(1.-0.65),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, plus, x = 0.65 = 0.65";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.7);
    refVal   = 66.1595464549260672152966155728;
    errBound = 66.1595464549260672152966155728 * std::max(eps/(1.-0.7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, plus, x = 0.7 = 0.7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.75);
    refVal   = 79.3914557459112806583559386873;
    errBound = 79.3914557459112806583559386873 * std::max(eps/(1.-0.75),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, plus, x = 0.75 = 0.75";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.8);
    refVal   = 99.2393196823891008229449233592;
    errBound = 99.2393196823891008229449233592 * std::max(eps/(1.-0.8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, plus, x = 0.8 = 0.8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.85);
    refVal   = 132.319092909852134430593231146;
    errBound = 132.319092909852134430593231146 * std::max(eps/(1.-0.85),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, plus, x = 0.85 = 0.85";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9);
    refVal   = 198.478639364778201645889846718;
    errBound = 198.478639364778201645889846718 * std::max(eps/(1.-0.9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, plus, x = 0.9 = 0.9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.95);
    refVal   = 396.957278729556403291779693437;
    errBound = 396.957278729556403291779693437 * std::max(eps/(1.-0.95),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, plus, x = 0.95 = 0.95";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99);
    refVal   = 1984.78639364778201645889846718;
    errBound = 1984.78639364778201645889846718 * std::max(eps/(1.-0.99),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, plus, x = 1-10^-2 = 0.99";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999);
    refVal   = 19847.8639364778201645889846718;
    errBound = 19847.8639364778201645889846718 * std::max(eps/(1.-0.999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, plus, x = 1-10^-3 = 0.999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999);
    refVal   = 198478.639364778201645889846718;
    errBound = 198478.639364778201645889846718 * std::max(eps/(1.-0.9999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, plus, x = 1-10^-4 = 0.9999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999);
    refVal   = 1.98478639364778201645889846718e6;
    errBound = 1.98478639364778201645889846718e6 * std::max(eps/(1.-0.99999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, plus, x = 1-10^-5 = 0.99999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999);
    refVal   = 1.98478639364778201645889846718e7;
    errBound = 1.98478639364778201645889846718e7 * std::max(eps/(1.-0.999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, plus, x = 1-10^-6 = 0.999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999);
    refVal   = 1.98478639364778201645889846718e8;
    errBound = 1.98478639364778201645889846718e8 * std::max(eps/(1.-0.9999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, plus, x = 1-10^-7 = 0.9999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999999);
    refVal   = 1.98478639364778201645889846718e9;
    errBound = 1.98478639364778201645889846718e9 * std::max(eps/(1.-0.99999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, plus, x = 1-10^-8 = 0.99999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999);
    refVal   = 1.98478639364778201645889846718e10;
    errBound = 1.98478639364778201645889846718e10 * std::max(eps/(1.-0.999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, plus, x = 1-10^-9 = 0.999999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999999);
    refVal   = 1.98478639364778201645889846718e11;
    errBound = 1.98478639364778201645889846718e11 * std::max(eps/(1.-0.9999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, FullNormalDoubles, plus, x = 1-10^-10 = 0.9999999999";
  }
}



TEST(AqqQNSEvenXspace,FullDeltaEvaluation)
{
  auto& eval = AqqQNSEven_delta;
  EXPECT_DOUBLE_EQ(eval(testAs,testLM,testNF),12.1501748418929015300593526677)
    << "^-- AqqQNSEven, Full, delta";
}







TEST(AqqQNSEvenXspace,TruncAs0NormalDoublesDeltaEvaluation)
{
  auto eval = AqqQNSEven_delta.truncate(0);
  EXPECT_DOUBLE_EQ(eval(testAs,testLM,testNF),1.)
    << "^-- AqqQNSEven, TruncAs0NormalDoubles, delta";
}







TEST(AqqQNSEvenXspace,TruncAs1NormalDoublesDeltaEvaluation)
{
  auto eval = AqqQNSEven_delta.truncate(1);
  EXPECT_DOUBLE_EQ(eval(testAs,testLM,testNF),1.)
    << "^-- AqqQNSEven, TruncAs1NormalDoubles, delta";
}



TEST(AqqQNSEvenXspace,TruncAs2PowersOf4RegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto eval = AqqQNSEven_reg.truncate(2);

  {
    testVal  = eval(testAs,testLM,testNF,9.31322574615478515625e-10);
    refVal   = 17.8217061035486357547375500369;
    errBound = 17.8217061035486357547375500369 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, reg, x = 4^-15 = 9.31322574615478515625e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.7252902984619140625e-9);
    refVal   = 15.7858095460317299789869394768;
    errBound = 15.7858095460317299789869394768 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, reg, x = 4^-14 = 3.7252902984619140625e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.490116119384765625e-8);
    refVal   = 13.856680476843573151637347734;
    errBound = 13.856680476843573151637347734 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, reg, x = 4^-13 = 1.490116119384765625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,5.9604644775390625e-8);
    refVal   = 12.0343192826499716016331048158;
    errBound = 12.0343192826499716016331048158 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, reg, x = 4^-12 = 5.9604644775390625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,2.384185791015625e-7);
    refVal   = 10.3187273173085245697509487745;
    errBound = 10.3187273173085245697509487745 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, reg, x = 4^-11 = 2.384185791015625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,9.5367431640625e-7);
    refVal   = 8.70990926799161824778922135492;
    errBound = 8.70990926799161824778922135492 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, reg, x = 4^-10 = 9.5367431640625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.814697265625e-6);
    refVal   = 7.20788114229739815587154947429;
    errBound = 7.20788114229739815587154947429 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, reg, x = 4^-9 = 3.814697265625e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0000152587890625);
    refVal   = 5.81269669510812129899603649672;
    errBound = 5.81269669510812129899603649672 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, reg, x = 4^-8 = 0.0000152587890625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00006103515625);
    refVal   = 4.52453261341313405656388449925;
    errBound = 4.52453261341313405656388449925 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, reg, x = 4^-7 = 0.00006103515625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.000244140625);
    refVal   = 3.34395356834756400818569948679;
    errBound = 3.34395356834756400818569948679 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, reg, x = 4^-6 = 0.000244140625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0009765625);
    refVal   = 2.27269703009075740014577399795;
    errBound = 2.27269703009075740014577399795 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, reg, x = 4^-5 = 0.0009765625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00390625);
    refVal   = 1.31583931756388828143992875279;
    errBound = 1.31583931756388828143992875279 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, reg, x = 4^-4 = 0.00390625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.015625);
    refVal   = 0.48726127075627485257793501304;
    errBound = 0.48726127075627485257793501304 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, reg, x = 4^-3 = 0.015625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0625);
    refVal   = -0.177135651603830542552739162964;
    errBound = 0.177135651603830542552739162964 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, reg, x = 4^-2 = 0.0625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9375);
    refVal   = -0.38499958453489546857703568513;
    errBound = 0.38499958453489546857703568513 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, reg, x = 1-4^-2 = 0.9375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.984375);
    refVal   = -0.3556611804860762279683890258;
    errBound = 0.3556611804860762279683890258 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, reg, x = 1-4^-3 = 0.984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99609375);
    refVal   = -0.34818370702466657859671556732;
    errBound = 0.34818370702466657859671556732 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, reg, x = 1-4^-4 = 0.99609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9990234375);
    refVal   = -0.346305752864286253367270487748;
    errBound = 0.346305752864286253367270487748 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, reg, x = 1-4^-5 = 0.9990234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999755859375);
    refVal   = -0.345835732823234028838374314547;
    errBound = 0.345835732823234028838374314547 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, reg, x = 1-4^-6 = 0.999755859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99993896484375);
    refVal   = -0.345718194673063083996800636007;
    errBound = 0.345718194673063083996800636007 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, reg, x = 1-4^-7 = 0.99993896484375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999847412109375);
    refVal   = -0.345688808065505242198483338781;
    errBound = 0.345688808065505242198483338781 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, reg, x = 1-4^-8 = 0.9999847412109375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999996185302734375);
    refVal   = -0.345681461284259028495949349392;
    errBound = 0.345681461284259028495949349392 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, reg, x = 1-4^-9 = 0.999996185302734375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999904632568359375);
    refVal   = -0.345679624580862977792671971337;
    errBound = 0.345679624580862977792671971337 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, reg, x = 1-4^-10 = 0.99999904632568359375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999997615814208984375);
    refVal   = -0.345679165404508688721165380759;
    errBound = 0.345679165404508688721165380759 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, reg, x = 1-4^-11 = 0.9999997615814208984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999940395355224609375);
    refVal   = -0.345679050610388536751747507152;
    errBound = 0.345679050610388536751747507152 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, reg, x = 1-4^-12 = 0.999999940395355224609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999998509883880615234375);
    refVal   = -0.345679021911856525029190290443;
    errBound = 0.345679021911856525029190290443 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, reg, x = 1-4^-13 = 0.99999998509883880615234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999962747097015380859375);
    refVal   = -0.345679014737223398740431182894;
    errBound = 0.345679014737223398740431182894 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, reg, x = 1-4^-14 = 0.9999999962747097015380859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999068677425384521484375);
    refVal   = -0.34567901294356510945835919749;
    errBound = 0.34567901294356510945835919749 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, reg, x = 1-4^-15 = 0.999999999068677425384521484375";
  }
}

TEST(AqqQNSEvenXspace,TruncAs2PowersOf4PlusEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;
  auto eval = AqqQNSEven_plus.truncate(2);
  {
    testVal  = eval(testAs,testLM,testNF,9.31322574615478515625e-10);
    refVal   = 1.27160493945587932437941415068;
    errBound = 1.27160493945587932437941415068 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, plus, x = 4^-15 = 9.31322574615478515625e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.7252902984619140625e-9);
    refVal   = 1.27160494300870249593813948108;
    errBound = 1.27160494300870249593813948108 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, plus, x = 4^-14 = 3.7252902984619140625e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.490116119384765625e-8);
    refVal   = 1.27160495721999538070250915835;
    errBound = 1.27160495721999538070250915835 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, plus, x = 4^-13 = 1.490116119384765625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,5.9604644775390625e-8);
    refVal   = 1.27160501406517009623166793273;
    errBound = 1.27160501406517009623166793273 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, plus, x = 4^-12 = 5.9604644775390625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,2.384185791015625e-7);
    refVal   = 1.27160524144591978190711202291;
    errBound = 1.27160524144591978190711202291 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, plus, x = 4^-11 = 2.384185791015625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,9.5367431640625e-7);
    refVal   = 1.2716061509697317023132216768;
    errBound = 1.2716061509697317023132216768 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, plus, x = 4^-10 = 9.5367431640625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.814697265625e-6);
    refVal   = 1.27160978907799027606410091034;
    errBound = 1.27160978907799027606410091034 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, plus, x = 4^-9 = 3.814697265625e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0000152587890625);
    refVal   = 1.27162434171920197199310141504;
    errBound = 1.27162434171920197199310141504 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, plus, x = 4^-8 = 0.0000152587890625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00006103515625);
    refVal   = 1.27168255561508730444008883041;
    errBound = 1.27168255561508730444008883041 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, plus, x = 4^-7 = 0.00006103515625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.000244140625);
    refVal   = 1.27191546450805710064969324229;
    errBound = 1.27191546450805710064969324229 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, plus, x = 4^-6 = 0.000244140625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0009765625);
    refVal   = 1.27284795385153808092876193235;
    errBound = 1.27284795385153808092876193235 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, plus, x = 4^-5 = 0.0009765625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00390625);
    refVal   = 1.27659162430404260469619946744;
    errBound = 1.27659162430404260469619946744 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, plus, x = 4^-4 = 0.00390625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.015625);
    refVal   = 1.2917891436409954928473446992;
    errBound = 1.2917891436409954928473446992 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, plus, x = 4^-3 = 0.015625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0625);
    refVal   = 1.35637860082304526748971193416;
    errBound = 1.35637860082304526748971193416 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, plus, x = 4^-2 = 0.0625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9375);
    refVal   = 20.3456790123456790123456790123;
    errBound = 20.3456790123456790123456790123 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, plus, x = 1-4^-2 = 0.9375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.984375);
    refVal   = 81.3827160493827160493827160494;
    errBound = 81.3827160493827160493827160494 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, plus, x = 1-4^-3 = 0.984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99609375);
    refVal   = 325.530864197530864197530864198;
    errBound = 325.530864197530864197530864198 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, plus, x = 1-4^-4 = 0.99609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9990234375);
    refVal   = 1302.12345679012345679012345679;
    errBound = 1302.12345679012345679012345679 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, plus, x = 1-4^-5 = 0.9990234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999755859375);
    refVal   = 5208.49382716049382716049382716;
    errBound = 5208.49382716049382716049382716 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, plus, x = 1-4^-6 = 0.999755859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99993896484375);
    refVal   = 20833.9753086419753086419753086;
    errBound = 20833.9753086419753086419753086 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, plus, x = 1-4^-7 = 0.99993896484375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999847412109375);
    refVal   = 83335.9012345679012345679012346;
    errBound = 83335.9012345679012345679012346 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, plus, x = 1-4^-8 = 0.9999847412109375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999996185302734375);
    refVal   = 333343.604938271604938271604938;
    errBound = 333343.604938271604938271604938 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, plus, x = 1-4^-9 = 0.999996185302734375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999904632568359375);
    refVal   = 1.33337441975308641975308641975e6;
    errBound = 1.33337441975308641975308641975e6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, plus, x = 1-4^-10 = 0.99999904632568359375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999997615814208984375);
    refVal   = 5.33349767901234567901234567901e6;
    errBound = 5.33349767901234567901234567901e6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, plus, x = 1-4^-11 = 0.9999997615814208984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999940395355224609375);
    refVal   = 2.1333990716049382716049382716e7;
    errBound = 2.1333990716049382716049382716e7 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, plus, x = 1-4^-12 = 0.999999940395355224609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999998509883880615234375);
    refVal   = 8.53359628641975308641975308642e7;
    errBound = 8.53359628641975308641975308642e7 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, plus, x = 1-4^-13 = 0.99999998509883880615234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999962747097015380859375);
    refVal   = 3.41343851456790123456790123457e8;
    errBound = 3.41343851456790123456790123457e8 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, plus, x = 1-4^-14 = 0.9999999962747097015380859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999068677425384521484375);
    refVal   = 1.36537540582716049382716049383e9;
    errBound = 1.36537540582716049382716049383e9 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2PowersOf4, plus, x = 1-4^-15 = 0.999999999068677425384521484375";
  }
}



TEST(AqqQNSEvenXspace,TruncAs2NormalDoublesRegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto eval = AqqQNSEven_reg.truncate(2);

  {
    testVal  = eval(testAs,testLM,testNF,1.e-10);
    refVal   = 21.3230102048050060083039198921;
    errBound = 21.3230102048050060083039198921 * std::max(eps/(1.-1.e-10),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, reg, x = 10^-10 = 1.e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-9);
    refVal   = 17.7146173942785623569748329156;
    errBound = 17.7146173942785623569748329156 * std::max(eps/(1.-1.e-9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, reg, x = 10^-9 = 1.e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-8);
    refVal   = 14.4007746358257017944876439667;
    errBound = 14.4007746358257017944876439667 * std::max(eps/(1.-1.e-8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, reg, x = 10^-8 = 1.e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-7);
    refVal   = 11.3814830621623545539010845103;
    errBound = 11.3814830621623545539010845103 * std::max(eps/(1.-1.e-7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, reg, x = 10^-7 = 1.e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-6);
    refVal   = 8.65675171268660515546140388389;
    errBound = 8.65675171268660515546140388389 * std::max(eps/(1.-1.e-6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, reg, x = 10^-6 = 1.e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-5);
    refVal   = 6.22665025234672232098097860793;
    errBound = 6.22665025234672232098097860793 * std::max(eps/(1.-1.e-5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, reg, x = 10^-5 = 1.e-5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0001);
    refVal   = 4.09168960766253686199028926281;
    errBound = 4.09168960766253686199028926281 * std::max(eps/(1.-0.0001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, reg, x = 10^-4 = 0.0001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.001);
    refVal   = 2.25534283688209607561653854713;
    errBound = 2.25534283688209607561653854713 * std::max(eps/(1.-0.001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, reg, x = 10^-3 = 0.001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.01);
    refVal   = 0.738488095420545443932353103742;
    errBound = 0.738488095420545443932353103742 * std::max(eps/(1.-0.01),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, reg, x = 10^-2 = 0.01";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.05);
    refVal   = -0.0840366145940149601139494039936;
    errBound = 0.0840366145940149601139494039936 * std::max(eps/(1.-0.05),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, reg, x = 0.05 = 0.05";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.1);
    refVal   = -0.351397413919272335738111721405;
    errBound = 0.351397413919272335738111721405 * std::max(eps/(1.-0.1),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, reg, x = 0.1 = 0.1";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.15);
    refVal   = -0.473358476701255463509266244707;
    errBound = 0.473358476701255463509266244707 * std::max(eps/(1.-0.15),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, reg, x = 0.15 = 0.15";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.2);
    refVal   = -0.540138954765073360749215426109;
    errBound = 0.540138954765073360749215426109 * std::max(eps/(1.-0.2),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, reg, x = 0.2 = 0.2";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.25);
    refVal   = -0.578207692291456844284520537334;
    errBound = 0.578207692291456844284520537334 * std::max(eps/(1.-0.25),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, reg, x = 0.25 = 0.25";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.3);
    refVal   = -0.598802628849339143129449818016;
    errBound = 0.598802628849339143129449818016 * std::max(eps/(1.-0.3),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, reg, x = 0.3 = 0.3";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.35);
    refVal   = -0.607691174761393407886918353833;
    errBound = 0.607691174761393407886918353833 * std::max(eps/(1.-0.35),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, reg, x = 0.35 = 0.35";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.4);
    refVal   = -0.608205938556899136474781397016;
    errBound = 0.608205938556899136474781397016 * std::max(eps/(1.-0.4),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, reg, x = 0.4 = 0.4";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.45);
    refVal   = -0.602437372263944781429130689028;
    errBound = 0.602437372263944781429130689028 * std::max(eps/(1.-0.45),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, reg, x = 0.45 = 0.45";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.5);
    refVal   = -0.591778591749814791614315820361;
    errBound = 0.591778591749814791614315820361 * std::max(eps/(1.-0.5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, reg, x = 0.5 = 0.5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.55);
    refVal   = -0.577202330444190885352981519144;
    errBound = 0.577202330444190885352981519144 * std::max(eps/(1.-0.55),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, reg, x = 0.55 = 0.55";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.6);
    refVal   = -0.559413489610987856045682073043;
    errBound = 0.559413489610987856045682073043 * std::max(eps/(1.-0.6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, reg, x = 0.6 = 0.6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.65);
    refVal   = -0.538938605094666549378838062134;
    errBound = 0.538938605094666549378838062134 * std::max(eps/(1.-0.65),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, reg, x = 0.65 = 0.65";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.7);
    refVal   = -0.516181035683198740514234446871;
    errBound = 0.516181035683198740514234446871 * std::max(eps/(1.-0.7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, reg, x = 0.7 = 0.7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.75);
    refVal   = -0.491456453381019672887120740997;
    errBound = 0.491456453381019672887120740997 * std::max(eps/(1.-0.75),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, reg, x = 0.75 = 0.75";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.8);
    refVal   = -0.465016477683309411427140319021;
    errBound = 0.465016477683309411427140319021 * std::max(eps/(1.-0.8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, reg, x = 0.8 = 0.8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.85);
    refVal   = -0.43706489002058434918027274715;
    errBound = 0.43706489002058434918027274715 * std::max(eps/(1.-0.85),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, reg, x = 0.85 = 0.85";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9);
    refVal   = -0.407769047200208012914042224406;
    errBound = 0.407769047200208012914042224406 * std::max(eps/(1.-0.9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, reg, x = 0.9 = 0.9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.95);
    refVal   = -0.377268097411200815373475106491;
    errBound = 0.377268097411200815373475106491 * std::max(eps/(1.-0.95),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, reg, x = 0.95 = 0.95";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99);
    refVal   = -0.352078889482817115652856466876;
    errBound = 0.352078889482817115652856466876 * std::max(eps/(1.-0.99),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, reg, x = 1-10^-2 = 0.99";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999);
    refVal   = -0.346320790001448995456553329598;
    errBound = 0.346320790001448995456553329598 * std::max(eps/(1.-0.999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, reg, x = 1-10^-3 = 0.999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999);
    refVal   = -0.345743207901112644906621845497;
    errBound = 0.345743207901112644906621845497 * std::max(eps/(1.-0.9999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, reg, x = 1-10^-4 = 0.9999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999);
    refVal   = -0.345685432079012223764490732886;
    errBound = 0.345685432079012223764490732886 * std::max(eps/(1.-0.99999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, reg, x = 1-10^-5 = 0.99999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999);
    refVal   = -0.345679654320790123334876449074;
    errBound = 0.345679654320790123334876449074 * std::max(eps/(1.-0.999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, reg, x = 1-10^-6 = 0.999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999);
    refVal   = -0.345679076543207901234445987645;
    errBound = 0.345679076543207901234445987645 * std::max(eps/(1.-0.9999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, reg, x = 1-10^-7 = 0.9999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999999);
    refVal   = -0.345679018765432079012345557099;
    errBound = 0.345679018765432079012345557099 * std::max(eps/(1.-0.99999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, reg, x = 1-10^-8 = 0.99999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999);
    refVal   = -0.345679012987654320790123456668;
    errBound = 0.345679012987654320790123456668 * std::max(eps/(1.-0.999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, reg, x = 1-10^-9 = 0.999999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999999);
    refVal   = -0.345679012409876543207901234568;
    errBound = 0.345679012409876543207901234568 * std::max(eps/(1.-0.9999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, reg, x = 1-10^-10 = 0.9999999999";
  }
}

TEST(AqqQNSEvenXspace,TruncAs2NormalDoublesPlusEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;
  auto eval = AqqQNSEven_plus.truncate(2);
  {
    testVal  = eval(testAs,testLM,testNF,1.e-10);
    refVal   = 1.27160493839876543211148148148;
    errBound = 1.27160493839876543211148148148 * std::max(eps/(1.-1.e-10),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, plus, x = 10^-10 = 1.e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-9);
    refVal   = 1.27160493954320987781481481609;
    errBound = 1.27160493954320987781481481609 * std::max(eps/(1.-1.e-9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, plus, x = 10^-9 = 1.e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-8);
    refVal   = 1.27160495098765444814814941975;
    errBound = 1.27160495098765444814814941975 * std::max(eps/(1.-1.e-8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, plus, x = 10^-8 = 1.e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-7);
    refVal   = 1.27160506543211148148275308655;
    errBound = 1.27160506543211148148275308655 * std::max(eps/(1.-1.e-7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, plus, x = 10^-7 = 1.e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-6);
    refVal   = 1.27160620987781481608642102469;
    errBound = 1.27160620987781481608642102469 * std::max(eps/(1.-1.e-6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, plus, x = 10^-6 = 1.e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-5);
    refVal   = 1.2716176544481494197658025963;
    errBound = 1.2716176544481494197658025963 * std::max(eps/(1.-1.e-5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, plus, x = 10^-5 = 1.e-5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0001);
    refVal   = 1.2717321114827532135929642347;
    errBound = 1.2717321114827532135929642347 * std::max(eps/(1.-0.0001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, plus, x = 10^-4 = 0.0001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.001);
    refVal   = 1.27287781608769263090250744572;
    errBound = 1.27287781608769263090250744572 * std::max(eps/(1.-0.001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, plus, x = 10^-3 = 0.001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.01);
    refVal   = 1.28444943259758074572889387704;
    errBound = 1.28444943259758074572889387704 * std::max(eps/(1.-0.01),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, plus, x = 10^-2 = 0.01";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.05);
    refVal   = 1.33853151397011046133853151397;
    errBound = 1.33853151397011046133853151397 * std::max(eps/(1.-0.05),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, plus, x = 0.05 = 0.05";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.1);
    refVal   = 1.41289437585733882030178326475;
    errBound = 1.41289437585733882030178326475 * std::max(eps/(1.-0.1),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, plus, x = 0.1 = 0.1";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.15);
    refVal   = 1.49600580973129992737835875091;
    errBound = 1.49600580973129992737835875091 * std::max(eps/(1.-0.15),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, plus, x = 0.15 = 0.15";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.2);
    refVal   = 1.58950617283950617283950617284;
    errBound = 1.58950617283950617283950617284 * std::max(eps/(1.-0.2),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, plus, x = 0.2 = 0.2";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.25);
    refVal   = 1.6954732510288065843621399177;
    errBound = 1.6954732510288065843621399177 * std::max(eps/(1.-0.25),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, plus, x = 0.25 = 0.25";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.3);
    refVal   = 1.81657848324514991181657848325;
    errBound = 1.81657848324514991181657848325 * std::max(eps/(1.-0.3),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, plus, x = 0.3 = 0.3";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.35);
    refVal   = 1.95631528964862298195631528965;
    errBound = 1.95631528964862298195631528965 * std::max(eps/(1.-0.35),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, plus, x = 0.35 = 0.35";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.4);
    refVal   = 2.11934156378600823045267489712;
    errBound = 2.11934156378600823045267489712 * std::max(eps/(1.-0.4),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, plus, x = 0.4 = 0.4";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.45);
    refVal   = 2.31200897867564534231200897868;
    errBound = 2.31200897867564534231200897868 * std::max(eps/(1.-0.45),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, plus, x = 0.45 = 0.45";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.5);
    refVal   = 2.54320987654320987654320987654;
    errBound = 2.54320987654320987654320987654 * std::max(eps/(1.-0.5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, plus, x = 0.5 = 0.5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.55);
    refVal   = 2.82578875171467764060356652949;
    errBound = 2.82578875171467764060356652949 * std::max(eps/(1.-0.55),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, plus, x = 0.55 = 0.55";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.6);
    refVal   = 3.17901234567901234567901234568;
    errBound = 3.17901234567901234567901234568 * std::max(eps/(1.-0.6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, plus, x = 0.6 = 0.6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.65);
    refVal   = 3.63315696649029982363315696649;
    errBound = 3.63315696649029982363315696649 * std::max(eps/(1.-0.65),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, plus, x = 0.65 = 0.65";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.7);
    refVal   = 4.23868312757201646090534979424;
    errBound = 4.23868312757201646090534979424 * std::max(eps/(1.-0.7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, plus, x = 0.7 = 0.7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.75);
    refVal   = 5.08641975308641975308641975309;
    errBound = 5.08641975308641975308641975309 * std::max(eps/(1.-0.75),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, plus, x = 0.75 = 0.75";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.8);
    refVal   = 6.35802469135802469135802469136;
    errBound = 6.35802469135802469135802469136 * std::max(eps/(1.-0.8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, plus, x = 0.8 = 0.8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.85);
    refVal   = 8.47736625514403292181069958848;
    errBound = 8.47736625514403292181069958848 * std::max(eps/(1.-0.85),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, plus, x = 0.85 = 0.85";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9);
    refVal   = 12.7160493827160493827160493827;
    errBound = 12.7160493827160493827160493827 * std::max(eps/(1.-0.9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, plus, x = 0.9 = 0.9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.95);
    refVal   = 25.4320987654320987654320987654;
    errBound = 25.4320987654320987654320987654 * std::max(eps/(1.-0.95),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, plus, x = 0.95 = 0.95";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99);
    refVal   = 127.160493827160493827160493827;
    errBound = 127.160493827160493827160493827 * std::max(eps/(1.-0.99),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, plus, x = 1-10^-2 = 0.99";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999);
    refVal   = 1271.60493827160493827160493827;
    errBound = 1271.60493827160493827160493827 * std::max(eps/(1.-0.999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, plus, x = 1-10^-3 = 0.999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999);
    refVal   = 12716.0493827160493827160493827;
    errBound = 12716.0493827160493827160493827 * std::max(eps/(1.-0.9999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, plus, x = 1-10^-4 = 0.9999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999);
    refVal   = 127160.493827160493827160493827;
    errBound = 127160.493827160493827160493827 * std::max(eps/(1.-0.99999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, plus, x = 1-10^-5 = 0.99999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999);
    refVal   = 1.27160493827160493827160493827e6;
    errBound = 1.27160493827160493827160493827e6 * std::max(eps/(1.-0.999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, plus, x = 1-10^-6 = 0.999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999);
    refVal   = 1.27160493827160493827160493827e7;
    errBound = 1.27160493827160493827160493827e7 * std::max(eps/(1.-0.9999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, plus, x = 1-10^-7 = 0.9999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999999);
    refVal   = 1.27160493827160493827160493827e8;
    errBound = 1.27160493827160493827160493827e8 * std::max(eps/(1.-0.99999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, plus, x = 1-10^-8 = 0.99999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999);
    refVal   = 1.27160493827160493827160493827e9;
    errBound = 1.27160493827160493827160493827e9 * std::max(eps/(1.-0.999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, plus, x = 1-10^-9 = 0.999999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999999);
    refVal   = 1.27160493827160493827160493827e10;
    errBound = 1.27160493827160493827160493827e10 * std::max(eps/(1.-0.9999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqqQNSEven, TruncAs2NormalDoubles, plus, x = 1-10^-10 = 0.9999999999";
  }
}



TEST(AqqQNSEvenXspace,TruncAs2NormalDoublesDeltaEvaluation)
{
  auto eval = AqqQNSEven_delta.truncate(2);
  EXPECT_DOUBLE_EQ(eval(testAs,testLM,testNF),1.45677694886353911970334838331)
    << "^-- AqqQNSEven, TruncAs2NormalDoubles, delta";
}



TEST(AqqQNSEvenNspace,FullMoments)
{
  double refVal = 0., errBound = 0.;
  auto& rpd = AqqQNSEven;
  integration_engine_gsl engine;
  auto mom = make_mellin_moment(rpd, engine, testAs, testLM, testNF);
  std::tuple<integration_status, double, double> res;

  {
    SCOPED_TRACE("AqqQNSEven Full Mellin moment N=2");
    refVal = -15.4108510022264770011147721835;
    errBound = 15.4108510022264770011147721835 * 3.e-11;
    res = mom.integrate(2, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AqqQNSEven Full Mellin moment N=4");
    refVal = -28.7584662796242410546737975434;
    errBound = 28.7584662796242410546737975434 * 3.e-11;
    res = mom.integrate(4, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AqqQNSEven Full Mellin moment N=6");
    refVal = -36.4195966475157942743467607082;
    errBound = 36.4195966475157942743467607082 * 3.e-11;
    res = mom.integrate(6, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AqqQNSEven Full Mellin moment N=8");
    refVal = -41.8758791904766359860067192509;
    errBound = 41.8758791904766359860067192509 * 3.e-11;
    res = mom.integrate(8, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AqqQNSEven Full Mellin moment N=10");
    refVal = -46.1282093835538857600877994342;
    errBound = 46.1282093835538857600877994342 * 3.e-11;
    res = mom.integrate(10, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

}


TEST(AqqQNSEvenNspace,TruncAs0Moments)
{
  double refVal = 0., errBound = 0.;
  auto rpd = AqqQNSEven.truncate(0);
  integration_engine_gsl engine;
  auto mom = make_mellin_moment(rpd, engine, testAs, testLM, testNF);
  std::tuple<integration_status, double, double> res;

  {
    SCOPED_TRACE("AqqQNSEven TruncAs0 Mellin moment N=2");
    refVal = 1.;
    errBound = 1. * 3.e-11;
    res = mom.integrate(2, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AqqQNSEven TruncAs0 Mellin moment N=4");
    refVal = 1.;
    errBound = 1. * 3.e-11;
    res = mom.integrate(4, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AqqQNSEven TruncAs0 Mellin moment N=6");
    refVal = 1.;
    errBound = 1. * 3.e-11;
    res = mom.integrate(6, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AqqQNSEven TruncAs0 Mellin moment N=8");
    refVal = 1.;
    errBound = 1. * 3.e-11;
    res = mom.integrate(8, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AqqQNSEven TruncAs0 Mellin moment N=10");
    refVal = 1.;
    errBound = 1. * 3.e-11;
    res = mom.integrate(10, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

}


TEST(AqqQNSEvenNspace,TruncAs1Moments)
{
  double refVal = 0., errBound = 0.;
  auto rpd = AqqQNSEven.truncate(1);
  integration_engine_gsl engine;
  auto mom = make_mellin_moment(rpd, engine, testAs, testLM, testNF);
  std::tuple<integration_status, double, double> res;

  {
    SCOPED_TRACE("AqqQNSEven TruncAs1 Mellin moment N=2");
    refVal = 1.;
    errBound = 1. * 3.e-11;
    res = mom.integrate(2, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AqqQNSEven TruncAs1 Mellin moment N=4");
    refVal = 1.;
    errBound = 1. * 3.e-11;
    res = mom.integrate(4, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AqqQNSEven TruncAs1 Mellin moment N=6");
    refVal = 1.;
    errBound = 1. * 3.e-11;
    res = mom.integrate(6, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AqqQNSEven TruncAs1 Mellin moment N=8");
    refVal = 1.;
    errBound = 1. * 3.e-11;
    res = mom.integrate(8, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AqqQNSEven TruncAs1 Mellin moment N=10");
    refVal = 1.;
    errBound = 1. * 3.e-11;
    res = mom.integrate(10, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

}


TEST(AqqQNSEvenNspace,TruncAs2Moments)
{
  double refVal = 0., errBound = 0.;
  auto rpd = AqqQNSEven.truncate(2);
  integration_engine_gsl engine;
  auto mom = make_mellin_moment(rpd, engine, testAs, testLM, testNF);
  std::tuple<integration_status, double, double> res;

  {
    SCOPED_TRACE("AqqQNSEven TruncAs2 Mellin moment N=2");
    refVal = -0.0617283950617283950617283950617;
    errBound = 0.0617283950617283950617283950617 * 3.e-11;
    res = mom.integrate(2, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AqqQNSEven TruncAs2 Mellin moment N=4");
    refVal = -0.988148919753086419753086419753;
    errBound = 0.988148919753086419753086419753 * 3.e-11;
    res = mom.integrate(4, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AqqQNSEven TruncAs2 Mellin moment N=6");
    refVal = -1.51807314574859926333849236344;
    errBound = 1.51807314574859926333849236344 * 3.e-11;
    res = mom.integrate(6, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AqqQNSEven TruncAs2 Mellin moment N=8");
    refVal = -1.89176165229200799138650822534;
    errBound = 1.89176165229200799138650822534 * 3.e-11;
    res = mom.integrate(8, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AqqQNSEven TruncAs2 Mellin moment N=10");
    refVal = -2.18060176847780770650295972924;
    errBound = 2.18060176847780770650295972924 * 3.e-11;
    res = mom.integrate(10, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

}
