/******************************************************************************
 * This file is part of libome                                                *
 * Copyright (C) 2025 Arnd Behring, Kay Schoenwald                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 ******************************************************************************/

#include <limits>
#include <gtest/gtest.h>
#include <ome/integration_engine_gsl.h>
#include <ome/mellin.h>
#include <ome/polAqqQNSEven.h>

using namespace ome;

const double eps = std::numeric_limits<double>::epsilon();
const double testAs = 0.25;
const double testLM = -5.;
const double testNF = 3.;



TEST(polAqqQNSEvenXspace,FullPowersOf4RegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto& eval = polAqqQNSEven_reg;

  {
    testVal  = eval(testAs,testLM,testNF,9.31322574615478515625e-10);
    refVal   = 3296.04537381517339333655888175;
    errBound = 3296.04537381517339333655888175 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, reg, x = 4^-15 = 9.31322574615478515625e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.7252902984619140625e-9);
    refVal   = 2626.51255671304019080229493741;
    errBound = 2626.51255671304019080229493741 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, reg, x = 4^-14 = 3.7252902984619140625e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.490116119384765625e-8);
    refVal   = 2062.95953460324253389856389005;
    errBound = 2062.95953460324253389856389005 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, reg, x = 4^-13 = 1.490116119384765625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,5.9604644775390625e-8);
    refVal   = 1593.84115551568814963638263612;
    errBound = 1593.84115551568814963638263612 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, reg, x = 4^-12 = 5.9604644775390625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,2.384185791015625e-7);
    refVal   = 1208.25066074472474306263839093;
    errBound = 1208.25066074472474306263839093 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, reg, x = 4^-11 = 2.384185791015625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,9.5367431640625e-7);
    refVal   = 895.91973923987743790628181377;
    errBound = 895.91973923987743790628181377 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, reg, x = 4^-10 = 9.5367431640625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.814697265625e-6);
    refVal   = 647.218652022395439497940051664;
    errBound = 647.218652022395439497940051664 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, reg, x = 4^-9 = 3.814697265625e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0000152587890625);
    refVal   = 453.156477520676356279561484685;
    errBound = 453.156477520676356279561484685 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, reg, x = 4^-8 = 0.0000152587890625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00006103515625);
    refVal   = 305.381456061735833845421369419;
    errBound = 305.381456061735833845421369419 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, reg, x = 4^-7 = 0.00006103515625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.000244140625);
    refVal   = 196.180925876489050568145818666;
    errBound = 196.180925876489050568145818666 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, reg, x = 4^-6 = 0.000244140625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0009765625);
    refVal   = 118.478387520361681482228903029;
    errBound = 118.478387520361681482228903029 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, reg, x = 4^-5 = 0.0009765625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00390625);
    refVal   = 65.8188784225058814958762683238;
    errBound = 65.8188784225058814958762683238 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, reg, x = 4^-4 = 0.00390625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.015625);
    refVal   = 32.3151241141230639313170540425;
    errBound = 32.3151241141230639313170540425 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, reg, x = 4^-3 = 0.015625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0625);
    refVal   = 12.4626999653935567448836053013;
    errBound = 12.4626999653935567448836053013 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, reg, x = 4^-2 = 0.0625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9375);
    refVal   = -19.533499707492819924166041588;
    errBound = 19.533499707492819924166041588 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, reg, x = 1-4^-2 = 0.9375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.984375);
    refVal   = -25.4738115673800742738168741014;
    errBound = 25.4738115673800742738168741014 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, reg, x = 1-4^-3 = 0.984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99609375);
    refVal   = -31.4114365190396091087493046279;
    errBound = 31.4114365190396091087493046279 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, reg, x = 1-4^-4 = 0.99609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9990234375);
    refVal   = -37.9881861315671808637795040215;
    errBound = 37.9881861315671808637795040215 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, reg, x = 1-4^-5 = 0.9990234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999755859375);
    refVal   = -45.5693036062218496499458888526;
    errBound = 45.5693036062218496499458888526 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, reg, x = 1-4^-6 = 0.999755859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99993896484375);
    refVal   = -54.4456195970858380912400625487;
    errBound = 54.4456195970858380912400625487 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, reg, x = 1-4^-7 = 0.99993896484375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999847412109375);
    refVal   = -64.8879313361013790381369055521;
    errBound = 64.8879313361013790381369055521 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, reg, x = 1-4^-8 = 0.9999847412109375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999996185302734375);
    refVal   = -77.1615433759623133029418607878;
    errBound = 77.1615433759623133029418607878 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, reg, x = 1-4^-9 = 0.999996185302734375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999904632568359375);
    refVal   = -91.5302163428994155674676009464;
    errBound = 91.5302163428994155674676009464 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, reg, x = 1-4^-10 = 0.99999904632568359375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999997615814208984375);
    refVal   = -108.25726644931167033232890761;
    errBound = 108.25726644931167033232890761 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, reg, x = 1-4^-11 = 0.9999997615814208984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999940395355224609375);
    refVal   = -127.60587955596877479006020191;
    errBound = 127.60587955596877479006020191 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, reg, x = 1-4^-12 = 0.999999940395355224609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999998509883880615234375);
    refVal   = -149.83920279915934867937229218;
    errBound = 149.83920279915934867937229218 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, reg, x = 1-4^-13 = 0.99999998509883880615234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999962747097015380859375);
    refVal   = -175.220371728624713668527837643;
    errBound = 175.220371728624713668527837643 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, reg, x = 1-4^-14 = 0.9999999962747097015380859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999068677425384521484375);
    refVal   = -204.012518418332054825472726641;
    errBound = 204.012518418332054825472726641 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, reg, x = 1-4^-15 = 0.999999999068677425384521484375";
  }
}

TEST(polAqqQNSEvenXspace,FullPowersOf4PlusEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;
  auto& eval = polAqqQNSEven_plus;
  {
    testVal  = eval(testAs,testLM,testNF,9.31322574615478515625e-10);
    refVal   = 19.8478639549625839237424922507;
    errBound = 19.8478639549625839237424922507 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, plus, x = 4^-15 = 9.31322574615478515625e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.7252902984619140625e-9);
    refVal   = 19.8478640104168754077863490609;
    errBound = 19.8478640104168754077863490609 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, plus, x = 4^-14 = 3.7252902984619140625e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.490116119384765625e-8);
    refVal   = 19.8478642322340444427118335806;
    errBound = 19.8478642322340444427118335806 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, plus, x = 4^-13 = 1.490116119384765625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,5.9604644775390625e-8);
    refVal   = 19.8478651195027701624175971435;
    errBound = 19.8478651195027701624175971435 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, plus, x = 4^-12 = 5.9604644775390625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,2.384185791015625e-7);
    refVal   = 19.8478686685784663214880366929;
    errBound = 19.8478686685784663214880366929 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, plus, x = 4^-11 = 2.384185791015625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,9.5367431640625e-7);
    refVal   = 19.847882864893943453643334231;
    errBound = 19.847882864893943453643334231 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, plus, x = 4^-10 = 9.5367431640625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.814697265625e-6);
    refVal   = 19.8479396503589326788280243906;
    errBound = 19.8479396503589326788280243906 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, plus, x = 4^-9 = 3.814697265625e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0000152587890625);
    refVal   = 19.8481667954682295308843167689;
    errBound = 19.8481667954682295308843167689 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, plus, x = 4^-8 = 0.0000152587890625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00006103515625);
    refVal   = 19.8490754278979799534045000832;
    errBound = 19.8490754278979799534045000832 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, plus, x = 4^-7 = 0.00006103515625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.000244140625);
    refVal   = 19.8527107896979612684142811272;
    errBound = 19.8527107896979612684142811272 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, plus, x = 4^-6 = 0.000244140625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0009765625);
    refVal   = 19.8672655630041914452972827996;
    errBound = 19.8672655630041914452972827996 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, plus, x = 4^-5 = 0.0009765625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00390625);
    refVal   = 19.9256986970130273024893336313;
    errBound = 19.9256986970130273024893336313 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, plus, x = 4^-4 = 0.00390625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.015625);
    refVal   = 20.1629093957869919132332542698;
    errBound = 20.1629093957869919132332542698 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, plus, x = 4^-3 = 0.015625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0625);
    refVal   = 21.1710548655763415088949169833;
    errBound = 21.1710548655763415088949169833 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, plus, x = 4^-2 = 0.0625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9375);
    refVal   = 317.565822983645122633423754749;
    errBound = 317.565822983645122633423754749 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, plus, x = 1-4^-2 = 0.9375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.984375);
    refVal   = 1270.263291934580490533695019;
    errBound = 1270.263291934580490533695019 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, plus, x = 1-4^-3 = 0.984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99609375);
    refVal   = 5081.05316773832196213478007599;
    errBound = 5081.05316773832196213478007599 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, plus, x = 1-4^-4 = 0.99609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9990234375);
    refVal   = 20324.212670953287848539120304;
    errBound = 20324.212670953287848539120304 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, plus, x = 1-4^-5 = 0.9990234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999755859375);
    refVal   = 81296.8506838131513941564812158;
    errBound = 81296.8506838131513941564812158 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, plus, x = 1-4^-6 = 0.999755859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99993896484375);
    refVal   = 325187.402735252605576625924863;
    errBound = 325187.402735252605576625924863 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, plus, x = 1-4^-7 = 0.99993896484375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999847412109375);
    refVal   = 1.30074961094101042230650369945e6;
    errBound = 1.30074961094101042230650369945e6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, plus, x = 1-4^-8 = 0.9999847412109375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999996185302734375);
    refVal   = 5.20299844376404168922601479781e6;
    errBound = 5.20299844376404168922601479781e6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, plus, x = 1-4^-9 = 0.999996185302734375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999904632568359375);
    refVal   = 2.08119937750561667569040591912e7;
    errBound = 2.08119937750561667569040591912e7 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, plus, x = 1-4^-10 = 0.99999904632568359375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999997615814208984375);
    refVal   = 8.3247975100224667027616236765e7;
    errBound = 8.3247975100224667027616236765e7 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, plus, x = 1-4^-11 = 0.9999997615814208984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999940395355224609375);
    refVal   = 3.3299190040089866811046494706e8;
    errBound = 3.3299190040089866811046494706e8 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, plus, x = 1-4^-12 = 0.999999940395355224609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999998509883880615234375);
    refVal   = 1.33196760160359467244185978824e9;
    errBound = 1.33196760160359467244185978824e9 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, plus, x = 1-4^-13 = 0.99999998509883880615234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999962747097015380859375);
    refVal   = 5.32787040641437868976743915296e9;
    errBound = 5.32787040641437868976743915296e9 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, plus, x = 1-4^-14 = 0.9999999962747097015380859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999068677425384521484375);
    refVal   = 2.13114816256575147590697566118e10;
    errBound = 2.13114816256575147590697566118e10 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullPowersOf4, plus, x = 1-4^-15 = 0.999999999068677425384521484375";
  }
}



TEST(polAqqQNSEvenXspace,FullNormalDoublesRegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto& eval = polAqqQNSEven_reg;

  {
    testVal  = eval(testAs,testLM,testNF,1.e-10);
    refVal   = 4627.3757936598036485792914346;
    errBound = 4627.3757936598036485792914346 * std::max(eps/(1.-1.e-10),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, reg, x = 10^-10 = 1.e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-9);
    refVal   = 3258.91249875239143576610693224;
    errBound = 3258.91249875239143576610693224 * std::max(eps/(1.-1.e-9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, reg, x = 10^-9 = 1.e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-8);
    refVal   = 2214.92851528303766407803176483;
    errBound = 2214.92851528303766407803176483 * std::max(eps/(1.-1.e-8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, reg, x = 10^-8 = 1.e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-7);
    refVal   = 1440.71250612257246852770020941;
    errBound = 1440.71250612257246852770020941 * std::max(eps/(1.-1.e-7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, reg, x = 10^-7 = 1.e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-6);
    refVal   = 886.411930043542775687533518323;
    errBound = 886.411930043542775687533518323 * std::max(eps/(1.-1.e-6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, reg, x = 10^-6 = 1.e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-5);
    refVal   = 507.03393146168643763336542804;
    errBound = 507.03393146168643763336542804 * std::max(eps/(1.-1.e-5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, reg, x = 10^-5 = 1.e-5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0001);
    refVal   = 262.447944337924715004325005455;
    errBound = 262.447944337924715004325005455 * std::max(eps/(1.-0.0001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, reg, x = 10^-4 = 0.0001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.001);
    refVal   = 117.384478301378045601730012454;
    errBound = 117.384478301378045601730012454 * std::max(eps/(1.-0.001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, reg, x = 10^-3 = 0.001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.01);
    refVal   = 41.3501687617931037093529611992;
    errBound = 41.3501687617931037093529611992 * std::max(eps/(1.-0.01),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, reg, x = 10^-2 = 0.01";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.05);
    refVal   = 14.9638654672620649628108337991;
    errBound = 14.9638654672620649628108337991 * std::max(eps/(1.-0.05),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, reg, x = 0.05 = 0.05";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.1);
    refVal   = 7.8389174993457382803695550162;
    errBound = 7.8389174993457382803695550162 * std::max(eps/(1.-0.1),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, reg, x = 0.1 = 0.1";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.15);
    refVal   = 4.37796692658074167521014437538;
    errBound = 4.37796692658074167521014437538 * std::max(eps/(1.-0.15),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, reg, x = 0.15 = 0.15";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.2);
    refVal   = 2.09277424866297614579396014491;
    errBound = 2.09277424866297614579396014491 * std::max(eps/(1.-0.2),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, reg, x = 0.2 = 0.2";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.25);
    refVal   = 0.342738744157654796021735368469;
    errBound = 0.342738744157654796021735368469 * std::max(eps/(1.-0.25),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, reg, x = 0.25 = 0.25";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.3);
    refVal   = -1.12204901138288742310449928615;
    errBound = 1.12204901138288742310449928615 * std::max(eps/(1.-0.3),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, reg, x = 0.3 = 0.3";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.35);
    refVal   = -2.42258722293067980697621017047;
    errBound = 2.42258722293067980697621017047 * std::max(eps/(1.-0.35),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, reg, x = 0.35 = 0.35";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.4);
    refVal   = -3.62647361268547902300482978342;
    errBound = 3.62647361268547902300482978342 * std::max(eps/(1.-0.4),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, reg, x = 0.4 = 0.4";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.45);
    refVal   = -4.77599609071370758655809954043;
    errBound = 4.77599609071370758655809954043 * std::max(eps/(1.-0.45),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, reg, x = 0.45 = 0.45";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.5);
    refVal   = -5.90055538725432210067313008755;
    errBound = 5.90055538725432210067313008755 * std::max(eps/(1.-0.5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, reg, x = 0.5 = 0.5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.55);
    refVal   = -7.02303033275694818277575834551;
    errBound = 7.02303033275694818277575834551 * std::max(eps/(1.-0.55),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, reg, x = 0.55 = 0.55";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.6);
    refVal   = -8.16364474289485451429605280997;
    errBound = 8.16364474289485451429605280997 * std::max(eps/(1.-0.6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, reg, x = 0.6 = 0.6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.65);
    refVal   = -9.3429728608690938637649500971;
    errBound = 9.3429728608690938637649500971 * std::max(eps/(1.-0.65),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, reg, x = 0.65 = 0.65";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.7);
    refVal   = -10.5851886866028829040909006686;
    errBound = 10.5851886866028829040909006686 * std::max(eps/(1.-0.7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, reg, x = 0.7 = 0.7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.75);
    refVal   = -11.9229391643875804851994872092;
    errBound = 11.9229391643875804851994872092 * std::max(eps/(1.-0.75),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, reg, x = 0.75 = 0.75";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.8);
    refVal   = -13.4067765366614478965360982304;
    errBound = 13.4067765366614478965360982304 * std::max(eps/(1.-0.8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, reg, x = 0.8 = 0.8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.85);
    refVal   = -15.1277941567597744513600342166;
    errBound = 15.1277941567597744513600342166 * std::max(eps/(1.-0.85),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, reg, x = 0.85 = 0.85";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9);
    refVal   = -17.2877502083643891205404526903;
    errBound = 17.2877502083643891205404526903 * std::max(eps/(1.-0.9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, reg, x = 0.9 = 0.9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.95);
    refVal   = -20.5362426960939354300841342341;
    errBound = 20.5362426960939354300841342341 * std::max(eps/(1.-0.95),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, reg, x = 0.95 = 0.95";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99);
    refVal   = -27.3449522535882835993822596828;
    errBound = 27.3449522535882835993822596828 * std::max(eps/(1.-0.99),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, reg, x = 1-10^-2 = 0.99";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999);
    refVal   = -37.8681119000460366046537987923;
    errBound = 37.8681119000460366046537987923 * std::max(eps/(1.-0.999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, reg, x = 1-10^-3 = 0.999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999);
    refVal   = -51.1186765072087164685341799815;
    errBound = 51.1186765072087164685341799815 * std::max(eps/(1.-0.9999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, reg, x = 1-10^-4 = 0.9999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999);
    refVal   = -68.4229761863136300472023204766;
    errBound = 68.4229761863136300472023204766 * std::max(eps/(1.-0.99999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, reg, x = 1-10^-5 = 0.99999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999);
    refVal   = -91.0011149317460472085185846922;
    errBound = 91.0011149317460472085185846922 * std::max(eps/(1.-0.999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, reg, x = 1-10^-6 = 0.999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999);
    refVal   = -120.060649312787370373306007226;
    errBound = 120.060649312787370373306007226 * std::max(eps/(1.-0.9999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, reg, x = 1-10^-7 = 0.9999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999999);
    refVal   = -156.807555772841716971706520241;
    errBound = 156.807555772841716971706520241 * std::max(eps/(1.-0.99999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, reg, x = 1-10^-8 = 0.99999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999);
    refVal   = -202.447602255404363621853154814;
    errBound = 202.447602255404363621853154814 * std::max(eps/(1.-0.999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, reg, x = 1-10^-9 = 0.999999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999999);
    refVal   = -258.186528602791732746182909573;
    errBound = 258.186528602791732746182909573 * std::max(eps/(1.-0.9999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, reg, x = 1-10^-10 = 0.9999999999";
  }
}

TEST(polAqqQNSEvenXspace,FullNormalDoublesPlusEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;
  auto& eval = polAqqQNSEven_plus;
  {
    testVal  = eval(testAs,testLM,testNF,1.e-10);
    refVal   = 19.8478639384626065584352453277;
    errBound = 19.8478639384626065584352453277 * std::max(eps/(1.-1.e-10),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, plus, x = 10^-10 = 1.e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-9);
    refVal   = 19.8478639563256841209146687927;
    errBound = 19.8478639563256841209146687927 * std::max(eps/(1.-1.e-9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, plus, x = 10^-9 = 1.e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-8);
    refVal   = 19.8478641349564615141535998134;
    errBound = 19.8478641349564615141535998134 * std::max(eps/(1.-1.e-8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, plus, x = 10^-8 = 1.e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-7);
    refVal   = 19.8478659212644122910302137749;
    errBound = 19.8478659212644122910302137749 * std::max(eps/(1.-1.e-7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, plus, x = 10^-7 = 1.e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-6);
    refVal   = 19.8478837843616045261935108653;
    errBound = 19.8478837843616045261935108653 * std::max(eps/(1.-1.e-6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, plus, x = 10^-6 = 1.e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-5);
    refVal   = 19.8480624171019911845008296801;
    errBound = 19.8480624171019911845008296801 * std::max(eps/(1.-1.e-5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, plus, x = 10^-5 = 1.e-5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0001);
    refVal   = 19.8498489213699571603050151733;
    errBound = 19.8498489213699571603050151733 * std::max(eps/(1.-0.0001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, plus, x = 10^-4 = 0.0001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.001);
    refVal   = 19.8677316681459661307197043762;
    errBound = 19.8677316681459661307197043762 * std::max(eps/(1.-0.001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, plus, x = 10^-3 = 0.001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.01);
    refVal   = 20.0483474105836567319080653251;
    errBound = 20.0483474105836567319080653251 * std::max(eps/(1.-0.01),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, plus, x = 10^-2 = 0.01";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.05);
    refVal   = 20.8924883541871791206199838651;
    errBound = 20.8924883541871791206199838651 * std::max(eps/(1.-0.05),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, plus, x = 0.05 = 0.05";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.1);
    refVal   = 22.0531821516420224050988718576;
    errBound = 22.0531821516420224050988718576 * std::max(eps/(1.-0.1),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, plus, x = 0.1 = 0.1";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.15);
    refVal   = 23.3504281605621413701046878492;
    errBound = 23.3504281605621413701046878492 * std::max(eps/(1.-0.15),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, plus, x = 0.15 = 0.15";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.2);
    refVal   = 24.8098299205972752057362308398;
    errBound = 24.8098299205972752057362308398 * std::max(eps/(1.-0.2),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, plus, x = 0.2 = 0.2";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.25);
    refVal   = 26.4638185819704268861186462291;
    errBound = 26.4638185819704268861186462291 * std::max(eps/(1.-0.25),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, plus, x = 0.25 = 0.25";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.3);
    refVal   = 28.3540913378254573779842638169;
    errBound = 28.3540913378254573779842638169 * std::max(eps/(1.-0.3),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, plus, x = 0.3 = 0.3";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.35);
    refVal   = 30.5351752868889540993676687259;
    errBound = 30.5351752868889540993676687259 * std::max(eps/(1.-0.35),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, plus, x = 0.35 = 0.35";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.4);
    refVal   = 33.0797732274630336076483077864;
    errBound = 33.0797732274630336076483077864 * std::max(eps/(1.-0.4),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, plus, x = 0.4 = 0.4";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.45);
    refVal   = 36.0870253390505821174345175851;
    errBound = 36.0870253390505821174345175851 * std::max(eps/(1.-0.45),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, plus, x = 0.45 = 0.45";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.5);
    refVal   = 39.6957278729556403291779693437;
    errBound = 39.6957278729556403291779693437 * std::max(eps/(1.-0.5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, plus, x = 0.5 = 0.5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.55);
    refVal   = 44.1063643032840448101977437152;
    errBound = 44.1063643032840448101977437152 * std::max(eps/(1.-0.55),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, plus, x = 0.55 = 0.55";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.6);
    refVal   = 49.6196598411945504114724616796;
    errBound = 49.6196598411945504114724616796 * std::max(eps/(1.-0.6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, plus, x = 0.6 = 0.6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.65);
    refVal   = 56.7081826756509147559685276338;
    errBound = 56.7081826756509147559685276338 * std::max(eps/(1.-0.65),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, plus, x = 0.65 = 0.65";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.7);
    refVal   = 66.1595464549260672152966155728;
    errBound = 66.1595464549260672152966155728 * std::max(eps/(1.-0.7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, plus, x = 0.7 = 0.7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.75);
    refVal   = 79.3914557459112806583559386873;
    errBound = 79.3914557459112806583559386873 * std::max(eps/(1.-0.75),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, plus, x = 0.75 = 0.75";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.8);
    refVal   = 99.2393196823891008229449233592;
    errBound = 99.2393196823891008229449233592 * std::max(eps/(1.-0.8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, plus, x = 0.8 = 0.8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.85);
    refVal   = 132.319092909852134430593231146;
    errBound = 132.319092909852134430593231146 * std::max(eps/(1.-0.85),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, plus, x = 0.85 = 0.85";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9);
    refVal   = 198.478639364778201645889846718;
    errBound = 198.478639364778201645889846718 * std::max(eps/(1.-0.9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, plus, x = 0.9 = 0.9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.95);
    refVal   = 396.957278729556403291779693437;
    errBound = 396.957278729556403291779693437 * std::max(eps/(1.-0.95),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, plus, x = 0.95 = 0.95";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99);
    refVal   = 1984.78639364778201645889846718;
    errBound = 1984.78639364778201645889846718 * std::max(eps/(1.-0.99),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, plus, x = 1-10^-2 = 0.99";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999);
    refVal   = 19847.8639364778201645889846718;
    errBound = 19847.8639364778201645889846718 * std::max(eps/(1.-0.999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, plus, x = 1-10^-3 = 0.999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999);
    refVal   = 198478.639364778201645889846718;
    errBound = 198478.639364778201645889846718 * std::max(eps/(1.-0.9999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, plus, x = 1-10^-4 = 0.9999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999);
    refVal   = 1.98478639364778201645889846718e6;
    errBound = 1.98478639364778201645889846718e6 * std::max(eps/(1.-0.99999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, plus, x = 1-10^-5 = 0.99999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999);
    refVal   = 1.98478639364778201645889846718e7;
    errBound = 1.98478639364778201645889846718e7 * std::max(eps/(1.-0.999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, plus, x = 1-10^-6 = 0.999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999);
    refVal   = 1.98478639364778201645889846718e8;
    errBound = 1.98478639364778201645889846718e8 * std::max(eps/(1.-0.9999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, plus, x = 1-10^-7 = 0.9999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999999);
    refVal   = 1.98478639364778201645889846718e9;
    errBound = 1.98478639364778201645889846718e9 * std::max(eps/(1.-0.99999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, plus, x = 1-10^-8 = 0.99999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999);
    refVal   = 1.98478639364778201645889846718e10;
    errBound = 1.98478639364778201645889846718e10 * std::max(eps/(1.-0.999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, plus, x = 1-10^-9 = 0.999999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999999);
    refVal   = 1.98478639364778201645889846718e11;
    errBound = 1.98478639364778201645889846718e11 * std::max(eps/(1.-0.9999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, FullNormalDoubles, plus, x = 1-10^-10 = 0.9999999999";
  }
}



TEST(polAqqQNSEvenXspace,FullDeltaEvaluation)
{
  auto& eval = polAqqQNSEven_delta;
  EXPECT_DOUBLE_EQ(eval(testAs,testLM,testNF),12.1501748418929015300593526677)
    << "^-- polAqqQNSEven, Full, delta";
}







TEST(polAqqQNSEvenXspace,TruncAs0NormalDoublesDeltaEvaluation)
{
  auto eval = polAqqQNSEven_delta.truncate(0);
  EXPECT_DOUBLE_EQ(eval(testAs,testLM,testNF),1.)
    << "^-- polAqqQNSEven, TruncAs0NormalDoubles, delta";
}







TEST(polAqqQNSEvenXspace,TruncAs1NormalDoublesDeltaEvaluation)
{
  auto eval = polAqqQNSEven_delta.truncate(1);
  EXPECT_DOUBLE_EQ(eval(testAs,testLM,testNF),1.)
    << "^-- polAqqQNSEven, TruncAs1NormalDoubles, delta";
}



TEST(polAqqQNSEvenXspace,TruncAs2PowersOf4RegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto eval = polAqqQNSEven_reg.truncate(2);

  {
    testVal  = eval(testAs,testLM,testNF,9.31322574615478515625e-10);
    refVal   = 24.2945391531051607896541861718;
    errBound = 24.2945391531051607896541861718 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, reg, x = 4^-15 = 9.31322574615478515625e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.7252902984619140625e-9);
    refVal   = 21.950577165068803555966157488;
    errBound = 21.950577165068803555966157488 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, reg, x = 4^-14 = 3.7252902984619140625e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.490116119384765625e-8);
    refVal   = 19.7133826179923339322847932777;
    errBound = 19.7133826179923339322847932777 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, reg, x = 4^-13 = 1.490116119384765625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,5.9604644775390625e-8);
    refVal   = 17.5829557667636721009561321108;
    errBound = 17.5829557667636721009561321108 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, reg, x = 4^-12 = 5.9604644775390625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,2.384185791015625e-7);
    refVal   = 15.5592974691149663221180610747;
    errBound = 15.5592974691149663221180610747 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, reg, x = 4^-11 = 2.384185791015625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,9.5367431640625e-7);
    refVal   = 13.6424105516611732880307771909;
    errBound = 13.6424105516611732880307771909 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, reg, x = 4^-10 = 9.5367431640625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.814697265625e-6);
    refVal   = 11.8323040755482182127343650263;
    errBound = 11.8323040755482182127343650263 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, reg, x = 4^-9 = 3.814697265625e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0000152587890625);
    refVal   = 10.1290059929594676451689047194;
    errBound = 10.1290059929594676451689047194 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, reg, x = 4^-8 = 0.0000152587890625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00006103515625);
    refVal   = 8.5325977125286692136080901449;
    errBound = 8.5325977125286692136080901449 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, reg, x = 4^-7 = 0.00006103515625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.000244140625);
    refVal   = 7.04329452172840778557083330018;
    errBound = 7.04329452172840778557083330018 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, reg, x = 4^-6 = 0.000244140625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0009765625);
    refVal   = 5.66156327512726792864800589954;
    errBound = 5.66156327512726792864800589954 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, reg, x = 4^-5 = 0.0009765625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00390625);
    refVal   = 4.38790550536101365575516619009;
    errBound = 4.38790550536101365575516619009 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, reg, x = 4^-4 = 0.00390625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.015625);
    refVal   = 3.21993361190786973785471883912;
    errBound = 3.21993361190786973785471883912 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, reg, x = 4^-3 = 0.015625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0625);
    refVal   = 2.13659810997390165973939871603;
    errBound = 2.13659810997390165973939871603 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, reg, x = 4^-2 = 0.0625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9375);
    refVal   = -0.268362475445021794896393223338;
    errBound = 0.268362475445021794896393223338 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, reg, x = 1-4^-2 = 0.9375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.984375);
    refVal   = -0.3266713135058627817810942937;
    errBound = 0.3266713135058627817810942937 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, reg, x = 1-4^-3 = 0.984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99609375);
    refVal   = -0.340946513246320684792526634059;
    errBound = 0.340946513246320684792526634059 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, reg, x = 1-4^-4 = 0.99609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9990234375);
    refVal   = -0.34449709175903860580951890324;
    errBound = 0.34449709175903860580951890324 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, reg, x = 1-4^-5 = 0.9990234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999755859375);
    refVal   = -0.345383607307621750884392650659;
    errBound = 0.345383607307621750884392650659 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, reg, x = 1-4^-6 = 0.999755859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99993896484375);
    refVal   = -0.345605165778065901050331736072;
    errBound = 0.345605165778065901050331736072 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, reg, x = 1-4^-7 = 0.99993896484375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999847412109375);
    refVal   = -0.345660550996982297006902759782;
    errBound = 0.345660550996982297006902759782 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, reg, x = 1-4^-8 = 0.9999847412109375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999996185302734375);
    refVal   = -0.345674397026829661536539041158;
    errBound = 0.345674397026829661536539041158 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, reg, x = 1-4^-9 = 0.999996185302734375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999904632568359375);
    refVal   = -0.345677858517111967299608103272;
    errBound = 0.345677858517111967299608103272 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, reg, x = 1-4^-10 = 0.99999904632568359375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999997615814208984375);
    refVal   = -0.345678723888608831733060846091;
    errBound = 0.345678723888608831733060846091 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, reg, x = 1-4^-11 = 0.9999997615814208984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999940395355224609375);
    refVal   = -0.345678940231415940980860170939;
    errBound = 0.345678940231415940980860170939 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, reg, x = 1-4^-12 = 0.999999940395355224609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999998509883880615234375);
    refVal   = -0.345678994317113524116210587615;
    errBound = 0.345678994317113524116210587615 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, reg, x = 1-4^-13 = 0.99999998509883880615234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999962747097015380859375);
    refVal   = -0.345679007838537657764044881895;
    errBound = 0.345679007838537657764044881895 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, reg, x = 1-4^-14 = 0.9999999962747097015380859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999068677425384521484375);
    refVal   = -0.345679011218893674792503782246;
    errBound = 0.345679011218893674792503782246 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, reg, x = 1-4^-15 = 0.999999999068677425384521484375";
  }
}

TEST(polAqqQNSEvenXspace,TruncAs2PowersOf4PlusEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;
  auto eval = polAqqQNSEven_plus.truncate(2);
  {
    testVal  = eval(testAs,testLM,testNF,9.31322574615478515625e-10);
    refVal   = 1.27160493945587932437941415068;
    errBound = 1.27160493945587932437941415068 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, plus, x = 4^-15 = 9.31322574615478515625e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.7252902984619140625e-9);
    refVal   = 1.27160494300870249593813948108;
    errBound = 1.27160494300870249593813948108 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, plus, x = 4^-14 = 3.7252902984619140625e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.490116119384765625e-8);
    refVal   = 1.27160495721999538070250915835;
    errBound = 1.27160495721999538070250915835 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, plus, x = 4^-13 = 1.490116119384765625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,5.9604644775390625e-8);
    refVal   = 1.27160501406517009623166793273;
    errBound = 1.27160501406517009623166793273 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, plus, x = 4^-12 = 5.9604644775390625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,2.384185791015625e-7);
    refVal   = 1.27160524144591978190711202291;
    errBound = 1.27160524144591978190711202291 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, plus, x = 4^-11 = 2.384185791015625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,9.5367431640625e-7);
    refVal   = 1.2716061509697317023132216768;
    errBound = 1.2716061509697317023132216768 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, plus, x = 4^-10 = 9.5367431640625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.814697265625e-6);
    refVal   = 1.27160978907799027606410091034;
    errBound = 1.27160978907799027606410091034 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, plus, x = 4^-9 = 3.814697265625e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0000152587890625);
    refVal   = 1.27162434171920197199310141504;
    errBound = 1.27162434171920197199310141504 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, plus, x = 4^-8 = 0.0000152587890625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00006103515625);
    refVal   = 1.27168255561508730444008883041;
    errBound = 1.27168255561508730444008883041 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, plus, x = 4^-7 = 0.00006103515625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.000244140625);
    refVal   = 1.27191546450805710064969324229;
    errBound = 1.27191546450805710064969324229 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, plus, x = 4^-6 = 0.000244140625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0009765625);
    refVal   = 1.27284795385153808092876193235;
    errBound = 1.27284795385153808092876193235 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, plus, x = 4^-5 = 0.0009765625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00390625);
    refVal   = 1.27659162430404260469619946744;
    errBound = 1.27659162430404260469619946744 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, plus, x = 4^-4 = 0.00390625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.015625);
    refVal   = 1.2917891436409954928473446992;
    errBound = 1.2917891436409954928473446992 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, plus, x = 4^-3 = 0.015625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0625);
    refVal   = 1.35637860082304526748971193416;
    errBound = 1.35637860082304526748971193416 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, plus, x = 4^-2 = 0.0625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9375);
    refVal   = 20.3456790123456790123456790123;
    errBound = 20.3456790123456790123456790123 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, plus, x = 1-4^-2 = 0.9375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.984375);
    refVal   = 81.3827160493827160493827160494;
    errBound = 81.3827160493827160493827160494 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, plus, x = 1-4^-3 = 0.984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99609375);
    refVal   = 325.530864197530864197530864198;
    errBound = 325.530864197530864197530864198 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, plus, x = 1-4^-4 = 0.99609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9990234375);
    refVal   = 1302.12345679012345679012345679;
    errBound = 1302.12345679012345679012345679 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, plus, x = 1-4^-5 = 0.9990234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999755859375);
    refVal   = 5208.49382716049382716049382716;
    errBound = 5208.49382716049382716049382716 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, plus, x = 1-4^-6 = 0.999755859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99993896484375);
    refVal   = 20833.9753086419753086419753086;
    errBound = 20833.9753086419753086419753086 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, plus, x = 1-4^-7 = 0.99993896484375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999847412109375);
    refVal   = 83335.9012345679012345679012346;
    errBound = 83335.9012345679012345679012346 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, plus, x = 1-4^-8 = 0.9999847412109375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999996185302734375);
    refVal   = 333343.604938271604938271604938;
    errBound = 333343.604938271604938271604938 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, plus, x = 1-4^-9 = 0.999996185302734375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999904632568359375);
    refVal   = 1.33337441975308641975308641975e6;
    errBound = 1.33337441975308641975308641975e6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, plus, x = 1-4^-10 = 0.99999904632568359375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999997615814208984375);
    refVal   = 5.33349767901234567901234567901e6;
    errBound = 5.33349767901234567901234567901e6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, plus, x = 1-4^-11 = 0.9999997615814208984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999940395355224609375);
    refVal   = 2.1333990716049382716049382716e7;
    errBound = 2.1333990716049382716049382716e7 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, plus, x = 1-4^-12 = 0.999999940395355224609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999998509883880615234375);
    refVal   = 8.53359628641975308641975308642e7;
    errBound = 8.53359628641975308641975308642e7 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, plus, x = 1-4^-13 = 0.99999998509883880615234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999962747097015380859375);
    refVal   = 3.41343851456790123456790123457e8;
    errBound = 3.41343851456790123456790123457e8 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, plus, x = 1-4^-14 = 0.9999999962747097015380859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999068677425384521484375);
    refVal   = 1.36537540582716049382716049383e9;
    errBound = 1.36537540582716049382716049383e9 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2PowersOf4, plus, x = 1-4^-15 = 0.999999999068677425384521484375";
  }
}



TEST(polAqqQNSEvenXspace,TruncAs2NormalDoublesRegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto eval = polAqqQNSEven_reg.truncate(2);

  {
    testVal  = eval(testAs,testLM,testNF,1.e-10);
    refVal   = 28.2917178181689775076785574172;
    errBound = 28.2917178181689775076785574172 * std::max(eps/(1.-1.e-10),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, reg, x = 10^-10 = 1.e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-9);
    refVal   = 24.1716394256614835390227244569;
    errBound = 24.1716394256614835390227244569 * std::max(eps/(1.-1.e-9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, reg, x = 10^-9 = 1.e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-8);
    refVal   = 20.346111037991381357292149948;
    errBound = 20.346111037991381357292149948 * std::max(eps/(1.-1.e-8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, reg, x = 10^-8 = 1.e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-7);
    refVal   = 16.8151334041954111521886314814;
    errBound = 16.8151334041954111521886314814 * std::max(eps/(1.-1.e-7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, reg, x = 10^-7 = 1.e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-6);
    refVal   = 13.5787120998985420754244804661;
    errBound = 13.5787120998985420754244804661 * std::max(eps/(1.-1.e-6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, reg, x = 10^-6 = 1.e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-5);
    refVal   = 10.6368858825057398138267022466;
    errBound = 10.6368858825057398138267022466 * std::max(eps/(1.-1.e-5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, reg, x = 10^-5 = 1.e-5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0001);
    refVal   = 7.9898939049823113326121467344;
    errBound = 7.9898939049823113326121467344 * std::max(eps/(1.-0.0001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, reg, x = 10^-4 = 0.0001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.001);
    refVal   = 5.63886450881613050117252085595;
    errBound = 5.63886450881613050117252085595 * std::max(eps/(1.-0.001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, reg, x = 10^-3 = 0.001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.01);
    refVal   = 3.58495886967125887823360267714;
    errBound = 3.58495886967125887823360267714 * std::max(eps/(1.-0.01),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, reg, x = 10^-2 = 0.01";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.05);
    refVal   = 2.30765501352664239775941261023;
    errBound = 2.30765501352664239775941261023 * std::max(eps/(1.-0.05),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, reg, x = 0.05 = 0.05";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.1);
    refVal   = 1.7757862713462034677321532362;
    errBound = 1.7757862713462034677321532362 * std::max(eps/(1.-0.1),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, reg, x = 0.1 = 0.1";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.15);
    refVal   = 1.45906048340681841206124818233;
    errBound = 1.45906048340681841206124818233 * std::max(eps/(1.-0.15),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, reg, x = 0.15 = 0.15";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.2);
    refVal   = 1.2274648222602481873279566035;
    errBound = 1.2274648222602481873279566035 * std::max(eps/(1.-0.2),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, reg, x = 0.2 = 0.2";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.25);
    refVal   = 1.04173025678408048107677905871;
    errBound = 1.04173025678408048107677905871 * std::max(eps/(1.-0.25),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, reg, x = 0.25 = 0.25";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.3);
    refVal   = 0.884778325897658307574829223266;
    errBound = 0.884778325897658307574829223266 * std::max(eps/(1.-0.3),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, reg, x = 0.3 = 0.3";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.35);
    refVal   = 0.747653502481008184131100025705;
    errBound = 0.747653502481008184131100025705 * std::max(eps/(1.-0.35),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, reg, x = 0.35 = 0.35";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.4);
    refVal   = 0.625077270137432649994133342331;
    errBound = 0.625077270137432649994133342331 * std::max(eps/(1.-0.4),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, reg, x = 0.4 = 0.4";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.45);
    refVal   = 0.513676531347856933945966319772;
    errBound = 0.513676531347856933945966319772 * std::max(eps/(1.-0.45),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, reg, x = 0.45 = 0.45";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.5);
    refVal   = 0.411163687571660613135747007949;
    errBound = 0.411163687571660613135747007949 * std::max(eps/(1.-0.5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, reg, x = 0.5 = 0.5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.55);
    refVal   = 0.315914702964704492917679814008;
    errBound = 0.315914702964704492917679814008 * std::max(eps/(1.-0.55),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, reg, x = 0.55 = 0.55";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.6);
    refVal   = 0.226733973242285389868882142924;
    errBound = 0.226733973242285389868882142924 * std::max(eps/(1.-0.6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, reg, x = 0.6 = 0.6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.65);
    refVal   = 0.142714880971783596565667340925;
    errBound = 0.142714880971783596565667340925 * std::max(eps/(1.-0.65),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, reg, x = 0.65 = 0.65";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.7);
    refVal   = 0.0631528494682723069688303561007;
    errBound = 0.0631528494682723069688303561007 * std::max(eps/(1.-0.7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, reg, x = 0.7 = 0.7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.75);
    refVal   = -0.0125111530596244361775344999229;
    errBound = 0.0125111530596244361775344999229 * std::max(eps/(1.-0.75),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, reg, x = 0.75 = 0.75";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.8);
    refVal   = -0.0847286161434186074671568335257;
    errBound = 0.0847286161434186074671568335257 * std::max(eps/(1.-0.8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, reg, x = 0.8 = 0.8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.85);
    refVal   = -0.153869814592880740962972004097;
    errBound = 0.153869814592880740962972004097 * std::max(eps/(1.-0.85),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, reg, x = 0.85 = 0.85";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9);
    refVal   = -0.220242517222626687701579239647;
    errBound = 0.220242517222626687701579239647 * std::max(eps/(1.-0.9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, reg, x = 0.9 = 0.9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.95);
    refVal   = -0.284105579325413216853924778963;
    errBound = 0.284105579325413216853924778963 * std::max(eps/(1.-0.95),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, reg, x = 0.95 = 0.95";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99);
    refVal   = -0.33353803688462414948726339534;
    errBound = 0.33353803688462414948726339534 * std::max(eps/(1.-0.99),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, reg, x = 1-10^-2 = 0.99";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999);
    refVal   = -0.344468715816189680597257001528;
    errBound = 0.344468715816189680597257001528 * std::max(eps/(1.-0.999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, reg, x = 1-10^-3 = 0.999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999);
    refVal   = -0.345558020493594118980140319567;
    errBound = 0.345558020493594118980140319567 * std::max(eps/(1.-0.9999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, reg, x = 1-10^-4 = 0.9999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999);
    refVal   = -0.345666913538271371911898134738;
    errBound = 0.345666913538271371911898134738 * std::max(eps/(1.-0.99999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, reg, x = 1-10^-5 = 0.99999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999);
    refVal   = -0.345677802468716049149691189815;
    errBound = 0.345677802468716049149691189815 * std::max(eps/(1.-0.999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, reg, x = 1-10^-6 = 0.999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999);
    refVal   = -0.345678891358020493826927469119;
    errBound = 0.345678891358020493826927469119 * std::max(eps/(1.-0.9999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, reg, x = 1-10^-7 = 0.9999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999999);
    refVal   = -0.345679000246913538271604705247;
    errBound = 0.345679000246913538271604705247 * std::max(eps/(1.-0.99999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, reg, x = 1-10^-8 = 0.99999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999);
    refVal   = -0.345679011135802468716049382483;
    errBound = 0.345679011135802468716049382483 * std::max(eps/(1.-0.999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, reg, x = 1-10^-9 = 0.999999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999999);
    refVal   = -0.34567901222469135802049382716;
    errBound = 0.34567901222469135802049382716 * std::max(eps/(1.-0.9999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, reg, x = 1-10^-10 = 0.9999999999";
  }
}

TEST(polAqqQNSEvenXspace,TruncAs2NormalDoublesPlusEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;
  auto eval = polAqqQNSEven_plus.truncate(2);
  {
    testVal  = eval(testAs,testLM,testNF,1.e-10);
    refVal   = 1.27160493839876543211148148148;
    errBound = 1.27160493839876543211148148148 * std::max(eps/(1.-1.e-10),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, plus, x = 10^-10 = 1.e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-9);
    refVal   = 1.27160493954320987781481481609;
    errBound = 1.27160493954320987781481481609 * std::max(eps/(1.-1.e-9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, plus, x = 10^-9 = 1.e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-8);
    refVal   = 1.27160495098765444814814941975;
    errBound = 1.27160495098765444814814941975 * std::max(eps/(1.-1.e-8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, plus, x = 10^-8 = 1.e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-7);
    refVal   = 1.27160506543211148148275308655;
    errBound = 1.27160506543211148148275308655 * std::max(eps/(1.-1.e-7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, plus, x = 10^-7 = 1.e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-6);
    refVal   = 1.27160620987781481608642102469;
    errBound = 1.27160620987781481608642102469 * std::max(eps/(1.-1.e-6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, plus, x = 10^-6 = 1.e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-5);
    refVal   = 1.2716176544481494197658025963;
    errBound = 1.2716176544481494197658025963 * std::max(eps/(1.-1.e-5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, plus, x = 10^-5 = 1.e-5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0001);
    refVal   = 1.2717321114827532135929642347;
    errBound = 1.2717321114827532135929642347 * std::max(eps/(1.-0.0001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, plus, x = 10^-4 = 0.0001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.001);
    refVal   = 1.27287781608769263090250744572;
    errBound = 1.27287781608769263090250744572 * std::max(eps/(1.-0.001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, plus, x = 10^-3 = 0.001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.01);
    refVal   = 1.28444943259758074572889387704;
    errBound = 1.28444943259758074572889387704 * std::max(eps/(1.-0.01),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, plus, x = 10^-2 = 0.01";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.05);
    refVal   = 1.33853151397011046133853151397;
    errBound = 1.33853151397011046133853151397 * std::max(eps/(1.-0.05),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, plus, x = 0.05 = 0.05";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.1);
    refVal   = 1.41289437585733882030178326475;
    errBound = 1.41289437585733882030178326475 * std::max(eps/(1.-0.1),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, plus, x = 0.1 = 0.1";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.15);
    refVal   = 1.49600580973129992737835875091;
    errBound = 1.49600580973129992737835875091 * std::max(eps/(1.-0.15),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, plus, x = 0.15 = 0.15";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.2);
    refVal   = 1.58950617283950617283950617284;
    errBound = 1.58950617283950617283950617284 * std::max(eps/(1.-0.2),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, plus, x = 0.2 = 0.2";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.25);
    refVal   = 1.6954732510288065843621399177;
    errBound = 1.6954732510288065843621399177 * std::max(eps/(1.-0.25),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, plus, x = 0.25 = 0.25";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.3);
    refVal   = 1.81657848324514991181657848325;
    errBound = 1.81657848324514991181657848325 * std::max(eps/(1.-0.3),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, plus, x = 0.3 = 0.3";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.35);
    refVal   = 1.95631528964862298195631528965;
    errBound = 1.95631528964862298195631528965 * std::max(eps/(1.-0.35),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, plus, x = 0.35 = 0.35";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.4);
    refVal   = 2.11934156378600823045267489712;
    errBound = 2.11934156378600823045267489712 * std::max(eps/(1.-0.4),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, plus, x = 0.4 = 0.4";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.45);
    refVal   = 2.31200897867564534231200897868;
    errBound = 2.31200897867564534231200897868 * std::max(eps/(1.-0.45),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, plus, x = 0.45 = 0.45";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.5);
    refVal   = 2.54320987654320987654320987654;
    errBound = 2.54320987654320987654320987654 * std::max(eps/(1.-0.5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, plus, x = 0.5 = 0.5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.55);
    refVal   = 2.82578875171467764060356652949;
    errBound = 2.82578875171467764060356652949 * std::max(eps/(1.-0.55),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, plus, x = 0.55 = 0.55";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.6);
    refVal   = 3.17901234567901234567901234568;
    errBound = 3.17901234567901234567901234568 * std::max(eps/(1.-0.6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, plus, x = 0.6 = 0.6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.65);
    refVal   = 3.63315696649029982363315696649;
    errBound = 3.63315696649029982363315696649 * std::max(eps/(1.-0.65),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, plus, x = 0.65 = 0.65";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.7);
    refVal   = 4.23868312757201646090534979424;
    errBound = 4.23868312757201646090534979424 * std::max(eps/(1.-0.7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, plus, x = 0.7 = 0.7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.75);
    refVal   = 5.08641975308641975308641975309;
    errBound = 5.08641975308641975308641975309 * std::max(eps/(1.-0.75),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, plus, x = 0.75 = 0.75";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.8);
    refVal   = 6.35802469135802469135802469136;
    errBound = 6.35802469135802469135802469136 * std::max(eps/(1.-0.8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, plus, x = 0.8 = 0.8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.85);
    refVal   = 8.47736625514403292181069958848;
    errBound = 8.47736625514403292181069958848 * std::max(eps/(1.-0.85),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, plus, x = 0.85 = 0.85";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9);
    refVal   = 12.7160493827160493827160493827;
    errBound = 12.7160493827160493827160493827 * std::max(eps/(1.-0.9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, plus, x = 0.9 = 0.9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.95);
    refVal   = 25.4320987654320987654320987654;
    errBound = 25.4320987654320987654320987654 * std::max(eps/(1.-0.95),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, plus, x = 0.95 = 0.95";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99);
    refVal   = 127.160493827160493827160493827;
    errBound = 127.160493827160493827160493827 * std::max(eps/(1.-0.99),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, plus, x = 1-10^-2 = 0.99";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999);
    refVal   = 1271.60493827160493827160493827;
    errBound = 1271.60493827160493827160493827 * std::max(eps/(1.-0.999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, plus, x = 1-10^-3 = 0.999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999);
    refVal   = 12716.0493827160493827160493827;
    errBound = 12716.0493827160493827160493827 * std::max(eps/(1.-0.9999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, plus, x = 1-10^-4 = 0.9999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999);
    refVal   = 127160.493827160493827160493827;
    errBound = 127160.493827160493827160493827 * std::max(eps/(1.-0.99999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, plus, x = 1-10^-5 = 0.99999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999);
    refVal   = 1.27160493827160493827160493827e6;
    errBound = 1.27160493827160493827160493827e6 * std::max(eps/(1.-0.999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, plus, x = 1-10^-6 = 0.999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999);
    refVal   = 1.27160493827160493827160493827e7;
    errBound = 1.27160493827160493827160493827e7 * std::max(eps/(1.-0.9999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, plus, x = 1-10^-7 = 0.9999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999999);
    refVal   = 1.27160493827160493827160493827e8;
    errBound = 1.27160493827160493827160493827e8 * std::max(eps/(1.-0.99999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, plus, x = 1-10^-8 = 0.99999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999);
    refVal   = 1.27160493827160493827160493827e9;
    errBound = 1.27160493827160493827160493827e9 * std::max(eps/(1.-0.999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, plus, x = 1-10^-9 = 0.999999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999999);
    refVal   = 1.27160493827160493827160493827e10;
    errBound = 1.27160493827160493827160493827e10 * std::max(eps/(1.-0.9999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSEven, TruncAs2NormalDoubles, plus, x = 1-10^-10 = 0.9999999999";
  }
}



TEST(polAqqQNSEvenXspace,TruncAs2NormalDoublesDeltaEvaluation)
{
  auto eval = polAqqQNSEven_delta.truncate(2);
  EXPECT_DOUBLE_EQ(eval(testAs,testLM,testNF),1.45677694886353911970334838331)
    << "^-- polAqqQNSEven, TruncAs2NormalDoubles, delta";
}



TEST(polAqqQNSEvenNspace,FullMoments)
{
  double refVal = 0., errBound = 0.;
  auto& rpd = polAqqQNSEven;
  integration_engine_gsl engine;
  auto mom = make_mellin_moment(rpd, engine, testAs, testLM, testNF);
  std::tuple<integration_status, double, double> res;

  {
    SCOPED_TRACE("polAqqQNSEven Full Mellin moment N=2");
    refVal = -13.1371375578672743446864514871;
    errBound = 13.1371375578672743446864514871 * 3.e-11;
    res = mom.integrate(2, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAqqQNSEven Full Mellin moment N=4");
    refVal = -28.0315757374001565266067553811;
    errBound = 28.0315757374001565266067553811 * 3.e-11;
    res = mom.integrate(4, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAqqQNSEven Full Mellin moment N=6");
    refVal = -36.060843394034104710099925415;
    errBound = 36.060843394034104710099925415 * 3.e-11;
    res = mom.integrate(6, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAqqQNSEven Full Mellin moment N=8");
    refVal = -41.661923747871398119379198654;
    errBound = 41.661923747871398119379198654 * 3.e-11;
    res = mom.integrate(8, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

}


TEST(polAqqQNSEvenNspace,TruncAs0Moments)
{
  double refVal = 0., errBound = 0.;
  auto rpd = polAqqQNSEven.truncate(0);
  integration_engine_gsl engine;
  auto mom = make_mellin_moment(rpd, engine, testAs, testLM, testNF);
  std::tuple<integration_status, double, double> res;

  {
    SCOPED_TRACE("polAqqQNSEven TruncAs0 Mellin moment N=2");
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
    SCOPED_TRACE("polAqqQNSEven TruncAs0 Mellin moment N=4");
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
    SCOPED_TRACE("polAqqQNSEven TruncAs0 Mellin moment N=6");
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
    SCOPED_TRACE("polAqqQNSEven TruncAs0 Mellin moment N=8");
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

}


TEST(polAqqQNSEvenNspace,TruncAs1Moments)
{
  double refVal = 0., errBound = 0.;
  auto rpd = polAqqQNSEven.truncate(1);
  integration_engine_gsl engine;
  auto mom = make_mellin_moment(rpd, engine, testAs, testLM, testNF);
  std::tuple<integration_status, double, double> res;

  {
    SCOPED_TRACE("polAqqQNSEven TruncAs1 Mellin moment N=2");
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
    SCOPED_TRACE("polAqqQNSEven TruncAs1 Mellin moment N=4");
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
    SCOPED_TRACE("polAqqQNSEven TruncAs1 Mellin moment N=6");
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
    SCOPED_TRACE("polAqqQNSEven TruncAs1 Mellin moment N=8");
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

}


TEST(polAqqQNSEvenNspace,TruncAs2Moments)
{
  double refVal = 0., errBound = 0.;
  auto rpd = polAqqQNSEven.truncate(2);
  integration_engine_gsl engine;
  auto mom = make_mellin_moment(rpd, engine, testAs, testLM, testNF);
  std::tuple<integration_status, double, double> res;

  {
    SCOPED_TRACE("polAqqQNSEven TruncAs2 Mellin moment N=2");
    refVal = 0.277777777777777777777777777778;
    errBound = 0.277777777777777777777777777778 * 3.e-11;
    res = mom.integrate(2, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAqqQNSEven TruncAs2 Mellin moment N=4");
    refVal = -0.89055632716049382716049382716;
    errBound = 0.89055632716049382716049382716 * 3.e-11;
    res = mom.integrate(4, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAqqQNSEven TruncAs2 Mellin moment N=6");
    refVal = -1.47234374287633924822133438913;
    errBound = 1.47234374287633924822133438913 * 3.e-11;
    res = mom.integrate(6, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAqqQNSEven TruncAs2 Mellin moment N=8");
    refVal = -1.86531274968569797766908710051;
    errBound = 1.86531274968569797766908710051 * 3.e-11;
    res = mom.integrate(8, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

}
