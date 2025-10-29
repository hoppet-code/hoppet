/******************************************************************************
 * This file is part of libome                                                *
 * Copyright (C) 2025 Arnd Behring, Kay Schoenwald                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 ******************************************************************************/

#include <limits>
#include <gtest/gtest.h>
#include <ome/integration_engine_gsl.h>
#include <ome/mellin.h>
#include <ome/polAqqQNSOdd.h>

using namespace ome;

const double eps = std::numeric_limits<double>::epsilon();
const double testAs = 0.25;
const double testLM = -5.;
const double testNF = 3.;



TEST(polAqqQNSOddXspace,FullPowersOf4RegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto& eval = polAqqQNSOdd_reg;

  {
    testVal  = eval(testAs,testLM,testNF,9.31322574615478515625e-10);
    refVal   = 4277.70101214639132947689711723;
    errBound = 4277.70101214639132947689711723 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, reg, x = 4^-15 = 9.31322574615478515625e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.7252902984619140625e-9);
    refVal   = 3410.43922779892380398777012758;
    errBound = 3410.43922779892380398777012758 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, reg, x = 4^-14 = 3.7252902984619140625e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.490116119384765625e-8);
    refVal   = 2679.96706947116633052388894681;
    errBound = 2679.96706947116633052388894681 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, reg, x = 4^-13 = 1.490116119384765625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,5.9604644775390625e-8);
    refVal   = 2071.44389170965514523000814093;
    errBound = 2071.44389170965514523000814093 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, reg, x = 4^-12 = 5.9604644775390625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,2.384185791015625e-7);
    refVal   = 1570.84982040592695899175172089;
    errBound = 1570.84982040592695899175172089 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, reg, x = 4^-11 = 2.384185791015625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,9.5367431640625e-7);
    refVal   = 1164.98579045707893735673763531;
    errBound = 1164.98579045707893735673763531 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, reg, x = 4^-10 = 9.5367431640625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.814697265625e-6);
    refVal   = 841.473630235022080583402478375;
    errBound = 841.473630235022080583402478375 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, reg, x = 4^-9 = 3.814697265625e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0000152587890625);
    refVal   = 588.756222541526337494710486275;
    errBound = 588.756222541526337494710486275 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, reg, x = 4^-8 = 0.0000152587890625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00006103515625);
    refVal   = 396.097709641153629387969604267;
    errBound = 396.097709641153629387969604267 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, reg, x = 4^-7 = 0.00006103515625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.000244140625);
    refVal   = 253.583344473775318771748595315;
    errBound = 253.583344473775318771748595315 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, reg, x = 4^-6 = 0.000244140625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0009765625);
    refVal   = 152.117169060496643428465976523;
    errBound = 152.117169060496643428465976523 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, reg, x = 4^-5 = 0.0009765625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00390625);
    refVal   = 83.4110578449723411075603983823;
    errBound = 83.4110578449723411075603983823 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, reg, x = 4^-4 = 0.00390625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.015625);
    refVal   = 39.943464366078806477615791985;
    errBound = 39.943464366078806477615791985 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, reg, x = 4^-3 = 0.015625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0625);
    refVal   = 14.7950009135025499124495707311;
    errBound = 14.7950009135025499124495707311 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, reg, x = 4^-2 = 0.0625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9375);
    refVal   = -19.5108624155838609072186356416;
    errBound = 19.5108624155838609072186356416 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, reg, x = 1-4^-2 = 0.9375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.984375);
    refVal   = -25.4673536173233760863614449979;
    errBound = 25.4673536173233760863614449979 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, reg, x = 1-4^-3 = 0.984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99609375);
    refVal   = -31.4096215293044683211221384916;
    errBound = 31.4096215293044683211221384916 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, reg, x = 1-4^-4 = 0.99609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9990234375);
    refVal   = -37.9876822443422242411982440131;
    errBound = 37.9876822443422242411982440131 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, reg, x = 1-4^-5 = 0.9990234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999755859375);
    refVal   = -45.5691650992202829111512419967;
    errBound = 45.5691650992202829111512419967 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, reg, x = 1-4^-6 = 0.999755859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99993896484375);
    refVal   = -54.4455818365323213187188189944;
    errBound = 54.4455818365323213187188189944 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, reg, x = 1-4^-7 = 0.99993896484375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999847412109375);
    refVal   = -64.887921112512144704580789027;
    errBound = 64.887921112512144704580789027 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, reg, x = 1-4^-8 = 0.9999847412109375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999996185302734375);
    refVal   = -77.1615406242022896735981046039;
    errBound = 77.1615406242022896735981046039 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, reg, x = 1-4^-9 = 0.999996185302734375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999904632568359375);
    refVal   = -91.5302156059937308772731748793;
    errBound = 91.5302156059937308772731748793 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, reg, x = 1-4^-10 = 0.99999904632568359375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999997615814208984375);
    refVal   = -108.257266252843829463705686112;
    errBound = 108.257266252843829463705686112 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, reg, x = 1-4^-11 = 0.9999997615814208984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999940395355224609375);
    refVal   = -127.605879503791459648879710192;
    errBound = 127.605879503791459648879710192 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, reg, x = 1-4^-12 = 0.999999940395355224609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999998509883880615234375);
    refVal   = -149.839202785349931163070897035;
    errBound = 149.839202785349931163070897035 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, reg, x = 1-4^-13 = 0.99999998509883880615234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999962747097015380859375);
    refVal   = -175.220371724981087106700919131;
    errBound = 175.220371724981087106700919131 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, reg, x = 1-4^-14 = 0.9999999962747097015380859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999068677425384521484375);
    refVal   = -204.012518417373330139328104554;
    errBound = 204.012518417373330139328104554 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, reg, x = 1-4^-15 = 0.999999999068677425384521484375";
  }
}

TEST(polAqqQNSOddXspace,FullPowersOf4PlusEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;
  auto& eval = polAqqQNSOdd_plus;
  {
    testVal  = eval(testAs,testLM,testNF,9.31322574615478515625e-10);
    refVal   = 19.8478639549625839237424922507;
    errBound = 19.8478639549625839237424922507 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, plus, x = 4^-15 = 9.31322574615478515625e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.7252902984619140625e-9);
    refVal   = 19.8478640104168754077863490609;
    errBound = 19.8478640104168754077863490609 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, plus, x = 4^-14 = 3.7252902984619140625e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.490116119384765625e-8);
    refVal   = 19.8478642322340444427118335806;
    errBound = 19.8478642322340444427118335806 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, plus, x = 4^-13 = 1.490116119384765625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,5.9604644775390625e-8);
    refVal   = 19.8478651195027701624175971435;
    errBound = 19.8478651195027701624175971435 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, plus, x = 4^-12 = 5.9604644775390625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,2.384185791015625e-7);
    refVal   = 19.8478686685784663214880366929;
    errBound = 19.8478686685784663214880366929 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, plus, x = 4^-11 = 2.384185791015625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,9.5367431640625e-7);
    refVal   = 19.847882864893943453643334231;
    errBound = 19.847882864893943453643334231 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, plus, x = 4^-10 = 9.5367431640625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.814697265625e-6);
    refVal   = 19.8479396503589326788280243906;
    errBound = 19.8479396503589326788280243906 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, plus, x = 4^-9 = 3.814697265625e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0000152587890625);
    refVal   = 19.8481667954682295308843167689;
    errBound = 19.8481667954682295308843167689 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, plus, x = 4^-8 = 0.0000152587890625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00006103515625);
    refVal   = 19.8490754278979799534045000832;
    errBound = 19.8490754278979799534045000832 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, plus, x = 4^-7 = 0.00006103515625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.000244140625);
    refVal   = 19.8527107896979612684142811272;
    errBound = 19.8527107896979612684142811272 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, plus, x = 4^-6 = 0.000244140625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0009765625);
    refVal   = 19.8672655630041914452972827996;
    errBound = 19.8672655630041914452972827996 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, plus, x = 4^-5 = 0.0009765625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00390625);
    refVal   = 19.9256986970130273024893336313;
    errBound = 19.9256986970130273024893336313 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, plus, x = 4^-4 = 0.00390625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.015625);
    refVal   = 20.1629093957869919132332542698;
    errBound = 20.1629093957869919132332542698 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, plus, x = 4^-3 = 0.015625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0625);
    refVal   = 21.1710548655763415088949169833;
    errBound = 21.1710548655763415088949169833 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, plus, x = 4^-2 = 0.0625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9375);
    refVal   = 317.565822983645122633423754749;
    errBound = 317.565822983645122633423754749 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, plus, x = 1-4^-2 = 0.9375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.984375);
    refVal   = 1270.263291934580490533695019;
    errBound = 1270.263291934580490533695019 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, plus, x = 1-4^-3 = 0.984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99609375);
    refVal   = 5081.05316773832196213478007599;
    errBound = 5081.05316773832196213478007599 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, plus, x = 1-4^-4 = 0.99609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9990234375);
    refVal   = 20324.212670953287848539120304;
    errBound = 20324.212670953287848539120304 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, plus, x = 1-4^-5 = 0.9990234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999755859375);
    refVal   = 81296.8506838131513941564812158;
    errBound = 81296.8506838131513941564812158 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, plus, x = 1-4^-6 = 0.999755859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99993896484375);
    refVal   = 325187.402735252605576625924863;
    errBound = 325187.402735252605576625924863 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, plus, x = 1-4^-7 = 0.99993896484375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999847412109375);
    refVal   = 1.30074961094101042230650369945e6;
    errBound = 1.30074961094101042230650369945e6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, plus, x = 1-4^-8 = 0.9999847412109375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999996185302734375);
    refVal   = 5.20299844376404168922601479781e6;
    errBound = 5.20299844376404168922601479781e6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, plus, x = 1-4^-9 = 0.999996185302734375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999904632568359375);
    refVal   = 2.08119937750561667569040591912e7;
    errBound = 2.08119937750561667569040591912e7 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, plus, x = 1-4^-10 = 0.99999904632568359375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999997615814208984375);
    refVal   = 8.3247975100224667027616236765e7;
    errBound = 8.3247975100224667027616236765e7 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, plus, x = 1-4^-11 = 0.9999997615814208984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999940395355224609375);
    refVal   = 3.3299190040089866811046494706e8;
    errBound = 3.3299190040089866811046494706e8 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, plus, x = 1-4^-12 = 0.999999940395355224609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999998509883880615234375);
    refVal   = 1.33196760160359467244185978824e9;
    errBound = 1.33196760160359467244185978824e9 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, plus, x = 1-4^-13 = 0.99999998509883880615234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999962747097015380859375);
    refVal   = 5.32787040641437868976743915296e9;
    errBound = 5.32787040641437868976743915296e9 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, plus, x = 1-4^-14 = 0.9999999962747097015380859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999068677425384521484375);
    refVal   = 2.13114816256575147590697566118e10;
    errBound = 2.13114816256575147590697566118e10 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullPowersOf4, plus, x = 1-4^-15 = 0.999999999068677425384521484375";
  }
}



TEST(polAqqQNSOddXspace,FullNormalDoublesRegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto& eval = polAqqQNSOdd_reg;

  {
    testVal  = eval(testAs,testLM,testNF,1.e-10);
    refVal   = 6000.87279535130638705135186729;
    errBound = 6000.87279535130638705135186729 * std::max(eps/(1.-1.e-10),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, reg, x = 10^-10 = 1.e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-9);
    refVal   = 4229.61568410987231584491526356;
    errBound = 4229.61568410987231584491526356 * std::max(eps/(1.-1.e-9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, reg, x = 10^-9 = 1.e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-8);
    refVal   = 2876.99993011825403606572779;
    errBound = 2876.99993011825403606572779 * std::max(eps/(1.-1.e-8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, reg, x = 10^-8 = 1.e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-7);
    refVal   = 1872.69670608742524211969716671;
    errBound = 1872.69670608742524211969716671 * std::max(eps/(1.-1.e-7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, reg, x = 10^-7 = 1.e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-6);
    refVal   = 1152.62405895190875574983638685;
    errBound = 1152.62405895190875574983638685 * std::max(eps/(1.-1.e-6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, reg, x = 10^-6 = 1.e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-5);
    refVal   = 658.947517407316462381489771447;
    errBound = 658.947517407316462381489771447 * std::max(eps/(1.-1.e-5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, reg, x = 10^-5 = 1.e-5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0001);
    refVal   = 340.081759907528114584944890211;
    errBound = 340.081759907528114584944890211 * std::max(eps/(1.-0.0001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, reg, x = 10^-4 = 0.0001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.001);
    refVal   = 150.688805802471530880671243372;
    errBound = 150.688805802471530880671243372 * std::max(eps/(1.-0.001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, reg, x = 10^-3 = 0.001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.01);
    refVal   = 51.6131688801579207335624982629;
    errBound = 51.6131688801579207335624982629 * std::max(eps/(1.-0.01),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, reg, x = 10^-2 = 0.01";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.05);
    refVal   = 17.886057646690013483387048357;
    errBound = 17.886057646690013483387048357 * std::max(eps/(1.-0.05),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, reg, x = 0.05 = 0.05";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.1);
    refVal   = 9.21066852414888580047439989416;
    errBound = 9.21066852414888580047439989416 * std::max(eps/(1.-0.1),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, reg, x = 0.1 = 0.1";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.15);
    refVal   = 5.19301917041247886458787738736;
    errBound = 5.19301917041247886458787738736 * std::max(eps/(1.-0.15),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, reg, x = 0.15 = 0.15";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.2);
    refVal   = 2.64090725323933430485693779708;
    errBound = 2.64090725323933430485693779708 * std::max(eps/(1.-0.2),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, reg, x = 0.2 = 0.2";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.25);
    refVal   = 0.74366858627870692023543875268;
    errBound = 0.74366858627870692023543875268 * std::max(eps/(1.-0.25),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, reg, x = 0.25 = 0.25";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.3);
    refVal   = -0.810080949028711172924740427451;
    errBound = 0.810080949028711172924740427451 * std::max(eps/(1.-0.3),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, reg, x = 0.3 = 0.3";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.35);
    refVal   = -2.16844262953297751297352004603;
    errBound = 2.16844262953297751297352004603 * std::max(eps/(1.-0.35),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, reg, x = 0.35 = 0.35";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.4);
    refVal   = -3.41247968960405945738296889479;
    errBound = 3.41247968960405945738296889479 * std::max(eps/(1.-0.4),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, reg, x = 0.4 = 0.4";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.45);
    refVal   = -4.5917226921782360501483253031;
    errBound = 4.5917226921782360501483253031 * std::max(eps/(1.-0.45),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, reg, x = 0.45 = 0.45";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.5);
    refVal   = -5.73969300929725450039665994846;
    errBound = 5.73969300929725450039665994846 * std::max(eps/(1.-0.5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, reg, x = 0.5 = 0.5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.55);
    refVal   = -6.88171406567811574525606552969;
    errBound = 6.88171406567811574525606552969 * std::max(eps/(1.-0.55),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, reg, x = 0.55 = 0.55";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.6);
    refVal   = -8.03950890693062758360103109962;
    errBound = 8.03950890693062758360103109962 * std::max(eps/(1.-0.6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, reg, x = 0.6 = 0.6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.65);
    refVal   = -9.23459860755594745524915754104;
    errBound = 9.23459860755594745524915754104 * std::max(eps/(1.-0.65),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, reg, x = 0.65 = 0.65";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.7);
    refVal   = -10.4917741627025124734072160184;
    errBound = 10.4917741627025124734072160184 * std::max(eps/(1.-0.7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, reg, x = 0.7 = 0.7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.75);
    refVal   = -11.8441022121010149284795270989;
    errBound = 11.8441022121010149284795270989 * std::max(eps/(1.-0.75),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, reg, x = 0.75 = 0.75";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.8);
    refVal   = -13.342441880580632405994285264;
    errBound = 13.342441880580632405994285264 * std::max(eps/(1.-0.8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, reg, x = 0.8 = 0.8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.85);
    refVal   = -15.0781426231572696352219676761;
    errBound = 15.0781426231572696352219676761 * std::max(eps/(1.-0.85),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, reg, x = 0.85 = 0.85";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9);
    refVal   = -17.2532323389491693477300844487;
    errBound = 17.2532323389491693477300844487 * std::max(eps/(1.-0.9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, reg, x = 0.9 = 0.9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.95);
    refVal   = -20.517724058116631500217617181;
    errBound = 20.517724058116631500217617181 * std::max(eps/(1.-0.95),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, reg, x = 0.95 = 0.95";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99);
    refVal   = -27.3406539710709365601867318856;
    errBound = 27.3406539710709365601867318856 * std::max(eps/(1.-0.99),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, reg, x = 1-10^-2 = 0.99";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999);
    refVal   = -37.8675967979137602224618389367;
    errBound = 37.8675967979137602224618389367 * std::max(eps/(1.-0.999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, reg, x = 1-10^-3 = 0.999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999);
    refVal   = -51.1186164689109434377416379702;
    errBound = 51.1186164689109434377416379702 * std::max(eps/(1.-0.9999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, reg, x = 1-10^-4 = 0.9999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999);
    refVal   = -68.4229693296745688338680768043;
    errBound = 68.4229693296745688338680768043 * std::max(eps/(1.-0.99999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, reg, x = 1-10^-5 = 0.99999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999);
    refVal   = -91.0011141608012117282328224133;
    errBound = 91.0011141608012117282328224133 * std::max(eps/(1.-0.999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, reg, x = 1-10^-6 = 0.999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999);
    refVal   = -120.060649227164793888275051521;
    errBound = 120.060649227164793888275051521 * std::max(eps/(1.-0.9999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, reg, x = 1-10^-7 = 0.9999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999999);
    refVal   = -156.807555763426650029501940281;
    errBound = 156.807555763426650029501940281 * std::max(eps/(1.-0.99999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, reg, x = 1-10^-8 = 0.99999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999);
    refVal   = -202.447602254377575998262546993;
    errBound = 202.447602254377575998262546993 * std::max(eps/(1.-0.999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, reg, x = 1-10^-9 = 0.999999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999999);
    refVal   = -258.186528602680525890886833807;
    errBound = 258.186528602680525890886833807 * std::max(eps/(1.-0.9999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, reg, x = 1-10^-10 = 0.9999999999";
  }
}

TEST(polAqqQNSOddXspace,FullNormalDoublesPlusEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;
  auto& eval = polAqqQNSOdd_plus;
  {
    testVal  = eval(testAs,testLM,testNF,1.e-10);
    refVal   = 19.8478639384626065584352453277;
    errBound = 19.8478639384626065584352453277 * std::max(eps/(1.-1.e-10),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, plus, x = 10^-10 = 1.e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-9);
    refVal   = 19.8478639563256841209146687927;
    errBound = 19.8478639563256841209146687927 * std::max(eps/(1.-1.e-9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, plus, x = 10^-9 = 1.e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-8);
    refVal   = 19.8478641349564615141535998134;
    errBound = 19.8478641349564615141535998134 * std::max(eps/(1.-1.e-8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, plus, x = 10^-8 = 1.e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-7);
    refVal   = 19.8478659212644122910302137749;
    errBound = 19.8478659212644122910302137749 * std::max(eps/(1.-1.e-7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, plus, x = 10^-7 = 1.e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-6);
    refVal   = 19.8478837843616045261935108653;
    errBound = 19.8478837843616045261935108653 * std::max(eps/(1.-1.e-6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, plus, x = 10^-6 = 1.e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-5);
    refVal   = 19.8480624171019911845008296801;
    errBound = 19.8480624171019911845008296801 * std::max(eps/(1.-1.e-5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, plus, x = 10^-5 = 1.e-5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0001);
    refVal   = 19.8498489213699571603050151733;
    errBound = 19.8498489213699571603050151733 * std::max(eps/(1.-0.0001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, plus, x = 10^-4 = 0.0001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.001);
    refVal   = 19.8677316681459661307197043762;
    errBound = 19.8677316681459661307197043762 * std::max(eps/(1.-0.001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, plus, x = 10^-3 = 0.001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.01);
    refVal   = 20.0483474105836567319080653251;
    errBound = 20.0483474105836567319080653251 * std::max(eps/(1.-0.01),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, plus, x = 10^-2 = 0.01";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.05);
    refVal   = 20.8924883541871791206199838651;
    errBound = 20.8924883541871791206199838651 * std::max(eps/(1.-0.05),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, plus, x = 0.05 = 0.05";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.1);
    refVal   = 22.0531821516420224050988718576;
    errBound = 22.0531821516420224050988718576 * std::max(eps/(1.-0.1),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, plus, x = 0.1 = 0.1";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.15);
    refVal   = 23.3504281605621413701046878492;
    errBound = 23.3504281605621413701046878492 * std::max(eps/(1.-0.15),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, plus, x = 0.15 = 0.15";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.2);
    refVal   = 24.8098299205972752057362308398;
    errBound = 24.8098299205972752057362308398 * std::max(eps/(1.-0.2),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, plus, x = 0.2 = 0.2";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.25);
    refVal   = 26.4638185819704268861186462291;
    errBound = 26.4638185819704268861186462291 * std::max(eps/(1.-0.25),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, plus, x = 0.25 = 0.25";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.3);
    refVal   = 28.3540913378254573779842638169;
    errBound = 28.3540913378254573779842638169 * std::max(eps/(1.-0.3),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, plus, x = 0.3 = 0.3";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.35);
    refVal   = 30.5351752868889540993676687259;
    errBound = 30.5351752868889540993676687259 * std::max(eps/(1.-0.35),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, plus, x = 0.35 = 0.35";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.4);
    refVal   = 33.0797732274630336076483077864;
    errBound = 33.0797732274630336076483077864 * std::max(eps/(1.-0.4),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, plus, x = 0.4 = 0.4";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.45);
    refVal   = 36.0870253390505821174345175851;
    errBound = 36.0870253390505821174345175851 * std::max(eps/(1.-0.45),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, plus, x = 0.45 = 0.45";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.5);
    refVal   = 39.6957278729556403291779693437;
    errBound = 39.6957278729556403291779693437 * std::max(eps/(1.-0.5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, plus, x = 0.5 = 0.5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.55);
    refVal   = 44.1063643032840448101977437152;
    errBound = 44.1063643032840448101977437152 * std::max(eps/(1.-0.55),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, plus, x = 0.55 = 0.55";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.6);
    refVal   = 49.6196598411945504114724616796;
    errBound = 49.6196598411945504114724616796 * std::max(eps/(1.-0.6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, plus, x = 0.6 = 0.6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.65);
    refVal   = 56.7081826756509147559685276338;
    errBound = 56.7081826756509147559685276338 * std::max(eps/(1.-0.65),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, plus, x = 0.65 = 0.65";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.7);
    refVal   = 66.1595464549260672152966155728;
    errBound = 66.1595464549260672152966155728 * std::max(eps/(1.-0.7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, plus, x = 0.7 = 0.7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.75);
    refVal   = 79.3914557459112806583559386873;
    errBound = 79.3914557459112806583559386873 * std::max(eps/(1.-0.75),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, plus, x = 0.75 = 0.75";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.8);
    refVal   = 99.2393196823891008229449233592;
    errBound = 99.2393196823891008229449233592 * std::max(eps/(1.-0.8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, plus, x = 0.8 = 0.8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.85);
    refVal   = 132.319092909852134430593231146;
    errBound = 132.319092909852134430593231146 * std::max(eps/(1.-0.85),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, plus, x = 0.85 = 0.85";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9);
    refVal   = 198.478639364778201645889846718;
    errBound = 198.478639364778201645889846718 * std::max(eps/(1.-0.9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, plus, x = 0.9 = 0.9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.95);
    refVal   = 396.957278729556403291779693437;
    errBound = 396.957278729556403291779693437 * std::max(eps/(1.-0.95),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, plus, x = 0.95 = 0.95";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99);
    refVal   = 1984.78639364778201645889846718;
    errBound = 1984.78639364778201645889846718 * std::max(eps/(1.-0.99),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, plus, x = 1-10^-2 = 0.99";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999);
    refVal   = 19847.8639364778201645889846718;
    errBound = 19847.8639364778201645889846718 * std::max(eps/(1.-0.999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, plus, x = 1-10^-3 = 0.999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999);
    refVal   = 198478.639364778201645889846718;
    errBound = 198478.639364778201645889846718 * std::max(eps/(1.-0.9999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, plus, x = 1-10^-4 = 0.9999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999);
    refVal   = 1.98478639364778201645889846718e6;
    errBound = 1.98478639364778201645889846718e6 * std::max(eps/(1.-0.99999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, plus, x = 1-10^-5 = 0.99999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999);
    refVal   = 1.98478639364778201645889846718e7;
    errBound = 1.98478639364778201645889846718e7 * std::max(eps/(1.-0.999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, plus, x = 1-10^-6 = 0.999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999);
    refVal   = 1.98478639364778201645889846718e8;
    errBound = 1.98478639364778201645889846718e8 * std::max(eps/(1.-0.9999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, plus, x = 1-10^-7 = 0.9999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999999);
    refVal   = 1.98478639364778201645889846718e9;
    errBound = 1.98478639364778201645889846718e9 * std::max(eps/(1.-0.99999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, plus, x = 1-10^-8 = 0.99999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999);
    refVal   = 1.98478639364778201645889846718e10;
    errBound = 1.98478639364778201645889846718e10 * std::max(eps/(1.-0.999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, plus, x = 1-10^-9 = 0.999999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999999);
    refVal   = 1.98478639364778201645889846718e11;
    errBound = 1.98478639364778201645889846718e11 * std::max(eps/(1.-0.9999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, FullNormalDoubles, plus, x = 1-10^-10 = 0.9999999999";
  }
}



TEST(polAqqQNSOddXspace,FullDeltaEvaluation)
{
  auto& eval = polAqqQNSOdd_delta;
  EXPECT_DOUBLE_EQ(eval(testAs,testLM,testNF),12.1501748418929015300593526677)
    << "^-- polAqqQNSOdd, Full, delta";
}







TEST(polAqqQNSOddXspace,TruncAs0NormalDoublesDeltaEvaluation)
{
  auto eval = polAqqQNSOdd_delta.truncate(0);
  EXPECT_DOUBLE_EQ(eval(testAs,testLM,testNF),1.)
    << "^-- polAqqQNSOdd, TruncAs0NormalDoubles, delta";
}







TEST(polAqqQNSOddXspace,TruncAs1NormalDoublesDeltaEvaluation)
{
  auto eval = polAqqQNSOdd_delta.truncate(1);
  EXPECT_DOUBLE_EQ(eval(testAs,testLM,testNF),1.)
    << "^-- polAqqQNSOdd, TruncAs1NormalDoubles, delta";
}



TEST(polAqqQNSOddXspace,TruncAs2PowersOf4RegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto eval = polAqqQNSOdd_reg.truncate(2);

  {
    testVal  = eval(testAs,testLM,testNF,9.31322574615478515625e-10);
    refVal   = 24.2945391531051607896541861718;
    errBound = 24.2945391531051607896541861718 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, reg, x = 4^-15 = 9.31322574615478515625e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.7252902984619140625e-9);
    refVal   = 21.950577165068803555966157488;
    errBound = 21.950577165068803555966157488 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, reg, x = 4^-14 = 3.7252902984619140625e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.490116119384765625e-8);
    refVal   = 19.7133826179923339322847932777;
    errBound = 19.7133826179923339322847932777 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, reg, x = 4^-13 = 1.490116119384765625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,5.9604644775390625e-8);
    refVal   = 17.5829557667636721009561321108;
    errBound = 17.5829557667636721009561321108 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, reg, x = 4^-12 = 5.9604644775390625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,2.384185791015625e-7);
    refVal   = 15.5592974691149663221180610747;
    errBound = 15.5592974691149663221180610747 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, reg, x = 4^-11 = 2.384185791015625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,9.5367431640625e-7);
    refVal   = 13.6424105516611732880307771909;
    errBound = 13.6424105516611732880307771909 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, reg, x = 4^-10 = 9.5367431640625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.814697265625e-6);
    refVal   = 11.8323040755482182127343650263;
    errBound = 11.8323040755482182127343650263 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, reg, x = 4^-9 = 3.814697265625e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0000152587890625);
    refVal   = 10.1290059929594676451689047194;
    errBound = 10.1290059929594676451689047194 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, reg, x = 4^-8 = 0.0000152587890625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00006103515625);
    refVal   = 8.5325977125286692136080901449;
    errBound = 8.5325977125286692136080901449 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, reg, x = 4^-7 = 0.00006103515625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.000244140625);
    refVal   = 7.04329452172840778557083330018;
    errBound = 7.04329452172840778557083330018 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, reg, x = 4^-6 = 0.000244140625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0009765625);
    refVal   = 5.66156327512726792864800589954;
    errBound = 5.66156327512726792864800589954 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, reg, x = 4^-5 = 0.0009765625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00390625);
    refVal   = 4.38790550536101365575516619009;
    errBound = 4.38790550536101365575516619009 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, reg, x = 4^-4 = 0.00390625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.015625);
    refVal   = 3.21993361190786973785471883912;
    errBound = 3.21993361190786973785471883912 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, reg, x = 4^-3 = 0.015625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0625);
    refVal   = 2.13659810997390165973939871603;
    errBound = 2.13659810997390165973939871603 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, reg, x = 4^-2 = 0.0625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9375);
    refVal   = -0.268362475445021794896393223338;
    errBound = 0.268362475445021794896393223338 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, reg, x = 1-4^-2 = 0.9375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.984375);
    refVal   = -0.3266713135058627817810942937;
    errBound = 0.3266713135058627817810942937 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, reg, x = 1-4^-3 = 0.984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99609375);
    refVal   = -0.340946513246320684792526634059;
    errBound = 0.340946513246320684792526634059 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, reg, x = 1-4^-4 = 0.99609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9990234375);
    refVal   = -0.34449709175903860580951890324;
    errBound = 0.34449709175903860580951890324 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, reg, x = 1-4^-5 = 0.9990234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999755859375);
    refVal   = -0.345383607307621750884392650659;
    errBound = 0.345383607307621750884392650659 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, reg, x = 1-4^-6 = 0.999755859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99993896484375);
    refVal   = -0.345605165778065901050331736072;
    errBound = 0.345605165778065901050331736072 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, reg, x = 1-4^-7 = 0.99993896484375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999847412109375);
    refVal   = -0.345660550996982297006902759782;
    errBound = 0.345660550996982297006902759782 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, reg, x = 1-4^-8 = 0.9999847412109375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999996185302734375);
    refVal   = -0.345674397026829661536539041158;
    errBound = 0.345674397026829661536539041158 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, reg, x = 1-4^-9 = 0.999996185302734375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999904632568359375);
    refVal   = -0.345677858517111967299608103272;
    errBound = 0.345677858517111967299608103272 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, reg, x = 1-4^-10 = 0.99999904632568359375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999997615814208984375);
    refVal   = -0.345678723888608831733060846091;
    errBound = 0.345678723888608831733060846091 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, reg, x = 1-4^-11 = 0.9999997615814208984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999940395355224609375);
    refVal   = -0.345678940231415940980860170939;
    errBound = 0.345678940231415940980860170939 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, reg, x = 1-4^-12 = 0.999999940395355224609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999998509883880615234375);
    refVal   = -0.345678994317113524116210587615;
    errBound = 0.345678994317113524116210587615 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, reg, x = 1-4^-13 = 0.99999998509883880615234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999962747097015380859375);
    refVal   = -0.345679007838537657764044881895;
    errBound = 0.345679007838537657764044881895 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, reg, x = 1-4^-14 = 0.9999999962747097015380859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999068677425384521484375);
    refVal   = -0.345679011218893674792503782246;
    errBound = 0.345679011218893674792503782246 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, reg, x = 1-4^-15 = 0.999999999068677425384521484375";
  }
}

TEST(polAqqQNSOddXspace,TruncAs2PowersOf4PlusEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;
  auto eval = polAqqQNSOdd_plus.truncate(2);
  {
    testVal  = eval(testAs,testLM,testNF,9.31322574615478515625e-10);
    refVal   = 1.27160493945587932437941415068;
    errBound = 1.27160493945587932437941415068 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, plus, x = 4^-15 = 9.31322574615478515625e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.7252902984619140625e-9);
    refVal   = 1.27160494300870249593813948108;
    errBound = 1.27160494300870249593813948108 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, plus, x = 4^-14 = 3.7252902984619140625e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.490116119384765625e-8);
    refVal   = 1.27160495721999538070250915835;
    errBound = 1.27160495721999538070250915835 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, plus, x = 4^-13 = 1.490116119384765625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,5.9604644775390625e-8);
    refVal   = 1.27160501406517009623166793273;
    errBound = 1.27160501406517009623166793273 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, plus, x = 4^-12 = 5.9604644775390625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,2.384185791015625e-7);
    refVal   = 1.27160524144591978190711202291;
    errBound = 1.27160524144591978190711202291 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, plus, x = 4^-11 = 2.384185791015625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,9.5367431640625e-7);
    refVal   = 1.2716061509697317023132216768;
    errBound = 1.2716061509697317023132216768 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, plus, x = 4^-10 = 9.5367431640625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.814697265625e-6);
    refVal   = 1.27160978907799027606410091034;
    errBound = 1.27160978907799027606410091034 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, plus, x = 4^-9 = 3.814697265625e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0000152587890625);
    refVal   = 1.27162434171920197199310141504;
    errBound = 1.27162434171920197199310141504 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, plus, x = 4^-8 = 0.0000152587890625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00006103515625);
    refVal   = 1.27168255561508730444008883041;
    errBound = 1.27168255561508730444008883041 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, plus, x = 4^-7 = 0.00006103515625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.000244140625);
    refVal   = 1.27191546450805710064969324229;
    errBound = 1.27191546450805710064969324229 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, plus, x = 4^-6 = 0.000244140625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0009765625);
    refVal   = 1.27284795385153808092876193235;
    errBound = 1.27284795385153808092876193235 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, plus, x = 4^-5 = 0.0009765625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00390625);
    refVal   = 1.27659162430404260469619946744;
    errBound = 1.27659162430404260469619946744 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, plus, x = 4^-4 = 0.00390625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.015625);
    refVal   = 1.2917891436409954928473446992;
    errBound = 1.2917891436409954928473446992 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, plus, x = 4^-3 = 0.015625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0625);
    refVal   = 1.35637860082304526748971193416;
    errBound = 1.35637860082304526748971193416 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, plus, x = 4^-2 = 0.0625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9375);
    refVal   = 20.3456790123456790123456790123;
    errBound = 20.3456790123456790123456790123 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, plus, x = 1-4^-2 = 0.9375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.984375);
    refVal   = 81.3827160493827160493827160494;
    errBound = 81.3827160493827160493827160494 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, plus, x = 1-4^-3 = 0.984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99609375);
    refVal   = 325.530864197530864197530864198;
    errBound = 325.530864197530864197530864198 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, plus, x = 1-4^-4 = 0.99609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9990234375);
    refVal   = 1302.12345679012345679012345679;
    errBound = 1302.12345679012345679012345679 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, plus, x = 1-4^-5 = 0.9990234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999755859375);
    refVal   = 5208.49382716049382716049382716;
    errBound = 5208.49382716049382716049382716 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, plus, x = 1-4^-6 = 0.999755859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99993896484375);
    refVal   = 20833.9753086419753086419753086;
    errBound = 20833.9753086419753086419753086 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, plus, x = 1-4^-7 = 0.99993896484375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999847412109375);
    refVal   = 83335.9012345679012345679012346;
    errBound = 83335.9012345679012345679012346 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, plus, x = 1-4^-8 = 0.9999847412109375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999996185302734375);
    refVal   = 333343.604938271604938271604938;
    errBound = 333343.604938271604938271604938 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, plus, x = 1-4^-9 = 0.999996185302734375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999904632568359375);
    refVal   = 1.33337441975308641975308641975e6;
    errBound = 1.33337441975308641975308641975e6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, plus, x = 1-4^-10 = 0.99999904632568359375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999997615814208984375);
    refVal   = 5.33349767901234567901234567901e6;
    errBound = 5.33349767901234567901234567901e6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, plus, x = 1-4^-11 = 0.9999997615814208984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999940395355224609375);
    refVal   = 2.1333990716049382716049382716e7;
    errBound = 2.1333990716049382716049382716e7 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, plus, x = 1-4^-12 = 0.999999940395355224609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999998509883880615234375);
    refVal   = 8.53359628641975308641975308642e7;
    errBound = 8.53359628641975308641975308642e7 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, plus, x = 1-4^-13 = 0.99999998509883880615234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999962747097015380859375);
    refVal   = 3.41343851456790123456790123457e8;
    errBound = 3.41343851456790123456790123457e8 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, plus, x = 1-4^-14 = 0.9999999962747097015380859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999068677425384521484375);
    refVal   = 1.36537540582716049382716049383e9;
    errBound = 1.36537540582716049382716049383e9 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2PowersOf4, plus, x = 1-4^-15 = 0.999999999068677425384521484375";
  }
}



TEST(polAqqQNSOddXspace,TruncAs2NormalDoublesRegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto eval = polAqqQNSOdd_reg.truncate(2);

  {
    testVal  = eval(testAs,testLM,testNF,1.e-10);
    refVal   = 28.2917178181689775076785574172;
    errBound = 28.2917178181689775076785574172 * std::max(eps/(1.-1.e-10),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, reg, x = 10^-10 = 1.e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-9);
    refVal   = 24.1716394256614835390227244569;
    errBound = 24.1716394256614835390227244569 * std::max(eps/(1.-1.e-9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, reg, x = 10^-9 = 1.e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-8);
    refVal   = 20.346111037991381357292149948;
    errBound = 20.346111037991381357292149948 * std::max(eps/(1.-1.e-8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, reg, x = 10^-8 = 1.e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-7);
    refVal   = 16.8151334041954111521886314814;
    errBound = 16.8151334041954111521886314814 * std::max(eps/(1.-1.e-7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, reg, x = 10^-7 = 1.e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-6);
    refVal   = 13.5787120998985420754244804661;
    errBound = 13.5787120998985420754244804661 * std::max(eps/(1.-1.e-6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, reg, x = 10^-6 = 1.e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-5);
    refVal   = 10.6368858825057398138267022466;
    errBound = 10.6368858825057398138267022466 * std::max(eps/(1.-1.e-5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, reg, x = 10^-5 = 1.e-5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0001);
    refVal   = 7.9898939049823113326121467344;
    errBound = 7.9898939049823113326121467344 * std::max(eps/(1.-0.0001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, reg, x = 10^-4 = 0.0001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.001);
    refVal   = 5.63886450881613050117252085595;
    errBound = 5.63886450881613050117252085595 * std::max(eps/(1.-0.001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, reg, x = 10^-3 = 0.001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.01);
    refVal   = 3.58495886967125887823360267714;
    errBound = 3.58495886967125887823360267714 * std::max(eps/(1.-0.01),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, reg, x = 10^-2 = 0.01";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.05);
    refVal   = 2.30765501352664239775941261023;
    errBound = 2.30765501352664239775941261023 * std::max(eps/(1.-0.05),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, reg, x = 0.05 = 0.05";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.1);
    refVal   = 1.7757862713462034677321532362;
    errBound = 1.7757862713462034677321532362 * std::max(eps/(1.-0.1),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, reg, x = 0.1 = 0.1";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.15);
    refVal   = 1.45906048340681841206124818233;
    errBound = 1.45906048340681841206124818233 * std::max(eps/(1.-0.15),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, reg, x = 0.15 = 0.15";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.2);
    refVal   = 1.2274648222602481873279566035;
    errBound = 1.2274648222602481873279566035 * std::max(eps/(1.-0.2),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, reg, x = 0.2 = 0.2";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.25);
    refVal   = 1.04173025678408048107677905871;
    errBound = 1.04173025678408048107677905871 * std::max(eps/(1.-0.25),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, reg, x = 0.25 = 0.25";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.3);
    refVal   = 0.884778325897658307574829223266;
    errBound = 0.884778325897658307574829223266 * std::max(eps/(1.-0.3),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, reg, x = 0.3 = 0.3";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.35);
    refVal   = 0.747653502481008184131100025705;
    errBound = 0.747653502481008184131100025705 * std::max(eps/(1.-0.35),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, reg, x = 0.35 = 0.35";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.4);
    refVal   = 0.625077270137432649994133342331;
    errBound = 0.625077270137432649994133342331 * std::max(eps/(1.-0.4),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, reg, x = 0.4 = 0.4";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.45);
    refVal   = 0.513676531347856933945966319772;
    errBound = 0.513676531347856933945966319772 * std::max(eps/(1.-0.45),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, reg, x = 0.45 = 0.45";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.5);
    refVal   = 0.411163687571660613135747007949;
    errBound = 0.411163687571660613135747007949 * std::max(eps/(1.-0.5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, reg, x = 0.5 = 0.5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.55);
    refVal   = 0.315914702964704492917679814008;
    errBound = 0.315914702964704492917679814008 * std::max(eps/(1.-0.55),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, reg, x = 0.55 = 0.55";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.6);
    refVal   = 0.226733973242285389868882142924;
    errBound = 0.226733973242285389868882142924 * std::max(eps/(1.-0.6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, reg, x = 0.6 = 0.6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.65);
    refVal   = 0.142714880971783596565667340925;
    errBound = 0.142714880971783596565667340925 * std::max(eps/(1.-0.65),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, reg, x = 0.65 = 0.65";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.7);
    refVal   = 0.0631528494682723069688303561007;
    errBound = 0.0631528494682723069688303561007 * std::max(eps/(1.-0.7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, reg, x = 0.7 = 0.7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.75);
    refVal   = -0.0125111530596244361775344999229;
    errBound = 0.0125111530596244361775344999229 * std::max(eps/(1.-0.75),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, reg, x = 0.75 = 0.75";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.8);
    refVal   = -0.0847286161434186074671568335257;
    errBound = 0.0847286161434186074671568335257 * std::max(eps/(1.-0.8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, reg, x = 0.8 = 0.8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.85);
    refVal   = -0.153869814592880740962972004097;
    errBound = 0.153869814592880740962972004097 * std::max(eps/(1.-0.85),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, reg, x = 0.85 = 0.85";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9);
    refVal   = -0.220242517222626687701579239647;
    errBound = 0.220242517222626687701579239647 * std::max(eps/(1.-0.9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, reg, x = 0.9 = 0.9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.95);
    refVal   = -0.284105579325413216853924778963;
    errBound = 0.284105579325413216853924778963 * std::max(eps/(1.-0.95),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, reg, x = 0.95 = 0.95";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99);
    refVal   = -0.33353803688462414948726339534;
    errBound = 0.33353803688462414948726339534 * std::max(eps/(1.-0.99),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, reg, x = 1-10^-2 = 0.99";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999);
    refVal   = -0.344468715816189680597257001528;
    errBound = 0.344468715816189680597257001528 * std::max(eps/(1.-0.999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, reg, x = 1-10^-3 = 0.999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999);
    refVal   = -0.345558020493594118980140319567;
    errBound = 0.345558020493594118980140319567 * std::max(eps/(1.-0.9999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, reg, x = 1-10^-4 = 0.9999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999);
    refVal   = -0.345666913538271371911898134738;
    errBound = 0.345666913538271371911898134738 * std::max(eps/(1.-0.99999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, reg, x = 1-10^-5 = 0.99999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999);
    refVal   = -0.345677802468716049149691189815;
    errBound = 0.345677802468716049149691189815 * std::max(eps/(1.-0.999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, reg, x = 1-10^-6 = 0.999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999);
    refVal   = -0.345678891358020493826927469119;
    errBound = 0.345678891358020493826927469119 * std::max(eps/(1.-0.9999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, reg, x = 1-10^-7 = 0.9999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999999);
    refVal   = -0.345679000246913538271604705247;
    errBound = 0.345679000246913538271604705247 * std::max(eps/(1.-0.99999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, reg, x = 1-10^-8 = 0.99999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999);
    refVal   = -0.345679011135802468716049382483;
    errBound = 0.345679011135802468716049382483 * std::max(eps/(1.-0.999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, reg, x = 1-10^-9 = 0.999999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999999);
    refVal   = -0.34567901222469135802049382716;
    errBound = 0.34567901222469135802049382716 * std::max(eps/(1.-0.9999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, reg, x = 1-10^-10 = 0.9999999999";
  }
}

TEST(polAqqQNSOddXspace,TruncAs2NormalDoublesPlusEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;
  auto eval = polAqqQNSOdd_plus.truncate(2);
  {
    testVal  = eval(testAs,testLM,testNF,1.e-10);
    refVal   = 1.27160493839876543211148148148;
    errBound = 1.27160493839876543211148148148 * std::max(eps/(1.-1.e-10),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, plus, x = 10^-10 = 1.e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-9);
    refVal   = 1.27160493954320987781481481609;
    errBound = 1.27160493954320987781481481609 * std::max(eps/(1.-1.e-9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, plus, x = 10^-9 = 1.e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-8);
    refVal   = 1.27160495098765444814814941975;
    errBound = 1.27160495098765444814814941975 * std::max(eps/(1.-1.e-8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, plus, x = 10^-8 = 1.e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-7);
    refVal   = 1.27160506543211148148275308655;
    errBound = 1.27160506543211148148275308655 * std::max(eps/(1.-1.e-7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, plus, x = 10^-7 = 1.e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-6);
    refVal   = 1.27160620987781481608642102469;
    errBound = 1.27160620987781481608642102469 * std::max(eps/(1.-1.e-6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, plus, x = 10^-6 = 1.e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-5);
    refVal   = 1.2716176544481494197658025963;
    errBound = 1.2716176544481494197658025963 * std::max(eps/(1.-1.e-5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, plus, x = 10^-5 = 1.e-5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0001);
    refVal   = 1.2717321114827532135929642347;
    errBound = 1.2717321114827532135929642347 * std::max(eps/(1.-0.0001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, plus, x = 10^-4 = 0.0001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.001);
    refVal   = 1.27287781608769263090250744572;
    errBound = 1.27287781608769263090250744572 * std::max(eps/(1.-0.001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, plus, x = 10^-3 = 0.001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.01);
    refVal   = 1.28444943259758074572889387704;
    errBound = 1.28444943259758074572889387704 * std::max(eps/(1.-0.01),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, plus, x = 10^-2 = 0.01";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.05);
    refVal   = 1.33853151397011046133853151397;
    errBound = 1.33853151397011046133853151397 * std::max(eps/(1.-0.05),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, plus, x = 0.05 = 0.05";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.1);
    refVal   = 1.41289437585733882030178326475;
    errBound = 1.41289437585733882030178326475 * std::max(eps/(1.-0.1),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, plus, x = 0.1 = 0.1";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.15);
    refVal   = 1.49600580973129992737835875091;
    errBound = 1.49600580973129992737835875091 * std::max(eps/(1.-0.15),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, plus, x = 0.15 = 0.15";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.2);
    refVal   = 1.58950617283950617283950617284;
    errBound = 1.58950617283950617283950617284 * std::max(eps/(1.-0.2),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, plus, x = 0.2 = 0.2";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.25);
    refVal   = 1.6954732510288065843621399177;
    errBound = 1.6954732510288065843621399177 * std::max(eps/(1.-0.25),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, plus, x = 0.25 = 0.25";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.3);
    refVal   = 1.81657848324514991181657848325;
    errBound = 1.81657848324514991181657848325 * std::max(eps/(1.-0.3),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, plus, x = 0.3 = 0.3";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.35);
    refVal   = 1.95631528964862298195631528965;
    errBound = 1.95631528964862298195631528965 * std::max(eps/(1.-0.35),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, plus, x = 0.35 = 0.35";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.4);
    refVal   = 2.11934156378600823045267489712;
    errBound = 2.11934156378600823045267489712 * std::max(eps/(1.-0.4),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, plus, x = 0.4 = 0.4";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.45);
    refVal   = 2.31200897867564534231200897868;
    errBound = 2.31200897867564534231200897868 * std::max(eps/(1.-0.45),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, plus, x = 0.45 = 0.45";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.5);
    refVal   = 2.54320987654320987654320987654;
    errBound = 2.54320987654320987654320987654 * std::max(eps/(1.-0.5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, plus, x = 0.5 = 0.5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.55);
    refVal   = 2.82578875171467764060356652949;
    errBound = 2.82578875171467764060356652949 * std::max(eps/(1.-0.55),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, plus, x = 0.55 = 0.55";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.6);
    refVal   = 3.17901234567901234567901234568;
    errBound = 3.17901234567901234567901234568 * std::max(eps/(1.-0.6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, plus, x = 0.6 = 0.6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.65);
    refVal   = 3.63315696649029982363315696649;
    errBound = 3.63315696649029982363315696649 * std::max(eps/(1.-0.65),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, plus, x = 0.65 = 0.65";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.7);
    refVal   = 4.23868312757201646090534979424;
    errBound = 4.23868312757201646090534979424 * std::max(eps/(1.-0.7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, plus, x = 0.7 = 0.7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.75);
    refVal   = 5.08641975308641975308641975309;
    errBound = 5.08641975308641975308641975309 * std::max(eps/(1.-0.75),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, plus, x = 0.75 = 0.75";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.8);
    refVal   = 6.35802469135802469135802469136;
    errBound = 6.35802469135802469135802469136 * std::max(eps/(1.-0.8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, plus, x = 0.8 = 0.8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.85);
    refVal   = 8.47736625514403292181069958848;
    errBound = 8.47736625514403292181069958848 * std::max(eps/(1.-0.85),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, plus, x = 0.85 = 0.85";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9);
    refVal   = 12.7160493827160493827160493827;
    errBound = 12.7160493827160493827160493827 * std::max(eps/(1.-0.9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, plus, x = 0.9 = 0.9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.95);
    refVal   = 25.4320987654320987654320987654;
    errBound = 25.4320987654320987654320987654 * std::max(eps/(1.-0.95),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, plus, x = 0.95 = 0.95";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99);
    refVal   = 127.160493827160493827160493827;
    errBound = 127.160493827160493827160493827 * std::max(eps/(1.-0.99),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, plus, x = 1-10^-2 = 0.99";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999);
    refVal   = 1271.60493827160493827160493827;
    errBound = 1271.60493827160493827160493827 * std::max(eps/(1.-0.999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, plus, x = 1-10^-3 = 0.999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999);
    refVal   = 12716.0493827160493827160493827;
    errBound = 12716.0493827160493827160493827 * std::max(eps/(1.-0.9999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, plus, x = 1-10^-4 = 0.9999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999);
    refVal   = 127160.493827160493827160493827;
    errBound = 127160.493827160493827160493827 * std::max(eps/(1.-0.99999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, plus, x = 1-10^-5 = 0.99999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999);
    refVal   = 1.27160493827160493827160493827e6;
    errBound = 1.27160493827160493827160493827e6 * std::max(eps/(1.-0.999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, plus, x = 1-10^-6 = 0.999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999);
    refVal   = 1.27160493827160493827160493827e7;
    errBound = 1.27160493827160493827160493827e7 * std::max(eps/(1.-0.9999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, plus, x = 1-10^-7 = 0.9999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999999);
    refVal   = 1.27160493827160493827160493827e8;
    errBound = 1.27160493827160493827160493827e8 * std::max(eps/(1.-0.99999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, plus, x = 1-10^-8 = 0.99999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999);
    refVal   = 1.27160493827160493827160493827e9;
    errBound = 1.27160493827160493827160493827e9 * std::max(eps/(1.-0.999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, plus, x = 1-10^-9 = 0.999999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999999);
    refVal   = 1.27160493827160493827160493827e10;
    errBound = 1.27160493827160493827160493827e10 * std::max(eps/(1.-0.9999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, plus, x = 1-10^-10 = 0.9999999999";
  }
}



TEST(polAqqQNSOddXspace,TruncAs2NormalDoublesDeltaEvaluation)
{
  auto eval = polAqqQNSOdd_delta.truncate(2);
  EXPECT_DOUBLE_EQ(eval(testAs,testLM,testNF),1.45677694886353911970334838331)
    << "^-- polAqqQNSOdd, TruncAs2NormalDoubles, delta";
}



TEST(polAqqQNSOddNspace,FullMoments)
{
  double refVal = 0., errBound = 0.;
  auto& rpd = polAqqQNSOdd;
  integration_engine_gsl engine;
  auto mom = make_mellin_moment(rpd, engine, testAs, testLM, testNF);
  std::tuple<integration_status, double, double> res;

  {
    SCOPED_TRACE("polAqqQNSOdd Full Mellin moment N=3");
    refVal = -22.0931532003857707128861290183;
    errBound = 22.0931532003857707128861290183 * 3.e-11;
    res = mom.integrate(3, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAqqQNSOdd Full Mellin moment N=5");
    refVal = -32.4675519992296219908912960781;
    errBound = 32.4675519992296219908912960781 * 3.e-11;
    res = mom.integrate(5, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAqqQNSOdd Full Mellin moment N=7");
    refVal = -39.0616711789071974294761166602;
    errBound = 39.0616711789071974294761166602 * 3.e-11;
    res = mom.integrate(7, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAqqQNSOdd Full Mellin moment N=9");
    refVal = -43.9417691634403222575059391201;
    errBound = 43.9417691634403222575059391201 * 3.e-11;
    res = mom.integrate(9, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

}


TEST(polAqqQNSOddNspace,TruncAs0Moments)
{
  double refVal = 0., errBound = 0.;
  auto rpd = polAqqQNSOdd.truncate(0);
  integration_engine_gsl engine;
  auto mom = make_mellin_moment(rpd, engine, testAs, testLM, testNF);
  std::tuple<integration_status, double, double> res;

  {
    SCOPED_TRACE("polAqqQNSOdd TruncAs0 Mellin moment N=3");
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
    SCOPED_TRACE("polAqqQNSOdd TruncAs0 Mellin moment N=5");
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
    SCOPED_TRACE("polAqqQNSOdd TruncAs0 Mellin moment N=7");
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
    SCOPED_TRACE("polAqqQNSOdd TruncAs0 Mellin moment N=9");
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


TEST(polAqqQNSOddNspace,TruncAs1Moments)
{
  double refVal = 0., errBound = 0.;
  auto rpd = polAqqQNSOdd.truncate(1);
  integration_engine_gsl engine;
  auto mom = make_mellin_moment(rpd, engine, testAs, testLM, testNF);
  std::tuple<integration_status, double, double> res;

  {
    SCOPED_TRACE("polAqqQNSOdd TruncAs1 Mellin moment N=3");
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
    SCOPED_TRACE("polAqqQNSOdd TruncAs1 Mellin moment N=5");
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
    SCOPED_TRACE("polAqqQNSOdd TruncAs1 Mellin moment N=7");
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
    SCOPED_TRACE("polAqqQNSOdd TruncAs1 Mellin moment N=9");
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


TEST(polAqqQNSOddNspace,TruncAs2Moments)
{
  double refVal = 0., errBound = 0.;
  auto rpd = polAqqQNSOdd.truncate(2);
  integration_engine_gsl engine;
  auto mom = make_mellin_moment(rpd, engine, testAs, testLM, testNF);
  std::tuple<integration_status, double, double> res;

  {
    SCOPED_TRACE("polAqqQNSOdd TruncAs2 Mellin moment N=3");
    refVal = -0.443383487654320987654320987654;
    errBound = 0.443383487654320987654320987654 * 3.e-11;
    res = mom.integrate(3, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAqqQNSOdd TruncAs2 Mellin moment N=5");
    refVal = -1.21589454732510288065843621399;
    errBound = 1.21589454732510288065843621399 * 3.e-11;
    res = mom.integrate(5, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAqqQNSOdd TruncAs2 Mellin moment N=7");
    refVal = -1.6843752520463652833267345739;
    errBound = 1.6843752520463652833267345739 * 3.e-11;
    res = mom.integrate(7, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAqqQNSOdd TruncAs2 Mellin moment N=9");
    refVal = -2.02322586911299564570749587966;
    errBound = 2.02322586911299564570749587966 * 3.e-11;
    res = mom.integrate(9, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

}
