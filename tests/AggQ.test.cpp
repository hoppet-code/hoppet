/******************************************************************************
 * This file is part of libome                                                *
 * Copyright (C) 2025 Arnd Behring, Kay Schoenwald                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 ******************************************************************************/

#include <limits>
#include <gtest/gtest.h>
#include <ome/integration_engine_gsl.h>
#include <ome/mellin.h>
#include <ome/AggQ.h>

using namespace ome;

const double eps = std::numeric_limits<double>::epsilon();
const double testAs = 0.25;
const double testLM = -5.;
const double testNF = 3.;



TEST(AggQXspace,FullPowersOf4RegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto& eval = AggQ_reg;

  {
    testVal  = eval(testAs,testLM,testNF,9.31322574615478515625e-10);
    refVal   = -3.18073222247051136702149781484e11;
    errBound = 3.18073222247051136702149781484e11 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, reg, x = 4^-15 = 9.31322574615478515625e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.7252902984619140625e-9);
    refVal   = -7.28774773316686520235330689017e10;
    errBound = 7.28774773316686520235330689017e10 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, reg, x = 4^-14 = 3.7252902984619140625e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.490116119384765625e-8);
    refVal   = -1.65590839834642301071592755218e10;
    errBound = 1.65590839834642301071592755218e10 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, reg, x = 4^-13 = 1.490116119384765625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,5.9604644775390625e-8);
    refVal   = -3.72464512427336701231177416472e9;
    errBound = 3.72464512427336701231177416472e9 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, reg, x = 4^-12 = 5.9604644775390625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,2.384185791015625e-7);
    refVal   = -8.27342998630197543690500731704e8;
    errBound = 8.27342998630197543690500731704e8 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, reg, x = 4^-11 = 2.384185791015625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,9.5367431640625e-7);
    refVal   = -1.80857249862532688257357832367e8;
    errBound = 1.80857249862532688257357832367e8 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, reg, x = 4^-10 = 9.5367431640625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.814697265625e-6);
    refVal   = -3.87048471364316009644684986669e7;
    errBound = 3.87048471364316009644684986669e7 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, reg, x = 4^-9 = 3.814697265625e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0000152587890625);
    refVal   = -8.04018075060200237915086119559e6;
    errBound = 8.04018075060200237915086119559e6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, reg, x = 4^-8 = 0.0000152587890625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00006103515625);
    refVal   = -1.59637586508102286660519155994e6;
    errBound = 1.59637586508102286660519155994e6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, reg, x = 4^-7 = 0.00006103515625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.000244140625);
    refVal   = -293454.895760814593724889740805;
    errBound = 293454.895760814593724889740805 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, reg, x = 4^-6 = 0.000244140625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0009765625);
    refVal   = -46097.9916694516288517306469225;
    errBound = 46097.9916694516288517306469225 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, reg, x = 4^-5 = 0.0009765625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00390625);
    refVal   = -4524.37759974140833844966541952;
    errBound = 4524.37759974140833844966541952 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, reg, x = 4^-4 = 0.00390625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.015625);
    refVal   = 540.560307910847184226974044704;
    errBound = 540.560307910847184226974044704 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, reg, x = 4^-3 = 0.015625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0625);
    refVal   = 423.813432837403576587263934077;
    errBound = 423.813432837403576587263934077 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, reg, x = 4^-2 = 0.0625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9375);
    refVal   = -79.3013156037107394090333314378;
    errBound = 79.3013156037107394090333314378 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, reg, x = 1-4^-2 = 0.9375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.984375);
    refVal   = -105.1855483551343758401036625932;
    errBound = 105.1855483551343758401036625932 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, reg, x = 1-4^-3 = 0.984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99609375);
    refVal   = -131.3291069367718749436262685271;
    errBound = 131.3291069367718749436262685271 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, reg, x = 1-4^-4 = 0.99609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9990234375);
    refVal   = -160.968772913091399763894043099;
    errBound = 160.968772913091399763894043099 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, reg, x = 1-4^-5 = 0.9990234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999755859375);
    refVal   = -195.950179975834752327419319079;
    errBound = 195.950179975834752327419319079 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, reg, x = 1-4^-6 = 0.999755859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99993896484375);
    refVal   = -237.724983921826997618146573082;
    errBound = 237.724983921826997618146573082 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, reg, x = 1-4^-7 = 0.99993896484375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999847412109375);
    refVal   = -287.647857411094997914197458036;
    errBound = 287.647857411094997914197458036 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, reg, x = 1-4^-8 = 0.9999847412109375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999996185302734375);
    refVal   = -347.053352874740774676294259052;
    errBound = 347.053352874740774676294259052 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, reg, x = 1-4^-9 = 0.999996185302734375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999904632568359375);
    refVal   = -417.273031909237309820805822467;
    errBound = 417.273031909237309820805822467 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, reg, x = 1-4^-10 = 0.99999904632568359375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999997615814208984375);
    refVal   = -499.638483087096496516487507857;
    errBound = 499.638483087096496516487507857 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, reg, x = 1-4^-11 = 0.9999997615814208984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999940395355224609375);
    refVal   = -595.481562539162476960480517572;
    errBound = 595.481562539162476960480517572 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, reg, x = 1-4^-12 = 0.999999940395355224609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999998509883880615234375);
    refVal   = -706.134274776820076751392063885;
    errBound = 706.134274776820076751392063885 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, reg, x = 1-4^-13 = 0.99999998509883880615234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999962747097015380859375);
    refVal   = -832.928685528836316015019009849;
    errBound = 832.928685528836316015019009849 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, reg, x = 1-4^-14 = 0.9999999962747097015380859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999068677425384521484375);
    refVal   = -977.196882664829020869825316776;
    errBound = 977.196882664829020869825316776 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, reg, x = 1-4^-15 = 0.999999999068677425384521484375";
  }
}

TEST(AggQXspace,FullPowersOf4PlusEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;
  auto& eval = AggQ_plus;
  {
    testVal  = eval(testAs,testLM,testNF,9.31322574615478515625e-10);
    refVal   = 42.2734346371860400952092060316;
    errBound = 42.2734346371860400952092060316 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, plus, x = 4^-15 = 9.31322574615478515625e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.7252902984619140625e-9);
    refVal   = 42.2734347552966524876352738601;
    errBound = 42.2734347552966524876352738601 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, plus, x = 4^-14 = 3.7252902984619140625e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.490116119384765625e-8);
    refVal   = 42.2734352277391086572844208845;
    errBound = 42.2734352277391086572844208845 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, plus, x = 4^-13 = 1.490116119384765625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,5.9604644775390625e-8);
    refVal   = 42.273437117509038935005216199;
    errBound = 42.273437117509038935005216199 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, plus, x = 4^-12 = 5.9604644775390625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,2.384185791015625e-7);
    refVal   = 42.273444676590449632272247516;
    errBound = 42.273444676590449632272247516 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, plus, x = 4^-11 = 2.384185791015625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,9.5367431640625e-7);
    refVal   = 42.2734749129431258288602113757;
    errBound = 42.2734749129431258288602113757 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, plus, x = 4^-10 = 9.5367431640625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.814697265625e-6);
    refVal   = 42.2735958587863667597428656719;
    errBound = 42.2735958587863667597428656719 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, plus, x = 4^-9 = 3.814697265625e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0000152587890625);
    refVal   = 42.2740796490800127470026475769;
    errBound = 42.2740796490800127470026475769 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, plus, x = 4^-8 = 0.0000152587890625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00006103515625);
    refVal   = 42.2760149209921661993349586302;
    errBound = 42.2760149209921661993349586302 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, plus, x = 4^-7 = 0.00006103515625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.000244140625);
    refVal   = 42.2837577808678057902139577069;
    errBound = 42.2837577808678057902139577069 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, plus, x = 4^-6 = 0.000244140625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0009765625);
    refVal   = 42.3147576032877968501774576759;
    errBound = 42.3147576032877968501774576759 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, plus, x = 4^-5 = 0.0009765625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00390625);
    refVal   = 42.439212772709231546795626669;
    errBound = 42.439212772709231546795626669 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, plus, x = 4^-4 = 0.00390625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.015625);
    refVal   = 42.944441496193865255686050796;
    errBound = 42.944441496193865255686050796 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, plus, x = 4^-3 = 0.015625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0625);
    refVal   = 45.0916635710035585184703533359;
    errBound = 45.0916635710035585184703533359 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, plus, x = 4^-2 = 0.0625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9375);
    refVal   = 676.374953565053377777055300038;
    errBound = 676.374953565053377777055300038 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, plus, x = 1-4^-2 = 0.9375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.984375);
    refVal   = 2705.49981426021351110822120015;
    errBound = 2705.49981426021351110822120015 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, plus, x = 1-4^-3 = 0.984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99609375);
    refVal   = 10821.9992570408540444328848006;
    errBound = 10821.9992570408540444328848006 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, plus, x = 1-4^-4 = 0.99609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9990234375);
    refVal   = 43287.9970281634161777315392024;
    errBound = 43287.9970281634161777315392024 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, plus, x = 1-4^-5 = 0.9990234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999755859375);
    refVal   = 173151.98811265366471092615681;
    errBound = 173151.98811265366471092615681 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, plus, x = 1-4^-6 = 0.999755859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99993896484375);
    refVal   = 692607.952450614658843704627239;
    errBound = 692607.952450614658843704627239 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, plus, x = 1-4^-7 = 0.99993896484375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999847412109375);
    refVal   = 2.77043180980245863537481850895e6;
    errBound = 2.77043180980245863537481850895e6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, plus, x = 1-4^-8 = 0.9999847412109375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999996185302734375);
    refVal   = 1.10817272392098345414992740358e7;
    errBound = 1.10817272392098345414992740358e7 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, plus, x = 1-4^-9 = 0.999996185302734375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999904632568359375);
    refVal   = 4.43269089568393381659970961433e7;
    errBound = 4.43269089568393381659970961433e7 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, plus, x = 1-4^-10 = 0.99999904632568359375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999997615814208984375);
    refVal   = 1.77307635827357352663988384573e8;
    errBound = 1.77307635827357352663988384573e8 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, plus, x = 1-4^-11 = 0.9999997615814208984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999940395355224609375);
    refVal   = 7.09230543309429410655953538292e8;
    errBound = 7.09230543309429410655953538292e8 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, plus, x = 1-4^-12 = 0.999999940395355224609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999998509883880615234375);
    refVal   = 2.83692217323771764262381415317e9;
    errBound = 2.83692217323771764262381415317e9 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, plus, x = 1-4^-13 = 0.99999998509883880615234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999962747097015380859375);
    refVal   = 1.13476886929508705704952566127e10;
    errBound = 1.13476886929508705704952566127e10 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, plus, x = 1-4^-14 = 0.9999999962747097015380859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999068677425384521484375);
    refVal   = 4.53907547718034822819810264507e10;
    errBound = 4.53907547718034822819810264507e10 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullPowersOf4, plus, x = 1-4^-15 = 0.999999999068677425384521484375";
  }
}



TEST(AggQXspace,FullNormalDoublesRegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto& eval = AggQ_reg;

  {
    testVal  = eval(testAs,testLM,testNF,1.e-10);
    refVal   = -3.36048956073637537615974763149e12;
    errBound = 3.36048956073637537615974763149e12 * std::max(eps/(1.-1.e-10),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, reg, x = 10^-10 = 1.e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-9);
    refVal   = -2.94959091267904158774139095788e11;
    errBound = 2.94959091267904158774139095788e11 * std::max(eps/(1.-1.e-9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, reg, x = 10^-9 = 1.e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-8);
    refVal   = -2.53867837854739144883424677841e10;
    errBound = 2.53867837854739144883424677841e10 * std::max(eps/(1.-1.e-8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, reg, x = 10^-8 = 1.e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-7);
    refVal   = -2.12768860961130313164381086901e9;
    errBound = 2.12768860961130313164381086901e9 * std::max(eps/(1.-1.e-7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, reg, x = 10^-7 = 1.e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-6);
    refVal   = -1.71630497775704030388882814193e8;
    errBound = 1.71630497775704030388882814193e8 * std::max(eps/(1.-1.e-6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, reg, x = 10^-6 = 1.e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-5);
    refVal   = -1.30314115150677855170342811574e7;
    errBound = 1.30314115150677855170342811574e7 * std::max(eps/(1.-1.e-5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, reg, x = 10^-5 = 1.e-5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0001);
    refVal   = -883283.759159587044308416900759;
    errBound = 883283.759159587044308416900759 * std::max(eps/(1.-0.0001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, reg, x = 10^-4 = 0.0001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.001);
    refVal   = -44554.0540324198346610872228876;
    errBound = 44554.0540324198346610872228876 * std::max(eps/(1.-0.001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, reg, x = 10^-3 = 0.001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.01);
    refVal   = 44.443447590686173568091974633;
    errBound = 44.443447590686173568091974633 * std::max(eps/(1.-0.01),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, reg, x = 10^-2 = 0.01";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.05);
    refVal   = 493.21562345876071400084535835;
    errBound = 493.21562345876071400084535835 * std::max(eps/(1.-0.05),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, reg, x = 0.05 = 0.05";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.1);
    refVal   = 282.287364107970722480820996055;
    errBound = 282.287364107970722480820996055 * std::max(eps/(1.-0.1),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, reg, x = 0.1 = 0.1";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.15);
    refVal   = 178.500736758262040260054085448;
    errBound = 178.500736758262040260054085448 * std::max(eps/(1.-0.15),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, reg, x = 0.15 = 0.15";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.2);
    refVal   = 118.0383316938327621759361310328;
    errBound = 118.0383316938327621759361310328 * std::max(eps/(1.-0.2),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, reg, x = 0.2 = 0.2";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.25);
    refVal   = 78.6465545745293819484226970249;
    errBound = 78.6465545745293819484226970249 * std::max(eps/(1.-0.25),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, reg, x = 0.25 = 0.25";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.3);
    refVal   = 51.0181330372177376970355036842;
    errBound = 51.0181330372177376970355036842 * std::max(eps/(1.-0.3),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, reg, x = 0.3 = 0.3";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.35);
    refVal   = 30.5854679968607024994670134455;
    errBound = 30.5854679968607024994670134455 * std::max(eps/(1.-0.35),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, reg, x = 0.35 = 0.35";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.4);
    refVal   = 14.8357562466563434485790865178;
    errBound = 14.8357562466563434485790865178 * std::max(eps/(1.-0.4),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, reg, x = 0.4 = 0.4";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.45);
    refVal   = 2.26213181633170537559483686196;
    errBound = 2.26213181633170537559483686196 * std::max(eps/(1.-0.45),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, reg, x = 0.45 = 0.45";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.5);
    refVal   = -8.10459201840268733707015752686;
    errBound = 8.10459201840268733707015752686 * std::max(eps/(1.-0.5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, reg, x = 0.5 = 0.5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.55);
    refVal   = -16.9277813465069506353078175011;
    errBound = 16.9277813465069506353078175011 * std::max(eps/(1.-0.55),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, reg, x = 0.55 = 0.55";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.6);
    refVal   = -24.6904389864237480581518532248;
    errBound = 24.6904389864237480581518532248 * std::max(eps/(1.-0.6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, reg, x = 0.6 = 0.6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.65);
    refVal   = -31.7697444055215097875401664358;
    errBound = 31.7697444055215097875401664358 * std::max(eps/(1.-0.65),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, reg, x = 0.65 = 0.65";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.7);
    refVal   = -38.4877501302640045648575855194;
    errBound = 38.4877501302640045648575855194 * std::max(eps/(1.-0.7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, reg, x = 0.7 = 0.7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.75);
    refVal   = -45.1551998538676991625825721868;
    errBound = 45.1551998538676991625825721868 * std::max(eps/(1.-0.75),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, reg, x = 0.75 = 0.75";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.8);
    refVal   = -52.1268461480542650462961117914;
    errBound = 52.1268461480542650462961117914 * std::max(eps/(1.-0.8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, reg, x = 0.8 = 0.8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.85);
    refVal   = -59.9085201112166743305392421668;
    errBound = 59.9085201112166743305392421668 * std::max(eps/(1.-0.85),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, reg, x = 0.85 = 0.85";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9);
    refVal   = -69.4649911752863262667979987405;
    errBound = 69.4649911752863262667979987405 * std::max(eps/(1.-0.9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, reg, x = 0.9 = 0.9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.95);
    refVal   = -83.6768563104991901809846302716;
    errBound = 83.6768563104991901809846302716 * std::max(eps/(1.-0.95),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, reg, x = 0.95 = 0.95";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99);
    refVal   = -113.3704830067264965284452997548;
    errBound = 113.3704830067264965284452997548 * std::max(eps/(1.-0.99),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, reg, x = 1-10^-2 = 0.99";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999);
    refVal   = -160.421213670264586027003845498;
    errBound = 160.421213670264586027003845498 * std::max(eps/(1.-0.999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, reg, x = 1-10^-3 = 0.999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999);
    refVal   = -221.982392128265874081141736211;
    errBound = 221.982392128265874081141736211 * std::max(eps/(1.-0.9999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, reg, x = 1-10^-4 = 0.9999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999);
    refVal   = -304.689776057899076631068507725;
    errBound = 304.689776057899076631068507725 * std::max(eps/(1.-0.99999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, reg, x = 1-10^-5 = 0.99999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999);
    refVal   = -414.677322360108083274457563806;
    errBound = 414.677322360108083274457563806 * std::max(eps/(1.-0.999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, reg, x = 1-10^-6 = 0.999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999);
    refVal   = -558.047048610332181519848843261;
    errBound = 558.047048610332181519848843261 * std::max(eps/(1.-0.9999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, reg, x = 1-10^-7 = 0.9999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999999);
    refVal   = -740.902016712486217932433829952;
    errBound = 740.902016712486217932433829952 * std::max(eps/(1.-0.99999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, reg, x = 1-10^-8 = 0.99999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999);
    refVal   = -969.346063816559281955543807006;
    errBound = 969.346063816559281955543807006 * std::max(eps/(1.-0.999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, reg, x = 1-10^-9 = 0.999999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999999);
    refVal   = -1249.483193986379054743603622781;
    errBound = 1249.483193986379054743603622781 * std::max(eps/(1.-0.9999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, reg, x = 1-10^-10 = 0.9999999999";
  }
}

TEST(AggQXspace,FullNormalDoublesPlusEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;
  auto& eval = AggQ_plus;
  {
    testVal  = eval(testAs,testLM,testNF,1.e-10);
    refVal   = 42.2734346020431795712702742095;
    errBound = 42.2734346020431795712702742095 * std::max(eps/(1.-1.e-10),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, plus, x = 10^-10 = 1.e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-9);
    refVal   = 42.2734346400892707511552270035;
    errBound = 42.2734346400892707511552270035 * std::max(eps/(1.-1.e-9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, plus, x = 10^-9 = 1.e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-8);
    refVal   = 42.273435020550186316567819418;
    errBound = 42.273435020550186316567819418 * std::max(eps/(1.-1.e-8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, plus, x = 10^-8 = 1.e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-7);
    refVal   = 42.2734388251597186270378189561;
    errBound = 42.2734388251597186270378189561 * std::max(eps/(1.-1.e-7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, plus, x = 10^-7 = 1.e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-6);
    refVal   = 42.2734768712927074037733600257;
    errBound = 42.2734768712927074037733600257 * std::max(eps/(1.-1.e-6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, plus, x = 10^-6 = 1.e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-5);
    refVal   = 42.2738573363892000030659869122;
    errBound = 42.2738573363892000030659869122 * std::max(eps/(1.-1.e-5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, plus, x = 10^-5 = 1.e-5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0001);
    refVal   = 42.2776623640522413351994762;
    errBound = 42.2776623640522413351994762 * std::max(eps/(1.-0.0001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, plus, x = 10^-4 = 0.0001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.001);
    refVal   = 42.3157503481640001111771333857;
    errBound = 42.3157503481640001111771333857 * std::max(eps/(1.-0.001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, plus, x = 10^-3 = 0.001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.01);
    refVal   = 42.700438987692763748551470962;
    errBound = 42.700438987692763748551470962 * std::max(eps/(1.-0.01),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, plus, x = 10^-2 = 0.01";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.05);
    refVal   = 44.4983522082271959063852171077;
    errBound = 44.4983522082271959063852171077 * std::max(eps/(1.-0.05),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, plus, x = 0.05 = 0.05";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.1);
    refVal   = 46.9704828864620401234066180582;
    errBound = 46.9704828864620401234066180582 * std::max(eps/(1.-0.1),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, plus, x = 0.1 = 0.1";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.15);
    refVal   = 49.7334524680186307189011250028;
    errBound = 49.7334524680186307189011250028 * std::max(eps/(1.-0.15),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, plus, x = 0.15 = 0.15";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.2);
    refVal   = 52.8417932472697951388324453155;
    errBound = 52.8417932472697951388324453155 * std::max(eps/(1.-0.2),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, plus, x = 0.2 = 0.2";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.25);
    refVal   = 56.3645794637544481480879416698;
    errBound = 56.3645794637544481480879416698 * std::max(eps/(1.-0.25),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, plus, x = 0.25 = 0.25";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.3);
    refVal   = 60.3906208540226230158085089319;
    errBound = 60.3906208540226230158085089319 * std::max(eps/(1.-0.3),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, plus, x = 0.3 = 0.3";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.35);
    refVal   = 65.0360532274089786324091634652;
    errBound = 65.0360532274089786324091634652 * std::max(eps/(1.-0.35),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, plus, x = 0.35 = 0.35";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.4);
    refVal   = 70.4557243296930601851099270873;
    errBound = 70.4557243296930601851099270873 * std::max(eps/(1.-0.4),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, plus, x = 0.4 = 0.4";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.45);
    refVal   = 76.8607901778469747473926477316;
    errBound = 76.8607901778469747473926477316 * std::max(eps/(1.-0.45),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, plus, x = 0.45 = 0.45";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.5);
    refVal   = 84.5468691956316722221319125047;
    errBound = 84.5468691956316722221319125047 * std::max(eps/(1.-0.5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, plus, x = 0.5 = 0.5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.55);
    refVal   = 93.9409657729240802468132361164;
    errBound = 93.9409657729240802468132361164 * std::max(eps/(1.-0.55),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, plus, x = 0.55 = 0.55";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.6);
    refVal   = 105.683586494539590277664890631;
    errBound = 105.683586494539590277664890631 * std::max(eps/(1.-0.6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, plus, x = 0.6 = 0.6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.65);
    refVal   = 120.781241708045246031617017864;
    errBound = 120.781241708045246031617017864 * std::max(eps/(1.-0.65),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, plus, x = 0.65 = 0.65";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.7);
    refVal   = 140.911448659386120370219854175;
    errBound = 140.911448659386120370219854175 * std::max(eps/(1.-0.7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, plus, x = 0.7 = 0.7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.75);
    refVal   = 169.093738391263344444263825009;
    errBound = 169.093738391263344444263825009 * std::max(eps/(1.-0.75),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, plus, x = 0.75 = 0.75";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.8);
    refVal   = 211.367172989079180555329781262;
    errBound = 211.367172989079180555329781262 * std::max(eps/(1.-0.8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, plus, x = 0.8 = 0.8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.85);
    refVal   = 281.822897318772240740439708349;
    errBound = 281.822897318772240740439708349 * std::max(eps/(1.-0.85),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, plus, x = 0.85 = 0.85";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9);
    refVal   = 422.734345978158361110659562524;
    errBound = 422.734345978158361110659562524 * std::max(eps/(1.-0.9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, plus, x = 0.9 = 0.9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.95);
    refVal   = 845.468691956316722221319125047;
    errBound = 845.468691956316722221319125047 * std::max(eps/(1.-0.95),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, plus, x = 0.95 = 0.95";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99);
    refVal   = 4227.34345978158361110659562524;
    errBound = 4227.34345978158361110659562524 * std::max(eps/(1.-0.99),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, plus, x = 1-10^-2 = 0.99";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999);
    refVal   = 42273.4345978158361110659562524;
    errBound = 42273.4345978158361110659562524 * std::max(eps/(1.-0.999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, plus, x = 1-10^-3 = 0.999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999);
    refVal   = 422734.345978158361110659562524;
    errBound = 422734.345978158361110659562524 * std::max(eps/(1.-0.9999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, plus, x = 1-10^-4 = 0.9999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999);
    refVal   = 4.22734345978158361110659562524e6;
    errBound = 4.22734345978158361110659562524e6 * std::max(eps/(1.-0.99999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, plus, x = 1-10^-5 = 0.99999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999);
    refVal   = 4.22734345978158361110659562524e7;
    errBound = 4.22734345978158361110659562524e7 * std::max(eps/(1.-0.999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, plus, x = 1-10^-6 = 0.999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999);
    refVal   = 4.22734345978158361110659562524e8;
    errBound = 4.22734345978158361110659562524e8 * std::max(eps/(1.-0.9999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, plus, x = 1-10^-7 = 0.9999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999999);
    refVal   = 4.22734345978158361110659562524e9;
    errBound = 4.22734345978158361110659562524e9 * std::max(eps/(1.-0.99999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, plus, x = 1-10^-8 = 0.99999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999);
    refVal   = 4.22734345978158361110659562524e10;
    errBound = 4.22734345978158361110659562524e10 * std::max(eps/(1.-0.999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, plus, x = 1-10^-9 = 0.999999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999999);
    refVal   = 4.22734345978158361110659562524e11;
    errBound = 4.22734345978158361110659562524e11 * std::max(eps/(1.-0.9999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, FullNormalDoubles, plus, x = 1-10^-10 = 0.9999999999";
  }
}



TEST(AggQXspace,FullDeltaEvaluation)
{
  auto& eval = AggQ_delta;
  EXPECT_DOUBLE_EQ(eval(testAs,testLM,testNF),-13.005021117466894651792463195)
    << "^-- AggQ, Full, delta";
}







TEST(AggQXspace,TruncAs0NormalDoublesDeltaEvaluation)
{
  auto eval = AggQ_delta.truncate(0);
  EXPECT_DOUBLE_EQ(eval(testAs,testLM,testNF),1.)
    << "^-- AggQ, TruncAs0NormalDoubles, delta";
}







TEST(AggQXspace,TruncAs1NormalDoublesDeltaEvaluation)
{
  auto eval = AggQ_delta.truncate(1);
  EXPECT_DOUBLE_EQ(eval(testAs,testLM,testNF),0.166666666666666666666666666667)
    << "^-- AggQ, TruncAs1NormalDoubles, delta";
}



TEST(AggQXspace,TruncAs2PowersOf4RegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto eval = AggQ_reg.truncate(2);

  {
    testVal  = eval(testAs,testLM,testNF,9.31322574615478515625e-10);
    refVal   = 5.29414259274318399053731796246e9;
    errBound = 5.29414259274318399053731796246e9 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, reg, x = 4^-15 = 9.31322574615478515625e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.7252902984619140625e-9);
    refVal   = 1.32353497521385245731012518388e9;
    errBound = 1.32353497521385245731012518388e9 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, reg, x = 4^-14 = 3.7252902984619140625e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.490116119384765625e-8);
    refVal   = 3.30883180286914360553381241384e8;
    errBound = 3.30883180286914360553381241384e8 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, reg, x = 4^-13 = 1.490116119384765625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,5.9604644775390625e-8);
    refVal   = 8.27203288504240101069948677175e7;
    errBound = 8.27203288504240101069948677175e7 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, reg, x = 4^-12 = 5.9604644775390625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,2.384185791015625e-7);
    refVal   = 2.06797017924109224289531081386e7;
    errBound = 2.06797017924109224289531081386e7 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, reg, x = 4^-11 = 2.384185791015625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,9.5367431640625e-7);
    refVal   = 5.16962000082860680023986806834e6;
    errBound = 5.16962000082860680023986806834e6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, reg, x = 4^-10 = 9.5367431640625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.814697265625e-6);
    refVal   = 1.29216436340470034741870406937e6;
    errBound = 1.29216436340470034741870406937e6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, reg, x = 4^-9 = 3.814697265625e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0000152587890625);
    refVal   = 322855.76722222995955354345078;
    errBound = 322855.76722222995955354345078 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, reg, x = 4^-8 = 0.0000152587890625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00006103515625);
    refVal   = 80575.0976234540345185036231447;
    errBound = 80575.0976234540345185036231447 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, reg, x = 4^-7 = 0.00006103515625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.000244140625);
    refVal   = 20043.235615484097645744702465;
    errBound = 20043.235615484097645744702465 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, reg, x = 4^-6 = 0.000244140625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0009765625);
    refVal   = 4941.0527599141165967776897082;
    errBound = 4941.0527599141165967776897082 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, reg, x = 4^-5 = 0.0009765625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00390625);
    refVal   = 1189.40504775026378431959087303;
    errBound = 1189.40504775026378431959087303 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, reg, x = 4^-4 = 0.00390625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.015625);
    refVal   = 269.141951246705427023882496576;
    errBound = 269.141951246705427023882496576 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, reg, x = 4^-3 = 0.015625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0625);
    refVal   = 51.1952023652574458171081067469;
    errBound = 51.1952023652574458171081067469 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, reg, x = 4^-2 = 0.0625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9375);
    refVal   = -2.45153281354499575288045860957;
    errBound = 2.45153281354499575288045860957 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, reg, x = 1-4^-2 = 0.9375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.984375);
    refVal   = -2.329304552662660586578788074;
    errBound = 2.329304552662660586578788074 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, reg, x = 1-4^-3 = 0.984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99609375);
    refVal   = -2.16571491805210376307697196242;
    errBound = 2.16571491805210376307697196242 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, reg, x = 1-4^-4 = 0.99609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9990234375);
    refVal   = -1.99428778137420063167726208063;
    errBound = 1.99428778137420063167726208063 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, reg, x = 1-4^-5 = 0.9990234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999755859375);
    refVal   = -1.82133553593995789882786003422;
    errBound = 1.82133553593995789882786003422 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, reg, x = 1-4^-6 = 0.999755859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99993896484375);
    refVal   = -1.64810043484687708684538840633;
    errBound = 1.64810043484687708684538840633 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, reg, x = 1-4^-7 = 0.99993896484375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999847412109375);
    refVal   = -1.4748186174546028450995511737;
    errBound = 1.4748186174546028450995511737 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, reg, x = 1-4^-8 = 0.9999847412109375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999996185302734375);
    refVal   = -1.30153108280922476297847223443;
    errBound = 1.30153108280922476297847223443 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, reg, x = 1-4^-9 = 0.999996185302734375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999904632568359375);
    refVal   = -1.12824360696337886804495006526;
    errBound = 1.12824360696337886804495006526 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, reg, x = 1-4^-10 = 0.99999904632568359375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999997615814208984375);
    refVal   = -0.954956517699302458254194401675;
    errBound = 0.954956517699302458254194401675 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, reg, x = 1-4^-11 = 0.9999997615814208984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999940395355224609375);
    refVal   = -0.781669618041996833183763238977;
    errBound = 0.781669618041996833183763238977 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, reg, x = 1-4^-12 = 0.999999940395355224609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999998509883880615234375);
    refVal   = -0.608382789026144581777576742322;
    errBound = 0.608382789026144581777576742322 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, reg, x = 1-4^-13 = 0.99999998509883880615234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999962747097015380859375);
    refVal   = -0.435095983480560129414532523947;
    errBound = 0.435095983480560129414532523947 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, reg, x = 1-4^-14 = 0.9999999962747097015380859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999068677425384521484375);
    refVal   = -0.261809185255016508615877830023;
    errBound = 0.261809185255016508615877830023 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, reg, x = 1-4^-15 = 0.999999999068677425384521484375";
  }
}

TEST(AggQXspace,TruncAs2PowersOf4PlusEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;
  auto eval = AggQ_plus.truncate(2);
  {
    testVal  = eval(testAs,testLM,testNF,9.31322574615478515625e-10);
    refVal   = 2.86111111377572847985368183904;
    errBound = 2.86111111377572847985368183904 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, plus, x = 4^-15 = 9.31322574615478515625e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.7252902984619140625e-9);
    refVal   = 2.86111112176958061586081383243;
    errBound = 2.86111112176958061586081383243 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, plus, x = 4^-14 = 3.7252902984619140625e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.490116119384765625e-8);
    refVal   = 2.8611111537449896065806456063;
    errBound = 2.8611111537449896065806456063 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, plus, x = 4^-13 = 1.490116119384765625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,5.9604644775390625e-8);
    refVal   = 2.86111128164663271652125284865;
    errBound = 2.86111128164663271652125284865 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, plus, x = 4^-12 = 5.9604644775390625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,2.384185791015625e-7);
    refVal   = 2.86111179325331950929100205154;
    errBound = 2.86111179325331950929100205154 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, plus, x = 4^-11 = 2.384185791015625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,9.5367431640625e-7);
    refVal   = 2.86111383968189633020474877281;
    errBound = 2.86111383968189633020474877281 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, plus, x = 4^-10 = 9.5367431640625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.814697265625e-6);
    refVal   = 2.86112202542547812114422704826;
    errBound = 2.86112202542547812114422704826 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, plus, x = 4^-9 = 3.814697265625e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0000152587890625);
    refVal   = 2.86115476886820443698447818384;
    errBound = 2.86115476886820443698447818384 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, plus, x = 4^-8 = 0.0000152587890625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00006103515625);
    refVal   = 2.86128575013394643499019986843;
    errBound = 2.86128575013394643499019986843 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, plus, x = 4^-7 = 0.00006103515625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.000244140625);
    refVal   = 2.86180979514312847646180979514;
    errBound = 2.86180979514312847646180979514 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, plus, x = 4^-6 = 0.000244140625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0009765625);
    refVal   = 2.86390789616596068208971434778;
    errBound = 2.86390789616596068208971434778 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, plus, x = 4^-5 = 0.0009765625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00390625);
    refVal   = 2.87233115468409586056644880174;
    errBound = 2.87233115468409586056644880174 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, plus, x = 4^-4 = 0.00390625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.015625);
    refVal   = 2.90652557319223985890652557319;
    errBound = 2.90652557319223985890652557319 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, plus, x = 4^-3 = 0.015625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0625);
    refVal   = 3.05185185185185185185185185185;
    errBound = 3.05185185185185185185185185185 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, plus, x = 4^-2 = 0.0625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9375);
    refVal   = 45.7777777777777777777777777778;
    errBound = 45.7777777777777777777777777778 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, plus, x = 1-4^-2 = 0.9375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.984375);
    refVal   = 183.111111111111111111111111111;
    errBound = 183.111111111111111111111111111 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, plus, x = 1-4^-3 = 0.984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99609375);
    refVal   = 732.444444444444444444444444444;
    errBound = 732.444444444444444444444444444 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, plus, x = 1-4^-4 = 0.99609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9990234375);
    refVal   = 2929.77777777777777777777777778;
    errBound = 2929.77777777777777777777777778 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, plus, x = 1-4^-5 = 0.9990234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999755859375);
    refVal   = 11719.1111111111111111111111111;
    errBound = 11719.1111111111111111111111111 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, plus, x = 1-4^-6 = 0.999755859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99993896484375);
    refVal   = 46876.4444444444444444444444444;
    errBound = 46876.4444444444444444444444444 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, plus, x = 1-4^-7 = 0.99993896484375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999847412109375);
    refVal   = 187505.777777777777777777777778;
    errBound = 187505.777777777777777777777778 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, plus, x = 1-4^-8 = 0.9999847412109375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999996185302734375);
    refVal   = 750023.111111111111111111111111;
    errBound = 750023.111111111111111111111111 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, plus, x = 1-4^-9 = 0.999996185302734375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999904632568359375);
    refVal   = 3.00009244444444444444444444444e6;
    errBound = 3.00009244444444444444444444444e6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, plus, x = 1-4^-10 = 0.99999904632568359375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999997615814208984375);
    refVal   = 1.20003697777777777777777777778e7;
    errBound = 1.20003697777777777777777777778e7 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, plus, x = 1-4^-11 = 0.9999997615814208984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999940395355224609375);
    refVal   = 4.80014791111111111111111111111e7;
    errBound = 4.80014791111111111111111111111e7 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, plus, x = 1-4^-12 = 0.999999940395355224609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999998509883880615234375);
    refVal   = 1.92005916444444444444444444444e8;
    errBound = 1.92005916444444444444444444444e8 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, plus, x = 1-4^-13 = 0.99999998509883880615234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999962747097015380859375);
    refVal   = 7.68023665777777777777777777778e8;
    errBound = 7.68023665777777777777777777778e8 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, plus, x = 1-4^-14 = 0.9999999962747097015380859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999068677425384521484375);
    refVal   = 3.07209466311111111111111111111e9;
    errBound = 3.07209466311111111111111111111e9 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2PowersOf4, plus, x = 1-4^-15 = 0.999999999068677425384521484375";
  }
}



TEST(AggQXspace,TruncAs2NormalDoublesRegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto eval = AggQ_reg.truncate(2);

  {
    testVal  = eval(testAs,testLM,testNF,1.e-10);
    refVal   = 4.93055541217086169017713479197e10;
    errBound = 4.93055541217086169017713479197e10 * std::max(eps/(1.-1.e-10),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, reg, x = 10^-10 = 1.e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-9);
    refVal   = 4.93055444186149963669208253487e9;
    errBound = 4.93055444186149963669208253487e9 * std::max(eps/(1.-1.e-9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, reg, x = 10^-9 = 1.e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-8);
    refVal   = 4.93054711693592074885008002798e8;
    errBound = 4.93054711693592074885008002798e8 * std::max(eps/(1.-1.e-8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, reg, x = 10^-8 = 1.e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-7);
    refVal   = 4.93049352742175628903846823955e7;
    errBound = 4.93049352742175628903846823955e7 * std::max(eps/(1.-1.e-7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, reg, x = 10^-7 = 1.e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-6);
    refVal   = 4.93011667250521238725087019037e6;
    errBound = 4.93011667250521238725087019037e6 * std::max(eps/(1.-1.e-6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, reg, x = 10^-6 = 1.e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-5);
    refVal   = 492759.956379484946611314159661;
    errBound = 492759.956379484946611314159661 * std::max(eps/(1.-1.e-5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, reg, x = 10^-5 = 1.e-5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0001);
    refVal   = 49119.1870963566944858589391231;
    errBound = 49119.1870963566944858589391231 * std::max(eps/(1.-0.0001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, reg, x = 10^-4 = 0.0001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.001);
    refVal   = 4823.39568733635666913719266506;
    errBound = 4823.39568733635666913719266506 * std::max(eps/(1.-0.001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, reg, x = 10^-3 = 0.001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.01);
    refVal   = 439.029046494628994408050743745;
    errBound = 439.029046494628994408050743745 * std::max(eps/(1.-0.01),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, reg, x = 10^-2 = 0.01";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.05);
    refVal   = 68.3765387743760604671858880117;
    errBound = 68.3765387743760604671858880117 * std::max(eps/(1.-0.05),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, reg, x = 0.05 = 0.05";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.1);
    refVal   = 26.4260494756299609412697194932;
    errBound = 26.4260494756299609412697194932 * std::max(eps/(1.-0.1),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, reg, x = 0.1 = 0.1";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.15);
    refVal   = 13.6021303462775014870743000624;
    errBound = 13.6021303462775014870743000624 * std::max(eps/(1.-0.15),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, reg, x = 0.15 = 0.15";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.2);
    refVal   = 7.6715854369708140715132152797;
    errBound = 7.6715854369708140715132152797 * std::max(eps/(1.-0.2),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, reg, x = 0.2 = 0.2";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.25);
    refVal   = 4.37179132155189211402510619226;
    errBound = 4.37179132155189211402510619226 * std::max(eps/(1.-0.25),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, reg, x = 0.25 = 0.25";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.3);
    refVal   = 2.33122899299787909504340529541;
    errBound = 2.33122899299787909504340529541 * std::max(eps/(1.-0.3),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, reg, x = 0.3 = 0.3";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.35);
    refVal   = 0.980634518861950321346751020208;
    errBound = 0.980634518861950321346751020208 * std::max(eps/(1.-0.35),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, reg, x = 0.35 = 0.35";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.4);
    refVal   = 0.0438089714550182727677507850283;
    errBound = 0.0438089714550182727677507850283 * std::max(eps/(1.-0.4),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, reg, x = 0.4 = 0.4";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.45);
    refVal   = -0.628311679536927586292383672361;
    errBound = 0.628311679536927586292383672361 * std::max(eps/(1.-0.45),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, reg, x = 0.45 = 0.45";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.5);
    refVal   = -1.12264403993169761328628485722;
    errBound = 1.12264403993169761328628485722 * std::max(eps/(1.-0.5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, reg, x = 0.5 = 0.5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.55);
    refVal   = -1.49290631524361361105591258342;
    errBound = 1.49290631524361361105591258342 * std::max(eps/(1.-0.55),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, reg, x = 0.55 = 0.55";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.6);
    refVal   = -1.77380306159488116972566010984;
    errBound = 1.77380306159488116972566010984 * std::max(eps/(1.-0.6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, reg, x = 0.6 = 0.6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.65);
    refVal   = -1.9884942941267986114896465105;
    errBound = 1.9884942941267986114896465105 * std::max(eps/(1.-0.65),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, reg, x = 0.65 = 0.65";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.7);
    refVal   = -2.1527338909732055764388350645;
    errBound = 2.1527338909732055764388350645 * std::max(eps/(1.-0.7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, reg, x = 0.7 = 0.7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.75);
    refVal   = -2.27720405980126644337415511514;
    errBound = 2.27720405980126644337415511514 * std::max(eps/(1.-0.75),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, reg, x = 0.75 = 0.75";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.8);
    refVal   = -2.36871981064039869977614929668;
    errBound = 2.36871981064039869977614929668 * std::max(eps/(1.-0.8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, reg, x = 0.8 = 0.8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.85);
    refVal   = -2.43040984994218161805525645894;
    errBound = 2.43040984994218161805525645894 * std::max(eps/(1.-0.85),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, reg, x = 0.85 = 0.85";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9);
    refVal   = -2.45993289902710457195002722186;
    errBound = 2.45993289902710457195002722186 * std::max(eps/(1.-0.9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, reg, x = 0.9 = 0.9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.95);
    refVal   = -2.43863918416962521639469861599;
    errBound = 2.43863918416962521639469861599 * std::max(eps/(1.-0.95),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, reg, x = 0.95 = 0.95";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99);
    refVal   = -2.2784526027252488173063417858;
    errBound = 2.2784526027252488173063417858 * std::max(eps/(1.-0.99),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, reg, x = 1-10^-2 = 0.99";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999);
    refVal   = -1.99724000822845995505772280423;
    errBound = 1.99724000822845995505772280423 * std::max(eps/(1.-0.999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, reg, x = 1-10^-3 = 0.999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999);
    refVal   = -1.70980729226882099600593085844;
    errBound = 1.70980729226882099600593085844 * std::max(eps/(1.-0.9999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, reg, x = 1-10^-4 = 0.9999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999);
    refVal   = -1.42199718013068885468785651056;
    errBound = 1.42199718013068885468785651056 * std::max(eps/(1.-0.99999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, reg, x = 1-10^-5 = 0.99999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999);
    refVal   = -1.13417275441423686871521078012;
    errBound = 1.13417275441423686871521078012 * std::max(eps/(1.-0.999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, reg, x = 1-10^-6 = 0.999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999);
    refVal   = -0.846349229828615177349449272861;
    errBound = 0.846349229828615177349449272861 * std::max(eps/(1.-0.9999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, reg, x = 1-10^-7 = 0.9999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999999);
    refVal   = -0.55852602850402734106135157994;
    errBound = 0.55852602850402734106135157994 * std::max(eps/(1.-0.99999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, reg, x = 1-10^-8 = 0.99999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999);
    refVal   = -0.270702882819329056049449894195;
    errBound = 0.270702882819329056049449894195 * std::max(eps/(1.-0.999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, reg, x = 1-10^-9 = 0.999999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999999);
    refVal   = 0.017120254970011746188719718967;
    errBound = 0.017120254970011746188719718967 * std::max(eps/(1.-0.9999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, reg, x = 1-10^-10 = 0.9999999999";
  }
}

TEST(AggQXspace,TruncAs2NormalDoublesPlusEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;
  auto eval = AggQ_plus.truncate(2);
  {
    testVal  = eval(testAs,testLM,testNF,1.e-10);
    refVal   = 2.86111111139722222225083333334;
    errBound = 2.86111111139722222225083333334 * std::max(eps/(1.-1.e-10),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, plus, x = 10^-10 = 1.e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-9);
    refVal   = 2.86111111397222222508333333619;
    errBound = 2.86111111397222222508333333619 * std::max(eps/(1.-1.e-9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, plus, x = 10^-9 = 1.e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-8);
    refVal   = 2.86111113972222250833333619444;
    errBound = 2.86111113972222250833333619444 * std::max(eps/(1.-1.e-8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, plus, x = 10^-8 = 1.e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-7);
    refVal   = 2.86111139722225083333619444473;
    errBound = 2.86111139722225083333619444473 * std::max(eps/(1.-1.e-7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, plus, x = 10^-7 = 1.e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-6);
    refVal   = 2.86111397222508333619444730556;
    errBound = 2.86111397222508333619444730556 * std::max(eps/(1.-1.e-6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, plus, x = 10^-6 = 1.e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-5);
    refVal   = 2.86113972250833619447305584167;
    errBound = 2.86113972250833619447305584167 * std::max(eps/(1.-1.e-5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, plus, x = 10^-5 = 1.e-5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0001);
    refVal   = 2.86139725083619473058416952806;
    errBound = 2.86139725083619473058416952806 * std::max(eps/(1.-0.0001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, plus, x = 10^-4 = 0.0001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.001);
    refVal   = 2.86397508619730841953064175286;
    errBound = 2.86397508619730841953064175286 * std::max(eps/(1.-0.001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, plus, x = 10^-3 = 0.001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.01);
    refVal   = 2.89001122334455667789001122334;
    errBound = 2.89001122334455667789001122334 * std::max(eps/(1.-0.01),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, plus, x = 10^-2 = 0.01";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.05);
    refVal   = 3.01169590643274853801169590643;
    errBound = 3.01169590643274853801169590643 * std::max(eps/(1.-0.05),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, plus, x = 0.05 = 0.05";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.1);
    refVal   = 3.17901234567901234567901234568;
    errBound = 3.17901234567901234567901234568 * std::max(eps/(1.-0.1),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, plus, x = 0.1 = 0.1";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.15);
    refVal   = 3.36601307189542483660130718954;
    errBound = 3.36601307189542483660130718954 * std::max(eps/(1.-0.15),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, plus, x = 0.15 = 0.15";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.2);
    refVal   = 3.57638888888888888888888888889;
    errBound = 3.57638888888888888888888888889 * std::max(eps/(1.-0.2),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, plus, x = 0.2 = 0.2";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.25);
    refVal   = 3.81481481481481481481481481481;
    errBound = 3.81481481481481481481481481481 * std::max(eps/(1.-0.25),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, plus, x = 0.25 = 0.25";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.3);
    refVal   = 4.0873015873015873015873015873;
    errBound = 4.0873015873015873015873015873 * std::max(eps/(1.-0.3),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, plus, x = 0.3 = 0.3";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.35);
    refVal   = 4.40170940170940170940170940171;
    errBound = 4.40170940170940170940170940171 * std::max(eps/(1.-0.35),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, plus, x = 0.35 = 0.35";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.4);
    refVal   = 4.76851851851851851851851851852;
    errBound = 4.76851851851851851851851851852 * std::max(eps/(1.-0.4),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, plus, x = 0.4 = 0.4";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.45);
    refVal   = 5.20202020202020202020202020202;
    errBound = 5.20202020202020202020202020202 * std::max(eps/(1.-0.45),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, plus, x = 0.45 = 0.45";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.5);
    refVal   = 5.72222222222222222222222222222;
    errBound = 5.72222222222222222222222222222 * std::max(eps/(1.-0.5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, plus, x = 0.5 = 0.5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.55);
    refVal   = 6.35802469135802469135802469136;
    errBound = 6.35802469135802469135802469136 * std::max(eps/(1.-0.55),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, plus, x = 0.55 = 0.55";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.6);
    refVal   = 7.15277777777777777777777777778;
    errBound = 7.15277777777777777777777777778 * std::max(eps/(1.-0.6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, plus, x = 0.6 = 0.6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.65);
    refVal   = 8.1746031746031746031746031746;
    errBound = 8.1746031746031746031746031746 * std::max(eps/(1.-0.65),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, plus, x = 0.65 = 0.65";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.7);
    refVal   = 9.53703703703703703703703703704;
    errBound = 9.53703703703703703703703703704 * std::max(eps/(1.-0.7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, plus, x = 0.7 = 0.7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.75);
    refVal   = 11.4444444444444444444444444444;
    errBound = 11.4444444444444444444444444444 * std::max(eps/(1.-0.75),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, plus, x = 0.75 = 0.75";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.8);
    refVal   = 14.3055555555555555555555555556;
    errBound = 14.3055555555555555555555555556 * std::max(eps/(1.-0.8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, plus, x = 0.8 = 0.8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.85);
    refVal   = 19.0740740740740740740740740741;
    errBound = 19.0740740740740740740740740741 * std::max(eps/(1.-0.85),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, plus, x = 0.85 = 0.85";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9);
    refVal   = 28.6111111111111111111111111111;
    errBound = 28.6111111111111111111111111111 * std::max(eps/(1.-0.9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, plus, x = 0.9 = 0.9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.95);
    refVal   = 57.2222222222222222222222222222;
    errBound = 57.2222222222222222222222222222 * std::max(eps/(1.-0.95),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, plus, x = 0.95 = 0.95";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99);
    refVal   = 286.111111111111111111111111111;
    errBound = 286.111111111111111111111111111 * std::max(eps/(1.-0.99),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, plus, x = 1-10^-2 = 0.99";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999);
    refVal   = 2861.11111111111111111111111111;
    errBound = 2861.11111111111111111111111111 * std::max(eps/(1.-0.999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, plus, x = 1-10^-3 = 0.999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999);
    refVal   = 28611.1111111111111111111111111;
    errBound = 28611.1111111111111111111111111 * std::max(eps/(1.-0.9999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, plus, x = 1-10^-4 = 0.9999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999);
    refVal   = 286111.111111111111111111111111;
    errBound = 286111.111111111111111111111111 * std::max(eps/(1.-0.99999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, plus, x = 1-10^-5 = 0.99999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999);
    refVal   = 2.86111111111111111111111111111e6;
    errBound = 2.86111111111111111111111111111e6 * std::max(eps/(1.-0.999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, plus, x = 1-10^-6 = 0.999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999);
    refVal   = 2.86111111111111111111111111111e7;
    errBound = 2.86111111111111111111111111111e7 * std::max(eps/(1.-0.9999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, plus, x = 1-10^-7 = 0.9999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999999);
    refVal   = 2.86111111111111111111111111111e8;
    errBound = 2.86111111111111111111111111111e8 * std::max(eps/(1.-0.99999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, plus, x = 1-10^-8 = 0.99999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999);
    refVal   = 2.86111111111111111111111111111e9;
    errBound = 2.86111111111111111111111111111e9 * std::max(eps/(1.-0.999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, plus, x = 1-10^-9 = 0.999999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999999);
    refVal   = 2.86111111111111111111111111111e10;
    errBound = 2.86111111111111111111111111111e10 * std::max(eps/(1.-0.9999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AggQ, TruncAs2NormalDoubles, plus, x = 1-10^-10 = 0.9999999999";
  }
}



TEST(AggQXspace,TruncAs2NormalDoublesDeltaEvaluation)
{
  auto eval = AggQ_delta.truncate(2);
  EXPECT_DOUBLE_EQ(eval(testAs,testLM,testNF),-2.99305555555555555555555555556)
    << "^-- AggQ, TruncAs2NormalDoubles, delta";
}



TEST(AggQNspace,FullMoments)
{
  double refVal = 0., errBound = 0.;
  auto& rpd = AggQ;
  integration_engine_gsl engine;
  auto mom = make_mellin_moment(rpd, engine, testAs, testLM, testNF);
  std::tuple<integration_status, double, double> res;

  {
    SCOPED_TRACE("AggQ Full Mellin moment N=2");
    refVal = -67.7252905373040287136420003452;
    errBound = 67.7252905373040287136420003452 * 3.e-11;
    res = mom.integrate(2, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AggQ Full Mellin moment N=4");
    refVal = -104.842284451535213240184896169;
    errBound = 104.842284451535213240184896169 * 3.e-11;
    res = mom.integrate(4, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AggQ Full Mellin moment N=6");
    refVal = -120.999566862102122724113388623;
    errBound = 120.999566862102122724113388623 * 3.e-11;
    res = mom.integrate(6, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AggQ Full Mellin moment N=8");
    refVal = -132.073899050232453003521924324;
    errBound = 132.073899050232453003521924324 * 3.e-11;
    res = mom.integrate(8, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

}


TEST(AggQNspace,TruncAs0Moments)
{
  double refVal = 0., errBound = 0.;
  auto rpd = AggQ.truncate(0);
  integration_engine_gsl engine;
  auto mom = make_mellin_moment(rpd, engine, testAs, testLM, testNF);
  std::tuple<integration_status, double, double> res;

  {
    SCOPED_TRACE("AggQ TruncAs0 Mellin moment N=2");
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
    SCOPED_TRACE("AggQ TruncAs0 Mellin moment N=4");
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
    SCOPED_TRACE("AggQ TruncAs0 Mellin moment N=6");
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
    SCOPED_TRACE("AggQ TruncAs0 Mellin moment N=8");
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


TEST(AggQNspace,TruncAs1Moments)
{
  double refVal = 0., errBound = 0.;
  auto rpd = AggQ.truncate(1);
  integration_engine_gsl engine;
  auto mom = make_mellin_moment(rpd, engine, testAs, testLM, testNF);
  std::tuple<integration_status, double, double> res;

  {
    SCOPED_TRACE("AggQ TruncAs1 Mellin moment N=2");
    refVal = 0.166666666666666666666666666667;
    errBound = 0.166666666666666666666666666667 * 3.e-11;
    res = mom.integrate(2, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AggQ TruncAs1 Mellin moment N=4");
    refVal = 0.166666666666666666666666666667;
    errBound = 0.166666666666666666666666666667 * 3.e-11;
    res = mom.integrate(4, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AggQ TruncAs1 Mellin moment N=6");
    refVal = 0.166666666666666666666666666667;
    errBound = 0.166666666666666666666666666667 * 3.e-11;
    res = mom.integrate(6, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AggQ TruncAs1 Mellin moment N=8");
    refVal = 0.166666666666666666666666666667;
    errBound = 0.166666666666666666666666666667 * 3.e-11;
    res = mom.integrate(8, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

}


TEST(AggQNspace,TruncAs2Moments)
{
  double refVal = 0., errBound = 0.;
  auto rpd = AggQ.truncate(2);
  integration_engine_gsl engine;
  auto mom = make_mellin_moment(rpd, engine, testAs, testLM, testNF);
  std::tuple<integration_status, double, double> res;

  {
    SCOPED_TRACE("AggQ TruncAs2 Mellin moment N=2");
    refVal = -5.99601337448559670781893004115;
    errBound = 5.99601337448559670781893004115 * 3.e-11;
    res = mom.integrate(2, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AggQ TruncAs2 Mellin moment N=4");
    refVal = -8.76058194444444444444444444444;
    errBound = 8.76058194444444444444444444444 * 3.e-11;
    res = mom.integrate(4, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AggQ TruncAs2 Mellin moment N=6");
    refVal = -9.90825099345608739842623872426;
    errBound = 9.90825099345608739842623872426 * 3.e-11;
    res = mom.integrate(6, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AggQ TruncAs2 Mellin moment N=8");
    refVal = -10.7053290939998303467096334037;
    errBound = 10.7053290939998303467096334037 * 3.e-11;
    res = mom.integrate(8, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

}
