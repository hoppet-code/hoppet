/******************************************************************************
 * This file is part of libome                                                *
 * Copyright (C) 2025 Arnd Behring, Kay Schoenwald                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 ******************************************************************************/

#include <limits>
#include <gtest/gtest.h>
#include <ome/integration_engine_gsl.h>
#include <ome/mellin.h>
#include <ome/AqgQ.h>

using namespace ome;

const double eps = std::numeric_limits<double>::epsilon();
const double testAs = 0.25;
const double testLM = -5.;
const double testNF = 3.;



TEST(AqgQXspace,FullPowersOf4RegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto& eval = AqgQ_reg;

  {
    testVal  = eval(testAs,testLM,testNF,9.31322574615478515625e-10);
    refVal   = 2.82481581072541794270818699696e10;
    errBound = 2.82481581072541794270818699696e10 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullPowersOf4, reg, x = 4^-15 = 9.31322574615478515625e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.7252902984619140625e-9);
    refVal   = 7.06196716840664015660596191866e9;
    errBound = 7.06196716840664015660596191866e9 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullPowersOf4, reg, x = 4^-14 = 3.7252902984619140625e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.490116119384765625e-8);
    refVal   = 1.76543916948360917451504085075e9;
    errBound = 1.76543916948360917451504085075e9 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullPowersOf4, reg, x = 4^-13 = 1.490116119384765625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,5.9604644775390625e-8);
    refVal   = 4.41322337696496851907904642452e8;
    errBound = 4.41322337696496851907904642452e8 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullPowersOf4, reg, x = 4^-12 = 5.9604644775390625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,2.384185791015625e-7);
    refVal   = 1.10304572981570094454666080852e8;
    errBound = 1.10304572981570094454666080852e8 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullPowersOf4, reg, x = 4^-11 = 2.384185791015625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,9.5367431640625e-7);
    refVal   = 2.75585830089805564636707508159e7;
    errBound = 2.75585830089805564636707508159e7 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullPowersOf4, reg, x = 4^-10 = 9.5367431640625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.814697265625e-6);
    refVal   = 6.87817462328847682794846571285e6;
    errBound = 6.87817462328847682794846571285e6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullPowersOf4, reg, x = 4^-9 = 3.814697265625e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0000152587890625);
    refVal   = 1.712334399492663334426254536e6;
    errBound = 1.712334399492663334426254536e6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullPowersOf4, reg, x = 4^-8 = 0.0000152587890625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00006103515625);
    refVal   = 423756.486847530753824124326988;
    errBound = 423756.486847530753824124326988 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullPowersOf4, reg, x = 4^-7 = 0.00006103515625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.000244140625);
    refVal   = 103482.271661838678816977714173;
    errBound = 103482.271661838678816977714173 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullPowersOf4, reg, x = 4^-6 = 0.000244140625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0009765625);
    refVal   = 24567.9408526892128827361053737;
    errBound = 24567.9408526892128827361053737 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullPowersOf4, reg, x = 4^-5 = 0.0009765625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00390625);
    refVal   = 5508.85704959744230188032003775;
    errBound = 5508.85704959744230188032003775 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullPowersOf4, reg, x = 4^-4 = 0.00390625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.015625);
    refVal   = 1102.95450547543514964766462084;
    errBound = 1102.95450547543514964766462084 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullPowersOf4, reg, x = 4^-3 = 0.015625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0625);
    refVal   = 174.550559615823031916830932261;
    errBound = 174.550559615823031916830932261 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullPowersOf4, reg, x = 4^-2 = 0.0625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9375);
    refVal   = -15.2449323029017029598183855054;
    errBound = 15.2449323029017029598183855054 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullPowersOf4, reg, x = 1-4^-2 = 0.9375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.984375);
    refVal   = -39.4074334179303944713710556933;
    errBound = 39.4074334179303944713710556933 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullPowersOf4, reg, x = 1-4^-3 = 0.984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99609375);
    refVal   = -74.4446716161665494656047105871;
    errBound = 74.4446716161665494656047105871 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullPowersOf4, reg, x = 1-4^-4 = 0.99609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9990234375);
    refVal   = -123.902176047606158531252524249;
    errBound = 123.902176047606158531252524249 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullPowersOf4, reg, x = 1-4^-5 = 0.9990234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999755859375);
    refVal   = -191.469962389733214703474714212;
    errBound = 191.469962389733214703474714212 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullPowersOf4, reg, x = 1-4^-6 = 0.999755859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99993896484375);
    refVal   = -280.907783145744864359440730493;
    errBound = 280.907783145744864359440730493 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullPowersOf4, reg, x = 1-4^-7 = 0.99993896484375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999847412109375);
    refVal   = -396.123838262943137179504315545;
    errBound = 396.123838262943137179504315545 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullPowersOf4, reg, x = 1-4^-8 = 0.9999847412109375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999996185302734375);
    refVal   = -541.233037869427348934646278245;
    errBound = 541.233037869427348934646278245 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullPowersOf4, reg, x = 1-4^-9 = 0.999996185302734375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999904632568359375);
    refVal   = -720.586258836041139451901968236;
    errBound = 720.586258836041139451901968236 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullPowersOf4, reg, x = 1-4^-10 = 0.99999904632568359375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999997615814208984375);
    refVal   = -938.78298279998796429018066858;
    errBound = 938.78298279998796429018066858 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullPowersOf4, reg, x = 1-4^-11 = 0.9999997615814208984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999940395355224609375);
    refVal   = -1200.67629953993066506834960779;
    errBound = 1200.67629953993066506834960779 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullPowersOf4, reg, x = 1-4^-12 = 0.999999940395355224609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999998509883880615234375);
    refVal   = -1511.37477310488079904867879697;
    errBound = 1511.37477310488079904867879697 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullPowersOf4, reg, x = 1-4^-13 = 0.99999998509883880615234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999962747097015380859375);
    refVal   = -1876.24310763554649242457082093;
    errBound = 1876.24310763554649242457082093 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullPowersOf4, reg, x = 1-4^-14 = 0.9999999962747097015380859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999068677425384521484375);
    refVal   = -2300.90237683626726558577268715;
    errBound = 2300.90237683626726558577268715 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullPowersOf4, reg, x = 1-4^-15 = 0.999999999068677425384521484375";
  }
}



TEST(AqgQXspace,FullNormalDoublesRegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto& eval = AqgQ_reg;

  {
    testVal  = eval(testAs,testLM,testNF,1.e-10);
    refVal   = 2.63082603791777098143297505151e11;
    errBound = 2.63082603791777098143297505151e11 * std::max(eps/(1.-1.e-10),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullNormalDoubles, reg, x = 10^-10 = 1.e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-9);
    refVal   = 2.63081394439389038485968567585e10;
    errBound = 2.63081394439389038485968567585e10 * std::max(eps/(1.-1.e-9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullNormalDoubles, reg, x = 10^-9 = 1.e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-8);
    refVal   = 2.63074087970705395188286194723e9;
    errBound = 2.63074087970705395188286194723e9 * std::max(eps/(1.-1.e-8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullNormalDoubles, reg, x = 10^-8 = 1.e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-7);
    refVal   = 2.63032352190006303596205291004e8;
    errBound = 2.63032352190006303596205291004e8 * std::max(eps/(1.-1.e-7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullNormalDoubles, reg, x = 10^-7 = 1.e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-6);
    refVal   = 2.62810011996653536830820687405e7;
    errBound = 2.62810011996653536830820687405e7 * std::max(eps/(1.-1.e-6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullNormalDoubles, reg, x = 10^-6 = 1.e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-5);
    refVal   = 2.61726225598597593020569042937e6;
    errBound = 2.61726225598597593020569042937e6 * std::max(eps/(1.-1.e-5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullNormalDoubles, reg, x = 10^-5 = 1.e-5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0001);
    refVal   = 257024.190353545405886337776397;
    errBound = 257024.190353545405886337776397 * std::max(eps/(1.-0.0001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullNormalDoubles, reg, x = 10^-4 = 0.0001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.001);
    refVal   = 23961.6796816816568702748178676;
    errBound = 23961.6796816816568702748178676 * std::max(eps/(1.-0.001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullNormalDoubles, reg, x = 10^-3 = 0.001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.01);
    refVal   = 1884.09418322918420335581288082;
    errBound = 1884.09418322918420335581288082 * std::max(eps/(1.-0.01),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullNormalDoubles, reg, x = 10^-2 = 0.01";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.05);
    refVal   = 240.852972898187802940871233045;
    errBound = 240.852972898187802940871233045 * std::max(eps/(1.-0.05),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullNormalDoubles, reg, x = 0.05 = 0.05";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.1);
    refVal   = 84.508430051282402421381348425;
    errBound = 84.508430051282402421381348425 * std::max(eps/(1.-0.1),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullNormalDoubles, reg, x = 0.1 = 0.1";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.15);
    refVal   = 42.3538889918790969329910661227;
    errBound = 42.3538889918790969329910661227 * std::max(eps/(1.-0.15),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullNormalDoubles, reg, x = 0.15 = 0.15";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.2);
    refVal   = 24.8844590679620286782490655898;
    errBound = 24.8844590679620286782490655898 * std::max(eps/(1.-0.2),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullNormalDoubles, reg, x = 0.2 = 0.2";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.25);
    refVal   = 16.169801524860865586644247638;
    errBound = 16.169801524860865586644247638 * std::max(eps/(1.-0.25),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullNormalDoubles, reg, x = 0.25 = 0.25";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.3);
    refVal   = 11.3583738296744914085707777989;
    errBound = 11.3583738296744914085707777989 * std::max(eps/(1.-0.3),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullNormalDoubles, reg, x = 0.3 = 0.3";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.35);
    refVal   = 8.5235073315059308608686626225;
    errBound = 8.5235073315059308608686626225 * std::max(eps/(1.-0.35),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullNormalDoubles, reg, x = 0.35 = 0.35";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.4);
    refVal   = 6.76222817022610150575227180004;
    errBound = 6.76222817022610150575227180004 * std::max(eps/(1.-0.4),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullNormalDoubles, reg, x = 0.4 = 0.4";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.45);
    refVal   = 5.59785565018806604945116821053;
    errBound = 5.59785565018806604945116821053 * std::max(eps/(1.-0.45),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullNormalDoubles, reg, x = 0.45 = 0.45";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.5);
    refVal   = 4.75269520094716912053202944183;
    errBound = 4.75269520094716912053202944183 * std::max(eps/(1.-0.5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullNormalDoubles, reg, x = 0.5 = 0.5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.55);
    refVal   = 4.04841174502991125774199439418;
    errBound = 4.04841174502991125774199439418 * std::max(eps/(1.-0.55),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullNormalDoubles, reg, x = 0.55 = 0.55";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.6);
    refVal   = 3.35643977808828948670345261336;
    errBound = 3.35643977808828948670345261336 * std::max(eps/(1.-0.6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullNormalDoubles, reg, x = 0.6 = 0.6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.65);
    refVal   = 2.56902327489450245173470049064;
    errBound = 2.56902327489450245173470049064 * std::max(eps/(1.-0.65),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullNormalDoubles, reg, x = 0.65 = 0.65";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.7);
    refVal   = 1.57716628358847373905871775411;
    errBound = 1.57716628358847373905871775411 * std::max(eps/(1.-0.7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullNormalDoubles, reg, x = 0.7 = 0.7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.75);
    refVal   = 0.245554829800327973656252024669;
    errBound = 0.245554829800327973656252024669 * std::max(eps/(1.-0.75),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullNormalDoubles, reg, x = 0.75 = 0.75";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.8);
    refVal   = -1.63058687362366185684608881825;
    errBound = 1.63058687362366185684608881825 * std::max(eps/(1.-0.8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullNormalDoubles, reg, x = 0.8 = 0.8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.85);
    refVal   = -4.42803595934598204645006938356;
    errBound = 4.42803595934598204645006938356 * std::max(eps/(1.-0.85),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullNormalDoubles, reg, x = 0.85 = 0.85";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9);
    refVal   = -9.02640589620356266098334919822;
    errBound = 9.02640589620356266098334919822 * std::max(eps/(1.-0.9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullNormalDoubles, reg, x = 0.9 = 0.9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.95);
    refVal   = -18.5273214987368645880930331905;
    errBound = 18.5273214987368645880930331905 * std::max(eps/(1.-0.95),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullNormalDoubles, reg, x = 0.95 = 0.95";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99);
    refVal   = -49.3364299109522939353626140649;
    errBound = 49.3364299109522939353626140649 * std::max(eps/(1.-0.99),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullNormalDoubles, reg, x = 1-10^-2 = 0.99";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999);
    refVal   = -122.91437465614291982536932247;
    errBound = 122.91437465614291982536932247 * std::max(eps/(1.-0.999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullNormalDoubles, reg, x = 1-10^-3 = 0.999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999);
    refVal   = -246.306040659118938709065534052;
    errBound = 246.306040659118938709065534052 * std::max(eps/(1.-0.9999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullNormalDoubles, reg, x = 1-10^-4 = 0.9999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999);
    refVal   = -436.993073706097945806746118527;
    errBound = 436.993073706097945806746118527 * std::max(eps/(1.-0.99999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullNormalDoubles, reg, x = 1-10^-5 = 0.99999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999);
    refVal   = -713.834673356106996712354513423;
    errBound = 713.834673356106996712354513423 * std::max(eps/(1.-0.999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullNormalDoubles, reg, x = 1-10^-6 = 0.999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999);
    refVal   = -1097.49466421117664655726530135;
    errBound = 1097.49466421117664655726530135 * std::max(eps/(1.-0.9999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullNormalDoubles, reg, x = 1-10^-7 = 0.9999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999999);
    refVal   = -1610.55888565636617005632973053;
    errBound = 1610.55888565636617005632973053 * std::max(eps/(1.-0.99999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullNormalDoubles, reg, x = 1-10^-8 = 0.99999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999);
    refVal   = -2277.55987429198865163788019632;
    errBound = 2277.55987429198865163788019632 * std::max(eps/(1.-0.999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullNormalDoubles, reg, x = 1-10^-9 = 0.999999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999999);
    refVal   = -3124.98137087834241174897599261;
    errBound = 3124.98137087834241174897599261 * std::max(eps/(1.-0.9999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- AqgQ, FullNormalDoubles, reg, x = 1-10^-10 = 0.9999999999";
  }
}



TEST(AqgQNspace,FullMoments)
{
  double refVal = 0., errBound = 0.;
  auto& rpd = AqgQ;
  integration_engine_gsl engine;
  auto mom = make_mellin_moment(rpd, engine, testAs, testLM, testNF);
  std::tuple<integration_status, double, double> res;

  {
    SCOPED_TRACE("AqgQ Full Mellin moment N=2");
    refVal = 0.370041006688864128977737475694;
    errBound = 0.370041006688864128977737475694 * 3.e-11;
    res = mom.integrate(2, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AqgQ Full Mellin moment N=4");
    refVal = -2.41999365942078596262046428092;
    errBound = 2.41999365942078596262046428092 * 3.e-11;
    res = mom.integrate(4, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AqgQ Full Mellin moment N=6");
    refVal = -2.41248281958931075363962699026;
    errBound = 2.41248281958931075363962699026 * 3.e-11;
    res = mom.integrate(6, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AqgQ Full Mellin moment N=8");
    refVal = -2.28185858612602839180936883132;
    errBound = 2.28185858612602839180936883132 * 3.e-11;
    res = mom.integrate(8, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("AqgQ Full Mellin moment N=10");
    refVal = -2.14572535946038682708224656616;
    errBound = 2.14572535946038682708224656616 * 3.e-11;
    res = mom.integrate(10, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

}
