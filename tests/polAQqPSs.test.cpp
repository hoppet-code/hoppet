/******************************************************************************
 * This file is part of libome                                                *
 * Copyright (C) 2025 Arnd Behring, Kay Schoenwald                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 ******************************************************************************/

#include <limits>
#include <gtest/gtest.h>
#include <ome/integration_engine_gsl.h>
#include <ome/mellin.h>
#include <ome/polAQqPSs.h>

using namespace ome;

const double eps = std::numeric_limits<double>::epsilon();
const double testAs = 0.25;
const double testLM = -5.;
const double testNF = 3.;



TEST(polAQqPSsXspace,FullPowersOf4RegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto& eval = polAQqPSs_reg;

  {
    testVal  = eval(testAs,testLM,testNF,9.31322574615478515625e-10);
    refVal   = 589.92745848515381058405166436;
    errBound = 589.92745848515381058405166436 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullPowersOf4, reg, x = 4^-15 = 9.31322574615478515625e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.7252902984619140625e-9);
    refVal   = 404.508666930815199629269112463;
    errBound = 404.508666930815199629269112463 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullPowersOf4, reg, x = 4^-14 = 3.7252902984619140625e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.490116119384765625e-8);
    refVal   = 255.729171546258045303061457322;
    errBound = 255.729171546258045303061457322 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullPowersOf4, reg, x = 4^-13 = 1.490116119384765625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,5.9604644775390625e-8);
    refVal   = 139.937053222465942531276814124;
    errBound = 139.937053222465942531276814124 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullPowersOf4, reg, x = 4^-12 = 5.9604644775390625e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,2.384185791015625e-7);
    refVal   = 53.4806277682027596631510918093;
    errBound = 53.4806277682027596631510918093 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullPowersOf4, reg, x = 4^-11 = 2.384185791015625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,9.5367431640625e-7);
    refVal   = -7.29108654475158827516095976797;
    errBound = 7.29108654475158827516095976797 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullPowersOf4, reg, x = 4^-10 = 9.5367431640625e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,3.814697265625e-6);
    refVal   = -46.0271078312026936726272106107;
    errBound = 46.0271078312026936726272106107 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullPowersOf4, reg, x = 4^-9 = 3.814697265625e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0000152587890625);
    refVal   = -66.3714186296231026937458519098;
    errBound = 66.3714186296231026937458519098 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullPowersOf4, reg, x = 4^-8 = 0.0000152587890625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00006103515625);
    refVal   = -71.9565828496139553164191724897;
    errBound = 71.9565828496139553164191724897 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullPowersOf4, reg, x = 4^-7 = 0.00006103515625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.000244140625);
    refVal   = -66.3941486251607640338672431222;
    errBound = 66.3941486251607640338672431222 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullPowersOf4, reg, x = 4^-6 = 0.000244140625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0009765625);
    refVal   = -53.2727866029634603883136784464;
    errBound = 53.2727866029634603883136784464 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullPowersOf4, reg, x = 4^-5 = 0.0009765625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.00390625);
    refVal   = -36.2095489149468225974281344526;
    errBound = 36.2095489149468225974281344526 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullPowersOf4, reg, x = 4^-4 = 0.00390625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.015625);
    refVal   = -19.0456198322635226040125618276;
    errBound = 19.0456198322635226040125618276 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullPowersOf4, reg, x = 4^-3 = 0.015625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0625);
    refVal   = -6.08584923577354661933984426047;
    errBound = 6.08584923577354661933984426047 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullPowersOf4, reg, x = 4^-2 = 0.0625";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9375);
    refVal   = -0.124890436737844573522467682269;
    errBound = 0.124890436737844573522467682269 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullPowersOf4, reg, x = 1-4^-2 = 0.9375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.984375);
    refVal   = -0.0340218471310461499688023521294;
    errBound = 0.0340218471310461499688023521294 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullPowersOf4, reg, x = 1-4^-3 = 0.984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99609375);
    refVal   = -0.00900400066694242494224478853687;
    errBound = 0.00900400066694242494224478853687 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullPowersOf4, reg, x = 1-4^-4 = 0.99609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9990234375);
    refVal   = -0.00234404877375438047735388940631;
    errBound = 0.00234404877375438047735388940631 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullPowersOf4, reg, x = 1-4^-5 = 0.9990234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999755859375);
    refVal   = -0.000602504966863507310058188297429;
    errBound = 0.000602504966863507310058188297429 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullPowersOf4, reg, x = 1-4^-6 = 0.999755859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99993896484375);
    refVal   = -0.00015313987596359173244313726157;
    errBound = 0.00015313987596359173244313726157 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullPowersOf4, reg, x = 1-4^-7 = 0.99993896484375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999847412109375);
    refVal   = -0.0000385172395263201569173556634501;
    errBound = 0.0000385172395263201569173556634501 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullPowersOf4, reg, x = 1-4^-8 = 0.9999847412109375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999996185302734375);
    refVal   = -9.58880699943087892071221064941e-6;
    errBound = 9.58880699943087892071221064941e-6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullPowersOf4, reg, x = 1-4^-9 = 0.999996185302734375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999904632568359375);
    refVal   = -2.3624670862248104569141647497e-6;
    errBound = 2.3624670862248104569141647497e-6 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullPowersOf4, reg, x = 1-4^-10 = 0.99999904632568359375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999997615814208984375);
    refVal   = -5.75783270046599156197543570465e-7;
    errBound = 5.75783270046599156197543570465e-7 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullPowersOf4, reg, x = 1-4^-11 = 0.9999997615814208984375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999940395355224609375);
    refVal   = -1.38700151702237182278230863405e-7;
    errBound = 1.38700151702237182278230863405e-7 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullPowersOf4, reg, x = 1-4^-12 = 0.999999940395355224609375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999998509883880615234375);
    refVal   = -3.29793105486082492818238771307e-8;
    errBound = 3.29793105486082492818238771307e-8 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullPowersOf4, reg, x = 1-4^-13 = 0.99999998509883880615234375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999962747097015380859375);
    refVal   = -7.72481886411978228298443483583e-9;
    errBound = 7.72481886411978228298443483583e-9 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullPowersOf4, reg, x = 1-4^-14 = 0.9999999962747097015380859375";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999068677425384521484375);
    refVal   = -1.777183345229454174682493864e-9;
    errBound = 1.777183345229454174682493864e-9 * 2048.*eps;
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullPowersOf4, reg, x = 1-4^-15 = 0.999999999068677425384521484375";
  }
}



TEST(polAQqPSsXspace,FullNormalDoublesRegEvaluation)
{
  double testVal = 0., refVal = 0., errBound = 0.;

  auto& eval = polAQqPSs_reg;

  {
    testVal  = eval(testAs,testLM,testNF,1.e-10);
    refVal   = 974.567208941906888222593840909;
    errBound = 974.567208941906888222593840909 * std::max(eps/(1.-1.e-10),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullNormalDoubles, reg, x = 10^-10 = 1.e-10";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-9);
    refVal   = 579.461379699703614867887742074;
    errBound = 579.461379699703614867887742074 * std::max(eps/(1.-1.e-9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullNormalDoubles, reg, x = 10^-9 = 1.e-9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-8);
    refVal   = 294.994138403650550754318855945;
    errBound = 294.994138403650550754318855945 * std::max(eps/(1.-1.e-8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullNormalDoubles, reg, x = 10^-8 = 1.e-8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-7);
    refVal   = 104.431231303151819278012088553;
    errBound = 104.431231303151819278012088553 * std::max(eps/(1.-1.e-7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullNormalDoubles, reg, x = 10^-7 = 1.e-7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-6);
    refVal   = -8.9597769294023099433934992139;
    errBound = 8.9597769294023099433934992139 * std::max(eps/(1.-1.e-6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullNormalDoubles, reg, x = 10^-6 = 1.e-6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,1.e-5);
    refVal   = -61.9009103422260706815572563641;
    errBound = 61.9009103422260706815572563641 * std::max(eps/(1.-1.e-5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullNormalDoubles, reg, x = 10^-5 = 1.e-5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.0001);
    refVal   = -71.067296345853365097899061449;
    errBound = 71.067296345853365097899061449 * std::max(eps/(1.-0.0001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullNormalDoubles, reg, x = 10^-4 = 0.0001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.001);
    refVal   = -53.0035491097247822502714408982;
    errBound = 53.0035491097247822502714408982 * std::max(eps/(1.-0.001),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullNormalDoubles, reg, x = 10^-3 = 0.001";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.01);
    refVal   = -24.3285654301047960850456564688;
    errBound = 24.3285654301047960850456564688 * std::max(eps/(1.-0.01),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullNormalDoubles, reg, x = 10^-2 = 0.01";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.05);
    refVal   = -7.70592156168025518488795236548;
    errBound = 7.70592156168025518488795236548 * std::max(eps/(1.-0.05),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullNormalDoubles, reg, x = 0.05 = 0.05";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.1);
    refVal   = -3.38362987098288417881951698037;
    errBound = 3.38362987098288417881951698037 * std::max(eps/(1.-0.1),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullNormalDoubles, reg, x = 0.1 = 0.1";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.15);
    refVal   = -1.84798152306432195993021509429;
    errBound = 1.84798152306432195993021509429 * std::max(eps/(1.-0.15),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullNormalDoubles, reg, x = 0.15 = 0.15";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.2);
    refVal   = -1.18023304130887052433534413223;
    errBound = 1.18023304130887052433534413223 * std::max(eps/(1.-0.2),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullNormalDoubles, reg, x = 0.2 = 0.2";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.25);
    refVal   = -0.86988391870649371005753475275;
    errBound = 0.86988391870649371005753475275 * std::max(eps/(1.-0.25),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullNormalDoubles, reg, x = 0.25 = 0.25";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.3);
    refVal   = -0.724334998005751105239585112817;
    errBound = 0.724334998005751105239585112817 * std::max(eps/(1.-0.3),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullNormalDoubles, reg, x = 0.3 = 0.3";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.35);
    refVal   = -0.657384994936130608051764644456;
    errBound = 0.657384994936130608051764644456 * std::max(eps/(1.-0.35),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullNormalDoubles, reg, x = 0.35 = 0.35";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.4);
    refVal   = -0.626498134999340616432511180528;
    errBound = 0.626498134999340616432511180528 * std::max(eps/(1.-0.4),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullNormalDoubles, reg, x = 0.4 = 0.4";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.45);
    refVal   = -0.609471019601062537073947410921;
    errBound = 0.609471019601062537073947410921 * std::max(eps/(1.-0.45),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullNormalDoubles, reg, x = 0.45 = 0.45";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.5);
    refVal   = -0.594404744417579761390279815127;
    errBound = 0.594404744417579761390279815127 * std::max(eps/(1.-0.5),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullNormalDoubles, reg, x = 0.5 = 0.5";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.55);
    refVal   = -0.574938201879723708158545432249;
    errBound = 0.574938201879723708158545432249 * std::max(eps/(1.-0.55),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullNormalDoubles, reg, x = 0.55 = 0.55";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.6);
    refVal   = -0.547807685949285700552385142599;
    errBound = 0.547807685949285700552385142599 * std::max(eps/(1.-0.6),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullNormalDoubles, reg, x = 0.6 = 0.6";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.65);
    refVal   = -0.511524150059708060203859217686;
    errBound = 0.511524150059708060203859217686 * std::max(eps/(1.-0.65),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullNormalDoubles, reg, x = 0.65 = 0.65";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.7);
    refVal   = -0.465622077075202725259844479054;
    errBound = 0.465622077075202725259844479054 * std::max(eps/(1.-0.7),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullNormalDoubles, reg, x = 0.7 = 0.7";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.75);
    refVal   = -0.410213599275956410465527636134;
    errBound = 0.410213599275956410465527636134 * std::max(eps/(1.-0.75),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullNormalDoubles, reg, x = 0.75 = 0.75";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.8);
    refVal   = -0.345707141453316609422431268331;
    errBound = 0.345707141453316609422431268331 * std::max(eps/(1.-0.8),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullNormalDoubles, reg, x = 0.8 = 0.8";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.85);
    refVal   = -0.272604410536640878328706289293;
    errBound = 0.272604410536640878328706289293 * std::max(eps/(1.-0.85),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullNormalDoubles, reg, x = 0.85 = 0.85";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9);
    refVal   = -0.191293802281915372447670058752;
    errBound = 0.191293802281915372447670058752 * std::max(eps/(1.-0.9),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullNormalDoubles, reg, x = 0.9 = 0.9";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.95);
    refVal   = -0.10163730882415640149085089018;
    errBound = 0.10163730882415640149085089018 * std::max(eps/(1.-0.95),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullNormalDoubles, reg, x = 0.95 = 0.95";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99);
    refVal   = -0.0222261270645372990468678699997;
    errBound = 0.0222261270645372990468678699997 * std::max(eps/(1.-0.99),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullNormalDoubles, reg, x = 1-10^-2 = 0.99";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999);
    refVal   = -0.00239891535253556149892065449652;
    errBound = 0.00239891535253556149892065449652 * std::max(eps/(1.-0.999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullNormalDoubles, reg, x = 1-10^-3 = 0.999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999);
    refVal   = -0.000249736468067421405440265118351;
    errBound = 0.000249736468067421405440265118351 * std::max(eps/(1.-0.9999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullNormalDoubles, reg, x = 1-10^-4 = 0.9999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999);
    refVal   = -0.0000252376465340800060130234752749;
    errBound = 0.0000252376465340800060130234752749 * std::max(eps/(1.-0.99999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullNormalDoubles, reg, x = 1-10^-5 = 0.99999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999);
    refVal   = -2.47889870328556386288655411611e-6;
    errBound = 2.47889870328556386288655411611e-6 * std::max(eps/(1.-0.999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullNormalDoubles, reg, x = 1-10^-6 = 0.999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999);
    refVal   = -2.36286810370102534243626205198e-7;
    errBound = 2.36286810370102534243626205198e-7 * std::max(eps/(1.-0.9999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullNormalDoubles, reg, x = 1-10^-7 = 0.9999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.99999999);
    refVal   = -2.17568530891433679868440748645e-8;
    errBound = 2.17568530891433679868440748645e-8 * std::max(eps/(1.-0.99999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullNormalDoubles, reg, x = 1-10^-8 = 0.99999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.999999999);
    refVal   = -1.91735179743997281663539099546e-9;
    errBound = 1.91735179743997281663539099546e-9 * std::max(eps/(1.-0.999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullNormalDoubles, reg, x = 1-10^-9 = 0.999999999";
  }
  {
    testVal  = eval(testAs,testLM,testNF,0.9999999999);
    refVal   = -1.58786774044255785643032460045e-10;
    errBound = 1.58786774044255785643032460045e-10 * std::max(eps/(1.-0.9999999999),2048.*eps);
    EXPECT_NEAR(testVal,refVal,errBound)
      << "^-- polAQqPSs, FullNormalDoubles, reg, x = 1-10^-10 = 0.9999999999";
  }
}



TEST(polAQqPSsNspace,FullMoments)
{
  double refVal = 0., errBound = 0.;
  auto& rpd = polAQqPSs;
  integration_engine_gsl engine;
  auto mom = make_mellin_moment(rpd, engine, testAs, testLM, testNF);
  std::tuple<integration_status, double, double> res;

  {
    SCOPED_TRACE("polAQqPSs Full Mellin moment N=2");
    refVal = -0.261888234416759729483466361737;
    errBound = 0.261888234416759729483466361737 * 3.e-11;
    res = mom.integrate(2, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAQqPSs Full Mellin moment N=4");
    refVal = -0.0759470470946909445137264805811;
    errBound = 0.0759470470946909445137264805811 * 3.e-11;
    res = mom.integrate(4, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAQqPSs Full Mellin moment N=6");
    refVal = -0.0392284865291078366221794416697;
    errBound = 0.0392284865291078366221794416697 * 3.e-11;
    res = mom.integrate(6, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAQqPSs Full Mellin moment N=8");
    refVal = -0.024128977611498422971485001866;
    errBound = 0.024128977611498422971485001866 * 3.e-11;
    res = mom.integrate(8, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

  {
    SCOPED_TRACE("polAQqPSs Full Mellin moment N=10");
    refVal = -0.0163608490693801339452428020083;
    errBound = 0.0163608490693801339452428020083 * 3.e-11;
    res = mom.integrate(10, 0., 3.e-11);
    EXPECT_EQ(std::get<0>(res),integration_status::success)
      << "Status code from gsl_integrate_mellin_moment didn't indicate success.";
    EXPECT_NEAR(std::get<1>(res),refVal,errBound)
      << "Integration result from gsl_integrate_mellin_moment didn't match the reference value.";
    EXPECT_LE(std::get<2>(res),errBound)
      << "Absolute integration error from gsl_integrate_mellin_moment didn't reach the target.";
  }

}
