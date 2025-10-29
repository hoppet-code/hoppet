/******************************************************************************
 * This file is part of libome                                                *
 * Copyright (C) 2025 Arnd Behring, Kay Schoenwald                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 ******************************************************************************/

#include <gtest/gtest.h>
#include <ome/piecewise.h>

using namespace ome;

template<typename Tnum>
class eval_id
{
  public:
    Tnum operator()(const Tnum x) const
    {
      return(x);
    };
};

template<typename Tnum>
class eval_shifted
{
  public:
    explicit
    eval_shifted(Tnum shift)
      : shift_(shift) {};

    Tnum operator()(const Tnum x) const
    {
      return(x+shift_);
    };

  private:
    Tnum shift_;
};


TEST(PiecewiseTest, DefaultConstruction)
{
  piecewise<double, eval_id<double>> pw;

  EXPECT_EQ(pw(1.),1.);
}

TEST(PiecewiseTest, InvalidConstructionThrows)
{
  using eid = eval_id<double>;
  using pw_type = piecewise<double, eid>;

  EXPECT_THROW(pw_type pw({},{}), inconsistent_piecewise)
    << "Inconsistency (empty pieces list) not detected";

  EXPECT_THROW(pw_type pw({0.},{eid(),eid(),eid()}), inconsistent_piecewise)
    << "Inconsistency (mismatch between # boundaries & # pieces) not detected";
  
  EXPECT_THROW(pw_type pw({0.,10.,20.},{eid(),eid(),eid()}), inconsistent_piecewise)
    << "Inconsistency (mismatch between # boundaries & # pieces) not detected";

  EXPECT_THROW(pw_type pw({0.,10.,20.},{eid()}), inconsistent_piecewise)
    << "Inconsistency (mismatch between # boundaries & # pieces) not detected";

  EXPECT_THROW(pw_type pw({1.,0.},{eid(),eid(),eid()}), inconsistent_piecewise)
    << "Inconsistency (unordered boundaries) not detected";
}

TEST(PiecewiseTest, EvaluateOnSingleInterval)
{
  piecewise<double, eval_shifted<double>> pw({},{eval_shifted(1.)});

  EXPECT_DOUBLE_EQ(pw(-10.),-9.);
  EXPECT_DOUBLE_EQ(pw(-1.),0.);
  EXPECT_DOUBLE_EQ(pw(0.),1.);
  EXPECT_DOUBLE_EQ(pw(1.),2.);
  EXPECT_DOUBLE_EQ(pw(10.),11.);
}

TEST(PiecewiseTest, EvaluateOnTwoIntervals)
{
  piecewise<double, eval_shifted<double>> pw({0.},{eval_shifted(-10.),eval_shifted(10.)});

  EXPECT_DOUBLE_EQ(pw(-10.),-20.);
  EXPECT_DOUBLE_EQ(pw(-1.),-11.);
  EXPECT_DOUBLE_EQ(pw(0.),10.);
  EXPECT_DOUBLE_EQ(pw(1.),11.);
  EXPECT_DOUBLE_EQ(pw(10.),20.);
}

TEST(PiecewiseTest, EvaluateOnThreeIntervals)
{
  piecewise<double, eval_shifted<double>>
    pw({-5.,5.},{eval_shifted(-10.),eval_shifted(0.),eval_shifted(10.)});

  EXPECT_DOUBLE_EQ(pw(-10.),-20.);
  EXPECT_DOUBLE_EQ(pw(-5.),-5);
  EXPECT_DOUBLE_EQ(pw(-1.),-1.);
  EXPECT_DOUBLE_EQ(pw(0.),0.);
  EXPECT_DOUBLE_EQ(pw(1.),1.);
  EXPECT_DOUBLE_EQ(pw(5.),15.);
  EXPECT_DOUBLE_EQ(pw(10.),20.);
}
