/******************************************************************************
 * This file is part of libome                                                *
 * Copyright (C) 2025 Arnd Behring, Kay Schoenwald                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 ******************************************************************************/

#include <gtest/gtest.h>
#include <vector>
#include <ome/rpd_distribution.h>
#include <ome/laurent_polynomial.h>

using namespace ome;

// Test cases for basic functionality

TEST(RPDDistributionTest, HasRegular)
{
  double d = 1.;
  rpd_distribution<double,double,double> empty(std::nullopt, std::nullopt, std::nullopt);
  rpd_distribution<double,double,double> with_reg(std::make_optional(d), std::nullopt, std::nullopt);

  EXPECT_FALSE(empty.has_regular());
  EXPECT_FALSE(empty.get_regular().has_value());

  EXPECT_TRUE(with_reg.has_regular());
  EXPECT_EQ(*with_reg.get_regular(), d);
}

TEST(RPDDistributionTest, HasPlus)
{
  double d = 1.;
  rpd_distribution<double,double,double> empty(std::nullopt, std::nullopt, std::nullopt);
  rpd_distribution<double,double,double> with_plus(std::nullopt, std::make_optional(d), std::nullopt);

  EXPECT_FALSE(empty.has_plus());
  EXPECT_FALSE(empty.get_plus().has_value());

  EXPECT_TRUE(with_plus.has_plus());
  EXPECT_EQ(*with_plus.get_plus(), d);
}

TEST(RPDDistributionTest, HasDelta)
{
  double d = 1.;
  rpd_distribution<double,double,double> empty(std::nullopt, std::nullopt, std::nullopt);
  rpd_distribution<double,double,double> with_delta(std::nullopt, std::nullopt, std::make_optional(d));

  EXPECT_FALSE(empty.has_delta());
  EXPECT_FALSE(empty.get_delta().has_value());

  EXPECT_TRUE(with_delta.has_delta());
  EXPECT_EQ(*with_delta.get_delta(), d);
}

// Test cases for get_coefficient on data type without a view

template<typename T>
class vec
{
  public:
    using coefficient_type = T;
    using coefficient_has_view = std::false_type;
    using coefficient_view_type = void;

    vec(std::vector<coefficient_type> coefficients)
      : coefficients_(coefficients) {};

    const coefficient_type& operator[](size_t pos) const
    {
      return(coefficients_[pos]);
    };

  private:
    std::vector<coefficient_type> coefficients_;
};

TEST(RPDDistributionTest, GetRPDCoefficientOnNonNested)
{
  vec<int> reg({1, 2, 3}),
           plus({11, 12, 13}),
           delta({21, 22, 23});

  rpd_distribution<vec<int>,vec<int>,vec<int>> rpd(reg, plus, delta);

  auto rpd_coeff0 = rpd.get_coefficient(0);
  ASSERT_TRUE(rpd_coeff0.has_regular());
  ASSERT_TRUE(rpd_coeff0.has_plus());
  ASSERT_TRUE(rpd_coeff0.has_delta());
  EXPECT_EQ(*(rpd_coeff0.get_regular()), 1);
  EXPECT_EQ(*(rpd_coeff0.get_plus()), 11);
  EXPECT_EQ(*(rpd_coeff0.get_delta()), 21);
}

TEST(RPDDistributionTest, GetRPDCoefficientOnPartiallyMissing)
{
  vec<int> reg({1, 2, 3}),
           plus({11, 12, 13}),
           delta({21, 22, 23});

  rpd_distribution<vec<int>, vec<int>, vec<int>>
    missing_reg{std::nullopt, plus, delta},
    missing_plus{reg, std::nullopt, delta},
    missing_delta{reg, plus, std::nullopt};

  auto missing_reg_coeff0 = missing_reg.get_coefficient(0);
  EXPECT_FALSE(missing_reg_coeff0.has_regular());
  ASSERT_TRUE (missing_reg_coeff0.has_plus());
  ASSERT_TRUE (missing_reg_coeff0.has_delta());
  EXPECT_FALSE(missing_reg_coeff0.get_regular().has_value());
  EXPECT_EQ(*(missing_reg_coeff0.get_plus()),    11);
  EXPECT_EQ(*(missing_reg_coeff0.get_delta()),   21);

  auto missing_plus_coeff0 = missing_plus.get_coefficient(0);
  ASSERT_TRUE (missing_plus_coeff0.has_regular());
  EXPECT_FALSE(missing_plus_coeff0.has_plus());
  ASSERT_TRUE (missing_plus_coeff0.has_delta());
  EXPECT_EQ(*(missing_plus_coeff0.get_regular()), 1);
  EXPECT_FALSE(missing_plus_coeff0.get_plus().has_value());
  EXPECT_EQ(*(missing_plus_coeff0.get_delta()),   21);

  auto missing_delta_coeff0 = missing_delta.get_coefficient(0);
  ASSERT_TRUE (missing_delta_coeff0.has_regular());
  ASSERT_TRUE (missing_delta_coeff0.has_plus());
  EXPECT_FALSE(missing_delta_coeff0.has_delta());
  EXPECT_EQ(*(missing_delta_coeff0.get_regular()), 1);
  EXPECT_EQ(*(missing_delta_coeff0.get_plus()),    11);
  EXPECT_FALSE(missing_delta_coeff0.get_delta().has_value());
}

TEST(RPDDistributionTest, GetRPDCoefficientOnNested)
{
  vec<vec<int>> reg  ({vec<int>({ 1,  2,  3}), vec<int>({11, 12, 13})}),
                plus ({vec<int>({21, 22, 23}), vec<int>({31, 32, 33})}),
                delta({vec<int>({41, 42, 43}), vec<int>({51, 52, 53})});

  rpd_distribution rpd(std::make_optional(reg),
                       std::make_optional(plus),
                       std::make_optional(delta));

  auto rpd_coeff1 = rpd.get_coefficient(1);
  EXPECT_TRUE(rpd_coeff1.has_regular());
  EXPECT_TRUE(rpd_coeff1.has_plus());
  EXPECT_TRUE(rpd_coeff1.has_delta());

  auto rpd_coeff1_coeff2 = rpd_coeff1.get_coefficient(2);
  ASSERT_TRUE(rpd_coeff1_coeff2.has_regular());
  ASSERT_TRUE(rpd_coeff1_coeff2.has_plus());
  ASSERT_TRUE(rpd_coeff1_coeff2.has_delta());
  EXPECT_EQ(*(rpd_coeff1_coeff2.get_regular()), 13);
  EXPECT_EQ(*(rpd_coeff1_coeff2.get_plus()),    33);
  EXPECT_EQ(*(rpd_coeff1_coeff2.get_delta()),   53);
}

// Test cases for get_coefficient on data type with a view

TEST(RPDDistributionTest, GetRPDCoefficientOnNestedWithView)
{
  using poly = laurent_polynomial<double, double>;
  using ppoly = laurent_polynomial<double, poly, double>;
  ppoly reg  ({poly({1., 2., 3.}), poly({9., 8., 7.})}),
        plus ({poly({4., 5., 6.}), poly({6., 5., 4.})}),
        delta({poly({7., 8., 9.}), poly({3., 2., 1.})});

  rpd_distribution<ppoly, ppoly, ppoly> rpd(reg, plus, delta);

  auto rpd_coeff1 = rpd.get_coefficient(1);
  EXPECT_TRUE(rpd_coeff1.has_regular());
  EXPECT_TRUE(rpd_coeff1.has_plus());
  EXPECT_TRUE(rpd_coeff1.has_delta());

  auto rpd_coeff1_coeff2 = rpd_coeff1.get_coefficient(2);
  ASSERT_TRUE(rpd_coeff1_coeff2.has_regular());
  ASSERT_TRUE(rpd_coeff1_coeff2.has_plus());
  ASSERT_TRUE(rpd_coeff1_coeff2.has_delta());
  EXPECT_EQ(*(rpd_coeff1_coeff2.get_regular()), 7.);
  EXPECT_EQ(*(rpd_coeff1_coeff2.get_plus()),    4.);
  EXPECT_EQ(*(rpd_coeff1_coeff2.get_delta()),   1.);
}

// Test cases for truncate operation

TEST(RPDDistributionTest, Truncate)
{
  using poly = laurent_polynomial<double, double>;
  poly reg  ({ 1.,  2.,  3.}),
       plus ({11., 12., 13.}),
       delta({21., 22., 23.});

  rpd_distribution<poly, poly, poly> rpd(reg, plus, delta);

  auto rpd_trunc0 = rpd.truncate(0);
  ASSERT_TRUE(rpd_trunc0.has_regular());
  ASSERT_TRUE(rpd_trunc0.has_plus());
  ASSERT_TRUE(rpd_trunc0.has_delta());
  EXPECT_EQ((*rpd_trunc0.get_regular())(10.),  1.);
  EXPECT_EQ((*rpd_trunc0.get_plus())(10.),    11.);
  EXPECT_EQ((*rpd_trunc0.get_delta())(10.),   21.);
}
