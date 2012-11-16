#include <limits.h>
#include <Helper.h>
#include <gtest/gtest.h>
#include <math.h>

TEST(HelperTest, Square) {
  double x = 5;
  EXPECT_EQ(pow(x, 2), std::SQR(x)) << "int var wasn't squared correctly";
}

TEST(HelperTest, Comparisons) {
  int x = 5;
  int y = 6;
  EXPECT_EQ(y, std::MAX(x, y)) << "smaller value was returned incorectly";
  EXPECT_EQ(x, std::MIN(x, y)) << "larger value was returned incorectly";
  EXPECT_EQ(x, Abs(x)) << "positive absolute value returned negative result";
  EXPECT_EQ(x, Abs(-x)) << "negative absolute value returned negative result";
  EXPECT_EQ(0.00000000000001, TOLERANCE) << "TOLERANCE set to 0.00000000000001";
  EXPECT_EQ(0.00000001, FTOLERANCE) << "FTOLERANCE set to 0.00000001";
}

TEST(HelperTest, Calculations) {
  double x = 5.0;
  double y = 6.0;
  EXPECT_EQ((Abs(x - y) / y), Methods::DoubDif(x, y)) << "DoubDif returned the wrong answer for positive values";
  EXPECT_EQ((float)(Abs(x - y) / y), Methods::FloatDif((float)x, (float)y)) << "FloatDif returned the wrong answer for positive values";
  EXPECT_EQ((Abs(-x - y) / y), Methods::DoubDif(-x, y)) << "DoubDif returned the wrong answer for one negative value";
  EXPECT_EQ((float)(Abs(-x - y) / y), Methods::FloatDif((float)-x, (float)y)) << "FloatDif returned the wrong answer for one negative value";
  EXPECT_EQ((Abs(-x + y) / y), Methods::DoubDif(-x, -y)) << "DoubDif returned the wrong answer for two negative value";
  EXPECT_EQ((float)(Abs(-x + y) / y), Methods::FloatDif((float)-x, (float)-y)) << "FloatDif returned the wrong answer for two negative value";

  EXPECT_TRUE(Methods::float_comp((float)y, (float)y));
  EXPECT_TRUE(Methods::fEquals((float)y, (float)y));
  EXPECT_FALSE(Methods::float_comp((float)x, (float)y));
  EXPECT_FALSE(Methods::fEquals((float)x, (float)y));

  EXPECT_TRUE(Methods::double_comp(y, y));
  EXPECT_TRUE(Methods::dEquals(y, y));
  EXPECT_FALSE(Methods::double_comp(x, y));
  EXPECT_FALSE(Methods::dEquals(x, y));
}