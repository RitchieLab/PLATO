#include <Marker.h>
#include <gtest/gtest.h>

namespace {

class MarkerTest : public ::testing::Test {
 protected:

   virtual void SetUp() {
     m2_.addAllele("x"); m2_.addAllele("y");
     // This is testing that alleles of size 3 are micro satellites
     m3_.addAllele("x"); m3_.addAllele("y"); m3_.addAllele("z");
   }

  Methods::Marker m1_;
  Methods::Marker m2_;
  Methods::Marker m3_;
};

TEST_F(MarkerTest, Initialize) {
  EXPECT_STREQ("", m1_.getProbeID().c_str()) << "probe_id not initialized";
  EXPECT_STREQ("", m1_.getRSID().c_str()) << "rsid not initialized";
  EXPECT_STREQ("", m1_.getEnzyme().c_str()) << "enzyme not initialized";
  EXPECT_STREQ("", m1_.getAllele1().c_str()) << "allele1 not initialized";
  EXPECT_STREQ("", m1_.getAllele2().c_str()) << "allele2 not initialized";
  EXPECT_STREQ("", m1_.getReferent().c_str()) << "referent_allele not initialized";
  
  EXPECT_EQ(-1, m1_.chrom) << "chrom not initialized";
  EXPECT_EQ(-1, m1_.bploc) << "bploc not initialized";
  EXPECT_EQ(-1, m1_.getLoc()) << "loc not initialized";
  EXPECT_EQ(-1, m1_.getMAF()) << "maf not initialized";
  EXPECT_EQ(-1, m1_.getReferentIndex()) << "ref_all_index not initialized";
  
  EXPECT_FALSE(m1_.isFlagged()) << "flag wasn't initialized correctly";
  EXPECT_FALSE(m1_.isEnabled()) << "enabled wasn't initialized correctly";
  EXPECT_FALSE(m1_.hasMAF()) << "freqflag wasn't initialized correctly";
  
}

TEST_F(MarkerTest, getAlleleLoc) {
  EXPECT_EQ(-1, m1_.getAlleleLoc("x"));
  EXPECT_EQ(0, m2_.getAlleleLoc("x"));
  EXPECT_EQ(1, m2_.getAlleleLoc("y"));
}

TEST_F(MarkerTest, isMicroSat) {
  EXPECT_FALSE(m1_.isMicroSat());
  EXPECT_TRUE(m3_.isMicroSat());
  
  // This is a strange test, but it is easy to see the error from the code
  m1_.setAllele1("x"); m1_.setAllele1("y");
  m1_.addAllele("y");
  EXPECT_FALSE(m1_.isMicroSat());
}

TEST_F(MarkerTest, setAllele1) {
  m1_.setAllele1("x");
  EXPECT_EQ("x", m1_.getAllele1());
  m1_.setAllele1("y");
  EXPECT_EQ("y", m1_.getAllele1());
}

TEST_F(MarkerTest, setAllele2) {
  m1_.setAllele2("x");
  EXPECT_EQ("x", m1_.getAllele2());
  m1_.setAllele2("y");
  EXPECT_EQ("y", m1_.getAllele2());
}

TEST_F(MarkerTest, getAllele) {
  m1_.setAllele1("x");
  EXPECT_EQ("x", m1_.getAllele(0));
  m1_.setAllele2("y");
  EXPECT_EQ("y", m1_.getAllele(1));
  m1_.setAllele2("z");
  EXPECT_EQ("z", m1_.getAllele(1));
  // EXPECT_ANY_THROW(m1_.getAllele(3));
}

} // namespace