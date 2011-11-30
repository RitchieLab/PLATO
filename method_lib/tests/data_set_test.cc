#include <limits.h>
#include <DataSet.h>
#include <gtest/gtest.h>
using namespace Methods;

TEST(DataSetTest, Initialize) {
  DataSet ds;
  EXPECT_EQ((unsigned int)3, ds.get_missing_value()) << "Missing value was not initialized to 3";
  EXPECT_EQ((unsigned int)2, ds.get_max_locus()) << "Max locus value was not initialized to 2";
  EXPECT_EQ((unsigned int)1, ds.get_max_allele()) << "Max allele was not initialized to 2";
  EXPECT_TRUE(ds.missing_data_present()) << "any_missing was not initialized to 'true'";
  EXPECT_EQ(-99999, ds.get_missing_covalue()) << "Missing covalue was not initialized to -99999";
  // EXPECT_FALSE(ds.matched_pairs) << "matched_pairs was not initialized to 'false'";
  //EXPECT_EQ(-1, ds.pheno_loc) << "Pheno loc was not initialized to -1";
}

TEST(DataSetTest, FindSnpIndexByName) {
  DataSet ds;
  vector<Marker*> markers;
  Marker* marker = new Marker();
  marker->setRSID("test");
  marker->setChrom(1);
  marker->setBPLOC(2);
  markers.push_back(marker);
  ds.set_markers(&markers);
  EXPECT_EQ(-1, ds.find_snp_index_by_name("bogus")) << "File was found when it doesn't exist";
  EXPECT_EQ(0, ds.find_snp_index_by_name("test")) << "Locus index was found when it doesn't exist";
}
