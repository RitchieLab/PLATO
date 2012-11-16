//BioFilter.h

#ifndef __BIOFILTER_H__
#define __BIOFILTER_H__

#include "Filter.h"

#ifdef MYSQLDEFINED
#include "DB.h"
#endif

#include <time.h>

#ifdef MYSQLDEFINED
#include "BioFilterCombinations.h"
#include "BioRecordResults.h"
#endif

namespace Filters{

struct table_info{
  unsigned int combo_size;
  string name;
};

/// Filter extracts SNP combinations from database
class BioFilter: public Filter{

  public:
    BioFilter();
    BioFilter(string filterName);
    void analyze(ResultSet & resultList, Methods::DataSet & dataset);
    void set_params(PARAMS params, Methods::DataSet* set);

    ProcessEstimate estimate_run_time(double num_models, Methods::DataSet& dataset,
      ResultSet& testList, int model_size=1);

  private:
    void initialize();

    void set_text_outname();

    void analyze_list(ResultSet & resultList, Methods::DataSet & dataset);

#ifdef MYSQLDEFINED
    void analyze_from_db(ResultSet & resultList, Methods::DataSet & dataset);
#endif

    void analyze_from_text(ResultSet & resultList, Methods::DataSet & dataset);

    void fill_index(vector<unsigned int>& index_combo, vector<string>& string_combo,
      Methods::DataSet& dataset);

    void fill_list(ResultSet& templist, Methods::DataSet& dataset);

    // add to list below for each additional parameter
    enum ConfigList{
      NoMatch,
      TableName,
      TableOut,
      InteractionIncluded,
      CovariateIncluded,
      OutSchema,
      InSchema,
      MaxModelSize,
      OutTextFile,
      InFileName
    };

    /// performs data extraction
    /// contains names and sizes of views/tables for extracting data
    vector<table_info> bio_tables;

#ifdef MYSQLDEFINED
    BioRecordResults recorder;
#endif

    map<string, ConfigList> ConfigMap;

    string output_table_name, in_schema, out_schema, in_filename, out_text_file;
    bool interact_included, include_cov;
    int max_model_size;

string current_time();

};

}

#endif
