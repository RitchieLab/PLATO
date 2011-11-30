 // Permutation.h

#ifndef __PERMUTATION_H__
#define __PERMUTATION_H__

#include "Filter.h"
#include <ComboGenerator.h>

///
/// Shuffles dataset statuses and runs permutation tests
/// by establishing distributions that can be used by filters
/// for assigning significance to results
///
namespace Filters{
/// Runs permutations for assigning significance to results
class Permutation{

  public:
  
    /// Constructor
    Permutation();
  
    /// Sets number of permutations to run
    void set_num_permutations(int n_perms);
    
    /// Adds filter to list for running permutation test
    void add_filter(Filter * permuted_filter);
    
    /// Returns p value for indicated value
    float get_p_value(float value, unsigned int filter_number=0);
    
    /// Runs permutations and creates distribution
    void run_permutations(Methods::DataSet & data);
    
    /// Runs indicated number of permutations and creates distribution
    void run_permutations(Methods::DataSet & data, unsigned int n_perms);
    
    /// Runs permutation tests optimized for use with contingency tables
    void run_permutations_contingency(Methods::DataSet & data);

    /// Sets combination sizes to run in generating models to test
    void set_model_sizes(unsigned int min, unsigned int max){min_model_size=min; max_model_size=max;}

    /// Permutes status 
    void permute_status(Methods::DataSet&, vector<unsigned int, std::allocator<unsigned int> >&, unsigned int);

    /// Sets status
    void set_status(vector<bool>& status, Methods::DataSet& data);

void output_set(Methods::DataSet& data, int perm_num);

  private:
    void initialize();
    
    void get_best_scores(Methods::DataSet & data);
    
    vector<Filter *> filters;
    int num_perms;
    vector<vector<float> > p_value_distribution;
    vector<Filter::BestScore> sort_priorities; 
    
    unsigned int min_model_size, max_model_size;

    double elapsed_time;
    
};
}

#endif
