//DataSet.h

#ifndef __DATASET_H__
#define __DATASET_H__
#include "Family.h"
#include "Sample.h"
#include "Marker.h"
#include "MethodException.h"

namespace Methods{
///
/// Stores dataset information.
/// Does not free memory allocated to the pointers stored in
/// the families, samples and markers vectors.
///


class DataSet{

  public:
    DataSet();
    DataSet(vector<Sample*>* samps, vector<Family*>* fams, vector<Marker*>* marks,
      vector<int>* mark_map);
    ~DataSet(){};
    /// Add data information to dataset
    void add_info(vector<Sample*>* samps, vector<Family*>* fams, vector<Marker*>* marks,
      vector<int>* mark_map);

    /// Set samples vector
    void set_samples(vector<Sample*>* samps){samples = *samps;}
    /// Set families vector
    void set_families(vector<Family*>* fams){families = *fams;}
    /// Set markers vector
    void set_markers(vector<Marker*>* loci){markers = *loci;}
    /// Sets the marker map
    void set_marker_map(vector<int>* mark_map){marker_map = *mark_map;}

    /// Get markers vector
    inline vector<Marker*>* get_markers(){return &markers;}
    /// Get marker map
    inline vector<int>* get_marker_map(){return &marker_map;}
    /// Get samples vector
    inline vector<Sample*>* get_samples(){return &samples;}
    /// Get family vector
    inline vector<Family*>* get_families(){return &families;}

	///Get loci by index
	Marker* get_locus(unsigned int index){return markers[index];}

	  ///Get sample by index
	Sample* get_sample(unsigned int index){return samples[index];}

    /// Overloaded to return sample at indicated index
    inline Sample* operator[](int indIndex){return samples[indIndex];}
    /// Total number of individuals (samples) in set
    int num_inds(){return samples.size();}

    /// Total number of enabled individuals
    int num_enabled_inds();

	///look for snp index
	int find_snp_index_by_name(string name);

    /// Number of affected individuals in dataset
    unsigned int num_affected(){return affected_inds.size();}
    /// Number of unaffected individuals in dataset
    unsigned int num_unaffected(){return unaffected_inds.size();}
	/// Number of Males in dataset
	unsigned int num_males();
	/// Number of females in dataset
	unsigned int num_females();
    /// Number of loci (markers) in set
    unsigned int num_loci(){return markers.size();}
    /// Number of families in set
    unsigned int num_pedigrees(){return families.size();}
    /// Returns pedigree identified by index
    Family* get_pedigree(unsigned int index){return families[index];}

    /// clears individuals from set
    void clear();

	/// clears all data from set
	void clear_all();

	/// converts 0,1,2,3 genotype to string representation (A_A, A_C, etc.)
	vector<string> convert_geno_tostring(vector<unsigned int> genos, vector<unsigned int> loci, vector<unsigned int> covs, vector<unsigned int> traits);

    /// adds individual to set
    void add_ind(Sample* ind);

    /// sets affected and unaffected vectors
    void set_affection_vectors();

    /// value that codes a missing genotype in the set
    unsigned int get_missing_value(){return missing_value;}
    void set_missing_value(unsigned int value){missing_value = value;}

    /// value that codes as missing for traits and covariates
    double get_missing_covalue(){return missing_covalue;}
    void set_missing_covalues(double value){missing_covalue = value;}

    /// returns true when any genotypes missing in set
    bool missing_data_present(){return any_missing;}
    void missing_data_present(bool any){any_missing = any;}

    /// maximum locus value in set
    unsigned int get_max_locus(){return max_locus_value;}
    void set_max_locus(unsigned int max){max_locus_value = max;}

    /// maximum allele value in set
    unsigned int get_max_allele(){return max_allele;}
    void set_max_allele(unsigned int max){max_allele = max;}

    /// get affected individual from affected individuals vector
    Sample* get_affected(unsigned int index){return affected_inds[index];}
    /// get unaffected individual from unaffected individuals vector
    Sample* get_unaffected(unsigned int index){return unaffected_inds[index];}

    /// returns affected samples in vector
    vector<Sample*> get_affected_vector(){return affected_inds;}
    /// returns unaffected samples in vector
    vector<Sample*> get_unaffected_vector(){return unaffected_inds;}

    ///return number of enabled affected
    int num_enabled_affected();
    ///return number of enabled unaffected
    int num_enabled_unaffected();

	/// returns covariates in vector
	inline vector<string>* get_covariates(){return &covariates;}
	inline void set_covariates(vector<string>* c){covariates = *c;}
	inline void add_covariate(string s){covariates.push_back(s);}
	inline vector<int>* get_covariate_map(){return &cov_map;}
	inline void set_covariate_map(vector<int>* m){cov_map = *m;}
	inline int num_covariates(){return covariates.size();}
	/// returns traits in vector
	inline vector<string>* get_traits(){return &traits;}
	inline void set_traits(vector<string>* t){traits = *t;}
	inline void add_trait(string s){traits.push_back(s);}
	inline vector<int>* get_trait_map(){return &trait_map;}
	inline void set_trait_map(vector<int>* m){trait_map = *m;}
	inline int num_traits(){return traits.size();}
	/// returns covariate by index
	string get_covariate_name(int i){
		if(i < (int)covariates.size()){
			return covariates[i];
		}
		else{
			return NULL;
		}
	}
	/// returns trait by index
	string get_trait_name(int i){
		if(i < (int)traits.size()){
			return traits[i];
		}
		else{
			return NULL;
		}
	}

	/// returns index of covariate name
	int get_covariate_index(string v);
	/// returns index of trait name
	int get_trait_index(string v);

	///returns vector of covariate values from samples for a specific index
	vector<double> get_covariates_by_index(int i);

	///returns vector of trait values from samples for a specific index
	vector<double> get_traits_by_index(int i);

    /// clears sample vector
    void clear_samples(){samples.clear();}

    /// shuffles samples
    void shuffle_inds();

    /// Returns index of marker in Marker vector -- throws MethodException when not found
    int get_locus_index(string name){
      if(marker_name_map.find(name) == marker_name_map.end())
        throw MethodException("No locus found matching " + name);
      return marker_name_map[name];
    }


    /// Sets up map with marker name as key and Marker* as value
    void create_marker_name_map(){
      for(unsigned int m=0; m<markers.size(); m++){
        marker_name_map[markers[m]->getRSID()] = m;
      }
    }

	map<int, string>* get_master_map(){return &master_map;}

	///creates the family vector if you only pass samples & markers and don't have family info available
	///creates temp families based on sample->famid
	void recreate_family_vector();

	bool get_binary_trait(){return binary_trait;}
	void set_binary_trait(bool b){binary_trait = b;}

	bool get_quant_trait(){return quant_trait;}
	void set_quant_trait(bool b){quant_trait = b;}

	DataSet* get_aff_unaff_dataset();

	vector<DataSet*> generate_case_control_subsets(int);

  private:

    void initialize();

    void check_for_missing();



    vector<Sample*> samples;
    vector<Family*> families;
    vector<Marker*> markers;
    vector<int> marker_map;
	vector<string> covariates, traits;
	vector<CovTrait*> covariates2, traits2;
	vector<int> cov_map; //used only for mapping location of covariate in ped file during reading
	vector<int> trait_map; //used only for mapping locatino of trait in ped file during reading
	map<int, string> master_map;
	vector<bool> orig_marker_good_bad;
	vector<bool> orig_sample_good_bad;
	vector<bool> orig_family_good_bad;

	bool binary_trait;
	bool quant_trait;
	bool alternate_pheno;
	int pheno_loc;

    map<string, int> marker_name_map;

    vector<Sample*> affected_inds, unaffected_inds;
    unsigned int missing_value, max_locus_value, max_allele;
    bool any_missing, matched_pairs; // true when affected and unaffected are matched pairs
    double missing_covalue;

};
};

#endif
