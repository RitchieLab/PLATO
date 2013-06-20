//DataSet.cpp
#include "Helpers.h"
#include "DataSet.h"
namespace Methods{

DataSet::~DataSet(){
	if(recreated_fams){
		for(vector<Family*>::iterator f_iter=families.begin(); f_iter != families.end(); ++f_iter){
			delete *f_iter;
		}
	}
}

///
/// returns number of enabled samples
///
int DataSet::num_enabled_inds(){
	int count = 0;
	for(int i = 0; i < (int)samples.size(); i++){
		if(samples.at(i)->isEnabled()){
			count++;
		}
	}
	return count;
}


///
/// returns number of enabled affected
///
int DataSet::num_enabled_affected(){
	int count = 0;
	for(int i = 0; i < (int)affected_inds.size(); i++){
		if(affected_inds.at(i)->isEnabled()){
			count++;
		}
	}
	return count;
}

///
/// returns number of enabled unaffected
///
int DataSet::num_enabled_unaffected(){
	int count = 0;
	for(int i = 0; i < (int)unaffected_inds.size(); i++){
		if(unaffected_inds.at(i)->isEnabled()){
			count++;
		}
	}
	return count;
}

///
/// clears all vectors
///
void DataSet::clear_all(){
	markers.clear();
	samples.clear();
	families.clear();
	marker_map.clear();
	covariates.clear();
	traits.clear();
	affected_inds.clear();
	unaffected_inds.clear();
	cov_map.clear();
	trait_map.clear();
	master_map.clear();
	marker_name_map.clear();
}

vector<DataSet*> DataSet::generate_case_control_subsets(int size){
	vector<DataSet*> subsets;

	int cases = this->num_affected();
	int controls = this->num_unaffected();

	vector<Sample*> case_pool;
	vector<Sample*> control_pool;
	for(int i = 0; i < cases; i++){
		if(affected_inds.at(i)->isEnabled()){
			case_pool.push_back(affected_inds.at(i));
		}
	}
	for(int i = 0; i < controls; i++){
		if(unaffected_inds.at(i)->isEnabled()){
			control_pool.push_back(unaffected_inds.at(i));
		}
	}
	unsigned int min_cases = case_pool.size() / 2;
	unsigned int min_controls = control_pool.size() / 2;

	if(min_cases == 0 && min_controls == 0){
		throw MethodException("Cannot generate case/control subsets since no cases or controls are defined!\n");
	}

	DataSet* test_set = new DataSet();
	DataSet* train_set = new DataSet();
	test_set->set_markers(this->get_markers());
	test_set->set_covariates(this->get_covariates());
	train_set->set_markers(this->get_markers());
	train_set->set_covariates(this->get_covariates());

	vector<Sample*> train_samples;
	vector<Sample*> test_samples;
	vector<Sample*> train_control_samples;
	vector<Sample*> test_control_samples;


	vector<int> used_samps;
	while(used_samps.size() != case_pool.size()){
		int random = int(rand() % case_pool.size());
		vector<int>::iterator used = find(used_samps.begin(), used_samps.end(), random);
		if(used == used_samps.end()){
			used_samps.push_back(random);
			if(train_samples.size() >= min_cases){
				test_samples.push_back(case_pool.at(random));
				for(int s = 0; s < (int)case_pool.size(); s++){
					used = find(used_samps.begin(), used_samps.end(), s);
					if(used == used_samps.end()){
						test_samples.push_back(case_pool.at(s));
					}
				}
				break;
			}
			else{
				train_samples.push_back(case_pool.at(random));
			}
		}
	}
	used_samps.clear();
	while(used_samps.size() != control_pool.size()){
		int random = int(rand() % case_pool.size());
		vector<int>::iterator used = find(used_samps.begin(), used_samps.end(), random);
		if(used == used_samps.end()){
			used_samps.push_back(random);
			if(train_control_samples.size() >= min_controls){
				test_control_samples.push_back(control_pool.at(random));
				for(int s = 0; s < (int)control_pool.size(); s++){
					used = find(used_samps.begin(), used_samps.end(), s);
					if(used == used_samps.end()){
						test_control_samples.push_back(control_pool.at(s));
					}
				}
				break;
			}
			else{
				train_control_samples.push_back(control_pool.at(random));
			}
		}
	}

	train_samples.insert(train_samples.end(), train_control_samples.begin(), train_control_samples.end());
	test_samples.insert(test_samples.end(), test_control_samples.begin(), test_control_samples.end());

	train_set->set_samples(&train_samples);
	test_set->set_samples(&test_samples);
	test_set->recreate_family_vector();
	train_set->recreate_family_vector();

	test_set->set_affection_vectors();
	train_set->set_affection_vectors();

	subsets.push_back(train_set);
	subsets.push_back(test_set);

	return subsets;
}

///
///returns copy of dataset object that has only the affecteds and unaffecteds
///
DataSet* DataSet::get_aff_unaff_dataset(){
	DataSet* newset = new DataSet();

    newset->set_markers(this->get_markers());
    newset->set_covariates(this->get_covariates());
    newset->set_traits(this->get_traits());
	for(int s = 0; s < this->num_inds(); s++){
		Sample* samp = this->get_sample(s);
		if(samp->isEnabled() && (samp->getPheno() == 1 || samp->getPheno() == 2)){
			newset->add_ind(samp);
		}
	}
	newset->recreate_family_vector();
	newset->set_affection_vectors();

	return newset;
}

///
///creates the family vector if it is not available to the dataset.
///uses the sample->famid to create families and organize the samples into families.
///
void DataSet::recreate_family_vector(){
	if(families.size() > 0){
		families.clear();
	}
	recreated_fams=true;
	vector<Family*>::iterator f_iter;
	for(unsigned int i = 0; i < samples.size(); i++){
		Sample* samp = samples.at(i);
		f_iter = find_if(families.begin(), families.end(), FindFamily(samp->getFamID()));
		if(f_iter == families.end()){
			Family* fam = new Family();
			fam->setFamID(samp->getFamID());
			fam->setEnabled(true);
			fam->AddInd(samp);
			families.push_back(fam);
			fam->setLoc((families.size() - 1));
		}
		else{
			(*f_iter)->AddInd(samp);
		}
	}
}

///
///
///
int DataSet::find_snp_index_by_name(string name){
	int found = -1;
	for(unsigned int i = 0; i < markers.size(); i++){
		if(markers.at(i)->getRSID() == name){
			return i;
		}
	}
	return found;
}

///
/// Sets affection and unaffection vectors
/// to hold split samples into two categories
/// that can be used by algorithms
///
void DataSet::set_affection_vectors(){
  affected_inds.clear();
  unaffected_inds.clear();

  unsigned int curr_ind;
  unsigned int num_samples = num_inds();
  for(curr_ind = 0; curr_ind < num_samples; curr_ind++){
    if(samples.at(curr_ind)->getAffected())
      affected_inds.push_back(samples.at(curr_ind));
    else
      unaffected_inds.push_back(samples.at(curr_ind));
  }

}

///
/// Constructor
///
DataSet::DataSet(){
  initialize();
}

///
/// Constructor that adds classes containing information for the set
/// @param samps vector<Sample*> *
/// @param fams
/// @param marks
/// @param mark_map
///
DataSet::DataSet(vector<Sample*>* samps, vector<Family*>* fams, vector<Marker*>* marks,
  vector<int>* mark_map){
  initialize();
  add_info(samps, fams, marks, mark_map);
}


///
/// Set initial state of set
///
void DataSet::initialize(){
  missing_value = 3;
  max_locus_value = 2;
  max_allele = 1;
  any_missing = true;
  matched_pairs = false;
  missing_covalue = -99999;
  alternate_pheno = false;
  pheno_loc = -1;
  recreated_fams=false;
}

///
/// Checks for presence of missing data in set
///
void DataSet::check_for_missing(){
  // check for missing data
  // if found, set missing in set to true and exit function
  unsigned int total_markers = num_loci();
  unsigned int total_inds = num_inds();

  unsigned int curr_marker, curr_ind;

  for(curr_ind=0; curr_ind < total_inds; curr_ind++){
    for(curr_marker=0; curr_marker < total_markers; curr_marker++){
      if(samples.at(curr_ind)->get_genotype(curr_marker) == missing_value){
        any_missing = true;
        return;
      }
    }
  }

  // if no missing have been found, set appropriate variable here
  any_missing = false;

}

///
/// returns a vector of covariates of sample values for a given index
///
vector<double> DataSet::get_covariates_by_index(int i){
	vector<double> values;
	for(unsigned int s = 0; s < samples.size(); s++){
		values.push_back(samples.at(s)->getCovariate(i));
	}
	return values;
}

///
/// returns a vector of traits of sample values for a given index
///
vector<double> DataSet::get_traits_by_index(int i){
	vector<double> values;
	for(unsigned int s = 0; s < samples.size(); s++){
		values.push_back(samples.at(s)->getTrait(i));
	}
	return values;
}

///
/// Converts numeric genotype sum (Sample.get_genotype) to string
/// based on alleles in the Loci being looked at.
/// Returns vector of genotypes being requested.
///
vector<string> DataSet::convert_geno_tostring(vector<unsigned int> genos, vector<unsigned int> loci, vector<unsigned int> covs, vector<unsigned int> traits){
	vector<string> strgenos(genos.size());

	for(unsigned int g = 0; g < genos.size(); g++){
		if(g < loci.size()){
			Marker* mark = markers[loci.at(g)];
			if(genos.at(g) == 0){
				strgenos.at(g) = mark->getAllele1() + "_" + mark->getAllele1();
			}
			else if(genos.at(g) == 1){
				strgenos.at(g) = mark->getAllele1() + "_" + mark->getAllele2();
			}
			else if(genos.at(g) == 2){
				strgenos.at(g) = mark->getAllele2() + "_" + mark->getAllele2();
			}
			else{
				strgenos.at(g) = "0_0";
			}
		}
		else{
			strgenos.at(g) = getString<int>(genos.at(g));//traits[traitcount++];
		}
	}
	return strgenos;
}

///
/// Constructor that adds classes containing information for the set
/// @param samps vector<Sample*> *  sample information
/// @param fams vector<Family*> * family information
/// @param marks vector<Marker*> * marker information
/// @param mark_map vector<int> * marker map
///
void DataSet::add_info(vector<Sample*>* samps, vector<Family*>* fams, vector<Marker*>* marks,
  vector<int>* mark_map){
    set_samples(samps);
    set_families(fams);
    set_markers(marks);
    set_marker_map(mark_map);
    set_affection_vectors();

    check_for_missing();

}


///
/// Shuffles individuals. <br>
/// When individual data are not matched pairs the
/// affected and unaffected individuals are sorted
/// independently <br>
/// When individual data are matched pairs the affected
/// and unaffected individuals maintain their relative
/// positions in the set
/// @return none
///
void DataSet::shuffle_inds(){
  // when matched data shuffle both affected
  // and unaffected at same time
  int randIndex, lastActive;
  unsigned int currSize, maxSize, currPos;
  Sample* tempHolder;

  if(matched_pairs){
    // when matched every other individual is affected (or unaffected)
    // starting ind can be either
    currSize=affected_inds.size();
    maxSize = currSize;
    for(currPos=0; currPos < maxSize;
      currPos++){

      // randomly select one pair
      randIndex = rand() % currSize;
      lastActive = currSize - 1;

      // swap last affected and unaffected
      tempHolder = affected_inds.at(lastActive);
      affected_inds.at(lastActive) = affected_inds.at(randIndex);
      affected_inds.at(randIndex) = tempHolder;

      tempHolder = unaffected_inds.at(lastActive);
      unaffected_inds.at(lastActive) = unaffected_inds.at(randIndex);
      unaffected_inds.at(randIndex) = tempHolder;

      currSize--;
    }
  }
  // shuffle affected and unaffected independently
  else{
    affected_inds.clear();
    unaffected_inds.clear();
    currSize=num_inds();
    maxSize = currSize;
    for(currPos=0; currPos < maxSize;
      currPos++){

      // randomly select one individual
      randIndex = rand() % currSize;
      lastActive = currSize-1;
      // swap this individual and one after with last positions

      // add to unaffected and affected
      if(samples.at(randIndex)->getAffected()){
        affected_inds.push_back(samples.at(randIndex));
      }
      else{
        unaffected_inds.push_back(samples.at(randIndex));
      }

      tempHolder = samples.at(randIndex);
      samples.at(randIndex) = samples.at(lastActive);
      samples.at(lastActive) = tempHolder;

      currSize--;
    }
  }

}


///
/// Clears all indivduals from the dataset
/// @return none
///
void DataSet::clear(){
  samples.clear();
  affected_inds.clear();
  unaffected_inds.clear();
}


///
/// Adds the individual <br>
/// Function can be used to create a set using
/// only a portion of the individuals in original
/// set
/// @param  ind  Sample* individual to add
/// @return none
///
void DataSet::add_ind(Sample* ind){
  samples.push_back(ind);
  if(ind->getAffected()){
    affected_inds.push_back(ind);
  }
  else{
    unaffected_inds.push_back(ind);
  }
}

///
/// Returns the Index of a covariate by name
/// @param v string value to find
/// @return int index
int DataSet::get_covariate_index(string v){
	for(int i = 0; i < (int)covariates.size(); i++){
		if(v == covariates.at(i)){
			return i;
		}
	}
	return -1;
}

///
/// Returns the Index of a trait by name
/// @param v string value to find
/// @return int index
int DataSet::get_trait_index(string v){
	for(int i = 0; i < (int)traits.size(); i++){
		if(v == traits.at(i)){
			return i;
		}
	}
	return -1;
}

///
/// Returns number of males in dataset
///
unsigned int DataSet::num_males(){
	unsigned int count = 0;
	for(unsigned int i = 0; i < samples.size(); i++){
		if(samples.at(i)->getSex()){
			count++;
		}
	}
	return count;
}

///
/// Returns number of females in dataset
///
unsigned int DataSet::num_females(){
	unsigned int count = 0;
	for(unsigned int i = 0; i < samples.size(); i++){
		if(!samples.at(i)->getSex()){
			count++;
		}
	}
	return count;
}

}
