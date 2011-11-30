#include "Helpers.h"
#include "InputFilter.h"
#include "MethodException.h"
namespace Methods{
/*
 * Function: ExcludeLocusFilter
 * Checks if Marker is in list provided.  If so, then return false meaning "don't use this marker".  If not found
 * then return true meaningn "use this marker"
 *
 * return: bool
 */
void InputFilter::ExcludeLocusFilter(vector<Marker*>* marks, vector<Marker*>* mlist){
	stable_sort(mlist->begin(), mlist->end(), greater<Methods::Marker*>());
	vector<Marker*> temp = *marks;
	stable_sort(temp.begin(), temp.end(), greater<Methods::Marker*>());

	unsigned int mi = 0;
	for(unsigned int i = 0; i < temp.size(); i++){
		if(mi >= mlist->size()){
			break;
		}
		while(temp[i]->getRSID() > (*mlist)[mi]->getRSID()){
			mi++;
			if(mi >= mlist->size()){
				break;
			}
		}
		if(mi >= mlist->size()){
			break;
		}
		if(temp[i]->getRSID() == (*mlist)[mi]->getRSID()){
			temp[i]->setEnabled(false);
		}
}

/*
 * Function: IncludeLocusFilter
 * Checks vector of markers against vector of rsids for inclusion (bulk compare instead of single compare
 *
 */
void InputFilter::IncludeLocusFilter(vector<Marker*>* marks, vector<Marker*>* mlist){
	stable_sort(mlist->begin(), mlist->end(), greater<Methods::Marker*>());
	vector<Marker*> temp = *marks;
	stable_sort(temp.begin(), temp.end(), greater<Methods::Marker*>());

	unsigned int mi = 0;
	for(unsigned int i = 0; i < temp.size(); i++){
		if(mi >= mlist->size()){
			temp[i]->setEnabled(false);
			continue;
		}
		while(temp[i]->getRSID() > (*mlist)[mi]->getRSID()){
			mi++;
			if(mi >= mlist->size()){
				break;
			}
		}
		if(mi >= mlist->size()){
			temp[i]->setEnabled(false);
			continue;
		}
		if(temp[i]->getRSID() == (*mlist)[mi]->getRSID()){
			temp[i]->setEnabled(true);
		}else{temp[i]->setEnabled(false);}
	}
}


bool InputFilter::IncludeSampleFilter(Sample* s, vector<Sample*>* list){
	for(unsigned int i = 0; i < list->size(); i++){
		if(s->getInd() == (*list)[i]->getInd() && s->getFamID() == (*list)[i]->getFamID()){
			return true;
		}
	}
	return false;
}

bool InputFilter::ExcludeSampleFilter(Sample* s, vector<Sample*>* list){
	for(unsigned int i = 0; i < list->size(); i++){
		if(s->getInd() == (*list)[i]->getInd() && s->getFamID() == (*list)[i]->getFamID()){
			return false;
		}
	}
	return true;
}

bool InputFilter::IncludeFamilyFilter(Family* f, vector<Family*>* list){
	for(unsigned int i = 0; i < list->size(); i++){
		if(f->getFamID() == (*list)[i]->getFamID()){
			return true;
		}
	}
	return false;
}

bool InputFilter::ExcludeFamilyFilter(Family* f, vector<Family*>* list){
	for(unsigned int i = 0; i < list->size(); i++){
		if(f->getFamID() == (*list)[i]->getFamID()){
			return false;
		}
	}
	return true;
}

void InputFilter::LocusChromFilter(vector<Marker*>* marks, vector<Marker*>* list){
	stable_sort(list->begin(), list->end(), less<Methods::Marker*>());
	vector<Marker*> temp = *marks;
	stable_sort(temp.begin(), temp.end(), less<Methods::Marker*>());

	unsigned int ml = 0;
	for(unsigned int i = 0; i < temp.size(); i++){
		if(ml >= list->size()){
			temp[i]->setEnabled(false);
			continue;
		}
		while((*list)[ml]->getChrom() < temp[i]->getChrom()){
			ml++;
			if(ml >= list->size()){
				break;
			}
		}

		if(ml >= list->size()){
			temp[i]->setEnabled(false);
			continue;
		}
		if((*list)[ml]->getChrom() < temp[i]->getChrom()){
			ml++;
		}
		if(ml >= list->size()){
			temp[i]->setEnabled(false);
			continue;
		}
		if(temp[i]->getChrom() == (*list)[ml]->getChrom()){
			temp[i]->setEnabled(true);
		}
		else{
			temp[i]->setEnabled(false);
		}
	}
}

/*
 * Function LocusBplocRangeFilter
 * First element of vector is assumed to be MIN value.  Second element is assumed to be the MAX value
 * If one or the other are NULL, then will use non-null value as min or max.
 */
void InputFilter::LocusBplocRangeFilter(vector<Marker*>* marks, vector<Marker*>* list){
	if(list == NULL){
		return;
	}
	if(list->size() != 2){
		throw MethodException("Only 2 elements may be passed in the vector to LocusBplocRangeFilter...\n");
	}
	Marker* min = (*list)[0];
	Marker* max = (*list)[1];

	vector<Marker*> temp = *marks;
	stable_sort(temp.begin(), temp.end(), less<Methods::Marker*>());

	for(unsigned int i = 0; i < temp.size(); i++){
		if(min != NULL){
			if(temp[i]->getBPLOC() < min->getBPLOC()){
				temp[i]->setEnabled(false);
			}
		}
		if(max != NULL){
			if(temp[i]->getBPLOC() > max->getBPLOC()){
				temp[i]->setEnabled(false);
			}
		}
	}
}

bool InputFilter::IncludeCovariateFilter(string cov, vector<string>* list){
	if(list == NULL){
		return true;
	}
	for(unsigned int i = 0; i < list->size(); i++){
		if(cov == (*list)[i]){
			return true;
		}
	}
	return false;
}

bool InputFilter::ExcludeCovariateFilter(string cov, vector<string>* list){
	if(list == NULL){
		return true;
	}
	for(unsigned int i = 0; i < list->size(); i++){
		if(cov == (*list)[i]){
			return false;
		}
	}
	return true;
}

bool InputFilter::IncludeTraitFilter(string trait, vector<string>* list){
	if(list == NULL){
		return true;
	}
	for(unsigned int i = 0; i < list->size(); i++){
		if(trait == (*list)[i]){
			return true;
		}
	}
	return false;
}

bool InputFilter::ExcludeTraitFilter(string trait, vector<string>* list){
	if(list == NULL){
		return true;
	}
	for(unsigned int i = 0; i < list->size(); i++){
		if(trait == (*list)[i]){
			return false;
		}
	}
	return true;
}

void InputFilter::run_locus_filter(int f, vector<Marker*>* marks){
	marker_filters[f](marks, marker_lists[f]);
}

bool InputFilter::run_sample_filter(int f, Sample* s){
	return sample_filters[f](s, sample_lists[f]);
}

bool InputFilter::run_family_filter(int f, Family* fam){
	return family_filters[f](fam, family_lists[f]);
}

bool InputFilter::run_covariate_filter(int f, string s){
	return cov_filters[f](s, cov_lists[f]);
}

bool InputFilter::run_trait_filter(int f, string s){
	return trait_filters[f](s, trait_lists[f]);
}

void InputFilter::run_sample_bprange_filter(vector<Sample*>* samples, vector<Marker*>* markers, vector<Sample*> filter_samples,
		vector<vector<Marker*> > filter_ranges){

	string fname = opts::_OUTPREFIX_ + "sample_bprange_filter_log.txt";
	ofstream eout(fname.c_str());
	if (!eout) {
		opts::printLog("Error opening " + fname + "!  Exiting!\n");
		throw MethodException("");
	}
	eout << "FID IID Chr bp SNP" << endl;

	int ssize = (int)filter_samples.size();
	int msize = (int)markers->size();

	for(int s = 0; s < ssize; s++){
		vector<Sample*>::iterator s_iter;
		s_iter = find_if(samples->begin(), samples->end(), FindSampleByFamAndID(filter_samples[s]->getFamID(), filter_samples[s]->getInd()));
		if(s_iter != samples->end()){
			Sample* found = (*s_iter);
			vector<Marker*> range = filter_ranges[s];
			Marker* min = range[0];
			Marker* max = range[1];

			for(int m = 0; m < msize; m++){
				Marker* mark = (*markers)[m];
				if(mark->getChrom() == min->getChrom() && mark->getBPLOC() >= min->getBPLOC() && mark->getBPLOC() <= max->getBPLOC()){
					int loc = mark->getLoc();
					found->addAone(loc, true);
					found->addAtwo(loc, true);
					found->addAmissing(loc, true);
					eout << found->getFamID() << " " << found->getInd() << " " << mark->getChrom() << " " << mark->getBPLOC() << " " << mark->getRSID() << endl;
				}
			}
		}
	}

	eout.close();
}

}
