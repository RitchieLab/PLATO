#include "Helpers.h"
#include "InputFilter.h"
#include "MethodException.h"
namespace Methods{
/*
 * Function: ExcludeLocusFilter
 * Checks if Marker is in list provided.  If so, then return false meaning "don't use this marker".  If not found
 * then return true meaning "use this marker"
 *
 * return: bool
 */
void InputFilter::ExcludeLocusFilter(vector<Marker*>* marks, vector<Marker*>* mlist){
	//stable_sort(mlist->begin(), mlist->end(), greater<Methods::Marker*>());
	stable_sort(mlist->begin(), mlist->end(), Helpers::markerGreater);
	vector<Marker*> temp = *marks;
	//stable_sort(temp.begin(), temp.end(), greater<Methods::Marker*>());
	stable_sort(temp.begin(), temp.end(), Helpers::markerGreater);

	unsigned int mi = 0;
	for(unsigned int i = 0; i < temp.size(); i++){
		if(mi >= mlist->size()){
			break;
		}
		while(temp.at(i)->getRSID() > (*mlist).at(mi)->getRSID()){
			i++;
			if(mi >= mlist->size()){
				break;
			}
		}
		if(mi >= mlist->size()){
			break;
		}
		if(temp.at(i)->getRSID() == (*mlist).at(mi)->getRSID()){
			temp.at(i)->setEnabled(false);
			mi++;
		}
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
			temp.at(i)->setEnabled(false);
			continue;
		}
		while(temp.at(i)->getRSID() > (*mlist).at(mi)->getRSID()){
			mi++;
			if(mi >= mlist->size()){
				break;
			}
		}
		if(mi >= mlist->size()){
			temp.at(i)->setEnabled(false);
			continue;
		}
		if(temp.at(i)->getRSID() == (*mlist).at(mi)->getRSID()){
			temp.at(i)->setEnabled(true);
		}else{temp.at(i)->setEnabled(false);}
	}
}


bool InputFilter::IncludeSampleFilter(Sample* s, vector<Sample*>* list){
	for(unsigned int i = 0; i < list->size(); i++){
		if(s->getInd() == (*list).at(i)->getInd() && s->getFamID() == (*list).at(i)->getFamID()){
			return true;
		}
	}
	return false;
}

bool InputFilter::ExcludeSampleFilter(Sample* s, vector<Sample*>* list){
	for(unsigned int i = 0; i < list->size(); i++){
		if(s->getInd() == (*list).at(i)->getInd() && s->getFamID() == (*list).at(i)->getFamID()){
			return false;
		}
	}
	return true;
}

bool InputFilter::IncludeFamilyFilter(Family* f, vector<Family*>* list){
	for(unsigned int i = 0; i < list->size(); i++){
		if(f->getFamID() == (*list).at(i)->getFamID()){
			return true;
		}
	}
	return false;
}

bool InputFilter::ExcludeFamilyFilter(Family* f, vector<Family*>* list){
	for(unsigned int i = 0; i < list->size(); i++){
		if(f->getFamID() == (*list).at(i)->getFamID()){
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
			temp.at(i)->setEnabled(false);
			continue;
		}
		while((*list).at(ml)->getChrom() < temp[i]->getChrom()){
			ml++;
			if(ml >= list->size()){
				break;
			}
		}

		if(ml >= list->size()){
			temp.at(i)->setEnabled(false);
			continue;
		}
		if((*list).at(ml)->getChrom() < temp.at(i)->getChrom()){
			ml++;
		}
		if(ml >= list->size()){
			temp.at(i)->setEnabled(false);
			continue;
		}
		if(temp.at(i)->getChrom() == (*list).at(ml)->getChrom()){
			temp.at(i)->setEnabled(true);
		}
		else{
			temp.at(i)->setEnabled(false);
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
	Marker* min = (*list).at(0);
	Marker* max = (*list).at(1);

	vector<Marker*> temp = *marks;
	stable_sort(temp.begin(), temp.end(), less<Methods::Marker*>());

	for(unsigned int i = 0; i < temp.size(); i++){
		if(min != NULL){
			if(temp.at(i)->getBPLOC() < min->getBPLOC()){
				temp.at(i)->setEnabled(false);
			}
		}
		if(max != NULL){
			if(temp.at(i)->getBPLOC() > max->getBPLOC()){
				temp.at(i)->setEnabled(false);
			}
		}
	}
}

bool InputFilter::IncludeCovariateFilter(string cov, vector<string>* list){
	if(list == NULL){
		return true;
	}
	for(unsigned int i = 0; i < list->size(); i++){
		if(cov == (*list).at(i)){
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
		if(cov == (*list).at(i)){
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
		if(trait == (*list).at(i)){
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
		if(trait == (*list).at(i)){
			return false;
		}
	}
	return true;
}

void InputFilter::run_locus_filter(int f, vector<Marker*>* marks){
	marker_filters.at(f)(marks, marker_lists.at(f));
}

bool InputFilter::run_sample_filter(int f, Sample* s){
	return sample_filters.at(f)(s, sample_lists.at(f));
}

bool InputFilter::run_family_filter(int f, Family* fam){
	return family_filters.at(f)(fam, family_lists.at(f));
}

bool InputFilter::run_covariate_filter(int f, string s){
	return cov_filters.at(f)(s, cov_lists.at(f));
}

bool InputFilter::run_trait_filter(int f, string s){
	return trait_filters.at(f)(s, trait_lists.at(f));
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
		s_iter = find_if(samples->begin(), samples->end(), FindSampleByFamAndID(filter_samples.at(s)->getFamID(), filter_samples.at(s)->getInd()));
		if(s_iter != samples->end()){
			Sample* found = (*s_iter);
			vector<Marker*> range = filter_ranges.at(s);
			Marker* min = range.at(0);
			Marker* max = range.at(1);

			for(int m = 0; m < msize; m++){
				Marker* mark = (*markers).at(m);
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
