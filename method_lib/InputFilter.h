#ifndef INPUTFILTER_H
#define INPUTFILTER_H

#include <vector>
#include "Family.h"
#include "Sample.h"
#include "Marker.h"
#include "MethodException.h"

using namespace std;
namespace Methods{
class InputFilter {
	private:
		vector<bool (*)(Sample*, vector<Sample*>*)> sample_filters;
		vector<bool (*)(Family*, vector<Family*>*)> family_filters;
//		vector<bool (*)(Marker*, vector<Marker*>*)> marker_filters;
		vector<void(*)(vector<Marker*>*, vector<Marker*>*)> marker_filters;
		vector<bool (*)(string, vector<string>*)> cov_filters;
		vector<bool (*)(string, vector<string>*)> trait_filters;

		vector<vector<Sample*>*> sample_lists;
		vector<vector<Family*>*> family_lists;
		vector<vector<Marker*>*> marker_lists;
		vector<vector<string>*> cov_lists;
		vector<vector<string>*> trait_lists;

	public:

		InputFilter(){}
		~InputFilter(){}

		///Predefined filters
		static void IncludeLocusFilter(vector<Marker*>* marks, vector<Marker*>* mlist);
//		static bool IncludeLocusFilter(Marker*, vector<Marker*>*);
		static void ExcludeLocusFilter(vector<Marker*>* marks, vector<Marker*>*);
		static bool IncludeSampleFilter(Sample*, vector<Sample*>*);
		static bool ExcludeSampleFilter(Sample*, vector<Sample*>*);
		static bool IncludeFamilyFilter(Family*, vector<Family*>*);
		static bool ExcludeFamilyFilter(Family*, vector<Family*>*);
		static void LocusChromFilter(vector<Marker*>*, vector<Marker*>*);
		static void LocusBplocRangeFilter(vector<Marker*>*, vector<Marker*>*);
		static bool IncludeCovariateFilter(string, vector<string>*);
		static bool ExcludeCovariateFilter(string, vector<string>*);
		static bool IncludeTraitFilter(string, vector<string>*);
		static bool ExcludeTraitFilter(string, vector<string>*);

		void add_sample_filter(bool (*func)(Sample*, vector<Sample*>*)){sample_filters.push_back(func);}
		void add_family_filter(bool (*func)(Family*, vector<Family*>*)){family_filters.push_back(func);}
//		void add_locus_filter(bool (*func)(Marker*, vector<Marker*>*)){marker_filters.push_back(func);}
		void add_locus_filter(void (*func)(vector<Marker*>*, vector<Marker*>*)){marker_filters.push_back(func);}
		void add_covariate_filter(bool (*func)(string, vector<string>*)){cov_filters.push_back(func);}
		void add_trait_filter(bool (*func)(string, vector<string>*)){trait_filters.push_back(func);}

		void add_sample_list(vector<Sample*>* l){sample_lists.push_back(l);}
		void add_family_list(vector<Family*>* l){family_lists.push_back(l);}
		void add_locus_list(vector<Marker*>* l){marker_lists.push_back(l);}
		void add_covariate_list(vector<string>* l){cov_lists.push_back(l);}
		void add_trait_list(vector<string>* l){trait_lists.push_back(l);}

		vector<bool (*)(Sample*, vector<Sample*>*)>* get_sample_filters(){return &sample_filters;}
		vector<bool (*)(Family*, vector<Family*>*)>* get_family_filters(){return &family_filters;}
//		vector<bool (*)(Marker*, vector<Marker*>*)>* get_marker_filters(){return &marker_filters;}
		vector<void (*)(vector<Marker*>*, vector<Marker*>*)>* get_marker_filters(){return &marker_filters;}
		vector<bool (*)(string, vector<string>*)>* get_covariate_filters(){return &cov_filters;}
		vector<bool (*)(string, vector<string>*)>* get_trait_filters(){return &trait_filters;}

		vector<vector<Sample*>*>* get_sample_lists(){return &sample_lists;}
		vector<vector<Family*>*>* get_family_lists(){return &family_lists;}
		vector<vector<Marker*>*>* get_marker_lists(){return &marker_lists;}
		vector<vector<string>*>* get_covariate_lists(){return &cov_lists;}
		vector<vector<string>*>* get_trait_lists(){return &trait_lists;}

		int num_locus_filters(){return marker_filters.size();}
		int num_sample_filters(){return sample_filters.size();}
		int num_family_filters(){return family_filters.size();}
		int num_covariate_filters(){return cov_filters.size();}
		int num_trait_filters(){return trait_filters.size();}

		void run_locus_filter(int, vector<Marker*>*);
//		bool run_locus_filter(int, Marker*);
		bool run_sample_filter(int, Sample*);
		bool run_family_filter(int, Family*);
		bool run_covariate_filter(int, string);
		bool run_trait_filter(int, string);

		void run_sample_bprange_filter(vector<Sample*>*, vector<Marker*>*, vector<Sample*>, vector<vector<Marker*> >);
};
};
#endif
