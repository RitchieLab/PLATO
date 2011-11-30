#ifndef HELPERS_H
#define HELPERS_H

#include <string>
#include <sstream>
#include <fstream>
#include <bitset>
#include <vector>
#include "Family.h"
#include "Sample.h"
#include "Marker.h"
#include "DataSet.h"
#include "StepOptions.h"
#include "InputFilter.h"
//#include <boost/dynamic_bitset.hpp>

//using namespace boost;
using namespace std;
namespace std{

template<>
class less<Methods::Marker*>{
	public:
		bool operator()(Methods::Marker const* p1, Methods::Marker const* p2) const{
			if(!p1){
				return true;
			}
			if(!p2){
				return false;
			}
			return (p1->chrom < p2->chrom || (p1->chrom == p2->chrom && p1->bploc < p2->bploc));
		}
};

template<>
class greater<Methods::Marker*>{
public:
	bool operator()(Methods::Marker const* p1, Methods::Marker const* p2) const{
		if(!p1){
			return true;
		}
		if(!p2){
			return false;
		}
		return (p1->rsid < p2->rsid);
	}
};


template <class T>
inline const std::string getString(const T& t){
	    std::stringstream os;
		os << t;
		return os.str();
};


template<class T>
inline const T SQR(const T a) {return a*a;}

template<class T>
inline const T MAX(const T &a, const T &b)
	        {return b > a ? (b) : (a);}

			template<class T>
			inline const T MIN(const T &a, const T &b)
	        {return b < a ? (b) : (a);}

			template<class T>
			inline const T SIGN(const T &a, const T &b)
	        {return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

			template<class T>
			inline void SWAP(T &a, T &b)
	        {T dum=a; a=b; b=dum;}
};//end namespace std

namespace Methods{

typedef vector<vector<int> > table_t;
typedef vector<vector<double> > matrix_t;
typedef vector<double> vector_t;
typedef vector<bool> boolvec_t;
typedef vector<int> intvec_t;
// typedef vector<Individual*>::iterator iIndividual;
// typedef vector<Locus*>::iterator iLocus;
// typedef vector<CSNP*>::iterator iSNP;
typedef vector<bool>::iterator iAllele;


class Family;
class Sample;
class Marker;

class FindFamily{
	public:
		FindFamily(string n) : value(n){}

		bool operator() (Methods::Family* element) const{
			return (element->getFamID() == value);
		}

	private:
		string value;
};

class FindSampleByID{
	public:
		FindSampleByID(string n) : value(n){}

		bool operator() (Methods::Sample* element) const{
			return(element->getInd() == value);
		}
	private:
		string value;
};

class FindSampleByFamAndID{
	public:
		FindSampleByFamAndID(string f, string i) : fam(f), ind(i){}

		bool operator() (Methods::Sample* element) const{
			return(element->getFamID() == fam && element->getInd() == ind);
		}
	private:
		string fam;
		string ind;
};

class FindMarker{
	public:
		FindMarker(string m) : marker(m){}

		bool operator() (Methods::Marker* mark) const{
			return(mark->getProbeID() == marker);
		}
	private:
		string marker;
};

class FindMarkerByChromBploc{
	public:
		FindMarkerByChromBploc(int chr, int bp) : chrom(chr), bploc(bp){}
		bool operator()(Methods::Marker* mark) const{
			return (mark->getChrom() == chrom && mark->getBPLOC() == bploc);
		}
	private:
		int chrom;
		int bploc;
};

class FindString
{
public:
	FindString(string s) : otherString(s){}

	bool operator() (string str) const
	{
		return str == otherString;
	}
private:
	string otherString;
};


#define Abs(x)    ((x) < 0 ? -(x) : (x))
#define Max(a, b) ((a) > (b) ? (a) : (b))
#define TOLERANCE 0.00000000000001
#define FTOLERANCE 0.00000001
#define LOW 0.02425
#define HIGH 0.97575

class Helpers{
	public:
		//stats
		static const double a[];/* =
		  {
			      -3.969683028665376e+01,
				      2.209460984245205e+02,
					      -2.759285104469687e+02,
						      1.383577518672690e+02,
							      -3.066479806614716e+01,
								       2.506628277459239e+00
										     };
*/
		static const double b[];/* =
		  {
			      -5.447609879822406e+01,
				      1.615858368580409e+02,
					      -1.556989798598866e+02,
						      6.680131188771972e+01,
							      -1.328068155288572e+01
									    };
*/
		static const double c[];/* =
		  {
			      -7.784894002430293e-03,
				      -3.223964580411365e-01,
					      -2.400758277161838e+00,
						      -2.549732539343734e+00,
							      4.374664141464968e+00,
								       2.938163982698783e+00
										     };
										     */
		static const double d[];/* =
		  {
			      7.784695709041462e-03,
				      3.224671290700398e-01,
					      2.445134137142996e+00,
						      3.754408661907416e+00
								    };
*/

static double DoubDif(double a, double b);

static float FloatDif(float a, float b);

static bool float_comp(float a, float b);

static bool double_comp(double a, double b);

static bool fEquals(float a, float b);

static bool dEquals(double a, double b);

static bool fLess(float a, float b);

static bool dLess(double a, double b);

static bool fLessOrEqual(float a, float b);

static bool dLessOrEqual(double a, double b);

static bool fGreater(float a, float b);

static bool dGreater(double a, double b);

static bool fGreaterOrEqual(float a, float b);

static bool dGreaterOrEqual(double a, double b);

static bool fileExists(const string& fileName);

static bool isAlphaNum( string s );

static double p_from_chi(double chi, double df);

static void printFamsToDigit(vector<Methods::Family*>* families, string name, StepOptions options);

static void remapFamsToDigit(vector<Methods::Family*>* families);

static bool isValidMarker(Methods::Marker* mark, StepOptions* options, int& prev_base, int& prev_chrom);

static vector<Methods::Marker*> findValidMarkers(vector<Methods::Marker*>* marks, StepOptions options);

static vector<Methods::Marker*> findValidMarkers(vector<Methods::Marker*>* marks, StepOptions* options);

static vector<int> findValidMarkersIndexes(vector<Methods::Marker*>* marks, StepOptions options);

static vector<int> findValidMarkersIndexes(vector<Methods::Marker*>* marks, StepOptions* options);

static vector<vector<Sample*> > generateSampleSets(DataSet* ds, StepOptions* options);

static vector<Methods::Marker*> findRandomMarkers(vector<Methods::Marker*> good_markers, vector<Methods::Marker*>* used_markers, StepOptions* options);

static void readCovariateFile(string file, DataSet* ds, StepOptions options, InputFilter* filters);

static void readTraitFile(string file, DataSet* ds, StepOptions options, InputFilter* filters);

//Added 02-23-2011 to support -update-ids option
static void readIDFile(string file, DataSet* ds);

static void readCovariates(string file, vector<Methods::Sample*>* samples, vector<string>* covheaders);

static bool readString(FILE* fp, string* s);

static void readZeroGenoFile(string file);

static void zeroSingleGenos(vector<Methods::Marker*>* markers, vector<Methods::Sample*>* samples);

static void assignLinks(vector<Methods::Family*>* families);

static void reorderAlleles(vector<Methods::Sample*>* samples, vector<Methods::Marker*>* markers);

static string map_allele(Methods::Marker* mark, string a, StepOptions *options);

static map<string, string> readCustomAlleles(string f);

static void remapSamples(vector<Methods::Sample*>* samples, vector<Methods::Marker*>* markers, vector<int>* marker_map, int loc);

static void readPedM(vector<Methods::Sample*>* samples, vector<Methods::Family*>* families, vector<Methods::Marker*>* markers, vector<int>* marker_map, StepOptions options);

static string removeBeginWhiteSpace(string l);

static void readSampleFile(string file, vector<Methods::Sample*>* mlist);

static void readFamilyFile(string file, vector<Methods::Family*>* mlist);

static void readCovTraitFile(string file, vector<string>* list);

static void readLocusFile(string file, vector<Methods::Marker*>* mlist);

static map<string, float> readFreqFile(string file);

static void readMapM(DataSet* ds, StepOptions options,  InputFilter* filters);

static void readMapM(DataSet* ds, StepOptions options);

static void readMapM(vector<Methods::Marker*>* markers, vector<int>* marker_map, StepOptions options);

static void readTPedM(DataSet* ds, StepOptions options, InputFilter* filters);

static void readTPedM(DataSet* ds, StepOptions options);

static void readTFamM(vector<Methods::Sample*>* samples, vector<Methods::Family*>* families, StepOptions options, InputFilter* filters);

static void readTFamM(vector<Methods::Sample*>* samples, vector<Methods::Family*>* families, StepOptions options);

static void readTFamM(DataSet* ds, StepOptions options, InputFilter* filters);

static void readTFamM(DataSet* ds, StepOptions options);

static void readBinM(vector<Methods::Sample*>* samples, vector<Methods::Family*>* families, vector<Methods::Marker*>* markers, vector<int>* marker_map, StepOptions options, InputFilter* filters);

static void readBinM(vector<Methods::Sample*>* samples, vector<Methods::Family*>* families, vector<Methods::Marker*>* markers, vector<int>* marker_map, StepOptions options);

static void readBinM(DataSet* ds, StepOptions options, InputFilter* filters);

static void readBinM(DataSet* ds, StepOptions options);

static void readPedInfo();

static void readPedM_3vec(vector<Methods::Sample*>* samples, vector<Methods::Family*>* families, vector<Methods::Marker*>* markers, vector<int>* marker_map, StepOptions options);

static void readMapMdr(DataSet* ds, StepOptions options,  InputFilter* filters);

static void readMdr(DataSet* set, StepOptions options, InputFilter* filters);

static void readLgenFile(DataSet* set, StepOptions options, InputFilter* filters);

static void readReferenceFile(DataSet* set, StepOptions options, InputFilter* filters);

static double ltqnorm(double p);

static void readPedM_3vec_set(DataSet* set, StepOptions options, InputFilter* filters);

static void readPedM_3vec_set(DataSet* set, StepOptions options);

static bool realnum(double d);

static double normdist(double z);

static double pT(double T, double df);

static double chiprobP(double chi, double df);

static double SQR(double a);

static double pythag(const double a, const double b);

static void svdcmp(vector<vector<double> > & a,
		        vector<double> & w,
				        vector<vector<double> > &v);

static vector< vector<double> > svd_inverse(vector< vector<double> > & u);

static void sizeMatrix(vector<vector<double> > &m, int r, int c);

static void svbksb(vector<vector<double> > &u, vector<double> &w, vector<vector<double> > &v,
		vector<double> &b, vector<double> &x);

static void multMatrix(vector<vector<double> > & a,
        vector<vector<double> > & b,
        vector<vector<double> > & c);

static string stringToLowerCase(string str);



};

};

#endif
