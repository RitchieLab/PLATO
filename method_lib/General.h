#ifndef GENERAL_H
#define GENERAL_H

#include <string>
#include <sstream>
#include <fstream>
#include <bitset>
#include <vector>
//#include <boost/dynamic_bitset.hpp>

//using namespace boost;
using namespace std;
namespace Methods{
class General{
	public:
		static void Tokenize(const string&, vector<string>&, const string&);
		static vector<string> ParseDelimitedLine(string);
		static vector<string> ParseDelimitedLine(string line, string delim);

		static string getStringInt(int);
		static string getStringFloat(float);
//		static string getString(int);
//		static string getString(float);
//		static string getString(double);
//		static string getString(unsigned int);
};
};
#endif
