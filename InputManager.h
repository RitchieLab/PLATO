#ifndef INPUT_MANAGER_H
#define INPUT_MANAGER_H

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <string>
#include <map>

/*!
 * Class that manages the input to PLATO
 */
class InputManager{

private:
	template<typename T>
	struct ci_less:std::binary_function<T,T,bool>
	  { bool operator() (const T& s1,const T& s2) const { return boost::ilexicographical_compare(s1,s2); }};

	typedef std::map<std::string, unsigned short, ci_less<std::string> > IChromMap;

public:

	/*!
	 * Converts a string into a chromosome number (i.e. X -> 23, Y -> 24, etc)
	 */
	static unsigned short chrStringToInt(const std::string& chr_str){
		IChromMap::const_iterator itr = s_chr_map.find(chr_str);
		return itr == s_chr_map.end() ? 0 : (*itr).second;
	}

	static const std::string& chrIntToString(unsigned int chr_int){
		std::map<unsigned short, std::string>::const_iterator itr = s_chrint_map.find(static_cast<unsigned short>(chr_int));
		return itr == s_chrint_map.end() ? s_missing_chr_str : (*itr).second;
	}

	/*!
	 * Adds some options relevant to parsing input
	 */
	static boost::program_options::options_description& addOptions(boost::program_options::options_description& opts);

	/*!
	 * Parse options related to the chromosome structure you're working with
	 */
	static void parseGlobalOptions(const boost::program_options::variables_map& vm);

private:

	//! A mapping of strings -> chromosome numbers
	static IChromMap s_chr_map;
	static std::map<unsigned short, std::string> s_chrint_map;
	static const std::string s_missing_chr_str;

};

#endif
