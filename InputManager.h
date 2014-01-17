#ifndef INPUT_MANAGER_H
#define INPUT_MANAGER_H

#include <string>
#include <sstream>
#include <map>

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>

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

	/*!
	 * Parse a vector of separated strings
	 */
	template <class I, class O>
	static void parseInput(const I& cont_in, O& cont_out, const std::string& sep=",");

private:

	//! A mapping of strings -> chromosome numbers
	static IChromMap s_chr_map;
	static std::map<unsigned short, std::string> s_chrint_map;
	static const std::string s_missing_chr_str;

};

template <class I, class O>
void InputManager::parseInput(const I& cont_in, O& cont_out, const std::string& sep){
	typename I::const_iterator itr = cont_in.begin();
	while(itr != cont_in.end()){
		std::stringstream out_s;
		out_s << (*itr);

		boost::char_separator<char> tok_sep(sep.c_str());

		boost::tokenizer<boost::char_separator<char> > tok(out_s.str());
		boost::tokenizer<boost::char_separator<char> >::iterator t_itr = tok.begin();
		while(t_itr != tok.end()){
			cont_out.insert(cont_out.end(), boost::lexical_cast<typename O::value_type>(*t_itr));
			++t_itr;
		}

		++itr;
	}
}

#endif
