#ifndef METHODS_UTIL_CONTAINER_H
#define METHODS_UTIL_CONTAINER_H

#include <sstream>
#include <iostream>

#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

namespace Methods{
namespace Util{

template<class C, class T>
class Container {
public:
	explicit Container<C, T> (const std::string& str) {
		std::stringstream ss(str);
		ss >> (*this);
	}
	Container<T> () {
	}

	operator C() const { return _data; }

	C::const_iterator begin(){ return _data.begin();}
	C::const_iterator end(){ return _data.end();}

	void push_back(const T& val) {
		_data.insert(_data.end(), val);
	}
private:
	C _data;
};

}
}

namespace std{

template <class C, class T>
std::ostream& operator<<(std::ostream& o, const Methods::Util::Container<C, T>& d){
	std::string sep = ",";
	typename C::const_iterator d_itr = d.begin();

	int i=0;
	while(d_itr != data.end()){
		if(i){
			o << sep;
		}
		o << *d_itr;

		++d_itr;
		++i;
	}

	return o;
}

template <class T>
std::istream& operator>>(std::istream& in, Methods::Util::Container<T>& d_out){
	std::string in_s;
	in >> in_s;
	// Now, split up the string
	boost::tokenizer<boost::escaped_list_separator<char> > tok(in_s);
	for(boost::tokenizer<boost::escaped_list_separator<char> >::iterator beg=tok.begin(); beg!=tok.end();++beg){
		try{
			T obj(boost::lexical_cast<T>(*beg));
			d_out.push_back(obj);
		}catch(boost::bad_lexical_cast& e){
	    	throw boost::program_options::validation_error(boost::program_options::validation_error::invalid_option_value);
		}
	}

	return in;
}


#endif
