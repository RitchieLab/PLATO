#ifndef UTIL_MPIUTILS_H
#define UTIL_MPIUTILS_H

#include <sstream>
#include <utility>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

namespace PLATO{
namespace Utility{

class MPIUtils{

public:
	template<class T>
	static std::pair<unsigned int, const char *> pack(const T& data);

	template<class T>
	static void unpack(unsigned int, const char*, T& data);

};

template <class T>
std::pair<unsigned int, const char*> MPIUtils::pack(const T& data){
	std::stringstream ss_out;
	boost::archive::binary_oarchive oa(ss_out);

	oa << data;
	unsigned int bufsz = ss_out.tellp();
	ss_out.seekg(0);
	char* buf_out = new char[bufsz];
	ss_out.read(buf_out,bufsz);

	return std::pair<unsigned int, const char*>(bufsz, buf_out);
}

template <class T>
void MPIUtils::unpack(unsigned int bufsz, const char* in_buf, T& data){
	std::stringstream ss_resp;
	ss_resp.write(in_buf, bufsz);
	ss_resp.seekg(0);
	boost::archive::binary_iarchive ia_resp(ss_resp);
	ia_resp >> data;
}

}
}


#endif
