#ifndef GENOTYPE_READER_H
#define GENOTYPE_READER_H

#include "Reader.h"

#include <string>
#include <map>
#include <boost/program_options.hpp>

#include "GenotypeReaderFactory.h"

#include "DataSet.h"

/*!
 * Class for reading genotypes from a given file
 */
class GenotypeReader : public Reader{

public:
	GenotypeReader(const boost::program_options::variables_map& vm);
	virtual ~GenotypeReader() = 0;

	virtual DataSet& read(DataSet& ds) = 0;
	virtual boost::program_options::options_description& addOptions(boost::program_options::options_description& opts) const = 0;

	/*!
	 * Adds options common to every genotype reader (such as start / stop
	 * position, chromosome, etc.)
	 */
	static boost::program_options::options_description& addOptions(boost::program_options::options_description& opts);

protected:

	// Some helpers used by more than one reader
	DataSet& readMap(const std::string& fn, DataSet& ds) const;
	DataSet& readFam(const std::string& fn, DataSet& ds) const;
	DataSet& readBim(const std::string& fn, DataSet& ds) const;

protected:

	//! A set of regions to include. Note: when used in the read() function,
	// this MUST be in sorted order (the data structure is naturally sorted)
	std::vector< std::pair<unsigned short, std::pair<unsigned int, unsigned int> > > _region_inc;
	//! family include/exclude lists
	std::set<std::string> _fam_inc;
	std::set<std::string> _fam_exc;
	//! sample include/exclue lists
	std::set<std::string> _sample_inc;
	std::set<std::string> _sample_exc;
	//! marker incldue/exclude lists
	std::set<std::string> _marker_inc;
	std::set<std::string> _marker_exc;

	// Vectors of booleans for if a sample or a marker should be read into memory
	// If empty, we can assume no filtering, and these will be in the same order
	// as the genotype file.
	std::vector<bool> _maker_read;
	std::vector<bool> _sample_read;

};


template <class T>
class GenotypeReaderImpl : public GenotypeReader {

public:
	GenotypeReaderImpl(const boost::program_options::variables_map& vm) : GenotypeReader(vm) {};
	virtual ~GenotypeReaderImpl(){}

	static GenotypeReader* create(const boost::program_options::variables_map& vm){return new T(vm);}

protected:
	static const boost::program_options::options_description& doRegister(const boost::program_options::options_description& opts);

};

template <class T>
const boost::program_options::options_description& GenotypeReaderImpl<T>::doRegister(const boost::program_options::options_description& opts){
	unsigned int n_options = opts.options().size();
	for (unsigned int i = 0; i < n_options; i++){
		GenotypeReaderFactory::getFactory().RegisterReader(opts.options()[i].key(), &T::create);
	}
	return opts;
}
#endif
