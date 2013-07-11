#ifndef GENOTYPE_READER_FACTORY_H
#define GENOTYPE_READER_FACTORY_H

#include <map>
#include <string>
#include <boost/program_options.hpp>

class GenotypeReader;

typedef GenotypeReader* (createGenoReader)(const boost::program_options::variables_map&);

class GenotypeReaderFactory{
public:
	typedef std::map<const std::string, createGenoReader*>::const_iterator const_iterator;

private:
	GenotypeReaderFactory(){}
	GenotypeReaderFactory(const GenotypeReaderFactory&);
	GenotypeReaderFactory& operator=(const GenotypeReaderFactory&);

public:
	void RegisterReader(const std::string& key, createFunc* ptr);
	GenotypeReader* Create(const boost::program_options::variables_map& vm);

	const_iterator begin() const{return creation_map.begin();}
	const_iterator end() const{return creation_map.end();}

	static GenotypeReaderFactory& getFactory(){static GenotypeReaderFactory f; return f;}
	static boost::program_options::options_description& addOptions(boost::program_options::options_description& opts);

private:
	std::map<const std::string, createGenoReader*> creation_map;
	boost::program_options::options_description combined_options;
};

#endif
