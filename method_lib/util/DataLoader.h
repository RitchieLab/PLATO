/*
 * DataLoader.h
 *
 *  Created on: Nov 11, 2013
 *      Author: jrw32
 */

#ifndef METHODS_DATA_LOADER_H
#define METHODS_DATA_LOADER_H

#include <string>
#include <sstream>

#include <boost/dynamic_bitset.hpp>
#include <boost/program_options.hpp>

namespace Methods {

class DataSet;
class Marker;
class Sample;

class DataLoader {

public:
	enum input_style{
		UNKNOWN,
		PED,
		BED,
		TPED
	};

	DataLoader();
	virtual ~DataLoader() {}

	boost::program_options::options_description& addOptions(boost::program_options::options_description& opts);
	void parseOptions(const boost::program_options::variables_map& vm);

protected:
	void setDataSet(DataSet& ds){ds_ptr = &ds;}
	void read();

private:
	void readPed(const std::string& fn);
	void readMap(const std::string& fn);
	void readBinPed(const std::string& fn);
	void readTPed(const std::string& fn);

	Marker* parseMap(std::stringstream& ss);
	void parseSample(Marker* m, Sample* samp, const std::string& s1, const std::string& s2);


	// does the map file contain a distance column?
	bool _map_no_distance;
	// does the map file contain a referent allele column?
	bool _map_ref;
	// does the map file contain an alternate allele column?
	bool _map_alt;

	// does the ped file contain genotype information (false when reading fam)?
	bool _ped_genotype;
	// does the ped file have fids?
	bool _ped_no_fid;
	// does the ped file have parental ids?
	bool _ped_no_parents;
	// does the ped file have gender?
	bool _ped_no_gender;
	// does the ped file have phenotype?
	bool _ped_no_pheno;
	// is the control value in the phenotype "1" (-> case is "2")
	// or is the control "0" (-> case is "1")
	bool _ped_control0;

	// bitset dictating the inclusion/exclusion of markers in the dataset
	boost::dynamic_bitset<> _marker_incl;

	// bitset indicating inclusion/exclusion of samples in the dataset
	boost::dynamic_bitset<> _sample_incl;

	std::string _ped_missing_geno;

	// the actual DataSet we're working with
	DataSet* ds_ptr;

	enum input_style input;

	//options
	std::string file_base;
	std::string ped_fn;
	std::string map_fn;
	std::string bfile_base;
	std::string bim_fn;
	std::string bed_fn;
	std::string fam_fn;
	std::string tfile_base;
	std::string tped_fn;
	std::string tfam_fn;

};

}

#endif /* DATALOADER_H_ */
