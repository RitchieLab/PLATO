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

namespace Methods {

class DataSet;
class Marker;
class Sample;

class DataLoader {

public:
	DataLoader();

	void read();

private:
	void readPed(const std::string& fn);
	void readMap(const std::string& fn);
	void readBinPed(const std::string& fn);
	void readTPed(const std::string& fn);

	Marker* parseMap(std::stringstream& ss);
	void parseSample(Marker* m, Sample* samp, const std::string& s1, const std::string& s2)


	// does the map file contain a distance column?
	bool _map_distance = true;
	// does the map file contain a referent allele column?
	bool _map_ref = false;
	// does the map file contain an alternate allele column?
	bool _map_alt = false;

	// does the ped file contain genotype information (false when reading fam)?
	bool _ped_genotype = true;
	// does the ped file have fids?
	bool _ped_fid = true;
	// does the ped file have parental ids?
	bool _ped_parents = true;
	// does the ped file have gender?
	bool _ped_gender = true;
	// does the ped file have phenotype?
	bool _ped_pheno = true;
	// is the control value in the phenotype "1" (-> case is "2")
	// or is the control "0" (-> case is "1")
	bool _ped_control1 = true

	// bitset dictating the inclusion/exclusion of markers in the dataset
	boost::dynamic_bitset _marker_incl;

	// bitset indicating inclusion/exclusion of samples in the dataset
	boost::dynamic_bitset _sample_incl;

	std::string _ped_missing_geno = "0";

	// the actual DataSet we're working with
	DataSet* ds_ptr;



};

}

#endif /* DATALOADER_H_ */
