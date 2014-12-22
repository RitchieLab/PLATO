/*
 * DataLoader.h
 *
 *  Created on: Nov 11, 2013
 *      Author: jrw32
 */

#ifndef UTILITY_DATA_LOADER_H
#define UTILITY_DATA_LOADER_H

#include <string>

#include <bitset>

#include <boost/dynamic_bitset.hpp>
#include <boost/program_options.hpp>

namespace PLATO {

namespace Data{

class DataSet;
class Marker;
class Sample;
}

namespace Utility{
class DataLoader {

public:
	enum input_style{
		UNKNOWN,
		PED,
		BED,
		TPED,
		LGEN,
		BEAGLE,
		VCF
	};

	DataLoader();
	virtual ~DataLoader() {}

	boost::program_options::options_description& addOptions(boost::program_options::options_description& opts);
	void parseOptions(const boost::program_options::variables_map& vm);

protected:
	//void setDataSet(PLATO::Data::DataSet& ds){ds_ptr = &ds;}
	void read(Data::DataSet& ds);

private:
	void readPed(const std::string& fn);
	void readMap(const std::string& fn);
	void readBinPed(const std::string& fn);
	void readTPed(const std::string& fn);
	void readLGen(const std::string& fn);

	void readBeagle();

	void readVCF();

	PLATO::Data::Marker* parseMap(std::stringstream& ss) const;
	void parseSample(PLATO::Data::Marker* m, PLATO::Data::Sample* samp,
			const std::string& s1, const std::string& s2,
			bool append = true, bool phased=false) const;

	/*!
	 * \brief Determines wheter to use the given marker
	 * Enables a pre-filter for loading the given marker based solely on the
	 * chromosome, base pair and ID.  Returns true if the given marker passes
	 * all pre-filter steps.
	 */
	bool filterMarker(const std::string& chrom, int bploc, const std::string& id) const;

	/*!
	 * \brief Determines wheter to use the given marker
	 * Enables a pre-filter for loading the given marker based solely on the
	 * chromosome, base pair and ID.  Returns true if the given marker passes
	 * all pre-filter steps.
	 */
	bool filterSample(const std::string& id, const std::string& fid="") const;

	static void readSampleList(const std::vector<std::string>& in_list, std::set<std::string>& out_set);
	static void readSampleFile(const std::vector<std::string>& in_list, std::set<std::string>& out_set);
	static void readMarkerFile(const std::vector<std::string>& fn_list, std::set<std::string>& out_set);
	static void addSampleToSet(const std::string& samp, std::set<std::string>& out_set);

	// does the map file contain a distance column?
	bool _map_no_distance;
	// does the map file contain a referent allele column?
	bool _map_ref;
	// does the map file contain an alternate allele column?
	bool _map_alt;
	// does the map file contain even more alleles?
	bool _map_others;

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
	// is the phenotype value quantitative (float)?
	bool _quant;

	// bitset dictating the inclusion/exclusion of markers in the dataset
	boost::dynamic_bitset<> _marker_incl;

	// bitset indicating inclusion/exclusion of samples in the dataset
	boost::dynamic_bitset<> _sample_incl;

	std::string _ped_missing_geno;

	// the actual DataSet we're working with
	PLATO::Data::DataSet* ds_ptr;

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
	std::string lfile_base;
	std::string lgen_fn;

	// BEAGLE options
	std::string beagle_prefix;
	std::string beagle_suffix;
	std::string marker_suffix;

	std::string bgl_missing_allele;

	std::vector<std::string> beagle_fns;
	std::vector<std::string> marker_fns;
	std::vector<std::string> chroms;

	std::bitset<3> _fns_provided;

	bool bgl_trio;
	bool bgl_pair;
	bool bgl_phased;
	bool bgl_poly;

	std::string vcf_fn;
	bool no_filter_marker;
	bool no_filter_geno;
	bool vcf_phased;
	bool vcf_poly;

	std::vector<std::string> chrom_list;
	std::vector<std::string> bp_window;
	std::vector<std::string> incl_marker_str;
	std::vector<std::string> excl_marker_str;
	std::vector<std::string> incl_marker_fns;
	std::vector<std::string> excl_marker_fns;
	std::vector<std::string> incl_sample_str;
	std::vector<std::string> excl_sample_str;
	std::vector<std::string> incl_sample_fns;
	std::vector<std::string> excl_sample_fns;

	std::set<unsigned short> chrom_ids;
	std::set<std::pair<unsigned int, unsigned int> > bp_window_set;
	std::set<std::string> incl_sample_set;
	std::set<std::string> excl_sample_set;
	std::set<std::string> incl_marker_set;
	std::set<std::string> excl_marker_set;

	static const std::string sampl_field_sep;
};

}
}

#endif /* DATALOADER_H_ */
