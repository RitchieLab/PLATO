#ifndef PROCESSLIB_OUTPUTBEAGLE_H
#define PROCESSLIB_OUTPUTBEAGLE_H

#include <string>

#include "Process.h"

#include "data/DataSet.h"
#include "data/Family.h"

namespace PLATO{

namespace ProcessLib{

class OutputBeagle : public ProcessImpl<OutputBeagle>{
private:
	const static std::string stepname;

	/*! A sample generator that will iterate over the samples in the correct
	 *  order (pair / trio data as appropriate)
	 *
	 */
	class SampleGenerator{

	public:
		SampleGenerator(const Data::DataSet& ds) : _ds(ds) {}
		virtual ~SampleGenerator() {}
		virtual const Data::Sample* next() = 0;
		virtual void reset() = 0;

	protected:
		const Data::DataSet& _ds;
	};

	class BasicSampleGenerator : public SampleGenerator{
	public:
		BasicSampleGenerator(const Data::DataSet& ds) : SampleGenerator(ds), _si(ds.beginSample()) {}
		virtual ~BasicSampleGenerator() {}

		virtual const Data::Sample* next() { return (_si == _ds.endSample()) ? 0 : (*(_si++));}
		virtual void reset() { _si = _ds.beginSample();}

	private:
		Data::DataSet::const_sample_iterator _si;
	};

	class FamilySampleGenerator : public SampleGenerator{
	public:
		FamilySampleGenerator(const Data::DataSet& ds, unsigned int n_found);
		virtual ~FamilySampleGenerator() {}

		virtual const Data::Sample* next();
		virtual void reset() {_fi = _ds.beginFamily(); _curr_founder = 0;}

	private:
		const unsigned int _FOUNDERS;
		// This is the current founder (if == _FOUNDERS, we're on children!)
		unsigned int _curr_founder;

		Data::DataSet::const_family_iterator _fi;


	};

public:
	OutputBeagle() : ProcessImpl<OutputBeagle>(stepname, "Example Command - please remove in released version"), _sep("\t") {};
	virtual ~OutputBeagle(){};

	virtual void parseOptions(const boost::program_options::variables_map& vm);


protected:
	virtual void process(Data::DataSet&);
	virtual boost::program_options::options_description& appendOptions(boost::program_options::options_description& opts);

private:
	SampleGenerator* getSampleGenerator(const Data::DataSet& ds);

private:
	std::string _prefix;
	std::string _suffix;
	std::string _marker_suff;
	std::string _sep;
	bool _incl_trait;
	bool _pair;
	bool _trio;

};

}
}

#endif
