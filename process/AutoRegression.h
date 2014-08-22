#ifndef PROCESSLIB_AUTOREGRESSION_H
#define PROCESSLIB_AUTOREGRESSION_H

#include "Process.h"
#include "analysis/Regression.h"
#include "MPIProcess.h"

#include <boost/array.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/serialization.hpp>

#include <utility>

namespace PLATO{
namespace ProcessLib{

class AutoRegression : public ProcessImpl<AutoRegression>,
	public Analysis::Regression, public MPIProcessImpl<AutoRegression>,
	private virtual LinearRegression, private virtual LogisticRegression {
private:
	const static std::string stepname;
	const static std::string MPIname;

public:
	AutoRegression() : ProcessImpl<AutoRegression>(stepname, "Run Auto Regression"),
		MPIProcessImpl<AutoRegression>(MPIname), class_data(0) {};
	virtual ~AutoRegression(){if(class_data){delete class_data;}};

	virtual void parseOptions(const boost::program_options::variables_map& vm);

protected:
	virtual void process(Data::DataSet&);
	virtual boost::program_options::options_description& appendOptions(boost::program_options::options_description& opts);

	virtual calc_fn& getCalcFn() const;
	virtual const Regression::ExtraData* getExtraData() const;

	virtual void printVarHeader(const std::string& var_name, std::ofstream& of) const;
	virtual bool initData();

	virtual void printExtraHeader(std::ofstream& of);
	virtual std::string printExtraResults(const Result& r);

private:
	static boost::array<double, 4> linkFunction(double v);
	static std::pair<float, float> calcPVal(Result* r, Result* submodel, unsigned int df, float null_ll);

	template<class Archive>
	static void* loadExtraData(const Archive& ar);

public:
	static void calculate_MPI(unsigned int bufsz, const char* buf, 
		std::deque<std::pair<unsigned int, const char*> >& result_queue, boost::mutex& result_mutex, boost::condition_variable& cv);

private:
	bool show_odds;
	unsigned int maxIterations;

public:
	class ExtraData : public LogisticRegression::ExtraData {
	public:

		unsigned short analysis_type;

		ExtraData(unsigned int n=0) : LogisticRegression::ExtraData(n) {}
		ExtraData(const Regression::ExtraData& o) : LogisticRegression::ExtraData(o) {}
		virtual ~ExtraData() {}

		template<class Archive>
		void serialize(Archive& ar, const unsigned int){
			ar & boost::serialization::base_object<LogisticRegression::ExtraData>(*this);
			ar & analysis_type;
		}
	};

private:
	mutable ExtraData* class_data;

};

}
}

#endif