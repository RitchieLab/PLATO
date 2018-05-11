#ifndef PROCESSLIB_AUTOREGRESSION_H
#define PROCESSLIB_AUTOREGRESSION_H

#include "Process.h"
#include "analysis/Regression.h"
#include "LinearRegression.h"
#include "LogisticRegression.h"
#include "MPIProcess.h"

#include <boost/array.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/export.hpp>

#include <utility>

namespace PLATO{
namespace ProcessLib{

class AutoRegression : public ProcessImpl<AutoRegression>,
	public virtual Analysis::Regression, public MPIProcessImpl<AutoRegression>,
	private virtual LinearRegression, private virtual LogisticRegression{

private:
	const static std::string stepname;
	const static std::string MPIname;

public:
	AutoRegression() : LinearRegression(), LogisticRegression(), ProcessImpl<AutoRegression>(stepname, "Run Auto Regression"),
		MPIProcessImpl<AutoRegression>(MPIname), class_data(0) {};
	virtual ~AutoRegression(){if(class_data){delete class_data;}};

	virtual void parseOptions(const boost::program_options::variables_map& vm);
	using ProcessImpl<AutoRegression>::create;

protected:

	virtual const std::string& getMPIName() const{ return MPIProcessImpl<AutoRegression>::getMPIName(); }

	using MPIProcessImpl<AutoRegression>::registerMPI;
	using ProcessImpl<AutoRegression>::doRegister;

	virtual const std::string& getName() const{ return ProcessImpl<AutoRegression>::getName(); }
	virtual const std::string& getDesc() const{ return ProcessImpl<AutoRegression>::getDesc(); }

	virtual void process(Data::DataSet&);
	virtual boost::program_options::options_description& appendOptions(boost::program_options::options_description& opts);
	virtual boost::program_options::options_description getExtraOptions();

	virtual calc_fn& getCalcFn() const;
	virtual const Regression::ExtraData* getExtraData() const;

	virtual void printVarHeader(const std::string& var_name, std::ofstream& of) const;
	virtual bool initData();

	static Result* calculate(const double* Y, const double* data,
			unsigned int n_cols, unsigned int n_rows, unsigned int offset,
			unsigned int n_covars, bool run_null, const Regression::ExtraData* other_data);

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

public:

	enum AnalysisType {LINEAR, LOGISTIC};

	class AutoData : public LogisticRegression::LogisticData {
	public:

		AnalysisType analysis_type;

		AutoData(unsigned int n=0) : LogisticRegression::LogisticData(n) {}
		AutoData(const Regression::ExtraData& o) : LogisticRegression::LogisticData(o) {}
		virtual ~AutoData() {}

		template<class Archive>
		void serialize(Archive& ar, const unsigned int){
			ar & boost::serialization::base_object<LogisticRegression::LogisticData>(*this);
			ar & analysis_type;
		}
	};

private:
	mutable AutoData* class_data;
	AnalysisType analysis_type;

	std::vector<std::string> linear_traits;
	std::vector<std::string> logistic_traits;

	std::set<std::string> linear_set;
	std::set<std::string> logistic_set;

};

}
}

BOOST_CLASS_EXPORT_KEY(PLATO::ProcessLib::AutoRegression::AutoData)

#endif
