#ifndef PROCESSLIB_LOGISTICREGRESSION_H
#define PROCESSLIB_LOGISTICREGRESSION_H

#include "Process.h"
#include "analysis/Regression.h"
#include "MPIProcess.h"

#include <boost/array.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/export.hpp>

#include <utility>

namespace PLATO{
namespace ProcessLib{

class LogisticRegression : public ProcessImpl<LogisticRegression>,
	public virtual Analysis::Regression, public MPIProcessImpl<LogisticRegression> {

	friend class AutoRegression;

private:
	const static std::string stepname;
	const static std::string MPIname;

public:
	LogisticRegression() : ProcessImpl<LogisticRegression>(stepname, "Run Logistic Regression"),
		MPIProcessImpl<LogisticRegression>(MPIname), class_data(0) {};
	virtual ~LogisticRegression(){if(class_data){delete class_data;}};

	virtual void parseOptions(const boost::program_options::variables_map& vm);

protected:
	virtual void process(Data::DataSet&);
	virtual boost::program_options::options_description& appendOptions(boost::program_options::options_description& opts);
	virtual boost::program_options::options_description getExtraOptions();

	virtual calc_fn& getCalcFn() const;
	virtual const Regression::ExtraData* getExtraData() const;

	static Result* calculate(const double* Y, const double* data,
			unsigned int n_cols, unsigned int n_rows, unsigned int offset,
			unsigned int n_covars, bool run_null, const Regression::ExtraData* other_data);

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
	class LogisticData : public Regression::ExtraData {
	public:
		bool show_odds;
		unsigned int maxIterations;

		LogisticData(unsigned int n=0) : Regression::ExtraData(n) {}
		LogisticData(const Regression::ExtraData& o) : Regression::ExtraData(o) {
			const LogisticData* lo = dynamic_cast<const LogisticData*>(&o);
			if(lo){
				show_odds = lo->show_odds;
				maxIterations = lo->maxIterations;
			}
		}
		virtual ~LogisticData() {}

		template<class Archive>
		void serialize(Archive& ar, const unsigned int){
			ar & boost::serialization::base_object<Regression::ExtraData>(*this);
			ar & show_odds;
			ar & maxIterations;
		}
	};

private:
	mutable LogisticData* class_data;

};

}
}

BOOST_CLASS_EXPORT_KEY(PLATO::ProcessLib::LogisticRegression::LogisticData)

#endif
