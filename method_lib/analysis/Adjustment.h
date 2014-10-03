#ifndef ANALYSIS_ADJUSTMENT_H
#define ANALYSIS_ADJUSTMENT_H

#include <vector>
#include <iostream>
#include <gsl/gsl_randist.h>

namespace PLATO{
namespace Analysis{

class AdjustmentModel;

class Adjustment {

protected:
	// Class used to erase the type of the container being passed in
	class Wrapper{
		struct WrapData{
			virtual ~WrapData() {}
			virtual double at(unsigned int) const = 0;
			virtual unsigned int size() const = 0;
		};

		template <typename T>
		struct WrapModel : WrapData {
			WrapModel(const T& t) : data(t) {}
			virtual ~WrapModel() {}
			virtual double at(unsigned int i) const{
				return data[i];
			}
			virtual unsigned int size() const{
				return data.size();
			}

		private:
			const T& data;
		};

		WrapData* wc;

	public:
		template<class T> Wrapper(const T& t) : wc( new WrapModel<T>(t) ) {}
		~Wrapper() {delete wc;}

		double operator[](unsigned int i) const{
			return wc->at(i);
		}

		unsigned int size() const{
			return wc->size();
		}
	};


public:
	enum adjust_ENUM { REGRESSION, GC };

protected:
	Adjustment(){}

public:
	virtual ~Adjustment(){};

	virtual double adjust(Wrapper cont) const = 0;
	virtual double correct(double pval, double gif_recip) const = 0;

	static Adjustment* getAdjustmentMethod(const AdjustmentModel& c);
	static std::string listAdjustmentMethods();
};

class AdjustmentModel{
public:
	AdjustmentModel(const std::string& s);
	AdjustmentModel(const char* s);
	AdjustmentModel() : _data(Adjustment::GC) {}
	AdjustmentModel(Adjustment::adjust_ENUM c) : _data(c) {}

	operator int() const{return _data;}

//	bool operator<(const CorrectionModel& o) {return static_cast<int>(_data) < static_cast<int>(o._data);}

private:
	Adjustment::adjust_ENUM _data;
};


class RegressionAdjustment : public Adjustment {

public:
	virtual double adjust(Wrapper cont) const;
	virtual double correct(double pval, double gif_recip) const;

};

class GCAdjustment : public Adjustment {
public:
	virtual double adjust(Wrapper cont) const;
	virtual double correct(double pval, double gif_recip) const;
};

}
}

namespace std{
istream& operator>>(istream& in, PLATO::Analysis::AdjustmentModel& model_out);
ostream& operator<<(ostream& o, const PLATO::Analysis::AdjustmentModel& m);
}


#endif
