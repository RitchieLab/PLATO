/*
 * MCAnnealParameters.h
 *
 *  Created on: Nov 3, 2009
 *      Author: gilesjt
 */

#ifndef MCANNEALPARAMETERS_H_
#define MCANNEALPARAMETERS_H_

namespace Methods {

class MCAnnealParameters {
public:
	MCAnnealParameters();
	MCAnnealParameters(int nb, int ni, int hp, int u, int o);
	virtual ~MCAnnealParameters();

	void setNburn(int s){nburn = s;}
	int getNburn(){return nburn;}
	void setNiter(int s){niter = s;}
	int getNiter(){return niter;}
	void setHyperPars(int s){hyperpars = s;}
	int getHyperPars(){return hyperpars;}
	void setUpdate(int s){update = s;}
	int getUpdate(){return update;}
	void setOutput(int s){output = s;}
	int getOutput(){return output;}

private:
	int nburn; //number of burnins
	int niter; //nuber of iters
	int hyperpars; //hyperparameters????????
	int update; //stdout update frequency
	int output; //output settings >1 bivariate, >2 trivariate, > 0, all, < 0 none
};

}

#endif /* MCANNEALPARAMETERS_H_ */
